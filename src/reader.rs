use flate2::read::MultiGzDecoder;
use log::error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;
use std::sync::mpsc;
use std::thread;

/// Detect whether a file is gzip-compressed by reading the magic bytes (0x1f, 0x8b).
pub fn is_gzipped(path: &Path) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let mut buf = [0u8; 2];
    match file.read_exact(&mut buf) {
        Ok(_) => Ok(buf[0] == 0x1f && buf[1] == 0x8b),
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(e),
    }
}

/// Open a file that may or may not be gzip-compressed, returning a buffered reader.
/// Uses magic bytes detection rather than file extension.
pub fn open_maybe_gzipped(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let gzipped = is_gzipped(path)?;
    let file = File::open(path)?;
    if gzipped {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Spawn a thread that reads a gzipped fragment file in chunks and sends
/// them over a sync channel. Returns the join handle and receiver.
pub fn spawn_fragment_reader(
    frag_file: &Path,
) -> (
    thread::JoinHandle<()>,
    mpsc::Receiver<Vec<u8>>,
    mpsc::Sender<Vec<u8>>,
) {
    let (tx, rx) = mpsc::sync_channel(100);
    let (pool_tx, pool_rx) = mpsc::channel();
    let frag_file = frag_file.to_path_buf();

    let handle = thread::spawn(move || {
        let mut reader = match open_maybe_gzipped(&frag_file) {
            Ok(r) => r,
            Err(e) => {
                error!("Failed to open fragment file: {}", e);
                return;
            }
        };

        const CHUNK_LINES: usize = 10_000;
        let mut buffer = pool_rx
            .try_recv()
            .unwrap_or_else(|_| Vec::with_capacity(CHUNK_LINES * 100));
        let mut total_fragments = 0;
        let mut line_count = 0;

        loop {
            let before = buffer.len();
            match reader.read_until(b'\n', &mut buffer) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    // Skip comment lines
                    if buffer.get(before) == Some(&b'#') {
                        buffer.truncate(before);
                        continue;
                    }

                    line_count += 1;
                    total_fragments += 1;

                    if line_count >= CHUNK_LINES {
                        if total_fragments % 1_000_000 == 0 {
                            eprint!("\rProcessed {} M fragments", total_fragments / 1_000_000);
                            std::io::stderr().flush().expect("Can't flush output");
                        }

                        let mut chunk_to_send = pool_rx
                            .try_recv()
                            .unwrap_or_else(|_| Vec::with_capacity(CHUNK_LINES * 100));
                        chunk_to_send.clear();
                        std::mem::swap(&mut buffer, &mut chunk_to_send);

                        if tx.send(chunk_to_send).is_err() {
                            break;
                        }
                        line_count = 0;
                    }
                }
                Err(e) => {
                    error!("Error reading fragment file: {}", e);
                    break;
                }
            }
        }

        // Send any remaining fragments
        if !buffer.is_empty() {
            let _ = tx.send(buffer);
        }

        eprintln!("\nFinished reading {} total fragments", total_fragments);
    });

    (handle, rx, pool_tx)
}
