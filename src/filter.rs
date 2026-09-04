use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashSet;
use smallvec::SmallVec;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub fn run(fragments: &str, cells: &str) -> std::io::Result<()> {
    // Load the cell barcodes into a FxHashSet for fast lookups
    let cell_barcodes = load_cells(cells)?;

    // Filter the fragment file based on the cell barcodes
    filter_fragments(fragments, &cell_barcodes)?;

    Ok(())
}

fn load_cells<P: AsRef<Path>>(path: P) -> std::io::Result<FxHashSet<Box<str>>> {
    let reader = crate::reader::open_maybe_gzipped(path.as_ref())?;

    let mut cell_barcodes = FxHashSet::default();

    for line in reader.lines() {
        let line = line?;
        cell_barcodes.insert(line.into_boxed_str());
    }

    Ok(cell_barcodes)
}

fn filter_fragments<P: AsRef<Path>>(
    fragments_path: P,
    cell_barcodes: &FxHashSet<Box<str>>,
) -> std::io::Result<()> {
    let fragments_file = File::open(fragments_path)?;
    let mut fragments_reader =
        BufReader::with_capacity(4 * 1024 * 1024, MultiGzDecoder::new(fragments_file));

    let stdout = std::io::stdout();
    let mut output_writer = BufWriter::with_capacity(4 * 1024 * 1024, stdout.lock());

    // Pre-allocate a buffer for collecting matching lines
    const BUFFER_SIZE: usize = 4 * 1024 * 1024;
    let mut output_buffer = Vec::with_capacity(BUFFER_SIZE);
    let mut line_count: u64 = 0;
    let mut matching_count: u64 = 0;
    let mut buffer = String::with_capacity(1024);

    loop {
        buffer.clear();
        match fragments_reader.read_line(&mut buffer) {
            Ok(0) => break, // End of file
            Ok(_) => {
                if buffer.ends_with('\n') {
                    buffer.pop();
                }
                if buffer.ends_with('\r') {
                    buffer.pop();
                }

                if buffer.starts_with('#') {
                    continue;
                }

                let fields = buffer.split('\t').collect::<SmallVec<[&str; 10]>>();
                if fields.len() < 4 {
                    continue;
                }

                let barcode_str = fields[3];
                if cell_barcodes.contains(barcode_str) {
                    output_buffer.extend_from_slice(buffer.as_bytes());
                    output_buffer.push(b'\n');
                    matching_count += 1;

                    if output_buffer.len() >= BUFFER_SIZE {
                        output_writer.write_all(&output_buffer)?;
                        output_buffer.clear();
                    }
                }

                line_count += 1;
                if line_count % 1_000_000 == 0 {
                    eprint!(
                        "\rProcessed {} M lines, matched {} M fragments",
                        line_count / 1_000_000,
                        matching_count / 1_000_000
                    );
                    std::io::stderr().flush().expect("Can't flush stderr");
                }
            }
            Err(e) => return Err(e),
        }
    }

    // Write any remaining data
    if !output_buffer.is_empty() {
        output_writer.write_all(&output_buffer)?;
    }

    // Explicitly flush the writer
    output_writer.flush()?;

    eprintln!(
        "\nTotal: processed {} fragments, matched {} fragments",
        line_count, matching_count
    );

    Ok(())
}
