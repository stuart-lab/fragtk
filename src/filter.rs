use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashSet;

pub fn run(
    fragments: &str,
    cells: &str
) -> std::io::Result<()> {
    // Load the cell barcodes into a FxHashSet for fast lookups
    let cell_barcodes = load_cells(cells)?;

    // Filter the fragment file based on the cell barcodes
    filter_fragments(fragments, &cell_barcodes)?;

    Ok(())
}

fn load_cells<P: AsRef<Path>>(path: P) -> std::io::Result<FxHashSet<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut cell_barcodes = FxHashSet::default();

    for line in reader.lines() {
        let line = line?;
        cell_barcodes.insert(line);
    }

    Ok(cell_barcodes)
}

fn filter_fragments<P: AsRef<Path>>(
    fragments_path: P,
    cell_barcodes: &FxHashSet<String>,
) -> std::io::Result<()> {
    let fragments_file = File::open(fragments_path)?;
    let mut fragments_reader = BufReader::with_capacity(1024 * 1024, MultiGzDecoder::new(fragments_file));

    let stdout = std::io::stdout();
    let mut output_writer = stdout.lock();

    let mut line_count: u64 = 0;
    let mut buffer = String::with_capacity(1024);

    loop {
        buffer.clear();
        match fragments_reader.read_line(&mut buffer) {
            Ok(0) => break, // End of file
            Ok(_) => {
                // Remove trailing newline
                if buffer.ends_with('\n') {
                    buffer.pop();
                }
                if buffer.ends_with('\r') {
                    buffer.pop();
                }

                // Skip comment lines
                if buffer.starts_with('#') {
                    continue;
                }

                if let Some(barcode) = buffer.split('\t').nth(3) {
                    if cell_barcodes.contains(barcode) {
                        writeln!(output_writer, "{}", buffer)?;
                    }
                }

                line_count += 1;
                if line_count % 1_000_000 == 0 {
                    eprint!("\rProcessed {} M lines", line_count / 1_000_000);
                    std::io::stderr().flush().expect("Can't flush stderr");
                }
            }
            Err(e) => return Err(e),
        }
    }

    Ok(())
}