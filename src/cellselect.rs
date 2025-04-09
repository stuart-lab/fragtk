use std::io;
use std::error::Error;
use std::path::Path;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Write;
use std::thread;
use std::sync::mpsc;
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashMap;
use log::info;


pub fn cellselect(
    fragments: &str,
    outfile: &str,
    threshold: &Option<usize>,
    ncells: &Option<usize>
) -> Result<(), Box<dyn Error>> {

    let frag_file = Path::new(fragments)
        .canonicalize()
        .expect("Can't find path to input fragment file");
    info!("Received fragment file: {:?}", frag_file);

    info!("Output file: {:?}", outfile);

    // Get either threshold or ncells
    match (threshold, ncells) {
        (Some(t), None) => {
            info!("Cell count cutoff: {:?}", t);
        },
        (None, Some(n)) => {
            info!("Cell number cutoff: {:?}", n);
        },
        (None, None) => {
            eprintln!("Either --threshold or --ncells must be specified");
            std::process::exit(1);
        },
        (Some(_), Some(_)) => {
            eprintln!("Cannot specify both --threshold and --ncells");
            std::process::exit(1);
        }
    };

    let bc_count = count_barcodes(&frag_file)?;
    let selected = if threshold.is_some() {
        select_barcodes(&bc_count, &threshold.unwrap())?
    } else {
        top_barcodes(&bc_count, &ncells.unwrap())?
    };

    // Output results to the specified file
    let mut writer = File::create(outfile)?;
    let mut output = String::new();

    for (barcode, count) in &bc_count {
        output.push_str(&format!("{}\t{}\n", barcode, count));
    }
    writer.write_all(output.as_bytes())?;

    // print selected cells to stdout
    for cell in selected {
        println!("{}", cell);
    }

    Ok(())
}

fn select_barcodes(
    barcodes: &FxHashMap<String, usize>,
    count_cutoff: &usize,
) -> io::Result<Vec<String>> {

    // iterate over key, value entries, adding cells if count is greater than threshold
    let mut filtered_cells = Vec::new();
    for (cell_barcode, &count) in barcodes.iter() {
        if count > *count_cutoff {
            filtered_cells.push(cell_barcode.clone());
        }
    }

    Ok(filtered_cells)
}

fn top_barcodes(
    barcodes: &FxHashMap<String, usize>,
    ncells: &usize,
) -> io::Result<Vec<String>> {
    // create vectors of cells and counts
    let mut cells: Vec<String> = Vec::new();
    let mut counts: Vec<usize> = Vec::new();

    // iterate over barcode hashmap, filling in the vectors
    for (cell_barcode, &count) in barcodes.iter() {
        cells.push(cell_barcode.clone());
        counts.push(count);
    }

    // create index vector and sort it based on counts
    let mut idx: Vec<usize> = (0..cells.len()).collect();
    idx.sort_unstable_by(|&a, &b| counts[b].cmp(&counts[a]));

    if cells.len() < *ncells {
        eprintln!("Warning: Only {} cells available, fewer than requested {}", cells.len(), ncells);
    }

    // Take the top n cells using the sorted indices
    let n = std::cmp::min(*ncells, cells.len());
    let selected: Vec<String> = idx
        .into_iter()
        .take(n)
        .map(|i| cells[i].clone())
        .collect();

    Ok(selected)
}

fn count_barcodes(frag_file: &Path) -> io::Result<FxHashMap<String, usize>> {
    // hashmap for cell barcode counts
    let mut cells: FxHashMap<String, usize> = FxHashMap::default();

    // Create a channel for communication between the decompression and processing threads
    let (tx, rx) = mpsc::sync_channel(100);

    // Spawn the decompression thread
    let frag_file = frag_file.to_path_buf();
    let decompress_handle = thread::spawn(move || {
        // Open the fragment file
        let file = match File::open(&frag_file) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Failed to open fragment file: {}", e);
                return;
            }
        };
        
        // Create buffered reader
        let reader = BufReader::with_capacity(4 * 1024 * 1024, MultiGzDecoder::new(file));
        
        // Number of lines to read in each chunk
        const CHUNK_SIZE: usize = 10_000;
        let mut lines = Vec::with_capacity(CHUNK_SIZE);
        let mut total_lines = 0;
        
        // Read the file line by line
        for line_result in reader.lines() {
            match line_result {
                Ok(line) => {
                    // Skip header lines
                    if !line.starts_with('#') {
                        lines.push(line);
                        
                        // When we have a full chunk, send it to the processing thread
                        if lines.len() >= CHUNK_SIZE {
                            total_lines += lines.len();
                            
                            // Report progress
                            if total_lines % 1_000_000 == 0 {
                                eprint!("\rProcessed {} M fragments", total_lines / 1_000_000);
                                std::io::stderr().flush().unwrap();
                            }
                            
                            // Create a new vector and send the current one
                            let chunk_to_send = std::mem::replace(&mut lines, Vec::with_capacity(CHUNK_SIZE));
                            if tx.send(chunk_to_send).is_err() {
                                // Channel closed, receiver dropped
                                break;
                            }
                        }
                    }
                },
                Err(e) => {
                    eprintln!("Error reading line: {}", e);
                    break;
                }
            }
        }
        
        // Send any remaining lines
        if !lines.is_empty() {
            total_lines += lines.len();
            let _ = tx.send(lines);
        }
        
        eprintln!("\nFinished reading {} total fragments", total_lines);
    });

    // Process chunks from the channel
    for chunk in rx {
        for line in chunk {
            // Count barcodes from each line
            let fields = line.split('\t').collect::<smallvec::SmallVec<[&str; 10]>>();
            // let fields: Vec<&str> = line.split('\t').collect();
            
            // Update count for cell barcode
            if let Some(cell_barcode) = fields.get(3) {
                let cell_barcode = cell_barcode.to_string();
                *cells.entry(cell_barcode).or_insert(0) += 1;
            }
        }
    }
    eprintln!();

    // Join thread to ensure it completes
    decompress_handle.join().expect("Failed to join decompression thread");

    eprintln!("Found {} unique cell barcodes", cells.len());
    Ok(cells)
}