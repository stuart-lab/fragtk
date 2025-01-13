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

pub fn cellselect(matches: &clap::ArgMatches) -> Result<(), Box<dyn Error>> {

    let frag_file = Path::new(matches.get_one::<String>("fragments").unwrap())
        .canonicalize()
        .expect("Can't find path to input fragment file");
    info!("Received fragment file: {:?}", frag_file);

    let output_file = matches.get_one::<String>("outfile").unwrap();
    info!("Output file: {:?}", output_file);

    // Get either threshold or ncells from matches
    let (threshold, ncells) = match (matches.get_one::<String>("threshold"), 
                                   matches.get_one::<String>("ncells")) {
        (Some(t), None) => {
            let threshold = t.parse().unwrap_or_else(|_| {
                eprintln!("Failed to parse threshold as usize");
                std::process::exit(1);
            });
            info!("Cell count cutoff: {:?}", threshold);
            (Some(threshold), None)
        },
        (None, Some(n)) => {
            let ncells = n.parse().unwrap_or_else(|_| {
                eprintln!("Failed to parse ncells as usize");
                std::process::exit(1);
            });
            info!("Cell number cutoff: {:?}", ncells);
            (None, Some(ncells))
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
    let mut writer = File::create(output_file)?;
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

fn count_barcodes(frag_file: &Path,) -> io::Result<FxHashMap<String, usize>> {

    // hashmap for cell barcode counts
    let mut cells: FxHashMap<String, usize> = FxHashMap::default();

    // Create a channel for communication between the decompression and processing threads
    let (tx, rx) = mpsc::sync_channel(500);

    // Spawn the decompression thread
    let frag_file = frag_file.to_path_buf();
    let decompress_handle = thread::spawn(move || {
        let reader = BufReader::new(MultiGzDecoder::new(File::open(frag_file).expect("Failed to open fragment file")));
        for line in reader.lines() {
            let line = line.expect("Failed to read line");
            if tx.send(line).is_err() {
                break;
            }
        }
    });

    // Progress counter
    let mut line_count: u64 = 0;
    let update_interval = 1_000_000;

    for line in rx {

        // Skip header lines that start with #
        if line.starts_with('#') {
            continue;
        }

        line_count += 1;
        if line_count % update_interval == 0 {
            eprint!("\rProcessed {} M fragments", line_count / 1_000_000 );
            std::io::stdout().flush().expect("Can't flush output");
        }

        // parse bed entry
        let fields: Vec<&str> = line.split('\t').collect();

        // update count for cell barcode
        if let Some(cell_barcode) = fields.get(3) {
            let cell_barcode = cell_barcode.to_string();
            *cells.entry(cell_barcode).or_insert(0) += 1;
        }
    }
    eprintln!();

    // Join thread to ensure it completes
    decompress_handle.join().expect("Failed to join decompression thread");

    Ok(cells)
}