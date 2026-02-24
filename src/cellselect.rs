use std::io;
use std::error::Error;
use std::path::Path;
use std::fs::File;
use std::io::Write;
use rustc_hash::FxHashMap;
use log::info;
use std::fs;

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
            return Err("Either --threshold or --ncells must be specified".into());
        },
        (Some(_), Some(_)) => {
            return Err("Cannot specify both --threshold and --ncells".into());
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
    barcodes: &FxHashMap<Box<str>, usize>,
    count_cutoff: &usize,
) -> io::Result<Vec<Box<str>>> {

    // iterate over key, value entries, adding cells if count is greater than threshold
    let mut filtered_cells = Vec::new();
    for (cell_barcode, &count) in barcodes.iter() {
        if count >= *count_cutoff {
            filtered_cells.push(cell_barcode.clone());
        }
    }

    Ok(filtered_cells)
}

fn top_barcodes(
    barcodes: &FxHashMap<Box<str>, usize>,
    ncells: &usize,
) -> io::Result<Vec<Box<str>>> {
    // create vectors of cells and counts
    let mut cells: Vec<Box<str>> = Vec::new();
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
    let selected: Vec<Box<str>> = idx
        .into_iter()
        .take(n)
        .map(|i| cells[i].clone())
        .collect();

    Ok(selected)
}

fn count_barcodes(frag_file: &Path) -> io::Result<FxHashMap<Box<str>, usize>> {

    let metadata = fs::metadata(&frag_file)?;
    let file_size = metadata.len() as usize;
    let estimated_lines: usize = file_size / 100;
    let estimated_cell_count: usize = (estimated_lines / 10_000).max(1000);

    let mut cells: FxHashMap<Box<str>, usize> = FxHashMap::with_capacity_and_hasher(
        estimated_cell_count, 
        Default::default()
    );

    // Spawn reader thread for decompression
    let (decompress_handle, rx, pool_tx) = crate::reader::spawn_fragment_reader(frag_file);

    // Process chunks from the channel
    for mut chunk in rx {
        for line in chunk.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }

            let mut iter = line.splitn(5, |&b| b == b'\t');
            let _ = match iter.next() { Some(b) => b, None => continue }; // chr
            let _ = match iter.next() { Some(b) => b, None => continue }; // start
            let _ = match iter.next() { Some(b) => b, None => continue }; // end
            let mut barcode_bytes = match iter.next() { Some(b) => b, None => continue };

            // Trim trailing \r
            if barcode_bytes.ends_with(b"\r") {
                barcode_bytes = &barcode_bytes[..barcode_bytes.len() - 1];
            }

            let barcode_str = unsafe { std::str::from_utf8_unchecked(barcode_bytes) };
            if let Some(count) = cells.get_mut(barcode_str) {
                *count += 1;
            } else {
                cells.insert(barcode_str.into(), 1);
            }
        }
        chunk.clear();
        let _ = pool_tx.send(chunk);
    }

    // Join thread to ensure it completes
    decompress_handle.join().expect("Failed to join decompression thread");
    eprintln!("Found {} unique cell barcodes", cells.len());
    Ok(cells)
}