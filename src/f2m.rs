use std::{
    io,
    fs,
    path::Path,
    error::Error,
    fs::File,
    io::BufReader,
    io::BufRead,
    io::Write,
    sync::mpsc,
    thread,
};
use std::fmt::Write as FmtWrite; // for write! on String
use rust_lapper::{Interval, Lapper};
use flate2::read::MultiGzDecoder;
use flate2::Compression;
use log::error;
use log::info;
use rustc_hash::FxHashMap;
use gzp::{
    deflate::Gzip,
    ZWriter,
    par::compress::{ParCompress, ParCompressBuilder},
};
use smallvec::SmallVec;
use lexical_core::parse;
use itoa;

pub fn f2m(
    fragments: &str,
    bed: &str,
    cells: &str,
    outdir: &str,
    num_threads: usize,
    group: bool,
    pic: bool
) -> Result<(), Box<dyn Error>> {

    let frag_file = Path::new(fragments)
        .canonicalize()
        .expect("Can't find path to input fragment file");
    info!("Received fragment file: {:?}", frag_file);

    let bed_file = Path::new(bed)
        .canonicalize()
        .expect("Can't find path to input BED file");
    info!("Received BED file: {:?}", bed_file);

    let cell_file = Path::new(cells)
        .canonicalize()
        .expect("Can't find path to input cell file");
    info!("Received cell file: {:?}", cell_file);

    info!("Received output directory: {:?}", outdir);
    info!("Grouping peaks: {:?}", group);

    let output_path = Path::new(outdir);

    // Create the directory if it does not exist
    if !output_path.exists() {
        if let Err(e) = fs::create_dir_all(output_path) {
            eprintln!("Failed to create output directory: {}", e);
            std::process::exit(1);
        }
    }

    // make sure output is a directory
    match fs::metadata(output_path) {
        Ok(metadata) => {
            if metadata.is_dir() {
                info!("{:?} is a directory.", output_path);
            } else {
                eprintln!("Provided output is not a directory: {}", output_path.display());
                std::process::exit(1);
            }
        }
        Err(e) => {
            eprintln!("Failed to get metadata for {:?}: {}", output_path, e);
            std::process::exit(1);
        }
    }

    fcount(&frag_file, &bed_file, &cell_file, output_path, group, pic, num_threads)?;
    
    Ok(())
}

fn fcount(
    frag_file: &Path,
    bed_file: &Path,
    cell_file: &Path,
    output: &Path,
    group: bool,
    pic: bool,
    num_threads: usize,
) -> io::Result<()> {
    info!(
        "Processing fragment file: {:?}, BED file: {:?}, Cell file: {:?}",
        frag_file, bed_file, cell_file
    );

    // create BED intervals for overlaps with fragment coordinates
    // returns hashmap with each key being chromosome name
    // each value is intervals for that chromosome
    // interval value gives the index of the feature
    // also writes features to output directory to avoid second iteration of file
    // write features
    let feature_path = output.join("features.tsv.gz");
    info!("Writing output feature file: {:?}", &feature_path);
    let (total_peaks, peaks) = match peak_intervals(bed_file, group, &feature_path, num_threads) {
        Ok(trees) => trees,
        Err(e) => {
            error!("Failed to read BED file: {}", e);
            return Err(e);
        }
    };
    
    // create hashmap for cell barcodes
    let cell_file_handle = File::open(cell_file)?;
    let cellreader: Box<dyn BufRead> = if cell_file.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(cell_file_handle)))
    } else {
        Box::new(BufReader::new(cell_file_handle))
    };
    
    let mut cells: FxHashMap<String, u32> = FxHashMap::default();
    for (index, line) in cellreader.lines().enumerate() {
        let line = line?;
        let index_u32 = index as u32;
        cells.insert(line, index_u32);
    }
    
    // Estimate the peak-cell count map size for preallocation
    let cell_count = cells.len();
    info!("Loaded {} cell barcodes", cell_count);
    
    // Estimate an average of cells per peak for preallocation 
    let avg_cells_per_peak = (cell_count / 5).min(1000);

    // vector of features
    // each element is hashmap of cell: count
    let mut peak_cell_counts: Vec<FxHashMap<u32, u32>> = Vec::with_capacity(total_peaks);
    for _ in 0..total_peaks {
        peak_cell_counts.push(FxHashMap::with_capacity_and_hasher(
            avg_cells_per_peak, 
            Default::default()
        ));
    }

    // Create a channel for communication between threads
    let (tx, rx) = mpsc::sync_channel(100);
    
    // Spawn thread for decompression and reading
    let frag_file = frag_file.to_path_buf();
    let reader_handle = thread::spawn(move || {
        let file = match File::open(&frag_file) {
            Ok(f) => f,
            Err(e) => {
                error!("Failed to open fragment file: {}", e);
                return;
            }
        };
        
        let reader = BufReader::with_capacity(4 * 1024 * 1024, MultiGzDecoder::new(file));
        
        // Process fragments in chunks for better performance
        const CHUNK_SIZE: usize = 10_000;
        let mut fragments = Vec::with_capacity(CHUNK_SIZE);
        let mut total_fragments = 0;
        
        for line_result in reader.lines() {
            match line_result {
                Ok(line) => {
                    // Skip header lines
                    if !line.starts_with('#') {
                        fragments.push(line);
                        
                        if fragments.len() >= CHUNK_SIZE {
                            total_fragments += fragments.len();
                            
                            // Report progress
                            if total_fragments % 1_000_000 == 0 {
                                print!("\rProcessed {} M fragments", total_fragments / 1_000_000);
                                std::io::stdout().flush().expect("Can't flush output");
                            }
                            
                            // Send chunks for processing
                            let chunk_to_send = std::mem::replace(&mut fragments, Vec::with_capacity(CHUNK_SIZE));
                            if tx.send(chunk_to_send).is_err() {
                                break;
                            }
                        }
                    }
                },
                Err(e) => {
                    error!("Error reading fragment file: {}", e);
                    break;
                }
            }
        }
        
        // Send any remaining fragments
        if !fragments.is_empty() {
            total_fragments += fragments.len();
            let _ = tx.send(fragments);
        }
        
        eprintln!("\nFinished reading {} total fragments", total_fragments);
    });
    
    // Cache for chromosome lookups to improve performance
    let mut current_chrom = String::new();
    let mut current_lapper: Option<&Lapper<u32, usize>> = None;
    let mut cursor = 0;
    let mut check_end: bool;

    let mut startpos: u32;
    let mut endpos: u32;
    
    // Process chunks from the channel
    for chunk in rx {
        for line in chunk {            
            // Parse BED entry
            let fields = line.split('\t').collect::<SmallVec<[&str; 10]>>();
            if fields.len() < 4 {
                continue;
            }
            
            // Check if cell is to be included
            let cell_barcode: &str = fields[3];
            if let Some(&cell_index) = cells.get(cell_barcode) {
                check_end = true;
                
                // Create intervals from fragment entry
                let seqname: &str = fields[0];
                
                // Update lapper if chromosome changed
                if seqname != current_chrom {
                    current_chrom = seqname.to_string();
                    current_lapper = peaks.get(&current_chrom);
                    cursor = 0;
                }
                
                let start_str = fields[1];
                let end_str = fields[2];
                startpos = match parse::<u32>(start_str.as_bytes()) {
                    Ok(num) => num,
                    Err(_) => continue,
                };
                
                endpos = match parse::<u32>(end_str.as_bytes()) {
                    Ok(num) => num,
                    Err(_) => continue,
                };
                
                if let Some(lapper) = &current_lapper {
                    // seems to be a problem with seek if lapper has one element
                    // set cursor to 0
                    if lapper.intervals.len() == 1 {
                        cursor = 0;
                    }

                    // Check for overlaps at start position
                    for interval in lapper.seek(startpos, startpos + 1, &mut cursor) {
                        let peak_index = interval.val;
                        let peak_end = interval.stop;
                        *peak_cell_counts[peak_index].entry(cell_index).or_insert(0) += 1;
                        
                        if endpos < peak_end {
                            // Check if fragment end is behind peak end (it overlaps)
                            check_end = false;
                            
                            // From Paired Insertion Counting paper
                            // https://www.nature.com/articles/s41592-023-02103-7
                            //
                            // In PIC, for a given chromosome interval, if the pair of insertions of an ATAC-seq fragment
                            // are both within the interval, they are counted as one (pair); if only one insertion is within
                            // the interval and the other is outside the interval, also count one (pair).
                            if !pic {
                                *peak_cell_counts[peak_index].entry(cell_index).or_insert(0) += 1;
                            }
                        }
                    }
                    
                    // Check for overlaps at end position if needed
                    if check_end {
                        for interval in lapper.seek(endpos, endpos + 1, &mut cursor) {
                            let peak_index = interval.val;
                            *peak_cell_counts[peak_index].entry(cell_index).or_insert(0) += 1;
                        }
                    }
                }
            }
        }
    }
    
    // Wait for reader thread to complete
    reader_handle.join().expect("Reader thread panicked");
    
    for counts in &mut peak_cell_counts {
        counts.shrink_to_fit();
    }

    // write count matrix
    let counts_path = output.join("matrix.mtx.gz");
    info!("Writing output counts file: {:?}", &counts_path);
    write_matrix_market(&counts_path, &peak_cell_counts, total_peaks, cells.len(), num_threads)
        .expect("Failed to write matrix"); // features stored as rows

    // write cells
    let cell_path = output.join("barcodes.tsv.gz");
    info!("Writing output cells file: {:?}", &cell_path);
    write_cells(&cell_path, cell_file, num_threads)
        .expect("Failed to write cells");

    Ok(())
}

fn write_cells(
    outfile: &Path,
    cells: &Path,
    num_threads: usize,
) -> io::Result<()> {

    // If input is gzipped, just copy the file
    if cells.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        fs::copy(cells, outfile)?;
        info!("Copied gzipped cell barcodes to {:?}", outfile);
        return Ok(());
    }

    // Otherwise proceed with reading, compressing, and writing
    let input = File::open(cells)?;
    let reader = BufReader::new(input);

    let output = File::create(outfile)?;
    let mut writer: ParCompress<Gzip> = ParCompressBuilder::new()
        .compression_level(Compression::default())
        .num_threads(num_threads)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?
        .from_writer(output);

    for line in reader.lines() {
        writeln!(writer, "{}", line?)?;
    }

    writer.finish().map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    
    info!("Successfully wrote compressed cell barcodes to {:?}", outfile);
    Ok(())
}

fn write_matrix_market(
    outfile: &Path,
    peak_cell_counts: &[FxHashMap<u32, u32>],
    nrow: usize,
    ncol: usize,
    num_threads: usize,
) -> io::Result<()> {

    // get nonzero value count
    let nonzero: usize = peak_cell_counts.iter().map(|map| map.len()).sum();

    // create output file
    let writer = File::create(outfile)?;
    let mut encoder: ParCompress<Gzip> = ParCompressBuilder::new()
        .compression_level(Compression::default())  // Set compression level
        .num_threads(num_threads)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))? 
        .from_writer(writer);

    // Create a string buffer to collect all lines
    let mut output = String::with_capacity(2 * 1024 * 1024);

    // Write the header for the Matrix Market format
    output.push_str("%%MatrixMarket matrix coordinate integer general\n");
    output.push_str(&format!("%%metadata json: {{\"software_version\": \"fragtk-{}\"}}\n", env!("CARGO_PKG_VERSION")));
    output.push_str(&itoa::Buffer::new().format(nrow));
    output.push(' ');
    output.push_str(&itoa::Buffer::new().format(ncol));
    output.push(' ');
    output.push_str(&itoa::Buffer::new().format(nonzero));
    output.push('\n');
    encoder.write_all(output.as_bytes())?;
    output.clear();

    const CHUNK_SIZE: usize = 50_000;
    let mut entries_in_chunk = 0;

    let mut row_buf = itoa::Buffer::new();
    let mut col_buf = itoa::Buffer::new();
    let mut val_buf = itoa::Buffer::new();

    for (index, hashmap) in peak_cell_counts.iter().enumerate() {
        for (key, value) in hashmap.iter() {

            write!(
                &mut output,
                "{} {} {}\n",
                row_buf.format(index + 1),
                col_buf.format(key + 1),
                val_buf.format(*value)
            ).unwrap();
    
            entries_in_chunk += 1;
    
            if entries_in_chunk >= CHUNK_SIZE {
                encoder.write_all(output.as_bytes())?;
                output.clear();
                entries_in_chunk = 0;
            }
        }
    }

    // Write the remaining string buffer
    if !output.is_empty() {
        encoder.write_all(output.as_bytes())?;
    }

    encoder.finish().map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

    Ok(())
}

fn peak_intervals(
    bed_file: &Path,
    group: bool,
    outfile: &Path,
    num_threads: usize,
) -> io::Result<(usize, FxHashMap<String, Lapper<u32, usize>>)> {

    // feature file
    let writer = File::create(outfile)?;
    let mut writer: ParCompress<Gzip> = ParCompressBuilder::new()
        .compression_level(Compression::default())
        .num_threads(num_threads)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?
        .from_writer(writer);
    
    // bed file reader
    let file = File::open(bed_file)?;
    let reader: Box<dyn BufRead> = if bed_file.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    
    // hashmap of peak intervals for each chromosome
    let mut chromosome_trees: FxHashMap<String, Vec<Interval<u32, usize>>> = FxHashMap::default();
    
    // Store peak group name and corresponding index
    let mut peak_group_index: FxHashMap<String, usize> = FxHashMap::default();
    
    // track total number of peaks
    let mut total_peaks: usize = 0;

    // index for peak groups
    let mut current_index: usize = 0;

    // skipped lines
    let mut skipped_lines: usize = 0;

    for (index, line) in reader.lines().enumerate() {

        match line {
            Ok(line) => {
                if line.starts_with('#') {
                    skipped_lines += 1;
                    continue;
                }
                let fields = line.split('\t').collect::<SmallVec<[&str; 10]>>();
                if fields.len() >= 3 {
                    let chromosome = fields[0].to_string();
                    let start: u32 = match parse::<u32>(fields[1].trim().as_bytes()) {
                        Ok(num) => num,
                        Err(_) => {
                            return Err(io::Error::new(io::ErrorKind::InvalidData, 
                                format!("Line {}: Failed to parse start position", index + 1)));
                        }
                    };
                    let end: u32 = match parse::<u32>(fields[2].trim().as_bytes()) {
                        Ok(num) => num,
                        Err(_) => {
                            return Err(io::Error::new(io::ErrorKind::InvalidData,
                                format!("Line {}: Failed to parse end position", index +1)));
                        }
                    };

                    let intervals = chromosome_trees.entry(chromosome.clone()).or_insert_with(Vec::new);

                    if group && (fields.len() >= 4) {
                        let peakgroup: String = match fields[3].parse() {
                            Ok(num) => num,
                            Err(_) => {
                                return Err(io::Error::new(io::ErrorKind::InvalidData,
                                    format!("Line {}: Failed to parse group information", index + 1)));
                            }
                        };

                        let group_index = peak_group_index.entry(peakgroup.clone()).or_insert_with(|| {
                            writeln!(writer, "{}", peakgroup).expect("Failed to write");
                            let idx: usize = current_index - skipped_lines;
                            current_index += 1;
                            idx
                        });

                        intervals.push(Interval { start, stop: end, val: *group_index });
                    } else {
                        intervals.push(Interval { start, stop: end, val: index - skipped_lines});
                        writeln!(writer, "{}-{}-{}", chromosome, start, end)?;
                    }
                    total_peaks += 1;
                } else {
                    return Err(io::Error::new(io::ErrorKind::InvalidData,
                        format!("Line {}: Less than three fields", index + 1)));
                }
            },
            Err(e) => {
                return Err(io::Error::new(io::ErrorKind::InvalidData,
                    format!("Error reading line {}: {}", index + 1, e)));
            }
        }
    }

    let lapper_map = chromosome_trees.into_iter()
        .map(|(chr, intervals)| (chr, Lapper::new(intervals)))
        .collect();

    if group {
        total_peaks = current_index;
    }

    // Finalize the compression, converting GzpError to io::Error
    writer.finish().map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

    Ok((total_peaks, lapper_map))
}