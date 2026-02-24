use std::{
    io,
    fs,
    path::Path,
    error::Error,
    fs::File,
    io::BufReader,
    io::BufWriter,
    io::BufRead,
    io::Write,
    sync::mpsc,
    thread,
    fs::OpenOptions,
};
use std::fmt::Write as FmtWrite;
use rust_lapper::{Interval, Lapper};
use flate2::Compression;
use flate2::write::GzEncoder;
use log::error;
use log::info;
use rustc_hash::{FxHashMap, FxHashSet};
use gzp::{
    deflate::Gzip,
    ZWriter,
    par::compress::{ParCompress, ParCompressBuilder},
};
use smallvec::SmallVec;
use lexical_core::parse;
use itoa;
use tempfile::NamedTempFile;

const MAX_COUNT: u32 = u16::MAX as u32;

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
        fs::create_dir_all(output_path)
            .map_err(|e| format!("Failed to create output directory: {}", e))?;
    }

    // make sure output is a directory
    let metadata = fs::metadata(output_path)
        .map_err(|e| format!("Failed to get metadata for {:?}: {}", output_path, e))?;
    if !metadata.is_dir() {
        return Err(format!("Provided output is not a directory: {}", output_path.display()).into());
    }
    info!("{:?} is a directory.", output_path);

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
    let cellreader = crate::reader::open_maybe_gzipped(cell_file)?;
    
    let mut cells: FxHashMap<Box<str>, u32> = FxHashMap::default();
    for (index, line) in cellreader.lines().enumerate() {
        let line = line?;
        let index_u32 = index as u32;
        cells.insert(line.into_boxed_str(), index_u32);
    }
    
    let cell_count = cells.len();
    info!("Loaded {} cell barcodes", cell_count);

    // Create channel for writer thread communication
    let (counts_tx, counts_rx) = mpsc::sync_channel::<(FxHashMap<(u32, u32), u16>, bool)>(3);

    let temp_file = NamedTempFile::new()?;
    let temp_path_clone = temp_file.path().to_str().unwrap().to_string();
    info!("temp_file: {:?}", temp_path_clone);

    let writer_threads = num_threads;

    // Spawn writer thread first
    let writer_handle = thread::spawn(move || -> io::Result<()> {
        let mut result = Ok(());
        while let Ok((counts, is_last)) = counts_rx.recv() {
            match write_matrix_market(&temp_path_clone, counts, writer_threads) {
                Ok(_) => {
                    if is_last {
                        break;
                    }
                }
                Err(e) => {
                    error!("Error writing matrix market: {}", e);
                    result = Err(e);
                    break;
                }
            }
        }
        result
    });
    
    // Spawn reader thread for decompression
    let (reader_handle, rx, pool_tx) = crate::reader::spawn_fragment_reader(&frag_file);

    let mut nonzero_counts: u64 = 0;
    let mut past_chromosomes: FxHashSet<Box<str>> = FxHashSet::default();
    let mut peak_cell_counts: FxHashMap<(u32, u32), u16> = FxHashMap::default();
    let mut current_chrom = String::new();
    let mut current_lapper: Option<&Lapper<u32, usize>> = None;
    let mut cursor = 0;
    let mut check_end: bool;
    
    for mut chunk in rx {
        for line in chunk.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            
            let mut iter = line.splitn(5, |&b| b == b'\t');
            let seqname_bytes = match iter.next() { Some(b) => b, None => continue };
            let start_bytes = match iter.next() { Some(b) => b, None => continue };
            let end_bytes = match iter.next() { Some(b) => b, None => continue };
            let mut barcode_bytes = match iter.next() { Some(b) => b, None => continue };
            
            // Trim trailing \r
            if barcode_bytes.ends_with(b"\r") {
                barcode_bytes = &barcode_bytes[..barcode_bytes.len() - 1];
            }

            let cell_barcode = unsafe { std::str::from_utf8_unchecked(barcode_bytes) };

            // Check if cell is to be included
            if let Some(&cell_index) = cells.get(cell_barcode) {
                check_end = true;
                
                // Create intervals from fragment entry
                let seqname = unsafe { std::str::from_utf8_unchecked(seqname_bytes) };
                
                // check if chromosome changed
                // if so: update lapper, write previous chromosome's counts to file, reset peak_cell_counts
                if seqname != current_chrom {
                    if past_chromosomes.contains(seqname) {
                        // return error, entries not sorted
                        let line_str = unsafe { std::str::from_utf8_unchecked(line) };
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Fragment file is not sorted by chromosome: {}", line_str)
                        ));
                    }

                    if !current_chrom.is_empty() {
                        info!("Writing counts for chromosome: {}", current_chrom);
                        past_chromosomes.insert(current_chrom.clone().into_boxed_str()); // remember what has been processed
                        nonzero_counts += peak_cell_counts.len() as u64;

                        // Send current counts to writer thread
                        if !peak_cell_counts.is_empty() {
                            let cap = peak_cell_counts.capacity();
                            // send counts to writer thread, replace with empty hashmap but keep capacity
                            let counts_to_send = std::mem::replace(&mut peak_cell_counts, FxHashMap::with_capacity_and_hasher(cap, Default::default()));
                            if let Err(e) = counts_tx.send((counts_to_send, false)) {
                                error!("Failed to send chromosome counts: {}", e);
                                return Err(io::Error::new(io::ErrorKind::Other, e));
                            }
                        }
                    }

                    current_chrom = seqname.to_string();
                    current_lapper = peaks.get(&current_chrom);
                    cursor = 0;
                }
                
                let startpos = match parse::<u32>(start_bytes) {
                    Ok(num) => num,
                    Err(_) => continue,
                };
                
                let endpos = match parse::<u32>(end_bytes) {
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
                        let peak_index = interval.val as u32;
                        let peak_end = interval.stop;
                        let count = peak_cell_counts.entry((peak_index, cell_index)).or_insert(0);
                        if *count < MAX_COUNT as u16 {
                            *count += 1;
                        }
                        
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
                                let count = peak_cell_counts.entry((peak_index, cell_index)).or_insert(0);
                                if *count < MAX_COUNT as u16 {
                                    *count += 1;
                                }
                            }
                        }
                    }
                    
                    // Check for overlaps at end position if needed
                    if check_end {
                        for interval in lapper.seek(endpos, endpos + 1, &mut cursor) {
                            let peak_index = interval.val as u32;
                            let count = peak_cell_counts.entry((peak_index, cell_index)).or_insert(0);
                            if *count < MAX_COUNT as u16 {
                                *count += 1;
                            }
                        }
                    }
                }
            }
        }
        chunk.clear();
        let _ = pool_tx.send(chunk);
    }

    // Wait for reader thread to complete
    reader_handle.join().expect("Reader thread panicked");

    // Send final counts and signal completion
    nonzero_counts += peak_cell_counts.len() as u64;
    let final_counts = std::mem::replace(&mut peak_cell_counts, FxHashMap::default());
    if let Err(e) = counts_tx.send((final_counts, true)) {
        error!("Failed to send final counts: {}", e);
        return Err(io::Error::new(io::ErrorKind::Other, e));
    }

    // Wait for writer thread to complete
    let _ = writer_handle.join().expect("Writer thread panicked");

    // write mtx header with proper gzip compression
    info!("Writing output counts file: {:?}", &output.join("matrix.mtx.gz"));
    let output_file = File::create(output.join("matrix.mtx.gz"))?;
    let mut header = Vec::new();
    {
        let mut header_writer = BufWriter::new(GzEncoder::new(&mut header, Compression::default()));
        writeln!(header_writer, "%%MatrixMarket matrix coordinate integer general")?;
        writeln!(header_writer, "%metadata json: {{\"software_version\": \"fragtk-{}\", \"command\": \"fragtk matrix\"}}", env!("CARGO_PKG_VERSION"))?;
        writeln!(header_writer, "{} {} {}", total_peaks, cells.len(), nonzero_counts)?;
        header_writer.flush()?;
    }
    
    // Write gzipped header and concatenate with gzipped counts
    info!("Copying counts to output file");
    let mut output_writer = BufWriter::new(output_file);
    output_writer.write_all(&header)?;
    let mut temp_reader = BufReader::new(File::open(temp_file.path().to_str().unwrap())?);
    io::copy(&mut temp_reader, &mut output_writer)?;
    output_writer.flush()?;

    // NamedTempFile handles cleanup on drop
    drop(temp_file);
    
    // write cells
    info!("Writing output cells file: {:?}", &output.join("barcodes.tsv.gz"));
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
    if crate::reader::is_gzipped(cells)? {
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
    file_path: &str,
    peak_cell_counts: FxHashMap<(u32, u32), u16>,
    num_threads: usize,
) -> io::Result<()> {

    // write count information only
    // header not written

    let file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(file_path)?;

    let mut encoder: ParCompress<Gzip> = ParCompressBuilder::new()
        .compression_level(Compression::default())
        .num_threads(num_threads)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?
        .from_writer(file);

    // Create a string buffer to collect all lines
    let mut output = String::with_capacity(2 * 1024 * 1024);
    const CHUNK_SIZE: usize = 50_000;
    let mut entries_in_chunk = 0;

    let mut row_buf = itoa::Buffer::new();
    let mut col_buf = itoa::Buffer::new();
    let mut val_buf = itoa::Buffer::new();

    let mut sorted_counts: Vec<_> = peak_cell_counts.into_iter().collect();
    // Sort by row index (peak position), then by column index (cell)
    sorted_counts.sort_unstable_by_key(|&(k, _)| k);

    for (key, value) in sorted_counts {

            write!(
                &mut output,
                "{} {} {}\n",
                row_buf.format(key.0 + 1),
                col_buf.format(key.1 + 1),
                val_buf.format(value)
            ).unwrap();
    
            entries_in_chunk += 1;
    
            if entries_in_chunk >= CHUNK_SIZE {
                encoder.write_all(output.as_bytes())?;
                output.clear();
                entries_in_chunk = 0;
        }
    }

    // Write the remaining string buffer
    if !output.is_empty() {
        encoder.write_all(output.as_bytes())?;
    }

    encoder.flush()?;
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
    let reader = crate::reader::open_maybe_gzipped(bed_file)?;
    
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
                            let idx = current_index;
                            current_index += 1;
                            idx
                        });

                        intervals.push(Interval { start, stop: end, val: *group_index });
                    } else {
                        intervals.push(Interval { start, stop: end, val: index - skipped_lines});
                        writeln!(writer, "{}:{}-{}", chromosome, start, end)?;
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