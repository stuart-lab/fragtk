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
use crate::intervals::{seek_position, SeekCursor};
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
use tempfile::NamedTempFile;

#[cfg(feature = "hdf5")]
use hdf5::File as H5File;
#[cfg(feature = "hdf5")]
use hdf5::types::VarLenUnicode;
#[cfg(feature = "hdf5")]
use ndarray::Array1;

const MAX_COUNT: u32 = u16::MAX as u32;

/// Whether an `--h5` output path names the file to write rather than a directory to
/// write `matrix.h5` into.
///
/// HDF5 output is a single file, so `-o counts.h5` should produce that file rather
/// than a directory called `counts.h5` containing `matrix.h5`.
#[cfg(feature = "hdf5")]
fn h5_output_is_file(path: &Path) -> bool {
    if path.is_dir() {
        return false;
    }
    if path.is_file() {
        return true;
    }
    matches!(
        path.extension().and_then(|e| e.to_str()),
        Some(ext) if ext.eq_ignore_ascii_case("h5") || ext.eq_ignore_ascii_case("hdf5")
    )
}

pub fn f2m(
    fragments: &str,
    bed: &str,
    cells: &str,
    outdir: &str,
    num_threads: usize,
    group: bool,
    pic: bool,
    #[cfg(feature = "hdf5")]
    h5: bool
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

    // With --h5 the output is a single file, and may name that file rather than a
    // directory to put it in. In that case only its parent has to exist.
    #[cfg(feature = "hdf5")]
    let output_is_h5_file = h5 && h5_output_is_file(output_path);
    #[cfg(not(feature = "hdf5"))]
    let output_is_h5_file = false;

    if output_is_h5_file {
        info!("Writing HDF5 output to file: {:?}", output_path);
        if let Some(parent) = output_path.parent() {
            if !parent.as_os_str().is_empty() && !parent.exists() {
                fs::create_dir_all(parent)
                    .map_err(|e| format!("Failed to create output directory: {}", e))?;
            }
        }
    } else {
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
    }

    #[cfg(feature = "hdf5")]
    fcount(&frag_file, &bed_file, &cell_file, output_path, group, pic, num_threads, h5)?;

    #[cfg(not(feature = "hdf5"))]
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
    #[cfg(feature = "hdf5")]
    h5: bool
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
    
    // Provide a Some(...) path only if not generating h5, though we could write it anyway
    #[cfg(feature = "hdf5")]
    let outfile_arg = if h5 { None } else { Some(feature_path.as_path()) };
    #[cfg(not(feature = "hdf5"))]
    let outfile_arg = Some(feature_path.as_path());
    
    let (total_peaks, peaks, _features_list) = match peak_intervals(bed_file, group, outfile_arg, num_threads) {
        Ok(trees) => trees,
        Err(e) => {
            error!("Failed to read BED file: {}", e);
            return Err(e);
        }
    };
    
    // create hashmap for cell barcodes
    let cellreader = crate::reader::open_maybe_gzipped(cell_file)?;
    
    let mut cells: FxHashMap<Box<str>, u32> = FxHashMap::default();
    
    #[cfg(feature = "hdf5")]
    let mut barcodes_list = Vec::new();
    
    for (index, line) in cellreader.lines().enumerate() {
        let line = line?;
        let index_u32 = index as u32;
        // Barcodes index the matrix columns, but barcodes.tsv.gz is a verbatim copy of
        // this file. A duplicate would collapse two columns into one in the map while
        // leaving both lines in the output, so the mtx header would declare fewer
        // columns than the largest column index written.
        if let Some(first) = cells.insert(line.clone().into_boxed_str(), index_u32) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Duplicate cell barcode {:?} in cell file (lines {} and {}). \
                     Cell barcodes must be unique.",
                    line,
                    first + 1,
                    index_u32 + 1
                ),
            ));
        }
        #[cfg(feature = "hdf5")]
        if h5 { barcodes_list.push(line); }
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
    let writer_handle = thread::spawn(move || -> io::Result<Option<Vec<(u32, u32, u16)>>> {
        let mut result: io::Result<Option<Vec<(u32, u32, u16)>>> = Ok(None);
        
        #[cfg(feature = "hdf5")]
        let mut all_h5_counts = if h5 { Some(Vec::new()) } else { None };
        
        while let Ok((counts, is_last)) = counts_rx.recv() {
            #[cfg(feature = "hdf5")]
            let is_h5 = all_h5_counts.is_some();
            #[cfg(not(feature = "hdf5"))]
            let is_h5 = false;
            
            if is_h5 {
                #[cfg(feature = "hdf5")]
                if let Some(ref mut all) = all_h5_counts {
                    for (k, v) in counts {
                        all.push((k.1, k.0, v)); // cell, peak, count
                    }
                }
            } else {
                match write_matrix_market(&temp_path_clone, counts, writer_threads) {
                    Ok(_) => {}
                    Err(e) => {
                        error!("Error writing matrix market: {}", e);
                        result = Err(e);
                        break;
                    }
                }
            }
            if is_last {
                break;
            }
        }
        
        #[cfg(feature = "hdf5")]
        if let Ok(_) = result {
            if h5 {
                return Ok(all_h5_counts);
            }
        }
        
        match result {
            Ok(_) => Ok(None),
            Err(e) => Err(e),
        }
    });
    
    // Spawn reader thread for decompression
    let (reader_handle, rx, pool_tx) = crate::reader::spawn_fragment_reader(&frag_file);

    let mut nonzero_counts: u64 = 0;
    let mut past_chromosomes: FxHashSet<Box<str>> = FxHashSet::default();
    let mut peak_cell_counts: FxHashMap<(u32, u32), u16> = FxHashMap::default();
    let mut current_chrom = String::new();
    let mut last_start: u32 = 0;
    let mut current_lapper: Option<&Lapper<u32, usize>> = None;
    // fragment starts and fragment ends are separate forward-moving series
    let mut start_cursor = SeekCursor::new();
    let mut end_cursor = SeekCursor::new();
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
                    start_cursor.reset();
                    end_cursor.reset();
                    last_start = 0;
                }
                
                let startpos = match parse::<u32>(start_bytes) {
                    Ok(num) => num,
                    Err(_) => continue,
                };

                // Overlaps are found with a cursor that only moves forward, so a
                // fragment that starts before its predecessor can silently miss peaks
                if startpos < last_start {
                    let line_str = unsafe { std::str::from_utf8_unchecked(line) };
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "Fragment file is not sorted by position ({} < {} on {}): {}",
                            startpos, last_start, current_chrom, line_str
                        )
                    ));
                }
                last_start = startpos;
                
                let endpos = match parse::<u32>(end_bytes) {
                    Ok(num) => num,
                    Err(_) => continue,
                };
                
                if let Some(lapper) = &current_lapper {
                    // Check for overlaps at start position
                    for interval in seek_position(lapper, startpos, &mut start_cursor) {
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
                        for interval in seek_position(lapper, endpos, &mut end_cursor) {
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
    let _writer_result = writer_handle.join().expect("Writer thread panicked")?;

    #[cfg(feature = "hdf5")]
    if h5 {
        // Output the HDF5 file, either at the path given or as matrix.h5 inside it
        let h5_path = if h5_output_is_file(output) {
            output.to_path_buf()
        } else {
            output.join("matrix.h5")
        };
        info!("Writing output HDF5 file: {:?}", &h5_path);
        if let Some(all_counts) = _writer_result {
            if let Some(features) = _features_list {
                write_hdf5(&h5_path, all_counts, total_peaks, cells.len(), features, barcodes_list)?;
            } else {
                return Err(io::Error::new(io::ErrorKind::Other, "Features list missing for HDF5 output"));
            }
        } else {
            return Err(io::Error::new(io::ErrorKind::Other, "Counts list missing for HDF5 output"));
        }
    } else {
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
    
        // write cells
        info!("Writing output cells file: {:?}", &output.join("barcodes.tsv.gz"));
        let cell_path = output.join("barcodes.tsv.gz");
        info!("Writing output cells file: {:?}", &cell_path);
        crate::f2m::write_cells(&cell_path, cell_file, num_threads)
            .expect("Failed to write cells");
    }

    #[cfg(not(feature = "hdf5"))]
    {
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
    
        // write cells
        info!("Writing output cells file: {:?}", &output.join("barcodes.tsv.gz"));
        let cell_path = output.join("barcodes.tsv.gz");
        info!("Writing output cells file: {:?}", &cell_path);
        crate::f2m::write_cells(&cell_path, cell_file, num_threads)
            .expect("Failed to write cells");
    }

    // NamedTempFile handles cleanup on drop
    drop(temp_file);

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
    outfile: Option<&Path>,
    num_threads: usize,
) -> io::Result<(usize, FxHashMap<String, Lapper<u32, usize>>, Option<Vec<String>>)> {

    // feature file
    let mut writer = if let Some(out) = outfile {
        let w = File::create(out)?;
        let pw: ParCompress<Gzip> = ParCompressBuilder::new()
            .compression_level(Compression::default())
            .num_threads(num_threads)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?
            .from_writer(w);
        Some(pw)
    } else {
        None
    };
    
    // bed file reader
    let reader = crate::reader::open_maybe_gzipped(bed_file)?;
    
    // hashmap of peak intervals for each chromosome
    let mut chromosome_trees: FxHashMap<String, Vec<Interval<u32, usize>>> = FxHashMap::default();
    
    // Store peak group name and corresponding index
    let mut peak_group_index: FxHashMap<String, usize> = FxHashMap::default();
    
    let mut features: Option<Vec<String>> = if outfile.is_none() { Some(Vec::new()) } else { None };
    
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
                            if let Some(ref mut w) = writer {
                                writeln!(w, "{}", peakgroup).expect("Failed to write");
                            }
                            if let Some(ref mut f) = features {
                                f.push(peakgroup.clone());
                            }
                            let idx = current_index;
                            current_index += 1;
                            idx
                        });

                        intervals.push(Interval { start, stop: end, val: *group_index });
                    } else {
                        intervals.push(Interval { start, stop: end, val: index - skipped_lines});
                        if let Some(ref mut w) = writer {
                            writeln!(w, "{}:{}-{}", chromosome, start, end)?;
                        }
                        if let Some(ref mut f) = features {
                            f.push(format!("{}:{}-{}", chromosome, start, end));
                        }
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
    if let Some(mut w) = writer {
        w.finish().map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    }

    Ok((total_peaks, lapper_map, features))
}

#[cfg(feature = "hdf5")]
fn write_hdf5(
    file_path: &Path,
    mut all_counts: Vec<(u32, u32, u16)>, // cell_idx, peak_idx, count
    total_peaks: usize,
    total_cells: usize,
    features: Vec<String>,
    barcodes: Vec<String>,
) -> io::Result<()> {
    // Sort primarily by cell index (column), secondarily by peak position (row) for CSC format
    all_counts.sort_unstable_by_key(|&(c, p, _)| (c, p));
    
    let mut data = Vec::with_capacity(all_counts.len());
    let mut indices = Vec::with_capacity(all_counts.len());
    let mut indptr = Vec::with_capacity(total_cells + 1);
    
    indptr.push(0);
    let mut current_col = 0;
    
    for &(cell, peak, value) in &all_counts {
        // Fill empty columns
        while current_col < cell {
            indptr.push(data.len() as u32);
            current_col += 1;
        }
        
        indices.push(peak);
        data.push(value);
    }
    
    // Fill remaining trailing columns
    while current_col < total_cells as u32 {
        indptr.push(data.len() as u32);
        current_col += 1;
    }
    
    let file = H5File::create(file_path).map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    let matrix_group = file.create_group("matrix").map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        
    let data_arr = Array1::from(data);
    matrix_group.new_dataset_builder().with_data(&data_arr).create("data").map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        
    let indices_arr = Array1::from(indices);
    matrix_group.new_dataset_builder().with_data(&indices_arr).create("indices").map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        
    let indptr_arr = Array1::from(indptr);
    matrix_group.new_dataset_builder().with_data(&indptr_arr).create("indptr").map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        
    let shape_arr = Array1::from(vec![total_peaks as u32, total_cells as u32]);
    matrix_group.new_dataset_builder().with_data(&shape_arr).create("shape").map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        
    // HDF5 string conversion natively uses Variable Length Strings over `&str`
    let barcodes_ref: Vec<VarLenUnicode> = barcodes.iter().map(|s| s.parse().unwrap()).collect();
    let barcodes_arr = Array1::from(barcodes_ref);
    matrix_group.new_dataset_builder().with_data(&barcodes_arr).create("barcodes").map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        
    let features_group = matrix_group.create_group("features").map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        
    let features_ref: Vec<VarLenUnicode> = features.iter().map(|s| s.parse().unwrap()).collect();
    let features_arr = Array1::from(features_ref);
    features_group.new_dataset_builder().with_data(&features_arr).create("id").map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    features_group.new_dataset_builder().with_data(&features_arr).create("name").map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        
    let feature_type: Vec<VarLenUnicode> = vec!["Peaks".parse().unwrap(); features.len()];
    let feature_type_arr = Array1::from(feature_type);
    features_group.new_dataset_builder().with_data(&feature_type_arr).create("feature_type").map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        
    let genome: Vec<VarLenUnicode> = vec!["GRCh38".parse().unwrap(); features.len()];
    let genome_arr = Array1::from(genome);
    features_group.new_dataset_builder().with_data(&genome_arr).create("genome").map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        
    Ok(())
}