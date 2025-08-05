use bio::io::gff;
use bio_types::strand::Strand;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use lexical_core::parse;
use log::error;
use log::info;
use rust_lapper::{Interval, Lapper};
use rustc_hash::FxHashMap;
use smallvec::SmallVec;
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::mpsc;
use std::sync::OnceLock;
use std::thread;

// https://www.encodeproject.org/data-standards/terms/
/* The reads around a reference set of TSSs are collected to form an aggregate
distribution of reads centered on the TSSs and extending to 2000 bp in either
direction (for a total of 4000bp). This distribution is then normalized by taking
the average read depth in the 100 bps at each of the end flanks of the
distribution (for a total of 200bp of averaged data) and calculating a fold change
at each position over that average read depth. This means that the flanks should
start at 1, and if there is high read signal at transcription start sites (highly
open regions of the genome) there should be an increase in signal up to a peak in the
middle. We take the signal value at the center of the distribution after this
normalization as our TSS enrichment metric. Used to evaluate ATAC-seq. */

// TSS enrichment parameters based on ENCODE standard
const TSS_WINDOW: u32 = 2000; // ±2000bp around TSS (total 4000bp)
const FLANK_SIZE: u32 = 100; // 100bp at each end for normalization
const CENTER_SIZE: u32 = 100; // 100bp at the center for score
const NUCLEOSOME: u32 = 147; // length of DNA wrapped around nucleosome

static MITO_CHROMS: OnceLock<HashSet<&'static str>> = OnceLock::new();
static CHRX_NAMES: OnceLock<HashSet<&'static str>> = OnceLock::new();
static CHRY_NAMES: OnceLock<HashSet<&'static str>> = OnceLock::new();

fn get_mito_chroms() -> &'static HashSet<&'static str> {
    MITO_CHROMS.get_or_init(|| {
        [
            "chrM",
            "MT",
            "chrMT",
            "M",
            "mitochondrion_genome",
            "chrMito",
            "Mito",
            "mtDNA",
        ]
        .iter()
        .copied()
        .collect()
    })
}

fn get_chr_x() -> &'static HashSet<&'static str> {
    CHRX_NAMES.get_or_init(|| ["chrX", "X"].iter().copied().collect())
}

fn get_chr_y() -> &'static HashSet<&'static str> {
    CHRY_NAMES.get_or_init(|| ["chrY", "Y"].iter().copied().collect())
}

pub fn is_mito_chr(chr: &str) -> bool {
    get_mito_chroms().contains(chr)
}

pub fn is_chr_x(chr: &str) -> bool {
    get_chr_x().contains(chr)
}

pub fn is_chr_y(chr: &str) -> bool {
    get_chr_y().contains(chr)
}

#[derive(Debug)]
pub struct FragmentCounts {
    pub tss_flank: u32,
    pub tss_center: u32,
    pub mononucleosome: u32,
    pub nucleosome_free: u32,
    pub total_fragments: u32,
    pub mitochondrial_fragments: u32,
    pub x_fragments: u32,
    pub y_fragments: u32,
}

pub fn tss_enrichment(
    fragments: &str,
    annotation: &str,
    annotation_is_gff: bool,
    outfile: &str
) -> Result<(), Box<dyn Error>> {
    // get TSS positions from the annotation file
    let annotation = Path::new(annotation)
        .canonicalize()
        .expect("Can't find path to input annotation file");
    info!("Received annotation file: {:?}", annotation);

    let frag_file = Path::new(fragments)
        .canonicalize()
        .expect("Can't find path to input fragment file");
    info!("Received fragment file: {:?}", frag_file);

    info!("Received output file path: {:?}", outfile);

    let annotation = annotation.to_path_buf();

    let tss_regions = if annotation_is_gff {
        extract_tss_gff(&annotation)?
    } else {
        extract_tss_bed(&annotation)?
    };

    // iterate over fragments
    // count overlaps with tss regions and flanking regions
    // store in hashmap of FragmentCounts
    let mut results_by_cell: FxHashMap<String, FragmentCounts> = FxHashMap::default();

    // Create channels for communication between threads
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
                            let chunk_to_send =
                                std::mem::replace(&mut fragments, Vec::with_capacity(CHUNK_SIZE));
                            if tx.send(chunk_to_send).is_err() {
                                break;
                            }
                        }
                    }
                }
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

    let mut current_chrom = String::new();
    let mut current_lapper: Option<&Lapper<u32, usize>> = None;
    let mut cursor = 0;

    // Process chunks from the channel
    for chunk in rx {
        for line in chunk {
            // Parse BED entry
            let fields = line.split('\t').collect::<SmallVec<[&str; 10]>>();
            if fields.len() < 4 {
                continue;
            }

            // TODO enable specifying list of cell barcodes to include
            let cell_barcode: &str = fields[3];
            results_by_cell
                .entry(cell_barcode.to_string())
                .or_insert(FragmentCounts {
                    tss_flank: 0,
                    tss_center: 0,
                    total_fragments: 0,
                    mononucleosome: 0,
                    nucleosome_free: 0,
                    mitochondrial_fragments: 0,
                    x_fragments: 0,
                    y_fragments: 0,
                });
            results_by_cell
                .get_mut(cell_barcode)
                .unwrap()
                .total_fragments += 1;

            // Create intervals from fragment entry
            let seqname: &str = fields[0];

            // check if chromosome changed
            // if so, update lapper
            if seqname != current_chrom {
                current_chrom = seqname.to_string();
                current_lapper = tss_regions.get(&current_chrom);
                cursor = 0;
            }

            // update mito count
            if is_mito_chr(seqname) {
                results_by_cell
                    .get_mut(cell_barcode)
                    .unwrap()
                    .mitochondrial_fragments += 1;
            } else if is_chr_x(seqname) {
                results_by_cell.get_mut(cell_barcode).unwrap().x_fragments += 1;
            } else if is_chr_y(seqname) {
                results_by_cell.get_mut(cell_barcode).unwrap().y_fragments += 1;
            }

            let start_str = fields[1];
            let end_str = fields[2];
            let startpos = match parse::<u32>(start_str.as_bytes()) {
                Ok(num) => num,
                Err(_) => continue,
            };

            let endpos = match parse::<u32>(end_str.as_bytes()) {
                Ok(num) => num,
                Err(_) => continue,
            };

            let fragment_width = endpos - startpos;

            // record nucleosome signal information
            if fragment_width < NUCLEOSOME {
                results_by_cell
                    .get_mut(cell_barcode)
                    .unwrap()
                    .nucleosome_free += 1;
            } else if fragment_width < (NUCLEOSOME * 2) {
                results_by_cell
                    .get_mut(cell_barcode)
                    .unwrap()
                    .mononucleosome += 1;
            }

            if let Some(lapper) = &current_lapper {
                if lapper.intervals.len() == 1 {
                    cursor = 0;
                }

                // Check for overlaps at start position
                for interval in lapper.seek(startpos, startpos + 1, &mut cursor) {
                    let val = interval.val as u32; // 0 for tss, 1 for flank
                    if val == 0 {
                        results_by_cell.get_mut(cell_barcode).unwrap().tss_center += 1;
                    } else {
                        results_by_cell.get_mut(cell_barcode).unwrap().tss_flank += 1;
                    }
                }

                // Check for overlaps at end position
                for interval in lapper.seek(endpos, endpos + 1, &mut cursor) {
                    let val = interval.val as u32; // 0 for tss, 1 for flank
                    if val == 0 {
                        results_by_cell.get_mut(cell_barcode).unwrap().tss_center += 1;
                    } else {
                        results_by_cell.get_mut(cell_barcode).unwrap().tss_flank += 1;
                    }
                }
            }
        }
    }

    // Wait for reader thread to complete
    reader_handle.join().expect("Reader thread panicked");

    // write to file
    let _ = write_results(&results_by_cell, &outfile);

    Ok(())
}

/// Extract TSS regions from a BED file
///
/// BED file contains TSS positions
/// returns a hashmap of chromosome names to vectors of Interval objects
/// each interval is a TSS region or flanking region
/// the value of the interval determines if it is a TSS region or flanking region
/// 0 indicates a TSS region, 1 indicates a flanking region
fn extract_tss_bed(
    bed_path: &Path,
) -> Result<FxHashMap<String, Lapper<u32, usize>>, Box<dyn std::error::Error>> {
    let file = File::open(bed_path)?;
    let reader: Box<dyn BufRead> =
        if bed_path.extension().and_then(|ext| ext.to_str()) == Some("gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

    let mut chromosome_trees: FxHashMap<String, Vec<Interval<u32, usize>>> = FxHashMap::default();

    eprintln!("Reading TSS positions");
    for line in reader.lines() {
        match line {
            Ok(line) => {
                if line.starts_with('#') {
                    continue;
                }
                let fields = line.split('\t').collect::<SmallVec<[&str; 10]>>();
                if fields.len() >= 3 {
                    let chromosome = fields[0].to_string();
                    let start: u32 = match parse::<u32>(fields[1].trim().as_bytes()) {
                        Ok(num) => num,
                        Err(_) => {
                            return Err(Box::new(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "Failed to parse start position",
                            )));
                        }
                    };

                    let intervals = chromosome_trees
                        .entry(chromosome.clone())
                        .or_insert_with(Vec::new);

                    // assume bed file is TSS base position only
                    // just take start position

                    let tss_start = start - (CENTER_SIZE / 2);
                    let tss_end = start + (CENTER_SIZE / 2);

                    let flank_upstream_start = start - TSS_WINDOW - (FLANK_SIZE / 2);
                    let flank_upstream_end = start - TSS_WINDOW + (FLANK_SIZE / 2);

                    let flank_downstream_start = start + TSS_WINDOW - (FLANK_SIZE / 2);
                    let flank_downstream_end = start + TSS_WINDOW + (FLANK_SIZE / 2);

                    intervals.push(Interval {
                        start: tss_start,
                        stop: tss_end,
                        val: 0,
                    });
                    intervals.push(Interval {
                        start: flank_upstream_start,
                        stop: flank_upstream_end,
                        val: 1,
                    });
                    intervals.push(Interval {
                        start: flank_downstream_start,
                        stop: flank_downstream_end,
                        val: 1,
                    });
                } else {
                    return Err(Box::new(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Less than three fields"
                    )));
                }
            }
            Err(_) => {
                return Err(Box::new(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Error reading line"
                )));
            }
        }
    }

    let lapper_map = chromosome_trees
        .into_iter()
        .map(|(chr, intervals)| (chr, Lapper::new(intervals)))
        .collect();

    Ok(lapper_map)
}

/// Extract TSS regions from a GFF file
///
/// returns a hashmap of chromosome names to vectors of Interval objects
/// each interval is a TSS region or flanking region
/// the value of the interval determines if it is a TSS region or flanking region
/// 0 indicates a TSS region, 1 indicates a flanking region
fn extract_tss_gff(
    gff_path: &Path,
) -> Result<FxHashMap<String, Lapper<u32, usize>>, Box<dyn std::error::Error>> {
    let file: Box<dyn io::Read> = if gff_path.extension().and_then(|ext| ext.to_str()) == Some("gz")
    {
        Box::new(MultiGzDecoder::new(File::open(gff_path)?))
    } else {
        Box::new(File::open(gff_path)?)
    };
    let mut reader = gff::Reader::new(file, gff::GffType::GFF3);

    // hashmap of peak intervals for each chromosome
    let mut chromosome_trees: FxHashMap<String, Vec<Interval<u32, usize>>> = FxHashMap::default();

    eprintln!("Reading gene annotations");
    for record in reader.records() {
        let rec = record?;

        // Filter for transcript or gene features
        if rec.feature_type() != "transcript" {
            continue;
        }

        let chromosome = rec.seqname().to_string();
        let intervals = chromosome_trees
            .entry(chromosome.clone())
            .or_insert_with(Vec::new);

        let mut start = *rec.start() as u32;
        let strand = rec.strand().unwrap_or(Strand::Unknown);

        // if strand is -, start = end
        if strand == Strand::Reverse {
            start = *rec.end() as u32;
        }

        let tss_start = start - (CENTER_SIZE / 2);
        let tss_end = start + (CENTER_SIZE / 2);

        let flank_upstream_start = start - TSS_WINDOW - (FLANK_SIZE / 2);
        let flank_upstream_end = start - TSS_WINDOW + (FLANK_SIZE / 2);

        let flank_downstream_start = start + TSS_WINDOW - (FLANK_SIZE / 2);
        let flank_downstream_end = start + TSS_WINDOW + (FLANK_SIZE / 2);

        intervals.push(Interval {
            start: tss_start,
            stop: tss_end,
            val: 0,
        });
        intervals.push(Interval {
            start: flank_upstream_start,
            stop: flank_upstream_end,
            val: 1,
        });
        intervals.push(Interval {
            start: flank_downstream_start,
            stop: flank_downstream_end,
            val: 1,
        });
    }

    let lapper_map = chromosome_trees
        .into_iter()
        .map(|(chr, intervals)| (chr, Lapper::new(intervals)))
        .collect();

    Ok(lapper_map)
}

/// Writes a TSV of fragment counts to a gzip-compressed file.
fn write_results(
    counts: &FxHashMap<String, FragmentCounts>,
    output_path: &str,
) -> std::io::Result<()> {
    let file = File::create(output_path)?;
    let gz = GzEncoder::new(file, Compression::default());
    let mut writer = BufWriter::new(gz);

    // Write header
    writeln!(
        writer,
        "cell_barcode\tTSS_flank\tTSS_center\ttotal_fragments\tTSS_enrichment\tFRiP\tNucleosome_signal\tMito_fragments\tMito_fraction\tchrX_fragments\tchrX_fraction\tchrY_fragments\tchrY_fraction"
    )?;

    // Write records
    for (barcode, fc) in counts {
        let mut tsse: f32 = 0.0;
        let mut nucleosome_signal: f32 = 0.0;
        if fc.tss_flank > 0 {
            // avoid division by zero
            // if there are zero counts in flank, must be extremely low total counts
            // TSSe = 0 ok to record
            tsse = fc.tss_center as f32 / fc.tss_flank as f32;
        }
        if fc.nucleosome_free > 0 {
            nucleosome_signal = fc.mononucleosome as f32 / fc.nucleosome_free as f32;
        }

        let frip: f32 = fc.tss_center as f32 / fc.total_fragments as f32;
        let mito_fraction: f32 = fc.mitochondrial_fragments as f32 / fc.total_fragments as f32;
        let x_fraction: f32 = fc.x_fragments as f32 / fc.total_fragments as f32;
        let y_fraction: f32 = fc.y_fragments as f32 / fc.total_fragments as f32;

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            barcode,
            fc.tss_flank,
            fc.tss_center,
            fc.total_fragments,
            format!("{:.4}", tsse),
            format!("{:.4}", frip),
            format!("{:.4}", nucleosome_signal),
            fc.mitochondrial_fragments,
            format!("{:.4}", mito_fraction),
            fc.x_fragments,
            format!("{:.4}", x_fraction),
            fc.y_fragments,
            format!("{:.4}", y_fraction),
        )?;
    }

    writer.flush()?;
    Ok(())
}
