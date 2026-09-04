use crate::intervals::{seek_position, SeekCursor};
use bio::io::gff;
use bio_types::strand::Strand;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use lexical_core::parse;
use log::info;
use rust_lapper::{Interval, Lapper};
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::SmallVec;
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;
use std::sync::OnceLock;

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

// TSS enrichment parameters, following the ENCODE definition quoted above: the
// distribution extends TSS_WINDOW either side of the TSS, the background is the mean
// depth over the FLANK_SIZE at each end of that distribution, and the score is the
// fold change of the centre over that background.
//
// ENCODE takes "the signal value at the center" -- a single position. That is written
// for a bulk aggregate profile; per single cell the centre base collects only a
// handful of insertions and the score bottoms out at zero for real cells, so the
// centre is averaged over CENTER_SIZE bases instead.
const TSS_WINDOW: u32 = 2000; // distribution extends 2000bp either side (4000bp total)
const FLANK_SIZE: u32 = 100; // 100bp at each end flank (200bp of averaged data)
const CENTER_SIZE: u32 = 100; // central 100bp; a single position is too sparse per cell
const PROMOTER_HALF: u32 = 1000; // promoter window is TSS +/-1000bp (2000bp total)
const NUCLEOSOME: u32 = 147; // length of DNA wrapped around nucleosome

// Total width feeding tss_flank, used with CENTER_SIZE to convert the raw counts into
// per-base depths before taking their ratio in write_results.
const FLANK_WIDTH: u32 = 2 * FLANK_SIZE; // upstream + downstream end flank

// Floor on the flank depth (insertions per base) used to normalize the score. A cell
// with only a handful of flank insertions otherwise divides by a near-zero background
// and scores arbitrarily high
const MIN_FLANK_DEPTH: f32 = 0.1;

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

/// Per-chromosome interval sets built from an annotation file: the TSS centre and
/// flank windows used for the enrichment score, and a separate flattened promoter set.
type TssRegions = (
    FxHashMap<String, Lapper<u32, usize>>,
    FxHashMap<String, Lapper<u32, usize>>,
);

fn build_lappers(
    trees: FxHashMap<String, Vec<Interval<u32, usize>>>,
) -> FxHashMap<String, Lapper<u32, usize>> {
    trees
        .into_iter()
        .map(|(chr, intervals)| (chr, Lapper::new(intervals)))
        .collect()
}

/// Build lappers with overlapping intervals merged into single spans.
///
/// Flattening makes each base covered only once, so an insertion can only be
/// counted a single time even if there are overlapping TSS sites
fn flatten_lappers(
    trees: FxHashMap<String, Vec<Interval<u32, usize>>>,
) -> FxHashMap<String, Lapper<u32, usize>> {
    trees
        .into_iter()
        .map(|(chr, intervals)| {
            let mut lapper = Lapper::new(intervals);
            lapper.merge_overlaps();
            (chr, lapper)
        })
        .collect()
}

#[derive(Debug)]
pub struct FragmentCounts {
    pub tss_flank: u32,
    pub tss_center: u32,
    pub mononucleosome: u32,
    pub nucleosome_free: u32,
    pub total_fragments: u32,
    pub total_insertions: u32,
    pub promoter_insertions: u32,
    pub mitochondrial_fragments: u32,
    pub x_fragments: u32,
    pub y_fragments: u32,
}

pub fn tss_enrichment(
    fragments: &str,
    annotation: &str,
    annotation_is_gff: bool,
    cells: Option<&str>,
    outfile: &str,
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

    // Load cell barcodes if provided
    let cell_filter: Option<FxHashSet<Box<str>>> = match cells {
        Some(cells_path) => {
            let cell_file = Path::new(cells_path)
                .canonicalize()
                .expect("Can't find path to input cell file");
            info!("Received cell file: {:?}", cell_file);
            let reader = crate::reader::open_maybe_gzipped(&cell_file)?;
            let mut set = FxHashSet::default();
            for line in reader.lines() {
                set.insert(line?.into_boxed_str());
            }
            info!("Loaded {} cell barcodes for filtering", set.len());
            Some(set)
        }
        None => {
            info!("No cell barcode filter provided, processing all cells");
            None
        }
    };

    let annotation = annotation.to_path_buf();

    let (tss_regions, promoter_regions) = if annotation_is_gff {
        extract_tss_gff(&annotation)?
    } else {
        extract_tss_bed(&annotation)?
    };

    // iterate over fragments
    // count overlaps with tss regions and flanking regions
    // store in hashmap of FragmentCounts
    let mut results_by_cell: FxHashMap<Box<str>, FragmentCounts> = FxHashMap::default();

    // Spawn reader thread for decompression
    let frag_file = frag_file.to_path_buf();
    let (reader_handle, rx, pool_tx) = crate::reader::spawn_fragment_reader(&frag_file);

    let mut current_chrom = String::new();
    let mut current_lapper: Option<&Lapper<u32, usize>> = None;
    let mut start_cursor = SeekCursor::new();
    let mut end_cursor = SeekCursor::new();
    let mut current_promoters: Option<&Lapper<u32, usize>> = None;
    let mut prom_start_cursor = SeekCursor::new();
    let mut prom_end_cursor = SeekCursor::new();

    // Process chunks from the channel
    for mut chunk in rx {
        for line in chunk.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }

            let mut iter = line.splitn(5, |&b| b == b'\t');
            let seqname_bytes = match iter.next() {
                Some(b) => b,
                None => continue,
            };
            let start_bytes = match iter.next() {
                Some(b) => b,
                None => continue,
            };
            let end_bytes = match iter.next() {
                Some(b) => b,
                None => continue,
            };
            let mut barcode_bytes = match iter.next() {
                Some(b) => b,
                None => continue,
            };

            // Trim trailing \r
            if barcode_bytes.ends_with(b"\r") {
                barcode_bytes = &barcode_bytes[..barcode_bytes.len() - 1];
            }

            let cell_barcode = unsafe { std::str::from_utf8_unchecked(barcode_bytes) };

            // Filter by cell barcodes if a filter was provided
            if let Some(ref filter) = cell_filter {
                if !filter.contains(cell_barcode) {
                    continue;
                }
            }

            let fc = if let Some(fc) = results_by_cell.get_mut(cell_barcode) {
                fc
            } else {
                results_by_cell.insert(
                    cell_barcode.into(),
                    FragmentCounts {
                        tss_flank: 0,
                        tss_center: 0,
                        total_fragments: 0,
                        total_insertions: 0,
                        promoter_insertions: 0,
                        mononucleosome: 0,
                        nucleosome_free: 0,
                        mitochondrial_fragments: 0,
                        x_fragments: 0,
                        y_fragments: 0,
                    },
                );
                results_by_cell.get_mut(cell_barcode).unwrap()
            };
            fc.total_fragments += 1;

            // Create intervals from fragment entry
            let seqname = unsafe { std::str::from_utf8_unchecked(seqname_bytes) };

            // check if chromosome changed
            // if so, update lapper
            if seqname != current_chrom {
                current_chrom = seqname.to_string();
                current_lapper = tss_regions.get(&current_chrom);
                start_cursor.reset();
                end_cursor.reset();
                current_promoters = promoter_regions.get(&current_chrom);
                prom_start_cursor.reset();
                prom_end_cursor.reset();
            }

            // update mito count
            if is_mito_chr(seqname) {
                fc.mitochondrial_fragments += 1;
            } else if is_chr_x(seqname) {
                fc.x_fragments += 1;
            } else if is_chr_y(seqname) {
                fc.y_fragments += 1;
            }

            let startpos = match parse::<u32>(start_bytes) {
                Ok(num) => num,
                Err(_) => continue,
            };

            let endpos = match parse::<u32>(end_bytes) {
                Ok(num) => num,
                Err(_) => continue,
            };

            let fragment_width = endpos - startpos;

            // record nucleosome signal information
            if fragment_width < NUCLEOSOME {
                fc.nucleosome_free += 1;
            } else if fragment_width < (NUCLEOSOME * 2) {
                fc.mononucleosome += 1;
            }

            // re-borrow fc mutably after the continue branches above
            let fc = results_by_cell.get_mut(cell_barcode).unwrap();

            // each fragment contributes two Tn5 insertions, at its start and its end
            fc.total_insertions += 2;

            // Promoter regions are flattened, so each insertion matches at most once
            if let Some(promoters) = &current_promoters {
                for (pos, cur) in [
                    (startpos, &mut prom_start_cursor),
                    (endpos, &mut prom_end_cursor),
                ] {
                    if seek_position(promoters, pos, cur).next().is_some() {
                        fc.promoter_insertions += 1;
                    }
                }
            }

            if let Some(lapper) = &current_lapper {
                // Check for overlaps at start position
                for interval in seek_position(lapper, startpos, &mut start_cursor) {
                    let val = interval.val as u32; // 0 for tss, 1 for flank
                    if val == 0 {
                        fc.tss_center += 1;
                    } else {
                        fc.tss_flank += 1;
                    }
                }

                // Check for overlaps at end position
                for interval in seek_position(lapper, endpos, &mut end_cursor) {
                    let val = interval.val as u32; // 0 for tss, 1 for flank
                    if val == 0 {
                        fc.tss_center += 1;
                    } else {
                        fc.tss_flank += 1;
                    }
                }
            }
        }
        chunk.clear();
        let _ = pool_tx.send(chunk);
    }

    // Wait for reader thread to complete
    reader_handle.join().expect("Reader thread panicked");

    // write to file
    write_results(&results_by_cell, &outfile)?;

    Ok(())
}

/// Extract TSS and promoter regions from a BED file of TSS positions.
///
/// Returns the enrichment-score intervals keyed by chromosome, where an interval value
/// of 0 marks the centre window and 1 a flanking window, plus a separate flattened set
/// of promoter intervals. See [`TssRegions`].
fn extract_tss_bed(bed_path: &Path) -> Result<TssRegions, Box<dyn std::error::Error>> {
    let reader = crate::reader::open_maybe_gzipped(bed_path)?;

    let mut chromosome_trees: FxHashMap<String, Vec<Interval<u32, usize>>> = FxHashMap::default();
    let mut promoter_trees: FxHashMap<String, Vec<Interval<u32, usize>>> = FxHashMap::default();

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

                    let intervals = chromosome_trees.entry(chromosome.clone()).or_default();

                    // assume bed file is TSS base position only
                    // just take start position

                    // Centre of the distribution: the central CENTER_SIZE bases, centred on
                    // the TSS. Deriving tss_end from tss_start keeps the width exactly
                    // CENTER_SIZE. saturating_sub guards TSSes near the chromosome start.
                    let tss_start = start.saturating_sub(CENTER_SIZE / 2);
                    let tss_end = tss_start + CENTER_SIZE;

                    // the FLANK_SIZE at each end of the +/-TSS_WINDOW distribution
                    let flank_upstream_start = start.saturating_sub(TSS_WINDOW);
                    let flank_upstream_end = start.saturating_sub(TSS_WINDOW - FLANK_SIZE);

                    let flank_downstream_start = start + TSS_WINDOW - FLANK_SIZE;
                    let flank_downstream_end = start + TSS_WINDOW;

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

                    promoter_trees
                        .entry(chromosome)
                        .or_default()
                        .push(Interval {
                            start: start.saturating_sub(PROMOTER_HALF),
                            stop: start + PROMOTER_HALF,
                            val: 0,
                        });
                } else {
                    return Err(Box::new(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Less than three fields",
                    )));
                }
            }
            Err(_) => {
                return Err(Box::new(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Error reading line",
                )));
            }
        }
    }

    Ok((
        build_lappers(chromosome_trees),
        flatten_lappers(promoter_trees),
    ))
}

/// Extract TSS and promoter regions from a GFF file of gene annotations.
///
/// Returns the enrichment-score intervals keyed by chromosome, where an interval value
/// of 0 marks the centre window and 1 a flanking window, plus a separate flattened set
/// of promoter intervals. See [`TssRegions`].
fn extract_tss_gff(gff_path: &Path) -> Result<TssRegions, Box<dyn std::error::Error>> {
    let file: Box<dyn io::Read> = if crate::reader::is_gzipped(gff_path)? {
        Box::new(MultiGzDecoder::new(File::open(gff_path)?))
    } else {
        Box::new(File::open(gff_path)?)
    };
    let mut reader = gff::Reader::new(file, gff::GffType::GFF3);

    // hashmap of TSS intervals for each chromosome
    let mut chromosome_trees: FxHashMap<String, Vec<Interval<u32, usize>>> = FxHashMap::default();
    let mut promoter_trees: FxHashMap<String, Vec<Interval<u32, usize>>> = FxHashMap::default();

    eprintln!("Reading gene annotations");
    for record in reader.records() {
        let rec = record?;

        // Filter for transcript or gene features
        if rec.feature_type() != "transcript" {
            continue;
        }

        let chromosome = rec.seqname().to_string();
        let intervals = chromosome_trees.entry(chromosome.clone()).or_default();

        let mut start = *rec.start() as u32;
        let strand = rec.strand().unwrap_or(Strand::Unknown);

        // if strand is -, start = end
        if strand == Strand::Reverse {
            start = *rec.end() as u32;
        }

        // Centre of the distribution: the central CENTER_SIZE bases, centred on the TSS.
        // Deriving tss_end from tss_start keeps the width exactly CENTER_SIZE.
        // saturating_sub guards TSSes near the chromosome start.
        let tss_start = start.saturating_sub(CENTER_SIZE / 2);
        let tss_end = tss_start + CENTER_SIZE;

        // the FLANK_SIZE at each end of the +/-TSS_WINDOW distribution
        let flank_upstream_start = start.saturating_sub(TSS_WINDOW);
        let flank_upstream_end = start.saturating_sub(TSS_WINDOW - FLANK_SIZE);

        let flank_downstream_start = start + TSS_WINDOW - FLANK_SIZE;
        let flank_downstream_end = start + TSS_WINDOW;

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

        promoter_trees
            .entry(chromosome)
            .or_default()
            .push(Interval {
                start: start.saturating_sub(PROMOTER_HALF),
                stop: start + PROMOTER_HALF,
                val: 0,
            });
    }

    Ok((
        build_lappers(chromosome_trees),
        flatten_lappers(promoter_trees),
    ))
}

/// Writes a TSV of fragment counts to a gzip-compressed file.
fn write_results(
    counts: &FxHashMap<Box<str>, FragmentCounts>,
    output_path: &str,
) -> std::io::Result<()> {
    let file = File::create(output_path)?;
    let gz = GzEncoder::new(file, Compression::default());
    let mut writer = BufWriter::new(gz);

    // Write header
    writeln!(
        writer,
        "cell_barcode\tTSS_flank\tTSS_center\ttotal_fragments\tTSS_enrichment\tFIP\tNucleosome_signal\tMito_fragments\tMito_fraction\tchrX_fragments\tchrX_fraction\tchrY_fragments\tchrY_fraction"
    )?;

    // Write records
    for (barcode, fc) in counts {
        let mut nucleosome_signal: f32 = 0.0;
        let mut fip: f32 = 0.0;
        let mut mito_fraction: f32 = 0.0;
        let mut x_fraction: f32 = 0.0;
        let mut y_fraction: f32 = 0.0;

        let center_depth = fc.tss_center as f32 / CENTER_SIZE as f32;
        let flank_depth = (fc.tss_flank as f32 / FLANK_WIDTH as f32).max(MIN_FLANK_DEPTH);
        let tsse = center_depth / flank_depth;
        if fc.nucleosome_free > 0 {
            nucleosome_signal = fc.mononucleosome as f32 / fc.nucleosome_free as f32;
        }
        if fc.total_insertions > 0 {
            fip = fc.promoter_insertions as f32 / fc.total_insertions as f32;
        }
        if fc.total_fragments > 0 {
            mito_fraction = fc.mitochondrial_fragments as f32 / fc.total_fragments as f32;
            x_fraction = fc.x_fragments as f32 / fc.total_fragments as f32;
            y_fraction = fc.y_fragments as f32 / fc.total_fragments as f32;
        }

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            barcode,
            fc.tss_flank,
            fc.tss_center,
            fc.total_fragments,
            format!("{:.4}", tsse),
            format!("{:.4}", fip),
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
