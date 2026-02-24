#[cfg(not(target_os = "windows"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: std::alloc::System = std::alloc::System;

use clap::{Parser, Subcommand, ArgGroup};
use std::error::Error;
use std::path::PathBuf;

mod f2m;
mod cellselect;
mod filter;
mod man;
mod qc;
mod reader;

#[derive(Parser)]
#[command(
    name = env!("CARGO_PKG_NAME"),
    author = env!("CARGO_PKG_AUTHORS"),
    about = env!("CARGO_PKG_DESCRIPTION"),
    version = env!("CARGO_PKG_VERSION")
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[command(
        about = "Create a feature x cell matrix",
        long_about = "Create a feature x cell matrix from fragments file and BED regions"
    )]
    Matrix {
        #[arg(
            short, 
            long,
            value_name = "FILE",
            help = "Path to the fragment file",
            long_help = "Path to the fragment file. The file should be gzipped and contain \
                         chromosome, start, end, and cell barcode columns tab-separated."
        )]
        fragments: String,

        #[arg(
            short,
            long,
            value_name = "FILE",
            help = "BED file containing non-overlapping genomic regions to quantify",
            long_help = "Input BED file containing regions to quantify. \
                         Should be a tab-delimited file containing chromosome, start, end \
                         as the first three columns (valid BED format). The file can be \
                         gzip-compressed."
        )]
        bed: String,

        #[arg(
            short,
            long,
            value_name = "FILE",
            help = "File containing cell barcodes to include",
            long_help = "Text file containing cell barcodes to include in the output matrix. \
                         This file should contain one cell barcode string per line, and may be \
                         gzip-compressed."
        )]
        cells: String,

        #[arg(
            short,
            long,
            value_name = "PATH",
            help = "Output directory name",
            long_help = "Output directory name. Directory will be created if it does not exist. \
                         The output directory will contain matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz"
        )]
        outdir: String,

        #[arg(
            short,
            long,
            value_name = "NUMBER",
            help = "Number of compression threads to use",
            default_value = "2",
            value_parser = clap::value_parser!(usize)
        )]
        threads: usize,

        #[arg(
            long,
            help = "Use paired insertion counting",
            long_help = "Use paired insertion counting, as described by Zhen Miao & Junhyong Kim \
                         (2023, Nature Methods; https://doi.org/10.1038/s41592-023-02103-7).",
            action = clap::ArgAction::SetTrue)]
        pic: bool,

        #[arg(
            long,
            help = "Group regions by variable in fourth BED column",
            long_help = "Sum counts across regions sharing a variable stored in the fourth \
                         column of the BED file. This enables quantification of non-contiguous \
                         regions as a single feature in the output matrix.",
            action = clap::ArgAction::SetTrue)]
        group: bool,

        #[arg(
            long,
            help = "Output the matrix in 10x Genomics HDF5 format",
            long_help = "Outputs matrix.h5 to the output directory instead of the standard \
                         Matrix Market (.mtx) and TSV files.",
            action = clap::ArgAction::SetTrue)]
        h5: bool,
    },

    #[command(
        about = "Count number of fragments per cell barcode",
        long_about = "Count the total number of fragments for each cell barcode in the fragments file"
    )]
    Count {
        #[arg(
            short, 
            long, 
            value_name = "FILE",
            help = "Path to the fragment file",
            long_help = "Path to the fragment file. The file should be gzipped and contain \
                         chromosome, start, end, and cell barcode columns tab-separated."
        )]
        fragments: String,

        #[arg(
            short, 
            long, 
            value_name = "FILE", 
            help = "Name of output file",
            long_help = "Name of output file. The file will contain each cell barcode and its total fragment count, tab-separated."
        )]
        outfile: String,

        #[arg(
            short, 
            long, 
            value_name = "NUMBER",
            help = "Minimum number of fragments for a cell to be included",
            long_help = "Sets the minimum number of fragments a cell must have to be included in the output.\
                         Cells with fewer fragments than this threshold will be filtered out. Cannot be used with --ncells",
            conflicts_with = "ncells",
            value_parser = clap::value_parser!(usize)
        )]
        threshold: Option<usize>,

        #[arg(
            short, 
            long, 
            value_name = "NUMBER",
            help = "Number of top cells to select",
            long_help = "Select this many cells with the highest fragment counts. Cannot be used together with --threshold.",
            conflicts_with = "threshold",
            value_parser = clap::value_parser!(usize)
        )]
        ncells: Option<usize>,
    },

    #[command(
        about = "Subset a fragment file to include only specified cell barcodes",
        long_about = "Filter a fragments file to keep only fragments from cells listed in a cell barcodes file. \n\
                      The filtered fragments are written to stdout as plain text and can be piped into a compression \
                      command, for example: \nfragtk filter -f <fragments> -c <cells> | bgzip -c > filtered.tsv.gz"
    )]
    Filter {
        #[arg(
            short, 
            long,
            value_name = "FILE",
            help = "Path to the fragment file",
            long_help = "Path to the fragment file. The file should be gzipped and contain \
                         chromosome, start, end, and cell barcode columns tab-separated."
        )]
        fragments: String,

        #[arg(
            short, 
            long,
            value_name = "FILE",
            help = "File containing cell barcodes to include",
            long_help = "Text file containing cell barcodes to include in the output matrix. \
                         This file should contain one cell barcode string per line, and may be \
                         gzip-compressed."
        )]
        cells: String,
    },

    #[command(
        about = "Compute scATAC-seq quality control metrics",
        long_about = "Compute TSS enrichment, nucleosome signal, total fragments, fraction of reads in promoters (FRiP), and fraction of fragments in mito, chrX, chrY chromosomes.",
        group(
            ArgGroup::new("annotation")
                .required(true)
                .args(&["gff", "bed"])
    ))]
    Qc {
        #[arg(
            short, 
            long,
            value_name = "FILE",
            help = "Path to the fragment file",
            long_help = "Path to the fragment file. The file should be gzipped and contain \
                         chromosome, start, end, and cell barcode columns tab-separated."
        )]
        fragments: String,

        #[arg(
            short,
            long,
            value_name = "FILE",
            help = "Path to a GFF file containing gene annotations. Supply either a GFF file or BED file, not both.",
        )]
        gff: Option<String>,

        #[arg(
            short,
            long,
            value_name = "FILE",
            help = "Path to a BED file containing TSS positions. Supply either a GFF file or BED file, not both.",
        )]
        bed: Option<String>,

        #[arg(
            short,
            long,
            value_name = "FILE",
            help = "File containing cell barcodes to include",
            long_help = "Optional text file containing cell barcodes to include. \
                         Only these cells will appear in the output. \
                         If not provided, all cell barcodes in the fragment file are included. \
                         The file may be gzip-compressed."
        )]
        cells: Option<String>,

        #[arg(
            short,
            long,
            value_name = "FILE",
            help = "Path to the output file",
            long_help = "Path to the output file. The file will contain the TSS enrichment for each cell barcode, tab-separated."
        )]
        outfile: String,
    },

    #[command(
        name = "generate-manpages",
        about = "Generate man pages",
        hide = true  // Hide from normal help output
    )]
    GenerateManPages {
        #[arg(
            short, 
            long,
            help = "Output directory for man pages"
        )]
        outdir: PathBuf,
    },
}

fn main() -> Result<(), Box<dyn Error>> {
    pretty_env_logger::init_timed();

    let cli = Cli::parse();

    match &cli.command {
        Commands::Matrix { fragments, bed, cells, outdir, threads, pic, group, h5 } => {
            f2m::f2m(fragments, bed, cells, outdir, *threads, *group, *pic, *h5)?
        },
        Commands::Count { fragments, outfile, threshold, ncells } => {
            cellselect::cellselect(fragments, outfile, threshold, ncells)?
        },
        Commands::Filter { fragments, cells } => {
            filter::run(fragments, cells)?
        },
        Commands::GenerateManPages { outdir } => {
            man::generate_manpages(outdir)?;
        },
        Commands::Qc { fragments, gff, bed, cells, outfile } => {
            let (annotation, annotation_is_gff) = match (gff.as_ref(), bed.as_ref()) {
                (Some(gff_path), None) => (gff_path, true),
                (None, Some(bed_path)) => (bed_path, false),
                _ => unreachable!("clap ArgGroup ensures exactly one of gff or bed is present"),
            };
            qc::tss_enrichment(fragments, annotation, annotation_is_gff, cells.as_deref(), outfile)?
        },
    }

    Ok(())
}