#[cfg(not(target_os = "windows"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: std::alloc::System = std::alloc::System;

use clap::{Parser, Subcommand};
use std::error::Error;

mod f2m;
mod cellselect;
mod filter;

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
    /// Create a feature x cell matrix from a fragment file
    Matrix {
        /// Path to the fragment file
        #[arg(short, long)]
        fragments: String,

        /// BED file containing non-overlapping genomic regions to quantify
        #[arg(short, long)]
        bed: String,

        /// File containing cell barcodes to include
        #[arg(short, long)]
        cells: String,

        /// Output directory name
        #[arg(short, long, help = "Output directory name", long_help = "Output directory name. Directory will be created if it does not exist. The output directory will contain matrix.mtx.gz, features.tsv, barcodes.tsv")]
        outdir: String,

        /// Number of compression threads to use
        #[arg(short, long, default_value = "4", value_parser = clap::value_parser!(usize))]
        threads: usize,

        /// Use paired insertion counting
        #[arg(long, action = clap::ArgAction::SetTrue)]
        pic: bool,

        /// Group peaks by variable in fourth BED column
        #[arg(long, action = clap::ArgAction::SetTrue)]
        group: bool,
    },

    /// Count number of fragments per cell barcode
    Count {
        /// Path to the fragment file
        #[arg(short, long, value_name = "FILE")]
        fragments: String,

        /// Name of output file
        #[arg(short, long, value_name = "FILE", help = "Name of output file", long_help = "Name of output file. The file will contain each cell barcode and its total fragment count, tab-separated.")]
        outfile: String,

        /// Minimum number of fragments for a cell to be included
        #[arg(short, long, value_name = "NUMBER", help = "Minimum number of fragments for a cell to be included", long_help = "Sets the minimum number of fragments a cell must have to be included in the output. Cells with fewer fragments than this threshold will be filtered out. Cannot be used with --ncells", conflicts_with = "ncells", value_parser = clap::value_parser!(usize))]
        threshold: Option<usize>,

        /// Number of top cells to select
        #[arg(short, long, value_name = "NUMBER", help = "Number of top cells to select", long_help = "Select this many cells with the highest fragment counts. Cannot be used together with --threshold.", conflicts_with = "threshold", value_parser = clap::value_parser!(usize))]
        ncells: Option<usize>,
    },

    /// Subset a fragment file to include only specified cell barcodes
    Filter {
        /// Path to the fragment file
        #[arg(short, long)]
        fragments: String,

        /// File containing cell barcodes to include
        #[arg(short, long)]
        cells: String,
    },
}

fn main() -> Result<(), Box<dyn Error>> {
    pretty_env_logger::init_timed();

    let cli = Cli::parse();

    match &cli.command {
        Commands::Matrix { fragments, bed, cells, outdir, threads, pic, group } => {
            f2m::f2m(fragments, bed, cells, outdir, *threads, *group, *pic)?
        },
        Commands::Count { fragments, outfile, threshold, ncells } => {
            cellselect::cellselect(fragments, outfile, threshold, ncells)?
        },
        Commands::Filter { fragments, cells } => {
            filter::run(fragments, cells)?
        },
    }

    Ok(())
}