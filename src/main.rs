#[cfg(not(target_os = "windows"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(target_os = "windows")]
#[global_allocator]
static GLOBAL: std::alloc::System = std::alloc::System;

use clap::{Command, Arg, ArgAction};
use std::error::Error;

mod f2m;
mod cellselect;
mod filter;


fn main() -> Result<(), Box<dyn Error>> {

    let matches = Command::new("fragtk")
        .version("1.1.0")
        .author("Tim Stuart")
        .about("Fragment file processing tools")
        .arg_required_else_help(true)
        .subcommand(
            Command::new("matrix")
                .about("Create a feature x cell matrix from a fragment file")
                .arg(
                    Arg::new("fragments")
                        .short('f')
                        .long("fragments")
                        .help("Path to the fragment file")
                        .required(true),
                )
                .arg(
                    Arg::new("bed")
                        .short('b')
                        .long("bed")
                        .help("BED file containing non-overlapping genomic regions to quantify")
                        .required(true),
                )
                .arg(
                    Arg::new("cells")
                        .short('c')
                        .long("cells")
                        .help("File containing cell barcodes to include")
                        .required(true),
                )
                .arg(
                    Arg::new("outdir")
                        .short('o')
                        .long("outdir")
                        .help("Output directory name")
                        .long_help("Output directory name. Directory will be created if it does not exist. \
                               The output directory will contain matrix.mtx.gz, features.tsv, barcodes.tsv")
                        .required(true),
                )
                .arg(
                    Arg::new("threads")
                        .short('t')
                        .long("threads")
                        .help("Number of compression threads to use")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("4")
                        .required(false),
                )
                .arg(
                    Arg::new("group")
                        .long("group")
                        .help("Group peaks by variable in fourth BED column")
                        .action(ArgAction::SetTrue),
                )
        )
        .subcommand(
            Command::new("count")
            .about("Count number of fragments per cell barcode")
            .arg(
                Arg::new("fragments")
                    .short('f')
                    .long("fragments")
                    .value_name("FILE")
                    .help("Path to the fragment file")
                    .required(true),
            )
            .arg(
                Arg::new("outfile")
                    .short('o')
                    .long("outfile")
                    .value_name("FILE")
                    .help("Name of output file")
                    .long_help(
                        "Name of output file. The file will contain each cell barcode \
                        and its total fragment count, tab-separated."
                    )
                    .required(true),
            )
            .arg(
                Arg::new("threshold")
                    .short('t')
                    .long("threshold")
                    .value_name("NUMBER")
                    .help("Minimum number of fragments for a cell to be included")
                    .long_help(
                        "Sets the minimum number of fragments a cell must have to be included in the output. \
                        Cells with fewer fragments than this threshold will be filtered out. Cannot be used with --ncells"
                    )
                    .conflicts_with("ncells")
                    .required(false),
            )
            .arg(
                Arg::new("ncells")
                    .short('n')
                    .long("ncells")
                    .value_name("NUMBER")
                    .help("Number of top cells to select")
                    .long_help(
                        "Select this many cells with the highest fragment counts. \
                        Cannot be used together with --threshold."
                    )
                    .conflicts_with("threshold")
                    .required(false),
            )
        )
        .subcommand(
            Command::new("filter")
                .about(
                    "Subset a fragment file to include only specified cell barcodes. \
                    Output is uncompressed data written to stdout"
                )
                .arg(
                    Arg::new("fragments")
                        .short('f')
                        .long("fragments")
                        .help("Path to the fragment file")
                        .required(true),
                )
                .arg(
                    Arg::new("cells")
                        .short('c')
                        .long("cells")
                        .help("File containing cell barcodes to include")
                        .required(true),
                )
        )
        .get_matches();

    pretty_env_logger::init_timed();

    match matches.subcommand() {
        Some(("matrix", sub_matches)) => f2m::f2m(sub_matches)?,
        Some(("count", sub_matches)) => cellselect::cellselect(sub_matches)?,
        Some(("filter", sub_matches)) => filter::run(sub_matches)?,
        _ => {

        }
    }

    Ok(())
}