use clap::CommandFactory;
use clap_mangen::Man;
use std::fs::File;
use std::io::Result;
use std::path::Path;

use crate::Cli;

pub fn generate_manpages<P: AsRef<Path>>(outdir: P) -> Result<()> {
    let outdir = outdir.as_ref();
    std::fs::create_dir_all(outdir)?;

    // Generate man page for main command
    let cmd = Cli::command();
    let mut out = File::create(outdir.join("fragtk.1"))?;
    Man::new(cmd.clone()).render(&mut out)?;

    // Generate man pages for subcommands
    for subcmd in cmd.get_subcommands() {
        let mut out = File::create(outdir.join(format!("fragtk-{}.1", subcmd.get_name())))?;
        Man::new(subcmd.clone()).render(&mut out)?;
    }

    Ok(())
}
