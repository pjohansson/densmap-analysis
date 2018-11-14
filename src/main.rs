mod analysis;
mod average;
mod densmap;

use structopt::StructOpt;

use std::io;
use std::path::{Path, PathBuf};

use crate::analysis::calc_density_per_radius;
use crate::average::smoothen_data_of_bins_within_radius;
use crate::densmap::{read_densmap, write_densmap};

#[derive(Debug, StructOpt)]
struct Args {
    #[structopt(parse(from_os_str))]
    filename: PathBuf,
}

fn main() -> Result<(), io::Error> {
    let args = Args::from_args();

    let (densmap, time) = read_densmap(&args.filename)?;

    let smoothed_densmap = smoothen_data_of_bins_within_radius(densmap, 0.5);
    let radial_density = calc_density_per_radius(&smoothed_densmap);

    radial_density
        .x
        .iter()
        .zip(density.y.iter())
        .for_each(|(x, y)| println!("{:12.5} {:12.5}", x, y));

    let out = Path::new("smooth.dat");
    write_densmap(&out, &smoothed_densmap, time)?;

    Ok(())
}
