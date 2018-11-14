mod analysis;
mod average;
mod densmap;

use structopt::StructOpt;

use std::io;
use std::path::{Path, PathBuf};

use crate::analysis::{get_radial_density_distribution, get_radius_from_distribution};
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
    let radial_density = get_radial_density_distribution(&smoothed_densmap);

    radial_density
        .x
        .iter()
        .zip(radial_density.y.iter())
        .for_each(|(x, y)| println!("{:12.5} {:12.5}", x, y));
    eprintln!("{:?}", get_radius_from_distribution(radial_density));

    let out = Path::new("smooth.dat");
    write_densmap(&out, &smoothed_densmap, time)?;

    Ok(())
}
