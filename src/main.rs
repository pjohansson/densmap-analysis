mod analysis;
mod average;
mod densmap;
mod graphdata;

use structopt::StructOpt;

use std::{
    io,
    path::{Path, PathBuf},
};

use self::{
    analysis::{
        radial_density::{get_radial_density_distribution, get_radius_from_distribution},
        sample_interface::sample_interface,
    },
    average::smoothen_data_of_bins_within_radius,
    densmap::{read_densmap, write_densmap},
    graphdata::{write_xvg, XYData},
};

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

    if let Ok(radius) = get_radius_from_distribution(radial_density) {
        let contact_line = sample_interface(&smoothed_densmap, radius);
        let new_angles = (0..2160).map(|n| n as f64 * 360.0 / 2160.0).collect::<Vec<_>>();
        let resampled = contact_line.resample(&new_angles);
        write_xvg(&resampled);
    }

    let out = Path::new("smooth.dat");
    write_densmap(&out, &smoothed_densmap, time)?;

    Ok(())
}
