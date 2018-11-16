// #![feature(impl_trait_in_bindings)]

mod analysis;
mod average;
mod densmap;
mod graphdata;

use regex::Regex;
use structopt::StructOpt;
use walkdir::{DirEntry, WalkDir};

use std::{
    ffi::{OsStr, OsString},
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
#[structopt(raw(setting = "structopt::clap::AppSettings::DeriveDisplayOrder"))]
struct Args {
    #[structopt(parse(from_os_str), conflicts_with = "base", required_unless = "base")]
    /// List of density map files to analyze
    filenames: Vec<PathBuf>,
    #[structopt(long = "base", parse(from_os_str))]
    /// Base file name for density maps
    base: Option<PathBuf>,
    #[structopt(long = "ext", default_value = "dat", parse(from_os_str))]
    /// Extension for file names
    ext: OsString,
    #[structopt(
        long = "time_sig",
        default_value = r"([0-9]{5}\.[0-9]{3})ps",
        value_name = "regex"
    )]
    /// Regular expression for time signature in file names
    time_regex: String,
    #[structopt(short = "b", long = "begin", requires = "base", value_name = "t0")]
    /// Only include times for which t >= t0
    begin: Option<f64>,
    #[structopt(short = "e", long = "end", requires = "base", value_name = "t1")]
    /// Only include times for which t <= t1
    end: Option<f64>,
    #[structopt(short = "d", long = "dt", requires = "base", value_name = "dt")]
    /// Only include times for which t % dt = 0
    dt: Option<f64>,
}

fn main() -> Result<(), io::Error> {
    let args = Args::from_args();

    let filenames = match args.base {
        Some(base) => construct_file_list(
            &base,
            &args.time_regex,
            &args.ext,
            args.begin,
            args.end,
            args.dt,
        ),
        None => args.filenames,
    };

    eprintln!("{:?}", filenames);

    let (densmap, time) = read_densmap(&filenames[0])?;
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

/// Use the given base filename path along with the regular expression for time signatures
/// and the file extension to identify file names and return them as an array.
///
/// Filter out file names with times which do not lie within the (optional) interval.
fn construct_file_list(
    base_path: &Path,
    time_regex: &str,
    ext: &OsStr,
    begin: Option<f64>,
    end: Option<f64>,
    dt: Option<f64>,
) -> Vec<PathBuf> {
    let dir = base_path.parent().unwrap_or(&Path::new("."));
    let base = base_path
        .file_name()
        .unwrap_or(&OsStr::new(""))
        .to_str()
        .unwrap();

    let regex_string = format!(r"{}{}\.{}", base, time_regex, ext.to_str().unwrap());
    let re = Regex::new(&regex_string).unwrap();

    WalkDir::new(dir)
        .min_depth(1)
        .max_depth(1)
        .sort_by(|a, b| a.file_name().cmp(b.file_name()))
        .into_iter()
        .filter_entry(|entry| {
            let file_name = entry.file_name().to_str().unwrap();

            match re.captures(&file_name) {
                Some(captures) => {
                    let time = captures.get(1).unwrap().as_str().parse::<f64>().unwrap();
                    begin.map(|b| time >= b).unwrap_or(true)
                        && end.map(|e| time <= e).unwrap_or(true)
                        && dt.map(|d| time % d == 0.0).unwrap_or(true)
                }
                None => false,
            }
        })
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.into_path())
        .collect()
}
