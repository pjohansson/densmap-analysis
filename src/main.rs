mod analysis;
mod average;
mod densmap;
mod graphdata;

use pbr::ProgressBar;
use regex::Regex;
use structopt::StructOpt;
use walkdir::WalkDir;

use std::{
    ffi::{OsStr, OsString},
    io,
    path::{Path, PathBuf},
};

use self::{
    analysis::{
        autocorrelation::calc_autocorrelation,
        radial_density::{get_radial_density_distribution, get_radius_from_distribution},
        sample_interface::sample_interface,
    },
    average::smoothen_data_of_bins_within_radius,
    densmap::{read_densmap, write_densmap},
    graphdata::{write_xvg, Graph, Histogram, XYData},
};

#[derive(Debug, StructOpt)]
#[structopt(raw(
    setting = "structopt::clap::AppSettings::ColoredHelp",
    setting = "structopt::clap::AppSettings::DeriveDisplayOrder"
))]
struct Args {
    #[structopt(parse(from_os_str), conflicts_with = "base", required_unless = "base")]
    /// List of density map files to analyze
    filenames: Vec<PathBuf>,

    #[structopt(long = "base", value_name = "path", parse(from_os_str))]
    /// Base file name for density maps
    base: Option<PathBuf>,

    #[structopt(short = "d", long = "densmap", value_name = "path", parse(from_os_str))]
    /// Base output file name for smoothed density maps
    smooth: Option<PathBuf>,
    /// Base output file name for contact line angular distributions
    #[structopt(
        long = "contact_line",
        value_name = "path",
        hidden_short_help = true,
        parse(from_os_str)
    )]
    contact_line: Option<PathBuf>,

    #[structopt(
        long = "interface",
        value_name = "path",
        hidden_short_help = true,
        parse(from_os_str)
    )]
    /// Base output file name for interface graphs
    interface: Option<PathBuf>,

    #[structopt(
        short = "r",
        long = "radius",
        default_value = "radius.xvg",
        value_name = "path",
        parse(from_os_str)
    )]
    /// Output file name for droplet radius time series
    radius: PathBuf,

    #[structopt(
        long = "rdd",
        value_name = "path",
        hidden_short_help = true,
        parse(from_os_str)
    )]
    /// Base output file name for radial density distributions
    radial_density: Option<PathBuf>,

    #[structopt(
        long = "ac",
        value_name = "path",
        hidden_short_help = true,
        parse(from_os_str)
    )]
    /// Output file name for contact line autocorrelation
    autocorrelation: Option<PathBuf>,

    #[structopt(long = "ext", default_value = "dat", parse(from_os_str))]
    /// Extension for density map file names
    ext: OsString,
    #[structopt(
        long = "time_sig",
        long_help = "Regular expression for time signature in file names. The expression must include a capture group around the time value.",
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
    #[structopt(long = "dt", requires = "base", value_name = "dt")]
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

    let mut radius_time_series = Vec::with_capacity(filenames.len());
    let mut times = Vec::with_capacity(filenames.len());

    // To calculate the autocorrelation of contact line fluctuations we need to save
    // the contact line for every time step.
    let mut contact_line_per_time = Vec::with_capacity(filenames.len());

    let mut pb = ProgressBar::new(filenames.len() as u64);
    pb.format("[=> ]");

    for (i, filename) in filenames.into_iter().enumerate() {
        pb.message(&format!("Processing '{}' ", &filename.to_str().unwrap()));
        pb.inc();

        let (densmap, time) = read_densmap(&filename)?;

        let dir = filename.parent().unwrap();
        let time_signature = read_time_signature_or_default(&filename, &args.time_regex, i);

        let smoothed_densmap = smoothen_data_of_bins_within_radius(densmap, 0.5);
        if let Some(base) = &args.smooth {
            let path = construct_file_name(&base, &time_signature, &args.ext, &dir);
            write_densmap(&path, &smoothed_densmap, time)?;
        }

        let radial_density = get_radial_density_distribution(&smoothed_densmap);
        if let Some(base) = &args.radial_density {
            let path = construct_file_name(&base, &time_signature, &args.ext, &dir);
            write_xvg(&path, &radial_density)?;
        }

        if let Ok(radius) = get_radius_from_distribution(radial_density) {
            radius_time_series.push(radius);
            times.push(time);

            let contact_line = sample_interface(&smoothed_densmap, radius);
            let interface = contact_line.to_carthesian();
            if let Some(base) = &args.interface {
                let path = construct_file_name(&base, &time_signature, &args.ext, &dir);
                write_xvg(&path, &interface)?;
            }

            let relative_contact_line = Graph::Polar {
                angles: contact_line.x().to_vec(),
                radius: contact_line.y().iter().map(|r| r - radius).collect(),
            };

            if let Some(base) = &args.contact_line {
                let path = construct_file_name(&base, &time_signature, &args.ext, &dir);
                write_xvg(&path, &relative_contact_line)?;
            }

            contact_line_per_time.push(relative_contact_line);
        }
    }

    pb.finish_print("Processed all density maps.");
    eprintln!("");

    if let Some(filename) = args.autocorrelation {
        let mut pb = ProgressBar::new(contact_line_per_time.len() as u64);
        pb.message("Calculating autocorrelation of contact line ");

        // To calculate the autocorrelation for the contact line we resample the data onto
        // a common set of angles. We use the largest number of sample points as the base.
        let resample_xvals = contact_line_per_time
            .iter()
            .max_by(|&a, &b| a.x().len().cmp(&b.x().len()))
            .unwrap()
            .x();
        let resampled_contact_lines = contact_line_per_time
            .iter()
            .map(|contact_line| contact_line.resample(&resample_xvals))
            .collect::<Vec<_>>();

        let autocorrelation_yvals = calc_autocorrelation(&resampled_contact_lines);
        let autocorrelation = Histogram {
            x: times.clone(),
            y: autocorrelation_yvals,
        };

        write_xvg(&filename, &autocorrelation)?;
        pb.finish_print("Finished autocorrelation calculation.");
    }

    let radius_per_time = Graph::Carthesian {
        x: times,
        y: radius_time_series,
    };
    write_xvg(&args.radius, &radius_per_time)?;

    Ok(())
}

fn construct_file_name(base: &Path, time_sig: &str, ext: &OsStr, dir: &Path) -> PathBuf {
    let file_name =
        PathBuf::from(base.to_str().unwrap().to_string() + time_sig + "." + ext.to_str().unwrap());
    dir.join(&file_name)
}

fn read_time_signature_or_default(path: &Path, time_regex: &str, index: usize) -> String {
    // For now we can recompile the regular expression for every time step, even though
    // it does not change. The compilation time is marginal at best compared to all the analysis.
    // It may be poor practice, but eh.
    let re = Regex::new(&time_regex).unwrap();

    match re.captures(&path.to_str().unwrap()) {
        Some(capture) => String::from(capture.get(0).unwrap().as_str()),
        None => format!("{:05}", index + 1),
    }
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
