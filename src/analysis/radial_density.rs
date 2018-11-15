use crate::{
    densmap::{index2tuple, DensMap},
    graphdata::Histogram,
};

/// Compute the radial density distribution function p(r) for the density map, using
/// the center point of the droplet as the origin.
///
/// The distribution is scaled to have units of mass / nm of the circumference at the radius.
pub fn get_radial_density_distribution(densmap: &DensMap) -> Histogram {
    let (rmin, dr, radius) = get_radius_values_for_histogram(&densmap);
    let histogram = get_radial_mass_sum_of_densmap(&densmap, rmin, dr, radius.len());
    let scaled_histogram = scale_histogram_to_per_unit_length(&histogram, &radius);

    Histogram {
        x: radius,
        y: scaled_histogram,
    }
}

/// Get the droplet radius from the radial density distribution by taking the midpoint
/// between the 10th and 90th percentile values.
///
/// Empty histogram bins are cut from the distribution before taking the percentiles.
/// The criteria for being empty is to have a value lower than 1% of the maximum.
///
/// # Error
/// If the radial density distribution is empty the percentiles cannot be calculated.
///
/// # Notes
/// Makes the assumption that the array of radius values match the density values
/// after invalid values (inf and NaN) have been removed from it. These shouldn't be
/// there in the first place unless something has gone *very* wrong when calculating
/// the density distribution in the first place.
pub fn get_radius_from_distribution(radial_density: Histogram) -> Result<f64, String> {
    // Ensure that we only have good numbers, no NaN or infs.
    let density = radial_density
        .y
        .into_iter()
        .filter(|v| v.is_finite())
        .collect::<Vec<_>>();

    let density_nonzero = cut_bins_below_percentage_of_max(&density, 1.0);
    let (lower_density, upper_density) = get_percentile_values(&density_nonzero, 10.0, 90.0)?;
    let mid_density = 0.5 * (lower_density + upper_density);

    // Find the vector index where the density value is reached by sweeping the histogram
    // in reverse. We do it in reverse because density fluctuations in the bulk may give
    // us an early hit in some edge cases, while this is less likely to happen if we look
    // move from the zero-density values to the bulk.
    //
    // We can unwrap since `get_percentile_values` would have returned an error if we
    // have no values in the array, which is the only way we would not be able to find
    // a mid point between two values.
    let i = density.iter().rposition(|&v| v >= mid_density).unwrap();

    Ok(radial_density.x[i])
}

/// Return bins which have values larger than or equal to a cutoff, determined by the maximum.
///
/// # Notes
/// Assumes that all values are positive.
fn cut_bins_below_percentage_of_max(values: &[f64], perc: f64) -> Vec<f64> {
    let max = values.iter().fold(0.0, |acc: f64, &v| acc.max(v));
    let cutoff = 0.01 * perc * max;

    values
        .into_iter()
        .cloned()
        .filter(|&v| v >= cutoff)
        .collect()
}

/// # Notes
/// Assumes that all input values are valid for a comparison, eg. that floats are not NaN or inf.
fn get_percentile_values<T: Copy + PartialOrd>(
    values: &[T],
    lower: f64,
    upper: f64,
) -> Result<(T, T), String> {
    if values.is_empty() {
        return Err(String::from(
            "cannot compute percentile values from an empty array",
        ));
    }
    if lower < 0.0 || lower > 100.0 {
        return Err(format!(
            "lower percentile value must be between 0 and 100, was {}",
            lower
        ));
    }
    if upper < 0.0 || upper > 100.0 {
        return Err(format!(
            "upper percentile value must be between 0 and 100, was {}",
            upper
        ));
    }

    let mut sorted_values = values.to_vec();
    sorted_values.sort_by(|a, b| {
        a.partial_cmp(b)
            .expect("could not compare two values when calculating percentile values")
    });

    let ilower = (0.01 * lower * sorted_values.len() as f64).round() as usize;
    let iupper = (0.01 * upper * sorted_values.len() as f64).round() as usize;

    Ok((sorted_values[ilower], sorted_values[iupper]))
}

fn get_radius_values_for_histogram(densmap: &DensMap) -> (f64, f64, Vec<f64>) {
    let [dx, dy, _] = densmap.bin_size;

    let dr = 0.5 * (dx + dy);
    let rmin = 1.0;
    let rmax = calc_maximum_radius(&densmap);

    let num_values = ((rmax - rmin) / dr) as usize;
    let values = (0..=num_values).map(|n| rmin + dr * n as f64).collect();

    (rmin, dr, values)
}

fn get_radial_mass_sum_of_densmap(
    densmap: &DensMap,
    rmin: f64,
    dr: f64,
    num_bins: usize,
) -> Vec<f64> {
    // Adjust the center coordinates to be relative to the bins, instead of adjusting
    // the bin coordinates. Those coordinates would be shifted every iteration, which
    // is unnecessary since the relative distance to the center is all we are interested in.
    let [dx, dy, _] = densmap.bin_size;
    let [xmin, ymin] = densmap.origin;
    let x0 = densmap.center[0] - xmin;
    let y0 = densmap.center[1] - ymin;

    let mut histogram = vec![0.0; num_bins];

    densmap
        .data
        .iter()
        .enumerate()
        // Get the 2D position of the bin from its 1D index.
        .map(|(i, v)| (index2tuple(i, densmap.shape).unwrap(), v))
        // Convert to system coordinates.
        .map(|((ix, iy), v)| ((dx * ix as f64, dy * iy as f64), v))
        // Calculate distance to center.
        .map(|((x, y), v)| (((x0 - x).powi(2) + (y0 - y).powi(2)).sqrt(), v))
        // Exclude points that are too close to the center, they're noisy.
        .filter(|(r, _)| r >= &rmin)
        // Add the value to the histogram at the radius.
        .for_each(|(r, v)| {
            let n = ((r - rmin) / dr) as usize;
            histogram[n] += v;
        });

    histogram
}

fn scale_histogram_to_per_unit_length(histogram: &[f64], radius: &[f64]) -> Vec<f64> {
    histogram
        .iter()
        .zip(radius.iter())
        .map(|(v, r)| v / (2.0 * std::f64::consts::PI * r))
        .collect()
}

/// Calculate the distance from the fitted droplet to the furthest away bin in the system.
fn calc_maximum_radius(densmap: &DensMap) -> f64 {
    let [dx, dy, _] = densmap.bin_size;
    let [nx, ny] = densmap.shape;

    let [xmin, ymin] = densmap.origin;
    let xmax = xmin + dx * nx as f64;
    let ymax = ymin + dy * ny as f64;

    let [x0, y0] = densmap.center;

    let mut rmax2 = (x0 - xmin).powi(2) + (y0 - ymin).powi(2);
    rmax2 = rmax2.max((x0 - xmin).powi(2) + (y0 - ymax).powi(2));
    rmax2 = rmax2.max((x0 - xmax).powi(2) + (y0 - ymin).powi(2));
    rmax2 = rmax2.max((x0 - xmax).powi(2) + (y0 - ymax).powi(2));

    rmax2.sqrt()
}

#[test]
fn test_cut_bins_below_50_percent_of_max() {
    assert_eq!(
        vec![4.1, 5.0, 4.5, 8.0],
        cut_bins_below_percentage_of_max(&vec![0.0, 4.1, 3.9, 5.0, 4.5, 8.0, 3.5], 50.0)
    );
}

#[test]
fn test_get_20th_and_80th_percentile_values_of_reversed_array() {
    // In this array, 20% of the values lie below (or at) 1.0 and 80% below (or at) 2.0.
    let values = vec![3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0, 0.0];
    assert_eq!(Ok((1.0, 2.0)), get_percentile_values(&values, 20.0, 80.0));
}

#[test]
fn test_getting_percentiles_from_empty_array_returns_error() {
    assert!(get_percentile_values(&Vec::<f64>::new(), 0.0, 0.0).is_err());
}

#[test]
fn test_getting_invalid_percentiles_returns_error() {
    let values = vec![3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0, 0.0];
    assert!(get_percentile_values(&values, -1.0, 0.0).is_err());
    assert!(get_percentile_values(&values, 0.0, -1.0).is_err());
    assert!(get_percentile_values(&values, 101.0, 0.0).is_err());
    assert!(get_percentile_values(&values, 0.0, 101.0).is_err());
}
