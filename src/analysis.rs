use crate::densmap::{index2tuple, DensMap};

#[derive(Debug)]
pub struct GraphData {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

pub fn calc_density_per_radius(densmap: &DensMap) -> GraphData {
    let [dx, dy, _] = densmap.bin_size;
    let [xmin, ymin] = densmap.origin;

    let dr = 0.5 * (dx + dy);
    let rmin = 1.0;
    let rmax = calc_maximum_radius(&densmap);

    let num_values = ((rmax - rmin) / dr) as usize;

    let radius = (0..=num_values)
        .map(|n| rmin + dr * n as f64)
        .collect::<Vec<_>>();
    let mut histogram = vec![0.0; radius.len()];

    // Adjust the center coordinates to be relative to the bins, instead of adjusting
    // the bin coordinates. Those coordinates would be shifted every iteration, which
    // is unnecessary since the relative distance to the center is all we are interested in.
    let x0 = densmap.center[0] - xmin;
    let y0 = densmap.center[1] - ymin;

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

    // Scale the histogram by the circumference at the radius to get the density per unit length.
    use std::f64::consts::PI;
    let scaled_histogram = histogram
        .iter()
        .zip(radius.iter())
        .map(|(v, r)| v / (2.0 * PI * r))
        .collect();

    GraphData {
        x: radius,
        y: scaled_histogram,
    }
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
