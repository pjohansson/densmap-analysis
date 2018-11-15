use crate::{
    densmap::{coord2index, DensMap},
    graphdata::Graph,
};

/// Sample the contact line interface per angle from the droplet center point.
pub fn sample_interface(densmap: &DensMap, base_radius: f64) -> Graph {
    // For now, use an angular precision that gives a 0.1nm resolution at the droplet radius.
    let num_values = (2.0 * std::f64::consts::PI * base_radius / 0.1).ceil() as usize;

    // Create the angles array in degrees.
    let da = 360.0 / num_values as f64;
    let angles = (0..num_values).map(|n| da * n as f64).collect::<Vec<_>>();

    let cutoff = get_density_cutoff(&densmap);

    let radius = angles
        .iter()
        .map(|&a| sample_interface_at_angle(&densmap, a, base_radius, cutoff))
        .collect();

    Graph::Polar { angles, radius }
}

#[derive(Clone, Copy, Debug)]
/// Whether we are increasing or decreasing the radius from the initial guess.
enum Direction {
    Increasing,
    Decreasing,
}

/// Find the radius where the density crosses a cutoff for an input angle in degrees.
fn sample_interface_at_angle(densmap: &DensMap, angle: f64, base_radius: f64, cutoff: f64) -> f64 {
    // TODO: Consider sampling the interface using several bins within range of the radius instead
    // of a single bin. Supersampling!

    // Work in the density map relative coordinate space by adjusting the center coordinates,
    // which lie in system absolute space.
    let x0 = densmap.center[0] - densmap.origin[0];
    let y0 = densmap.center[1] - densmap.origin[1];

    // Change the radius by the largest possible amount that cannot skip any bins:
    // The smallest length of a bin.
    let [dx_bin, dy_bin, _] = densmap.bin_size;
    let dr_abs = dx_bin.min(dy_bin);

    let (mut radius, direction) =
        get_initial_radius_and_direction(&densmap, base_radius, angle, cutoff, dr_abs);

    let dr = match direction {
        Direction::Increasing => dr_abs,
        Direction::Decreasing => -dr_abs,
    };

    let (dy, dx) = angle.to_radians().sin_cos();

    while radius >= dr {
        radius += dr;
        let x = x0 + radius * dx;
        let y = y0 + radius * dy;

        match coord2index(x, y, densmap.bin_size, densmap.shape) {
            Some(i) => {
                let is_filled = densmap.data[i] >= cutoff;

                match (direction, is_filled) {
                    // If we are increasing the radius and find an empty bin, use the last radius.
                    (Direction::Increasing, false) => {
                        radius -= dr;
                        break;
                    }
                    // If we are decreasing the radius and find the first full bin, just break.
                    (Direction::Decreasing, true) => {
                        break;
                    }
                    // Otherwise continue the search until we get a match.
                    _ => (),
                }
            }
            // We are outside the system: use the final good value as a radius.
            None => {
                radius -= dr;
                break;
            }
        }
    }

    radius
}

/// Get the density cutoff for contact line determination as half of the maximum density value.
fn get_density_cutoff(densmap: &DensMap) -> f64 {
    0.5 * densmap.data.iter().fold(0.0, |acc: f64, &v| acc.max(v))
}

/// From the initial base radius, get the first valid radius and search direction.
///
/// The first valid radius may be smaller than the input base radius, since it may lie outside
/// of the system.
fn get_initial_radius_and_direction(
    densmap: &DensMap,
    base_radius: f64,
    angle: f64,
    cutoff: f64,
    dr: f64,
) -> (f64, Direction) {
    let x0 = densmap.center[0] - densmap.origin[0];
    let y0 = densmap.center[1] - densmap.origin[1];
    let (dy, dx) = angle.to_radians().sin_cos();

    let mut radius = base_radius;

    let direction = loop {
        let x = x0 + radius * dx;
        let y = y0 + radius * dy;

        if let Some(index) = coord2index(x, y, densmap.bin_size, densmap.shape) {
            if densmap.data[index] >= cutoff {
                break Direction::Increasing;
            } else {
                break Direction::Decreasing;
            }
        }

        // We reach this point if and only if the initial guess is outside the system.
        radius -= dr;
    };

    (radius, direction)
}
