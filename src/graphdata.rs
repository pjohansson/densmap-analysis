/// Write input data as an .xvg formatted file.
pub fn write_xvg<T: XYData>(data: &T) {
    data.x()
        .iter()
        .zip(data.y().iter())
        .for_each(|(x, y)| println!("{:12.5} {:12.5}", x, y));
}

/// Trait for data that has values corresponding to x and y axes.
pub trait XYData: PartialEq {
    /// Resample the data onto a new set of x values.
    ///
    /// # Notes
    /// Assumes that the current data is sorted along x.
    fn resample(&self, xs: &[f64]) -> Self;
    fn x(&self) -> &[f64];
    fn y(&self) -> &[f64];
}

#[derive(Clone, Debug, PartialEq)]
/// Histogram data with bin center points and values.
pub struct Histogram {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

impl XYData for Histogram {
    fn resample(&self, xs: &[f64]) -> Self {
        Histogram {
            x: xs.to_vec(),
            y: interpolate_data(&self.x, &self.y, &xs),
        }
    }

    fn x(&self) -> &[f64] {
        &self.x
    }

    fn y(&self) -> &[f64] {
        &self.y
    }
}

#[derive(Clone, Debug, PartialEq)]
/// Two dimensional graph types.
pub enum Graph {
    Carthesian {
        x: Vec<f64>,
        y: Vec<f64>,
    },
    /// Polar coordinates have angles in degrees.
    Polar {
        angles: Vec<f64>,
        radius: Vec<f64>,
    },
}

/// For carthesian coordinates, x and y values correspond directly to the variables.
/// For polar coordinates the angles become the x values and the radius the y values.
impl XYData for Graph {
    fn resample(&self, xs: &[f64]) -> Self {
        match self {
            Graph::Carthesian { x, y } => Graph::Carthesian {
                x: xs.to_vec(),
                y: interpolate_data(&x, &y, &xs),
            },
            Graph::Polar { angles, radius } => Graph::Polar {
                angles: xs.to_vec(),
                radius: interpolate_data(&angles, &radius, &xs),
            },
        }
    }

    fn x(&self) -> &[f64] {
        match self {
            Graph::Carthesian { x, y: _ } => &x,
            Graph::Polar { angles, radius: _ } => &angles,
        }
    }

    fn y(&self) -> &[f64] {
        match self {
            Graph::Carthesian { x: _, y } => &y,
            Graph::Polar { angles: _, radius } => &radius,
        }
    }
}

impl Graph {
    pub fn to_carthesian(&self) -> Self {
        match self {
            Graph::Polar { angles, radius } => {
                let mut x = Vec::with_capacity(angles.len());
                let mut y = Vec::with_capacity(angles.len());

                angles.iter().zip(radius.iter()).for_each(|(a, r)| {
                    let (dy, dx) = a.to_radians().sin_cos();
                    x.push(r * dx);
                    y.push(r * dy);
                });

                Graph::Carthesian { x, y }
            }
            Graph::Carthesian { x: _, y: _ } => self.clone(),
        }
    }

    pub fn to_polar(&self) -> Self {
        match self {
            Graph::Carthesian { x: xs, y: ys } => {
                let mut angles = Vec::with_capacity(xs.len());
                let mut radius = Vec::with_capacity(xs.len());

                xs.iter().zip(ys.iter()).for_each(|(&x, &y)| {
                    angles.push(y.atan2(x).to_degrees());
                    radius.push((x.powi(2) + y.powi(2)).sqrt());
                });

                Graph::Polar { angles, radius }
            }
            Graph::Polar {
                angles: _,
                radius: _,
            } => self.clone(),
        }
    }
}

/// Resample data from a set of input x values onto another using linear interpolation.
///
/// # Notes
/// Assumes that the x values of the input data and the final x values are sorted.
fn interpolate_data(from_xs: &[f64], ys: &[f64], onto_xs: &[f64]) -> Vec<f64> {
    onto_xs
        .iter()
        .map(|x| (x, from_xs.binary_search_by(|x1| x1.partial_cmp(x).unwrap())))
        .map(|(x, r)| {
            match r {
                Ok(i) => ys[i],
                Err(i) => {
                    let (i0, i1) = if i == 0 {
                        (0, 1)
                    } else if i == from_xs.len() {
                        (i - 2, i - 1)
                    } else {
                        (i - 1, i)
                    };

                    let x0 = from_xs[i0];
                    let x1 = from_xs[i1];
                    let y0 = ys[i0];
                    let y1 = ys[i1];

                    (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0)
                },
            }
        })
        .collect()
}

#[test]
fn test_interpolate_onto_midpoint_values() {
    let from_xs = vec![0.0, 1.0, 2.0, 3.0];
    let ys =      vec![5.0, 1.0, 3.0, 0.0];

    let onto_xs = vec![0.5, 1.5, 2.5];

    assert_eq!(
        vec![3.0, 2.0, 1.5],
        interpolate_data(&from_xs, &ys, &onto_xs)
    );
}

#[test]
fn test_interpolating_values_outside_of_initial_range_is_linear() {
    let from_xs = vec![0.0, 1.0, 2.0, 3.0];
    let ys =      vec![5.0, 1.0, 3.0, 0.0];

    let onto_xs = vec![-0.5, 4.0];

    assert_eq!(
        vec![7.0, -3.0],
        interpolate_data(&from_xs, &ys, &onto_xs)
    );
}
