/// Write input data as an .xvg formatted file.
pub fn write_xvg<T: XYData>(data: &T) {
    data.x()
        .iter()
        .zip(data.y().iter())
        .for_each(|(x, y)| println!("{:12.5} {:12.5}", x, y));
}

/// Trait for data that has values corresponding to x and y axes.
pub trait XYData {
    fn x(&self) -> &[f64];
    fn y(&self) -> &[f64];
}

#[derive(Clone, Debug)]
/// Histogram data with bin center points and values.
pub struct Histogram {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

impl XYData for Histogram {
    fn x(&self) -> &[f64] {
        &self.x
    }

    fn y(&self) -> &[f64] {
        &self.y
    }
}

#[derive(Clone, Debug)]
/// Two dimensional graph types.
pub enum Graph {
    Carthesian { x: Vec<f64>, y: Vec<f64> },
    /// Polar coordinates have angles in degrees.
    Polar { angles: Vec<f64>, radius: Vec<f64> },
}

/// For carthesian coordinates, x and y values correspond directly to the variables.
/// For polar coordinates the angles become the x values and the radius the y values.
impl XYData for Graph {
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
