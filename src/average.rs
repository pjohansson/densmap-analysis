use rayon::prelude::*;
use crate::densmap::{index2tuple, tuple2index, DensMap, Shape, Vec3};

pub fn smoothen_data_of_bins_within_radius(densmap: DensMap, radius: f64) -> DensMap {
    let neighbours = get_system_bin_neighbours(radius, densmap.bin_size, densmap.shape);

    DensMap {
        data: get_averaged_system(densmap.data, neighbours),
        ..densmap
    }
}

/// In order and for all bins in the system, get a list of its neighbours and return them all.
pub fn get_system_bin_neighbours(radius: f64, bin_size: Vec3, shape: Shape) -> Vec<Vec<usize>> {
    let sieve = get_averaging_bin_sieve(radius, bin_size);
    let [nx, ny] = shape;
    (0..(nx * ny) as usize)
        .into_par_iter()
        .map(|i| get_bin_neighbours(i, shape, &sieve))
        .collect()
}

/// For every bin in the system, average the value using its neighbours and return as a new system.
///
/// # Notes
/// Assumes that the index order of the input `data` and `neighbours` are identical and that
/// the vectors are of equal size.
fn get_averaged_system(data: Vec<f64>, neighbours: Vec<Vec<usize>>) -> Vec<f64> {
    neighbours
        .into_par_iter()
        .map(|bins| average_value_of_bins(&data, &bins))
        .collect()
}

/// Average the data of bins with input indices.
///
/// # Notes
/// Assumes that all input indices are valid, ie. lie within bounds of the `data` array.
fn average_value_of_bins(data: &Vec<f64>, bins: &Vec<usize>) -> f64 {
    if bins.is_empty() {
        0.0
    } else {
        bins.iter().map(|&i| data[i]).sum::<f64>() / bins.len() as f64
    }
}

/// Get the indices of neighbouring bins to the input bin, using the neighbour "sieve".
/// Each candidate in the sieve is used along with the input bin to assert that it lies
/// within the system. If so, the 1D index of that candidate is calculated and stored.
/// The list of all valid candidates is returned.
fn get_bin_neighbours(i: usize, shape: Shape, sieve: &Vec<(isize, isize)>) -> Vec<usize> {
    let (ix, iy) = index2tuple(i, shape)
        .map(|(i, j)| (i as isize, j as isize))
        .unwrap();

    sieve
        .iter()
        .filter_map(|(ix_add, iy_add)| tuple2index(ix + ix_add, iy + iy_add, shape))
        .collect()
}

/// Use an input radius and the bin sizes to get a general "sieve" of neighbour candidates
/// for a bin. These candidates have values (ix, iy) that should be *added* to a current bin's
/// 2D coordinate, resulting in a new 2D coordinate that can be asserted to lie within
/// the system.
///
/// # Bugs
/// Currently assumes that the bin size is equal along x and y.
fn get_averaging_bin_sieve(radius: f64, [dx, dy, _]: Vec3) -> Vec<(isize, isize)> {
    let nx = (radius / dx).ceil() as isize;
    let ny = (radius / dy).ceil() as isize;

    let r2 = radius.powi(2);

    let mut bins = Vec::new();

    for ix in -nx..=nx {
        let x = ix as f64 * dx;
        for iy in -ny..=ny {
            let y = iy as f64 * dy;

            if x.powi(2) + y.powi(2) <= r2 {
                bins.push((ix, iy));
            }
        }
    }

    bins
}

mod tests {
    use super::*;

    #[test]
    fn test_average_is_calculated_from_correct_bin_indices() {
        let data = vec![10.0, 20.0, 30.0];

        assert_eq!(0.0, average_value_of_bins(&data, &vec![]));
        assert_eq!(10.0, average_value_of_bins(&data, &vec![0]));
        assert_eq!(15.0, average_value_of_bins(&data, &vec![0, 1]));
        assert_eq!(20.0, average_value_of_bins(&data, &vec![0, 2]));
        assert_eq!(20.0, average_value_of_bins(&data, &vec![0, 1, 2]));
    }

    #[test]
    fn test_getting_bin_neighbours_only_includes_bins_within_the_system() {
        let sieve = vec![(-1, -1), (0, 0), (1, 1)];
        let shape = [2, 2];

        assert_eq!(vec![0, 3], get_bin_neighbours(0, shape, &sieve));
        assert_eq!(vec![1], get_bin_neighbours(1, shape, &sieve));
        assert_eq!(vec![2], get_bin_neighbours(2, shape, &sieve));
        assert_eq!(vec![0, 3], get_bin_neighbours(3, shape, &sieve));
    }

    #[test]
    fn test_bin_sizes_of_1_with_radius_1_returns_itself_and_four_neighbours() {
        let bin_size = [1.0, 1.0, 0.0];
        let radius = 1.0;

        let bins = get_averaging_bin_sieve(radius, bin_size);

        assert_eq!(5, bins.len());
        assert!(bins.contains(&(0, 0)));
        assert!(bins.contains(&(-1, 0)));
        assert!(bins.contains(&(1, 0)));
        assert!(bins.contains(&(0, -1)));
        assert!(bins.contains(&(0, 1)));
    }

    #[test]
    fn test_bin_sizes_of_half_with_larger_radius_returns_itself_and_eight_neighbours() {
        let bin_size = [0.5, 0.5, 0.0];
        let radius = 0.75;

        let bins = get_averaging_bin_sieve(radius, bin_size);

        assert_eq!(9, bins.len());
        assert!(bins.contains(&(-1, -1)));
        assert!(bins.contains(&(-1, 1)));
        assert!(bins.contains(&(1, -1)));
        assert!(bins.contains(&(1, 1)));
    }
}
