use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use flate2::{bufread::GzDecoder, write::GzEncoder, Compression};

use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Read, Write},
    path::Path,
};

pub type Vec2 = [f64; 2];
pub type Vec3 = [f64; 3];
pub type Shape = [u64; 2];

pub struct DensMap {
    /// Bin size in all directions.
    pub bin_size: Vec3,
    /// Origin of system along x and y.
    pub origin: Vec2,
    /// Shape of system along x and y.
    pub shape: Shape,
    /// Center of fitted spherical cap along x and y.
    pub center: Vec2,
    /// Density map data as a 1D vector, in order of x changing every index.
    pub data: Vec<f64>,
}

/// Read a density map from the input path.
///
/// If the path extension ends with '.gz' the file is assumed to be compressed with gzip
/// and is decompressed during reading.
pub fn read_densmap(path: &Path) -> Result<(DensMap, f64), io::Error> {
    let fp = File::open(&path)?;
    let mut reader = BufReader::new(fp);

    match path.extension().map(|p| p.to_str().unwrap()) {
        Some("gz") => read_densmap_from_reader(&mut GzDecoder::new(reader)),
        _ => read_densmap_from_reader(&mut reader),
    }
}

fn read_densmap_from_reader<R: Read>(reader: &mut R) -> Result<(DensMap, f64), io::Error> {
    let bin_size = [
        reader.read_f64::<LittleEndian>()?,
        reader.read_f64::<LittleEndian>()?,
        reader.read_f64::<LittleEndian>()?,
    ];

    let origin = [
        reader.read_f64::<LittleEndian>()?,
        reader.read_f64::<LittleEndian>()?,
    ];

    let shape = [
        reader.read_u64::<LittleEndian>()?,
        reader.read_u64::<LittleEndian>()?,
    ];

    let center = [
        reader.read_f64::<LittleEndian>()?,
        reader.read_f64::<LittleEndian>()?,
    ];

    let time = reader.read_f64::<LittleEndian>()?;

    let [nx, ny] = shape;
    let num_bins = nx * ny;

    let mut data: Vec<f64> = Vec::with_capacity(num_bins as usize);

    for _ in 0..num_bins {
        data.push(reader.read_f64::<LittleEndian>()?);
    }

    Ok((
        DensMap {
            bin_size,
            origin,
            shape,
            center,
            data,
        },
        time,
    ))
}

/// Write a density map to the input path.
///
/// If the extension ends with '.gz' the file will be encoded as a gzipped file.
pub fn write_densmap(path: &Path, densmap: &DensMap, time: f64) -> Result<(), io::Error> {
    let fp = File::create(&path)?;
    let mut writer = BufWriter::new(fp);

    match path.extension().map(|p| p.to_str().unwrap()) {
        Some("gz") => write_densmap_to_writer(
            &mut GzEncoder::new(writer, Compression::fast()),
            &densmap,
            time,
        ),
        _ => write_densmap_to_writer(&mut writer, &densmap, time),
    }
}

fn write_densmap_to_writer<W: Write>(
    writer: &mut W,
    densmap: &DensMap,
    time: f64,
) -> Result<(), io::Error> {
    writer.write_f64::<LittleEndian>(densmap.bin_size[0])?;
    writer.write_f64::<LittleEndian>(densmap.bin_size[1])?;
    writer.write_f64::<LittleEndian>(densmap.bin_size[2])?;

    writer.write_f64::<LittleEndian>(densmap.origin[0])?;
    writer.write_f64::<LittleEndian>(densmap.origin[1])?;

    writer.write_u64::<LittleEndian>(densmap.shape[0])?;
    writer.write_u64::<LittleEndian>(densmap.shape[1])?;

    writer.write_f64::<LittleEndian>(densmap.center[0])?;
    writer.write_f64::<LittleEndian>(densmap.center[1])?;

    writer.write_f64::<LittleEndian>(time)?;

    for v in densmap.data.iter().cloned() {
        writer.write_f64::<LittleEndian>(v)?;
    }

    Ok(())
}

/// Get the 2D grid position from a 1D index in the array.
///
/// Return `None` if the input index lies outside of the system.
pub fn index2tuple(i: usize, [nx, ny]: Shape) -> Option<(usize, usize)> {
    if i < (nx * ny) as usize {
        Some((i % nx as usize, i / nx as usize))
    } else {
        None
    }
}

/// Get the 1D array index from a 2D grid position.
///
/// Return `None` if the 2D position lies outside of the system.
pub fn tuple2index(ix: isize, iy: isize, [nx, ny]: Shape) -> Option<usize> {
    if ix >= 0 && ix < nx as isize && iy >= 0 && iy < ny as isize {
        Some((iy * nx as isize + ix) as usize)
    } else {
        None
    }
}

/// Get the 1D array index from a 2D system coordinate.
///
/// Return `None` if the 2D coordinate lies outside of the system.
pub fn coord2index(x: f64, y: f64, [dx, dy, _]: Vec3, shape: Shape) -> Option<usize> {
    let ix = (x / dx).floor() as isize;
    let iy = (y / dy).floor() as isize;
    tuple2index(ix, iy, shape)
}

#[test]
fn test_correct_index_from_coordinate_in_uniform_bin_sized_system() {
    let shape = [6, 9];
    let bin_size = [0.5, 0.5, 0.0];

    assert_eq!(Some(0), coord2index(0.0, 0.0, bin_size, shape));
    assert_eq!(Some(1), coord2index(0.5, 0.0, bin_size, shape));
    assert_eq!(Some(5), coord2index(2.5, 0.0, bin_size, shape));
    assert_eq!(Some(6), coord2index(0.0, 0.5, bin_size, shape));
    assert_eq!(Some(7), coord2index(0.5, 0.5, bin_size, shape));
    assert_eq!(Some(53), coord2index(2.9, 4.4, bin_size, shape));
    assert_eq!(None, coord2index(-0.1, 0.0, bin_size, shape));
    assert_eq!(None, coord2index(0.0, -0.1, bin_size, shape));
    assert_eq!(None, coord2index(3.0, 0.0, bin_size, shape));
    assert_eq!(None, coord2index(0.0, 4.5, bin_size, shape));
}

#[test]
fn test_correct_ix_values_from_index() {
    let shape = [6, 9];
    assert_eq!(0, index2tuple(0, shape).unwrap().0);
    assert_eq!(1, index2tuple(1, shape).unwrap().0);
    assert_eq!(5, index2tuple(5, shape).unwrap().0);
    assert_eq!(0, index2tuple(6, shape).unwrap().0);
    assert_eq!(1, index2tuple(7, shape).unwrap().0);
    assert_eq!(5, index2tuple(53, shape).unwrap().0);
    assert_eq!(None, index2tuple(54, shape));
}

#[test]
fn test_correct_iy_values_from_index() {
    let shape = [6, 9];
    assert_eq!(0, index2tuple(0, shape).unwrap().1);
    assert_eq!(0, index2tuple(1, shape).unwrap().1);
    assert_eq!(0, index2tuple(5, shape).unwrap().1);
    assert_eq!(1, index2tuple(6, shape).unwrap().1);
    assert_eq!(1, index2tuple(7, shape).unwrap().1);
    assert_eq!(8, index2tuple(53, shape).unwrap().1);
    assert_eq!(None, index2tuple(54, shape));
}

#[test]
fn test_correct_index_from_tuple() {
    let shape = [6, 9];
    assert_eq!(Some(0), tuple2index(0, 0, shape));
    assert_eq!(Some(1), tuple2index(1, 0, shape));
    assert_eq!(Some(5), tuple2index(5, 0, shape));
    assert_eq!(Some(6), tuple2index(0, 1, shape));
    assert_eq!(Some(7), tuple2index(1, 1, shape));
    assert_eq!(Some(53), tuple2index(5, 8, shape));
    assert_eq!(None, tuple2index(-1, 0, shape));
    assert_eq!(None, tuple2index(0, -1, shape));
    assert_eq!(None, tuple2index(6, 0, shape));
    assert_eq!(None, tuple2index(0, 9, shape));
}
