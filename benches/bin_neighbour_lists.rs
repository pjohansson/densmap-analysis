use criterion::*;
use densmap::{
    average::get_system_bin_neighbours,
    densmap::{Shape, Vec3},
};

fn bin_neighbour_lists(c: &mut Criterion) {
    c.bench_function("bin_neighbours 100 200", |b| {
        let shape: Shape = [100, 200];
        let bin_size: Vec3 = [0.1, 0.1, 0.0];
        let radius = 0.3_f64;
        b.iter(|| get_system_bin_neighbours(radius, bin_size, shape))
    });
}

criterion_group!(benches, bin_neighbour_lists);
criterion_main!(benches);
