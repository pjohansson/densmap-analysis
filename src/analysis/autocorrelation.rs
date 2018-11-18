use crate::graphdata::XYData;

pub fn calc_autocorrelation<T: XYData>(data: &[T]) -> Vec<f64> {
    let mut values = vec![0.0; data.len()];

    for i in 0..data.len() {
        for j in i..data.len() {
            let n = j - i;

            let a = data[i].y();
            let b = data[j].y();

            for (val0, val1) in a.iter().zip(b.iter()) {
                values[n] += val0 * val1;
            }
        }
    }

    // let max = values[0];
    // values.into_iter().map(|v| v / max).collect()
    values
}
