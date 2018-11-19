use crate::graphdata::XYData;

pub fn calc_autocorrelation<T: XYData>(data: &[T]) -> Vec<f64> {
    let mut values = vec![0.0; data.len()];

    for i in 0..data.len() {
        for j in i..data.len() {
            let n = j - i;

            let y0 = data[i].y();
            let y1 = data[j].y();

            for (a, b) in y0.iter().zip(y1.iter()) {
                values[n] += a * b;
            }
        }
    }

    // Rescale the values by dividing with the number of measurement points for the time lag.
    let rescaled_values = values
        .into_iter()
        .enumerate()
        .map(|(i, v)| v / (data.len() - i) as f64)
        .collect::<Vec<_>>();

    let max = rescaled_values[0];
    rescaled_values.into_iter().map(|v| v / max).collect()
}
