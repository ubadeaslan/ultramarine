//! Rank-correlation statistics (deterministic, std-only).
//!
//! Spearman's rho (tie-aware via average ranks) and Kendall's tau-b
//! (tie-corrected). These measure *ranking* agreement, which is the
//! primary metric for the screening engine: we care whether the engine
//! orders emitters like experiment, not absolute values.

/// Average ranks (1-based), assigning the mean rank to tied values.
pub fn average_ranks(xs: &[f64]) -> Vec<f64> {
    let n = xs.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| xs[a].partial_cmp(&xs[b]).unwrap());
    let mut ranks = vec![0.0; n];
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && xs[idx[j]] == xs[idx[i]] {
            j += 1;
        }
        // positions i..j are tied; average of 1-based ranks
        let avg = ((i + 1 + j) as f64) / 2.0; // ((i+1)+(j))/2 over inclusive 1-based
        for k in i..j {
            ranks[idx[k]] = avg;
        }
        i = j;
    }
    ranks
}

fn pearson(a: &[f64], b: &[f64]) -> Option<f64> {
    let n = a.len();
    if n < 2 {
        return None;
    }
    let ma = a.iter().sum::<f64>() / n as f64;
    let mb = b.iter().sum::<f64>() / n as f64;
    let mut num = 0.0;
    let mut da = 0.0;
    let mut db = 0.0;
    for i in 0..n {
        let xa = a[i] - ma;
        let xb = b[i] - mb;
        num += xa * xb;
        da += xa * xa;
        db += xb * xb;
    }
    if da == 0.0 || db == 0.0 {
        return None;
    }
    Some(num / (da.sqrt() * db.sqrt()))
}

/// Spearman's rho between two equal-length samples.
pub fn spearman(x: &[f64], y: &[f64]) -> Option<f64> {
    if x.len() != y.len() || x.len() < 2 {
        return None;
    }
    pearson(&average_ranks(x), &average_ranks(y))
}

/// Kendall's tau-b between two equal-length samples, with a reference
/// noise floor.
///
/// `tie_x` encodes the experimental noise floor on x: any pair with
/// `|x_i - x_j| <= tie_x` is *unrankable by reference* and is dropped
/// entirely (it neither rewards nor penalizes the predictor). Pass
/// `tie_x = 0.0` to drop only exact reference ties. Among the remaining
/// rankable pairs, standard tau-b tie-correction is applied for ties in
/// the predictor (y).
pub fn kendall_tau_b(x: &[f64], y: &[f64], tie_x: f64) -> Option<f64> {
    let n = x.len();
    if n != y.len() || n < 2 {
        return None;
    }
    let mut concordant = 0i64;
    let mut discordant = 0i64;
    let mut ty = 0i64; // predictor ties among rankable pairs
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = x[i] - x[j];
            if dx.abs() <= tie_x {
                continue; // unrankable by reference -> excluded
            }
            let dy = y[i] - y[j];
            if dy == 0.0 {
                ty += 1;
                continue;
            }
            if dx.signum() == dy.signum() {
                concordant += 1;
            } else {
                discordant += 1;
            }
        }
    }
    let n0 = concordant + discordant + ty; // rankable reference pairs
    if n0 == 0 {
        return None;
    }
    let denom = (n0 as f64 * (n0 - ty) as f64).sqrt();
    if denom == 0.0 {
        return None;
    }
    Some((concordant - discordant) as f64 / denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_and_reversed() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((spearman(&x, &y).unwrap() - 1.0).abs() < 1e-9);
        let yr = [5.0, 4.0, 3.0, 2.0, 1.0];
        assert!((spearman(&x, &yr).unwrap() + 1.0).abs() < 1e-9);
        assert!((kendall_tau_b(&x, &y, 0.0).unwrap() - 1.0).abs() < 1e-9);
        assert!((kendall_tau_b(&x, &yr, 0.0).unwrap() + 1.0).abs() < 1e-9);
    }

    #[test]
    fn tie_floor_collapses_close_pairs() {
        // x values within 0.05 are treated as tied at threshold 0.1
        let x = [0.10, 0.12, 0.50];
        let y = [1.0, 2.0, 3.0];
        // pair (0,1) is an x-tie; (0,2) and (1,2) concordant -> tau-b = 1.0
        assert!((kendall_tau_b(&x, &y, 0.1).unwrap() - 1.0).abs() < 1e-9);
    }
}
