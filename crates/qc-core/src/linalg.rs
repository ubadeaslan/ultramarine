//! Minimal dense linear algebra for the Stage-0 spike (std-only).
//! Production replaces this with `faer`.

pub fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let (n, m) = (a.len(), a[0].len());
    let mut t = vec![vec![0.0; n]; m];
    for i in 0..n {
        for j in 0..m {
            t[j][i] = a[i][j];
        }
    }
    t
}

pub fn matmul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let (n, k, m) = (a.len(), b.len(), b[0].len());
    let mut c = vec![vec![0.0; m]; n];
    for i in 0..n {
        for p in 0..k {
            let aip = a[i][p];
            if aip == 0.0 {
                continue;
            }
            for j in 0..m {
                c[i][j] += aip * b[p][j];
            }
        }
    }
    c
}

/// Solve A x = b for a small dense system (Gaussian elimination, partial pivot).
pub fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut m = a.to_vec();
    let mut x = b.to_vec();
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        m.swap(col, piv);
        x.swap(col, piv);
        let d = m[col][col];
        if d.abs() < 1e-14 {
            continue;
        }
        for r in 0..n {
            if r == col {
                continue;
            }
            let f = m[r][col] / d;
            for c in col..n {
                m[r][c] -= f * m[col][c];
            }
            x[r] -= f * x[col];
        }
    }
    for i in 0..n {
        if m[i][i].abs() > 1e-14 {
            x[i] /= m[i][i];
        }
    }
    x
}

/// Symmetric eigensolver via cyclic Jacobi rotations.
/// Returns (eigenvalues ascending, eigenvectors as COLUMNS): A = V Λ Vᵀ.
pub fn eigh(a_in: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a_in.len();
    let mut a: Vec<Vec<f64>> = a_in.to_vec();
    let mut v = vec![vec![0.0; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }
    if n == 1 {
        return (vec![a[0][0]], v);
    }

    for _sweep in 0..100 {
        // off-diagonal magnitude
        let mut off = 0.0;
        for p in 0..n {
            for q in (p + 1)..n {
                off += a[p][q] * a[p][q];
            }
        }
        if off.sqrt() < 1e-14 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                // rotate rows/cols p,q
                for k in 0..n {
                    let akp = a[k][p];
                    let akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[p][k];
                    let aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let vkp = v[k][p];
                    let vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }

    let mut vals: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&x, &y| vals[x].partial_cmp(&vals[y]).unwrap());
    let sorted_vals: Vec<f64> = idx.iter().map(|&i| vals[i]).collect();
    let mut sorted_vecs = vec![vec![0.0; n]; n];
    for (newc, &oldc) in idx.iter().enumerate() {
        for r in 0..n {
            sorted_vecs[r][newc] = v[r][oldc];
        }
    }
    vals = sorted_vals;
    (vals, sorted_vecs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diag_2x2() {
        // [[2,0],[0,3]] -> eigenvalues 2,3
        let a = vec![vec![2.0, 0.0], vec![0.0, 3.0]];
        let (vals, _) = eigh(&a);
        assert!((vals[0] - 2.0).abs() < 1e-12 && (vals[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn symmetric_reconstruct() {
        let a = vec![vec![4.0, 1.0], vec![1.0, 3.0]];
        let (vals, vecs) = eigh(&a);
        // V Λ Vᵀ == A
        let n = 2;
        let mut d = vec![vec![0.0; n]; n];
        for i in 0..n {
            d[i][i] = vals[i];
        }
        let recon = matmul(&matmul(&vecs, &d), &transpose(&vecs));
        for i in 0..n {
            for j in 0..n {
                assert!((recon[i][j] - a[i][j]).abs() < 1e-10);
            }
        }
    }
}
