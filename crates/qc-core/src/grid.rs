//! Numerical integration grid (Becke fuzzy-cell molecular grid) for KS-DFT.
//!
//! Radial: Mura-Knowles. Angular: Gauss-Legendre(cosθ) × uniform(φ) product
//! (simple, table-free; swap to Lebedev later for efficiency). Multi-center
//! partitioning: Becke 1988 smoothed Voronoi weights.
//!
//! Validated independently of the basis by integrating analytic Gaussians
//! (∫ exp(-α r²) dr = (π/α)^{3/2}), so grid correctness is established before
//! AO/XC evaluation is layered on top.

use crate::basis::Molecule;

const PI: f64 = std::f64::consts::PI;

#[derive(Clone, Copy)]
pub struct GridPoint {
    pub r: [f64; 3],
    pub w: f64,
}

/// Gauss-Legendre nodes & weights on [-1, 1] (Newton iteration).
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    let nf = n as f64;
    let legendre = |xi: f64| -> (f64, f64) {
        // returns (P_n, P_{n-1})
        let (mut p0, mut p1) = (1.0_f64, xi);
        if n == 0 {
            return (1.0, 0.0);
        }
        for k in 2..=n {
            let kf = k as f64;
            let p2 = ((2.0 * kf - 1.0) * xi * p1 - (kf - 1.0) * p0) / kf;
            p0 = p1;
            p1 = p2;
        }
        (p1, p0)
    };
    for i in 0..n {
        let mut xi = (PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
        for _ in 0..100 {
            let (pn, pn1) = legendre(xi);
            let dp = nf * (xi * pn - pn1) / (xi * xi - 1.0);
            let dx = -pn / dp;
            xi += dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        let (_pn, pn1) = legendre(xi);
        let dp = nf * (xi * legendre(xi).0 - pn1) / (xi * xi - 1.0);
        x[i] = xi;
        w[i] = 2.0 / ((1.0 - xi * xi) * dp * dp);
    }
    (x, w)
}

/// Mura-Knowles radial quadrature: ∫_0^∞ f(r) r² dr ≈ Σ wᵢ f(rᵢ).
pub fn radial_mura_knowles(n: usize, alpha: f64) -> (Vec<f64>, Vec<f64>) {
    let mut r = vec![0.0; n];
    let mut w = vec![0.0; n];
    let np1 = n as f64 + 1.0;
    for i in 1..=n {
        let x = i as f64 / np1;
        let one_m = 1.0 - x * x * x;
        let ri = -alpha * one_m.ln();
        let drdx = alpha * 3.0 * x * x / one_m;
        r[i - 1] = ri;
        w[i - 1] = (1.0 / np1) * drdx * ri * ri; // includes r² volume factor
    }
    (r, w)
}

/// Product angular grid on the unit sphere: Σ w = 4π. (Reference; superseded by
/// the Lebedev grid which needs ~4× fewer points for the same accuracy.)
pub fn angular_product(n_theta: usize, n_phi: usize) -> Vec<([f64; 3], f64)> {
    let (ct, wt) = gauss_legendre(n_theta);
    let mut pts = Vec::with_capacity(n_theta * n_phi);
    for j in 0..n_theta {
        let cth = ct[j];
        let sth = (1.0 - cth * cth).max(0.0).sqrt();
        for k in 0..n_phi {
            let phi = 2.0 * PI * k as f64 / n_phi as f64;
            let dir = [sth * phi.cos(), sth * phi.sin(), cth];
            pts.push((dir, wt[j] * (2.0 * PI / n_phi as f64)));
        }
    }
    pts
}

/// Lebedev-302 angular grid (embedded, dumped from PySCF). Weights scaled so
/// Σ w = 4π (the file stores them normalized to 1).
pub fn angular_lebedev() -> Vec<([f64; 3], f64)> {
    const DATA: &str = include_str!("../../../data/lebedev-302.txt");
    let four_pi = 4.0 * PI;
    DATA.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let v: Vec<f64> = l.split_whitespace().map(|t| t.parse().unwrap()).collect();
            ([v[0], v[1], v[2]], v[3] * four_pi)
        })
        .collect()
}

/// Radial scale for Mura-Knowles (loosely element-dependent).
fn radial_alpha(z: i32) -> f64 {
    match z {
        1 => 5.0,      // H
        2..=10 => 5.0, // first row
        _ => 7.0,
    }
}

fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Becke (1988) fuzzy-Voronoi weight of `target` atom at point `r`.
fn becke_weight(mol: &Molecule, target: usize, ra: &[f64], rij: &[Vec<f64>]) -> f64 {
    let n = mol.atoms.len();
    if n == 1 {
        return 1.0;
    }
    let mut p = vec![1.0_f64; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let mu = (ra[i] - ra[j]) / rij[i][j];
            let mut f = mu;
            for _ in 0..3 {
                f = 1.5 * f - 0.5 * f * f * f; // Becke smoothing
            }
            p[i] *= 0.5 * (1.0 - f);
        }
    }
    let sum: f64 = p.iter().sum();
    if sum <= 0.0 {
        0.0
    } else {
        p[target] / sum
    }
}

/// Build the full molecular grid (Lebedev-302 angular × Mura-Knowles radial,
/// Becke partitioning).
pub fn molecular_grid(mol: &Molecule, n_rad: usize) -> Vec<GridPoint> {
    let ang = angular_lebedev();
    let n = mol.atoms.len();
    // precompute inter-atomic distances
    let mut rij = vec![vec![1.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            if i != j {
                rij[i][j] = dist(&mol.atoms[i].coord, &mol.atoms[j].coord);
            }
        }
    }

    let mut grid = Vec::new();
    for (ai, atom) in mol.atoms.iter().enumerate() {
        let (rr, rw) = radial_mura_knowles(n_rad, radial_alpha(atom.z));
        for ir in 0..n_rad {
            for (dir, aw) in &ang {
                let r = [
                    atom.coord[0] + rr[ir] * dir[0],
                    atom.coord[1] + rr[ir] * dir[1],
                    atom.coord[2] + rr[ir] * dir[2],
                ];
                let w_atom = rw[ir] * aw;
                // Becke partition weight
                let ra: Vec<f64> = mol.atoms.iter().map(|a| dist(&r, &a.coord)).collect();
                let p = becke_weight(mol, ai, &ra, &rij);
                let w = w_atom * p;
                if w > 1e-15 {
                    grid.push(GridPoint { r, w });
                }
            }
        }
    }
    grid
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::{Atom, Molecule};

    #[test]
    fn gauss_legendre_integrates() {
        // ∫_{-1}^{1} x² dx = 2/3
        let (x, w) = gauss_legendre(8);
        let s: f64 = x.iter().zip(&w).map(|(xi, wi)| wi * xi * xi).sum();
        assert!((s - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn single_gaussian() {
        // ∫ exp(-α r²) dV = (π/α)^{3/2}, single center.
        let mol = Molecule { atoms: vec![Atom { z: 8, coord: [0.0, 0.0, 0.0] }] };
        let grid = molecular_grid(&mol, 96);
        let alpha = 1.3;
        let integral: f64 = grid.iter().map(|g| g.w * (-alpha * dot(&g.r, &g.r)).exp()).sum();
        let exact = (PI / alpha).powf(1.5);
        assert!((integral - exact).abs() / exact < 1e-4, "got {integral}, exact {exact}");
    }

    #[test]
    fn two_center_partition() {
        // Two Gaussians on two atoms; Becke partition must reproduce the sum.
        let mol = Molecule {
            atoms: vec![
                Atom { z: 1, coord: [0.0, 0.0, 0.0] },
                Atom { z: 1, coord: [0.0, 0.0, 1.4] },
            ],
        };
        let grid = molecular_grid(&mol, 96);
        let a = 1.0;
        let f = |r: &[f64; 3]| {
            (-a * (r[0].powi(2) + r[1].powi(2) + r[2].powi(2))).exp()
                + (-a * (r[0].powi(2) + r[1].powi(2) + (r[2] - 1.4).powi(2))).exp()
        };
        let integral: f64 = grid.iter().map(|g| g.w * f(&g.r)).sum();
        let exact = 2.0 * (PI / a).powf(1.5);
        assert!((integral - exact).abs() / exact < 1e-4, "got {integral}, exact {exact}");
    }

    fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }
}
