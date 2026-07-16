//! Restricted Hartree-Fock driver — integral-source-agnostic.
//!
//! Consumes precomputed S, Hcore (= T + V) and the ERI tensor, so the SAME
//! code runs on hand-coded integrals (Windows) or libcint (WSL). This is the
//! first real layer of the owned engine; KS-DFT and TDDFT build on it.

use crate::integrals::Eri;
use crate::linalg::{eigh, matmul, solve_linear, transpose};

/// Pulay DIIS accelerator: extrapolates the Fock matrix from a history of
/// (F, error=FPS−SPF) pairs so SCF converges in ~10 iters instead of ~30.
pub struct Diis {
    focks: Vec<Vec<Vec<f64>>>,
    errors: Vec<Vec<Vec<f64>>>,
    max: usize,
}

impl Diis {
    pub fn new(max: usize) -> Self {
        Diis { focks: Vec::new(), errors: Vec::new(), max }
    }

    /// Store (F, error) and return the DIIS-extrapolated Fock matrix.
    pub fn update(&mut self, f: &[Vec<f64>], p: &[Vec<f64>], s: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let nbf = f.len();
        let fps = matmul(&matmul(f, p), s);
        let spf = matmul(&matmul(s, p), f);
        let mut err = vec![vec![0.0; nbf]; nbf];
        for i in 0..nbf {
            for j in 0..nbf {
                err[i][j] = fps[i][j] - spf[i][j];
            }
        }
        self.focks.push(f.to_vec());
        self.errors.push(err);
        while self.focks.len() > self.max {
            self.focks.remove(0);
            self.errors.remove(0);
        }
        let n = self.focks.len();
        if n < 2 {
            return f.to_vec();
        }
        // B matrix with Lagrange constraint
        let mut b = vec![vec![0.0; n + 1]; n + 1];
        for i in 0..n {
            for j in 0..n {
                let mut acc = 0.0;
                for u in 0..nbf {
                    for v in 0..nbf {
                        acc += self.errors[i][u][v] * self.errors[j][u][v];
                    }
                }
                b[i][j] = acc;
            }
            b[i][n] = -1.0;
            b[n][i] = -1.0;
        }
        let mut rhs = vec![0.0; n + 1];
        rhs[n] = -1.0;
        let c = solve_linear(&b, &rhs);
        let mut fnew = vec![vec![0.0; nbf]; nbf];
        for (i, fi) in self.focks.iter().enumerate() {
            for u in 0..nbf {
                for v in 0..nbf {
                    fnew[u][v] += c[i] * fi[u][v];
                }
            }
        }
        fnew
    }

    /// Largest absolute error-vector element (SCF convergence gauge).
    pub fn max_error(&self) -> f64 {
        self.errors
            .last()
            .map(|e| e.iter().flatten().fold(0.0_f64, |m, &x| m.max(x.abs())))
            .unwrap_or(1.0)
    }
}

pub struct Scf {
    pub e_total: f64,
    pub e_elec: f64,
    pub iters: usize,
    /// Final density matrix.
    pub density: Vec<Vec<f64>>,
    /// MO coefficients (columns = MOs).
    pub mo_coeff: Vec<Vec<f64>>,
    /// MO energies (ascending).
    pub mo_energy: Vec<f64>,
}

/// S^{-1/2} via symmetric orthogonalization.
pub fn s_inv_sqrt(s: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let (vals, vecs) = eigh(s);
    let n = s.len();
    let mut d = vec![vec![0.0; n]; n];
    for i in 0..n {
        d[i][i] = 1.0 / vals[i].sqrt();
    }
    matmul(&matmul(&vecs, &d), &transpose(&vecs))
}

/// Diagonalize a Fock matrix in the orthonormal basis: returns (MO energies
/// ascending, MO coefficients with columns = MOs). Reusable by HF and KS-DFT.
pub fn diagonalize_fock(f: &[Vec<f64>], x: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let xt = transpose(x);
    let fp = matmul(&matmul(&xt, f), x);
    let (eps, cp) = eigh(&fp);
    (eps, matmul(x, &cp))
}

/// Closed-shell density matrix from MO coefficients.
pub fn build_density(c: &[Vec<f64>], nocc: usize) -> Vec<Vec<f64>> {
    let n = c.len();
    let mut p = vec![vec![0.0; n]; n];
    for u in 0..n {
        for v in 0..n {
            let mut acc = 0.0;
            for a in 0..nocc {
                acc += c[u][a] * c[v][a];
            }
            p[u][v] = 2.0 * acc;
        }
    }
    p
}

/// Closed-shell RHF SCF. `nocc` = doubly-occupied orbitals (= n_elec/2).
pub fn rhf(
    s: &[Vec<f64>],
    hcore: &[Vec<f64>],
    eri: &Eri,
    nocc: usize,
    e_nuc: f64,
) -> Scf {
    let n = s.len();
    let x = s_inv_sqrt(s);
    let xt = transpose(&x);

    let mut p = vec![vec![0.0; n]; n]; // density (core guess)
    let mut e_elec = 0.0;
    let mut iters = 0;
    let mut mo_coeff = vec![vec![0.0; n]; n];
    let mut mo_energy = vec![0.0; n];

    for _ in 0..200 {
        iters += 1;
        // Fock: F = Hcore + G, G_uv = Σ_λσ P_λσ [ (uv|σλ) - 0.5 (uλ|σv) ]
        let mut f = hcore.to_vec();
        for u in 0..n {
            for v in 0..n {
                let mut g = 0.0;
                for l in 0..n {
                    for sg in 0..n {
                        g += p[l][sg] * (eri[u][v][sg][l] - 0.5 * eri[u][l][sg][v]);
                    }
                }
                f[u][v] += g;
            }
        }
        let fp = matmul(&matmul(&xt, &f), &x);
        let (eps, cp) = eigh(&fp);
        let c = matmul(&x, &cp);

        let mut p_new = vec![vec![0.0; n]; n];
        for u in 0..n {
            for v in 0..n {
                let mut acc = 0.0;
                for a in 0..nocc {
                    acc += c[u][a] * c[v][a];
                }
                p_new[u][v] = 2.0 * acc;
            }
        }

        let mut e = 0.0;
        for u in 0..n {
            for v in 0..n {
                e += 0.5 * p_new[u][v] * (hcore[u][v] + f[u][v]);
            }
        }

        let de = (e - e_elec).abs();
        p = p_new;
        e_elec = e;
        mo_coeff = c;
        mo_energy = eps;
        if de < 1e-10 {
            break;
        }
    }

    Scf {
        e_total: e_elec + e_nuc,
        e_elec,
        iters,
        density: p,
        mo_coeff,
        mo_energy,
    }
}
