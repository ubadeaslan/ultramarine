//! Excited states via the Tamm-Dancoff approximation.
//!
//! Stage 3a: CIS (= TDA on a Hartree-Fock reference, no XC kernel). Builds the
//! spin-adapted CIS matrices and returns the lowest singlet/triplet excitation
//! energies — i.e. S1, T1, and ΔEST = E(S1) − E(T1), THE quantity the whole
//! TADF screen ranks on. TDA-DFT (with the XC kernel, for real accuracy on CT
//! states) layers on top of this.
//!
//! Spin-adapted CIS (Szabo & Ostlund):
//!   ¹A_{ia,jb} = (ε_a − ε_i)δ_ij δ_ab + 2(ia|jb) − (ij|ab)   [singlet]
//!   ³A_{ia,jb} = (ε_a − ε_i)δ_ij δ_ab − (ij|ab)              [triplet]

use crate::integrals::Eri;
use crate::linalg::eigh;

/// Transform the first index of a 4-index tensor: O[p..] = Σ_μ C[μ][p] E[μ..].
fn trans_first(e: &Eri, c: &[Vec<f64>]) -> Eri {
    let n = c.len();
    let mut o = vec![vec![vec![vec![0.0; n]; n]; n]; n];
    for p in 0..n {
        for q in 0..n {
            for r in 0..n {
                for s in 0..n {
                    let mut acc = 0.0;
                    for mu in 0..n {
                        acc += c[mu][p] * e[mu][q][r][s];
                    }
                    o[p][q][r][s] = acc;
                }
            }
        }
    }
    o
}

/// Cyclic axis rotation [p][q][r][s] -> [q][r][s][p].
fn rotate(e: &Eri) -> Eri {
    let n = e.len();
    let mut o = vec![vec![vec![vec![0.0; n]; n]; n]; n];
    for p in 0..n {
        for q in 0..n {
            for r in 0..n {
                for s in 0..n {
                    o[q][r][s][p] = e[p][q][r][s];
                }
            }
        }
    }
    o
}

/// AO→MO four-index transform (chemist notation preserved).
fn ao2mo(eri_ao: &Eri, c: &[Vec<f64>]) -> Eri {
    let mut t = eri_ao.clone();
    for _ in 0..4 {
        t = trans_first(&t, c);
        t = rotate(&t);
    }
    t
}

pub struct CisResult {
    pub s1: f64,         // lowest singlet excitation (Hartree)
    pub t1: f64,         // lowest triplet excitation (Hartree)
    pub singlets: Vec<f64>,
    pub triplets: Vec<f64>,
    /// Lowest-singlet amplitude vector X_{ia} (length nocc·nvir) for the
    /// transition dipole / oscillator strength.
    pub s1_vec: Vec<f64>,
}

impl CisResult {
    pub fn delta_est(&self) -> f64 {
        self.s1 - self.t1
    }
}

/// General spin-adapted TDA. `c_x` = exact-exchange fraction (1 for HF/CIS,
/// 0 for pure DFT, 0.25 for PBE0). `kernel_s`/`kernel_t` are the precomputed
/// XC-kernel blocks (ia|f_xc|jb) for singlet/triplet (zero for CIS).
///
///   ¹A = (ε_a−ε_i)δ + 2(ia|jb) − c_x(ij|ab) + kernel_s
///   ³A = (ε_a−ε_i)δ            − c_x(ij|ab) + kernel_t
#[allow(clippy::too_many_arguments)]
pub fn tda_general(
    mo_coeff: &[Vec<f64>],
    mo_energy: &[f64],
    eri_ao: &Eri,
    nocc: usize,
    c_x: f64,
    kernel_s: &[Vec<f64>],
    kernel_t: &[Vec<f64>],
) -> CisResult {
    let n = mo_energy.len();
    let nvir = n - nocc;
    let mo = ao2mo(eri_ao, mo_coeff);
    let dim = nocc * nvir;
    let idx = |i: usize, a: usize| i * nvir + a;

    let mut a_s = vec![vec![0.0; dim]; dim];
    let mut a_t = vec![vec![0.0; dim]; dim];
    for i in 0..nocc {
        for a in 0..nvir {
            let av = nocc + a;
            for j in 0..nocc {
                for b in 0..nvir {
                    let bv = nocc + b;
                    let coul = mo[i][av][j][bv]; // (ia|jb)
                    let exch = mo[i][j][av][bv]; // (ij|ab)
                    let diag = if i == j && a == b { mo_energy[av] - mo_energy[i] } else { 0.0 };
                    let (p, q) = (idx(i, a), idx(j, b));
                    a_s[p][q] = diag + 2.0 * coul - c_x * exch + kernel_s[p][q];
                    a_t[p][q] = diag - c_x * exch + kernel_t[p][q];
                }
            }
        }
    }

    let (es, vs) = eigh(&a_s);
    let (et, _) = eigh(&a_t);
    let s1_vec: Vec<f64> = (0..dim).map(|k| vs[k][0]).collect();
    CisResult { s1: es[0], t1: et[0], singlets: es, triplets: et, s1_vec }
}

/// CIS / TDA-HF (= tda_general with full exact exchange, no XC kernel).
pub fn cis(mo_coeff: &[Vec<f64>], mo_energy: &[f64], eri_ao: &Eri, nocc: usize) -> CisResult {
    let dim = nocc * (mo_energy.len() - nocc);
    let z = vec![vec![0.0; dim]; dim];
    tda_general(mo_coeff, mo_energy, eri_ao, nocc, 1.0, &z, &z)
}
