//! Stage-0 from-scratch s-type Gaussian integrals (pure Rust, no C).
//! Production replaces these with libcint (see crates/cint-sys); the SCF
//! driver in `scf.rs` is integral-source-agnostic.

const PI: f64 = std::f64::consts::PI;

/// Two-electron integral tensor (uv|ws), chemist notation.
pub type Eri = Vec<Vec<Vec<Vec<f64>>>>;

/// A contracted s-type basis function (coeffs normalized so <φ|φ>=1).
#[derive(Clone)]
pub struct BasisS {
    pub center: [f64; 3],
    pub prims: Vec<(f64, f64)>, // (exponent α, contraction coeff d)
}

fn dist2(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    (0..3).map(|i| (a[i] - b[i]).powi(2)).sum()
}

fn norm_s(alpha: f64) -> f64 {
    (2.0 * alpha / PI).powf(0.75)
}

/// Error function (Abramowitz & Stegun 7.1.26).
fn erf(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let y = 1.0
        - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
            + 0.254829592)
            * t
            * (-x * x).exp();
    sign * y
}

fn boys0(t: f64) -> f64 {
    if t < 1e-12 {
        1.0 - t / 3.0
    } else {
        0.5 * (PI / t).sqrt() * erf(t.sqrt())
    }
}

fn s_prim(a: f64, ca: &[f64; 3], b: f64, cb: &[f64; 3]) -> f64 {
    let p = a + b;
    let mu = a * b / p;
    norm_s(a) * norm_s(b) * (PI / p).powf(1.5) * (-mu * dist2(ca, cb)).exp()
}

pub fn overlap(u: &BasisS, v: &BasisS) -> f64 {
    let mut s = 0.0;
    for &(a, da) in &u.prims {
        for &(b, db) in &v.prims {
            s += da * db * s_prim(a, &u.center, b, &v.center);
        }
    }
    s
}

fn kinetic(u: &BasisS, v: &BasisS) -> f64 {
    let mut t = 0.0;
    for &(a, da) in &u.prims {
        for &(b, db) in &v.prims {
            let p = a + b;
            let mu = a * b / p;
            let sp = s_prim(a, &u.center, b, &v.center);
            t += da * db * mu * (3.0 - 2.0 * mu * dist2(&u.center, &v.center)) * sp;
        }
    }
    t
}

fn nuclear(u: &BasisS, v: &BasisS, c: &[f64; 3], z: f64) -> f64 {
    let mut val = 0.0;
    for &(a, da) in &u.prims {
        for &(b, db) in &v.prims {
            let p = a + b;
            let mu = a * b / p;
            let mut pc = [0.0; 3];
            for i in 0..3 {
                pc[i] = (a * u.center[i] + b * v.center[i]) / p;
            }
            let pre = norm_s(a) * norm_s(b) * (-z) * (2.0 * PI / p)
                * (-mu * dist2(&u.center, &v.center)).exp();
            val += da * db * pre * boys0(p * dist2(&pc, c));
        }
    }
    val
}

fn eri_4(u: &BasisS, v: &BasisS, w: &BasisS, s: &BasisS) -> f64 {
    let mut val = 0.0;
    for &(a, da) in &u.prims {
        for &(b, db) in &v.prims {
            let p = a + b;
            let mut pc = [0.0; 3];
            for i in 0..3 {
                pc[i] = (a * u.center[i] + b * v.center[i]) / p;
            }
            let kab = (-a * b / p * dist2(&u.center, &v.center)).exp();
            for &(g, dg) in &w.prims {
                for &(h, dh) in &s.prims {
                    let q = g + h;
                    let mut qc = [0.0; 3];
                    for i in 0..3 {
                        qc[i] = (g * w.center[i] + h * s.center[i]) / q;
                    }
                    let kcd = (-g * h / q * dist2(&w.center, &s.center)).exp();
                    let pre = norm_s(a) * norm_s(b) * norm_s(g) * norm_s(h) * 2.0
                        * PI.powf(2.5)
                        / (p * q * (p + q).sqrt());
                    val += da * db * dg * dh * pre * kab * kcd
                        * boys0(p * q / (p + q) * dist2(&pc, &qc));
                }
            }
        }
    }
    val
}

pub fn normalize(mut bf: BasisS) -> BasisS {
    let scale = 1.0 / overlap(&bf, &bf).sqrt();
    for p in &mut bf.prims {
        p.1 *= scale;
    }
    bf
}

pub fn nuclear_repulsion(nuclei: &[([f64; 3], f64)]) -> f64 {
    let mut e = 0.0;
    for a in 0..nuclei.len() {
        for b in (a + 1)..nuclei.len() {
            e += nuclei[a].1 * nuclei[b].1 / dist2(&nuclei[a].0, &nuclei[b].0).sqrt();
        }
    }
    e
}

/// Build (S, Hcore=T+V, ERI) from hand-coded integrals.
pub fn build_matrices(
    bfs: &[BasisS],
    nuclei: &[([f64; 3], f64)],
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Eri) {
    let n = bfs.len();
    let mut s = vec![vec![0.0; n]; n];
    let mut hcore = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            s[i][j] = overlap(&bfs[i], &bfs[j]);
            let mut h = kinetic(&bfs[i], &bfs[j]);
            for (c, z) in nuclei {
                h += nuclear(&bfs[i], &bfs[j], c, *z);
            }
            hcore[i][j] = h;
        }
    }
    let mut eri = vec![vec![vec![vec![0.0; n]; n]; n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                for l in 0..n {
                    eri[i][j][k][l] = eri_4(&bfs[i], &bfs[j], &bfs[k], &bfs[l]);
                }
            }
        }
    }
    (s, hcore, eri)
}

/// STO-3G hydrogen 1s.
pub fn h_sto3g(center: [f64; 3]) -> BasisS {
    normalize(BasisS {
        center,
        prims: vec![
            (3.42525091, 0.15432897),
            (0.62391373, 0.53532814),
            (0.16885540, 0.44463454),
        ],
    })
}
