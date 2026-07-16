//! Atomic-orbital evaluation on real-space points (for KS-DFT grid integration).
//!
//! AO = (normalized contracted radial) × (real solid harmonic). The radial
//! coefficients must be normalized identically to what libcint uses for the
//! integrals/density matrix, so callers pass already-normalized primitives.
//! The intra-shell component ORDER must also match libcint's `_sph` order;
//! this is verified by checking S_numerical ≈ S_analytic (see cint-sys).
//!
//! Implemented: l = 0 (s), l = 1 (p). l = 2 (d) is added with the per-component
//! normalization once needed for def2-SVP.

/// A shell with fully-normalized contracted primitives.
#[derive(Clone)]
pub struct AoShell {
    pub center: [f64; 3],
    pub l: i32,
    pub prims: Vec<(f64, f64)>, // (exponent, normalized coeff)
}

/// Number of spherical functions emitted for this shell.
pub fn shell_nf(l: i32) -> usize {
    (2 * l + 1) as usize
}

use std::f64::consts::PI;

/// Append AO values at point `r` for one shell to `out`.
///
/// φ = [Σ cᵢ exp(-αᵢ r²)] · S_lm, where S_lm is the real solid harmonic
/// normalized to ∫|Y_lm|² dΩ = 1 (libcint's `_sph` convention). The radial
/// coeffs already carry libcint's CINTgto_norm, so the angular factor here is
/// what makes the evaluated AO equal libcint's basis function.
pub fn eval_shell(sh: &AoShell, r: &[f64; 3], out: &mut Vec<f64>) {
    let dx = r[0] - sh.center[0];
    let dy = r[1] - sh.center[1];
    let dz = r[2] - sh.center[2];
    let r2 = dx * dx + dy * dy + dz * dz;
    let radial: f64 = sh.prims.iter().map(|&(a, c)| c * (-a * r2).exp()).sum();
    match sh.l {
        0 => {
            let c0 = (1.0 / (4.0 * PI)).sqrt();
            out.push(radial * c0);
        }
        // libcint sph order for l=1 is (x, y, z); factor √(3/4π).
        1 => {
            let c1 = (3.0 / (4.0 * PI)).sqrt();
            out.push(radial * c1 * dx);
            out.push(radial * c1 * dy);
            out.push(radial * c1 * dz);
        }
        // libcint sph order m=-2,-1,0,1,2: xy, yz, (2z²-x²-y²), xz, x²-y²
        2 => {
            let ca = (15.0 / (4.0 * PI)).sqrt();
            let c0 = (5.0 / (16.0 * PI)).sqrt();
            let cb = (15.0 / (16.0 * PI)).sqrt();
            out.push(radial * ca * dx * dy);
            out.push(radial * ca * dy * dz);
            out.push(radial * c0 * (2.0 * dz * dz - dx * dx - dy * dy));
            out.push(radial * ca * dx * dz);
            out.push(radial * cb * (dx * dx - dy * dy));
        }
        _ => panic!("aoeval: l={} not yet implemented", sh.l),
    }
}

/// Evaluate all AOs (in shell order) at point `r`.
pub fn eval_ao(shells: &[AoShell], r: &[f64; 3]) -> Vec<f64> {
    let mut out = Vec::new();
    for sh in shells {
        eval_shell(sh, r, &mut out);
    }
    out
}

/// Evaluate AO values AND their Cartesian gradients ∇φ at point `r` (for GGA).
pub fn eval_ao_grad(shells: &[AoShell], r: &[f64; 3]) -> (Vec<f64>, Vec<[f64; 3]>) {
    let mut val = Vec::new();
    let mut grad = Vec::new();
    for sh in shells {
        let dx = r[0] - sh.center[0];
        let dy = r[1] - sh.center[1];
        let dz = r[2] - sh.center[2];
        let r2 = dx * dx + dy * dy + dz * dz;
        // radial Σ cᵢ exp(-αᵢ r²) and its gradient (= Σ cᵢ(-2αᵢ)exp · (dx,dy,dz))
        let mut radial = 0.0;
        let mut gr = [0.0; 3];
        for &(a, c) in &sh.prims {
            let e = c * (-a * r2).exp();
            radial += e;
            let g = -2.0 * a * e;
            gr[0] += g * dx;
            gr[1] += g * dy;
            gr[2] += g * dz;
        }
        match sh.l {
            0 => {
                let c0 = (1.0 / (4.0 * PI)).sqrt();
                val.push(c0 * radial);
                grad.push([c0 * gr[0], c0 * gr[1], c0 * gr[2]]);
            }
            1 => {
                let c1 = (3.0 / (4.0 * PI)).sqrt();
                // components (x, y, z): φ_q = c1·d_q·radial
                let comps = [(dx, [1.0, 0.0, 0.0]), (dy, [0.0, 1.0, 0.0]), (dz, [0.0, 0.0, 1.0])];
                for (dq, e_unit) in comps {
                    val.push(c1 * dq * radial);
                    let mut g = [0.0; 3];
                    for k in 0..3 {
                        g[k] = c1 * (dq * gr[k] + radial * e_unit[k]);
                    }
                    grad.push(g);
                }
            }
            2 => {
                let ca = (15.0 / (4.0 * PI)).sqrt();
                let c0 = (5.0 / (16.0 * PI)).sqrt();
                let cb = (15.0 / (16.0 * PI)).sqrt();
                // (constant, P, [∂P/∂x, ∂P/∂y, ∂P/∂z])
                let terms = [
                    (ca, dx * dy, [dy, dx, 0.0]),
                    (ca, dy * dz, [0.0, dz, dy]),
                    (c0, 2.0 * dz * dz - dx * dx - dy * dy, [-2.0 * dx, -2.0 * dy, 4.0 * dz]),
                    (ca, dx * dz, [dz, 0.0, dx]),
                    (cb, dx * dx - dy * dy, [2.0 * dx, -2.0 * dy, 0.0]),
                ];
                for (c, poly, dpoly) in terms {
                    val.push(c * poly * radial);
                    let mut g = [0.0; 3];
                    for k in 0..3 {
                        g[k] = c * (poly * gr[k] + radial * dpoly[k]);
                    }
                    grad.push(g);
                }
            }
            _ => panic!("aoeval grad: l={} not yet implemented", sh.l),
        }
    }
    (val, grad)
}
