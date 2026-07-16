//! RHF via libcint + owned SCF, driven by a parsed basis set.
//!
//! Validates H2 and H2O / STO-3G. The molecule→libcint (atm/bas/env) builder
//! is general: any element/basis the parser knows can be fed to the engine.

use std::os::raw::{c_double, c_int, c_void};

use qc_core::aoeval::{eval_ao, AoShell};
use qc_core::basis::{Atom, BasisSet, Molecule};
use qc_core::grid::molecular_grid;
use qc_core::integrals::Eri;
use qc_core::scf::{build_density, diagonalize_fock, rhf, s_inv_sqrt, Diis};
use qc_core::tda::{cis, tda_general};

const HARTREE_EV: f64 = 27.211386245988;

const PTR_ENV_START: usize = 20;
const ATM_SLOTS: usize = 6;
const BAS_SLOTS: usize = 8;

extern "C" {
    fn CINTgto_norm(n: c_int, a: c_double) -> c_double;
    fn cint1e_ovlp_sph(buf: *mut c_double, shls: *const c_int, atm: *const c_int, natm: c_int, bas: *const c_int, nbas: c_int, env: *const c_double) -> c_int;
    fn cint1e_kin_sph(buf: *mut c_double, shls: *const c_int, atm: *const c_int, natm: c_int, bas: *const c_int, nbas: c_int, env: *const c_double) -> c_int;
    fn cint1e_nuc_sph(buf: *mut c_double, shls: *const c_int, atm: *const c_int, natm: c_int, bas: *const c_int, nbas: c_int, env: *const c_double) -> c_int;
    fn cint1e_r_sph(buf: *mut c_double, shls: *const c_int, atm: *const c_int, natm: c_int, bas: *const c_int, nbas: c_int, env: *const c_double) -> c_int;
    fn cint2e_sph(buf: *mut c_double, shls: *const c_int, atm: *const c_int, natm: c_int, bas: *const c_int, nbas: c_int, env: *const c_double, opt: *mut c_void) -> c_int;
}

// libxc (via C shim).
extern "C" {
    fn qc_xc_init(id: c_int) -> *mut c_void;
    fn qc_xc_lda_exc_vxc(f: *mut c_void, n: c_int, rho: *const c_double, exc: *mut c_double, vxc: *mut c_double);
    fn qc_xc_gga_exc_vxc(f: *mut c_void, n: c_int, rho: *const c_double, sigma: *const c_double, exc: *mut c_double, vrho: *mut c_double, vsigma: *mut c_double);
    fn qc_xc_init_pol(id: c_int) -> *mut c_void;
    fn qc_xc_lda_fxc(f: *mut c_void, n: c_int, rho: *const c_double, v2rho2: *mut c_double);
    fn qc_xc_cam_coef(f: *mut c_void, omega: *mut c_double, alpha: *mut c_double, beta: *mut c_double);
    fn qc_xc_free(f: *mut c_void);
}

/// Validate the libxc binding: LDA exchange exc(ρ) vs the analytic Slater value
/// ε_x = -(3/4)(3/π)^{1/3} ρ^{1/3}.
fn xc_smoke() {
    const XC_LDA_X: c_int = 1;
    let rho = [1.0_f64, 0.5, 0.1];
    let mut exc = [0.0_f64; 3];
    let mut vxc = [0.0_f64; 3];
    unsafe {
        let f = qc_xc_init(XC_LDA_X);
        assert!(!f.is_null(), "xc_func_init failed");
        qc_xc_lda_exc_vxc(f, 3, rho.as_ptr(), exc.as_mut_ptr(), vxc.as_mut_ptr());
        qc_xc_free(f);
    }
    let cx = -0.75 * (3.0 / std::f64::consts::PI).powf(1.0 / 3.0);
    let analytic = |r: f64| cx * r.powf(1.0 / 3.0);
    println!("\n=== libxc LDA_X (Slater exchange) smoke test ===");
    let mut ok = true;
    for i in 0..3 {
        let a = analytic(rho[i]);
        let d = (exc[i] - a).abs();
        if d > 1e-9 {
            ok = false;
        }
        println!("  ρ={:<4} exc_libxc={:.8}  analytic={:.8}  Δ={:.1e}", rho[i], exc[i], a, d);
    }
    println!("  {}", if ok { "✅ libxc FFI validated." } else { "❌ mismatch." });
}

struct Mol {
    atm: Vec<c_int>,
    bas: Vec<c_int>,
    env: Vec<f64>,
    natm: c_int,
    nbas: c_int,
}

impl Mol {
    fn shell_l(&self, i: usize) -> i32 { self.bas[i * BAS_SLOTS + 1] }
    fn shell_nctr(&self, i: usize) -> i32 { self.bas[i * BAS_SLOTS + 3] }
    fn shell_nf(&self, i: usize) -> usize { ((2 * self.shell_l(i) + 1) * self.shell_nctr(i)) as usize }
    fn ao_offsets(&self) -> (Vec<usize>, usize) {
        let mut off = Vec::new();
        let mut acc = 0;
        for i in 0..self.nbas as usize { off.push(acc); acc += self.shell_nf(i); }
        (off, acc)
    }
    fn nuclear_repulsion(&self) -> f64 {
        let mut e = 0.0;
        let na = self.natm as usize;
        for i in 0..na {
            let zi = self.atm[i * ATM_SLOTS] as f64;
            let pi = self.atm[i * ATM_SLOTS + 1] as usize;
            for j in (i + 1)..na {
                let zj = self.atm[j * ATM_SLOTS] as f64;
                let pj = self.atm[j * ATM_SLOTS + 1] as usize;
                let r2: f64 = (0..3).map(|k| (self.env[pi + k] - self.env[pj + k]).powi(2)).sum();
                e += zi * zj / r2.sqrt();
            }
        }
        e
    }
}

/// Fully normalized contracted coefficients (libcint primitive norm via
/// CINTgto_norm + analytic contraction renorm to <φ|φ>=1). Shared by the
/// libcint env builder and the grid AO-shell builder so both use IDENTICAL
/// normalization.
fn normalized_coeffs(l: c_int, exps: &[f64], coeffs: &[f64]) -> Vec<f64> {
    let n = exps.len();
    let mut c: Vec<f64> = exps
        .iter()
        .zip(coeffs)
        .map(|(&a, &d)| d * unsafe { CINTgto_norm(l, a) })
        .collect();
    let mut self_ov = 0.0;
    for i in 0..n {
        for j in 0..n {
            let (ai_, aj) = (exps[i], exps[j]);
            self_ov += c[i] * c[j] * (2.0 * (ai_ * aj).sqrt() / (ai_ + aj)).powf(l as f64 + 1.5);
        }
    }
    let renorm = 1.0 / self_ov.sqrt();
    for ci in c.iter_mut() {
        *ci *= renorm;
    }
    c
}

/// Build libcint arrays for a molecule + basis.
fn build_mol(mol: &Molecule, bs: &BasisSet) -> Mol {
    let mut env = vec![0.0_f64; PTR_ENV_START];
    let mut atm: Vec<c_int> = Vec::new();
    for atom in &mol.atoms {
        let cptr = env.len() as c_int;
        env.extend_from_slice(&atom.coord);
        atm.extend_from_slice(&[atom.z, cptr, 0, 0, 0, 0]);
    }
    let mut bas: Vec<c_int> = Vec::new();
    for (ai, atom) in mol.atoms.iter().enumerate() {
        let eb = bs.elements.get(&atom.z).unwrap_or_else(|| panic!("no basis for Z={}", atom.z));
        for sh in &eb.shells {
            let exp_ptr = env.len() as c_int;
            env.extend_from_slice(&sh.exps);
            let coeff_ptr = env.len() as c_int;
            env.extend_from_slice(&normalized_coeffs(sh.l, &sh.exps, &sh.coeffs));
            bas.extend_from_slice(&[ai as c_int, sh.l, sh.exps.len() as c_int, 1, 0, exp_ptr, coeff_ptr, 0]);
        }
    }
    let natm = mol.atoms.len() as c_int;
    let nbas = (bas.len() / BAS_SLOTS) as c_int;
    Mol { atm, bas, env, natm, nbas }
}

/// AO shells (same order + normalization as libcint) for grid evaluation.
fn build_ao_shells(mol: &Molecule, bs: &BasisSet) -> Vec<AoShell> {
    let mut shells = Vec::new();
    for atom in &mol.atoms {
        let eb = bs.elements.get(&atom.z).unwrap();
        for sh in &eb.shells {
            let c = normalized_coeffs(sh.l, &sh.exps, &sh.coeffs);
            shells.push(AoShell {
                center: atom.coord,
                l: sh.l,
                prims: sh.exps.iter().cloned().zip(c).collect(),
            });
        }
    }
    shells
}

/// Validate grid + AO evaluation by reproducing the analytic overlap matrix:
/// S_num[μν] = Σ_grid w φ_μ φ_ν  ≈  S_libcint.
fn validate_grid_ao(name: &str, mol: &Molecule, bs: &BasisSet) {
    let m = build_mol(mol, bs);
    let s_ana = one_electron(&m, cint1e_ovlp_sph);
    let nbf = s_ana.len();
    let shells = build_ao_shells(mol, bs);
    let grid = molecular_grid(mol, 96);
    let mut s_num = vec![vec![0.0_f64; nbf]; nbf];
    for g in &grid {
        let ao = eval_ao(&shells, &g.r);
        for i in 0..nbf {
            if ao[i] != 0.0 {
                for j in 0..nbf {
                    s_num[i][j] += g.w * ao[i] * ao[j];
                }
            }
        }
    }
    let mut maxd = 0.0;
    let mut worst = (0, 0);
    for i in 0..nbf {
        for j in 0..nbf {
            let d = (s_num[i][j] - s_ana[i][j]).abs();
            if d > maxd {
                maxd = d;
                worst = (i, j);
            }
        }
    }
    println!(
        "{:<6} grid+AO: {} pts, max|ΔS|={:.2e} at {:?}  S_num[0][0]={:.4} S_ana[0][0]={:.4} ratio={:.4}  {}",
        name,
        grid.len(),
        maxd,
        worst,
        s_num[0][0],
        s_ana[0][0],
        s_num[0][0] / s_ana[0][0],
        if maxd < 5e-4 { "✅" } else { "❌" }
    );
}

type Int1e = unsafe extern "C" fn(*mut c_double, *const c_int, *const c_int, c_int, *const c_int, c_int, *const c_double) -> c_int;

fn one_electron(mol: &Mol, func: Int1e) -> Vec<Vec<f64>> {
    let (off, nbf) = mol.ao_offsets();
    let mut m = vec![vec![0.0; nbf]; nbf];
    for ish in 0..mol.nbas as usize {
        let di = mol.shell_nf(ish);
        for jsh in 0..mol.nbas as usize {
            let dj = mol.shell_nf(jsh);
            let shls = [ish as c_int, jsh as c_int];
            let mut buf = vec![0.0_f64; di * dj];
            unsafe { func(buf.as_mut_ptr(), shls.as_ptr(), mol.atm.as_ptr(), mol.natm, mol.bas.as_ptr(), mol.nbas, mol.env.as_ptr()); }
            for p in 0..di { for q in 0..dj { m[off[ish] + p][off[jsh] + q] = buf[p + di * q]; } }
        }
    }
    m
}

fn two_electron(mol: &Mol) -> Eri {
    let (off, nbf) = mol.ao_offsets();
    let mut eri = vec![vec![vec![vec![0.0; nbf]; nbf]; nbf]; nbf];
    let ns = mol.nbas as usize;
    for ish in 0..ns {
        let di = mol.shell_nf(ish);
        for jsh in 0..ns {
            let dj = mol.shell_nf(jsh);
            for ksh in 0..ns {
                let dk = mol.shell_nf(ksh);
                for lsh in 0..ns {
                    let dl = mol.shell_nf(lsh);
                    let shls = [ish as c_int, jsh as c_int, ksh as c_int, lsh as c_int];
                    let mut buf = vec![0.0_f64; di * dj * dk * dl];
                    unsafe { cint2e_sph(buf.as_mut_ptr(), shls.as_ptr(), mol.atm.as_ptr(), mol.natm, mol.bas.as_ptr(), mol.nbas, mol.env.as_ptr(), std::ptr::null_mut()); }
                    for p in 0..di { for q in 0..dj { for r in 0..dk { for s in 0..dl {
                        let idx = p + di * (q + dj * (r + dk * s));
                        eri[off[ish] + p][off[jsh] + q][off[ksh] + r][off[lsh] + s] = buf[idx];
                    }}}}
                }
            }
        }
    }
    eri
}

fn run(name: &str, mol: &Molecule, bs: &BasisSet, reference: f64, tol: f64) {
    let m = build_mol(mol, bs);
    let s = one_electron(&m, cint1e_ovlp_sph);
    let t = one_electron(&m, cint1e_kin_sph);
    let v = one_electron(&m, cint1e_nuc_sph);
    let n = s.len();
    let mut hcore = vec![vec![0.0; n]; n];
    for i in 0..n { for j in 0..n { hcore[i][j] = t[i][j] + v[i][j]; } }
    let eri = two_electron(&m);
    let nocc = (mol.n_electrons() / 2) as usize;
    let scf = rhf(&s, &hcore, &eri, nocc, m.nuclear_repulsion());
    let ok = (scf.e_total - reference).abs() < tol;
    println!(
        "{:<10} nbf={:<3} iters={:<3} E={:.6}  ref={:.6}  Δ={:.2e}  {}",
        name, n, scf.iters, scf.e_total, reference, (scf.e_total - reference).abs(),
        if ok { "✅" } else { "❌" }
    );
}

/// Coulomb matrix J[μν] = Σ_λσ P_λσ (μν|λσ).
fn coulomb(p: &[Vec<f64>], eri: &Eri) -> Vec<Vec<f64>> {
    let n = p.len();
    let mut j = vec![vec![0.0; n]; n];
    for u in 0..n {
        for v in 0..n {
            let mut s = 0.0;
            for l in 0..n {
                for sg in 0..n {
                    s += p[l][sg] * eri[u][v][l][sg];
                }
            }
            j[u][v] = s;
        }
    }
    j
}

/// Density on the grid from the density matrix: ρ_g = Σ_μν P_μν φ_μ φ_ν.
fn density_on_grid(p: &[Vec<f64>], ao_grid: &[Vec<f64>]) -> Vec<f64> {
    let nbf = p.len();
    ao_grid
        .iter()
        .map(|ao| {
            let mut r = 0.0;
            for mu in 0..nbf {
                let mut t = 0.0;
                for nu in 0..nbf {
                    t += p[mu][nu] * ao[nu];
                }
                r += ao[mu] * t;
            }
            r.max(0.0)
        })
        .collect()
}

/// Exchange matrix K[μν] = Σ_λσ P_λσ (μλ|σν) (for hybrid functionals).
fn exchange(p: &[Vec<f64>], eri: &Eri) -> Vec<Vec<f64>> {
    let n = p.len();
    let mut k = vec![vec![0.0; n]; n];
    for u in 0..n {
        for v in 0..n {
            let mut s = 0.0;
            for l in 0..n {
                for sg in 0..n {
                    s += p[l][sg] * eri[u][l][sg][v];
                }
            }
            k[u][v] = s;
        }
    }
    k
}

/// Functional selector.
#[derive(Clone, Copy)]
enum Xc {
    Lda,
    Pbe,
    Pbe0,
    CamB3lyp,
}

impl Xc {
    fn is_gga(self) -> bool {
        !matches!(self, Xc::Lda)
    }
    fn func_ids(self) -> Vec<c_int> {
        match self {
            Xc::Lda => vec![1, 7],     // LDA_X + LDA_C_VWN
            Xc::Pbe => vec![101, 130], // GGA_X_PBE + GGA_C_PBE
            Xc::Pbe0 => vec![406],     // HYB_GGA_XC_PBEH
            Xc::CamB3lyp => vec![433], // HYB_GGA_XC_CAM_B3LYP
        }
    }
    fn exx(self) -> f64 {
        match self {
            Xc::Pbe0 => 0.25,
            _ => 0.0, // CAM exchange handled via range separation
        }
    }
    fn label(self) -> &'static str {
        match self {
            Xc::Lda => "SVWN5",
            Xc::Pbe => "PBE",
            Xc::Pbe0 => "PBE0",
            Xc::CamB3lyp => "CAM-B3LYP",
        }
    }
}

/// XC energy + potential matrix from the density matrix (LDA or GGA).
fn xc_eval(
    p: &[Vec<f64>],
    ao_grid: &[Vec<f64>],
    dao_grid: &[Vec<[f64; 3]>],
    w: &[f64],
    nbf: usize,
    xc: Xc,
) -> (f64, Vec<Vec<f64>>) {
    let ng = ao_grid.len();
    let gga = xc.is_gga();
    use rayon::prelude::*;
    // step 1: density ρ (+ ∇ρ for GGA) on each grid point — parallel.
    let dens: Vec<(f64, [f64; 3])> = (0..ng)
        .into_par_iter()
        .map(|g| {
            let ao = &ao_grid[g];
            // AO screening: only basis funcs that are non-negligible at this point.
            let sig: Vec<usize> = (0..nbf).filter(|&mu| ao[mu].abs() > 1e-12).collect();
            let mut pf = vec![0.0; nbf]; // pf[μ] = Σ_ν P_μν φ_ν
            for &mu in &sig {
                let mut s = 0.0;
                for &nu in &sig {
                    s += p[mu][nu] * ao[nu];
                }
                pf[mu] = s;
            }
            let mut r = 0.0;
            for &mu in &sig {
                r += ao[mu] * pf[mu];
            }
            let mut gr = [0.0; 3];
            if gga {
                let dao = &dao_grid[g];
                for k in 0..3 {
                    let mut s = 0.0;
                    for &mu in &sig {
                        s += pf[mu] * dao[mu][k];
                    }
                    gr[k] = 2.0 * s; // ∇ρ = 2 Σ_μ pf[μ] ∇φ_μ
                }
            }
            (r.max(0.0), gr)
        })
        .collect();
    let rho: Vec<f64> = dens.iter().map(|d| d.0).collect();
    let grad_rho: Vec<[f64; 3]> = dens.iter().map(|d| d.1).collect();
    let sigma: Vec<f64> = grad_rho.iter().map(|g| g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).collect();

    // libxc
    let mut exc = vec![0.0; ng];
    let mut vrho = vec![0.0; ng];
    let mut vsigma = vec![0.0; ng];
    for &id in &xc.func_ids() {
        let mut e = vec![0.0; ng];
        unsafe {
            let f = qc_xc_init(id);
            if gga {
                let mut vr = vec![0.0; ng];
                let mut vs = vec![0.0; ng];
                qc_xc_gga_exc_vxc(f, ng as c_int, rho.as_ptr(), sigma.as_ptr(), e.as_mut_ptr(), vr.as_mut_ptr(), vs.as_mut_ptr());
                for g in 0..ng {
                    exc[g] += e[g];
                    vrho[g] += vr[g];
                    vsigma[g] += vs[g];
                }
            } else {
                let mut vr = vec![0.0; ng];
                qc_xc_lda_exc_vxc(f, ng as c_int, rho.as_ptr(), e.as_mut_ptr(), vr.as_mut_ptr());
                for g in 0..ng {
                    exc[g] += e[g];
                    vrho[g] += vr[g];
                }
            }
            qc_xc_free(f);
        }
    }

    // step 3: E_xc (scalar) + V_xc (matrix) — parallel fold/reduce over points.
    (0..ng)
        .into_par_iter()
        .fold(
            || (0.0_f64, vec![vec![0.0_f64; nbf]; nbf]),
            |mut acc, g| {
                acc.0 += w[g] * rho[g] * exc[g];
                let ao = &ao_grid[g];
                let sig: Vec<usize> = (0..nbf).filter(|&mu| ao[mu].abs() > 1e-12).collect();
                let wv = w[g] * vrho[g];
                let vmat = &mut acc.1;
                if gga {
                    let dao = &dao_grid[g];
                    let gr = grad_rho[g];
                    let twovs = 2.0 * w[g] * vsigma[g];
                    let mut gd = vec![0.0; nbf]; // ∇ρ · ∇φ_μ
                    for &mu in &sig {
                        gd[mu] = gr[0] * dao[mu][0] + gr[1] * dao[mu][1] + gr[2] * dao[mu][2];
                    }
                    for &mu in &sig {
                        let am = ao[mu];
                        let gm = gd[mu];
                        for &nu in &sig {
                            vmat[mu][nu] += wv * am * ao[nu] + twovs * (gm * ao[nu] + gd[nu] * am);
                        }
                    }
                } else if wv != 0.0 {
                    for &mu in &sig {
                        let wva = wv * ao[mu];
                        for &nu in &sig {
                            vmat[mu][nu] += wva * ao[nu];
                        }
                    }
                }
                acc
            },
        )
        .reduce(
            || (0.0_f64, vec![vec![0.0_f64; nbf]; nbf]),
            |mut a, b| {
                a.0 += b.0;
                for u in 0..nbf {
                    for v in 0..nbf {
                        a.1[u][v] += b.1[u][v];
                    }
                }
                a
            },
        )
}

/// Old LDA-only helper kept for the standalone density path (unused now).
#[allow(dead_code)]
fn xc_lda(rho: &[f64], ao_grid: &[Vec<f64>], w: &[f64], nbf: usize) -> (f64, Vec<Vec<f64>>) {
    let ng = rho.len();
    let mut exc = vec![0.0; ng];
    let mut vxc = vec![0.0; ng];
    for id in [1, 7] {
        // XC_LDA_X, XC_LDA_C_VWN
        let (mut e, mut v) = (vec![0.0; ng], vec![0.0; ng]);
        unsafe {
            let f = qc_xc_init(id);
            qc_xc_lda_exc_vxc(f, ng as c_int, rho.as_ptr(), e.as_mut_ptr(), v.as_mut_ptr());
            qc_xc_free(f);
        }
        for g in 0..ng {
            exc[g] += e[g];
            vxc[g] += v[g];
        }
    }
    let mut e_xc = 0.0;
    let mut vmat = vec![vec![0.0; nbf]; nbf];
    for g in 0..ng {
        e_xc += w[g] * rho[g] * exc[g];
        let wv = w[g] * vxc[g];
        if wv == 0.0 {
            continue;
        }
        let ao = &ao_grid[g];
        for mu in 0..nbf {
            let a = ao[mu];
            if a != 0.0 {
                let wva = wv * a;
                for nu in 0..nbf {
                    vmat[mu][nu] += wva * ao[nu];
                }
            }
        }
    }
    (e_xc, vmat)
}

/// Schwarz bound Q[ish][jsh] = sqrt(max|(ij|ij)|) for integral screening.
fn schwarz_matrix(mol: &Mol) -> Vec<Vec<f64>> {
    let ns = mol.nbas as usize;
    let mut q = vec![vec![0.0; ns]; ns];
    let maxnf = (0..ns).map(|i| mol.shell_nf(i)).max().unwrap_or(1);
    let mut buf = vec![0.0_f64; maxnf.pow(4)];
    for ish in 0..ns {
        let di = mol.shell_nf(ish);
        for jsh in 0..ns {
            let dj = mol.shell_nf(jsh);
            let shls = [ish as c_int, jsh as c_int, ish as c_int, jsh as c_int];
            unsafe {
                cint2e_sph(buf.as_mut_ptr(), shls.as_ptr(), mol.atm.as_ptr(), mol.natm, mol.bas.as_ptr(), mol.nbas, mol.env.as_ptr(), std::ptr::null_mut());
            }
            let mut mx = 0.0_f64;
            for p in 0..di {
                for qq in 0..dj {
                    let idx = p + di * (qq + dj * (p + di * qq));
                    mx = mx.max(buf[idx].abs());
                }
            }
            q[ish][jsh] = mx.sqrt();
        }
    }
    q
}

/// Direct J/K build from libcint shell quartets with Schwarz screening,
/// 8-fold permutational symmetry, and rayon parallelism — no O(N⁴) storage.
/// For integral g=(ab|cd): J[a][b]+=P[c][d]g, K[a][d]+=P[b][c]g (validated rule).
/// Unique shell quartets (ish>=jsh, ksh>=lsh, ij>=kl) with canonical element
/// ranges; each unique integral scattered to its (deduped) 8 permutations.
fn jk_direct(mol: &Mol, p: &[Vec<f64>], schwarz: &[Vec<f64>], thresh: f64, want_k: bool) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    use rayon::prelude::*;
    let (off, nbf) = mol.ao_offsets();
    let ns = mol.nbas as usize;
    let maxnf = (0..ns).map(|i| mol.shell_nf(i)).max().unwrap_or(1);
    let sp = |a: usize, b: usize| a * (a + 1) / 2 + b; // shell-pair index, a>=b

    let zeros = || (vec![vec![0.0_f64; nbf]; nbf], vec![vec![0.0_f64; nbf]; nbf]);
    (0..ns)
        .into_par_iter()
        .fold(zeros, |mut acc, ish| {
            let (ref mut j, ref mut k) = acc;
            let di = mol.shell_nf(ish);
            let mut buf = vec![0.0_f64; maxnf.pow(4)];
            for jsh in 0..=ish {
                let dj = mol.shell_nf(jsh);
                let ij = sp(ish, jsh);
                let qij = schwarz[ish][jsh];
                let eq_ij = ish == jsh;
                for ksh in 0..ns {
                    let dk = mol.shell_nf(ksh);
                    for lsh in 0..=ksh {
                        let kl = sp(ksh, lsh);
                        if kl > ij {
                            continue;
                        }
                        if qij * schwarz[ksh][lsh] < thresh {
                            continue;
                        }
                        let dl = mol.shell_nf(lsh);
                        let eq_kl = ksh == lsh;
                        let eq_braket = ij == kl;
                        let shls = [ish as c_int, jsh as c_int, ksh as c_int, lsh as c_int];
                        unsafe {
                            cint2e_sph(buf.as_mut_ptr(), shls.as_ptr(), mol.atm.as_ptr(), mol.natm, mol.bas.as_ptr(), mol.nbas, mol.env.as_ptr(), std::ptr::null_mut());
                        }
                        for ip in 0..di {
                            let mu = off[ish] + ip;
                            let q_hi = if eq_ij { ip } else { dj - 1 };
                            for iq in 0..=q_hi {
                                let nu = off[jsh] + iq;
                                for ir in 0..dk {
                                    let la = off[ksh] + ir;
                                    let s_hi = if eq_kl { ir } else { dl - 1 };
                                    for is in 0..=s_hi {
                                        // bra==ket: keep only compound (ip,iq) >= (ir,is)
                                        if eq_braket && (ir * dl + is) > (ip * dj + iq) {
                                            continue;
                                        }
                                        let si = off[lsh] + is;
                                        let g = buf[ip + di * (iq + dj * (ir + dk * is))];
                                        let perms = [
                                            [mu, nu, la, si], [nu, mu, la, si],
                                            [mu, nu, si, la], [nu, mu, si, la],
                                            [la, si, mu, nu], [si, la, mu, nu],
                                            [la, si, nu, mu], [si, la, nu, mu],
                                        ];
                                        let mut seen = [[usize::MAX; 4]; 8];
                                        let mut nseen = 0;
                                        for pm in perms {
                                            if seen[..nseen].iter().any(|x| *x == pm) {
                                                continue;
                                            }
                                            seen[nseen] = pm;
                                            nseen += 1;
                                            let [a, b, c, d] = pm;
                                            j[a][b] += p[c][d] * g;
                                            if want_k {
                                                k[a][d] += p[b][c] * g;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            acc
        })
        .reduce(zeros, |mut a, b| {
            for u in 0..nbf {
                for v in 0..nbf {
                    a.0[u][v] += b.0[u][v];
                    if want_k {
                        a.1[u][v] += b.1[u][v];
                    }
                }
            }
            a
        })
}

/// Direct-SCF KS-DFT: no stored ERI tensor (scales to large molecules).
fn ks_direct(name: &str, mol: &Molecule, bs: &BasisSet, xc: Xc) {
    let t0 = std::time::Instant::now();
    let m = build_mol(mol, bs);
    let s = one_electron(&m, cint1e_ovlp_sph);
    let t = one_electron(&m, cint1e_kin_sph);
    let v = one_electron(&m, cint1e_nuc_sph);
    let nbf = s.len();
    let mut hcore = vec![vec![0.0; nbf]; nbf];
    for i in 0..nbf {
        for j in 0..nbf {
            hcore[i][j] = t[i][j] + v[i][j];
        }
    }
    let e_nuc = m.nuclear_repulsion();
    let nocc = (mol.n_electrons() / 2) as usize;
    let exx = xc.exx();
    let schwarz = schwarz_matrix(&m);

    let shells = build_ao_shells(mol, bs);
    let grid = molecular_grid(mol, 96);
    let w: Vec<f64> = grid.iter().map(|g| g.w).collect();
    let ao_grid: Vec<Vec<f64>> = grid.iter().map(|g| eval_ao(&shells, &g.r)).collect();
    let dao_grid: Vec<Vec<[f64; 3]>> = if xc.is_gga() {
        grid.iter().map(|g| qc_core::aoeval::eval_ao_grad(&shells, &g.r).1).collect()
    } else {
        Vec::new()
    };

    let x = s_inv_sqrt(&s);
    let mut p = vec![vec![0.0; nbf]; nbf];
    let mut diis = Diis::new(8);
    let mut e_old = 0.0;
    let mut iters = 0;
    let verbose = nbf > 50; // print per-iteration progress for big systems
    for _ in 0..200 {
        iters += 1;
        let it0 = std::time::Instant::now();
        let (jm, km) = jk_direct(&m, &p, &schwarz, 1e-10, exx > 0.0);
        let (e_xc, vxc) = xc_eval(&p, &ao_grid, &dao_grid, &w, nbf, xc);
        let mut f = vec![vec![0.0; nbf]; nbf];
        let mut e1 = 0.0;
        let mut ej = 0.0;
        let mut ek = 0.0;
        for u in 0..nbf {
            for vv in 0..nbf {
                let mut fuv = hcore[u][vv] + jm[u][vv] + vxc[u][vv];
                if exx > 0.0 {
                    fuv -= exx * 0.5 * km[u][vv];
                    ek += 0.25 * p[u][vv] * km[u][vv];
                }
                f[u][vv] = fuv;
                e1 += p[u][vv] * hcore[u][vv];
                ej += 0.5 * p[u][vv] * jm[u][vv];
            }
        }
        let e = e1 + ej - exx * ek + e_xc + e_nuc;
        let f_ext = diis.update(&f, &p, &s);
        let (_eps, c) = diagonalize_fock(&f_ext, &x);
        p = build_density(&c, nocc);
        if verbose {
            eprintln!(
                "    iter {:3}  E={:.6}  dE={:.2e}  diis_err={:.2e}  ({:.1}s/iter)",
                iters,
                e,
                (e - e_old).abs(),
                diis.max_error(),
                it0.elapsed().as_secs_f64()
            );
        }
        if iters > 1 && (e - e_old).abs() < 1e-8 && diis.max_error() < 1e-5 {
            e_old = e;
            break;
        }
        e_old = e;
    }
    println!(
        "{:<6} {:<6} [direct+DIIS]  nbf={}  iters={}  E={:.6} Ha   ({:.1}s)",
        name, xc.label(), nbf, iters, e_old, t0.elapsed().as_secs_f64()
    );
}

/// Closed-shell KS-DFT (LDA / PBE / PBE0) via libcint + libxc + owned SCF.
fn ks(name: &str, mol: &Molecule, bs: &BasisSet, xc: Xc) {
    let m = build_mol(mol, bs);
    let s = one_electron(&m, cint1e_ovlp_sph);
    let t = one_electron(&m, cint1e_kin_sph);
    let v = one_electron(&m, cint1e_nuc_sph);
    let nbf = s.len();
    let mut hcore = vec![vec![0.0; nbf]; nbf];
    for i in 0..nbf {
        for j in 0..nbf {
            hcore[i][j] = t[i][j] + v[i][j];
        }
    }
    let eri = two_electron(&m);
    let e_nuc = m.nuclear_repulsion();
    let nocc = (mol.n_electrons() / 2) as usize;
    let exx = xc.exx();

    let shells = build_ao_shells(mol, bs);
    let grid = molecular_grid(mol, 96);
    let w: Vec<f64> = grid.iter().map(|g| g.w).collect();
    let ao_grid: Vec<Vec<f64>> = grid.iter().map(|g| eval_ao(&shells, &g.r)).collect();
    let dao_grid: Vec<Vec<[f64; 3]>> = if xc.is_gga() {
        grid.iter().map(|g| qc_core::aoeval::eval_ao_grad(&shells, &g.r).1).collect()
    } else {
        Vec::new()
    };

    let x = s_inv_sqrt(&s);
    let mut p = vec![vec![0.0; nbf]; nbf];
    let mut e_old = 0.0;
    let mut iters = 0;
    let mix = 0.5;
    for _ in 0..300 {
        iters += 1;
        let j = coulomb(&p, &eri);
        let (e_xc, vxc) = xc_eval(&p, &ao_grid, &dao_grid, &w, nbf, xc);
        let k = if exx > 0.0 { exchange(&p, &eri) } else { Vec::new() };

        let mut f = vec![vec![0.0; nbf]; nbf];
        let mut e1 = 0.0;
        let mut ej = 0.0;
        let mut ek = 0.0;
        for u in 0..nbf {
            for vv in 0..nbf {
                let mut fuv = hcore[u][vv] + j[u][vv] + vxc[u][vv];
                if exx > 0.0 {
                    fuv -= exx * 0.5 * k[u][vv];
                    ek += 0.25 * p[u][vv] * k[u][vv];
                }
                f[u][vv] = fuv;
                e1 += p[u][vv] * hcore[u][vv];
                ej += 0.5 * p[u][vv] * j[u][vv];
            }
        }
        let e = e1 + ej - exx * ek + e_xc + e_nuc;
        let (_eps, c) = diagonalize_fock(&f, &x);
        let p_new = build_density(&c, nocc);
        for u in 0..nbf {
            for vv in 0..nbf {
                p[u][vv] = (1.0 - mix) * p[u][vv] + mix * p_new[u][vv];
            }
        }
        if iters > 1 && (e - e_old).abs() < 1e-8 {
            e_old = e;
            break;
        }
        e_old = e;
    }
    let rho = density_on_grid(&p, &ao_grid);
    let n_elec: f64 = rho.iter().zip(&w).map(|(r, wi)| r * wi).sum();
    println!(
        "{:<6} {:<6} nbf={}  iters={}  E={:.6} Ha   ∫ρ={:.4} (N={})",
        name, xc.label(), nbf, iters, e_old, n_elec, mol.n_electrons()
    );
}

/// Dipole integral matrices (x, y, z) in the AO basis.
fn dipole_ao(mol: &Mol) -> [Vec<Vec<f64>>; 3] {
    let (off, nbf) = mol.ao_offsets();
    let mut r = [
        vec![vec![0.0; nbf]; nbf],
        vec![vec![0.0; nbf]; nbf],
        vec![vec![0.0; nbf]; nbf],
    ];
    for ish in 0..mol.nbas as usize {
        let di = mol.shell_nf(ish);
        for jsh in 0..mol.nbas as usize {
            let dj = mol.shell_nf(jsh);
            let shls = [ish as c_int, jsh as c_int];
            let mut buf = vec![0.0_f64; 3 * di * dj];
            unsafe {
                cint1e_r_sph(buf.as_mut_ptr(), shls.as_ptr(), mol.atm.as_ptr(), mol.natm, mol.bas.as_ptr(), mol.nbas, mol.env.as_ptr());
            }
            for c in 0..3 {
                for i in 0..di {
                    for j in 0..dj {
                        r[c][off[ish] + i][off[jsh] + j] = buf[i + di * (j + dj * c)];
                    }
                }
            }
        }
    }
    r
}

/// Oscillator strength of the lowest singlet: f = (2/3) ω |μ|²,
/// μ_c = √2 Σ_ia X_ia ⟨i|r_c|a⟩.
fn oscillator_strength(mo_c: &[Vec<f64>], r: &[Vec<Vec<f64>>; 3], x: &[f64], nocc: usize, nvir: usize, omega: f64) -> f64 {
    let n = mo_c.len();
    let mut mu = [0.0_f64; 3];
    for (c, mu_c) in mu.iter_mut().enumerate() {
        let mut acc = 0.0;
        for i in 0..nocc {
            for a in 0..nvir {
                let av = nocc + a;
                let mut d = 0.0;
                for muu in 0..n {
                    for nu in 0..n {
                        d += mo_c[muu][i] * mo_c[nu][av] * r[c][muu][nu];
                    }
                }
                acc += x[i * nvir + a] * d;
            }
        }
        *mu_c = 2.0_f64.sqrt() * acc;
    }
    let mu2 = mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2];
    (2.0 / 3.0) * omega * mu2
}

/// CIS (TDA-HF) excited states: S1, T1, ΔEST, oscillator strength, emission λ.
fn cis_run(name: &str, mol: &Molecule, bs: &BasisSet) {
    let m = build_mol(mol, bs);
    let s = one_electron(&m, cint1e_ovlp_sph);
    let t = one_electron(&m, cint1e_kin_sph);
    let v = one_electron(&m, cint1e_nuc_sph);
    let nbf = s.len();
    let mut hcore = vec![vec![0.0; nbf]; nbf];
    for i in 0..nbf {
        for j in 0..nbf {
            hcore[i][j] = t[i][j] + v[i][j];
        }
    }
    let eri = two_electron(&m);
    let nocc = (mol.n_electrons() / 2) as usize;
    let scf = rhf(&s, &hcore, &eri, nocc, m.nuclear_repulsion());
    let r = cis(&scf.mo_coeff, &scf.mo_energy, &eri, nocc);
    let dip = dipole_ao(&m);
    let nvir = nbf - nocc;
    let f = oscillator_strength(&scf.mo_coeff, &dip, &r.s1_vec, nocc, nvir, r.s1);
    let lam = 1239.841984 / (r.s1 * HARTREE_EV);
    println!(
        "{:<6} CIS   S1={:.4} eV   T1={:.4} eV   ΔEST={:.4} eV   f={:.4}   λ={:.1} nm",
        name,
        r.s1 * HARTREE_EV,
        r.t1 * HARTREE_EV,
        r.delta_est() * HARTREE_EV,
        f,
        lam
    );
}

/// LDA XC kernel on the grid: singlet v_s = f_uu+f_ud, triplet v_t = f_uu−f_ud,
/// summed over X (LDA_X) and C (LDA_C_VWN).
fn lda_fxc(rho: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let ng = rho.len();
    let mut rho_pol = vec![0.0; 2 * ng];
    for g in 0..ng {
        rho_pol[2 * g] = rho[g] / 2.0;
        rho_pol[2 * g + 1] = rho[g] / 2.0;
    }
    let mut vs = vec![0.0; ng];
    let mut vt = vec![0.0; ng];
    for id in [1, 7] {
        let mut v2 = vec![0.0; 3 * ng];
        unsafe {
            let f = qc_xc_init_pol(id);
            qc_xc_lda_fxc(f, ng as c_int, rho_pol.as_ptr(), v2.as_mut_ptr());
            qc_xc_free(f);
        }
        for g in 0..ng {
            let fuu = v2[3 * g];
            let fud = v2[3 * g + 1];
            vs[g] += fuu + fud;
            vt[g] += fuu - fud;
        }
    }
    (vs, vt)
}

/// Build the singlet/triplet XC-kernel blocks (ia|f_xc|jb) on the grid.
fn build_kernel(
    mo_c: &[Vec<f64>],
    ao_grid: &[Vec<f64>],
    w: &[f64],
    vs: &[f64],
    vt: &[f64],
    nocc: usize,
    nvir: usize,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = mo_c.len();
    let dim = nocc * nvir;
    let mut ks = vec![vec![0.0; dim]; dim];
    let mut kt = vec![vec![0.0; dim]; dim];
    for g in 0..ao_grid.len() {
        let ao = &ao_grid[g];
        // MO values at this point: φ_p = Σ_μ C_μp ao_μ
        let mut phi = vec![0.0; n];
        for p in 0..n {
            let mut s = 0.0;
            for mu in 0..n {
                s += mo_c[mu][p] * ao[mu];
            }
            phi[p] = s;
        }
        // occ-virt products ρ_ia = φ_i φ_a
        let mut pia = vec![0.0; dim];
        for i in 0..nocc {
            for a in 0..nvir {
                pia[i * nvir + a] = phi[i] * phi[nocc + a];
            }
        }
        let ws = w[g] * vs[g];
        let wt = w[g] * vt[g];
        for x in 0..dim {
            let px = pia[x];
            if px != 0.0 {
                let (as_, at_) = (ws * px, wt * px);
                for y in 0..dim {
                    ks[x][y] += as_ * pia[y];
                    kt[x][y] += at_ * pia[y];
                }
            }
        }
    }
    (ks, kt)
}

/// TDA-DFT (LDA/SVWN5 kernel) excited states: S1, T1, ΔEST.
fn tda_dft_lda(name: &str, mol: &Molecule, bs: &BasisSet) {
    let m = build_mol(mol, bs);
    let s = one_electron(&m, cint1e_ovlp_sph);
    let t = one_electron(&m, cint1e_kin_sph);
    let v = one_electron(&m, cint1e_nuc_sph);
    let nbf = s.len();
    let mut hcore = vec![vec![0.0; nbf]; nbf];
    for i in 0..nbf {
        for j in 0..nbf {
            hcore[i][j] = t[i][j] + v[i][j];
        }
    }
    let eri = two_electron(&m);
    let e_nuc = m.nuclear_repulsion();
    let nocc = (mol.n_electrons() / 2) as usize;
    let nvir = nbf - nocc;

    let shells = build_ao_shells(mol, bs);
    let grid = molecular_grid(mol, 96);
    let w: Vec<f64> = grid.iter().map(|g| g.w).collect();
    let ao_grid: Vec<Vec<f64>> = grid.iter().map(|g| eval_ao(&shells, &g.r)).collect();
    let no_dao: Vec<Vec<[f64; 3]>> = Vec::new();

    // LDA SCF, capturing converged MOs.
    let x = s_inv_sqrt(&s);
    let mut p = vec![vec![0.0; nbf]; nbf];
    let mut mo_c = vec![vec![0.0; nbf]; nbf];
    let mut mo_e = vec![0.0; nbf];
    let mut e_old = 0.0;
    for it in 0..300 {
        let j = coulomb(&p, &eri);
        let (e_xc, vxc) = xc_eval(&p, &ao_grid, &no_dao, &w, nbf, Xc::Lda);
        let mut f = vec![vec![0.0; nbf]; nbf];
        let mut e1 = 0.0;
        let mut ej = 0.0;
        for u in 0..nbf {
            for vv in 0..nbf {
                f[u][vv] = hcore[u][vv] + j[u][vv] + vxc[u][vv];
                e1 += p[u][vv] * hcore[u][vv];
                ej += 0.5 * p[u][vv] * j[u][vv];
            }
        }
        let e = e1 + ej + e_xc + e_nuc;
        let (eps, c) = diagonalize_fock(&f, &x);
        mo_e = eps;
        mo_c = c;
        let p_new = build_density(&mo_c, nocc);
        for u in 0..nbf {
            for vv in 0..nbf {
                p[u][vv] = 0.5 * p[u][vv] + 0.5 * p_new[u][vv];
            }
        }
        if it > 1 && (e - e_old).abs() < 1e-9 {
            e_old = e;
            break;
        }
        e_old = e;
    }

    // XC kernel on the final density, then TDA.
    let rho = density_on_grid(&p, &ao_grid);
    let (vs, vt) = lda_fxc(&rho);
    let (ks, kt) = build_kernel(&mo_c, &ao_grid, &w, &vs, &vt, nocc, nvir);
    let r = tda_general(&mo_c, &mo_e, &eri, nocc, 0.0, &ks, &kt);
    println!(
        "{:<6} TDA-SVWN5  S1={:.4} eV  T1={:.4} eV  ΔEST={:.4} eV",
        name,
        r.s1 * HARTREE_EV,
        r.t1 * HARTREE_EV,
        r.delta_est() * HARTREE_EV
    );
}

/// ERI tensor with the range-separated erf(ωr)/r operator (ω=0 → full 1/r).
fn eri_with_omega(mol: &Mol, omega: f64) -> Eri {
    let mm = Mol {
        atm: mol.atm.clone(),
        bas: mol.bas.clone(),
        env: {
            let mut e = mol.env.clone();
            e[8] = omega; // PTR_RANGE_OMEGA
            e
        },
        natm: mol.natm,
        nbas: mol.nbas,
    };
    two_electron(&mm)
}

/// Range-separated hybrid (CAM-B3LYP) ground-state SCF.
/// Exact exchange = α·K_full + β·K_longrange, with (ω,α,β) from libxc.
fn rsh_scf(name: &str, mol: &Molecule, bs: &BasisSet) {
    let m = build_mol(mol, bs);
    let s = one_electron(&m, cint1e_ovlp_sph);
    let t = one_electron(&m, cint1e_kin_sph);
    let v = one_electron(&m, cint1e_nuc_sph);
    let nbf = s.len();
    let mut hcore = vec![vec![0.0; nbf]; nbf];
    for i in 0..nbf {
        for j in 0..nbf {
            hcore[i][j] = t[i][j] + v[i][j];
        }
    }
    let eri = two_electron(&m);
    let e_nuc = m.nuclear_repulsion();
    let nocc = (mol.n_electrons() / 2) as usize;
    let xc = Xc::CamB3lyp;

    let (mut omega, mut alpha, mut beta) = (0.0_f64, 0.0_f64, 0.0_f64);
    unsafe {
        let f = qc_xc_init(xc.func_ids()[0]);
        qc_xc_cam_coef(f, &mut omega, &mut alpha, &mut beta);
        qc_xc_free(f);
    }
    let eri_lr = eri_with_omega(&m, omega);

    let shells = build_ao_shells(mol, bs);
    let grid = molecular_grid(mol, 96);
    let w: Vec<f64> = grid.iter().map(|g| g.w).collect();
    let ao_grid: Vec<Vec<f64>> = grid.iter().map(|g| eval_ao(&shells, &g.r)).collect();
    let dao_grid: Vec<Vec<[f64; 3]>> =
        grid.iter().map(|g| qc_core::aoeval::eval_ao_grad(&shells, &g.r).1).collect();

    let x = s_inv_sqrt(&s);
    let mut p = vec![vec![0.0; nbf]; nbf];
    let mut e_old = 0.0;
    let mut iters = 0;
    for _ in 0..300 {
        iters += 1;
        let j = coulomb(&p, &eri);
        let (e_xc, vxc) = xc_eval(&p, &ao_grid, &dao_grid, &w, nbf, xc);
        let kf = exchange(&p, &eri);
        let kl = exchange(&p, &eri_lr);
        let mut f = vec![vec![0.0; nbf]; nbf];
        let mut e1 = 0.0;
        let mut ej = 0.0;
        let mut ek = 0.0;
        for u in 0..nbf {
            for vv in 0..nbf {
                // libxc: EXX = α/r + β·erfc(ωr)/r; erfc = 1/r − erf, so with
                // K_LR using erf(ωr): K = (α+β)·K_full − β·K_LR.
                let kx = (alpha + beta) * kf[u][vv] - beta * kl[u][vv];
                f[u][vv] = hcore[u][vv] + j[u][vv] + vxc[u][vv] - 0.5 * kx;
                e1 += p[u][vv] * hcore[u][vv];
                ej += 0.5 * p[u][vv] * j[u][vv];
                ek += 0.25 * p[u][vv] * kx;
            }
        }
        let e = e1 + ej - ek + e_xc + e_nuc;
        let (_eps, c) = diagonalize_fock(&f, &x);
        let p_new = build_density(&c, nocc);
        for u in 0..nbf {
            for vv in 0..nbf {
                p[u][vv] = 0.5 * p[u][vv] + 0.5 * p_new[u][vv];
            }
        }
        if iters > 1 && (e - e_old).abs() < 1e-9 {
            e_old = e;
            break;
        }
        e_old = e;
    }
    println!(
        "{:<6} {:<10} iters={}  E={:.6} Ha   (ω={:.3} α={:.3} β={:.3})",
        name, xc.label(), iters, e_old, omega, alpha, beta
    );
}

fn main() {
    // Emitter mode: `cint-sys <molecule.xyz> <basis.nw>` -> direct-SCF PBE0 on it.
    let args: Vec<String> = std::env::args().collect();
    if args.len() >= 3 {
        let xyz = std::fs::read_to_string(&args[1]).expect("read xyz");
        let mol = Molecule::from_xyz_str(&xyz).expect("parse xyz");
        let basis_txt = std::fs::read_to_string(&args[2]).expect("read basis");
        let bs = BasisSet::parse_nwchem(&basis_txt).expect("parse basis");
        let name: String = std::path::Path::new(&args[1])
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("mol")
            .to_string();
        println!(
            "EMİTÖR: {} — {} atom, {} elektron",
            name,
            mol.atoms.len(),
            mol.n_electrons()
        );
        ks_direct(&name, &mol, &bs, Xc::Pbe0);
        return;
    }

    let bs = BasisSet::parse_nwchem(STO3G).expect("parse STO-3G");
    println!("=== RHF/STO-3G via libcint + owned SCF ===");

    // H2 at R = 1.4 bohr
    let h2 = Molecule { atoms: vec![
        Atom { z: 1, coord: [0.0, 0.0, 0.0] },
        Atom { z: 1, coord: [0.0, 0.0, 1.4] },
    ]};
    run("H2", &h2, &bs, -1.1167, 1e-3);

    // H2O: R(OH)=0.9584 Å, angle 104.45°  (coords in bohr)
    let h2o = Molecule { atoms: vec![
        Atom { z: 8, coord: [0.0, 0.0, 0.0] },
        Atom { z: 1, coord: [1.43200, 0.0, 1.10891] },
        Atom { z: 1, coord: [-1.43200, 0.0, 1.10891] },
    ]};
    run("H2O", &h2o, &bs, -74.9659, 2e-2);

    xc_smoke();

    println!("\n=== grid + AO evaluation vs analytic overlap ===");
    validate_grid_ao("H2", &h2, &bs);
    validate_grid_ao("H2O", &h2o, &bs);

    println!("\n=== KS-DFT (owned engine: LDA / PBE / PBE0) ===");
    for xc in [Xc::Lda, Xc::Pbe, Xc::Pbe0] {
        ks("H2", &h2, &bs, xc);
    }
    for xc in [Xc::Lda, Xc::Pbe, Xc::Pbe0] {
        ks("H2O", &h2o, &bs, xc);
    }

    println!("\n=== Excited states: CIS (TDA-HF) → S1, T1, ΔEST ===");
    cis_run("H2", &h2, &bs);
    cis_run("H2O", &h2o, &bs);

    println!("\n=== Excited states: TDA-DFT (SVWN5 kernel) ===");
    tda_dft_lda("H2", &h2, &bs);
    tda_dft_lda("H2O", &h2o, &bs);

    println!("\n=== Range-separated hybrid (CAM-B3LYP) ground state ===");
    rsh_scf("H2", &h2, &bs);
    rsh_scf("H2O", &h2o, &bs);

    println!("\n=== def2-SVP (d functions) — water ===");
    match std::fs::read_to_string("../../data/def2-svp.nw") {
        Ok(txt) => {
            let bs2 = BasisSet::parse_nwchem(&txt).expect("parse def2-svp");
            validate_grid_ao("H2O", &h2o, &bs2);
            ks("H2O", &h2o, &bs2, Xc::Pbe0);
            ks_direct("H2O", &h2o, &bs2, Xc::Pbe0); // must equal stored-ERI result
        }
        Err(e) => println!("  could not read ../../data/def2-svp.nw: {e}"),
    }
}

/// STO-3G basis (NWChem format) for H and O.
const STO3G: &str = r#"
H   S
   3.42525091   0.15432897
   0.62391373   0.53532814
   0.16885540   0.44463454
O   S
 130.70932000   0.15432897
  23.80886100   0.53532814
   6.44360830   0.44463454
O   SP
   5.03315130  -0.09996723   0.15591628
   1.16959610   0.39951283   0.60768372
   0.38038900   0.70011547   0.39195739
"#;
