//! qc-core Stage 0 driver: RHF/STO-3G on H2, validated against the textbook
//! reference (Szabo & Ostlund: -1.1167 Hartree at R = 1.4 bohr).

use qc_core::integrals::{build_matrices, h_sto3g, nuclear_repulsion};
use qc_core::scf::rhf;

fn main() {
    let r = 1.4;
    let bfs = vec![h_sto3g([0.0, 0.0, 0.0]), h_sto3g([0.0, 0.0, r])];
    let nuclei = vec![([0.0, 0.0, 0.0], 1.0), ([0.0, 0.0, r], 1.0)];

    let (s, hcore, eri) = build_matrices(&bfs, &nuclei);
    let scf = rhf(&s, &hcore, &eri, 1, nuclear_repulsion(&nuclei));
    let reference = -1.1167;

    println!("=== qc-core Stage 0: RHF/STO-3G H2 (R = {r} bohr) ===");
    println!("  SCF iterations : {}", scf.iters);
    println!("  E_total        : {:.6} Hartree", scf.e_total);
    println!("  reference      : {reference:.6} Hartree (Szabo & Ostlund)");
    println!("  Δ              : {:.2e} Hartree", (scf.e_total - reference).abs());
    println!("  MO energies    : {:?}", scf.mo_energy.iter().map(|e| (e * 1e4).round() / 1e4).collect::<Vec<_>>());
    if (scf.e_total - reference).abs() < 1e-3 {
        println!("  ✅ MATCH — owned SCF core validated.");
    } else {
        println!("  ❌ mismatch.");
    }
}

#[cfg(test)]
mod tests {
    use qc_core::integrals::{build_matrices, h_sto3g, nuclear_repulsion};
    use qc_core::scf::rhf;

    #[test]
    fn h2_sto3g_energy() {
        let r = 1.4;
        let bfs = vec![h_sto3g([0.0, 0.0, 0.0]), h_sto3g([0.0, 0.0, r])];
        let nuclei = vec![([0.0, 0.0, 0.0], 1.0), ([0.0, 0.0, r], 1.0)];
        let (s, hcore, eri) = build_matrices(&bfs, &nuclei);
        let scf = rhf(&s, &hcore, &eri, 1, nuclear_repulsion(&nuclei));
        assert!((scf.e_total - (-1.1167)).abs() < 1e-3, "got {}", scf.e_total);
    }
}
