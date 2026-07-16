#!/usr/bin/env python3
"""PySCF physics kernel for the deep-blue TADF screening engine (Seçenek A).

Pure-physics, deterministic. Given ONE molecular geometry, computes the
descriptors the screening engine ranks on, and emits a JSON record with full
method PROVENANCE. No ML, no fitting — every number comes from an equation.

First-pass model (cheap, defensible — see docs/REFERENCE_SET.md):
  * Ground-state (S0) restricted KS-DFT  -> HOMO, LUMO, gap.
  * VERTICAL ΔEST at the S0 geometry via TDA-TDDFT:
      run TDA for singlets -> E(S1), oscillator strength f(S1)
      run TDA for triplets -> E(T1)
      ΔEST = E(S1) - E(T1);  emission λ ≈ 1239.84 / E(S1[eV])
    The Bonn/Herbert SI shows the vertical approximation (S0 geometry for both
    S1 and T1) only marginally worsens ΔEST vs. relaxed geometries — so this is
    a legitimate, much cheaper first pass. Relaxed-S1 emission is a later step.

Recommended methods (confirmed on STGABS27):
  * functional: PBE0 (MAD ~0.06 eV w/ state-specific PCM) or OT-ωB97M-V (best).
  * basis: def2-SVP is sufficient for these CT states (<0.1 eV vs def2-TZVPP).
  * solvent: optional IEF-PCM (ε for the measurement medium).

MR-TADF (DABNA family): TDA-TDDFT FAILS (multireference). For those, ΔDFT/ΔSCF
is required — flagged here, implemented in a later milestone.

Usage:
  python engine/pyscf_kernel.py --xyz mol.xyz --name MCz-XT \
      --functional pbe0 --basis def2-svp [--solvent-eps 2.38] [--tda/--full]
Outputs one JSON object to stdout.
"""
from __future__ import annotations

import argparse
import json
import sys
import time

HARTREE_EV = 27.211386245988
HC_EV_NM = 1239.841984


def read_xyz(path: str):
    """Read a .xyz file -> PySCF atom string. Deterministic (preserves order)."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{path}: not a valid xyz")
    try:
        natom = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"{path}: first line must be atom count") from exc
    atoms = lines[2 : 2 + natom]
    return "\n".join(atoms), natom


def _scf(atom, args):
    """Ground-state RKS at a geometry. Returns (mf, e_scf_hartree)."""
    from pyscf import dft, gto

    mol = gto.M(atom=atom, basis=args.basis, charge=args.charge, spin=args.spin, verbose=0)
    mf = dft.RKS(mol)
    mf.xc = args.functional
    if args.solvent_eps is not None:
        mf = mf.PCM()
        mf.with_solvent.eps = args.solvent_eps
        mf.with_solvent.method = "IEF-PCM"
    e = mf.kernel()
    return mf, float(e)


def _excitation(mf, singlet, nstates, tda):
    """Lowest excitation energy (eV) of given spin; osc strength for singlets."""
    from pyscf import tddft

    td = tddft.TDA(mf) if tda else tddft.TDDFT(mf)
    td.singlet = singlet
    td.nstates = nstates
    td.kernel()
    e_ev = float(td.e[0]) * HARTREE_EV
    osc = float(td.oscillator_strength()[0]) if singlet else None
    return e_ev, osc


def _frontier(mf):
    mo_e, occ = mf.mo_energy, mf.mo_occ
    h = max(i for i, o in enumerate(occ) if o > 0)
    return mo_e[h] * HARTREE_EV, mo_e[h + 1] * HARTREE_EV


def run(args) -> dict:
    t0 = time.time()
    atom_s1, natom = read_xyz(args.xyz)

    # S1 state: evaluated at the S1-optimized geometry.
    mf_s1, e_scf_s1 = _scf(atom_s1, args)
    exc_s1, f_s1 = _excitation(mf_s1, True, args.nstates, args.tda)
    homo, lumo = _frontier(mf_s1)
    e_total_s1 = e_scf_s1 * HARTREE_EV + exc_s1
    lam_em = HC_EV_NM / exc_s1 if exc_s1 > 0 else None

    if args.xyz_t1:
        # ADIABATIC ΔEST: T1 evaluated at its own optimized geometry.
        atom_t1, _ = read_xyz(args.xyz_t1)
        mf_t1, e_scf_t1 = _scf(atom_t1, args)
        exc_t1, _ = _excitation(mf_t1, False, args.nstates, args.tda)
        e_total_t1 = e_scf_t1 * HARTREE_EV + exc_t1
        delta_est = e_total_s1 - e_total_t1
        approximation = "adiabatic (S1@S1geom - T1@T1geom)"
    else:
        # VERTICAL ΔEST: both states at the single supplied geometry.
        exc_t1, _ = _excitation(mf_s1, False, args.nstates, args.tda)
        delta_est = exc_s1 - exc_t1
        approximation = "vertical@input-geom"

    return {
        "name": args.name,
        "natom": natom,
        "deltaEST_eV": round(delta_est, 4),
        "E_S1_vert_eV": round(exc_s1, 4),
        "E_T1_vert_eV": round(exc_t1, 4),
        "f_osc_S1": (round(f_s1, 5) if f_s1 is not None else None),
        "lambda_em_nm": (round(lam_em, 1) if lam_em else None),
        "HOMO_eV": round(homo, 4),
        "LUMO_eV": round(lumo, 4),
        "HL_gap_eV": round(lumo - homo, 4),
        "provenance": {
            "method": "TDA-TDDFT" if args.tda else "TDDFT",
            "approximation": approximation,
            "functional": args.functional,
            "basis": args.basis,
            "solvent": (
                {"model": "IEF-PCM", "eps": args.solvent_eps}
                if args.solvent_eps is not None
                else None
            ),
            "geom_s1": args.xyz,
            "geom_t1": args.xyz_t1,
            "code": "pyscf",
            "note": "MR-TADF (DABNA) needs ΔDFT; TDA-TDDFT not reliable there",
        },
        "wall_s": round(time.time() - t0, 1),
    }


def main() -> int:
    p = argparse.ArgumentParser(description="PySCF TADF descriptor kernel")
    p.add_argument("--xyz", required=True, help="S1-optimized geometry (or single geom for vertical)")
    p.add_argument("--xyz-t1", default=None, help="T1-optimized geometry; enables adiabatic ΔEST")
    p.add_argument("--name", default="unnamed")
    p.add_argument("--functional", default="pbe0")
    p.add_argument("--basis", default="def2-svp")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--spin", type=int, default=0)
    p.add_argument("--nstates", type=int, default=3)
    p.add_argument("--solvent-eps", type=float, default=None)
    grp = p.add_mutually_exclusive_group()
    grp.add_argument("--tda", dest="tda", action="store_true", default=True)
    grp.add_argument("--full", dest="tda", action="store_false")
    args = p.parse_args()

    try:
        rec = run(args)
    except ImportError as exc:
        print(json.dumps({"name": args.name, "error": f"pyscf not available: {exc}"}))
        return 3
    except Exception as exc:  # noqa: BLE001 - report any kernel failure as JSON
        print(json.dumps({"name": args.name, "error": str(exc)}))
        return 1

    print(json.dumps(rec))
    return 0


if __name__ == "__main__":
    sys.exit(main())
