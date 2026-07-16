"""PySCF validation oracle: same geometry, basis, and libxc functionals as the
owned engine, so the numbers are directly comparable."""
from pyscf import gto, dft


def run(atom):
    mol = gto.M(atom=atom, basis="sto-3g", unit="Bohr", verbose=0)
    out = {}
    for label, xc in [("SVWN5", "LDA_X,LDA_C_VWN"), ("PBE", "GGA_X_PBE,GGA_C_PBE"), ("PBE0", "PBE0")]:
        mf = dft.RKS(mol)
        mf.xc = xc
        out[label] = mf.kernel()
    return out


def cis(atom):
    from pyscf import tdscf
    au = 27.211386245988
    mol = gto.M(atom=atom, basis="sto-3g", unit="Bohr", verbose=0)
    mf = mol.RHF().run()
    ts = tdscf.TDA(mf); ts.singlet = True; ts.nstates = 1; ts.kernel()
    tt = tdscf.TDA(mf); tt.singlet = False; tt.nstates = 1; tt.kernel()
    s1 = ts.e[0] * au
    t1 = tt.e[0] * au
    f = ts.oscillator_strength()[0]
    lam = 1239.841984 / s1
    return s1, t1, s1 - t1, f, lam


for name, atom in [("H2", "H 0 0 0; H 0 0 1.4"),
                   ("H2O", "O 0 0 0; H 1.432 0 1.10891; H -1.432 0 1.10891")]:
    r = run(atom)
    print(f"{name:<4} SVWN5={r['SVWN5']:.6f}  PBE={r['PBE']:.6f}  PBE0={r['PBE0']:.6f}")

print("--- CIS (TDA-HF) ---")
for name, atom in [("H2", "H 0 0 0; H 0 0 1.4"),
                   ("H2O", "O 0 0 0; H 1.432 0 1.10891; H -1.432 0 1.10891")]:
    s1, t1, d, f, lam = cis(atom)
    print(f"{name:<4} S1={s1:.4f} eV  T1={t1:.4f} eV  dEST={d:.4f} eV  f={f:.4f}  lam={lam:.1f} nm")


def tda_dft(atom, xc):
    from pyscf import tdscf
    au = 27.211386245988
    mol = gto.M(atom=atom, basis="sto-3g", unit="Bohr", verbose=0)
    mf = dft.RKS(mol); mf.xc = xc; mf.kernel()
    ts = mf.TDA(); ts.singlet = True; ts.nstates = 1; ts.kernel()
    tt = mf.TDA(); tt.singlet = False; tt.nstates = 1; tt.kernel()
    s1 = ts.e[0] * au; t1 = tt.e[0] * au
    return s1, t1, s1 - t1


print("--- TDA-DFT (SVWN5 kernel) ---")
for name, atom in [("H2", "H 0 0 0; H 0 0 1.4"),
                   ("H2O", "O 0 0 0; H 1.432 0 1.10891; H -1.432 0 1.10891")]:
    s1, t1, d = tda_dft(atom, "LDA_X,LDA_C_VWN")
    print(f"{name:<4} S1={s1:.4f} eV  T1={t1:.4f} eV  dEST={d:.4f} eV")


def rsh(atom):
    mol = gto.M(atom=atom, basis="sto-3g", unit="Bohr", verbose=0)
    mf = dft.RKS(mol); mf.xc = "CAM_B3LYP"
    return mf.kernel()


print("--- CAM-B3LYP ground state ---")
for name, atom in [("H2", "H 0 0 0; H 0 0 1.4"),
                   ("H2O", "O 0 0 0; H 1.432 0 1.10891; H -1.432 0 1.10891")]:
    print(f"{name:<4} E={rsh(atom):.6f} Ha")

print("--- def2-SVP water PBE0 ---")
mol = gto.M(atom="O 0 0 0; H 1.432 0 1.10891; H -1.432 0 1.10891",
            basis="def2svp", unit="Bohr", verbose=0)
mf = dft.RKS(mol); mf.xc = "PBE0"
print(f"H2O  PBE0/def2-SVP  nbf={mol.nao}  E={mf.kernel():.6f} Ha")
