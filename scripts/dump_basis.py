"""Dump a basis set (default def2-SVP) in NWChem format for the TADF elements,
so the Rust parser reads identical data (no transcription errors). General
contractions are split into separate single-contraction shells."""
import sys
from pyscf import gto

name = sys.argv[1] if len(sys.argv) > 1 else "def2svp"
out = sys.argv[2] if len(sys.argv) > 2 else "data/def2-svp.nw"
elements = ["H", "B", "C", "N", "O", "F", "S", "Si"]
angl = {0: "S", 1: "P", 2: "D", 3: "F", 4: "G"}

lines = [f"# {name} (dumped from PySCF)"]
for el in elements:
    try:
        b = gto.basis.load(name, el)
    except Exception as e:
        print(f"skip {el}: {e}", file=sys.stderr)
        continue
    for shell in b:
        l = shell[0]
        prims = shell[1:]  # rows: [exp, c1, c2, ...]
        ncontr = len(prims[0]) - 1
        for ci in range(ncontr):
            lines.append(f"{el}   {angl[l]}")
            for p in prims:
                lines.append(f"  {p[0]:.10f}   {p[1 + ci]:.10f}")

with open(out, "w") as fh:
    fh.write("\n".join(lines) + "\n")
print(f"wrote {out} ({len(lines)} lines)")
