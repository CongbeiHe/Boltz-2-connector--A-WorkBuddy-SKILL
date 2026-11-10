#!/usr/bin/env python3
"""Convert a Boltz predicted mmCIF into a plain PDB (B-factor = pLDDT*100).

Usage: python3 cif2pdb.py <input.cif> <output.pdb>
Handles Boltz pLDDT stored as 0-1 in B_iso column; keeps chain IDs as-is.
Compatible with PyMOL/ChimeraX. No third-party dependencies.
"""
import sys

def read_atoms(path):
    rows, cols, i = [], [], 0
    lines = open(path).read().splitlines()
    while i < len(lines):
        if lines[i].startswith("_atom_site."):
            while i < len(lines) and lines[i].startswith("_atom_site."):
                cols.append(lines[i].strip().split(".", 1)[1])
                i += 1
            idx = {c: k for k, c in enumerate(cols)}
            while i < len(lines):
                ln = lines[i]
                if ln.startswith(("#", "loop_", "_")):
                    break
                if ln.strip():
                    p = ln.split()
                    if p[0] in ("ATOM", "HETATM"):
                        rows.append((p, idx))
                i += 1
            break
        i += 1
    return rows

def main(src, dst):
    rows = read_atoms(src)
    out = []
    for p, idx in rows:
        name = p[idx["label_atom_id"]].strip('"')
        resn = p[idx["label_comp_id"]].strip('"')
        chain = p[idx["label_asym_id"]]
        resi = p[idx["label_seq_id"]]
        if resi in (".", "?"):
            resi = "1"
        x, y, z = (float(p[idx["Cartn_x"]]), float(p[idx["Cartn_y"]]),
                   float(p[idx["Cartn_z"]]))
        b = float(p[idx["B_iso_or_equiv"]])
        if b <= 1.5:      # Boltz stores 0-1 pLDDT
            b *= 100.0
        elem = p[idx["type_symbol"]]
        aname = name if len(name) == 4 else f" {name:<3s}"
        out.append(
            f"ATOM  {len(out)+1:5d} {aname:<4s}{' ' if len(resn)<4 else ''}{resn:>3s} {chain}"
            f"{int(float(resi)):4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{b:6.2f}"
            f"          {elem:>2s}"
        )
    out.append("END")
    open(dst, "w").write("\n".join(out) + "\n")
    print(f"wrote {dst} ({len(out)-1} atoms)")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
