#!/usr/bin/env python3
"""Analyze Boltz prediction CIF(s) against the constraints defined in a Boltz input JSON.

Usage:
  python3 analyze_structure.py --input <boltz_input.json> [--cif DIR_OR_FILE] [--bond "B:7:C,B:1:N"]

- Constraints are read from the input JSON (type=contact, 0-based residue_index,
  max_distance_angstrom) and checked as min heavy-atom distance in the CIF (1-based).
- Also reports: chain composition, per-chain mean pLDDT, pairwise interfaces (<4 A),
  optional cyclization / custom bond check.
- If --cif points at a directory, every sample_*.cif plus the best-sample
  *_predicted.cif is analyzed; distances are aggregated (median/range).

Pure standard library; works without numpy/gemmi/BioPython.
"""
import sys, os, json, math, glob, argparse, statistics, collections

def read_cif_atoms(path):
    lines = open(path).read().splitlines()
    atoms, cols, i = [], [], 0
    while i < len(lines):
        if lines[i].startswith("_atom_site."):
            while i < len(lines) and lines[i].startswith("_atom_site."):
                cols.append(lines[i].strip().split(".", 1)[1]); i += 1
            idx = {c: k for k, c in enumerate(cols)}
            while i < len(lines):
                ln = lines[i]
                if ln.startswith(("#", "loop_", "_")):
                    break
                if ln.strip():
                    p = ln.split()
                    if p[0] in ("ATOM", "HETATM"):
                        atoms.append(dict(
                            chain=p[idx["label_asym_id"]],
                            resname=p[idx["label_comp_id"]],
                            resid=int(p[idx["label_seq_id"]]),
                            atom=p[idx["label_atom_id"]].strip('"'),
                            element=p[idx["type_symbol"]],
                            x=float(p[idx["Cartn_x"]]), y=float(p[idx["Cartn_y"]]),
                            z=float(p[idx["Cartn_z"]]), b=float(p[idx["B_iso_or_equiv"]])))
                i += 1
            break
        i += 1
    return atoms

def res_atoms(atoms, chain, resid, heavy=True):
    return [a for a in atoms if a["chain"] == chain and a["resid"] == resid
            and (a["element"] != "H" or not heavy)]

def min_dist(a1, a2):
    return min(math.dist((p["x"], p["y"], p["z"]), (q["x"], q["y"], q["z"]))
               for p in a1 for q in a2)

def interface(atoms, c1, c2, cutoff=4.0):
    A = [a for a in atoms if a["chain"] == c1 and a["element"] != "H"]
    B = [a for a in atoms if a["chain"] == c2 and a["element"] != "H"]
    c2s = cutoff ** 2
    s1, s2 = collections.Counter(), collections.Counter()
    for p in A:
        for q in B:
            if (p["x"]-q["x"])**2 + (p["y"]-q["y"])**2 + (p["z"]-q["z"])**2 <= c2s:
                s1[(p["resid"], p["resname"])] += 1
                s2[(q["resid"], q["resname"])] += 1
    return s1, s2

def fmt(counter):
    return " ".join(f"{r}{i}" for (i, r) in sorted(counter))

def load_constraints(input_json):
    """Return list of (label, chain1, res1_1based, chain2, res2_1based, limit)."""
    d = json.load(open(input_json))
    out = []
    for c in d.get("constraints", []):
        t1, t2 = c["token1"], c["token2"]
        out.append((f"{t1['chain_id']}{t1['residue_index']+1}-{t2['chain_id']}{t2['residue_index']+1}",
                    t1["chain_id"], t1["residue_index"] + 1,
                    t2["chain_id"], t2["residue_index"] + 1,
                    float(c["max_distance_angstrom"])))
    return out

def analyze_cif(path, cons, bond_spec):
    atoms = read_cif_atoms(path)
    chains = collections.Counter(a["chain"] for a in atoms)
    print(f"\n=== {os.path.basename(path)} ===")
    print("chains:", {c: f"{len({a['resid'] for a in atoms if a['chain']==c})} res / {n} atoms"
                      for c, n in sorted(chains.items())})
    for c in sorted(chains):
        ca = [a for a in atoms if a["chain"] == c and a["atom"] == "CA"]
        if ca:
            pl = sum(a["b"] for a in ca) / len(ca)
            if pl <= 1.5:
                pl *= 100
            print(f"  chain {c}: mean pLDDT(CA) = {pl:.1f}")

    if cons:
        print("[constraints] (0-based input index -> 1-based residue)")
        ok = 0
        for label, c1, r1, c2, r2, lim in cons:
            a1, a2 = res_atoms(atoms, c1, r1), res_atoms(atoms, c2, r2)
            if not a1 or not a2:
                print(f"  {label:28s} MISSING residue")
                continue
            dmin = min_dist(a1, a2)
            good = dmin <= lim
            ok += good
            print(f"  {label:28s} min heavy-atom {dmin:6.2f} A  limit {lim:.0f} A  "
                  f"{'OK' if good else 'VIOLATION'}")
        print(f"  satisfied {ok}/{len(cons)}")

    if bond_spec:
        print(f"[bond {bond_spec}]")
        a, b = bond_spec.split(",")
        (ca, ri, at), (cb, rj, bt) = a.split(":"), b.split(":")
        p = [x for x in atoms if x["chain"] == ca and x["resid"] == int(ri) and x["atom"] == at]
        q = [x for x in atoms if x["chain"] == cb and x["resid"] == int(rj) and x["atom"] == bt]
        if p and q:
            d = math.dist((p[0]["x"], p[0]["y"], p[0]["z"]), (q[0]["x"], q[0]["y"], q[0]["z"]))
            print(f"  {at}{ri}({ca}) - {bt}{rj}({cb}) = {d:.2f} A  "
                  f"({'bonded' if d < 1.8 else 'NOT bonded'})")

    print("[interfaces <4 A]")
    for c1, c2 in sorted({(x, y) for x in chains for y in chains if x < y}):
        s1, s2 = interface(atoms, c1, c2)
        print(f"  {c1}-{c2}: {len(s1)} res / {len(s2)} res")
        print(f"    {c1}: {fmt(s1)}")
        print(f"    {c2}: {fmt(s2)}")
    return atoms

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Boltz input JSON (constraints source)")
    ap.add_argument("--cif", required=True, help="CIF file, or prediction dir")
    ap.add_argument("--bond", default=None, help='e.g. "B:7:C,B:1:N" (chain:res:atom pairs)')
    args = ap.parse_args()

    cons = load_constraints(args.input)
    if os.path.isdir(args.cif):
        files = sorted(glob.glob(os.path.join(args.cif, "sample_*_predicted_structure.cif")))
        best = [f for f in glob.glob(os.path.join(args.cif, "*_predicted.cif"))]
        if best:
            files = best + files
        if not files:
            sys.exit(f"no CIF found in {args.cif}")
        rows = []
        for f in files:
            tag = "best" if f in best else os.path.basename(f).split("_")[1]
            atoms = analyze_cif(f, cons, args.bond)
            if cons:
                ds = [min_dist(res_atoms(atoms, c1, r1), res_atoms(atoms, c2, r2))
                      for _, c1, r1, c2, r2, _ in cons]
                rows.append((tag, ds))
        if rows and cons:
            print("\n[aggregate over samples]")
            for i, (_, c1, r1, c2, r2, lim) in enumerate(cons):
                col = [r[1][i] for r in rows]
                sat = sum(1 for d in col if d <= lim)
                print(f"  {c1}{r1}-{c2}{r2}: median {statistics.median(col):5.2f}  "
                      f"range {min(col):.2f}-{max(col):.2f}  limit {lim:.0f}  满足 {sat}/{len(col)}")
    else:
        analyze_cif(args.cif, cons, args.bond)

if __name__ == "__main__":
    main()
