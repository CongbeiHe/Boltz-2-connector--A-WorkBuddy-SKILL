#!/usr/bin/env python3
"""Compare two interfaces of a binder (e.g. cyclic peptide) against its partners.

Usage:
  python3 compare_ifaces.py --input <boltz_input.json> --binder B \
      --cif <sample.cif> --pae <sample_pae.npz>

Metrics per interface (binder-partnerX):
  - contact heavy-atom pairs <4.5 A
  - interface residue count + mean pLDDT on both sides
  - interface-restricted mean PAE (position certainty)
  - whole-chain mean PAE (binder vs entire partner)

Token order/base is derived automatically from entity lengths in the input JSON
(entities are tokenized one-token-per-residue in declaration order).
PAE npz/npy is parsed with the stdlib (zipfile + struct) — numpy not required.
"""
import sys, os, json, glob, argparse, zipfile, struct, ast, math, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_structure import read_cif_atoms, interface, fmt

# ---------- npy parsing (stdlib only) ----------
def load_npy(fp):
    hdr = fp.read(10)
    assert hdr[:6] == b"\x93NUMPY", "not a .npy stream"
    ver = hdr[6]
    if ver == 1:
        hlen = struct.unpack("<H", hdr[8:10])[0]; off = 10
    else:
        hlen = struct.unpack("<I", fp.read(4))[0]; off = 12
    meta = ast.literal_eval(fp.read(hlen).decode("latin1"))
    descr, shape = meta["descr"], meta["shape"]
    fp.seek(off + hlen)
    n = 1
    for s in shape: n *= s
    item = 4 if descr.strip("<>|") == "f4" else 8
    raw = fp.read(n * item)
    vals = struct.unpack("<%d%s" % (n, "f" if item == 4 else "d"), raw)
    return vals, shape

def load_pae(path):
    z = zipfile.ZipFile(path)
    name = [n for n in z.namelist() if "pae" in n.lower()][0]
    with z.open(name) as fp:
        return load_npy(fp)

def pae_matrix(path, expect):
    vals, shape = load_pae(path)
    n = shape[0]
    if n != expect:
        print(f"  (PAE matrix is {n}x{n}, expected {expect} tokens — check entity order)")
    return lambda i, j: vals[i * n + j]

# ---------- metrics ----------
def pairs_under(atoms, c1, c2, cut):
    sq = cut * cut
    a = [x for x in atoms if x["chain"] == c1 and x["element"] != "H"]
    b = [x for x in atoms if x["chain"] == c2 and x["element"] != "H"]
    return sum(1 for p in a for q in b
               if (p["x"]-q["x"])**2 + (p["y"]-q["y"])**2 + (p["z"]-q["z"])**2 <= sq)

def iface_plddt(atoms, res_set, chain):
    resids = {i for (i, _) in res_set}
    bs = [a["b"] for a in atoms if a["chain"] == chain and a["resid"] in resids and a["atom"] == "CA"]
    if not bs: return float("nan")
    v = sum(bs) / len(bs)
    return v * 100 if v <= 1.5 else v

def token_bases(entities):
    """entities -> {chain_id: (base, length)} in token space (0-based)."""
    base = 0
    out = {}
    for e in entities:
        for cid in e["chain_ids"]:
            out[cid] = (base, len(e["value"]))
        base += len(e["value"])
    return out, base

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Boltz input JSON (entity order)")
    ap.add_argument("--binder", required=True, help="binder chain id, e.g. B")
    ap.add_argument("--cif", required=True)
    ap.add_argument("--pae", default=None)
    args = ap.parse_args()

    d = json.load(open(args.input))
    entities = d["entities"]
    bases, ntok = token_bases(entities)
    binder = args.binder
    partners = [c for c in bases if c != binder]

    atoms = read_cif_atoms(args.cif)
    pae = pae_matrix(args.pae, ntok) if args.pae else None
    bb, blen = bases[binder]

    for partner in partners:
        sb, sp = interface(atoms, partner, binder, 4.5)
        bf = interface(atoms, binder, partner, 4.5)[0]
        ncon = pairs_under(atoms, binder, partner, 4.5)
        print(f"[{binder}-{partner}] 接触原子对 {ncon:4d} | {partner}侧界面残基 {len(sb):2d} 个 "
              f"pLDDT {iface_plddt(atoms, sb, partner):.1f} | {binder}侧 {fmt(bf)} "
              f"pLDDT {iface_plddt(atoms, bf, binder):.1f}")
        if pae:
            pb, plen = bases[partner]
            ib = sorted(i for (i, _) in bf)
            ip = sorted(i for (i, _) in sb)
            if ib and ip:
                m = sum(pae(bb + i, pb + j) for i in ib for j in ip) / (len(ib) * len(ip))
                print(f"   界面局部 PAE 均值: {m:.2f} A")
            f = sum(pae(bb + i, pb + j) for i in range(blen) for j in range(plen)) / (blen * plen)
            print(f"   全链 PAE 均值: {f:.2f} A")

if __name__ == "__main__":
    main()
