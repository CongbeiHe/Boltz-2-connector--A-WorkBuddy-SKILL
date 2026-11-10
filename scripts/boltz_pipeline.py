#!/usr/bin/env python3
"""End-to-end Boltz-2 structure prediction pipeline (wraps the boltz-api CLI).

Usage:
  python3 boltz_pipeline.py --input complex.json --name my_complex \
      [--root-dir boltz_runs] [--analyze]

Steps:
 1. auth check (auto device-code re-login if the session expired)
 2. estimate-cost (input validation + price, aborts on schema errors)
 3. predictions:structure-and-binding run  (submit -> wait -> download)
    - if `run` dies mid-flight (auth expiry etc.), recovers the prediction id
      from `list` and re-downloads with `download-results`
 4. best sample CIF -> PDB (B-factor = pLDDT*100)
 5. metrics summary (pTM / ipTM / pLDDT per sample, binding confidence)
 6. optional --analyze: constraint satisfaction + interface report

NOTE: always runs boltz-api with BOLTZ_API_TEST_DISABLE_KEYRING=1 (file-backed
credentials; the macOS keychain write is blocked in sandboxed environments and
silently loses the refresh token otherwise).
"""
import sys, os, json, glob, argparse, subprocess, time, re

ENV = dict(os.environ, BOLTZ_API_TEST_DISABLE_KEYRING="1")

def run_cli(cmd, timeout=None, print_output=True):
    r = subprocess.run(cmd, capture_output=True, text=True, env=ENV, timeout=timeout)
    if print_output and r.stdout.strip():
        print(r.stdout.strip()[:3000])
    return r

def auth_ok():
    r = run_cli(["boltz-api", "auth", "validate"], print_output=False)
    try:
        return json.loads(r.stdout).get("valid", False)
    except Exception:
        return False

def relogin(wait=180):
    print("[auth] session expired — starting device-code login ...")
    p = subprocess.Popen(["boltz-api", "auth", "login", "--device-code", "--json-events"],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=ENV)
    t0 = time.time()
    url = code = None
    while time.time() - t0 < wait:
        line = p.stdout.readline()
        if not line:
            if p.poll() is not None: break
            time.sleep(0.5); continue
        try: ev = json.loads(line)
        except Exception: continue
        if ev.get("event") == "auth_url":
            url, code = ev.get("verification_uri_complete"), ev.get("user_code")
            print(f"[auth] open this URL and confirm (code {code}):\n        {url}")
        elif ev.get("event") == "success":
            print("[auth] success")
            return True
    print("[auth] not confirmed in time — rerun `boltz-api auth login` manually")
    return False

def find_pred_id(root, name):
    # from the run dir's .boltz-run.json / run.json, else from the remote list
    for f in glob.glob(os.path.join(root, name, ".boltz-run.json")) + \
             glob.glob(os.path.join(root, name, "run.json")):
        try:
            d = json.load(open(f))
            for k in ("prediction_id", "id", "pred_id"):
                if d.get(k): return d[k]
            s = json.dumps(d)
            m = re.search(r"sab_pred_\w+", s)
            if m: return m.group(0)
        except Exception: pass
    r = run_cli(["boltz-api", "predictions:structure-and-binding", "list"], print_output=False)
    m = re.findall(r"sab_pred_\w+", r.stdout)
    return m[0] if m else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--root-dir", default="boltz_runs")
    ap.add_argument("--analyze", action="store_true",
                    help="run analyze_structure.py on the results")
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))

    if not auth_ok():
        if not relogin():
            sys.exit("ABORT: no valid Boltz session")

    print("\n[1/4] estimate-cost")
    r = run_cli(["boltz-api", "predictions:structure-and-binding", "estimate-cost",
                 "--input", f"@json://{args.input}"])
    if r.returncode != 0 or "estimated_cost_usd" not in r.stdout:
        sys.exit(f"ABORT: input rejected:\n{r.stdout}{r.stderr}")
    cost = re.search(r'"estimated_cost_usd":\s*"([^"]+)"', r.stdout)
    print(f"      cost: ${cost.group(1) if cost else '?'}")

    print("\n[2/4] submit + wait + download")
    r = run_cli(["boltz-api", "predictions:structure-and-binding", "run",
                 "--input", f"@json://{args.input}", "--name", args.name,
                 "--root-dir", args.root_dir, "--progress-format", "text"], timeout=3600)
    pred_dir = os.path.join(args.root_dir, args.name, "outputs", "files", "prediction")
    if not os.path.isdir(pred_dir):
        # auth probably died mid-run — recover from the remote
        if not auth_ok() and not relogin():
            sys.exit("ABORT: session lost mid-run")
        pid = find_pred_id(args.root_dir, args.name)
        if not pid:
            sys.exit(f"ABORT: run failed and no prediction id found\n{r.stdout}{r.stderr}")
        print(f"      recovering prediction {pid} via download-results")
        run_cli(["boltz-api", "download-results", "--id", pid, "--name", args.name,
                 "--root-dir", args.root_dir, "--progress-format", "text"], timeout=1800)
    if not os.path.isdir(pred_dir):
        sys.exit("ABORT: results not downloaded")

    print("\n[3/4] convert best sample to PDB")
    best = [f for f in glob.glob(os.path.join(pred_dir, "*_predicted.cif"))]
    if best:
        pdb = os.path.join(pred_dir, f"{args.name}_best.pdb")
        run_cli([sys.executable, os.path.join(here, "cif2pdb.py"), best[0], pdb])
    else:
        print("      (no best-sample cif found)")

    print("\n[4/4] metrics")
    mj = os.path.join(pred_dir, "metrics.json")
    if os.path.exists(mj):
        d = json.load(open(mj))
        for i, s in enumerate(d.get("all_sample_results", [])):
            m = s["metrics"]
            print(f"  sample {i}: ptm {m['ptm']:.3f}  iptm {m['iptm']:.3f}  "
                  f"plddt {m['complex_plddt']:.3f}")
        b = d.get("best_sample", {}).get("metrics", {})
        if b:
            print(f"  BEST: ptm {b['ptm']:.3f}  iptm {b['iptm']:.3f}  "
                  f"plddt {b['complex_plddt']:.3f}")
        print(f"  binding_confidence: {d['binding_metrics']['binding_confidence']:.4f}")
    if args.analyze:
        print("\n[extra] constraint / interface analysis")
        subprocess.run([sys.executable, os.path.join(here, "analyze_structure.py"),
                        "--input", args.input, "--cif", pred_dir])

if __name__ == "__main__":
    main()
