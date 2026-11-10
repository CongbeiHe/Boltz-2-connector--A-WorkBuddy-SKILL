---
name: boltz2-structure-prediction
description: Use this skill when the task involves Boltz / Boltz-2 / boltz-api / Boltz Bio cloud structure prediction: submitting protein-protein / protein-peptide / protein-ligand complex predictions, cyclic peptide (head-to-tail backbone cyclization) modeling, distance/binding constraints, running multi-sample predictions, interpreting metrics (pTM, ipTM, pLDDT, PAE, binding_confidence), downloading and converting predicted mmCIF to PDB, or analyzing constraint satisfaction and interface quality. Trigger words: boltz, boltz-api, structure-and-binding, 复合物结构预测, 结构预测任务, 环肽建模, 界面约束, ipTM 分析.
agent_created: true
---

# Boltz-2 Structure Prediction (boltz-api CLI)

## Overview

This skill drives the Boltz Bio cloud service (`boltz-api` CLI) end-to-end: authentication, input JSON construction (entities / constraints / bonds / binding), cost estimation, multi-sample submission, result download, CIF→PDB conversion, and post-run analysis (constraint satisfaction, interface composition, per-interface quality with PAE). All bundled scripts are pure Python standard library — no numpy/gemmi/BioPython needed.

## Prerequisites

1. Install the CLI: `curl -fsSL https://install.boltz.bio/boltz-api/install.sh | sh`
2. Log in once: `boltz-api auth login` (browser device-code flow).
3. **Always invoke boltz-api with `BOLTZ_API_TEST_DISABLE_KEYRING=1`** — file-backed credentials. The macOS keychain write is blocked in sandboxed environments and the refresh token is then silently lost. See `references/pitfalls.md` for the full failure signature.
4. Sessions expire after ~1h and token refresh can fail; if a command returns `auth_required`, re-login via device code. The pipeline script automates this.

## Workflow

### 1. Build the input JSON

Full schema (field names, enums, hard-won from server-side validation errors — the official docs only document the UI) is in `references/api-schema.md`. Core shape:

```json
{
  "entities": [
    {"type": "protein", "chain_ids": ["A"], "value": "<SEQ1>"},
    {"type": "protein", "chain_ids": ["B"], "value": "<SEQ2>"}
  ],
  "constraints": [
    {"type": "contact",
     "token1": {"type": "polymer_contact", "chain_id": "A", "residue_index": 25},
     "token2": {"type": "polymer_contact", "chain_id": "B", "residue_index": 1},
     "max_distance_angstrom": 4.0, "force": true}
  ],
  "bonds": [
    {"atom1": {"type": "polymer_atom", "chain_id": "B", "residue_index": 6, "atom_name": "C"},
     "atom2": {"type": "polymer_atom", "chain_id": "B", "residue_index": 0, "atom_name": "N"}}
  ],
  "binding": {"type": "protein_protein_binding", "binder_chain_ids": ["B"]},
  "num_samples": 5,
  "model_options": {"seed": 42, "recycling_steps": 3, "sampling_steps": 200}
}
```

Critical details:
- `residue_index` in constraints/bonds is **0-based**; CIF residue numbers are 1-based.
- `force: false` (soft) constraints are almost always ignored — use `force: true` when geometry matters, but note it steers diffusion and biases ipTM (cite as "constrained").
- Cyclic peptide head-to-tail cyclization = a `bonds` entry joining last-residue C to first-residue N.
- MSA is generated server-side automatically; no input needed.
- `--model boltz-2.2` is accepted by `estimate-cost` but rejected by `start`/`run` (`Unsupported parameter: body.model`) — effectively only the default `boltz-2.1` runs.

### 2. Run the pipeline

```bash
python3 scripts/boltz_pipeline.py --input complex.json --name my_complex \
    --root-dir boltz_runs --analyze
```

Automates: auth check + device-code re-login (prints the confirmation URL when manual approval is needed) → `estimate-cost` (aborts on schema errors with the server message) → `run` → mid-run auth-expiry recovery via `download-results --id` → best-sample CIF→PDB → metrics table (per-sample pTM/ipTM/pLDDT, binding_confidence) → constraint/interface report.

Manual equivalent:

```bash
export BOLTZ_API_TEST_DISABLE_KEYRING=1
boltz-api predictions:structure-and-binding estimate-cost --input @json://complex.json
boltz-api predictions:structure-and-binding run --input @json://complex.json \
    --name my_complex --root-dir boltz_runs --progress-format text
# if run dies mid-flight: find id via `predictions:structure-and-binding list`,
# then `download-results --id sab_pred_XXX --name my_complex --root-dir boltz_runs`
```

Results land in `<root-dir>/<name>/outputs/files/prediction/`: best `*_predicted.cif`, `sample_N_predicted_structure.cif`, `sample_N_pae.npz`, `metrics.json`.

### 3. Analyze

```bash
# constraint satisfaction (reads constraints from the input JSON, 0->1-based aware),
# chain composition + pLDDT, pairwise interfaces, optional bond/cyclization check;
# directory mode aggregates median/range over all samples:
python3 scripts/analyze_structure.py --input complex.json \
    --cif boltz_runs/my_complex/outputs/files/prediction --bond "B:7:C,B:1:N"

# compare a binder's interfaces against each partner (contacts, interface pLDDT,
# interface-restricted and whole-chain PAE; token bases derived from entity order):
python3 scripts/compare_ifaces.py --input complex.json --binder B \
    --cif boltz_runs/my_complex/outputs/files/prediction/sample_0_predicted_structure.cif \
    --pae boltz_runs/my_complex/outputs/files/prediction/sample_0_pae.npz

# CIF -> PDB (B-factor = pLDDT*100), opens in PyMOL/ChimeraX:
python3 scripts/cif2pdb.py sample.cif out.pdb
```

## Interpreting metrics

- ipTM (interface) is the headline number for complex plausibility; pTM is dominated by large well-folded chains and overstates ternary confidence.
- `binding_confidence` can contradict ipTM by 2-3 orders of magnitude — treat as unreliable for ranking.
- Single-sample constraint hits can be seed luck: always check the median over ≥5 samples before drawing conclusions.
- Narrow distance distributions across samples that violate constraints indicate systematic geometric incompatibility of the anchor residues, not sampling noise — force cannot fix it.
- A binder's better interface is identified by larger contact counts AND lower whole-chain PAE vs that partner (see compare_ifaces.py).

## Resources

- `scripts/boltz_pipeline.py` — end-to-end submission wrapper (auth, cost, run, recovery, PDB, metrics).
- `scripts/analyze_structure.py` — constraint/interface/pLDDT analysis, multi-sample aggregation.
- `scripts/compare_ifaces.py` — binder-vs-partner interface quality with PAE (stdlib npz/npy parser).
- `scripts/cif2pdb.py` — mmCIF→PDB conversion.
- `references/api-schema.md` — exact input JSON schema, field enums, cost table, CLI subcommands.
- `references/pitfalls.md` — auth/session/sandbox failure modes and fixes.
- `examples/` — smoke test, constrained binary complex, constrained ternary complex with cyclic peptide (all previously executed successfully; ~$0.025-0.04 per sample).
