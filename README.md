# boltz2-structure-prediction

End-to-end workflow for [Boltz-2](https://boltz.bio) cloud structure prediction: authentication → input construction (entities / interface constraints / cyclic-peptide covalent cyclization) → cost estimation → multi-sample submission → result download → CIF→PDB conversion → constraint-satisfaction and interface-quality analysis.

This repository doubles as a **WorkBuddy (CodeBuddy-family) Skill**: drop it into `~/.workbuddy/skills/` and simply say "run a Boltz structure prediction" in chat — the full workflow triggers automatically.

**Zero third-party dependencies**: every script is pure Python standard library — mmCIF and PAE (npz/npy) are parsed by hand with `zipfile` + `struct`; no numpy, gemmi, or BioPython required.

## What's inside

| Component | Purpose |
|---|---|
| `scripts/boltz_pipeline.py` | One command, full pipeline: auth check + automatic device-code re-login on session expiry (prints the confirmation URL) → `estimate-cost` (schema errors surfaced verbatim from the server) → `run` → automatic mid-run recovery via `download-results --id` if the CLI dies → best-sample CIF→PDB → metrics summary |
| `scripts/analyze_structure.py` | Constraint satisfaction (constraints read automatically from the input JSON, 0→1-based index conversion), per-chain pLDDT, pairwise interface residues, optional cyclization/covalent-bond check; directory mode aggregates median/range over all samples |
| `scripts/compare_ifaces.py` | Compares a binder's (e.g. cyclic peptide) interfaces against each partner: contact atom pairs, interface-residue pLDDT, interface-restricted PAE, whole-chain PAE (token bases derived automatically from entity order) |
| `scripts/cif2pdb.py` | mmCIF→PDB with B-factor = pLDDT×100; opens directly in PyMOL/ChimeraX |
| `references/api-schema.md` | The input JSON schema as measured by experiment — the official docs only cover the UI; every field name here was reverse-engineered from server-side validation errors |
| `references/pitfalls.md` | Every auth / sandbox / session / constraint pitfall encountered, with fixes |
| `examples/` | Three inputs that have all been executed successfully: a smoke test, a binary complex with 4 Å hard constraints, and a ternary complex (cyclic peptide cyclization + hard constraints + 5 samples) |

## Installation

### 1. Install the Boltz CLI and log in

```bash
curl -fsSL https://install.boltz.bio/boltz-api/install.sh | sh
boltz-api auth login     # browser device-code flow
```

### 2. Install as a WorkBuddy skill

```bash
git clone https://github.com/<you>/boltz2-structure-prediction.git \
    ~/.workbuddy/skills/boltz2-structure-prediction
```

Or use the bundled script: `./install.sh --skill`.

### 3. Works without WorkBuddy too

This is just a set of Python scripts + reference docs. Run `python3 scripts/boltz_pipeline.py ...` directly.

## Quick start

```bash
# Smoke test (~$0.025)
python3 scripts/boltz_pipeline.py --input examples/smoke_test.json \
    --name smoke --root-dir boltz_runs --analyze

# Real complex: constraints + 5 samples + cyclic peptide
python3 scripts/boltz_pipeline.py --input examples/ternary_cyclic_peptide_constrained.json \
    --name my_ternary --root-dir boltz_runs --analyze
```

Analyze existing results (no new submission):

```bash
python3 scripts/analyze_structure.py \
    --input my_input.json --cif boltz_runs/my_ternary/outputs/files/prediction \
    --bond "B:7:C,B:1:N"      # cyclization check: C7–N1 of peptide chain B

python3 scripts/compare_ifaces.py --input my_input.json --binder B \
    --cif .../sample_0_predicted_structure.cif --pae .../sample_0_pae.npz
```

## Input JSON essentials (full details in references/api-schema.md)

```jsonc
{
  "entities": [                                    // declaration order = token order
    {"type": "protein", "chain_ids": ["A"], "value": "..."},
    {"type": "protein", "chain_ids": ["B"], "value": "..."}   // cyclic peptide is also "protein"
  ],
  "constraints": [{                                // distance constraints
    "type": "contact",
    "token1": {"type": "polymer_contact", "chain_id": "B", "residue_index": 1},  // 0-based!
    "token2": {"type": "polymer_contact", "chain_id": "A", "residue_index": 25},
    "max_distance_angstrom": 4.0,
    "force": true                                  // soft constraints measured as ineffective
  }],
  "bonds": [{                                      // head-to-tail peptide cyclization
    "atom1": {"type": "polymer_atom", "chain_id": "B", "residue_index": 6, "atom_name": "C"},
    "atom2": {"type": "polymer_atom", "chain_id": "B", "residue_index": 0, "atom_name": "N"}
  }],
  "binding": {"type": "protein_protein_binding", "binder_chain_ids": ["B"]},
  "num_samples": 5,
  "model_options": {"seed": 42, "recycling_steps": 3, "sampling_steps": 200}
}
```

**Three things to remember**: `residue_index` is 0-based (CIF residue numbers are 1-based); always invoke the CLI with `BOLTZ_API_TEST_DISABLE_KEYRING=1` (built into the pipeline — the macOS keychain write is blocked in sandboxed environments and silently loses the refresh token otherwise); `force:true` steers diffusion and biases ipTM — cite results as "constrained".

## Metrics cheat sheet

- **ipTM** is the headline number for interface plausibility; pTM is inflated by large well-folded chains.
- `binding_confidence` can contradict ipTM by 2–3 orders of magnitude in the current version — do not use it for ranking.
- A single-sample constraint hit can be seed luck — **check the median over ≥5 samples before drawing conclusions**.
- Narrow distance distributions that still violate constraints = the anchor geometry is systematically incompatible with the model's preferred conformation; `force` cannot fix it.

## License

MIT
