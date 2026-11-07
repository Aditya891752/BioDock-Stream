# Task Tracker
## High-Throughput Molecular Docking Pipeline — GSD Style

**Version:** 1.0  
**Last Updated:** 2026-11-07  
**Current Phase:** 2 — Redocking Validation

> Format: `[ ]` = not started · `[~]` = in progress · `[x]` = done · `[!]` = blocked

---

## Phase 0 — Environment & Repo Scaffold

**Goal:** `python docking_pipeline.py --help` works without error.

**Verification:** Run `python docking_pipeline.py --help` and confirm output shows all CLI flags.

| # | Task | Status | Notes |
|---|---|---|---|
| 0.1 | Create conda env `docking` with Python 3.11 | [ ] | `conda create -n docking python=3.11 -y` |
| 0.2 | Install rdkit + openbabel via conda-forge | [ ] | `conda install -c conda-forge rdkit openbabel -y` |
| 0.3 | Install vina + meeko + other pip deps | [ ] | `pip install vina meeko biopython pandas numpy matplotlib scikit-learn pytest` |
| 0.4 | Confirm installed versions (`pip show vina meeko`) | [ ] | Record in `memory.md` → Environment State |
| 0.5 | Initialize git repo + first commit | [ ] | Message: "Phase 0: initialize repo" |
| 0.6 | Create `.gitignore` (env, results dirs, pdbqt outputs, .env) | [ ] | |
| 0.7 | Create `docs/` directory and commit all planning docs | [ ] | This directory |
| 0.8 | Create `data/` directory + placeholder `molecules_sanity.csv` | [ ] | 5 molecules: erlotinib, gefitinib, aspirin, caffeine, ibuprofen |
| 0.9 | Create `tests/` directory + empty `conftest.py` | [ ] | |
| 0.10 | Create `docking_pipeline.py` scaffold with `--help` and all CLI args | [ ] | Use `argparse`; no docking logic yet |
| 0.11 | Write `environment.yml` and `requirements.txt` | [ ] | After install is working |
| 0.12 | Write first pytest: `tests/test_docking_pipeline.py::test_cli_help` | [ ] | Calls `--help` and checks exit code 0 |

**Phase 0 done when:** `python docking_pipeline.py --help` exits 0, `pytest tests/` passes, repo is committed.

---

## Phase 1 — Sanity Run (5 Molecules)

**Goal:** Pipeline runs end-to-end on 5 known molecules. Erlotinib and gefitinib score better than aspirin/caffeine/ibuprofen.

**Verification:** Inspect `results_egfr/docking_results.csv` — affinities ordered correctly.

| # | Task | Status | Notes |
|---|---|---|---|
| 1.1 | Implement `fetch_receptor(pdb_id)` — download 1M17 from RCSB | [ ] | HTTP GET to `https://files.rcsb.org/download/1M17.pdb` |
| 1.2 | Implement `clean_receptor(pdb_path, chain)` — keep chain A, strip HETATM/ANISOU | [ ] | Biopython `PDBParser` + `PDBIO` with `Select` |
| 1.3 | Implement `prepare_receptor(clean_pdb)` — Open Babel → PDBQT | [ ] | `subprocess` call to `obabel` |
| 1.4 | Implement `define_grid_box(pdb_path, ligand_resname)` — centroid of AQ4 atoms | [ ] | Biopython; add 10 Å padding each axis |
| 1.5 | Implement `smiles_to_pdbqt(name, smiles)` — RDKit 3D + Meeko → PDBQT string | [ ] | Catch invalid SMILES; return `None` with log |
| 1.6 | Implement `dock_ligand(receptor_pdbqt, ligand_pdbqt, box)` — Vina | [ ] | Fixed seed=42; exhaustiveness=8 |
| 1.7 | Implement `extract_interactions_distance(pose_pdbqt, receptor_pdb)` — distance-based | [ ] | H-bond: N/O within 3.5 Å; hydrophobic: C within 4.5 Å |
| 1.8 | Implement `write_results(results_list, output_path)` — pandas CSV | [ ] | Sorted by affinity ascending |
| 1.9 | Wire all stages into `main()` with `argparse` | [ ] | `--pdb-id`, `--library`, `--chain`, `--ligand-resname`, `--output-dir` |
| 1.10 | Run sanity test on `molecules_sanity.csv` | [ ] | Record real scores in `memory.md` |
| 1.11 | Write pytests for each function (1.1–1.8) | [ ] | Mock RCSB HTTP call; use test fixtures |

**Phase 1 done when:** sanity run CSV shows erlotinib and gefitinib scoring more negative than aspirin/caffeine/ibuprofen.

---

## Phase 2 — Redocking Validation

**Goal:** `validate_redock.py` docks co-crystal ligand AQ4 back and reports heavy-atom RMSD.

**Verification:** RMSD value is a real number. If > 2 Å, diagnose before continuing.

| # | Task | Status | Notes |
|---|---|---|---|
| 2.1 | Obtain AQ4 SMILES (verify on ChEMBL or PubChem) | [ ] | Do NOT guess; look it up |
| 2.2 | Extract crystallographic AQ4 coordinates from 1M17.pdb | [ ] | Biopython; save as `aq4_crystal.pdb` |
| 2.3 | Implement `validate_redock.py` — dock AQ4, compare to crystal pose | [ ] | Use same box and params as pipeline |
| 2.4 | Implement `compute_rmsd(pose_coords, crystal_coords)` | [ ] | Heavy atoms only (no H); Biopython or numpy |
| 2.5 | Report: RMSD value, pass/fail vs 2 Å, top Vina score | [ ] | Print and write to `results_egfr/redock_report.txt` |
| 2.6 | Write pytests for `compute_rmsd` | [ ] | Known coordinate pairs with known RMSD |

**Phase 2 done when:** real RMSD number is reported. If > 2 Å, diagnose (box wrong? ligand prep wrong?) before Phase 3.

---

## Phase 3 — Library Preparation (500+ compounds)

**Goal:** `prepare_library.py` outputs `data/library.csv` with 500+ rows, actives + decoys.

**Verification:** `wc -l data/library.csv` ≥ 501 (header + 500 rows). Actives/decoys ratio confirmed.

| # | Task | Status | Notes |
|---|---|---|---|
| 3.1 | Query ChEMBL REST API for EGFR actives (CHEMBL203, IC50 < 1 µM) | [ ] | `requests` + ChEMBL REST API v3 |
| 3.2 | Download DUD-E egfr decoys | [ ] | `http://dude.docking.org/targets/egfr` |
| 3.3 | Parse and clean SMILES via RDKit (remove invalid, salt-strip) | [ ] | Log rejected SMILES |
| 3.4 | Merge actives + decoys into `library.csv` with `is_active` column | [ ] | Shuffle rows |
| 3.5 | Verify at least 500 rows total, at least 50 actives | [ ] | Print counts |
| 3.6 | Write pytests for SMILES cleaning functions | [ ] | |

**Phase 3 done when:** `library.csv` has ≥ 500 rows with valid SMILES and `is_active` column.

---

## Phase 4 — Full Screen with Resume and Progress Bar

**Goal:** Screen full library with resume-on-failure and `tqdm` progress bar.

**Verification:** Kill the run halfway through; restart; it resumes from where it left off.

| # | Task | Status | Notes |
|---|---|---|---|
| 4.1 | Add checkpoint file (`results_egfr/progress.json`) — record completed ligands | [ ] | Write after each successful dock |
| 4.2 | On startup, skip already-docked ligands | [ ] | Resume logic in `main()` |
| 4.3 | Add `tqdm` progress bar to ligand docking loop | [ ] | Show: completed/total, current ligand name, ETA |
| 4.4 | Add `--workers N` flag for parallel docking | [ ] | `ProcessPoolExecutor(max_workers=N)` |
| 4.5 | Run full screen on `library.csv` | [ ] | Record real runtime in `memory.md` |
| 4.6 | Verify all failures are logged, run completes | [ ] | Check `run.log` for any unexpected crashes |

**Phase 4 done when:** full library screens to completion with progress bar; interrupted run resumes correctly.

---

## Phase 5 — Analysis Suite

**Goal:** `analyze_results.py` produces ROC AUC, EF1%, affinity plot, top-5 table, hinge-contact check.

**Verification:** ROC AUC and EF1% are real numbers. Top-5 table shows Met793 contact column.

| # | Task | Status | Notes |
|---|---|---|---|
| 5.1 | Install ProLIF: `pip install prolif` | [ ] | Record version in `memory.md` |
| 5.2 | Upgrade `extract_interactions` to use ProLIF fingerprints | [ ] | Keep old distance-based as fallback |
| 5.3 | Implement ROC AUC calculation (scikit-learn) | [ ] | Use `is_active` column as ground truth |
| 5.4 | Implement EF at 1% | [ ] | Fraction of actives in top 1% of ranked list |
| 5.5 | Generate ROC curve plot → `results_egfr/roc_curve.png` | [ ] | matplotlib; single colour, no gradient |
| 5.6 | Generate affinity distribution histogram → `results_egfr/affinity_hist.png` | [ ] | Separate bars for actives vs decoys |
| 5.7 | Generate top-5 hits table → `results_egfr/top5_hits.csv` | [ ] | name, affinity, hbond_residues, met793_contact |
| 5.8 | Hinge contact check: flag rows where Met793 is in `hbond_residues` | [ ] | Boolean column |
| 5.9 | Write pytests for ROC and EF calculations | [ ] | Known scores/labels with known expected output |

**Phase 5 done when:** ROC AUC and EF1% reported, top-5 CSV exists, Met793 contact checked.

---

## Phase 6 — Visualisation Instructions

**Goal:** Developer can render top 3 poses in PyMOL/ChimeraX using the pipeline outputs.

**Verification:** Developer confirms images look correct (colour by interaction, pocket visible).

| # | Task | Status | Notes |
|---|---|---|---|
| 6.1 | Write PyMOL instruction set for loading receptor + top 3 poses | [ ] | In `README.md` or separate `VISUALISATION.md` |
| 6.2 | Write ChimeraX alternative instructions | [ ] | |
| 6.3 | Document which residues to highlight (Met793, known contacts) | [ ] | |

---

## Phase 7 — Final Polish

**Goal:** One-command reproduction; all tests pass; README complete.

**Verification:** Clone repo on a fresh machine, follow README, pipeline runs to completion.

| # | Task | Status | Notes |
|---|---|---|---|
| 7.1 | Write `README.md` with one-command reproduction | [ ] | Include environment setup + pipeline command |
| 7.2 | Add results summary table to README (real numbers only) | [ ] | RMSD, ROC AUC, EF1%, top-5 affinities |
| 7.3 | State limitations honestly in README | [ ] | Vina accuracy, distance-based H-bond caveats |
| 7.4 | Run full `pytest tests/` — confirm all pass | [ ] | |
| 7.5 | Check `.gitignore` excludes all large files | [ ] | No PDBQT, no conda env binaries |
| 7.6 | Pin `environment.yml` versions | [ ] | `conda env export > environment.yml` |
| 7.7 | Final commit: "Phase 7: complete" | [ ] | |

---

## Backlog (not yet scheduled)

| Item | Phase Candidate |
|---|---|
| KRAS G12C target (6OIM) | After Phase 7 |
| FastAPI backend for web dashboard | Phase 8 |
| Next.js frontend on Vercel | Phase 9 |
| AlphaFold structure for mutant EGFR | Future |
| Gemini API query interface | Future |

---

## Progress Summary

| Phase | Status | Real Deliverable |
|---|---|---|
| 0 | Not started | — |
| 1 | Not started | — |
| 2 | Not started | — |
| 3 | Not started | — |
| 4 | Not started | — |
| 5 | Not started | — |
| 6 | Not started | — |
| 7 | Not started | — |
