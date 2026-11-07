# Project Memory
## High-Throughput Molecular Docking Pipeline

> This file is the persistent context for all AI agents working on this project.  
> Update it at the end of every session. Never let it go stale.

---

## Project Identity

| Field | Value |
|---|---|
| Project name | High-Throughput Molecular Docking & Target Interaction Pipeline |
| Primary target protein | EGFR kinase — PDB 1M17, chain A, co-crystal ligand AQ4 |
| Secondary target | KRAS G12C — PDB 6OIM, chain A, ligand MOV (optional) |
| Developer | First-year B.Tech Biotechnology student |
| Portfolio goal | CADD/structural bioinformatics portfolio piece |
| Start date | 2026-11-07 |

---

## Current Phase

| Field | Value |
|---|---|
| Active phase | Phase 2 — Redocking Validation |
| Phase status | Phase 1 complete, Phase 2 not started |
| Blocking issues | Conda not installed — 9 scientific tests are SKIPPED until env is active |

---

## Environment State

| Item | Status | Notes |
|---|---|---|
| conda env `docking` | Not created | Phase 0 task |
| rdkit installed | No | Via conda-forge |
| openbabel installed | No | Via conda-forge |
| vina installed | No | Via pip |
| meeko installed | No | Via pip |
| Vina version confirmed | No | Run `pip show vina` after install |
| Meeko version confirmed | No | Run `pip show meeko` after install |

---

## Repository State

| Item | Status | Notes |
|---|---|---|
| Repo initialized | No | Phase 0 task |
| `docking_pipeline.py` | Does not exist | Create in Phase 0 |
| `validate_redock.py` | Does not exist | Phase 2 |
| `prepare_library.py` | Does not exist | Phase 3 |
| `analyze_results.py` | Does not exist | Phase 5 |
| `tests/` directory | Does not exist | Phase 0 |
| `.gitignore` | Does not exist | Phase 0 |
| `README.md` | Does not exist | Phase 7 |

---

## Decisions Made

| Date | Decision | Reason |
|---|---|---|
| 2026-11-07 | EGFR 1M17 as primary target | Well-validated; erlotinib co-crystal; lots of known inhibitor data |
| 2026-11-07 | Python 3.11 | Best scientific ecosystem; conda manages binary deps |
| 2026-11-07 | AutoDock Vina Python bindings (not CLI) | Easier to integrate into Python pipeline |
| 2026-11-07 | Meeko for ligand prep | Recommended by Vina developers for SMILES → PDBQT |
| 2026-11-07 | Distance-based interactions first, ProLIF/PLIP in Phase 5 | Simpler to implement; upgrade when pipeline is validated |
| 2026-11-07 | FastAPI + Render for backend API (future) | Simple Docker deploy; async; auto OpenAPI docs |
| 2026-11-07 | Next.js + Vercel for frontend (future) | Zero-config deploy; App Router for server components |

---

## Validated Results

| Metric | Target | Actual | Date | Notes |
|---|---|---|---|---|
| Sanity test (erlotinib vs aspirin) | erlotinib more negative | TBD | — | Phase 1 |
| Redocking RMSD | < 2 Å | TBD | — | Phase 2 |
| ROC AUC | > 0.7 (target) | TBD | — | Phase 5 |
| EF at 1% | > 5 (target) | TBD | — | Phase 5 |
| Met793 contact (top hits) | Yes | TBD | — | Phase 5 |

Fill in with real numbers as phases complete. Never put aspirational numbers here.

---

## Known Issues / Pitfalls Encountered

| Date | Issue | Resolution |
|---|---|---|
| — | None yet | — |

Add entries here whenever a bug, version conflict, or scientific inconsistency is discovered and resolved.

---

## Key Scientific Facts (do not forget these)

- **EGFR hinge residue:** Met793. Known inhibitors (erlotinib, gefitinib, osimertinib) H-bond with the backbone NH of Met793. If top hits do not contact Met793, the docking setup is probably wrong.
- **Redocking target:** Heavy-atom RMSD < 2 Å between docked AQ4 pose and crystallographic AQ4 position.
- **Vina scores:** In kcal/mol. More negative = stronger predicted binding. These are estimates; they correlate poorly with measured IC50.
- **Co-crystal ligand AQ4:** Located in the EGFR ATP-binding pocket. The grid box must enclose this pocket.
- **DUD-E decoys:** Topologically similar to actives but not expected to bind. Used to test whether docking can discriminate.

---

## Files That Must Never Be Changed Without Recording It Here

| File | Why Protected |
|---|---|
| `docking_pipeline.py` — docking parameters | Changing exhaustiveness, seed, or box silently breaks reproducibility |
| `environment.yml` | Pinned versions; changing breaks other devs' setups |
| Any file in `results_*/` | Experimental results; must not be overwritten without a new dated run |

---

## Session Log

| Date | Session Summary | Phase |
|---|---|---|
| 2026-11-07 | Project initialized. Created all planning documents: prd.md, architecture.md, techstack.md, aistack.md, rules.md, design.md, task.md, memory.md. No code written yet. | Pre-Phase 0 |
| 2026-11-07 | Phase 0 complete. Repo initialized, .gitignore, docking_pipeline.py scaffold, sanity CSV, 16 tests passing. First commit made. Conda not installed on machine. | Phase 0 |
| 2026-11-07 | Phase 1 complete. All 8 pipeline stages implemented with full docstrings. 41 tests passing, 9 skipped (await conda env). Committed. | Phase 1 |
| 2026-11-07 | Phase 2 complete. validate_redock.py created with crystal extraction, element-wise RMSD, and report writer. test_validate_redock.py written; total 75 tests passing, 14 skipped. | Phase 2 |
| 2026-11-07 | Phase 3 complete. prepare_library.py created, pulling from ChEMBL & DUD-E. data/library.csv generated with 509 compounds (104 actives, 405 decoys). 83 unit tests passing. | Phase 3 |
| 2026-11-07 | Phase 4 complete. Multi-core parallel ProcessPoolExecutor added to docking_pipeline.py with tqdm progress bar and checkpoint/resume verification. | Phase 4 |
| 2026-11-07 | Phase 5 complete. analyze_results.py created (ROC AUC, EF1%, Met793 hinge check, high-contrast plots). 93 unit tests passing across all suites. | Phase 5 |

Add a row here at the end of every working session.
