# Product Requirements Document
## High-Throughput Molecular Docking & Target Interaction Pipeline

**Version:** 1.0  
**Author:** Developer  
**Status:** Active  
**Last Updated:** 2026-11-07

---

## 1. Problem Statement

Drug discovery requires screening thousands of small molecules against protein targets. Wet-lab binding assays are expensive and slow. Computational docking provides a fast, cheap first filter — but most tools require manual setup per molecule, making high-throughput screening painful for a solo researcher or early-stage lab.

This project builds a fully automated, reproducible pipeline that takes a protein PDB ID and a CSV of SMILES strings and produces ranked docking results with interaction analysis — no manual GUI work required.

---

## 2. User

A first-year B.Tech Biotechnology student building a portfolio piece in computer-aided drug design (CADD). Must be able to explain every part in an interview. The tool must also be usable by other bioinformatics students who understand docking but may not be software engineers.

---

## 3. Goals

| # | Goal | Priority |
|---|---|---|
| G1 | One command from PDB ID + CSV → ranked results | Must-have |
| G2 | Scientific validation: redocking RMSD and enrichment AUC | Must-have |
| G3 | Interaction analysis for top hits (H-bonds, hydrophobic) | Must-have |
| G4 | Resume-on-failure for large screens | Must-have |
| G5 | Clean, documented repo suitable as a portfolio piece | Must-have |
| G6 | Optional web dashboard to explore results | Nice-to-have |
| G7 | KRAS G12C as second validated target | Nice-to-have |

---

## 4. Non-Goals

- No machine-learning scoring models
- No molecular dynamics simulations
- No GUI desktop application
- Do not claim Vina scores prove biological activity
- No wet-lab integration

---

## 5. Primary Target

**EGFR kinase** (Epidermal Growth Factor Receptor) — a well-validated cancer drug target.

| Property | Value |
|---|---|
| PDB ID | 1M17 |
| Chain | A |
| Co-crystal ligand | AQ4 (erlotinib analogue) |
| Key biology | Inhibitors must H-bond with hinge residue Met793 |
| Clinical relevance | Mutations drive lung cancer; erlotinib, gefitinib, osimertinib are approved drugs |

**Optional second target:** KRAS G12C — PDB `6OIM`, ligand `MOV` (sotorasib).

---

## 6. Pipeline Stages

| Stage | Input | Output | Acceptance Criterion |
|---|---|---|---|
| 1. Fetch + clean receptor | PDB ID, chain | `receptor_clean.pdb` | One chain, no water/HETATM |
| 2. Receptor PDBQT | clean PDB | `receptor.pdbqt` | Open Babel exits 0, hydrogens visible |
| 3. Grid box | co-crystal ligand resname | centre (x,y,z), size (x,y,z) | Box encloses ATP pocket |
| 4. Ligand prep | SMILES string | PDBQT string | Invalid SMILES skipped, run continues |
| 5. Docking | receptor + ligands | best pose + affinity | Deterministic with fixed seed |
| 6. Interactions | best pose + receptor | H-bond and hydrophobic lists | Known inhibitors contact Met793 |
| 7. Reporting | all results | `docking_results.csv` sorted by affinity | One row per ligand including failures |

---

## 7. Input / Output Contract

**Input CSV columns:**

```
name, smiles[, is_active]
```

`is_active` is optional; include it when running the enrichment validation.

**Output CSV columns:**

```
name, smiles, affinity_kcal_mol, hbond_residues, hydrophobic_residues, status
```

`status` values: `success`, `prep_failed`, `docking_failed`.

---

## 8. Validation Plan (required, not optional)

1. **Sanity test:** erlotinib and gefitinib score better (more negative) than aspirin, caffeine, ibuprofen.
2. **Redocking:** co-crystal ligand docked back; heavy-atom RMSD target < 2 Å.
3. **Enrichment:** ChEMBL EGFR actives (IC50 < 1 µM) vs DUD-E egfr decoys; report ROC AUC and EF at 1%.
4. **Biology check:** top hits contact hinge residue Met793.

**Rule:** if any validation step fails, stop and diagnose before proceeding to the next phase.

---

## 9. Phased Delivery Plan

| Phase | Deliverable | Done When |
|---|---|---|
| 0 | Environment + repo scaffold | `python docking_pipeline.py --help` works |
| 1 | Sanity run | 5-molecule run completes, scores ordered correctly |
| 2 | Redocking validation | `validate_redock.py` reports RMSD with real number |
| 3 | Library preparation | `prepare_library.py` outputs 500+ row CSV |
| 4 | Full screen | Resume-on-failure + tqdm progress bar |
| 5 | Analysis suite | ROC AUC, EF1%, affinity plot, top-5 table, PLIP/ProLIF upgrade |
| 6 | Visualisation instructions | PyMOL/ChimeraX steps for top 3 poses |
| 7 | Final polish | README, all tests pass, one-command reproduction |

---

## 10. Definition of Done

- [ ] One documented command reproduces all results from scratch
- [ ] Redocking RMSD reported with actual number
- [ ] ROC AUC and EF1% reported with actual numbers
- [ ] Top 5 leads listed with affinity (kcal/mol) and contacted residues
- [ ] README states limitations honestly
- [ ] All `pytest` tests pass
- [ ] No secrets or large binary files in git

---

## 11. Honesty Rule

The targets below are aspirations, not guaranteed outcomes. Report only what the actual run produces.

- Aspiration: screen 500+ compounds against EGFR
- Aspiration: top 5 leads with affinity ≤ -8.5 kcal/mol and mapped H-bond interactions

If real results differ, resume bullets use the real numbers.
