# PROJECT BRIEF: High-Throughput Molecular Docking & Target Interaction Pipeline

> Read this file fully before writing any code. It is the source of truth for scope, rules and definition of done.

## 1. One-paragraph summary

Build a Python command-line pipeline that takes (a) a protein drug target from the RCSB Protein Data Bank and (b) a CSV of small molecules (SMILES strings), docks every molecule into the protein's binding pocket with AutoDock Vina, ranks them by predicted binding energy, and reports which protein residues each top molecule contacts (hydrogen bonds, hydrophobic contacts). Primary target: **EGFR kinase** (cancer drug target). Optional second target: **KRAS G12C**. The project is a portfolio piece for computer-aided drug design (CADD) and structural bioinformatics.

## 2. About the developer (adjust your explanations)

- First-year B.Tech Biotechnology student, strong on biology, still learning software engineering and docking.
- Explain non-obvious code and every scientific choice in plain language. Prefer readable code over clever code.
- The developer must be able to explain every line in an interview. Never add code you cannot justify in a comment or in your reply.

## 3. Background (so you do not need to guess)

| Term | Meaning here |
|---|---|
| Docking | Computationally placing a small molecule (ligand) into a protein pocket and scoring how well it fits |
| Ligand | The small molecule being tested (a candidate drug) |
| Receptor | The protein target (from a PDB file) |
| SMILES | Text encoding of a molecule, e.g. `CC(=O)Oc1ccccc1C(=O)O` is aspirin |
| PDBQT | File format AutoDock Vina needs: 3D coordinates + partial charges + atom types (+ rotatable bonds for ligands) |
| Grid box | 3D box around the binding pocket; Vina only searches inside it |
| Binding affinity | Vina's score in kcal/mol. More negative = predicted tighter binding. It is an estimate, not a measurement |
| Co-crystal ligand | The inhibitor already bound in the PDB structure; used to locate the pocket and to validate docking |
| Redocking | Docking the co-crystal ligand back and measuring RMSD to the real pose. RMSD < 2 Å means the setup is trustworthy |
| Decoys | Molecules that look like drugs but are not true binders; used to test whether docking can tell actives from inactives |
| Hinge (EGFR Met793) | Key kinase residue that known EGFR inhibitors form a hydrogen bond with |

## 4. Goals and non-goals

**Goals**
1. Fully automated, reproducible pipeline: one command from PDB ID + CSV to ranked results.
2. Scientific validation: redocking RMSD and actives-vs-decoys enrichment (ROC AUC).
3. Interaction analysis for top hits (H-bonds, hydrophobic contacts), compared against known EGFR biology.
4. Clean repo: README, tests, figures, requirements.

**Non-goals (do not build unless asked)**
- No web app, GUI or database.
- No machine-learning scoring models.
- No molecular dynamics (may be mentioned as future work).
- Do not claim Vina scores prove a molecule is a drug.

## 5. Tech stack

- Python 3.11, conda environment
- `rdkit` and `openbabel` (install via conda-forge), `vina` (AutoDock Vina Python bindings), `meeko`, `biopython`, `pandas`, `numpy`, `matplotlib`, `scikit-learn` (ROC AUC), `pytest`
- Visualisation (manual, outside the pipeline): PyMOL or ChimeraX
- Later upgrade: `prolif` or `plip` for proper interaction fingerprints

Install:
```bash
conda create -n docking python=3.11 -y && conda activate docking
conda install -c conda-forge rdkit openbabel -y
pip install vina meeko biopython pandas numpy matplotlib scikit-learn pytest
```

## 6. Existing starter code

`docking_pipeline.py` already implements the core flow (fetch PDB, clean receptor, box from co-crystal ligand, SMILES to PDBQT, parallel Vina docking, distance-based interaction extraction, ranked CSV). It has **not been tested against live data**, so version issues with Meeko, Vina or Open Babel are expected. Fix them minimally; do not rewrite the file from scratch.

## 7. Target repository layout

```
docking-pipeline/
  PROJECT_BRIEF.md
  README.md
  requirements / environment.yml
  docking_pipeline.py          # main CLI
  validate_redock.py           # redocking RMSD
  prepare_library.py           # clean ChEMBL/DUD-E data -> library.csv
  analyze_results.py           # ROC/enrichment, plots, top-hit table
  data/                        # library CSVs (small files only)
  results_<target>/            # outputs (git-ignored except summaries)
  tests/
```

## 8. Pipeline stages and contracts

| # | Stage | Input | Output | Acceptance check |
|---|---|---|---|---|
| 1 | Fetch + clean receptor | PDB ID, chain | `receptor_clean.pdb` (one chain, no water/hetero) | Opens in PyMOL, only protein remains |
| 2 | Receptor PDBQT | clean PDB | `receptor.pdbqt` | Open Babel exits 0, hydrogens present |
| 3 | Grid box | co-crystal ligand resname | centre (x,y,z), size (x,y,z) | Box visibly encloses the ATP pocket |
| 4 | Ligand prep | SMILES | ligand PDBQT string | Invalid SMILES skipped and logged, not crashing the run |
| 5 | Docking | receptor + ligands | best pose + affinity per ligand | Seeded, deterministic enough to re-run |
| 6 | Interactions | best pose + receptor | H-bond and hydrophobic residue lists | Known inhibitors show hinge contact |
| 7 | Reporting | all results | `docking_results.csv` sorted by affinity | One row per ligand incl. failures with status |

Input CSV format: `name,smiles` (plus optional `is_active` column for validation).

## 9. Targets (verify on rcsb.org before use)

| Target | PDB | Chain | Co-crystal ligand |
|---|---|---|---|
| EGFR wild-type with erlotinib | 1M17 | A | AQ4 |
| KRAS G12C with sotorasib | 6OIM | A | MOV |

## 10. Validation plan (required, not optional)

1. **Sanity test:** erlotinib and gefitinib must score better than aspirin, caffeine, ibuprofen.
2. **Redocking:** dock the co-crystal ligand back; report heavy-atom RMSD; target < 2 Å.
3. **Enrichment:** ChEMBL EGFR actives (CHEMBL203, IC50 < 1 µM) vs DUD-E `egfr` decoys; report ROC AUC and enrichment factor at 1%.
4. **Biology check:** do top hits contact EGFR hinge residue Met793?

If a validation step fails, stop and diagnose before building further.

## 11. Rules for the AI agent

1. Work **one phase at a time** (see Section 12). Stop after each phase and summarise what changed and how to verify it.
2. Never change docking parameters, scoring function, box size or random seeds silently. State every change and why.
3. Run the sanity test after any change to docking, ligand prep or receptor prep.
4. Never fabricate results, scores, RMSD values or PDB/ligand IDs. If unsure, say so and tell the developer to verify.
5. Keep dependencies to the list in Section 5. Ask before adding more.
6. Handle failures per ligand (log, mark status, continue). One bad molecule must never kill a 500-compound run.
7. Add type hints and short docstrings; add a pytest for each new function that has testable logic.
8. Do not delete or overwrite result files without asking.

## 12. Phased plan (agent works through these in order)

| Phase | Deliverable |
|---|---|
| 0 | Environment, repo, `--help` works |
| 1 | Sanity run on 5 molecules passes |
| 2 | `validate_redock.py` with RMSD < 2 Å |
| 3 | `prepare_library.py` produces `library.csv` (500+ rows, actives + decoys) |
| 4 | Full screen completes with resume-on-failure and progress bar |
| 5 | `analyze_results.py`: ROC AUC, EF1%, affinity plot, top-5 table, hinge-contact check; upgrade interactions to PLIP/ProLIF |
| 6 | PyMOL/ChimeraX images of top 3 poses (instructions only; developer renders) |
| 7 | README with results and limitations, tests, one-command reproduction |

## 13. Definition of done

- One documented command reproduces the full results from scratch.
- Redocking RMSD, ROC AUC and EF1% are reported with real numbers.
- Top 5 leads are listed with affinity and contacted residues.
- README states limitations honestly (Vina is a rough estimate; distance-only H-bond detection if not upgraded).
- All tests pass.

## 14. Honesty rule for resume claims

The targets below are **aspirations, not facts**. Report only what the actual run produces.

- Aspiration: screen 500+ compounds against EGFR/KRAS.
- Aspiration: top 5 leads with affinity <= -8.5 kcal/mol and mapped H-bond interactions.

If real results differ (fewer compounds, weaker scores), the resume bullets use the real numbers.

## 15. Known pitfalls

- Meeko and Vina Python APIs change between versions; check installed versions before assuming function signatures.
- Open Babel receptor conversion can mis-assign atom types; inspect the receptor in PyMOL if scores look odd.
- Vina scores correlate poorly with measured affinity; enrichment and redocking are the honest metrics.
- Mutant EGFR PDB IDs and ligand codes must be checked on rcsb.org; do not guess.
