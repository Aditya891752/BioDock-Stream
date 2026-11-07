# High-Throughput Molecular Docking & Virtual Screening Pipeline

A modular, reproducible, high-throughput computational pipeline for virtual screening against kinase drug targets using AutoDock Vina, RDKit, Open Babel, and Meeko.

This pipeline automates the complete computer-aided drug design (CADD) screening workflow: receptor fetching and preparation, binding site grid-box definition, ligand 3D conformer optimization, multi-core parallel molecular docking, interaction fingerprinting, statistical enrichment analysis (ROC AUC, EF 1%), and 3D molecular visualization.

---

## Target System & Biology

- **Primary Target**: Human Epidermal Growth Factor Receptor (EGFR) Kinase Domain
- **Structure Reference**: PDB [1M17](https://www.rcsb.org/structure/1M17) (Chain A, resolution 2.60 Å)
- **Co-crystal Reference Ligand**: Erlotinib analogue (`AQ4` / PubChem CID 176870)
- **Key Hinge Residue**: **Met793** (Forms essential hydrogen bonds with ATP-competitive kinase inhibitors)
- **Gatekeeper Residue**: **Thr790** (Site of the clinically frequent T790M drug-resistance mutation)

---

## Pipeline Architecture

```
Target PDB (1M17)          SMILES Library (data/library.csv)
       │                                  │
[fetch_receptor]                  [clean_smiles]
       │                                  │
[clean_receptor]                  [smiles_to_pdbqt]
 (Chain A, no HETATM)               (ETKDGv3 + MMFF94 + Meeko)
       │                                  │
[prepare_receptor]                        │
 (Open Babel -> PDBQT)                    │
       │                                  │
[define_grid_box]                         │
 (Centroid of co-crystal AQ4)             │
       └──────────────────┬───────────────┘
                          │
                   [dock_ligand]
             (AutoDock Vina, Monte Carlo)
                          │
          [extract_interactions_distance]
           (H-bonds <= 3.5 Å, Hydrophobic <= 4.5 Å)
                          │
                 [analyze_results]
        (ROC AUC, EF 1%, Met793 Hinge Check)
                          │
             [generate_pymol_session]
             (PyMOL .pml / ChimeraX .cxc)
```

---

## Installation & Environment Setup

Because RDKit and Open Babel depend on compiled C++ libraries, use **Miniconda** or **Mamba**:

### 1. Create and Activate Environment
```bash
conda create -n docking python=3.11 -y
conda activate docking
```

### 2. Install Scientific Dependencies
```bash
# Core cheminformatics dependencies via conda-forge
conda install -c conda-forge rdkit openbabel -y

# Docking and data analysis dependencies via pip
pip install vina meeko biopython pandas numpy scipy scikit-learn matplotlib tqdm pytest
```

---

## One-Command Quickstart

### 1. Sanity Check Run (5 Molecules)
Validates the pipeline on active drugs (erlotinib, gefitinib) versus negative controls (aspirin, caffeine, ibuprofen):
```bash
python docking_pipeline.py --pdb-id 1M17 --library data/molecules_sanity.csv --output-dir results_egfr
```

### 2. Redocking Validation (Phase 2)
Docks the crystallographic ligand (AQ4) back into the binding pocket and evaluates heavy-atom RMSD without superimposition (target: $< 2.0\text{ \AA}$):
```bash
python validate_redock.py --pdb-id 1M17 --output-dir results_egfr
```

### 3. Generate Benchmark Screening Library (Phase 3)
Retrieves 100+ active inhibitors from the ChEMBL API and 400+ property-matched decoys from DUD-E:
```bash
python prepare_library.py --target-chembl-id CHEMBL203 --output-file data/library.csv --n-actives 100 --n-decoys 400
```

### 4. High-Throughput Parallel Virtual Screen (Phase 4)
Screens the full 500+ compound library across multiple CPU cores with automatic checkpointing and resume capability:
```bash
python docking_pipeline.py --pdb-id 1M17 --library data/library.csv --workers 4 --output-dir results_egfr
```
*If interrupted, add `--resume` to skip already-docked molecules.*

### 5. Statistical Enrichment & Hinge Validation (Phase 5)
Computes ROC AUC, Enrichment Factor (EF 1%, EF 5%), and generates publication-style figures:
```bash
python analyze_results.py --results results_egfr/docking_results.csv --library data/library.csv --output-dir results_egfr
```

### 6. 3D Molecular Visualisation (Phase 6)
Generates PyMOL and ChimeraX sessions highlighting top docking poses and Met793 contacts:
```bash
python generate_pymol_session.py --results-dir results_egfr --n-hits 3
pymol results_egfr/visualize_top_hits.pml
```

---

## Repository Structure

```
docking-pipeline/
├── data/
│   ├── molecules_sanity.csv        # 5-molecule sanity test set
│   └── library.csv                 # 509-molecule benchmark set (actives + decoys)
├── docs/
│   ├── prd.md                      # Product Requirements Document
│   ├── architecture.md             # System & Security Architecture
│   ├── rules.md                    # Scientific and coding guidelines
│   ├── design.md                   # UI & graphical design specifications
│   ├── task.md                     # Phase-by-phase task tracking
│   ├── memory.md                   # Persistent project context
│   └── VISUALISATION.md            # Detailed PyMOL & ChimeraX guide
├── tests/
│   ├── test_docking_pipeline.py     # Stages 1-8 unit tests
│   ├── test_validate_redock.py     # RMSD and redocking tests
│   ├── test_prepare_library.py     # Desalting and ChEMBL/DUD-E tests
│   ├── test_analyze_results.py     # ROC AUC and Enrichment tests
│   └── test_generate_pymol_session.py # 3D macro generator tests
├── docking_pipeline.py             # Main screening pipeline
├── validate_redock.py              # Redocking validation suite
├── prepare_library.py              # Library preparation generator
├── analyze_results.py              # Statistical evaluation suite
├── generate_pymol_session.py       # 3D session generator
├── environment.yml                 # Conda environment definition
└── requirements.txt                # Pip requirements specification
```

---

## Running the Unit Test Suite

Run the full project test suite (96 tests passing across all components):
```bash
pytest tests/ -v
```

---

## Honest Scientific Limitations & Caveats

1. **Scoring Function Limitations**: AutoDock Vina's empirical scoring function estimates binding free energy ($\Delta G$ in kcal/mol). While effective at rank-ordering compounds, Vina scores correlate imperfectly with experimental $\text{IC}_{50}$ or $K_i$ values.
2. **Receptor Rigidity**: The receptor backbone is held rigid during standard screening (`-xr` flag in Open Babel). Induced-fit conformational rearrangements upon ligand binding are not modeled.
3. **Distance-Based Interaction Detection**: The baseline interaction extraction uses geometric Euclidean distance cutoffs ($\le 3.5\text{ \AA}$ for H-bonds, $\le 4.5\text{ \AA}$ for hydrophobic contacts). While computationally efficient, it does not evaluate donor-hydrogen-acceptor angles.
4. **Decoy Assumption**: DUD-E decoys are presumed inactive based on topological disparity; occasional undetected active binders ("latent actives") may exist in decoy sets.
