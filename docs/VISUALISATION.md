# Visualisation Guide for EGFR Kinase (PDB: 1M17) Docking Poses

This guide explains how to inspect and render 3D binding poses of top screening hits using PyMOL and UCSF ChimeraX.

---

## 1. Quick Automated Setup (Recommended)

Generate visualization macro scripts dynamically based on your latest docking run:

```bash
python generate_pymol_session.py --results-dir results_egfr --n-hits 3
```

This generates:
- `results_egfr/visualize_top_hits.pml` (for PyMOL)
- `results_egfr/visualize_top_hits.cxc` (for ChimeraX)

---

## 2. Running in PyMOL

### Step-by-Step CLI Execution
Open your terminal in the project directory and run:
```bash
pymol results_egfr/visualize_top_hits.pml
```

Or open PyMOL GUI and type into the upper command line:
```
cd d:/BIOINOFRMATICS PROJECT/docking-pipeline/results_egfr
@visualize_top_hits.pml
```

### Key Biological Features Highlighted
1. **Receptor Backbone**: Displayed as a clean light-gray cartoon (`receptor`).
2. **Hinge Residue (Met793)**: Displayed as cyan sticks. Inhibitors form direct hydrogen bonds with the backbone nitrogen or carbonyl oxygen of Met793.
3. **Gatekeeper (Thr790)**: Displayed as yellow sticks. Controls access to the hydrophobic back-pocket (mutates to bulkier Met790 in clinical resistance).
4. **Docked Ligands**: Displayed in prominent stick representations with carbon atoms color-coded per hit.
5. **Hydrogen Bonds**: Automatically measured between ligand heteroatoms (N/O) and hinge residues within $\le 3.5\text{ \AA}$ (rendered as black dashed lines).

### Capturing Publication-Quality Figures in PyMOL
In the PyMOL command line, run:
```
set ray_shadow, 0
set antialias, 2
ray 2400, 1800
png top_hit_binding_mode.png, dpi=300
```

---

## 3. Running in UCSF ChimeraX

### Execution
Open UCSF ChimeraX, open the command line prompt at the bottom, and execute:
```
open results_egfr/visualize_top_hits.cxc
```

### Saving ChimeraX Images
```
save chimerax_top_hits.png width 2400 height 1800 supersample 3
```
