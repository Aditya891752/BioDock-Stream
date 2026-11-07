"""
generate_pymol_session.py — Automated PyMOL & ChimeraX Visualisation Script Generator
=====================================================================================
Phase 6 deliverable.

Scientific objective
--------------------
To visually validate docking results, inspect ligand-receptor binding modes,
and verify critical hydrogen bonds with the hinge residue (Met793) of EGFR kinase,
researchers use 3D molecular graphics tools like PyMOL or ChimeraX.

This script automatically generates:
  1. `visualize_top_hits.pml`: A PyMOL macro script that loads:
     - The clean receptor structure (`receptor_clean.pdb` or raw PDB)
     - The crystallographic reference ligand (AQ4 / erlotinib)
     - The top docked poses from the screening run (`poses/<name>_best.pdbqt`)
     - Pre-configured representation: cartoon protein, stick ligands, pocket surface,
       and dashed lines representing hydrogen bonds to Met793.
  2. `visualize_top_hits.cxc`: A ChimeraX command script providing the equivalent
     setup for UCSF ChimeraX users.

Usage
-----
python generate_pymol_session.py --results-dir results_egfr --n-hits 3
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger("generate_pymol_session")


def configure_logging() -> None:
    if logger.handlers:
        return
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("%(levelname)-8s %(message)s"))
    logger.addHandler(handler)


def generate_pymol_script(
    results_dir: Path,
    n_hits: int = 3,
    output_path: Optional[Path] = None,
) -> Path:
    """
    Generate a PyMOL macro script (.pml).
    """
    if output_path is None:
        output_path = results_dir / "visualize_top_hits.pml"

    # Identify top hits from top5_hits.csv or docking_results.csv
    top_hits = []
    top5_csv = results_dir / "top5_hits.csv"
    docking_csv = results_dir / "docking_results.csv"

    if top5_csv.exists():
        df = pd.read_csv(top5_csv)
        top_hits = df["name"].head(n_hits).tolist()
    elif docking_csv.exists():
        df = pd.read_csv(docking_csv)
        succ = df[df["status"] == "success"].sort_values("affinity_kcal_mol")
        top_hits = succ["name"].head(n_hits).tolist()
    else:
        # Placeholder names if run before docking completes
        top_hits = ["erlotinib", "gefitinib", "hit_3"]

    lines = [
        "# ===========================================================================",
        "# PyMOL Visualization Script for EGFR (PDB 1M17) Docking Hits",
        "# Generated automatically by Bioinformatics Docking Pipeline (Phase 6)",
        "# ===========================================================================",
        "",
        "# 1. Initialize environment",
        "reinitialize",
        "set bg_rgb, [1, 1, 1]              # Clean white background for clarity",
        "set ray_shadow, 0                  # Clean publication view without dark shadows",
        "set orthoscopic, on                # True perspective without distortion",
        "",
        "# 2. Load receptor",
        "load receptor_clean.pdb, receptor",
        "hide everything, receptor",
        "show cartoon, receptor",
        "color gray80, receptor             # Neutral light gray for receptor backbone",
        "",
        "# 3. Highlight key binding pocket residues",
        "# Met793: Key hinge residue forming canonical H-bonds with active inhibitors",
        "# Thr790: Gatekeeper residue (mutates to Met790 in drug resistance)",
        "# Lys745: Catalytic lysine in beta3 strand",
        "select hinge, receptor and resi 793 and resn MET",
        "select gatekeeper, receptor and resi 790 and resn THR",
        "select catalytic_lys, receptor and resi 745 and resn LYS",
        "select pocket_res, receptor and (resi 718-726 or resi 743-747 or resi 790-798 or resi 844-846)",
        "",
        "show sticks, hinge",
        "show sticks, gatekeeper",
        "color cyan, hinge and elem C",
        "color yellow, gatekeeper and elem C",
        "label hinge and name CA, '\"Met793 (Hinge)\"'",
        "label gatekeeper and name CA, '\"Thr790 (Gatekeeper)\"'",
        "set label_color, black",
        "set label_size, 14",
        "",
        "# 4. Load co-crystal reference ligand if present",
        "load 1M17_raw.pdb, raw_pdb",
        "select crystal_aq4, raw_pdb and resn AQ4",
        "show sticks, crystal_aq4",
        "color magenta, crystal_aq4 and elem C",
        "set stick_radius, 0.22, crystal_aq4",
        "hide everything, raw_pdb and not crystal_aq4",
        "",
    ]

    # Load top docked poses
    colors = ["tv_blue", "tv_green", "orange", "marine", "firebrick"]
    lines.append("# 5. Load top docked poses")

    for i, name in enumerate(top_hits):
        col = colors[i % len(colors)]
        lines.extend([
            f"load poses/{name}_best.pdbqt, pose_{i+1}_{name}",
            f"show sticks, pose_{i+1}_{name}",
            f"color {col}, pose_{i+1}_{name} and elem C",
            f"set stick_radius, 0.28, pose_{i+1}_{name}",
            f"# Compute hydrogen bonds to hinge (distance <= 3.5 Å)",
            f"distance hb_{i+1}, pose_{i+1}_{name} and (elem N or elem O), hinge and (elem N or elem O), 3.5",
            f"color dash, black",
            f"set dash_gap, 0.2",
            f"set dash_radius, 0.05",
            "",
        ])

    lines.extend([
        "# 6. Orient camera on the ATP-binding pocket",
        "zoom hinge, 12",
        "clip near, 2",
        "",
        "# 7. Render high-resolution PNG image",
        "# png results_egfr/top_hits_binding_mode.png, width=1800, height=1400, dpi=300, ray=1",
        "print 'Visualization loaded successfully. Inspect pocket and hinge H-bonds.'",
    ])

    script_text = "\n".join(lines)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(script_text, encoding="utf-8")
    logger.info("Generated PyMOL script: %s", output_path)
    return output_path


def generate_chimerax_script(
    results_dir: Path,
    n_hits: int = 3,
    output_path: Optional[Path] = None,
) -> Path:
    """
    Generate a UCSF ChimeraX command script (.cxc).
    """
    if output_path is None:
        output_path = results_dir / "visualize_top_hits.cxc"

    top_hits = []
    top5_csv = results_dir / "top5_hits.csv"
    if top5_csv.exists():
        df = pd.read_csv(top5_csv)
        top_hits = df["name"].head(n_hits).tolist()
    else:
        top_hits = ["erlotinib", "gefitinib", "hit_3"]

    lines = [
        "# ===========================================================================",
        "# ChimeraX Visualization Script for EGFR (PDB 1M17) Docking Hits",
        "# ===========================================================================",
        "",
        "# 1. Initialize environment",
        "close session",
        "set bgColor white",
        "",
        "# 2. Open receptor and set representation",
        "open receptor_clean.pdb",
        "hide #1 atoms",
        "show #1 cartoon",
        "color #1 lightgray",
        "",
        "# 3. Highlight key residues",
        "show #1:793,790 atoms,bonds",
        "color #1:793 cyan",
        "color #1:790 goldenrod",
        "label #1:793 text 'Met793' color black",
        "label #1:790 text 'Thr790' color black",
        "",
    ]

    for i, name in enumerate(top_hits, start=2):
        lines.extend([
            f"open poses/{name}_best.pdbqt",
            f"show #{i} sticks",
            f"# Calculate H-bonds between ligand #{i} and receptor #1",
            f"hbonds #{i} makePseudobonds true color black",
            "",
        ])

    lines.extend([
        "# View orientation",
        "view #1:793 distance 20",
        "# save results_egfr/chimerax_top_hits.png width 1800 height 1400",
    ])

    script_text = "\n".join(lines)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(script_text, encoding="utf-8")
    logger.info("Generated ChimeraX script: %s", output_path)
    return output_path


def main() -> None:
    configure_logging()
    parser = argparse.ArgumentParser(
        prog="generate_pymol_session",
        description="Generate PyMOL (.pml) and ChimeraX (.cxc) scripts to visualize top docking hits.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results_egfr"),
        help="Path to results directory containing receptor and poses.",
    )
    parser.add_argument(
        "--n-hits",
        type=int,
        default=3,
        help="Number of top hits to include in the visual session (default: 3).",
    )
    args = parser.parse_args()
    generate_pymol_script(args.results_dir, args.n_hits)
    generate_chimerax_script(args.results_dir, args.n_hits)


if __name__ == "__main__":
    main()
