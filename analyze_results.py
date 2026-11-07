"""
analyze_results.py — Post-Docking Analysis, Enrichment & Interaction Suite
==========================================================================
Phase 5 deliverable.

Scientific objective
--------------------
A raw list of AutoDock Vina affinities (kcal/mol) is insufficient to judge
whether a virtual screen succeeded. This analysis suite evaluates the results
using rigorous statistical, graphical, and biological metrics:

1. Enrichment & Statistical Discrimination:
   - ROC AUC (Receiver Operating Characteristic Area Under the Curve):
     Measures how effectively the docking scores rank true active inhibitors
     above decoys (0.5 = random chance, 1.0 = perfect ranking).
     *Note on sign*: More negative docking scores indicate stronger predicted
     binding, so we rank by ascending affinity (or negate score for ROC).
   - EF 1% (Enrichment Factor at top 1%):
     Measures how much the top 1% of the ranked library is enriched with active
     molecules compared to random selection:
       EF_1% = (Actives in top 1% / Total compounds in top 1%) /
               (Total Actives in library / Total compounds in library)

2. Biological Hinge-Residue Validation:
   - Known ATP-competitive EGFR kinase inhibitors (erlotinib, gefitinib,
     osimertinib) depend on a key hydrogen bond to the hinge backbone nitrogen
     of Met793. We inspect whether the highest-ranked molecules form contacts
     with Met793.

3. Visualizations (following design rules: no gradients, clean scientific style):
   - ROC curve (single accent color, reference diagonal line).
   - Affinity distribution histogram (contrasting actives vs decoys).

4. Reporting:
   - top5_hits.csv: Top 5 lead compounds with scores, SMILES, and contacts.
   - analysis_report.txt: Complete audit report with all numbers.

Usage
-----
python analyze_results.py --results results_egfr/docking_results.csv --library data/library.csv --output-dir results_egfr
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for headless generation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

logger = logging.getLogger("analyze_results")


# ===========================================================================
# LOGGING SETUP
# ===========================================================================

def configure_logging(verbose: bool = False) -> None:
    if logger.handlers:
        return
    level = logging.DEBUG if verbose else logging.INFO
    logger.setLevel(level)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter("%(levelname)-8s %(message)s"))
    logger.addHandler(handler)


# ===========================================================================
# 1. ENRICHMENT & ROC METRICS
# ===========================================================================

def compute_roc_auc(df: pd.DataFrame) -> float:
    """
    Calculate ROC AUC for ranking actives versus decoys by docking affinity.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing 'affinity_kcal_mol' and 'is_active'.
        Only successfully docked rows with non-null affinity are evaluated.

    Returns
    -------
    float
        ROC AUC value between 0.0 and 1.0.
    """
    valid = df.dropna(subset=["affinity_kcal_mol", "is_active"]).copy()
    if len(valid) == 0:
        raise ValueError("No valid rows with both affinity and is_active found.")

    y_true = valid["is_active"].astype(int).values
    if len(np.unique(y_true)) < 2:
        raise ValueError("ROC AUC requires both active (1) and decoy (0) samples.")

    # More negative score = stronger binding = higher rank.
    # sklearn assumes higher numerical score = positive class, so we negate.
    y_scores = -valid["affinity_kcal_mol"].values
    return float(roc_auc_score(y_true, y_scores))


def compute_enrichment_factor(df: pd.DataFrame, fraction: float = 0.01) -> float:
    """
    Compute Enrichment Factor at a specified top percentile (default: 1% / 0.01).

    Formula:
      EF_x% = (a / n) / (A / N)
      where:
        a = count of active molecules in top x%
        n = total compounds in top x% (n = ceil(fraction * N))
        A = total count of active molecules in library
        N = total compounds in library

    Parameters
    ----------
    df : pd.DataFrame
        Ranked or unranked DataFrame with 'affinity_kcal_mol' and 'is_active'.
    fraction : float
        Screening cutoff fraction (e.g. 0.01 for 1%).

    Returns
    -------
    float
        Enrichment factor score.
    """
    valid = df.dropna(subset=["affinity_kcal_mol", "is_active"]).copy()
    if len(valid) == 0:
        return 0.0

    # Sort ascending: most negative affinity first
    ranked = valid.sort_values(by="affinity_kcal_mol", ascending=True).reset_index(drop=True)

    N = len(ranked)
    A = int(ranked["is_active"].sum())
    if A == 0 or N == 0:
        return 0.0

    n = max(1, int(np.ceil(fraction * N)))
    top_slice = ranked.iloc[:n]
    a = int(top_slice["is_active"].sum())

    base_rate = A / N
    top_rate = a / n
    return float(top_rate / base_rate)


# ===========================================================================
# 2. PLOTTING (Scientific, High-Contrast Style — Rule-Compliant)
# ===========================================================================

def plot_roc_curve(df: pd.DataFrame, output_path: Path) -> None:
    """
    Plot clean ROC curve using solid accent color without gradients or gimmicks.
    """
    valid = df.dropna(subset=["affinity_kcal_mol", "is_active"]).copy()
    y_true = valid["is_active"].astype(int).values
    y_scores = -valid["affinity_kcal_mol"].values

    fpr, tpr, _ = roc_curve(y_true, y_scores)
    auc_val = roc_auc_score(y_true, y_scores)

    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)

    # Styling per design.md (clean, high contrast, #2563EB accent)
    ax.plot(fpr, tpr, color="#2563EB", lw=2.0, label=f"AutoDock Vina (AUC = {auc_val:.3f})")
    ax.plot([0, 1], [0, 1], color="#9C9C94", lw=1.2, linestyle="--", label="Random Selection (AUC = 0.500)")

    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.02])
    ax.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=11, fontweight="medium")
    ax.set_ylabel("True Positive Rate (Sensitivity)", fontsize=11, fontweight="medium")
    ax.set_title("EGFR Virtual Screening — ROC Curve", fontsize=13, fontweight="bold", pad=12)
    ax.legend(loc="lower right", frameon=True, framealpha=0.9, edgecolor="#E4E4E0")

    # Clean borders
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#6B6B63")
    ax.spines["bottom"].set_color("#6B6B63")
    ax.grid(color="#E4E4E0", linestyle=":", linewidth=0.8, alpha=0.7)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    logger.info("Saved ROC curve: %s", output_path)


def plot_affinity_distribution(df: pd.DataFrame, output_path: Path) -> None:
    """
    Plot histogram of docking affinities comparing actives vs decoys.
    """
    valid = df.dropna(subset=["affinity_kcal_mol", "is_active"]).copy()

    actives = valid[valid["is_active"] == 1]["affinity_kcal_mol"].values
    decoys = valid[valid["is_active"] == 0]["affinity_kcal_mol"].values

    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)

    bins = np.linspace(
        min(valid["affinity_kcal_mol"].min(), -12.0),
        max(valid["affinity_kcal_mol"].max(), -2.0),
        30,
    )

    ax.hist(
        decoys,
        bins=bins,
        alpha=0.65,
        color="#6B6B63",
        label=f"Decoys (N={len(decoys)})",
        edgecolor="#1A1A18",
        linewidth=0.5,
    )
    ax.hist(
        actives,
        bins=bins,
        alpha=0.85,
        color="#2563EB",
        label=f"Actives (N={len(actives)})",
        edgecolor="#1D4ED8",
        linewidth=0.5,
    )

    ax.set_xlabel("Predicted Binding Affinity (kcal/mol)", fontsize=11, fontweight="medium")
    ax.set_ylabel("Molecule Count", fontsize=11, fontweight="medium")
    ax.set_title("Affinity Distribution: Actives vs Decoys", fontsize=13, fontweight="bold", pad=12)
    ax.legend(loc="upper left", frameon=True, edgecolor="#E4E4E0")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#6B6B63")
    ax.spines["bottom"].set_color("#6B6B63")
    ax.grid(color="#E4E4E0", linestyle=":", linewidth=0.8, alpha=0.7)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    logger.info("Saved affinity distribution: %s", output_path)


# ===========================================================================
# 3. INTERACTION & HINGE CHECK
# ===========================================================================

def evaluate_met793_hinge_contacts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add boolean column 'contacts_met793' checking for H-bonds with Met793.
    """
    df_eval = df.copy()

    def has_met793(val) -> bool:
        if not val or not isinstance(val, str):
            return False
        return "MET793" in val.upper()

    df_eval["contacts_met793"] = df_eval["hbond_residues"].apply(has_met793)
    return df_eval


# ===========================================================================
# 4. REPORT GENERATOR
# ===========================================================================

def run_analysis(
    results_path: Path,
    library_path: Optional[Path],
    output_dir: Path,
) -> None:
    """
    Execute full analysis pipeline:
      - Merge results with active/decoy ground-truth labels if separate.
      - Calculate ROC AUC and EF 1%.
      - Check hinge contacts.
      - Write top5_hits.csv.
      - Save ROC curve and histogram plots.
      - Generate analysis_report.txt.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    configure_logging()

    logger.info("Loading docking results from: %s", results_path)
    res_df = pd.read_csv(results_path)

    # Merge is_active label from library CSV if not present in docking results
    if "is_active" not in res_df.columns and library_path and library_path.exists():
        logger.info("Merging ground truth is_active column from %s", library_path)
        lib_df = pd.read_csv(library_path)
        if "is_active" in lib_df.columns:
            res_df = res_df.merge(lib_df[["name", "is_active"]], on="name", how="left")

    # Filter successful rows
    success_df = res_df[res_df["status"] == "success"].copy()
    if len(success_df) == 0:
        raise ValueError("Results file contains zero successful docking entries.")

    success_df = evaluate_met793_hinge_contacts(success_df)
    sorted_df = success_df.sort_values(by="affinity_kcal_mol", ascending=True).reset_index(drop=True)

    # Top 5 hits table
    top5_cols = [c for c in ["name", "affinity_kcal_mol", "is_active", "contacts_met793", "hbond_residues", "smiles"] if c in sorted_df.columns]
    top5_df = sorted_df.head(5)[top5_cols]
    top5_path = output_dir / "top5_hits.csv"
    top5_df.to_csv(top5_path, index=False)
    logger.info("Saved top 5 hits table: %s", top5_path)

    # Calculate statistics if is_active is available
    has_labels = "is_active" in sorted_df.columns and sorted_df["is_active"].dropna().nunique() > 1

    if has_labels:
        auc_score = compute_roc_auc(sorted_df)
        ef_1pct = compute_enrichment_factor(sorted_df, fraction=0.01)
        ef_5pct = compute_enrichment_factor(sorted_df, fraction=0.05)

        roc_plot_path = output_dir / "roc_curve.png"
        plot_roc_curve(sorted_df, roc_plot_path)

        hist_plot_path = output_dir / "affinity_hist.png"
        plot_affinity_distribution(sorted_df, hist_plot_path)
    else:
        auc_score = None
        ef_1pct = None
        ef_5pct = None
        logger.warning("No ground-truth 'is_active' labels found; ROC and EF calculation omitted.")

    # Generate comprehensive text report
    n_total = len(res_df)
    n_success = len(sorted_df)
    n_met793 = sorted_df["contacts_met793"].sum()
    pct_met793 = (n_met793 / n_success) * 100 if n_success > 0 else 0.0

    lines = [
        "=" * 60,
        "VIRTUAL SCREENING ANALYSIS & ENRICHMENT REPORT",
        "=" * 60,
        f"Input Results File : {results_path}",
        f"Total Molecules    : {n_total}",
        f"Successful Docks   : {n_success} ({n_success/n_total*100:.1f}%)",
        "",
        "Enrichment & Statistical Performance",
        "-" * 40,
    ]

    if has_labels:
        lines += [
            f"  ROC AUC            : {auc_score:.3f} (Random = 0.500, Ideal = 1.000)",
            f"  Enrichment @ 1%    : {ef_1pct:.2f}x",
            f"  Enrichment @ 5%    : {ef_5pct:.2f}x",
            "",
            "  Interpretation:",
            (
                "  * High discriminatory performance (AUC > 0.70)."
                if auc_score >= 0.70
                else "  * Moderate/Low enrichment (AUC < 0.70). Review grid box bounds."
            ),
        ]
    else:
        lines += [
            "  ROC AUC            : N/A (unlabeled screening library)",
            "  Enrichment Factor  : N/A",
        ]

    lines += [
        "",
        "Biological Hinge (Met793) Interaction Validation",
        "-" * 40,
        f"  Poses contacting Met793 : {n_met793} / {n_success} ({pct_met793:.1f}%)",
    ]

    lines += [
        "",
        "Top 5 Ranked Hits",
        "-" * 40,
    ]
    for i, row in top5_df.iterrows():
        act_str = f" [Active={row['is_active']}]" if "is_active" in row else ""
        lines.append(
            f"  {i+1}. {row['name']:<24} {row['affinity_kcal_mol']:+6.2f} kcal/mol | "
            f"Met793: {'YES' if row.get('contacts_met793') else 'no ':3s}{act_str}"
        )

    lines += [
        "",
        "=" * 60,
        "DELIVERABLES GENERATED",
        "=" * 60,
        f"  1. Top 5 hits table : {top5_path.name}",
    ]
    if has_labels:
        lines += [
            f"  2. ROC curve plot   : roc_curve.png",
            f"  3. Affinity histogram: affinity_hist.png",
        ]
    lines.append("=" * 60)

    report_text = "\n".join(lines)
    report_path = output_dir / "analysis_report.txt"
    report_path.write_text(report_text, encoding="utf-8")
    logger.info("Saved analysis report: %s", report_path)
    print("\n" + report_text)


# ===========================================================================
# CLI ARGUMENT PARSER
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="analyze_results",
        description="Compute ROC AUC, EF 1%, Met793 hinge contacts, and plots from docking results.",
    )
    parser.add_argument(
        "--results",
        type=Path,
        required=True,
        help="Path to docking_results.csv output from docking_pipeline.py",
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=None,
        help="Optional path to library.csv containing ground truth is_active column.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results_egfr"),
        help="Directory to save report, plots, and top hits CSV (default: results_egfr).",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    run_analysis(args.results, args.library, args.output_dir)


if __name__ == "__main__":
    main()
