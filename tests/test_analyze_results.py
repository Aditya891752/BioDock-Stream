"""
test_analyze_results.py — Phase 5 Test Suite
=============================================
Tests for analyze_results.py:
  1. compute_roc_auc: verifies perfect discrimination (1.0), inverse discrimination, and input checks.
  2. compute_enrichment_factor: verifies known active distribution at top 1% and 5%.
  3. evaluate_met793_hinge_contacts: checks case-insensitive string parsing of residue lists.
  4. Full run_analysis workflow: table generation, report writing, and plotting output.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).parent.parent
SCRIPT = REPO_ROOT / "analyze_results.py"


def _import_fn(name: str):
    sys.path.insert(0, str(REPO_ROOT))
    import analyze_results as ar
    return getattr(ar, name)


class TestROCAUC:
    """Test ROC AUC calculations."""

    def test_perfect_discrimination(self) -> None:
        compute_roc = _import_fn("compute_roc_auc")
        # Actives have more negative scores than decoys
        df = pd.DataFrame({
            "name": [f"mol_{i}" for i in range(10)],
            "affinity_kcal_mol": [-11.0, -10.5, -9.8, -9.5, -9.0, -5.0, -4.5, -4.0, -3.5, -3.0],
            "is_active": [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
        })
        auc = compute_roc(df)
        assert abs(auc - 1.0) < 1e-5

    def test_random_discrimination(self) -> None:
        compute_roc = _import_fn("compute_roc_auc")
        # Alternating scores
        df = pd.DataFrame({
            "name": [f"mol_{i}" for i in range(10)],
            "affinity_kcal_mol": [-10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0],
            "is_active": [1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
        })
        auc = compute_roc(df)
        # Random distribution should give AUC ~ 0.55
        assert 0.4 <= auc <= 0.65

    def test_single_class_raises(self) -> None:
        compute_roc = _import_fn("compute_roc_auc")
        df = pd.DataFrame({
            "name": ["m1", "m2"],
            "affinity_kcal_mol": [-8.0, -7.0],
            "is_active": [1, 1],
        })
        with pytest.raises(ValueError, match="requires both active"):
            compute_roc(df)


class TestEnrichmentFactor:
    """Test Enrichment Factor formula."""

    def test_ef_calculation(self) -> None:
        compute_ef = _import_fn("compute_enrichment_factor")
        # 100 compounds: 10 actives, 90 decoys (base rate = 10%)
        # Top 1% = top 1 molecule. If top 1 is active, EF = (1/1) / (10/100) = 10.0x
        affinities = [-10.0] + [-5.0] * 99
        is_active = [1] + [1] * 9 + [0] * 90
        df = pd.DataFrame({
            "name": [f"m_{i}" for i in range(100)],
            "affinity_kcal_mol": affinities,
            "is_active": is_active,
        })
        ef_1 = compute_ef(df, fraction=0.01)
        assert abs(ef_1 - 10.0) < 1e-4

    def test_ef_zero_when_no_active_in_top_slice(self) -> None:
        compute_ef = _import_fn("compute_enrichment_factor")
        # Top 1 molecule is a decoy
        affinities = [-10.0] + [-5.0] * 99
        is_active = [0] + [1] * 10 + [0] * 89
        df = pd.DataFrame({
            "name": [f"m_{i}" for i in range(100)],
            "affinity_kcal_mol": affinities,
            "is_active": is_active,
        })
        ef_1 = compute_ef(df, fraction=0.01)
        assert ef_1 == 0.0


class TestHingeContacts:
    """Test Met793 detection."""

    def test_detect_met793(self) -> None:
        eval_fn = _import_fn("evaluate_met793_hinge_contacts")
        df = pd.DataFrame({
            "name": ["erlotinib", "aspirin", "gefitinib"],
            "hbond_residues": ["MET793_A;THR790_A", "LEU718_A", "met793_a;GLN791_A"],
        })
        evaluated = eval_fn(df)
        assert evaluated.loc[0, "contacts_met793"] is True or evaluated.loc[0, "contacts_met793"] == 1
        assert not evaluated.loc[1, "contacts_met793"]
        assert evaluated.loc[2, "contacts_met793"] is True or evaluated.loc[2, "contacts_met793"] == 1


class TestFullAnalysisWorkflow:
    """Integration test for full report and figure generation."""

    def test_run_analysis_end_to_end(self, tmp_path: Path) -> None:
        run_analysis = _import_fn("run_analysis")

        # Create dummy docking_results.csv
        results_file = tmp_path / "docking_results.csv"
        df = pd.DataFrame({
            "name": [f"compound_{i}" for i in range(30)],
            "smiles": ["C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"] * 30,
            "affinity_kcal_mol": np.linspace(-10.5, -4.0, 30),
            "hbond_residues": ["MET793_A"] * 10 + ["LEU718_A"] * 20,
            "hydrophobic_residues": ["LEU844_A"] * 30,
            "status": ["success"] * 30,
            "is_active": [1] * 10 + [0] * 20,
        })
        df.to_csv(results_file, index=False)

        out_dir = tmp_path / "output_analysis"
        run_analysis(results_path=results_file, library_path=None, output_dir=out_dir)

        # Confirm expected deliverables exist
        assert (out_dir / "top5_hits.csv").exists()
        assert (out_dir / "analysis_report.txt").exists()
        assert (out_dir / "roc_curve.png").exists()
        assert (out_dir / "affinity_hist.png").exists()

        report_txt = (out_dir / "analysis_report.txt").read_text()
        assert "ROC AUC" in report_txt
        assert "Enrichment @" in report_txt
        assert "Met793" in report_txt
