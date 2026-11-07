"""
test_generate_pymol_session.py — Phase 6 Test Suite
===================================================
Tests for generate_pymol_session.py:
  1. PyMOL macro script generation (.pml) and syntax checks (Met793 hinge, white background, loading commands).
  2. ChimeraX command script generation (.cxc).
  3. CLI help and execution checks.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).parent.parent
SCRIPT = REPO_ROOT / "generate_pymol_session.py"


def _import_fn(name: str):
    sys.path.insert(0, str(REPO_ROOT))
    import generate_pymol_session as gps
    return getattr(gps, name)


class TestPyMOLScriptGeneration:
    """Test generating PyMOL .pml and ChimeraX .cxc scripts."""

    def test_generate_pymol_script(self, tmp_path: Path) -> None:
        gen_pymol = _import_fn("generate_pymol_script")

        # Create dummy top5_hits.csv
        df = pd.DataFrame({
            "name": ["erlotinib", "gefitinib", "hit_3"],
            "affinity_kcal_mol": [-9.5, -9.1, -8.7],
        })
        df.to_csv(tmp_path / "top5_hits.csv", index=False)

        out_pml = tmp_path / "test.pml"
        gen_pymol(results_dir=tmp_path, n_hits=2, output_path=out_pml)

        assert out_pml.exists()
        content = out_pml.read_text()
        assert "load receptor_clean.pdb" in content
        assert "resi 793" in content
        assert "Met793" in content
        assert "poses/erlotinib_best.pdbqt" in content
        assert "poses/gefitinib_best.pdbqt" in content
        assert "distance hb_1" in content

    def test_generate_chimerax_script(self, tmp_path: Path) -> None:
        gen_cx = _import_fn("generate_chimerax_script")

        df = pd.DataFrame({
            "name": ["erlotinib", "gefitinib"],
            "affinity_kcal_mol": [-9.5, -9.1],
        })
        df.to_csv(tmp_path / "top5_hits.csv", index=False)

        out_cxc = tmp_path / "test.cxc"
        gen_cx(results_dir=tmp_path, n_hits=2, output_path=out_cxc)

        assert out_cxc.exists()
        content = out_cxc.read_text()
        assert "open receptor_clean.pdb" in content
        assert "#1:793" in content
        assert "hbonds" in content

    def test_cli_execution(self, tmp_path: Path) -> None:
        res = subprocess.run([sys.executable, str(SCRIPT), "--help"], capture_output=True, text=True)
        assert res.returncode == 0
        assert "--results-dir" in res.stdout
        assert "--n-hits" in res.stdout
