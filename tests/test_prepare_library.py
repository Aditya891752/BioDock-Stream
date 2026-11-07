"""
test_prepare_library.py — Phase 3 Test Suite
=============================================
Tests for prepare_library.py:
  1. clean_smiles: Desalting, canonicalization, rejection of malformed inputs.
  2. Fallback actives: Validation that curated drugs parse cleanly.
  3. Decoy generation: Property-matched, non-kinase diverse decoys.
  4. build_screening_library: Output schema, row counts, deterministic shuffling.
  5. CLI help and argument parsing.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).parent.parent
SCRIPT = REPO_ROOT / "prepare_library.py"


def _import_fn(name: str):
    sys.path.insert(0, str(REPO_ROOT))
    import prepare_library as pl
    return getattr(pl, name)


class TestCleanSMILES:
    """Test chemical structure cleaning and desalting."""

    def test_clean_standard_molecule(self) -> None:
        clean = _import_fn("clean_smiles")
        # Aspirin
        res = clean("CC(=O)Oc1ccccc1C(=O)O")
        assert res is not None
        assert "c1ccccc1" in res or "C1=CC=CC=C1" in res

    def test_salt_stripping(self) -> None:
        clean = _import_fn("clean_smiles")
        # Molecule with sodium or chloride salt fragment: CC(=O)O.[Na+]
        res = clean("CC(=O)O.[Na+]")
        assert res is not None
        assert "." not in res  # Salt dot separated component stripped
        assert "CC(=O)O" in res or "C(=O)O" in res

    def test_invalid_smiles_rejected(self) -> None:
        clean = _import_fn("clean_smiles")
        assert clean("NOT_A_SMILES") is None
        assert clean("C1(=O") is None  # Unbalanced parentheses
        assert clean("") is None
        assert clean(None) is None


class TestActivesAndDecoys:
    """Test retrieval and assembly logic."""

    def test_fallback_actives_structure(self) -> None:
        get_fallbacks = _import_fn("_get_fallback_actives")
        actives = get_fallbacks(limit=25)
        assert len(actives) == 25
        assert all("name" in a and "smiles" in a for a in actives)
        assert any("erlotinib" in a["name"] for a in actives)

    def test_synthetic_decoys(self) -> None:
        gen_decoys = _import_fn("_generate_property_matched_decoys")
        decoys = gen_decoys(count=30)
        assert len(decoys) == 30
        assert all(d["name"].startswith("DUD_E_decoy_") for d in decoys)
        # Verify valid SMILES string produced
        clean = _import_fn("clean_smiles")
        assert all(clean(d["smiles"]) is not None for d in decoys)


class TestBuildLibrary:
    """Test full library compilation."""

    def test_build_library_schema_and_counts(self, tmp_path: Path, monkeypatch) -> None:
        import prepare_library as pl
        # Mock network functions with deterministic fallbacks to guarantee fast, offline execution
        monkeypatch.setattr(pl, "fetch_chembl_actives", lambda target_chembl_id, max_records: pl._get_fallback_actives(max_records))
        monkeypatch.setattr(pl, "fetch_dude_decoys", lambda target, max_records: pl._generate_property_matched_decoys(max_records))

        build_fn = pl.build_screening_library
        out_csv = tmp_path / "test_lib.csv"

        df = build_fn(
            n_actives=20,
            n_decoys=50,
            output_path=out_csv,
            seed=42,
        )

        assert out_csv.exists()
        assert len(df) == 70
        assert set(df.columns) == {"name", "smiles", "is_active"}

        # Check binary active labels
        n_act = (df["is_active"] == 1).sum()
        n_dec = (df["is_active"] == 0).sum()
        assert n_act == 20
        assert n_dec == 50

    def test_shuffling_is_deterministic(self, monkeypatch) -> None:
        import prepare_library as pl
        monkeypatch.setattr(pl, "fetch_chembl_actives", lambda target_chembl_id, max_records: pl._get_fallback_actives(max_records))
        monkeypatch.setattr(pl, "fetch_dude_decoys", lambda target, max_records: pl._generate_property_matched_decoys(max_records))

        build_fn = pl.build_screening_library
        df1 = build_fn(n_actives=10, n_decoys=20, seed=123)
        df2 = build_fn(n_actives=10, n_decoys=20, seed=123)
        pd.testing.assert_frame_equal(df1, df2)


class TestCLI:
    """Test CLI commands."""

    def test_cli_help(self) -> None:
        res = subprocess.run([sys.executable, str(SCRIPT), "--help"], capture_output=True, text=True)
        assert res.returncode == 0
        assert "--target-chembl-id" in res.stdout
        assert "--n-actives" in res.stdout
        assert "--n-decoys" in res.stdout
        assert "--output-file" in res.stdout
