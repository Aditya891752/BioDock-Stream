"""
test_docking_pipeline.py — Phase 0 + Phase 1 tests
====================================================
Tests are split into two classes:

  TestCLI*       — test the CLI parser and argument handling.
                   These need NO scientific libraries.
                   They run on any Python installation.

  TestParsing*   — test helper functions that parse PDB/PDBQT text.
                   These need NO scientific libraries either.

  TestGridBox*   — test define_grid_box() with synthetic PDB content.
                   Needs: biopython, numpy.

  TestInteractions — test extract_interactions_distance() end-to-end.
                   Needs: numpy, biopython.

  TestWriteResults — test write_results() CSV output.
                   Needs: pandas.

  TestSMILES*    — test smiles_to_pdbqt() with real molecules.
                   Needs: rdkit, meeko. SKIPPED if not installed.

  TestDocking*   — test dock_ligand() with a real receptor.
                   Needs: vina, rdkit, meeko, biopython, openbabel.
                   SKIPPED if not installed — these require conda env.

Run with:
    pytest tests/test_docking_pipeline.py -v
    pytest tests/test_docking_pipeline.py -v -k "not Docking"   # skip slow tests
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).parent.parent
PIPELINE  = REPO_ROOT / "docking_pipeline.py"
SANITY_CSV = REPO_ROOT / "data" / "molecules_sanity.csv"


# ===========================================================================
# Phase 0 tests — CLI, --help, argument parsing (no scientific deps)
# ===========================================================================

class TestCLIHelp:
    """The --help flag must exit 0 and document all key arguments."""

    def test_help_exits_zero(self) -> None:
        """
        --help must exit with code 0.
        A non-zero exit means a syntax error or import crash at startup.
        """
        result = subprocess.run(
            [sys.executable, str(PIPELINE), "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            f"--help exited {result.returncode}.\n"
            f"stderr: {result.stderr}\nstdout: {result.stdout}"
        )

    def test_help_mentions_pdb_id(self) -> None:
        result = subprocess.run(
            [sys.executable, str(PIPELINE), "--help"],
            capture_output=True, text=True,
        )
        assert "--pdb-id" in result.stdout

    def test_help_mentions_library(self) -> None:
        result = subprocess.run(
            [sys.executable, str(PIPELINE), "--help"],
            capture_output=True, text=True,
        )
        assert "--library" in result.stdout

    def test_help_mentions_exhaustiveness(self) -> None:
        result = subprocess.run(
            [sys.executable, str(PIPELINE), "--help"],
            capture_output=True, text=True,
        )
        assert "--exhaustiveness" in result.stdout

    def test_help_mentions_seed(self) -> None:
        result = subprocess.run(
            [sys.executable, str(PIPELINE), "--help"],
            capture_output=True, text=True,
        )
        assert "--seed" in result.stdout

    def test_help_mentions_resume(self) -> None:
        result = subprocess.run(
            [sys.executable, str(PIPELINE), "--help"],
            capture_output=True, text=True,
        )
        assert "--resume" in result.stdout

    def test_no_args_exits_2(self) -> None:
        """Missing required args → argparse exits 2."""
        result = subprocess.run(
            [sys.executable, str(PIPELINE)],
            capture_output=True, text=True,
        )
        assert result.returncode == 2

    def test_missing_library_exits_2(self) -> None:
        result = subprocess.run(
            [sys.executable, str(PIPELINE), "--pdb-id", "1M17"],
            capture_output=True, text=True,
        )
        assert result.returncode == 2


class TestCLIArgumentParsing:
    """Verify that arguments are parsed and reach run_pipeline()."""

    def _run(self, *extra_args: str, tmp_path: Path) -> subprocess.CompletedProcess:
        """Helper to run the pipeline with minimal required args."""
        return subprocess.run(
            [
                sys.executable, str(PIPELINE),
                "--pdb-id", "1M17",
                "--library", str(SANITY_CSV),
                "--output-dir", str(tmp_path / "results"),
                *extra_args,
            ],
            capture_output=True, text=True,
        )

    def test_pdb_id_logged(self, tmp_path: Path) -> None:
        """PDB ID passed on CLI must appear in stdout."""
        r = self._run(tmp_path=tmp_path)
        assert "1M17" in r.stdout, f"stdout: {r.stdout[:500]}"

    def test_default_seed_42(self, tmp_path: Path) -> None:
        """Default seed must be 42 as documented in rules.md."""
        r = self._run(tmp_path=tmp_path)
        assert "42" in r.stdout

    def test_custom_chain(self, tmp_path: Path) -> None:
        """--chain B must be reflected in the log."""
        r = self._run("--chain", "B", tmp_path=tmp_path)
        assert "B" in r.stdout

    def test_output_dir_created(self, tmp_path: Path) -> None:
        """Output directory must be created automatically."""
        out = tmp_path / "auto_created_dir"
        assert not out.exists()
        self._run("--output-dir", str(out), tmp_path=tmp_path)
        assert out.exists(), "run_pipeline() must create the output directory"

    def test_log_file_created(self, tmp_path: Path) -> None:
        """run.log must be written inside the output directory."""
        out = tmp_path / "results"
        self._run("--output-dir", str(out), tmp_path=tmp_path)
        assert (out / "run.log").exists()


# ===========================================================================
# Sanity CSV tests (no scientific deps)
# ===========================================================================

class TestSanityCSV:
    """The five-molecule sanity CSV must have the correct structure."""

    def test_exists(self) -> None:
        assert SANITY_CSV.exists(), f"Missing: {SANITY_CSV}"

    def test_header(self) -> None:
        lines = SANITY_CSV.read_text().strip().splitlines()
        assert lines[0] == "name,smiles,is_active"

    def test_five_data_rows(self) -> None:
        lines = [
            l for l in SANITY_CSV.read_text().strip().splitlines()
            if l and not l.startswith("name")
        ]
        assert len(lines) == 5, f"Expected 5 molecules, got {len(lines)}"

    def test_contains_erlotinib(self) -> None:
        assert "erlotinib" in SANITY_CSV.read_text()

    def test_contains_aspirin(self) -> None:
        assert "aspirin" in SANITY_CSV.read_text()

    def test_erlotinib_is_active(self) -> None:
        """erlotinib must have is_active=1 — it's the positive control."""
        for line in SANITY_CSV.read_text().splitlines():
            if line.startswith("erlotinib"):
                assert line.endswith(",1"), (
                    f"erlotinib row must end with ',1' (is_active=1). Got: {line}"
                )

    def test_aspirin_is_inactive(self) -> None:
        """aspirin must have is_active=0 — it's the negative control."""
        for line in SANITY_CSV.read_text().splitlines():
            if line.startswith("aspirin"):
                assert line.endswith(",0"), (
                    f"aspirin row must end with ',0' (is_active=0). Got: {line}"
                )


# ===========================================================================
# _parse_atom_records tests (no scientific deps — pure string parsing)
# ===========================================================================

class TestParseAtomRecords:
    """
    Test _parse_atom_records() with synthetic PDB text.
    This function has no external dependencies.
    """

    def _import(self):
        """Import the function under test from the pipeline module."""
        sys.path.insert(0, str(REPO_ROOT))
        from docking_pipeline import _parse_atom_records
        return _parse_atom_records

    # Minimal valid PDB ATOM line (80 characters)
    # Columns: record(1-6) serial(7-11) name(13-16) resname(18-20) chain(22)
    #          resnum(23-26) x(31-38) y(39-46) z(47-54) element(77-78)
    SAMPLE_PDB = textwrap.dedent("""\
        ATOM      1  N   MET A   1      11.104  13.207  12.994  1.00  0.00           N
        ATOM      2  CA  MET A   1      11.639  13.051  11.651  1.00  0.00           C
        ATOM      3  O   MET A   1      10.819  12.398  10.869  1.00  0.00           O
        HETATM  101  C1  LIG A 200       5.000   6.000   7.000  1.00  0.00           C
        REMARK This line should be ignored
    """)

    def test_counts_atom_and_hetatm(self) -> None:
        parse = self._import()
        atoms = parse(self.SAMPLE_PDB)
        # 3 ATOM + 1 HETATM = 4 total; REMARK is ignored
        assert len(atoms) == 4

    def test_coordinates_parsed_correctly(self) -> None:
        parse = self._import()
        atoms = parse(self.SAMPLE_PDB)
        n_atom = atoms[0]
        assert abs(n_atom["x"] - 11.104) < 0.001
        assert abs(n_atom["y"] - 13.207) < 0.001
        assert abs(n_atom["z"] - 12.994) < 0.001

    def test_residue_name_parsed(self) -> None:
        parse = self._import()
        atoms = parse(self.SAMPLE_PDB)
        assert atoms[0]["resname"] == "MET"

    def test_chain_parsed(self) -> None:
        parse = self._import()
        atoms = parse(self.SAMPLE_PDB)
        assert atoms[0]["chain"] == "A"

    def test_resnum_parsed(self) -> None:
        parse = self._import()
        atoms = parse(self.SAMPLE_PDB)
        assert atoms[0]["resnum"] == 1
        assert atoms[3]["resnum"] == 200

    def test_element_from_column(self) -> None:
        """Element should be read from column 77-78 when present."""
        parse = self._import()
        atoms = parse(self.SAMPLE_PDB)
        assert atoms[0]["element"] == "N"
        assert atoms[1]["element"] == "C"
        assert atoms[2]["element"] == "O"

    def test_element_fallback_from_atom_name(self) -> None:
        """When element column is absent, fall back to atom name."""
        parse = self._import()
        # Line without element column (only 54 characters)
        short_line = "ATOM      1  CA  ALA A   1      1.000   2.000   3.000\n"
        atoms = parse(short_line)
        assert len(atoms) == 1
        # "CA" → first alphabetic char is "C"
        assert atoms[0]["element"] == "C"

    def test_empty_text_returns_empty_list(self) -> None:
        parse = self._import()
        assert parse("") == []
        assert parse("REMARK just a comment\nEND\n") == []

    def test_malformed_coordinate_is_skipped(self) -> None:
        """Lines with non-numeric coordinates must be silently skipped."""
        parse = self._import()
        bad = "ATOM      1  N   MET A   1      XX.XX  13.207  12.994\n"
        atoms = parse(bad)
        assert atoms == []


# ===========================================================================
# extract_interactions_distance tests (needs numpy — stdlib adjacent)
# ===========================================================================

class TestExtractInteractions:
    """
    Test distance-based interaction extraction with synthetic data.
    Requires: numpy (should always be available once docking env is set up).
    """

    @pytest.fixture(autouse=True)
    def _require_numpy(self):
        pytest.importorskip("numpy")

    def _import(self):
        sys.path.insert(0, str(REPO_ROOT))
        from docking_pipeline import extract_interactions_distance
        return extract_interactions_distance

    def _write_pdb(self, tmp_path: Path, content: str) -> Path:
        p = tmp_path / "receptor.pdb"
        p.write_text(content)
        return p

    def test_hbond_detected_within_cutoff(self, tmp_path: Path) -> None:
        """
        An O–O pair at 3.0 Å (< 3.5 Å cutoff) must be reported as an H-bond.
        """
        extract = self._import()

        # Ligand: O atom at (0, 0, 0)
        pose_pdbqt = (
            "ATOM      1  O   LIG A   1       0.000   0.000   0.000"
            "  1.00  0.00           O  \n"
        )
        # Receptor: O atom of MET793 at (0, 0, 3.0) — within 3.5 Å cutoff
        receptor_pdb = (
            "ATOM      1  O   MET A 793       0.000   0.000   3.000"
            "  1.00  0.00           O  \n"
        )
        rec_path = self._write_pdb(tmp_path, receptor_pdb)
        result = extract(pose_pdbqt, rec_path)
        assert "MET793_A" in result["hbond_residues"], (
            "O–O pair at 3.0 Å should be detected as H-bond"
        )

    def test_hbond_not_detected_beyond_cutoff(self, tmp_path: Path) -> None:
        """O–O pair at 4.0 Å (> 3.5 Å cutoff) must NOT be reported."""
        extract = self._import()

        pose_pdbqt = (
            "ATOM      1  O   LIG A   1       0.000   0.000   0.000"
            "  1.00  0.00           O  \n"
        )
        receptor_pdb = (
            "ATOM      1  O   MET A 793       0.000   0.000   4.000"
            "  1.00  0.00           O  \n"
        )
        rec_path = self._write_pdb(tmp_path, receptor_pdb)
        result = extract(pose_pdbqt, rec_path)
        assert "MET793_A" not in result["hbond_residues"]

    def test_hydrophobic_detected_within_cutoff(self, tmp_path: Path) -> None:
        """C–C pair at 4.0 Å (< 4.5 Å cutoff) must be a hydrophobic contact."""
        extract = self._import()

        pose_pdbqt = (
            "ATOM      1  C   LIG A   1       0.000   0.000   0.000"
            "  1.00  0.00           C  \n"
        )
        receptor_pdb = (
            "ATOM      1  CB  LEU A 844       0.000   0.000   4.000"
            "  1.00  0.00           C  \n"
        )
        rec_path = self._write_pdb(tmp_path, receptor_pdb)
        result = extract(pose_pdbqt, rec_path)
        assert "LEU844_A" in result["hydrophobic_residues"]

    def test_polar_atom_not_in_hydrophobic(self, tmp_path: Path) -> None:
        """N atoms must never appear in the hydrophobic contact list."""
        extract = self._import()

        pose_pdbqt = (
            "ATOM      1  N   LIG A   1       0.000   0.000   0.000"
            "  1.00  0.00           N  \n"
        )
        receptor_pdb = (
            "ATOM      1  N   MET A 793       0.000   0.000   1.000"
            "  1.00  0.00           N  \n"
        )
        rec_path = self._write_pdb(tmp_path, receptor_pdb)
        result = extract(pose_pdbqt, rec_path)
        assert result["hydrophobic_residues"] == []

    def test_empty_pose_returns_empty_lists(self, tmp_path: Path) -> None:
        """Empty ligand PDBQT → no interactions."""
        extract = self._import()
        rec_path = self._write_pdb(tmp_path, "ATOM      1  N   MET A 793  0.0 0.0 0.0\n")
        result = extract("", rec_path)
        assert result["hbond_residues"] == []
        assert result["hydrophobic_residues"] == []

    def test_multiple_residues_detected(self, tmp_path: Path) -> None:
        """Multiple residues within cutoff must all appear in the output."""
        extract = self._import()

        pose_pdbqt = (
            "ATOM      1  O   LIG A   1       0.000   0.000   0.000"
            "  1.00  0.00           O  \n"
        )
        receptor_pdb = (
            "ATOM      1  O   MET A 793       0.000   0.000   2.000  1.00 0.00  O  \n"
            "ATOM      2  N   THR A 790       0.000   0.000   3.000  1.00 0.00  N  \n"
        )
        rec_path = self._write_pdb(tmp_path, receptor_pdb)
        result = extract(pose_pdbqt, rec_path)
        assert "MET793_A" in result["hbond_residues"]
        assert "THR790_A" in result["hbond_residues"]


# ===========================================================================
# write_results tests (needs pandas)
# ===========================================================================

class TestWriteResults:
    """Test that write_results() produces correct CSV output."""

    @pytest.fixture(autouse=True)
    def _require_pandas(self):
        pytest.importorskip("pandas")

    def _import(self):
        sys.path.insert(0, str(REPO_ROOT))
        from docking_pipeline import write_results
        return write_results

    def _make_results(self) -> list[dict]:
        return [
            {
                "name": "erlotinib", "smiles": "C#Cc1cccc",
                "affinity": -9.5,
                "hbond_residues": ["MET793_A"],
                "hydrophobic_residues": ["LEU844_A"],
                "status": "success",
            },
            {
                "name": "aspirin", "smiles": "CC(=O)Oc1ccccc1",
                "affinity": -5.1,
                "hbond_residues": [],
                "hydrophobic_residues": [],
                "status": "success",
            },
            {
                "name": "badmol", "smiles": "INVALID",
                "affinity": None,
                "hbond_residues": [],
                "hydrophobic_residues": [],
                "status": "prep_failed",
            },
        ]

    def test_csv_created(self, tmp_path: Path) -> None:
        import pandas as pd
        write = self._import()
        out = tmp_path / "results.csv"
        write(self._make_results(), out)
        assert out.exists()

    def test_all_rows_present(self, tmp_path: Path) -> None:
        import pandas as pd
        write = self._import()
        out = tmp_path / "results.csv"
        write(self._make_results(), out)
        df = pd.read_csv(out)
        assert len(df) == 3

    def test_sorted_by_affinity(self, tmp_path: Path) -> None:
        """
        Successful rows must be sorted by affinity ascending
        (most negative = best binding = first row).
        """
        import pandas as pd
        write = self._import()
        out = tmp_path / "results.csv"
        write(self._make_results(), out)
        df = pd.read_csv(out)
        success = df[df["status"] == "success"]
        affinities = success["affinity_kcal_mol"].tolist()
        assert affinities == sorted(affinities), (
            f"Affinities not sorted ascending: {affinities}"
        )

    def test_failures_at_bottom(self, tmp_path: Path) -> None:
        """Failed molecules must appear after all successes."""
        import pandas as pd
        write = self._import()
        out = tmp_path / "results.csv"
        write(self._make_results(), out)
        df = pd.read_csv(out)
        # Last row should be the failure
        assert df.iloc[-1]["status"] == "prep_failed"
        assert df.iloc[-1]["name"] == "badmol"

    def test_hbond_residues_semicolon_separated(self, tmp_path: Path) -> None:
        """H-bond residues must be a single semicolon-separated string per row."""
        import pandas as pd
        write = self._import()
        # Add a molecule with two H-bond residues
        results = self._make_results()
        results[0]["hbond_residues"] = ["MET793_A", "THR790_A"]
        out = tmp_path / "results.csv"
        write(results, out)
        df = pd.read_csv(out)
        erlotinib_row = df[df["name"] == "erlotinib"].iloc[0]
        assert "MET793_A" in erlotinib_row["hbond_residues"]
        assert "THR790_A" in erlotinib_row["hbond_residues"]
        assert ";" in erlotinib_row["hbond_residues"]

    def test_required_columns_present(self, tmp_path: Path) -> None:
        """CSV must have all required columns."""
        import pandas as pd
        write = self._import()
        out = tmp_path / "results.csv"
        write(self._make_results(), out)
        df = pd.read_csv(out)
        required = {
            "name", "smiles", "affinity_kcal_mol",
            "hbond_residues", "hydrophobic_residues", "status",
        }
        assert required.issubset(set(df.columns)), (
            f"Missing columns: {required - set(df.columns)}"
        )


# ===========================================================================
# define_grid_box tests (needs biopython + numpy)
# ===========================================================================

class TestDefineGridBox:
    """Test grid box computation with synthetic PDB containing a known ligand."""

    @pytest.fixture(autouse=True)
    def _require_deps(self):
        pytest.importorskip("Bio")
        pytest.importorskip("numpy")

    def _make_pdb_with_ligand(self, tmp_path: Path, coords: list[tuple]) -> Path:
        """
        Write a minimal PDB file containing a HETATM ligand with given coordinates.
        """
        lines = []
        for i, (x, y, z) in enumerate(coords, start=1):
            lines.append(
                f"HETATM{i:5d}  C1  TST A{200:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  "
            )
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text("\n".join(lines) + "\n")
        return pdb_path

    def test_centroid_is_correct(self, tmp_path: Path) -> None:
        """
        Box centre must equal the geometric centroid of ligand heavy atoms.
        For atoms at (0,0,0) and (2,0,0), centroid is (1,0,0).
        """
        sys.path.insert(0, str(REPO_ROOT))
        from docking_pipeline import define_grid_box

        pdb_path = self._make_pdb_with_ligand(tmp_path, [(0, 0, 0), (2, 0, 0)])
        box = define_grid_box(pdb_path, ligand_resname="TST", padding=0.0)
        assert abs(box["center_x"] - 1.0) < 0.01
        assert abs(box["center_y"] - 0.0) < 0.01
        assert abs(box["center_z"] - 0.0) < 0.01

    def test_box_size_includes_padding(self, tmp_path: Path) -> None:
        """
        Box size = ligand extent + 2 × padding.
        Atoms at x=0 and x=2 → extent=2 → size_x = 2 + 2×5 = 12 with padding=5.
        """
        sys.path.insert(0, str(REPO_ROOT))
        from docking_pipeline import define_grid_box

        pdb_path = self._make_pdb_with_ligand(tmp_path, [(0, 0, 0), (2, 0, 0)])
        box = define_grid_box(pdb_path, ligand_resname="TST", padding=5.0)
        assert abs(box["size_x"] - 12.0) < 0.01

    def test_wrong_resname_raises(self, tmp_path: Path) -> None:
        """Wrong resname must raise ValueError with a helpful message."""
        sys.path.insert(0, str(REPO_ROOT))
        from docking_pipeline import define_grid_box

        pdb_path = self._make_pdb_with_ligand(tmp_path, [(0, 0, 0)])
        with pytest.raises(ValueError, match="No heavy atoms found"):
            define_grid_box(pdb_path, ligand_resname="XXX")

    def test_returns_six_keys(self, tmp_path: Path) -> None:
        """Box dict must have exactly the six expected keys."""
        sys.path.insert(0, str(REPO_ROOT))
        from docking_pipeline import define_grid_box

        pdb_path = self._make_pdb_with_ligand(tmp_path, [(1, 2, 3)])
        box = define_grid_box(pdb_path, ligand_resname="TST", padding=10.0)
        expected_keys = {"center_x", "center_y", "center_z", "size_x", "size_y", "size_z"}
        assert set(box.keys()) == expected_keys


# ===========================================================================
# smiles_to_pdbqt tests — SKIPPED if rdkit/meeko not installed
# ===========================================================================

class TestSMILEStoPDBQT:
    """Test ligand preparation. Skipped automatically if rdkit/meeko absent."""

    @pytest.fixture(autouse=True)
    def _require_deps(self):
        pytest.importorskip("rdkit")
        pytest.importorskip("meeko")

    def _convert(self, name: str, smiles: str) -> str | None:
        sys.path.insert(0, str(REPO_ROOT))
        from docking_pipeline import smiles_to_pdbqt
        return smiles_to_pdbqt(name, smiles)

    def test_aspirin_produces_pdbqt(self) -> None:
        """Aspirin (simple, well-parameterised) must produce a valid PDBQT string."""
        pdbqt = self._convert("aspirin", "CC(=O)Oc1ccccc1C(=O)O")
        assert pdbqt is not None, "smiles_to_pdbqt() returned None for aspirin"
        assert "ATOM" in pdbqt or "HETATM" in pdbqt, (
            "PDBQT output must contain ATOM or HETATM records"
        )

    def test_erlotinib_produces_pdbqt(self) -> None:
        """Erlotinib (the primary validation molecule) must produce a PDBQT string."""
        pdbqt = self._convert(
            "erlotinib",
            "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",
        )
        assert pdbqt is not None, "smiles_to_pdbqt() returned None for erlotinib"

    def test_invalid_smiles_returns_none(self) -> None:
        """
        An invalid SMILES must return None, not raise an exception.
        The pipeline marks this molecule as prep_failed and continues.
        """
        pdbqt = self._convert("bad_mol", "THISISNOTSMILES!!!")
        assert pdbqt is None, "Invalid SMILES must return None"

    def test_empty_smiles_returns_none(self) -> None:
        pdbqt = self._convert("empty", "")
        assert pdbqt is None

    def test_pdbqt_contains_root_line(self) -> None:
        """
        Valid Meeko PDBQT output always starts with a ROOT record for rigid molecules.
        This confirms the output is properly formatted for Vina.
        """
        pdbqt = self._convert("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C")
        if pdbqt is None:
            pytest.skip("Meeko returned None for caffeine — check meeko version")
        assert "ROOT" in pdbqt, (
            "Meeko PDBQT must contain ROOT record. "
            f"Got first 200 chars: {pdbqt[:200]}"
        )
