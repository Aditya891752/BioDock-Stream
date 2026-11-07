"""
test_validate_redock.py — Phase 2 test suite
=============================================
Tests for every function in validate_redock.py that has testable logic.

Test groups
-----------
  TestComputeRMSD          — the core maths. No external deps. Always runs.
  TestExtractCrystalCoords — PDB parsing with synthetic text. Needs: biopython.
  TestExtractDockedCoords  — PDBQT parsing.      No external deps. Always runs.
  TestWriteReport          — report file content. No external deps. Always runs.
  TestCLIHelp              — --help and arg parsing. No external deps.

Run with:
    pytest tests/test_validate_redock.py -v
    pytest tests/test_validate_redock.py -v -k "RMSD"   # just the maths tests
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Optional

import pytest
import numpy as np

# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------
REPO_ROOT   = Path(__file__).parent.parent
SCRIPT      = REPO_ROOT / "validate_redock.py"


def _import_fn(name: str):
    """Import a function from validate_redock.py by name."""
    sys.path.insert(0, str(REPO_ROOT))
    import validate_redock as vr
    return getattr(vr, name)


# ===========================================================================
# TestComputeRMSD — core RMSD mathematics
# (needs numpy only — always runs once docking env is active)
# ===========================================================================

class TestComputeRMSD:
    """
    Tests for compute_rmsd().

    These cover the most important function in Phase 2.
    If compute_rmsd() is wrong, everything downstream is wrong.
    No biopython, rdkit, or vina needed.
    """

    @pytest.fixture(autouse=True)
    def _require_numpy(self):
        pytest.importorskip("numpy")

    def _rmsd(self, c_coords, c_elems, d_coords, d_elems) -> Optional[float]:
        fn = _import_fn("compute_rmsd")
        return fn(
            np.array(c_coords, dtype=float),
            c_elems,
            np.array(d_coords, dtype=float),
            d_elems,
        )

    # ── Basic correctness ─────────────────────────────────────────────────

    def test_identical_poses_rmsd_zero(self) -> None:
        """
        When crystal and docked poses are identical, RMSD must be 0.
        This is the simplest sanity check for the formula.
        """
        coords = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
        elems  = ["C", "N", "O"]
        rmsd = self._rmsd(coords, elems, coords, elems)
        assert rmsd is not None
        assert rmsd < 1e-6, f"Expected ~0, got {rmsd}"

    def test_known_rmsd_single_atom(self) -> None:
        """
        Single C atom displaced 3.0 Å along z.
        RMSD = sqrt((0² + 0² + 3²) / 1) = 3.0 Å.
        """
        crystal = [[0.0, 0.0, 0.0]]
        docked  = [[0.0, 0.0, 3.0]]
        rmsd = self._rmsd(crystal, ["C"], docked, ["C"])
        assert abs(rmsd - 3.0) < 1e-6, f"Expected 3.0, got {rmsd}"

    def test_known_rmsd_two_atoms(self) -> None:
        """
        Two C atoms: first displaced 0 Å, second displaced 2 Å in y.
        RMSD = sqrt((0² + 2²) / 2) = sqrt(2) ≈ 1.4142 Å.
        """
        crystal = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
        docked  = [[0.0, 0.0, 0.0], [3.0, 2.0, 0.0]]
        rmsd = self._rmsd(crystal, ["C", "C"], docked, ["C", "C"])
        expected = float(np.sqrt(2.0))
        assert abs(rmsd - expected) < 1e-5, f"Expected {expected:.4f}, got {rmsd}"

    def test_known_rmsd_mixed_elements(self) -> None:
        """
        One C (displaced 0 Å) and one N (displaced 2 Å).
        RMSD = sqrt((0 + 4) / 2) = sqrt(2) ≈ 1.4142 Å.
        C and N are matched separately within their element groups.
        """
        crystal = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]]
        docked  = [[0.0, 0.0, 0.0], [5.0, 2.0, 0.0]]
        rmsd = self._rmsd(crystal, ["C", "N"], docked, ["C", "N"])
        expected = float(np.sqrt(2.0))
        assert abs(rmsd - expected) < 1e-5, f"Expected {expected:.4f}, got {rmsd}"

    def test_rmsd_below_2a_threshold(self) -> None:
        """
        A good redocking result should give RMSD < 2 Å.
        Test with a displacement of 1.5 Å per atom.
        """
        crystal = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
        docked  = [[1.5, 0.0, 0.0], [4.5, 0.0, 0.0]]  # shifted 1.5 Å in x
        rmsd = self._rmsd(crystal, ["C", "C"], docked, ["C", "C"])
        assert rmsd < 2.0, (
            f"Expected RMSD < 2.0 Å for 1.5 Å displacement, got {rmsd:.3f}"
        )

    def test_rmsd_above_2a_threshold(self) -> None:
        """
        A bad redocking result gives RMSD > 2 Å.
        Test with a displacement of 3 Å per atom.
        """
        crystal = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
        docked  = [[0.0, 0.0, 3.0], [3.0, 0.0, 3.0]]  # shifted 3 Å in z
        rmsd = self._rmsd(crystal, ["C", "C"], docked, ["C", "C"])
        assert rmsd > 2.0, (
            f"Expected RMSD > 2.0 Å for 3.0 Å displacement, got {rmsd:.3f}"
        )

    # ── Element matching ──────────────────────────────────────────────────

    def test_element_matching_does_not_cross_elements(self) -> None:
        """
        A C atom must never be matched to an N atom.
        Crystal: C at (0,0,0), N at (10,0,0)
        Docked:  C at (10,0,0), N at (0,0,0)
        If elements are matched correctly: C→C (10 Å) and N→N (10 Å), RMSD=10.
        If elements are crossed: C→N (0 Å) and N→C (0 Å), RMSD=0 (WRONG).
        """
        crystal = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
        # In docked, atom 0 is N at (10,0,0) and atom 1 is C at (0,0,0) -> wait:
        # If crystal C is at (0,0,0) and docked C is at (10,0,0), then C->C distance is 10.
        # If crystal N is at (10,0,0) and docked N is at (0,0,0), then N->N distance is 10.
        # In docked: coords=[[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]], elements=["N", "C"]
        # Then docked N is at (0,0,0) and docked C is at (10,0,0).
        docked  = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
        # Scenario A: same elements (C, C) — greedy matching may flip
        # Scenario B: different elements (C, N) — must NOT cross
        rmsd_mixed = self._rmsd(crystal, ["C", "N"], docked, ["N", "C"])
        # C at (0,0,0) matches docked C at (10,0,0) => dist = 10
        # N at (10,0,0) matches docked N at (0,0,0) => dist = 10
        # RMSD = sqrt((100 + 100) / 2) = 10 Å
        assert rmsd_mixed > 5.0, (
            "Elements must not be matched across types. "
            f"Expected RMSD ~10.0 Å, got {rmsd_mixed:.3f}"
        )

    def test_element_mismatch_returns_partial_or_none(self) -> None:
        """
        When crystal has 2 C atoms but docked has only 1 C, the function
        must NOT crash. It may return a partial result or None.
        """
        fn = _import_fn("compute_rmsd")
        result = fn(
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            ["C", "C"],
            np.array([[0.0, 0.0, 0.0]]),
            ["C"],
        )
        # Must not raise; return value is None or a float
        assert result is None or isinstance(result, float)

    def test_empty_coords_returns_none(self) -> None:
        """Empty coordinate arrays must return None, not crash."""
        fn = _import_fn("compute_rmsd")
        result = fn(
            np.array([]).reshape(0, 3), [],
            np.array([]).reshape(0, 3), [],
        )
        assert result is None

    # ── RMSD is symmetric ────────────────────────────────────────────────

    def test_rmsd_is_symmetric(self) -> None:
        """
        RMSD(A, B) == RMSD(B, A).
        Swapping crystal and docked must give the same result.
        """
        crystal = [[0.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 3.0, 1.0]]
        docked  = [[0.5, 0.0, 0.0], [2.5, 1.0, 0.5], [0.0, 3.5, 1.0]]
        elems   = ["C", "C", "N"]

        rmsd_ab = self._rmsd(crystal, elems, docked, elems)
        rmsd_ba = self._rmsd(docked, elems, crystal, elems)

        assert rmsd_ab is not None and rmsd_ba is not None
        assert abs(rmsd_ab - rmsd_ba) < 1e-5, (
            f"RMSD should be symmetric: {rmsd_ab:.5f} vs {rmsd_ba:.5f}"
        )

    # ── Realistic AQ4 scale ───────────────────────────────────────────────

    def test_realistic_atom_count(self) -> None:
        """
        AQ4 (erlotinib) has 32 heavy atoms (23 C, 4 N, 5 O).
        A correctly working RMSD function must handle this scale.
        """
        rng = np.random.default_rng(seed=0)
        # Simulate a 'good' docking: small random displacements (< 1 Å)
        crystal = rng.uniform(0, 20, size=(32, 3))
        docked  = crystal + rng.uniform(-0.5, 0.5, size=(32, 3))
        elements = (["C"] * 23) + (["N"] * 4) + (["O"] * 5)

        fn = _import_fn("compute_rmsd")
        rmsd = fn(
            crystal, elements,
            docked, elements,
        )
        assert rmsd is not None
        # With uniform ±0.5 Å noise, RMSD should be small (< 1 Å)
        assert rmsd < 1.5, (
            f"Expected small RMSD with ±0.5 Å noise, got {rmsd:.3f} Å"
        )


# ===========================================================================
# TestExtractDockedCoords — PDBQT string parsing (no scientific deps)
# ===========================================================================

class TestExtractDockedCoords:
    """Test extract_docked_coords() with synthetic PDBQT strings."""

    # Minimal synthetic PDBQT output (Vina format)
    # Only the ATOM/HETATM lines matter; ROOT/BRANCH lines are ignored by our parser
    SAMPLE_PDBQT = textwrap.dedent("""\
        MODEL 1
        REMARK VINA RESULT:     -9.5      0.000      0.000
        ROOT
        ATOM      1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00    0.000 C
        ATOM      2  N1  LIG A   1       4.000   5.000   6.000  1.00  0.00    0.000 N
        ATOM      3  O1  LIG A   1       7.000   8.000   9.000  1.00  0.00    0.000 O
        ATOM      4  H1  LIG A   1       0.500   0.500   0.500  1.00  0.00    0.000 H
        ENDROOT
        TORSDOF 0
        ENDMDL
    """)

    def test_heavy_atoms_only(self) -> None:
        """H atoms must be excluded from the parsed output."""
        fn = _import_fn("extract_docked_coords")
        coords, elems = fn(self.SAMPLE_PDBQT)
        assert "H" not in elems, f"Hydrogen found in elements: {elems}"

    def test_correct_atom_count(self) -> None:
        """3 heavy atoms (C, N, O) must be parsed; H is skipped."""
        fn = _import_fn("extract_docked_coords")
        coords, elems = fn(self.SAMPLE_PDBQT)
        assert len(coords) == 3
        assert len(elems) == 3

    def test_coordinates_parsed_correctly(self) -> None:
        """Coordinates of the first atom must match the PDBQT line."""
        fn = _import_fn("extract_docked_coords")
        coords, elems = fn(self.SAMPLE_PDBQT)
        assert abs(coords[0][0] - 1.0) < 1e-5   # x
        assert abs(coords[0][1] - 2.0) < 1e-5   # y
        assert abs(coords[0][2] - 3.0) < 1e-5   # z

    def test_elements_in_order(self) -> None:
        """Elements must be returned in ATOM-record order."""
        fn = _import_fn("extract_docked_coords")
        _, elems = fn(self.SAMPLE_PDBQT)
        assert elems == ["C", "N", "O"]

    def test_empty_string_raises(self) -> None:
        """Empty PDBQT must raise ValueError with an informative message."""
        fn = _import_fn("extract_docked_coords")
        with pytest.raises(ValueError, match="No heavy atoms"):
            fn("")

    def test_only_hydrogen_raises(self) -> None:
        """PDBQT with only H atoms must raise ValueError."""
        fn = _import_fn("extract_docked_coords")
        h_only = (
            "ATOM      1  H1  LIG A   1       0.0   0.0   0.0  1.00  0.00    0.000 H\n"
        )
        with pytest.raises(ValueError, match="No heavy atoms"):
            fn(h_only)


# ===========================================================================
# TestExtractCrystalCoords — PDB parsing (needs biopython)
# ===========================================================================

class TestExtractCrystalCoords:
    """
    Test extract_crystal_coords() with a synthetic PDB containing a HETATM ligand.
    Skipped if biopython is not installed.
    """

    @pytest.fixture(autouse=True)
    def _require_biopython(self):
        pytest.importorskip("Bio")

    def _write_pdb(self, tmp_path: Path, hetatm_lines: list[str]) -> Path:
        """Write a minimal PDB file with given HETATM lines."""
        content = "\n".join(hetatm_lines) + "\nEND\n"
        p = tmp_path / "test.pdb"
        p.write_text(content)
        return p

    def _make_hetatm(
        self, serial: int, atom_name: str, resname: str,
        x: float, y: float, z: float, element: str,
    ) -> str:
        """Format a HETATM line in standard PDB fixed-width format."""
        return (
            f"HETATM{serial:5d} {atom_name:<4s} {resname:3s} A{200:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  "
        )

    def test_extracts_heavy_atoms(self, tmp_path: Path) -> None:
        """Crystal coord extraction must return heavy atoms, not H."""
        fn = _import_fn("extract_crystal_coords")
        pdb = self._write_pdb(tmp_path, [
            self._make_hetatm(1, "C1 ", "AQ4", 1.0, 2.0, 3.0, "C"),
            self._make_hetatm(2, "N1 ", "AQ4", 4.0, 5.0, 6.0, "N"),
            self._make_hetatm(3, "H1 ", "AQ4", 0.1, 0.1, 0.1, "H"),  # should be excluded
        ])
        coords, elems = fn(pdb, "AQ4")
        assert len(coords) == 2, f"Expected 2 heavy atoms, got {len(coords)}"
        assert "H" not in elems

    def test_coordinates_correct(self, tmp_path: Path) -> None:
        """Coordinates must match the PDB HETATM record exactly."""
        fn = _import_fn("extract_crystal_coords")
        pdb = self._write_pdb(tmp_path, [
            self._make_hetatm(1, "C1 ", "TST", 7.5, -3.2, 11.0, "C"),
        ])
        coords, elems = fn(pdb, "TST")
        assert abs(coords[0][0] - 7.5) < 0.01
        assert abs(coords[0][1] - (-3.2)) < 0.01
        assert abs(coords[0][2] - 11.0) < 0.01

    def test_wrong_resname_raises(self, tmp_path: Path) -> None:
        """Wrong resname must raise ValueError with a helpful message."""
        fn = _import_fn("extract_crystal_coords")
        pdb = self._write_pdb(tmp_path, [
            self._make_hetatm(1, "C1 ", "AQ4", 1.0, 2.0, 3.0, "C"),
        ])
        with pytest.raises(ValueError, match="No heavy atoms found"):
            fn(pdb, "XXX")

    def test_multiple_atoms_returned(self, tmp_path: Path) -> None:
        """All heavy atoms of the ligand must be returned."""
        fn = _import_fn("extract_crystal_coords")
        pdb = self._write_pdb(tmp_path, [
            self._make_hetatm(i, f"C{i} ", "AQ4", float(i), 0.0, 0.0, "C")
            for i in range(1, 6)
        ])
        coords, elems = fn(pdb, "AQ4")
        assert len(coords) == 5

    def test_elements_list_matches_coords(self, tmp_path: Path) -> None:
        """len(elements) must equal len(coords)."""
        fn = _import_fn("extract_crystal_coords")
        pdb = self._write_pdb(tmp_path, [
            self._make_hetatm(1, "C1 ", "AQ4", 1.0, 0.0, 0.0, "C"),
            self._make_hetatm(2, "N1 ", "AQ4", 2.0, 0.0, 0.0, "N"),
            self._make_hetatm(3, "O1 ", "AQ4", 3.0, 0.0, 0.0, "O"),
        ])
        coords, elems = fn(pdb, "AQ4")
        assert len(coords) == len(elems)


# ===========================================================================
# TestWriteReport — report file content (no scientific deps)
# ===========================================================================

class TestWriteReport:
    """Test that write_report() creates a file with expected content."""

    def _write(
        self, tmp_path: Path,
        rmsd: Optional[float] = 1.5,
        score: Optional[float] = -9.5,
    ) -> Path:
        fn = _import_fn("write_report")
        report_path = tmp_path / "redock_report.txt"
        fn(
            report_path,
            pdb_id="1M17",
            ligand_resname="AQ4",
            n_crystal_atoms=32,
            n_docked_atoms=32,
            best_vina_score=score,
            rmsd=rmsd,
            box={"center_x": 1.0, "center_y": 2.0, "center_z": 3.0,
                 "size_x": 20.0, "size_y": 20.0, "size_z": 20.0},
        )
        return report_path

    def test_report_file_created(self, tmp_path: Path) -> None:
        """Report file must be created."""
        path = self._write(tmp_path)
        assert path.exists()

    def test_report_contains_pdb_id(self, tmp_path: Path) -> None:
        path = self._write(tmp_path)
        assert "1M17" in path.read_text()

    def test_report_contains_rmsd_value(self, tmp_path: Path) -> None:
        path = self._write(tmp_path, rmsd=1.234)
        assert "1.234" in path.read_text()

    def test_report_contains_vina_score(self, tmp_path: Path) -> None:
        path = self._write(tmp_path, score=-9.5)
        assert "-9.50" in path.read_text()

    def test_pass_verdict_when_rmsd_below_threshold(self, tmp_path: Path) -> None:
        """RMSD 1.5 Å < 2.0 Å → verdict must be PASS."""
        path = self._write(tmp_path, rmsd=1.5)
        assert "PASS" in path.read_text()

    def test_fail_verdict_when_rmsd_above_threshold(self, tmp_path: Path) -> None:
        """RMSD 2.5 Å > 2.0 Å → verdict must be FAIL."""
        path = self._write(tmp_path, rmsd=2.5)
        assert "FAIL" in path.read_text()

    def test_fail_verdict_when_rmsd_is_none(self, tmp_path: Path) -> None:
        """When RMSD is None (docking failed), report must still say FAIL."""
        path = self._write(tmp_path, rmsd=None)
        assert "FAIL" in path.read_text()

    def test_report_contains_threshold(self, tmp_path: Path) -> None:
        """Report must state the 2.0 Å threshold."""
        path = self._write(tmp_path)
        text = path.read_text()
        assert "2.0" in text

    def test_report_contains_atom_counts(self, tmp_path: Path) -> None:
        """Crystal and docked atom counts must appear in the report."""
        path = self._write(tmp_path)
        text = path.read_text()
        assert "32" in text   # n_crystal_atoms and n_docked_atoms

    def test_report_contains_box_centre(self, tmp_path: Path) -> None:
        """Grid box centre coordinates must appear in the report."""
        path = self._write(tmp_path)
        text = path.read_text()
        assert "1.00" in text   # center_x
        assert "2.00" in text   # center_y


# ===========================================================================
# TestCLIHelp — validate_redock.py CLI (no scientific deps)
# ===========================================================================

class TestCLIHelp:
    """The --help flag must work and show key arguments."""

    def _run_help(self) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, str(SCRIPT), "--help"],
            capture_output=True, text=True,
        )

    def test_help_exits_zero(self) -> None:
        """--help must exit 0."""
        assert self._run_help().returncode == 0

    def test_help_mentions_pdb_id(self) -> None:
        assert "--pdb-id" in self._run_help().stdout

    def test_help_mentions_exhaustiveness(self) -> None:
        assert "--exhaustiveness" in self._run_help().stdout

    def test_help_mentions_output_dir(self) -> None:
        assert "--output-dir" in self._run_help().stdout

    def test_help_mentions_rmsd_threshold(self) -> None:
        """Help text must mention the 2 Å threshold so users know the target."""
        assert "2" in self._run_help().stdout

    def test_default_pdb_id_is_1m17(self) -> None:
        """Default PDB ID must be 1M17 (the EGFR + erlotinib structure)."""
        assert "1M17" in self._run_help().stdout

    def test_default_ligand_resname_is_aq4(self) -> None:
        """Default ligand resname must be AQ4."""
        assert "AQ4" in self._run_help().stdout
