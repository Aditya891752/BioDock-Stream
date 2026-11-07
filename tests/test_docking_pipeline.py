"""
test_docking_pipeline.py — Phase 0 tests for docking_pipeline.py
=================================================================
Tests that can be run RIGHT NOW without any conda environment or Vina installed.
They test the CLI argument parser and logging setup only.

Phase 1 will add tests for the scientific stages (receptor fetch, clean, dock).

Run with:
    pytest tests/test_docking_pipeline.py -v
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Helper: project root path
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).parent.parent
PIPELINE_SCRIPT = REPO_ROOT / "docking_pipeline.py"


# ---------------------------------------------------------------------------
# CLI tests (no scientific deps needed)
# ---------------------------------------------------------------------------

class TestCLIHelp:
    """The --help flag must exit 0 and print key argument names."""

    def test_help_exits_zero(self) -> None:
        """
        Running --help must exit with code 0 (success).

        argparse exits with 0 on --help by design. If this test fails it means
        the script has a syntax error or an import that crashes at startup.
        """
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT), "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"--help exited with code {result.returncode}.\n"
            f"stdout: {result.stdout}\n"
            f"stderr: {result.stderr}"
        )

    def test_help_contains_pdb_id(self) -> None:
        """--help output must mention --pdb-id so the user knows it's required."""
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT), "--help"],
            capture_output=True,
            text=True,
        )
        assert "--pdb-id" in result.stdout

    def test_help_contains_library(self) -> None:
        """--help output must mention --library."""
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT), "--help"],
            capture_output=True,
            text=True,
        )
        assert "--library" in result.stdout

    def test_help_contains_output_dir(self) -> None:
        """--help output must mention --output-dir."""
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT), "--help"],
            capture_output=True,
            text=True,
        )
        assert "--output-dir" in result.stdout

    def test_missing_required_args_exits_nonzero(self) -> None:
        """
        Running with no arguments must exit with code 2 (argparse error).

        --pdb-id and --library are required; omitting them should fail fast
        with a clear error message, not a Python exception.
        """
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 2, (
            "Expected exit code 2 (missing required args) "
            f"but got {result.returncode}"
        )

    def test_missing_library_exits_nonzero(self) -> None:
        """Running with --pdb-id but no --library must also exit 2."""
        result = subprocess.run(
            [sys.executable, str(PIPELINE_SCRIPT), "--pdb-id", "1M17"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 2


class TestCLIArgumentParsing:
    """
    Verify that CLI arguments are correctly parsed and reach run_pipeline().

    These tests run the full script with valid arguments and check that the
    scaffold logs the parsed values correctly. They do NOT perform docking.
    """

    def test_scaffold_logs_pdb_id(self, tmp_path: Path) -> None:
        """
        The scaffold run_pipeline() must log the PDB ID passed on the CLI.

        This confirms argument parsing works end-to-end: CLI → argparse →
        Namespace → run_pipeline() → logger.
        """
        result = subprocess.run(
            [
                sys.executable, str(PIPELINE_SCRIPT),
                "--pdb-id", "1M17",
                "--library", str(REPO_ROOT / "data" / "molecules_sanity.csv"),
                "--output-dir", str(tmp_path / "results"),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert "1M17" in result.stdout

    def test_scaffold_logs_chain(self, tmp_path: Path) -> None:
        """The scaffold must log the --chain value."""
        result = subprocess.run(
            [
                sys.executable, str(PIPELINE_SCRIPT),
                "--pdb-id", "1M17",
                "--library", str(REPO_ROOT / "data" / "molecules_sanity.csv"),
                "--output-dir", str(tmp_path / "results"),
                "--chain", "B",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "B" in result.stdout

    def test_default_seed_is_42(self, tmp_path: Path) -> None:
        """
        The default random seed must be 42.

        The seed controls Vina's random search; the default must match what
        is documented in rules.md so results are reproducible.
        """
        result = subprocess.run(
            [
                sys.executable, str(PIPELINE_SCRIPT),
                "--pdb-id", "1M17",
                "--library", str(REPO_ROOT / "data" / "molecules_sanity.csv"),
                "--output-dir", str(tmp_path / "results"),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "42" in result.stdout  # seed=42 is logged

    def test_output_dir_is_created(self, tmp_path: Path) -> None:
        """
        The output directory must be created automatically if it does not exist.

        A researcher should not need to manually mkdir before running the pipeline.
        """
        output_dir = tmp_path / "new_results_dir"
        assert not output_dir.exists(), "Directory should not exist before the run"

        subprocess.run(
            [
                sys.executable, str(PIPELINE_SCRIPT),
                "--pdb-id", "1M17",
                "--library", str(REPO_ROOT / "data" / "molecules_sanity.csv"),
                "--output-dir", str(output_dir),
            ],
            capture_output=True,
            text=True,
        )
        assert output_dir.exists(), "run_pipeline() must create the output directory"

    def test_log_file_is_created(self, tmp_path: Path) -> None:
        """
        A run.log file must be created inside the output directory.

        This log file is the primary debugging tool when docking fails on a
        specific molecule. If it doesn't exist, diagnosing failures is hard.
        """
        output_dir = tmp_path / "results"
        subprocess.run(
            [
                sys.executable, str(PIPELINE_SCRIPT),
                "--pdb-id", "1M17",
                "--library", str(REPO_ROOT / "data" / "molecules_sanity.csv"),
                "--output-dir", str(output_dir),
            ],
            capture_output=True,
            text=True,
        )
        assert (output_dir / "run.log").exists(), "run.log must be written to output-dir"


class TestSanityCSV:
    """Verify the sanity test CSV has the right structure."""

    def test_sanity_csv_exists(self) -> None:
        """data/molecules_sanity.csv must exist in the repo."""
        assert (REPO_ROOT / "data" / "molecules_sanity.csv").exists()

    def test_sanity_csv_has_header(self) -> None:
        """The CSV must start with name,smiles,is_active."""
        content = (REPO_ROOT / "data" / "molecules_sanity.csv").read_text()
        first_line = content.strip().splitlines()[0]
        assert first_line == "name,smiles,is_active"

    def test_sanity_csv_has_five_molecules(self) -> None:
        """There must be exactly 5 molecules (header + 5 data rows = 6 lines)."""
        content = (REPO_ROOT / "data" / "molecules_sanity.csv").read_text()
        data_lines = [l for l in content.strip().splitlines() if l and not l.startswith("name")]
        assert len(data_lines) == 5, f"Expected 5 molecules, got {len(data_lines)}"

    def test_sanity_csv_contains_erlotinib(self) -> None:
        """erlotinib must be in the sanity CSV — it's the primary validation molecule."""
        content = (REPO_ROOT / "data" / "molecules_sanity.csv").read_text()
        assert "erlotinib" in content

    def test_sanity_csv_contains_aspirin(self) -> None:
        """aspirin must be in the sanity CSV — used as a negative control."""
        content = (REPO_ROOT / "data" / "molecules_sanity.csv").read_text()
        assert "aspirin" in content
