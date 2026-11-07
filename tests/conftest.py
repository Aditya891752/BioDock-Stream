"""
conftest.py — pytest configuration and shared fixtures
=======================================================
This file is loaded automatically by pytest before any tests run.
Add shared fixtures and configuration here so all test modules can use them
without importing from each other.

Currently empty in Phase 0 — fixtures will be added in Phase 1 when
test coverage for the pipeline stages begins.
"""

from pathlib import Path
import pytest


# ---------------------------------------------------------------------------
# Path fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def data_dir() -> Path:
    """
    Return the path to the project data/ directory.

    Use this fixture in tests that need to load the sanity CSV or library CSV
    without hardcoding absolute paths.
    """
    return Path(__file__).parent.parent / "data"


@pytest.fixture
def results_dir(tmp_path: Path) -> Path:
    """
    Return a temporary directory for test outputs.

    pytest's built-in tmp_path fixture creates a unique directory per test
    invocation that is automatically cleaned up after the test session.
    Use this for any test that writes files (logs, CSVs, PDBQT strings).
    """
    return tmp_path / "results"
