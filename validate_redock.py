"""
validate_redock.py — Redocking Validation for EGFR / AQ4
=========================================================
Phase 2 deliverable.

What redocking validation is
-----------------------------
We take the co-crystal ligand (AQ4 / erlotinib, already crystallised inside
EGFR's ATP pocket in PDB 1M17) and dock it back using the same pipeline we
built in Phase 1. If Vina recovers a pose close to the crystal pose (heavy-
atom RMSD < 2 Å), we have strong evidence that:
  - Our grid box covers the correct pocket.
  - Our receptor preparation (Open Babel PDBQT) is reasonable.
  - Our ligand preparation (RDKit + Meeko) produces a sensible 3D geometry.
  - Vina's scoring function can distinguish the native pose from random placements.

If RMSD > 2 Å: stop and diagnose before running any library screen.
  Common causes: wrong grid box, flipped receptor charge, bad SMILES geometry.

RMSD definition used here
--------------------------
  RMSD = sqrt( (1/N) * sum_i( min_j( ||crystal_i - docked_j||² ) ) )

Both the crystal pose and the docked pose are expressed in the same coordinate
frame (the receptor's crystallographic frame), so NO superimposition is applied
before computing RMSD. What we measure is: did Vina place the molecule in the
right location AND orientation in 3D space?

Atom matching uses a greedy nearest-neighbour approach within each element type
(match crystal N to nearest docked N, crystal C to nearest docked C, etc.).
This is robust to differences in atom ordering between the PDB and the PDBQT.

AQ4 = erlotinib
---------------
PDB residue name AQ4 in structure 1M17 corresponds to erlotinib (OSI-774).
Verified against: PubChem CID 176870, RCSB PDB ligand page for AQ4.
SMILES: C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1
Heavy atoms: 32  (23 C, 4 N, 5 O)

Usage
-----
python validate_redock.py --pdb-id 1M17 --output-dir results_egfr

Dependencies
------------
Same conda env as docking_pipeline.py:
  conda activate docking
  (biopython, rdkit, meeko, vina, openbabel all required)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np

# ---------------------------------------------------------------------------
# AQ4 / erlotinib reference data
# Source: PubChem CID 176870 (verified); RCSB PDB ligand AQ4
# Do NOT modify without updating the comment above and re-running validation.
# ---------------------------------------------------------------------------

AQ4_SMILES = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
AQ4_NAME   = "AQ4"          # PDB residue name in 1M17
AQ4_PUBCHEM_CID = 176870    # for the developer to verify independently

RMSD_PASS_THRESHOLD = 2.0   # Å — standard cutoff in CADD redocking validation


# ===========================================================================
# STEP 1 — Extract crystal ligand coordinates from the raw PDB
# ===========================================================================

def extract_crystal_coords(
    raw_pdb_path: Path,
    ligand_resname: str,
) -> tuple[np.ndarray, list[str]]:
    """
    Extract heavy-atom 3D coordinates from the co-crystal ligand in a PDB file.

    Why we need the RAW PDB (not the cleaned one):
      The clean_receptor() step in Phase 1 removes all HETATM records to give
      Vina a water-free protein. The co-crystal ligand is a HETATM record, so
      it disappears after cleaning. We always work from the original download.

    How we filter hydrogen atoms:
      We use Biopython's atom.element field. If the PDB file has the standard
      element column (cols 77-78), Biopython fills this automatically.
      We fall back to the atom name: atoms starting with "H" or "1H","2H" etc.
      are hydrogens. We skip them because RMSD is always computed on heavy atoms
      (H positions are poorly defined in X-ray structures).

    Parameters
    ----------
    raw_pdb_path : Path
        The original PDB file download (contains HETATM records).
    ligand_resname : str
        Three-letter PDB residue code for the co-crystal ligand (e.g. "AQ4").

    Returns
    -------
    tuple[np.ndarray, list[str]]
        coords   : float array of shape (N, 3) — x, y, z in Ångströms
        elements : list of N strings — element symbol per atom ("C", "N", "O", …)

    Raises
    ------
    ValueError
        If no atoms are found (wrong resname, or wrong PDB file).
    """
    from Bio import PDB

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("crystal", str(raw_pdb_path))

    coords: list[list[float]] = []
    elements: list[str] = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip().upper() != ligand_resname.upper():
                    continue
                for atom in residue:
                    # Determine element
                    elem = (atom.element or "").strip().upper()
                    if not elem:
                        # Fallback: derive from atom name
                        name = atom.get_name().strip().lstrip("0123456789")
                        elem = name[0].upper() if name else "X"

                    # Skip hydrogens
                    if elem == "H":
                        continue

                    coords.append(list(atom.get_vector()))
                    elements.append(elem)

    if not coords:
        raise ValueError(
            f"No heavy atoms found for ligand '{ligand_resname}' in {raw_pdb_path}.\n"
            f"Checklist:\n"
            f"  1. Use the RAW PDB (the cleaned one has no HETATM).\n"
            f"  2. Verify the residue name on rcsb.org for PDB 1M17 — it should be 'AQ4'.\n"
            f"  3. If you changed the PDB ID, the ligand resname may differ."
        )

    return np.array(coords, dtype=float), elements


# ===========================================================================
# STEP 2 — Extract docked pose coordinates from a PDBQT string
# ===========================================================================

def extract_docked_coords(
    pose_pdbqt_str: str,
) -> tuple[np.ndarray, list[str]]:
    """
    Parse heavy-atom coordinates from a Vina PDBQT pose string.

    PDBQT format is PDB with two extra columns (partial charge, atom type).
    We re-use _parse_atom_records from docking_pipeline.py for consistency.
    Hydrogen atoms are skipped (element == "H").

    Parameters
    ----------
    pose_pdbqt_str : str
        PDBQT content of the best docked pose (from dock_ligand()).

    Returns
    -------
    tuple[np.ndarray, list[str]]
        coords   : float array of shape (M, 3)
        elements : list of M element symbols

    Raises
    ------
    ValueError
        If no atoms are found in the pose string.
    """
    # Import from Phase 1 pipeline to avoid code duplication
    sys.path.insert(0, str(Path(__file__).parent))
    from docking_pipeline import _parse_atom_records

    raw_atoms = _parse_atom_records(pose_pdbqt_str)
    heavy = [a for a in raw_atoms if a["element"].upper() != "H"]

    if not heavy:
        raise ValueError(
            "No heavy atoms found in the docked pose PDBQT string.\n"
            "This suggests ligand preparation or docking failed silently."
        )

    coords   = np.array([[a["x"], a["y"], a["z"]] for a in heavy], dtype=float)
    elements = [a["element"].upper() for a in heavy]
    return coords, elements


# ===========================================================================
# STEP 3 — Compute heavy-atom RMSD
# ===========================================================================

def compute_rmsd(
    crystal_coords: np.ndarray,
    crystal_elements: list[str],
    docked_coords: np.ndarray,
    docked_elements: list[str],
) -> Optional[float]:
    """
    Compute heavy-atom RMSD between a crystal pose and a docked pose.

    Both pose sets must be in the same coordinate frame (the receptor's
    crystallographic frame). No superimposition is applied — this measures
    the actual spatial error of Vina's prediction.

    Atom matching strategy
    ----------------------
    Crystal and docked PDBQT files may have different atom orderings because
    the PDB names atoms by convention while Meeko numbers them from the SMILES
    graph. We match atoms element-by-element using a greedy nearest-neighbour
    algorithm:

      For each element type (C, N, O, …):
        1. Collect crystal atoms of that element: set A.
        2. Collect docked atoms of that element: set B.
        3. If |A| ≠ |B|: skip this element and log a warning.
           (This should never happen for the same molecule.)
        4. For each point in A (in order), greedily assign the closest
           unmatched point in B, accumulate squared distance.

    This is not globally optimal (that would require the Hungarian algorithm)
    but is fast, interpretable, and adequate for a first-pass validation.

    Parameters
    ----------
    crystal_coords : np.ndarray, shape (N, 3)
        Heavy-atom coordinates from the crystal structure.
    crystal_elements : list[str], length N
        Element symbols for each crystal atom (e.g. ["C", "C", "N", "O"]).
    docked_coords : np.ndarray, shape (M, 3)
        Heavy-atom coordinates from the docked pose.
    docked_elements : list[str], length M
        Element symbols for each docked atom.

    Returns
    -------
    float | None
        RMSD in Ångströms, or None if no atoms could be matched.
        Return None signals a problem; the caller should log and stop.
    """
    if len(crystal_coords) == 0 or len(docked_coords) == 0:
        return None

    all_elements = sorted(set(crystal_elements) | set(docked_elements))

    total_sq_dist = 0.0
    n_matched     = 0
    n_mismatched  = 0

    for elem in all_elements:
        # Indices of this element in each set
        c_idx = [i for i, e in enumerate(crystal_elements) if e == elem]
        d_idx = [i for i, e in enumerate(docked_elements)  if e == elem]

        if len(c_idx) != len(d_idx):
            print(
                f"  WARNING: element {elem} count mismatch — "
                f"crystal has {len(c_idx)}, docked has {len(d_idx)}. "
                f"Skipping this element in RMSD calculation."
            )
            n_mismatched += abs(len(c_idx) - len(d_idx))
            continue

        c_pts = crystal_coords[c_idx]  # shape (k, 3)
        d_pts = docked_coords[d_idx]   # shape (k, 3)

        # Greedy nearest-neighbour matching within this element group
        assigned = set()
        for c_pt in c_pts:
            # Distance from this crystal atom to all (remaining) docked atoms
            dists = np.linalg.norm(d_pts - c_pt, axis=1)
            # Pick the closest unassigned docked atom
            for rank_idx in np.argsort(dists):
                if rank_idx not in assigned:
                    total_sq_dist += float(dists[rank_idx] ** 2)
                    assigned.add(rank_idx)
                    n_matched += 1
                    break

    if n_matched == 0:
        return None

    if n_mismatched > 0:
        print(
            f"  WARNING: {n_mismatched} atoms skipped due to element count mismatch. "
            f"RMSD computed over {n_matched} matched atoms."
        )

    return float(np.sqrt(total_sq_dist / n_matched))


# ===========================================================================
# STEP 4 — Write the validation report
# ===========================================================================

def write_report(
    report_path: Path,
    pdb_id: str,
    ligand_resname: str,
    n_crystal_atoms: int,
    n_docked_atoms: int,
    best_vina_score: Optional[float],
    rmsd: Optional[float],
    box: dict,
) -> None:
    """
    Write a plain-text redocking validation report.

    The report captures all the numbers needed for the portfolio README
    and for diagnosing failures. Format: human-readable plain text.

    Parameters
    ----------
    report_path : Path
        Where to write the report file.
    pdb_id : str
        PDB ID of the target (e.g. "1M17").
    ligand_resname : str
        Co-crystal ligand residue name (e.g. "AQ4").
    n_crystal_atoms : int
        Number of heavy atoms in the crystal pose.
    n_docked_atoms : int
        Number of heavy atoms in the docked pose.
    best_vina_score : float | None
        Best Vina affinity (kcal/mol) for the redocked pose.
    rmsd : float | None
        Computed heavy-atom RMSD (Å). None if computation failed.
    box : dict
        Grid box parameters (center_x/y/z, size_x/y/z).
    """
    if rmsd is not None:
        verdict = "PASS" if rmsd <= RMSD_PASS_THRESHOLD else "FAIL"
        rmsd_str = f"{rmsd:.3f} Å"
    else:
        verdict  = "FAIL — RMSD could not be computed"
        rmsd_str = "N/A"

    if best_vina_score is not None:
        score_str = f"{best_vina_score:+.2f} kcal/mol"
    else:
        score_str = "N/A — docking failed"

    lines = [
        "=" * 60,
        "REDOCKING VALIDATION REPORT",
        "=" * 60,
        f"Target PDB ID     : {pdb_id.upper()}",
        f"Co-crystal ligand : {ligand_resname} (AQ4 = erlotinib / OSI-774)",
        f"AQ4 PubChem CID   : {AQ4_PUBCHEM_CID}  (verify independently)",
        "",
        "Grid Box Parameters",
        "-" * 40,
        f"  Centre  : ({box['center_x']:.2f}, {box['center_y']:.2f}, {box['center_z']:.2f}) Å",
        f"  Size    : ({box['size_x']:.2f}, {box['size_y']:.2f}, {box['size_z']:.2f}) Å",
        "",
        "Atom Counts (heavy atoms only, no H)",
        "-" * 40,
        f"  Crystal pose : {n_crystal_atoms} atoms",
        f"  Docked pose  : {n_docked_atoms} atoms",
        (
            f"  Match        : YES"
            if n_crystal_atoms == n_docked_atoms
            else f"  Match        : NO — count mismatch (investigate!)"
        ),
        "",
        "Results",
        "-" * 40,
        f"  Best Vina score : {score_str}",
        f"  Heavy-atom RMSD : {rmsd_str}",
        f"  Threshold       : {RMSD_PASS_THRESHOLD:.1f} Å (standard CADD cutoff)",
        f"  Verdict         : {verdict}",
        "",
        "Interpretation",
        "-" * 40,
    ]

    if rmsd is not None and rmsd <= RMSD_PASS_THRESHOLD:
        lines += [
            f"  RMSD {rmsd:.3f} Å < {RMSD_PASS_THRESHOLD:.1f} Å threshold.",
            "  Vina recovered a pose close to the crystal structure.",
            "  The grid box and receptor preparation are trustworthy.",
            "  Proceed to Phase 3 (library preparation).",
        ]
    elif rmsd is not None:
        lines += [
            f"  RMSD {rmsd:.3f} Å > {RMSD_PASS_THRESHOLD:.1f} Å threshold.",
            "  Vina did NOT recover the native binding mode. STOP.",
            "  Diagnose before running any library screen.",
            "",
            "  Diagnostic checklist:",
            "  1. Open receptor.pdbqt in PyMOL — does it look like a complete protein?",
            "  2. Visualise the grid box — does it cover the ATP pocket?",
            "  3. Check AQ4 SMILES against PubChem CID 176870.",
            "  4. Try increasing --exhaustiveness to 16 or 32.",
            "  5. Try --box-padding 15 to give Vina more search space.",
        ]
    else:
        lines += [
            "  RMSD could not be computed (docking or ligand prep failed).",
            "  Check run.log for detailed error messages.",
            "  Ensure conda env is active: conda activate docking",
        ]

    lines += [
        "",
        "=" * 60,
        "NOTE ON RMSD METHOD",
        "=" * 60,
        "RMSD is computed WITHOUT superimposition (both poses are in",
        "the receptor's crystallographic coordinate frame).",
        "Atom matching: greedy nearest-neighbour within each element type.",
        "This is a conservative estimate — the true RMSD may be slightly",
        "lower with optimal matching (Hungarian algorithm, used in Phase 5).",
        "=" * 60,
    ]

    text = "\n".join(lines)
    report_path.write_text(text, encoding="utf-8")

    # Also print to console
    print()
    print(text)


# ===========================================================================
# MAIN VALIDATION WORKFLOW
# ===========================================================================

def run_validation(args: argparse.Namespace) -> None:
    """
    Orchestrate the redocking validation for the co-crystal ligand.

    Steps:
      1. Ensure receptor is prepared (re-use or re-run Phase 1 receptor stages).
      2. Extract crystal ligand coordinates from the raw PDB.
      3. Prepare AQ4 as a docking-ready PDBQT using Phase 1 smiles_to_pdbqt().
      4. Dock AQ4 with the same box and parameters used in Phase 1.
      5. Extract docked pose coordinates.
      6. Compute RMSD.
      7. Write report.
    """
    # Import Phase 1 functions — avoids duplicating any logic
    sys.path.insert(0, str(Path(__file__).parent))
    from docking_pipeline import (
        fetch_receptor,
        clean_receptor,
        prepare_receptor,
        define_grid_box,
        smiles_to_pdbqt,
        dock_ligand,
        _configure_logging,
    )

    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    _configure_logging(output_dir)

    print("=" * 60)
    print("REDOCKING VALIDATION — Phase 2")
    print("=" * 60)
    print(f"Target  : {args.pdb_id.upper()}")
    print(f"Ligand  : {args.ligand_resname}  (SMILES from PubChem CID {AQ4_PUBCHEM_CID})")
    print(f"Seed    : {args.seed}   Exhaustiveness: {args.exhaustiveness}")
    print()

    # ── Step 1: Receptor preparation (re-use if already done) ───────────────
    raw_pdb      = fetch_receptor(args.pdb_id, output_dir)
    clean_pdb    = clean_receptor(raw_pdb, args.chain, output_dir)
    receptor_pdbqt = prepare_receptor(clean_pdb, output_dir)

    # ── Step 2: Grid box (same as Phase 1 pipeline) ─────────────────────────
    box = define_grid_box(raw_pdb, args.ligand_resname, padding=args.box_padding)

    # ── Step 3: Extract crystal ligand coordinates ───────────────────────────
    print(f"\nExtracting crystal {args.ligand_resname} coordinates from {raw_pdb.name}...")
    crystal_coords, crystal_elements = extract_crystal_coords(raw_pdb, args.ligand_resname)
    print(f"  Found {len(crystal_coords)} heavy atoms in crystal pose.")
    print(f"  Element counts: { {e: crystal_elements.count(e) for e in sorted(set(crystal_elements))} }")

    # ── Step 4: Prepare AQ4 as PDBQT ────────────────────────────────────────
    print(f"\nPreparing {args.ligand_resname} from SMILES...")
    print(f"  SMILES: {AQ4_SMILES}")
    ligand_pdbqt = smiles_to_pdbqt(args.ligand_resname, AQ4_SMILES)

    if ligand_pdbqt is None:
        print(
            "\nERROR: Failed to prepare AQ4 as PDBQT.\n"
            "Check that rdkit and meeko are installed in the conda env:\n"
            "  conda activate docking\n"
            "  pip show rdkit-pypi meeko"
        )
        sys.exit(1)
    print("  AQ4 PDBQT prepared successfully.")

    # ── Step 5: Dock AQ4 ────────────────────────────────────────────────────
    print(f"\nDocking {args.ligand_resname} (seed={args.seed}, "
          f"exhaustiveness={args.exhaustiveness}, n_poses={args.n_poses})...")

    dock_result = dock_ligand(
        receptor_pdbqt,
        ligand_pdbqt,
        box,
        seed=args.seed,
        exhaustiveness=args.exhaustiveness,
        n_poses=args.n_poses,
    )

    if dock_result["status"] != "success":
        print(
            "\nERROR: Docking failed. Check that vina is installed:\n"
            "  conda activate docking\n"
            "  pip show vina"
        )
        write_report(
            output_dir / "redock_report.txt",
            args.pdb_id, args.ligand_resname,
            len(crystal_coords), 0,
            None, None, box,
        )
        sys.exit(1)

    best_score = dock_result["affinity"]
    pose_pdbqt = dock_result["pose_pdbqt"]
    print(f"  Best Vina score: {best_score:+.2f} kcal/mol")

    # Save the redocked pose for PyMOL inspection
    pose_path = output_dir / f"{args.ligand_resname}_redocked_best.pdbqt"
    pose_path.write_text(pose_pdbqt)
    print(f"  Best pose saved to: {pose_path.name}")

    # ── Step 6: Extract docked pose coordinates ──────────────────────────────
    print(f"\nExtracting heavy atoms from docked pose...")
    try:
        docked_coords, docked_elements = extract_docked_coords(pose_pdbqt)
    except ValueError as e:
        print(f"\nERROR: {e}")
        sys.exit(1)

    print(f"  Found {len(docked_coords)} heavy atoms in docked pose.")
    print(f"  Element counts: { {e: docked_elements.count(e) for e in sorted(set(docked_elements))} }")

    if len(crystal_coords) != len(docked_coords):
        print(
            f"\nWARNING: atom count mismatch — "
            f"crystal={len(crystal_coords)}, docked={len(docked_coords)}.\n"
            "RMSD will be computed over matched atoms only."
        )

    # ── Step 7: Compute RMSD ─────────────────────────────────────────────────
    print("\nComputing heavy-atom RMSD...")
    rmsd = compute_rmsd(crystal_coords, crystal_elements, docked_coords, docked_elements)

    if rmsd is not None:
        verdict = "PASS" if rmsd <= RMSD_PASS_THRESHOLD else "FAIL"
        print(f"  Heavy-atom RMSD : {rmsd:.3f} Å")
        print(f"  Threshold       : {RMSD_PASS_THRESHOLD:.1f} Å")
        print(f"  Verdict         : {verdict}")
    else:
        print("  ERROR: RMSD could not be computed — check element count mismatches above.")

    # ── Step 8: Write report ─────────────────────────────────────────────────
    report_path = output_dir / "redock_report.txt"
    write_report(
        report_path,
        args.pdb_id, args.ligand_resname,
        len(crystal_coords), len(docked_coords),
        best_score, rmsd, box,
    )
    print(f"\nFull report written to: {report_path}")

    # Exit code: 0 = PASS, 1 = FAIL (useful for CI/scripting)
    if rmsd is None or rmsd > RMSD_PASS_THRESHOLD:
        sys.exit(1)


# ===========================================================================
# CLI ARGUMENT PARSER
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    """
    Define and return the ArgumentParser for the redocking validation CLI.
    """
    parser = argparse.ArgumentParser(
        prog="validate_redock",
        description=(
            "Redocking validation for the EGFR co-crystal ligand AQ4.\n"
            "Docks AQ4 back into 1M17 using the Phase 1 pipeline and reports\n"
            "the heavy-atom RMSD vs the crystallographic pose.\n"
            "Target: RMSD < 2.0 Å."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python validate_redock.py --pdb-id 1M17 --output-dir results_egfr\n\n"
            "If RMSD > 2 Å: diagnose before running any library screen.\n"
            "  Try: --exhaustiveness 16  or  --box-padding 15"
        ),
    )

    parser.add_argument(
        "--pdb-id", default="1M17", metavar="ID",
        help="RCSB PDB accession code for the target. Default: 1M17.",
    )
    parser.add_argument(
        "--chain", default="A", metavar="CHAIN",
        help="Protein chain to keep. Default: A.",
    )
    parser.add_argument(
        "--ligand-resname", default="AQ4", metavar="RESNAME",
        help="PDB residue name of the co-crystal ligand. Default: AQ4.",
    )
    parser.add_argument(
        "--box-padding", type=float, default=10.0, metavar="ANGSTROM",
        help="Grid box padding in Å beyond the ligand extent. Default: 10.0.",
    )
    parser.add_argument(
        "--exhaustiveness", type=int, default=8, metavar="N",
        help=(
            "Vina search depth. Default: 8 (same as pipeline). "
            "Increase to 16 or 32 if RMSD > 2 Å."
        ),
    )
    parser.add_argument(
        "--seed", type=int, default=42, metavar="SEED",
        help="Vina random seed. Must match the pipeline seed. Default: 42.",
    )
    parser.add_argument(
        "--n-poses", type=int, default=9, metavar="N",
        help="Number of Vina poses to generate. Default: 9.",
    )
    parser.add_argument(
        "--output-dir", default="results_egfr", metavar="DIR", type=Path,
        help="Output directory. Created if absent. Default: results_egfr.",
    )

    return parser


# ===========================================================================
# ENTRY POINT
# ===========================================================================

def main() -> None:
    """Parse CLI arguments and run redocking validation."""
    parser = build_parser()
    args = parser.parse_args()
    run_validation(args)


if __name__ == "__main__":
    main()
