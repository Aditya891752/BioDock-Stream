"""
docking_pipeline.py — High-Throughput Molecular Docking Pipeline
=================================================================
Phase 1: Full implementation of all 8 pipeline stages.

Scientific summary
------------------
This pipeline takes a protein drug target (from RCSB PDB) and a CSV of
small molecules encoded as SMILES strings, docks each molecule into the
protein's binding pocket using AutoDock Vina, and reports which protein
residues the top molecules contact.

Stage overview
--------------
1. fetch_receptor      — download PDB file from RCSB
2. clean_receptor      — keep one chain, remove water and heteroatoms
3. prepare_receptor    — convert PDB → PDBQT (Vina's input format)
4. define_grid_box     — locate the binding pocket from the co-crystal ligand
5. smiles_to_pdbqt     — convert each SMILES to a 3D PDBQT ligand file
6. dock_ligand         — run AutoDock Vina for each ligand
7. extract_interactions— identify H-bonds and hydrophobic contacts (distance-based)
8. write_results       — produce a ranked CSV

Usage
-----
python docking_pipeline.py --pdb-id 1M17 --library data/molecules_sanity.csv

Dependencies (must be installed in the docking conda env)
---------------------------------------------------------
conda install -c conda-forge rdkit openbabel
pip install vina meeko biopython pandas numpy tqdm
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Module-level logger
# All functions in this file use this logger.
# Handlers (console + file) are attached in _configure_logging().
# ---------------------------------------------------------------------------
logger = logging.getLogger("docking")


# ===========================================================================
# LOGGING SETUP
# ===========================================================================

def _configure_logging(output_dir: Path) -> None:
    """
    Attach a console handler (INFO) and a file handler (DEBUG) to the logger.

    Why two handlers:
      - Console shows the researcher what is happening in real time.
      - run.log captures DEBUG-level detail for diagnosing failures later.

    This function is safe to call multiple times — it checks if handlers
    are already attached before adding new ones.
    """
    if logger.handlers:
        return  # already configured; do not add duplicate handlers

    logger.setLevel(logging.DEBUG)

    # Console: show INFO and above in a clean format
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter("%(levelname)-8s %(message)s"))

    # File: capture everything for post-run diagnosis
    log_path = output_dir / "run.log"
    file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)-8s %(name)s — %(message)s")
    )

    logger.addHandler(console)
    logger.addHandler(file_handler)


# ===========================================================================
# STAGE 1 — Fetch receptor
# ===========================================================================

def fetch_receptor(pdb_id: str, output_dir: Path) -> Path:
    """
    Download a PDB structure file from RCSB and save it locally.

    RCSB (the Protein Data Bank) hosts all public protein structures.
    Every structure has a 4-character ID (e.g. "1M17" = EGFR with erlotinib).
    The download URL is always: https://files.rcsb.org/download/<ID>.pdb

    We use Python's built-in urllib so we do not need the `requests` library.
    If the file already exists (from a previous run), we skip the download.

    Parameters
    ----------
    pdb_id : str
        Four-character RCSB PDB accession code, case-insensitive.
    output_dir : Path
        Directory where the downloaded file is saved.

    Returns
    -------
    Path
        Path to the saved PDB file: output_dir/<PDB_ID>_raw.pdb

    Raises
    ------
    urllib.error.URLError
        If the download fails (network error or invalid PDB ID).
        Check the PDB ID on rcsb.org before concluding it is a network issue.
    """
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    out_path = output_dir / f"{pdb_id}_raw.pdb"

    if out_path.exists():
        logger.info("Receptor already downloaded: %s (skipping)", out_path.name)
        return out_path

    logger.info("Downloading %s from RCSB PDB...", url)
    urllib.request.urlretrieve(url, out_path)
    size_kb = out_path.stat().st_size // 1024
    logger.info("Saved %s (%d KB)", out_path.name, size_kb)
    return out_path


# ===========================================================================
# STAGE 2 — Clean receptor
# ===========================================================================

def clean_receptor(raw_pdb_path: Path, chain: str, output_dir: Path) -> Path:
    """
    Keep only one protein chain and remove water/heteroatoms.

    Why cleaning is needed:
      - PDB files often contain multiple protein chains, which creates
        ambiguity about which surface to dock into.
      - Water molecules (HOH) and heteroatoms (ligands, ions) interfere
        with receptor preparation and inflate the grid box.
      - ANISOU records (anisotropic displacement parameters) are not
        needed for docking and can confuse Open Babel.

    What we keep: ATOM records for the specified chain (standard amino acids).
    What we remove: HETATM records (all heteroatoms including water),
                    all other chains.

    We use Biopython's PDB parser and Select class rather than text
    manipulation — fixed-width PDB columns make string slicing fragile.

    Parameters
    ----------
    raw_pdb_path : Path
        The downloaded PDB file (may contain many chains and HETATM records).
    chain : str
        Single chain letter to keep, e.g. "A".
    output_dir : Path
        Directory for the cleaned file.

    Returns
    -------
    Path
        Path to receptor_clean.pdb — one chain, protein atoms only.
    """
    from Bio import PDB  # Biopython; imported here for a clear missing-dep error

    out_path = output_dir / "receptor_clean.pdb"
    if out_path.exists():
        logger.info("Clean receptor already exists (skipping): %s", out_path.name)
        return out_path

    logger.info("Cleaning receptor (chain=%s, removing water and heteroatoms)...", chain)

    # PDBParser: QUIET=True suppresses harmless "discontinuous chain" warnings
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(raw_pdb_path))

    class ProteinChainSelect(PDB.Select):
        """
        Biopython Select subclass: keep only protein residues from one chain.

        PDBIO calls accept_chain() and accept_residue() on every entity.
        Return 1 to keep, 0 to discard.
        """
        def accept_chain(self, c) -> int:
            return 1 if c.id == chain else 0

        def accept_residue(self, residue) -> int:
            # Biopython encodes the residue type in the hetflag:
            #   " "  = standard amino acid (keep)
            #   "W"  = water molecule (discard)
            #   "H_" = other heteroatom (discard)
            hetflag, _, _ = residue.get_id()
            return 1 if hetflag == " " else 0

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), ProteinChainSelect())

    n_atoms = sum(1 for line in out_path.read_text().splitlines()
                  if line.startswith("ATOM"))
    logger.info("Clean receptor saved: %s (%d ATOM records)", out_path.name, n_atoms)
    return out_path


# ===========================================================================
# STAGE 3 — Prepare receptor (PDB → PDBQT)
# ===========================================================================

def prepare_receptor(clean_pdb_path: Path, output_dir: Path) -> Path:
    """
    Convert the cleaned PDB receptor to PDBQT format using Open Babel.

    What PDBQT is:
      PDBQT is a PDB extension that AutoDock Vina requires. It adds two
      columns to each atom record: (1) a Gasteiger partial charge and
      (2) an AutoDock atom type (e.g., "C.3", "N.am", "OA").
      These are used by Vina's scoring function.

    Open Babel flags used:
      -xr   → rigid receptor (no rotatable bonds are assigned to the receptor;
               only the ligand is flexible during docking).
      --partialcharge gasteiger → assign Gasteiger electrostatic charges.

    Parameters
    ----------
    clean_pdb_path : Path
        Single-chain, water-free PDB file from clean_receptor().
    output_dir : Path
        Directory where receptor.pdbqt is written.

    Returns
    -------
    Path
        Path to receptor.pdbqt.

    Raises
    ------
    RuntimeError
        If obabel is not found or exits with a non-zero code.
        Make sure Open Babel is installed: conda install -c conda-forge openbabel
    """
    out_path = output_dir / "receptor.pdbqt"
    if out_path.exists():
        logger.info("Receptor PDBQT already exists (skipping): %s", out_path.name)
        return out_path

    logger.info("Converting receptor to PDBQT with Open Babel...")

    cmd = [
        "obabel",
        str(clean_pdb_path),
        "-O", str(out_path),
        "-xr",
        "--partialcharge", "gasteiger",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"Open Babel failed (exit code {result.returncode}).\n"
            f"Command: {' '.join(cmd)}\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}\n\n"
            "Fix: conda install -c conda-forge openbabel"
        )

    if not out_path.exists() or out_path.stat().st_size == 0:
        raise RuntimeError(
            f"Open Babel ran (exit 0) but produced no output at {out_path}.\n"
            f"stderr: {result.stderr}\n"
            "This often means the PDB file had no readable ATOM records."
        )

    logger.info(
        "Receptor PDBQT saved: %s (%d KB)",
        out_path.name,
        out_path.stat().st_size // 1024,
    )
    return out_path


# ===========================================================================
# STAGE 4 — Define grid box
# ===========================================================================

def define_grid_box(
    raw_pdb_path: Path,
    ligand_resname: str,
    padding: float = 10.0,
) -> dict[str, float]:
    """
    Define the Vina docking grid box from the co-crystal ligand position.

    What the grid box is:
      Vina only searches for poses *inside* a user-defined 3D box. If the
      true binding pocket is outside the box, Vina will never find it.
      The standard approach is to centre the box on the co-crystal ligand —
      the molecule that was crystallised inside the pocket. This guarantees
      the box covers the right region.

    How we compute it:
      1. Find all heavy atoms of the specified ligand in the RAW PDB file
         (the cleaned file has no HETATM, so we must use the raw one).
      2. Compute the centroid (geometric centre) of those atoms.
      3. Measure the extent of the ligand along each axis.
      4. Add padding on both sides so the box is larger than the ligand
         and there is room for binding poses to rotate.

    Parameters
    ----------
    raw_pdb_path : Path
        The original downloaded PDB file (still contains HETATM records).
    ligand_resname : str
        Three-letter PDB residue name of the co-crystal ligand, e.g. "AQ4".
    padding : float
        Extra space (Å) added on each side of the ligand extent. Default 10 Å.

    Returns
    -------
    dict[str, float]
        Keys: center_x, center_y, center_z, size_x, size_y, size_z (all in Å).

    Raises
    ------
    ValueError
        If no atoms are found for the given resname. This means the resname
        is wrong or you accidentally passed the cleaned PDB (no HETATM).
    """
    from Bio import PDB

    logger.info(
        "Defining grid box from co-crystal ligand '%s' (padding=%.1f Å)...",
        ligand_resname, padding,
    )

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("raw", str(raw_pdb_path))

    # Collect heavy-atom coordinates for the specified ligand residue
    coords: list[list[float]] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip().upper() == ligand_resname.strip().upper():
                    for atom in residue:
                        if atom.element and atom.element.strip().upper() != "H":
                            coords.append(list(atom.get_vector()))
                        elif not atom.element:
                            # element field absent in some PDB files; skip H by name prefix
                            if not atom.get_name().strip().startswith("H"):
                                coords.append(list(atom.get_vector()))

    if not coords:
        raise ValueError(
            f"No heavy atoms found for ligand resname '{ligand_resname}' in {raw_pdb_path}.\n"
            "Checklist:\n"
            "  1. Use the RAW PDB file, not the cleaned one (clean has no HETATM).\n"
            "  2. Verify the exact 3-letter resname on rcsb.org for your PDB ID.\n"
            f"  3. For 1M17 the correct resname is 'AQ4'. You passed: '{ligand_resname}'."
        )

    arr = np.array(coords)   # shape: (n_heavy_atoms, 3)
    centroid = arr.mean(axis=0)
    extent = arr.max(axis=0) - arr.min(axis=0)
    box_size = extent + 2.0 * padding   # padding on both sides of each axis

    box = {
        "center_x": float(centroid[0]),
        "center_y": float(centroid[1]),
        "center_z": float(centroid[2]),
        "size_x":   float(box_size[0]),
        "size_y":   float(box_size[1]),
        "size_z":   float(box_size[2]),
    }

    logger.info(
        "Box centre: (%.2f, %.2f, %.2f) Å",
        box["center_x"], box["center_y"], box["center_z"],
    )
    logger.info(
        "Box size:   (%.2f, %.2f, %.2f) Å  [ligand extent + 2×%.1f Å padding]",
        box["size_x"], box["size_y"], box["size_z"], padding,
    )
    return box


# ===========================================================================
# STAGE 5 — Ligand preparation (SMILES → PDBQT)
# ===========================================================================

def smiles_to_pdbqt(name: str, smiles: str) -> Optional[str]:
    """
    Convert a SMILES string to a 3D PDBQT string ready for Vina.

    Steps
    -----
    1. Parse the SMILES with RDKit and sanitize it (checks valence, aromaticity).
       If the SMILES is invalid, return None and log a warning — do not crash.
    2. Add explicit hydrogens. Meeko and Vina need H atoms to assign
       hydrogen-bond donor/acceptor types and partial charges correctly.
    3. Generate a 3D conformer with ETKDGv3.
       ETKDG (Experimental Torsion Knowledge Distance Geometry) uses dihedral
       angle preferences from the Cambridge Structural Database, so the
       starting geometry is chemically sensible (not random).
    4. Minimize energy with MMFF94 force field. This removes bad bond lengths
       and angles that could make Vina score the pose badly from the start.
    5. Use Meeko to assign AutoDock atom types and produce the PDBQT string.
       Meeko handles both rigid and flexible torsion bonds in the ligand.

    Meeko API note:
    Meeko changed its public API between v0.3 and v0.4. This function tries
    the new API first and falls back to the old one automatically.

    Parameters
    ----------
    name : str
        Molecule name — used only for logging messages.
    smiles : str
        SMILES string for the molecule.

    Returns
    -------
    str | None
        PDBQT string on success; None if any step fails.
        Failure is logged as a WARNING, not raised as an exception.
        The pipeline will mark this molecule as status=prep_failed and continue.
    """
    # Import inside the function so the module still loads even if rdkit/meeko
    # are not installed — the caller handles the None return value.
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        logger.error(
            "rdkit not found. Install with: conda install -c conda-forge rdkit"
        )
        return None

    try:
        import meeko
    except ImportError:
        logger.error("meeko not found. Install with: pip install meeko")
        return None

    # ── Step 1: Parse and sanitize SMILES ───────────────────────────────────
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning(
            "[%s] Invalid SMILES — could not parse: %s", name, smiles[:60]
        )
        return None

    # ── Step 2: Add explicit hydrogens ──────────────────────────────────────
    mol = Chem.AddHs(mol)

    # ── Step 3: Generate 3D conformer with ETKDGv3 ──────────────────────────
    params = AllChem.ETKDGv3()
    params.randomSeed = 42   # fixed seed = reproducible geometry
    embed_result = AllChem.EmbedMolecule(mol, params)
    if embed_result == -1:
        logger.warning(
            "[%s] 3D conformer generation failed (ETKDGv3 returned -1). "
            "This can happen with very rigid or unusual ring systems.",
            name,
        )
        return None

    # ── Step 4: Minimize with MMFF94 (fallback to UFF) ──────────────────────
    # MMFFOptimizeMolecule returns 0=success, 1=not converged, -1=MMFF unavailable
    ff_result = AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    if ff_result == -1:
        # MMFF94 not parameterised for this molecule (rare; some exotic atoms)
        logger.debug("[%s] MMFF94 unavailable — trying UFF force field.", name)
        uff_result = AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
        if uff_result == -1:
            logger.warning(
                "[%s] No force field available for minimization. "
                "Using unoptimized geometry (docking quality may be reduced).",
                name,
            )

    # ── Step 5: Meeko preparation ────────────────────────────────────────────
    return _meeko_to_pdbqt(name, mol)


def _meeko_to_pdbqt(name: str, mol) -> Optional[str]:
    """
    Convert an RDKit Mol with 3D coordinates to a PDBQT string via Meeko.

    Meeko changed its API significantly in v0.4. We try the new API first
    and fall back to the old one so the pipeline works with either version.

    New API (Meeko >= 0.4):
      MoleculePreparation().prepare(mol) → list of mol_setups
      PDBQTWriterLegacy.write_string(mol_setups[0]) → (pdbqt_str, is_ok, error_msg)

    Old API (Meeko <= 0.3):
      MoleculePreparation().prepare(mol)
      preparator.write_pdbqt_string() → pdbqt_str

    Parameters
    ----------
    name : str
        Molecule name for logging.
    mol : rdkit.Chem.Mol
        RDKit molecule with 3D coordinates and explicit hydrogens.

    Returns
    -------
    str | None
        PDBQT string or None if Meeko fails.
    """
    from meeko import MoleculePreparation

    try:
        # ── Try new Meeko >= 0.4 API ─────────────────────────────────────────
        from meeko import PDBQTWriterLegacy

        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)

        if not mol_setups:
            logger.warning("[%s] Meeko returned no molecule setups.", name)
            return None

        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setups[0])

        if not is_ok:
            logger.warning("[%s] Meeko PDBQT write error: %s", name, error_msg)
            return None

        return pdbqt_string

    except (ImportError, AttributeError):
        # PDBQTWriterLegacy doesn't exist → old Meeko API
        pass
    except Exception as e:
        logger.warning("[%s] Meeko new API failed (%s), trying old API.", name, e)

    try:
        # ── Fallback: old Meeko <= 0.3 API ──────────────────────────────────
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        return preparator.write_pdbqt_string()

    except Exception as e:
        logger.warning("[%s] Meeko old API also failed: %s", name, e)
        return None


# ===========================================================================
# STAGE 6 — Docking
# ===========================================================================

def dock_ligand(
    receptor_pdbqt_path: Path,
    ligand_pdbqt_str: str,
    box: dict[str, float],
    seed: int = 42,
    exhaustiveness: int = 8,
    n_poses: int = 9,
) -> dict:
    """
    Dock a single ligand into the receptor using AutoDock Vina.

    How Vina works
    --------------
    Vina uses a Monte Carlo stochastic search algorithm. It randomly moves
    and rotates the ligand inside the grid box, and evaluates each pose with
    a knowledge-based scoring function. The score (kcal/mol) estimates
    binding free energy — more negative = predicted tighter binding.

    Important caveats:
    - Vina scores are predictions, not measurements. They correlate poorly
      with measured IC50 values in absolute terms, but rank-ordering is
      useful for comparing molecules against the same target.
    - `exhaustiveness=8` is Vina's default and adequate for screening.
      Increase to 16–32 for final validation of top hits.
    - Fixing the seed makes docking reproducible. Never change the seed
      without documenting it in the commit message (rule A2 in rules.md).

    Parameters
    ----------
    receptor_pdbqt_path : Path
        Path to the prepared receptor PDBQT file.
    ligand_pdbqt_str : str
        PDBQT content of the ligand (from smiles_to_pdbqt).
    box : dict[str, float]
        Grid box from define_grid_box(): center_x/y/z and size_x/y/z.
    seed : int
        Random seed for Vina. Default: 42.
    exhaustiveness : int
        Search depth. Default: 8.
    n_poses : int
        Number of poses Vina returns. We use only the best (rank 1). Default: 9.

    Returns
    -------
    dict with keys:
        "affinity"   : float | None — best pose energy in kcal/mol
        "pose_pdbqt" : str | None   — PDBQT text of the best pose
        "status"     : str          — "success" or "docking_failed"
    """
    try:
        from vina import Vina
    except ImportError:
        logger.error("vina not installed. Run: pip install vina")
        return {"affinity": None, "pose_pdbqt": None, "status": "docking_failed"}

    try:
        # verbosity=0 silences Vina's own output — we log what we need ourselves
        v = Vina(sf_name="vina", seed=seed, verbosity=0)

        # Load the rigid receptor PDBQT
        v.set_receptor(rigid_pdbqt_filename=str(receptor_pdbqt_path))

        # Load the ligand from a PDBQT string (no temp file needed)
        # Vina expects a list of strings even for a single ligand
        v.set_ligand_from_string([ligand_pdbqt_str])

        # Compute the affinity map (the grid of pre-computed energies inside the box)
        v.compute_vina_maps(
            center=[box["center_x"], box["center_y"], box["center_z"]],
            box_size=[box["size_x"], box["size_y"], box["size_z"]],
        )

        # Run the actual docking search
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

        # Retrieve results: energies() returns an array of shape (n_poses, 6)
        # Column 0 is the total binding energy in kcal/mol
        energies = v.energies(n_poses=1)
        best_affinity = float(energies[0][0])

        # Retrieve the best pose as a PDBQT string
        best_pose_pdbqt = v.poses(n_poses=1)

        return {
            "affinity": best_affinity,
            "pose_pdbqt": best_pose_pdbqt,
            "status": "success",
        }

    except Exception as e:
        logger.debug("Vina docking exception: %s", e, exc_info=True)
        return {"affinity": None, "pose_pdbqt": None, "status": "docking_failed"}


# ===========================================================================
# WORKER FUNCTION FOR SCREENING (Stages 5, 6, 7 per ligand)
# ===========================================================================

def process_single_ligand(
    name: str,
    smiles: str,
    receptor_pdbqt_path: Path,
    clean_pdb_path: Path,
    box: dict[str, float],
    seed: int = 42,
    exhaustiveness: int = 8,
    n_poses: int = 9,
    poses_dir: Optional[Path] = None,
) -> dict:
    """
    Execute Stages 5, 6, and 7 on a single ligand in a self-contained manner.

    This function is designed to be executed either serially or in parallel
    across multiple worker processes via concurrent.futures.ProcessPoolExecutor.

    Parameters
    ----------
    name : str
        Molecule identifier.
    smiles : str
        SMILES string.
    receptor_pdbqt_path : Path
        Prepared receptor PDBQT.
    clean_pdb_path : Path
        Receptor clean PDB (for distance calculations).
    box : dict
        Grid box bounds and center.
    seed : int
        Deterministic random seed.
    exhaustiveness : int
        Vina search exhaustiveness.
    n_poses : int
        Number of output poses.
    poses_dir : Optional[Path]
        Directory to save the best pose PDBQT file.

    Returns
    -------
    dict
        Completed result record ready for checkpointing and CSV reporting.
    """
    res: dict = {"name": name, "smiles": smiles}

    # Stage 5: Ligand preparation
    lig_pdbqt = smiles_to_pdbqt(name, smiles)
    if lig_pdbqt is None:
        res.update({
            "affinity": None,
            "pose_pdbqt": None,
            "hbond_residues": [],
            "hydrophobic_residues": [],
            "status": "prep_failed",
        })
        return res

    # Stage 6: Docking
    dock_res = dock_ligand(
        receptor_pdbqt_path=receptor_pdbqt_path,
        ligand_pdbqt_str=lig_pdbqt,
        box=box,
        seed=seed,
        exhaustiveness=exhaustiveness,
        n_poses=n_poses,
    )
    res.update(dock_res)

    if dock_res["status"] != "success":
        res.update({"hbond_residues": [], "hydrophobic_residues": []})
        return res

    # Save best pose file if directory provided
    if poses_dir is not None:
        try:
            poses_dir.mkdir(parents=True, exist_ok=True)
            pose_file = poses_dir / f"{name}_best.pdbqt"
            pose_file.write_text(dock_res["pose_pdbqt"], encoding="utf-8")
        except Exception as e:
            logger.debug("[%s] Failed to save pose file: %s", name, e)

    # Stage 7: Interaction detection
    interactions = extract_interactions_distance(dock_res["pose_pdbqt"], clean_pdb_path)
    res.update(interactions)
    return res


# ===========================================================================
# STAGE 7 — Interaction extraction (distance-based)
# ===========================================================================

def _parse_atom_records(text: str) -> list[dict]:
    """
    Parse ATOM and HETATM lines from PDB or PDBQT text.

    PDB format columns (1-indexed, fixed-width):
      1–6   : Record type ("ATOM  " or "HETATM")
      13–16 : Atom name
      18–20 : Residue name
      22    : Chain ID
      23–26 : Residue sequence number
      31–38 : X coordinate
      39–46 : Y coordinate
      47–54 : Z coordinate
      77–78 : Element symbol (may be absent in older PDB files)

    PDBQT adds two more columns after position 78 (charge, atom type),
    but we only read up to position 78.

    Parameters
    ----------
    text : str
        Contents of a PDB or PDBQT file.

    Returns
    -------
    list[dict]
        One dict per ATOM/HETATM line with keys:
        name, resname, chain, resnum, x, y, z, element.
    """
    atoms = []
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            atom_name = line[12:16].strip()
            resname   = line[17:20].strip()
            chain     = line[21].strip() if len(line) > 21 else ""
            resnum    = int(line[22:26].strip())
            x         = float(line[30:38])
            y         = float(line[38:46])
            z         = float(line[46:54])

            # Derive element from column 77-78 if available; otherwise from atom name.
            # Atom names like "CA", "CB", "NZ" — first alphabetic character is the element.
            if len(line) >= 78 and line[76:78].strip():
                element = line[76:78].strip()[0].upper()
            else:
                # Strip leading digits (e.g. "1HB" → "HB" → "H")
                stripped = atom_name.lstrip("0123456789")
                element = stripped[0].upper() if stripped else "X"

            atoms.append({
                "name":    atom_name,
                "resname": resname,
                "chain":   chain,
                "resnum":  resnum,
                "x": x, "y": y, "z": z,
                "element": element,
            })
        except (ValueError, IndexError):
            # Malformed line — skip silently; log at DEBUG only
            logger.debug("Could not parse PDB/PDBQT line: %r", line[:60])
            continue
    return atoms


def extract_interactions_distance(
    pose_pdbqt_str: str,
    receptor_pdb_path: Path,
    hbond_cutoff: float = 3.5,
    hydrophobic_cutoff: float = 4.5,
) -> dict[str, list[str]]:
    """
    Identify protein residues that contact the docked ligand pose.

    This is a simplified distance-based method adequate for Phase 1–4.
    Phase 5 upgrades this to ProLIF fingerprints (more accurate because
    ProLIF also checks donor/acceptor roles and geometric angles).

    H-bond rule (simplified):
      Any N or O atom in the ligand within hbond_cutoff Å of any N or O
      atom in a receptor residue. This *over-detects* (ignores angle and
      H-atom position) but is a valid first-pass screen.

    Hydrophobic contact rule:
      Any C or S atom in the ligand within hydrophobic_cutoff Å of any
      C or S atom in a receptor residue.

    Biology check (EGFR specific):
      Met793 is the hinge residue. Known inhibitors (erlotinib, gefitinib,
      osimertinib) always H-bond with Met793. If your top hits miss Met793,
      the docking setup is likely wrong.

    Parameters
    ----------
    pose_pdbqt_str : str
        PDBQT text of the best docked pose (from dock_ligand()).
    receptor_pdb_path : Path
        Path to receptor_clean.pdb (the cleaned single-chain protein).
    hbond_cutoff : float
        Max N/O–N/O distance for a putative H-bond in Å. Default: 3.5.
    hydrophobic_cutoff : float
        Max C/S–C/S distance for a hydrophobic contact in Å. Default: 4.5.

    Returns
    -------
    dict[str, list[str]]
        {
          "hbond_residues":       ["MET793_A", "THR790_A", ...],
          "hydrophobic_residues": ["LEU844_A", "VAL726_A", ...]
        }
        Residue labels are formatted as "<resname><resnum>_<chain>".
    """
    # Parse ligand atoms from the pose PDBQT
    ligand_atoms = _parse_atom_records(pose_pdbqt_str)
    if not ligand_atoms:
        logger.debug("No atoms in pose PDBQT — cannot extract interactions.")
        return {"hbond_residues": [], "hydrophobic_residues": []}

    # Parse receptor atoms from the cleaned PDB
    receptor_text = receptor_pdb_path.read_text(encoding="utf-8", errors="ignore")
    receptor_atoms = _parse_atom_records(receptor_text)

    # Ligand atom coordinate arrays, split by element for vectorised distance calc
    lig_polar = np.array([
        [a["x"], a["y"], a["z"]]
        for a in ligand_atoms if a["element"] in ("N", "O")
    ], dtype=float)

    lig_hydrophobic = np.array([
        [a["x"], a["y"], a["z"]]
        for a in ligand_atoms if a["element"] in ("C", "S")
    ], dtype=float)

    hbond_hits: set[str] = set()
    hydrophobic_hits: set[str] = set()

    for rec in receptor_atoms:
        rec_xyz = np.array([rec["x"], rec["y"], rec["z"]], dtype=float)
        label = f"{rec['resname']}{rec['resnum']}_{rec['chain']}"

        # H-bond: receptor polar atom vs all ligand polar atoms
        if rec["element"] in ("N", "O") and lig_polar.size > 0:
            dists = np.linalg.norm(lig_polar - rec_xyz, axis=1)
            if dists.min() <= hbond_cutoff:
                hbond_hits.add(label)

        # Hydrophobic: receptor C/S atom vs all ligand C/S atoms
        if rec["element"] in ("C", "S") and lig_hydrophobic.size > 0:
            dists = np.linalg.norm(lig_hydrophobic - rec_xyz, axis=1)
            if dists.min() <= hydrophobic_cutoff:
                hydrophobic_hits.add(label)

    return {
        "hbond_residues":       sorted(hbond_hits),
        "hydrophobic_residues": sorted(hydrophobic_hits),
    }


# ===========================================================================
# STAGE 8 — Write results CSV
# ===========================================================================

def write_results(results: list[dict], output_path: Path) -> None:
    """
    Write all docking results to a CSV sorted by affinity (best first).

    Every molecule appears in the output — including failures.
    Failure rows have NaN for affinity and empty strings for interactions.
    This ensures the CSV is always a complete record of the run, not just
    the successes, so the developer can see the full picture.

    Output columns:
      name                 : molecule name from input CSV
      smiles               : input SMILES string
      affinity_kcal_mol    : best Vina score (NaN for failures)
      hbond_residues       : semicolon-separated list (e.g. "MET793_A;THR790_A")
      hydrophobic_residues : semicolon-separated list
      status               : "success", "prep_failed", or "docking_failed"

    Parameters
    ----------
    results : list[dict]
        One result dict per molecule (from the docking loop in run_pipeline).
    output_path : Path
        Where to write the CSV file.
    """
    rows = []
    for r in results:
        rows.append({
            "name":                 r.get("name", ""),
            "smiles":               r.get("smiles", ""),
            "affinity_kcal_mol":    r.get("affinity"),          # None → NaN
            "hbond_residues":       ";".join(r.get("hbond_residues", [])),
            "hydrophobic_residues": ";".join(r.get("hydrophobic_residues", [])),
            "status":               r.get("status", "unknown"),
        })

    df = pd.DataFrame(rows)

    # Sort: successes by affinity (ascending = most negative = best first),
    # then all failures appended below
    successes = df[df["status"] == "success"].sort_values("affinity_kcal_mol")
    failures  = df[df["status"] != "success"]
    df_out = pd.concat([successes, failures], ignore_index=True)

    df_out.to_csv(output_path, index=False)
    logger.info(
        "Results written: %s  (%d succeeded, %d failed)",
        output_path,
        len(successes),
        len(failures),
    )


# ===========================================================================
# CHECKPOINT HELPERS (Phase 4 resume feature — basic support from Phase 1)
# ===========================================================================

def _load_checkpoint(output_dir: Path) -> dict:
    """
    Load a JSON checkpoint of molecules already docked in a previous run.

    The checkpoint maps molecule_name → result_dict. If no checkpoint
    file exists yet, return an empty dict (fresh run).

    Why we do this: a 500-molecule screen can take hours. If the run is
    interrupted (power cut, error), the researcher should not have to
    start over from molecule 1.
    """
    path = output_dir / "progress.json"
    if path.exists():
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return {}


def _save_checkpoint(output_dir: Path, progress: dict) -> None:
    """
    Overwrite progress.json with the current progress dict.

    Called after every molecule (success or failure) so a re-run with
    --resume can skip anything already done.
    """
    path = output_dir / "progress.json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(progress, f, indent=2)


# ===========================================================================
# MAIN PIPELINE ORCHESTRATION
# ===========================================================================

def run_pipeline(args: argparse.Namespace) -> None:
    """
    Run all pipeline stages in order using the parsed CLI arguments.

    Receptor stages (1–4) run once and stop the whole pipeline on failure —
    there is no point docking ligands if the receptor is broken.

    Ligand stages (5–7) run once per molecule. Failures are caught,
    logged as warnings, stored with the appropriate status code, and skipped.
    One bad SMILES never kills the full run.
    """
    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    _configure_logging(output_dir)

    logger.info("=" * 60)
    logger.info("Molecular Docking Pipeline — Phase 1")
    logger.info("=" * 60)
    logger.info("Target PDB     : %s", args.pdb_id)
    logger.info("Library        : %s", args.library)
    logger.info("Chain          : %s  |  Co-crystal lig: %s", args.chain, args.ligand_resname)
    logger.info("Exhaustiveness : %d  |  Seed: %d  |  Poses: %d",
                args.exhaustiveness, args.seed, args.n_poses)
    logger.info("Box padding    : %.1f Å  |  Workers: %d", args.box_padding, args.workers)
    logger.info("Output dir     : %s", output_dir)
    logger.info("=" * 60)

    # ── Receptor preparation (stages 1–4) ───────────────────────────────────

    if args.skip_receptor_prep:
        raw_pdb      = output_dir / f"{args.pdb_id.upper()}_raw.pdb"
        clean_pdb    = output_dir / "receptor_clean.pdb"
        receptor_pdbqt = output_dir / "receptor.pdbqt"
        if not receptor_pdbqt.exists():
            raise FileNotFoundError(
                f"--skip-receptor-prep was set but {receptor_pdbqt} does not exist.\n"
                "Remove the flag to run receptor preparation."
            )
        logger.info("Skipping receptor prep — using existing %s", receptor_pdbqt)
    else:
        raw_pdb      = fetch_receptor(args.pdb_id, output_dir)        # Stage 1
        clean_pdb    = clean_receptor(raw_pdb, args.chain, output_dir) # Stage 2
        receptor_pdbqt = prepare_receptor(clean_pdb, output_dir)       # Stage 3

    box = define_grid_box(raw_pdb, args.ligand_resname, args.box_padding)  # Stage 4

    # ── Load ligand library ─────────────────────────────────────────────────

    library_df = pd.read_csv(args.library)
    required_cols = {"name", "smiles"}
    missing = required_cols - set(library_df.columns)
    if missing:
        raise ValueError(
            f"Library CSV is missing required columns: {missing}\n"
            "Required: name, smiles (optional: is_active)"
        )
    total = len(library_df)
    logger.info("Library loaded: %d molecules", total)

    # ── Resume checkpoint ───────────────────────────────────────────────────

    progress = _load_checkpoint(output_dir) if args.resume else {}
    if progress:
        logger.info("Resuming: %d molecules already done, skipping them.", len(progress))

    # ── Per-ligand docking loop (stages 5–7) ───────────────────────────────

    results: list[dict] = list(progress.values())
    pending_rows = [row for row in library_df.itertuples(index=False) if str(row.name) not in progress]
    logger.info("Molecules to process: %d (already completed: %d)", len(pending_rows), len(progress))

    poses_dir = output_dir / "poses"
    poses_dir.mkdir(parents=True, exist_ok=True)

    # Use tqdm progress bar if available; fall back to plain progress
    try:
        from tqdm import tqdm
        has_tqdm = True
    except ImportError:
        has_tqdm = False
        logger.info("(Install tqdm for a live progress bar: pip install tqdm)")

    if args.workers > 1 and len(pending_rows) > 0:
        logger.info("Executing parallel virtual screen across %d worker processes...", args.workers)
        import concurrent.futures

        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
            # Map futures to molecule name
            future_to_mol = {
                executor.submit(
                    process_single_ligand,
                    name=str(row.name),
                    smiles=str(row.smiles),
                    receptor_pdbqt_path=receptor_pdbqt,
                    clean_pdb_path=clean_pdb,
                    box=box,
                    seed=args.seed,
                    exhaustiveness=args.exhaustiveness,
                    n_poses=args.n_poses,
                    poses_dir=poses_dir,
                ): str(row.name)
                for row in pending_rows
            }

            pbar = tqdm(total=len(pending_rows), desc="Screening", unit="mol") if has_tqdm else None

            for future in concurrent.futures.as_completed(future_to_mol):
                name = future_to_mol[future]
                try:
                    result = future.result()
                except Exception as e:
                    logger.error("[%s] Worker exception: %s", name, e)
                    result = {
                        "name": name,
                        "smiles": "",
                        "affinity": None,
                        "pose_pdbqt": None,
                        "hbond_residues": [],
                        "hydrophobic_residues": [],
                        "status": "docking_failed",
                    }

                results.append(result)
                progress[name] = result
                _save_checkpoint(output_dir, progress)

                if pbar:
                    pbar.update(1)
                else:
                    if len(results) % 25 == 0 or len(results) == total:
                        logger.info("Progress: %d / %d completed", len(results), total)

            if pbar:
                pbar.close()

    else:
        # Serial execution (single worker or debugging)
        iterator = tqdm(pending_rows, desc="Screening", unit="mol") if has_tqdm else iter(pending_rows)

        for row in iterator:
            name   = str(row.name)
            smiles = str(row.smiles)

            result = process_single_ligand(
                name=name,
                smiles=smiles,
                receptor_pdbqt_path=receptor_pdbqt,
                clean_pdb_path=clean_pdb,
                box=box,
                seed=args.seed,
                exhaustiveness=args.exhaustiveness,
                n_poses=args.n_poses,
                poses_dir=poses_dir,
            )

            status = result.get("status")
            if status == "success":
                affinity = result.get("affinity")
                hbonds = result.get("hbond_residues", [])
                met793 = any("MET793" in r for r in hbonds)
                logger.info(
                    "%-32s  %+6.2f kcal/mol  H-bonds: %-3s  Met793: %s",
                    name, affinity,
                    str(len(hbonds)),
                    "YES" if met793 else "no",
                )
            else:
                logger.warning("%-32s  %s", name, status)

            results.append(result)
            progress[name] = result
            _save_checkpoint(output_dir, progress)

    # ── Stage 8: Write results CSV ──────────────────────────────────────────
    output_csv = output_dir / "docking_results.csv"
    write_results(results, output_csv)

    # ── Summary ─────────────────────────────────────────────────────────────
    successful = [r for r in results if r.get("status") == "success"]
    failed     = [r for r in results if r.get("status") != "success"]

    logger.info("")
    logger.info("=" * 60)
    logger.info("SUMMARY")
    logger.info("=" * 60)
    logger.info("Total molecules : %d", len(results))
    logger.info("Successful docks: %d", len(successful))
    logger.info("Failed          : %d", len(failed))

    if successful:
        ranked = sorted(successful, key=lambda r: r["affinity"])
        logger.info("")
        logger.info("Top 5 hits:")
        for i, r in enumerate(ranked[:5], 1):
            met793 = any("MET793" in res for res in r.get("hbond_residues", []))
            logger.info(
                "  %d. %-30s  %+.2f kcal/mol  Met793: %s",
                i, r["name"], r["affinity"],
                "YES" if met793 else "no",
            )

        # Sanity check: at least one known inhibitor should be near the top
        known_inhibitors = {"erlotinib", "gefitinib"}
        top_names = {r["name"] for r in ranked[:3]}
        if known_inhibitors & top_names:
            logger.info("")
            logger.info(
                "SANITY CHECK PASSED: known inhibitor(s) %s are in top 3.",
                known_inhibitors & top_names,
            )
        else:
            logger.warning(
                "SANITY CHECK: known inhibitors NOT in top 3 of this run. "
                "Top 3: %s. Review box definition and receptor prep.",
                [r["name"] for r in ranked[:3]],
            )

    logger.info("")
    logger.info("Results CSV : %s", output_csv)
    logger.info("Log file    : %s", output_dir / "run.log")
    logger.info("=" * 60)


# ===========================================================================
# CLI ARGUMENT PARSER
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    """
    Define and return the ArgumentParser for the pipeline CLI.

    Every argument has a detailed help string explaining what it does
    and why it exists — required since the developer must explain each
    one in an interview.
    """
    parser = argparse.ArgumentParser(
        prog="docking_pipeline",
        description=(
            "High-throughput molecular docking pipeline.\n"
            "Fetches a protein from RCSB PDB, prepares a ligand library,\n"
            "docks every molecule with AutoDock Vina, and reports ranked results."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Sanity run (5 molecules)\n"
            "  python docking_pipeline.py --pdb-id 1M17 --library data/molecules_sanity.csv\n\n"
            "  # Full library screen with 4 CPU cores\n"
            "  python docking_pipeline.py --pdb-id 1M17 --library data/library.csv --workers 4\n\n"
            "  # Resume an interrupted run\n"
            "  python docking_pipeline.py --pdb-id 1M17 --library data/library.csv --resume\n"
        ),
    )

    # Required
    parser.add_argument(
        "--pdb-id", required=True, metavar="ID",
        help="Four-character RCSB PDB accession code. Example: 1M17 (EGFR + erlotinib).",
    )
    parser.add_argument(
        "--library", required=True, metavar="CSV", type=Path,
        help="Ligand library CSV. Required columns: name, smiles. Optional: is_active.",
    )

    # Receptor options
    parser.add_argument(
        "--chain", default="A", metavar="CHAIN",
        help="PDB chain to keep during receptor cleaning. All others are removed. Default: A.",
    )
    parser.add_argument(
        "--ligand-resname", default="AQ4", metavar="RESNAME",
        help=(
            "PDB residue name of the co-crystal ligand — used to locate the binding pocket. "
            "Default: AQ4 (erlotinib analogue in 1M17)."
        ),
    )

    # Docking parameters
    parser.add_argument(
        "--exhaustiveness", type=int, default=8, metavar="N",
        help="Vina search depth. Higher = slower but more thorough. Default: 8.",
    )
    parser.add_argument(
        "--seed", type=int, default=42, metavar="SEED",
        help="Random seed for Vina. Fixing it makes docking reproducible. Default: 42.",
    )
    parser.add_argument(
        "--n-poses", type=int, default=9, metavar="N",
        help="Number of poses Vina generates per ligand. Only the best is reported. Default: 9.",
    )
    parser.add_argument(
        "--box-padding", type=float, default=10.0, metavar="ANGSTROM",
        help="Extra space (Å) added to each side of the grid box. Default: 10.0.",
    )

    # Execution options
    parser.add_argument(
        "--workers", type=int, default=1, metavar="N",
        help="Parallel docking processes (Phase 4). Default: 1 (serial).",
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume interrupted run — skips molecules already in progress.json.",
    )
    parser.add_argument(
        "--output-dir", default="results_egfr", metavar="DIR", type=Path,
        help="Directory for all outputs. Created if absent. Default: results_egfr.",
    )
    parser.add_argument(
        "--skip-receptor-prep", action="store_true",
        help="Skip stages 1–3 if receptor.pdbqt already exists in output-dir.",
    )

    return parser


# ===========================================================================
# ENTRY POINT
# ===========================================================================

def main() -> None:
    """
    Parse CLI arguments and run the pipeline.

    argparse handles --help and argument validation automatically.
    Scientific errors (bad PDB ID, missing ligand) raise with descriptive
    messages so the developer knows exactly what to fix.
    """
    parser = build_parser()
    args = parser.parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()
