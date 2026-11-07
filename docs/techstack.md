# Tech Stack
## High-Throughput Molecular Docking Pipeline

**Version:** 1.0  
**Last Updated:** 2026-11-07

---

## Runtime Environment

| Component | Choice | Reason |
|---|---|---|
| Language | Python 3.11 | Best scientific ecosystem; conda manages complex binary deps |
| Package manager | conda + pip | conda for RDKit/OpenBabel (compiled libs); pip for everything else |
| Environment name | `docking` | Isolated; reproducible across machines |

---

## Core Scientific Libraries

| Library | Version (pin after first working run) | Role |
|---|---|---|
| `rdkit` | conda-forge latest | SMILES parsing, 3D conformer generation, structure sanitization |
| `openbabel` | conda-forge latest | PDB → PDBQT receptor conversion; hydrogen addition |
| `vina` | PyPI (`pip install vina`) | AutoDock Vina Python bindings — the docking engine |
| `meeko` | PyPI | SMILES → PDBQT ligand prep (Vina-compatible) |
| `biopython` | PyPI | PDB file parsing, residue extraction, RMSD calculation |

**Version conflict risk:** Meeko and Vina APIs change between releases. Always check installed versions with `pip show vina meeko` before assuming function signatures. See `rules.md` for handling policy.

---

## Data & Analysis Libraries

| Library | Role |
|---|---|
| `pandas` | Results CSV, library CSV manipulation |
| `numpy` | Coordinate math, distance calculations |
| `matplotlib` | Affinity distribution plots, ROC curves |
| `scikit-learn` | ROC AUC, enrichment factor calculation |

---

## Interaction Analysis (Phase 5 upgrade)

| Library | Role | Install |
|---|---|---|
| `prolif` | Protein-ligand interaction fingerprints | `pip install prolif` |
| `plip` | H-bond, hydrophobic, salt bridge detection | `pip install plip` |

Start with distance-based detection (Phase 1–4). Upgrade to ProLIF/PLIP in Phase 5. Never mix both in the same reported result — pick one and stick with it.

---

## Testing

| Tool | Role |
|---|---|
| `pytest` | Unit and integration tests |
| `pytest-cov` | Coverage reporting |

Every function with testable logic gets a `tests/test_<module>.py` test. See `rules.md` for testing policy.

---

## Development Tools

| Tool | Role |
|---|---|
| `black` | Code formatter (line length 88) |
| `ruff` | Linter |
| `mypy` | Type checking |
| `pre-commit` | Runs black + ruff before each commit |

---

## Visualisation (manual, outside pipeline)

| Tool | Use Case |
|---|---|
| PyMOL | Render top poses, colour by interaction type, export publication images |
| ChimeraX | Alternative; better for large complexes |

These are run manually by the developer after the pipeline outputs poses. The pipeline does not call PyMOL programmatically.

---

## Frontend (Phase 6+, optional)

| Component | Choice | Reason |
|---|---|---|
| Framework | Next.js 14 (App Router) | File-based routing, SSR, deploys to Vercel |
| Styling | Tailwind CSS | Utility-first; avoids component library lock-in |
| Data tables | TanStack Table | Headless; full control over markup |
| Charts | Recharts | React-native; works with Tailwind |
| Hosting | Vercel | Zero-config Next.js deploy |

---

## Backend API (Phase 6+, optional)

| Component | Choice | Reason |
|---|---|---|
| Framework | FastAPI | Async, typed, auto OpenAPI docs |
| Server | Uvicorn | ASGI; pairs with FastAPI |
| Background jobs | Celery + Redis | Long-running docking jobs need a queue |
| Hosting | Render | Simple Docker deploy; free tier available |
| Container | Docker | Reproducible environment for Render |

---

## Installation Commands

```bash
# Create and activate environment
conda create -n docking python=3.11 -y
conda activate docking

# Scientific core (conda for compiled binaries)
conda install -c conda-forge rdkit openbabel -y

# Everything else via pip
pip install vina meeko biopython pandas numpy matplotlib scikit-learn pytest pytest-cov

# Dev tools
pip install black ruff mypy pre-commit

# Phase 5 interaction upgrade (add when ready)
pip install prolif plip
```

After a working run, freeze versions:

```bash
conda env export > environment.yml
pip freeze > requirements.txt
```

---

## Dependency Policy

- Do not add new dependencies without justification in a commit message.
- Ask before adding anything not in this document.
- Never install GUI frameworks (tkinter, PyQt) — the pipeline is CLI-only.
