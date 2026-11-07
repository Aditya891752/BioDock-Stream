# System Architecture
## High-Throughput Molecular Docking Pipeline

**Version:** 1.0  
**Last Updated:** 2026-11-07

---

## 1. System Overview

The project has two distinct layers that are built sequentially and can be deployed independently:

| Layer | What it is | When built |
|---|---|---|
| **Core pipeline** | Python CLI — the scientific computing engine | Phases 0–5 |
| **Web layer** | FastAPI backend + Next.js frontend — browses results | Phase 6+ |

The CLI pipeline is the primary deliverable. The web layer is optional and reads pipeline outputs; it does not duplicate pipeline logic.

---

## 2. CLI Pipeline Architecture

```
INPUT
  PDB ID (e.g. "1M17")
  molecules.csv (name, smiles[, is_active])
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│                  docking_pipeline.py                    │
│                                                         │
│  Stage 1: fetch_receptor()                              │
│    └─ RCSB HTTP API → receptor_raw.pdb                  │
│                                                         │
│  Stage 2: clean_receptor()                              │
│    └─ Biopython → receptor_clean.pdb (1 chain, no HETATM)│
│                                                         │
│  Stage 3: prepare_receptor()                            │
│    └─ Open Babel → receptor.pdbqt                       │
│                                                         │
│  Stage 4: define_grid_box()                             │
│    └─ Biopython → centre(x,y,z), size(x,y,z)           │
│       (derived from co-crystal ligand atom centroids)   │
│                                                         │
│  Stage 5: prepare_ligands()   [parallel, per molecule]  │
│    └─ RDKit + Meeko → ligand_<name>.pdbqt               │
│       invalid SMILES: logged, status=prep_failed        │
│                                                         │
│  Stage 6: dock_ligands()      [parallel, per molecule]  │
│    └─ Vina Python API → best pose, affinity (kcal/mol)  │
│       docking error: logged, status=docking_failed      │
│                                                         │
│  Stage 7: extract_interactions()  [per docked pose]     │
│    └─ Distance-based (Phase 1–4)                        │
│       ProLIF/PLIP fingerprints (Phase 5+)               │
│       → hbond_residues, hydrophobic_residues            │
│                                                         │
│  Stage 8: write_results()                               │
│    └─ pandas → docking_results.csv                      │
│                                                         │
└─────────────────────────────────────────────────────────┘
        │
        ▼
OUTPUT
  results_<target>/
    receptor_clean.pdb
    receptor.pdbqt
    poses/          ← best pose PDBQT per ligand
    docking_results.csv
    run.log
```

### Key Design Choices

| Choice | Reason |
|---|---|
| Per-ligand error isolation | One bad molecule must never kill a 500-compound run |
| Fixed random seed | Reproducible results across re-runs |
| Distance-based interactions first | Simpler; upgrade to ProLIF in Phase 5 |
| Grid box from co-crystal ligand | Standard practice; ensures box covers known binding pocket |
| Parallel ligand prep and docking | Speed; use `concurrent.futures.ProcessPoolExecutor` |

---

## 3. Validation Scripts Architecture

```
validate_redock.py
  └─ Uses: receptor.pdbqt, AQ4 SMILES
  └─ Runs: Vina docking with same parameters as pipeline
  └─ Computes: heavy-atom RMSD (Biopython) vs crystallographic AQ4
  └─ Outputs: RMSD value + pass/fail vs 2 Å threshold

prepare_library.py
  └─ Sources: ChEMBL API (actives, IC50 < 1 µM) + DUD-E decoys download
  └─ Cleans: SMILES validity via RDKit
  └─ Outputs: library.csv (name, smiles, is_active)

analyze_results.py
  └─ Reads: docking_results.csv + library.csv (is_active column)
  └─ Computes: ROC AUC, EF at 1%, Met793 contact check
  └─ Outputs: roc_curve.png, affinity_hist.png, top5_hits.csv
```

---

## 4. Web Layer Architecture (Phase 6+)

```
                        INTERNET
                            │
              ┌─────────────┴─────────────┐
              │                           │
        ┌─────▼─────┐             ┌───────▼──────┐
        │  Vercel   │             │    Render     │
        │ (Frontend)│             │  (Backend)    │
        │           │             │               │
        │ Next.js   │  REST API   │  FastAPI      │
        │ App Router│────────────▶│  (Python)     │
        │           │  JSON       │               │
        └─────┬─────┘             └───────┬───────┘
              │                           │
              │                     ┌─────▼──────┐
              │                     │ File store │
              │                     │ results/   │
              │                     │ (Render    │
              │                     │  disk or   │
              │                     │  S3-compat)│
              │                     └────────────┘
              │
        (Browser)
```

### Frontend (Vercel — Next.js)

| Route | Purpose |
|---|---|
| `/` | Landing: what this tool is |
| `/results` | Table of all docking results, sortable/filterable |
| `/results/[run_id]` | Single run detail: top hits, charts, interaction summary |
| `/api/...` | Next.js API routes for BFF — never expose backend URL to client |

### Backend (Render — FastAPI)

| Endpoint | Method | Purpose |
|---|---|---|
| `/health` | GET | Health check for Render deploy |
| `/runs` | GET | List completed docking runs |
| `/runs/{run_id}` | GET | Results for a specific run |
| `/runs/{run_id}/csv` | GET | Download raw CSV |
| `/dock` | POST | (Future) Submit a new docking job |

**Rule:** The backend URL is an environment variable on Vercel, server-side only. It is never in client JS.

---

## 5. Data Flow for Web Layer

```
1. Developer runs CLI pipeline locally (or on a server)
   → produces results_egfr/docking_results.csv

2. Developer uploads CSV to backend (or backend reads from shared volume)

3. Frontend calls Next.js API route (server-side)
   → Next.js API route calls FastAPI backend
   → FastAPI returns JSON

4. Next.js renders table/charts server-side (SSR) or client-side (CSR)
   → TanStack Table + Recharts display data

5. User interacts: sorts, filters, downloads CSV
```

No real-time updates needed. Polling is acceptable if future job submission is added.

---

## 6. File System Layout

```
docking-pipeline/           ← git repo root
├── docs/                   ← planning documents (this directory)
│   ├── prd.md
│   ├── architecture.md
│   ├── techstack.md
│   ├── aistack.md
│   ├── rules.md
│   ├── design.md
│   ├── task.md
│   ├── memory.md
│   └── architecture_diagram.md
├── docking_pipeline.py     ← main CLI
├── validate_redock.py
├── prepare_library.py
├── analyze_results.py
├── data/
│   ├── molecules_sanity.csv    ← 5-molecule sanity test set
│   └── library.csv             ← full 500+ compound library
├── results_egfr/           ← git-ignored (large binary + CSV)
│   ├── receptor_clean.pdb
│   ├── receptor.pdbqt
│   ├── poses/
│   ├── docking_results.csv
│   └── run.log
├── tests/
│   ├── test_docking_pipeline.py
│   ├── test_validate_redock.py
│   └── test_prepare_library.py
├── frontend/               ← Next.js app (Phase 6+)
│   ├── app/
│   ├── components/
│   └── package.json
├── backend/                ← FastAPI app (Phase 6+)
│   ├── main.py
│   ├── Dockerfile
│   └── requirements.txt
├── environment.yml
├── requirements.txt
├── .gitignore
└── README.md
```

---

## 7. Security Architecture

| Concern | Mitigation |
|---|---|
| Secrets in git | `.env` in `.gitignore`; `.env.example` with placeholders only |
| Backend URL exposed | Server-side environment variable on Vercel only |
| API abuse | Rate limiting on FastAPI endpoints via `slowapi` |
| File path traversal | Validate `run_id` against whitelist of known run directories |
| Unauthenticated admin | FastAPI dependency injection for auth check on any write endpoint |
