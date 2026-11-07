# AI Stack
## Tools, Models & Agents Used in This Project

**Version:** 1.0  
**Last Updated:** 2026-11-07

---

## Overview

This document tracks every AI tool, model, and agent involved in building, maintaining, or extending this project. It is distinct from the scientific algorithms (AutoDock Vina is not AI — it is a physics-based scoring function).

---

## 1. Primary Coding Agent

| Property | Value |
|---|---|
| Tool | Antigravity (Google DeepMind) |
| Model | Claude Sonnet 4.6 (Thinking) |
| Role | Primary pair-programmer; writes code, documents, fixes bugs |
| Scope | All phases 0–7 |
| Rules file | `rules.md` |
| Session memory | `memory.md` |

**What it must do:**
- Work one phase at a time
- State every scientific choice explicitly
- Never fabricate docking scores, RMSD values, or PDB IDs
- Run sanity tests after any change to docking, ligand prep, or receptor prep

**What it must not do:**
- Silently change docking parameters or box sizes
- Add dependencies not in `techstack.md` without asking
- Delete or overwrite result files without confirmation

---

## 2. Scientific Algorithms (not AI, but tracked here for completeness)

| Tool | Type | Role |
|---|---|---|
| AutoDock Vina | Physics-based docking | Scores ligand poses in protein pocket |
| RDKit conformer engine | ETKDG algorithm | Generates 3D coordinates from SMILES |
| Open Babel | Cheminformatics | PDB → PDBQT format conversion |
| Meeko | Ligand prep | SMILES → Vina-ready PDBQT |
| ProLIF / PLIP (Phase 5) | Fingerprint-based interaction | H-bonds, hydrophobic contacts, salt bridges |

These produce deterministic outputs given the same inputs and seed. They are not stochastic models.

---

## 3. Data Sources (AI-adjacent — LLM should not hallucinate these)

| Source | Access Method | Used For |
|---|---|---|
| RCSB PDB | HTTP REST (`https://files.rcsb.org/download/`) | Fetching protein structure files |
| ChEMBL | REST API or downloaded SDF | EGFR actives for enrichment validation |
| DUD-E | Downloaded from dude.docking.org | Decoy compounds for enrichment validation |

**Rule:** Never hardcode PDB IDs or ligand names without verifying them on rcsb.org. Never guess ChEMBL IDs.

---

## 4. Future AI Integrations (not yet built)

| Tool | Potential Use | Phase |
|---|---|---|
| Gemini API | Natural-language query over docking results in web dashboard | 6+ |
| AlphaFold (via EBI API) | Predict structure of mutant EGFR targets | Future |
| AlphaGenome | Variant impact scoring if genomic context added | Future |

**Rule:** These are not currently in scope. Do not build toward them without explicit instruction.

---

## 5. Version Tracking

| AI Tool | Current Version | Date Confirmed |
|---|---|---|
| AutoDock Vina (Python) | TBD — run `pip show vina` | — |
| Meeko | TBD — run `pip show meeko` | — |
| RDKit | TBD — run `conda list rdkit` | — |

Fill in after Phase 0 environment setup.

---

## 6. Failure Policy for AI-Generated Code

If the AI agent generates code that:

1. Produces scientifically wrong results (scores, coordinates) → stop, diagnose, correct before continuing
2. Introduces a dependency not in `techstack.md` → flag it; do not silently install
3. Silently changes a parameter (box size, seed, exhaustiveness) → revert and re-run
4. Fabricates a result that cannot be reproduced → treat as a bug, not a warning

The developer is responsible for understanding every line well enough to explain it in an interview.
