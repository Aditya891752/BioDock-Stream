# Rules for AI Agents
## High-Throughput Molecular Docking Pipeline

**Version:** 1.0  
**Last Updated:** 2026-11-07  
**Applies to:** All AI coding agents (Antigravity, Claude, GPT-4, Gemini, or any future model) working on this project.

> Read this file before writing any code, any config, or any documentation.

---

## Part A: Scientific Integrity Rules

**A1. Work one phase at a time.**  
Complete the current phase, summarise what changed, and explain how to verify it before starting the next phase. Do not combine phases in a single session.

**A2. Never silently change scientific parameters.**  
If you change docking exhaustiveness, box size, box centre, random seed, or scoring function — state the old value, the new value, and the scientific reason. Log the change in the commit message.

**A3. Run the sanity test after any change to:**
- Docking parameters
- Ligand preparation (Meeko settings)
- Receptor preparation (Open Babel flags)

**A4. Never fabricate scientific results.**  
Do not invent docking scores, RMSD values, PDB IDs, ligand names, or ChEMBL IDs. If unsure, say so and tell the developer to verify on rcsb.org or ChEMBL directly.

**A5. If a validation step fails, stop.**  
Do not build Phase N+1 if Phase N validation has failed. Diagnose the failure first.

**A6. Handle per-ligand failures gracefully.**  
One bad molecule (invalid SMILES, prep failure) must never kill a 500-compound run. Log the failure, set `status=prep_failed` or `status=docking_failed`, continue to the next molecule.

---

## Part B: Code Quality Rules

**B1. Add type hints to every function signature.**

**B2. Add a short docstring (one sentence minimum) to every function.**  
Explain what it does and why, not just what the parameters are.

**B3. Write a pytest test for every function with testable logic.**  
Tests go in `tests/test_<module>.py`. Each test must be runnable with `pytest` from the repo root.

**B4. Keep dependencies to `techstack.md`.**  
Ask before adding anything not already listed. Do not silently `pip install` new libraries.

**B5. Prefer readable code over clever code.**  
The developer is a first-year student who must explain every line in an interview. No clever one-liners that require a Python expert to decode.

**B6. Never delete or overwrite result files without asking the developer.**

---

## Part C: Security Rules (apply to web frontend/backend if built)

### Highest-Risk (fix these first)

**C1-HIGH. No API keys or secrets in frontend or client-side JavaScript.**  
Environment variables must only be accessed server-side. Never expose them to the browser bundle.

**C2-HIGH. No `.env` files committed to git.**  
Add `.env`, `.env.local`, `.env.*` to `.gitignore` before the first commit. Use `.env.example` with placeholder values only.

**C3-HIGH. No secrets anywhere in git history.**  
If a secret was committed, rotate it immediately. Use `git filter-repo` to scrub history. Treat the old secret as compromised.

### High Risk

**C4. No overly permissive database rules.**  
Supabase RLS must be enabled on every table. Firebase security rules must deny by default. Never ship `allow read, write: if true`.

**C5. No public storage buckets.**  
S3, Firebase Storage, and Supabase Storage buckets must require authentication. Never set a bucket to world-readable.

**C6. Rate-limiting on all API endpoints and auth routes.**  
Implement rate limiting on login, registration, and any endpoint that calls the docking engine. Use middleware (e.g., `slowapi` for FastAPI).

**C7. Debug mode and dev tools disabled in production.**  
`DEBUG=False` in Django/FastAPI. No `console.log` with internal data in browser. Remove all debug routes before deploy.

**C8. All admin routes require authentication.**  
Admin endpoints must check for a valid, privileged session token server-side. Client-side auth checks alone are not sufficient.

**C9. Sanitize all user input.**  
Any input used in file paths, database queries, or subprocess calls must be validated and sanitized. Use parameterized queries; never concatenate user input into SQL or shell commands.

**C10. Passwords must use strong hashing.**  
Use bcrypt or argon2 (via `passlib`). Never store plain text. Never use MD5 or SHA-1 for passwords.

**C11. Never expose stack traces to users.**  
Catch exceptions at the API boundary. Return a generic error message to the client. Log the full trace server-side only.

**C12. Auth checks must exist server-side.**  
Never rely solely on client-side redirects or hidden UI elements to protect routes. The server must verify the token on every protected request.

---

## Part D: Design Anti-Pattern Rules (UI / Frontend)

The following patterns are **banned**. Do not use them under any circumstances.

### Visual Anti-Patterns

| Banned Pattern | Why |
|---|---|
| Purple-to-blue gradient backgrounds | Overused; looks like a template |
| Gradient hero text | Illegible at small sizes; screams "template" |
| Glassmorphism cards | Low contrast; fails accessibility |
| Coloured border cards | Decorative noise with no semantic meaning |
| 3 icon boxes in a row | Cliché marketing layout |
| Badge above the headline ("New", "Beta") | Cheap attention-grab |
| Fade-in on scroll | Adds no UX value; annoying on revisit |
| Cursor-following beam effect | Gimmicky; hurts performance |
| Button colour fades on hover | Insufficient contrast signal; use solid state change |
| Grain over a gradient | Visual noise with no function |
| Lucide icons used as content (not UI) | Icons are UI chrome, not data |

### Typography Anti-Patterns

| Banned Pattern | Why |
|---|---|
| Inter font everywhere | Bland; no visual hierarchy without pairing |
| Space Grotesk + Instrument Serif pairing | Overused SaaS combo; lacks personality |
| Serif italic accent text | Generic "startup quote" styling |
| Em dashes everywhere (—) | Overused for "premium" tone; use sparingly |

### Code Anti-Patterns

| Banned Pattern | Why |
|---|---|
| Inconsistent spacing | Breaks visual rhythm; use a spacing scale (4px base) |
| Untouched shadcn/ui defaults | Shows no design decision was made |
| Generic placeholder passwords in docs | Creates false sense of security |
| Emojis in headings | Unprofessional in a scientific/technical context |

---

## Part E: Version Control Rules

**E1.** Every commit message must state: what changed, why it changed.  
**E2.** Do not commit generated result files (`results_*/`, `*.pdbqt` outputs).  
**E3.** Do not commit the conda environment binary files. Commit `environment.yml` only.  
**E4.** Use meaningful branch names: `phase-0-setup`, `phase-1-sanity`, etc.
