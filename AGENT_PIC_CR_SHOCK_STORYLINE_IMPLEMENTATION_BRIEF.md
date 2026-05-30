> **HISTORICAL ONLY.** Do not execute commands or infer qualification,
> publication, branch, or Frontier-submission authority from this file. Use
> `tst/publication/PIC_PRODUCTION_READINESS_PLAN.md` as the sole controlling plan.

# CR-Shock MHD-PIC Storyline Implementation Brief (Parallel Agent)

## Mission
Build a publication-grade CR-shock MHD-PIC storyline starting from
`pic_amr_shock_lb_smoke`, without regressing existing proxy/smoke gates.

Target outcome:
1. Longer-duration and higher-dynamic-range shock runs.
2. Rich CR diagnostics for `f(x,p)` and spectra.
3. Controlled scan campaign (Mach, CR loading, resolution).
4. Multi-panel figure pipeline covering fluid, CR phase-space, spectra, and
   AMR/load-balance context.

## Parallel Work Contract
This task is intended to run in parallel with ongoing F05/F07/F08 work.

Do:
1. Work in an isolated branch from current `dev/PIC`.
2. Restrict edits to shock-storyline files listed below.
3. Commit incrementally after each step and keep run IDs in commit bodies.

Do not:
1. Modify Bell/CRSI/CRPAI publication scripts or thresholds.
2. Change shared figure logic for F05/F07/F08.
3. Force-merge unrelated files.

## Allowed Edit Surface
Primary:
1. `inputs/tests/pic_amr_shock_lb_smoke.athinput`
2. `inputs/tests/pic_amr_shock_lb_publication_local.athinput`
3. `inputs/tests/pic_amr_shock_lb_publication_hpc.athinput`
4. `tst/scripts/particles/pic_amr_shock_lb_smoke.py`
5. `tst/publication/run_pic_shock_scan.py`
6. `tst/publication/pvtk_particles.py`
7. New shock-only analysis/plot scripts under `tst/publication/`
8. Shock-specific docs under repo root or `tst/.codex/.../reviews/`

Secondary (only if needed):
1. `src/pgen/tests/orszag_tang.cpp` and `src/pgen/AGENTS.md`
2. `inputs/AGENTS.md`, `tst/publication/AGENTS.md`

## Required Deliverables
1. New publication-local shock deck (CPU-workstation feasible).
2. New publication-hpc shock deck (cluster/GPU intended).
3. Scan driver with explicit grids:
   - Mach number scan,
   - CR loading scan,
   - resolution scan.
4. Post-processing script that emits:
   - shock profile panels (`rho`, `|B|`, `|J|`),
   - CR phase-space proxy (`x-p` from particle VTK),
   - CR spectrum panel (`dN/dp` and optional segmented upstream/downstream),
   - AMR/load-balance context panel (MeshBlock counts/cost/rank stats).
5. One concise report with:
   - run table (deck, ranks, walltime, status),
   - selected canonical run and rationale,
   - remaining risks for HPC campaign.

## Step-by-Step Plan (Mandatory)
Complete each step fully, then review, stage, and commit before moving on.

### Step 1: Baseline capture
1. Record current head and dirty-state summary.
2. Run smoke baseline:
   - `python3 -m scripts.particles.pic_amr_shock_lb_smoke` path via test harness style.
3. Save baseline metrics/log snippets under:
   - `tst/.codex/pic_publication_runs/<run_id>/reviews/shock_step1/`.

Gate:
1. Smoke script remains passing.
2. No runtime fatal/warning regressions.

### Step 2: Deck uplift (local + hpc)
1. Tune local deck for richer diagnostics and moderate runtime.
2. Tune hpc deck for longer evolution and larger dynamic range.
3. Ensure output cadence supports:
   - shock front tracking,
   - CR spectra time series,
   - AMR/load-balance trend extraction.
4. Keep clear comments on intended compute tier.

Gate:
1. Local deck runs to completion on `np2` without fatal.
2. Output files include required fields and particle dumps.

### Step 3: Diagnostic extraction
1. Build parser utilities for:
   - shock location (e.g., peak `|d rho / dx|` proxy),
   - downstream/upstream windows,
   - CR momentum magnitude and histogram bins.
2. Produce reproducible `f(x,p)` proxy arrays and spectra artifacts.

Gate:
1. Diagnostic script runs end-to-end on at least one local publication run.
2. Generated data products are versioned in run artifact directory.

### Step 4: Controlled scan runner
1. Implement explicit scan matrix in `run_pic_shock_scan.py`:
   - Mach: at least 3 points,
   - CR loading: at least 3 points (e.g., ppc or injection proxy),
   - resolution: at least 2 points.
2. Support `--tier local` and `--tier hpc` profiles.
3. Emit machine-readable summary (JSON/CSV).

Gate:
1. Local mini-scan executes successfully.
2. Summary table contains run IDs, key params, pass/fail.

### Step 5: Figure pipeline
1. Generate multi-panel figure template:
   - Panel A: density + B/current structure,
   - Panel B: CR phase-space proxy,
   - Panel C: CR spectrum (`dN/dp`) with optional temporal overlays,
   - Panel D: AMR/load-balance context (block counts/cost or rank histograms).
2. Add CLI for selecting run ID and output directory.

Gate:
1. Figure script reproduces same outputs when rerun.
2. Panels are legible and physically interpretable.

### Step 6: Canonical run selection
1. Pick one local canonical run and one hpc-target candidate.
2. Document why selected (signal quality, runtime feasibility, stability).
3. Note limitations and required HPC follow-up.

Gate:
1. Selection rationale documented and defensible with metrics.

### Step 7: Final review and handoff
1. Re-run smoke + local publication canonical + one scan subset.
2. Confirm no regression to existing particle smoke gate.
3. Produce final handoff report with exact file paths and commands.

Gate:
1. All step outputs reproducible from clean checkout instructions.

## Testing Requirements
At each step touching runtime behavior:
1. `python3 -m py_compile` on changed Python files.
2. `python3 -m flake8` on changed Python files.
3. Re-run affected tests:
   - `particles.pic_amr_shock_lb_smoke` minimum.
4. For changed decks, include at least:
   - serial (`np1`) sanity,
   - MPI (`np2`) parity sanity.

## Commit Discipline
Use one commit per completed step:
1. `PIC-SHOCK Step 1: baseline snapshot`
2. `PIC-SHOCK Step 2: publication deck uplift`
3. `PIC-SHOCK Step 3: CR diagnostics extraction`
4. `PIC-SHOCK Step 4: parameter scan runner`
5. `PIC-SHOCK Step 5: multi-panel figure pipeline`
6. `PIC-SHOCK Step 6: canonical run selection`
7. `PIC-SHOCK Step 7: final validation + handoff`

Each commit body must include:
1. run ID(s),
2. exact commands,
3. key metrics,
4. known caveats.

## Reporting Format (Final)
Deliver one markdown report:
1. Executive summary.
2. What changed (files + purpose).
3. Scan matrix and outcomes.
4. Canonical runs and why.
5. Figure pack path.
6. Residual risks and next HPC actions.

Use absolute paths in all report references.
