# CR-Shock Storyline Execution Plan and Live Log

## Scope
- Repository: `/Users/dbf75/Work/Research/AthenaK/athenak-DF`
- Driver brief: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_CR_SHOCK_STORYLINE_IMPLEMENTATION_BRIEF.md`
- Execution mode: implement and validate continuously, updating this file at each step.

## Operator Constraints
- Local machine is a laptop; testing is constrained to **8 CPU ranks** (`np=8`) for runtime runs.
- Local runs are designed to validate correctness and pipeline behavior, not final publication scale.
- HPC portability is mandatory: all scan/analysis tooling must support scaling to large MPI/GPU launches without code changes.

## Allowed Edit Surface
- Primary:
  - `inputs/tests/pic_amr_shock_lb_smoke.athinput`
  - `inputs/tests/pic_amr_shock_lb_publication_local.athinput`
  - `inputs/tests/pic_amr_shock_lb_publication_hpc.athinput`
  - `tst/scripts/particles/pic_amr_shock_lb_smoke.py`
  - `tst/publication/run_pic_shock_scan.py`
  - `tst/publication/pvtk_particles.py`
  - new shock-only scripts under `tst/publication/`
  - shock-specific docs under repo root or `tst/.codex/.../reviews/`
- Secondary only if needed:
  - `src/pgen/tests/orszag_tang.cpp`
  - `src/pgen/AGENTS.md`
  - `inputs/AGENTS.md`
  - `tst/publication/AGENTS.md`

## Branch and Baseline
- Start branch: `dev/PIC`
- Start commit: `a771feadd23d2489d8d3b80daa270983844068e5`
- Working branch for this execution: `codex/pic-shock-storyline`

## Standard Commands
- Build (MPI expected):
  - `cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst`
  - `cmake -S .. -B build -DAthena_ENABLE_MPI=ON`
  - `cmake --build build -j8`
- Smoke anchor module:
  - `cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst`
  - `python3 -m scripts.particles.pic_amr_shock_lb_smoke`
- 8-rank deck sanity pattern:
  - `cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/build/src`
  - `mpiexec -n 8 ./athena -i <deck> job/basename=<run_id>`

## Step Plan (Implementation + Test Gate)

### Step 1: Baseline capture
- [x] Record git head + dirty-state summary.
- [x] Run smoke baseline and archive logs/metrics under:
  - `tst/.codex/pic_publication_runs/<run_id>/reviews/shock_step1/`
- Gate:
  - [x] Smoke remains passing.
  - [x] No new fatal/warning regressions.

### Step 2: Deck uplift (local + hpc)
- [x] Tune local deck for richer diagnostics and laptop-feasible runtime.
- [x] Tune hpc deck for longer evolution and larger dynamic range.
- [x] Ensure cadence enables shock tracking, spectra series, and AMR/LB trends.
- Gate:
  - [x] Local deck completes on `np=8`.
  - [x] HPC deck plumbing sanity run on `np=8` (reduced overrides).
  - [x] Required mesh + particle outputs present.

### Step 3: Diagnostic extraction
- [x] Implement shock diagnostics utility for:
  - shock location proxy (peak `|d rho / dx|`),
  - upstream/downstream windows,
  - CR momentum magnitude + `dN/dp`,
  - `f(x,p)` proxy arrays.
- Gate:
  - [x] End-to-end diagnostics run on local publication run artifacts.
  - [x] Artifacts versioned under run directory.

### Step 4: Controlled scan runner
- [x] Extend scan matrix to explicit Mach / CR loading / resolution grid.
- [x] Add `--tier local|hpc` profiles and machine-readable summaries.
- [x] Ensure launcher abstraction supports HPC migration (MPI and scheduler wrappers).
- Gate:
  - [x] Local mini-scan executes on `np=8`.
  - [x] JSON/CSV summary contains run IDs, parameters, status.

### Step 5: Figure pipeline
- [x] Add shock-only multi-panel figure script:
  - Panel A: `rho`, `|B|`, `|J|`
  - Panel B: CR `x-p` proxy
  - Panel C: `dN/dp` (with optional overlays/segment split)
  - Panel D: AMR/load-balance context
- Gate:
  - [x] Re-running script reproduces outputs.
  - [x] Panels are legible and physically interpretable.

### Step 6: Canonical run selection
- [x] Select one local canonical run.
- [x] Select one HPC candidate run profile.
- [x] Document rationale and limitations.
- Gate:
  - [x] Selection backed by measured metrics.

### Step 7: Final review and handoff
- [x] Re-run smoke anchor.
- [x] Re-run local canonical and one scan subset on `np=8`.
- [x] Produce final report with exact commands and absolute paths.
- Gate:
  - [x] Reproducible from clean checkout instructions.

## Live Progress Log
- 2026-02-08: Initialized execution plan and logging document.
- 2026-02-08: Step 1 baseline captured.
  - Run ID: `pic_shock_storyline_20260208_124724`
  - Artifacts: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step1/`
  - Smoke module result: pass (`serial`, `serial_alt`, `mpi2`, `mpi4`).
  - Added explicit local-policy sanity: `np=8` smoke-deck run (`b1_std_amp=1.0072494745`).
- 2026-02-08: Step 2 deck uplift + validation complete.
  - Local deck: laptop-stable staging profile, rich diagnostics, `np=8` complete run.
  - HPC deck: long-duration/high-dynamic-range profile retained for cluster campaigns.
  - Step 2 artifacts: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step2/`
  - Output presence summary confirms all required mesh fields, particle moments, `pvtk`, and `rst`.
  - AMR churn is deferred in local/plumbing validation (`ncycle_check/refinement_interval=1000`) to avoid known orphan bursts on long laptop runs.
- 2026-02-08: Step 3 diagnostics extraction complete.
  - Added `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/publication/analyze_pic_shock_storyline.py`.
  - Generated `shock_series`, latest profiles, `f(x,p)` proxy, and segmented spectra artifacts.
  - Deterministic rerun check passed (`checksum_diff.txt` empty).
  - Step 3 artifacts: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step3/`
- 2026-02-08: Step 4 controlled scan runner complete.
  - Updated `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/publication/run_pic_shock_scan.py` with tier profiles, launcher templating, `np=8` local default, and JSON/CSV summaries.
  - Local mini-scan (`3` cases) passed on `np=8`:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step4/scan_summary.json`
  - HPC-ready emission with 4096-rank GPU template:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724_hpc_emit/reviews/shock_step4/scan_commands.sh`
- 2026-02-08: Step 5 figure pipeline complete.
  - Added `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/publication/plot_pic_shock_storyline.py`.
  - Generated figure pack:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step5/shock_storyline_multipanel.png`
  - Deterministic figure rerun verified (`figure_sha_diff.txt` empty).
- 2026-02-08: Step 6 canonical selection complete.
  - Selection doc:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step6/canonical_selection.md`
- 2026-02-08: Step 7 final validation complete.
  - Smoke rerun passed:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step7/smoke_rerun_metrics.json`
  - Local canonical rerun passed on `np=8` with clean error scan.
  - One-case scan subset rerun passed:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724_step7_subset/reviews/shock_step4/scan_summary.json`
  - Final handoff report written:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_CR_SHOCK_STORYLINE_FINAL_REPORT.md`
- 2026-02-08: Follow-up figure request complete.
  - Added Athena++-style panel script:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/publication/plot_pic_shock_paper_style.py`
  - Generated additional figure:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step5/shock_storyline_paper_style.png`
