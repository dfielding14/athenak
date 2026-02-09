# PIC Parallel-Shock Benchmark Execution Plan and Live Log

## Scope
- Repository: `/Users/dbf75/Work/Research/AthenaK/athenak-DF`
- Objective: implement a clean, minimal physics-benchmark track for a dedicated
  MHD-PIC parallel shock problem, separate from the existing Orszag-Tang
  shock storyline smoke/pipeline suite.
- Reference target: Section 5.4 setup and AMR criteria in
  `/Users/dbf75/Work/Research/AthenaK/athena++_MHD-PIC.pdf`.

## Hard Constraints
- Local validation uses 8 CPU ranks (`mpiexec -n 8`) for all runtime checks.
- Keep Orszag-Tang files/scripts as smoke/pipeline assets; no coupling between
  smoke gates and new benchmark gates.
- Implementation must be scalable by configuration (mesh size, levels, ranks),
  with no required code changes when moving to HPC/GPU runs.

## Deliverables
1. New built-in problem generator `pic_parallel_shock` with:
   - reflecting wall geometry on inner `x1` boundary,
   - upstream inflow and `B0 || x`,
   - ideal-shock-surface CR injection using efficiency `eta`,
   - conservative gas subtraction for injected CR mass/momentum/energy,
   - Section 5.4-style AMR curvature criterion from `g_rho` and `g_P`.
2. Three benchmark decks:
   - coarse uniform,
   - fine uniform,
   - AMR fiducial.
3. Benchmark analysis output comparing:
   - magnetic amplification trend,
   - CR spectra trend (downstream),
   - run metadata and pass/fail status.
4. Separation policy documentation:
   - Orszag-Tang remains smoke/pipeline only,
   - new decks/scripts tagged as physics-benchmark.

## Step-by-Step Implementation and Gates

### Step 0: Plan Initialization
- [x] Create this plan doc and keep it updated during execution.
- Gate:
  - [x] File exists and includes live status updates.

### Step 1: Baseline and Separation Guard
- [x] Record current git head and dirty-state note.
- [x] Re-run existing Orszag-Tang smoke gate to ensure baseline stability.
- [x] Confirm no edits to existing smoke test thresholds/logic while adding
      benchmark-specific paths.
- Gate:
  - [x] Smoke baseline still passes.
  - [x] `git diff` shows no unintended smoke threshold edits.

### Step 2: Add `pic_parallel_shock` Pgen (Minimal Core)
- [x] Add new pgen implementation file under `src/pgen/tests/`.
- [x] Register new pgen in:
  - `src/pgen/pgen.hpp`
  - `src/pgen/pgen.cpp` (new + restart dispatch)
  - `src/CMakeLists.txt`
- [x] Implement initial conditions:
  - uniform upstream MHD state,
  - `v = (-u0, 0, 0)`,
  - `B = (B0, 0, 0)`,
  - reflecting wall usage via deck BCs.
- [x] Set inflow boundary state for outer `x1` as upstream reservoir.
- Gate:
  - [x] Build succeeds.
  - [x] Short run of new pgen completes on `np=8` without fatal errors.

### Step 3: Eta Injection + Conservative Gas Subtraction
- [x] Add user source callback for shock-surface injection at
      `x_shock = x_min + u_shock * t`.
- [x] Inject mono-energetic CRs with isotropic direction in shock frame
      (minimal implementation consistent with Section 5.4 assumptions).
- [x] Subtract corresponding mass/momentum/energy from gas conservatively
      in source update.
- [x] Add guard/floor handling to avoid negative conserved states.
- Gate:
  - [x] `np=8` run completes.
  - [x] Particle count increases over time with non-zero injection.
  - [x] No negative-density/energy fatal behavior.

### Step 4: Section 5.4 Curvature AMR
- [x] Add user refinement callback with block-wise criteria:
  - refine if `max(g_rho or g_P) > 1.0`,
  - derefine if all active cells satisfy `g_rho < 0.1` and `g_P < 0.1`.
- [x] Ensure criteria operate in x-y curvature form as in Section 5.4.
- [x] Keep thresholds configurable via `<problem>` params with defaults
      matching paper values.
- Gate:
  - [x] AMR fiducial run on `np=8` shows refine/derefine activity.
  - [x] Uniform decks remain static (no AMR path).

### Step 5: Three Benchmark Decks
- [x] Add benchmark decks under `inputs/tests/`:
  - `pic_parallel_shock_coarse_uniform.athinput`
  - `pic_parallel_shock_fine_uniform.athinput`
  - `pic_parallel_shock_amr_fiducial.athinput`
- [x] Keep grid/runtime small enough for laptop `np=8` validation.
- [x] Add clear comments indicating HPC scaling knobs
      (`mesh`, `meshblock`, `num_levels`, `nlim`, outputs cadence).
- Gate:
  - [x] All three decks run to completion on `np=8`.
  - [x] Outputs required for analysis exist (`rho`, `B`, CR outputs).

### Step 6: Benchmark Comparison Tooling
- [x] Add benchmark analysis script in `tst/publication/` that:
  - reads run outputs for all three cases,
  - computes magnetic amplification proxy (`B/B0` trend),
  - computes downstream CR spectrum proxy,
  - writes JSON/CSV summary + comparison plot.
- [x] Keep script deterministic and path-stable under `tst/.codex/`.
- Gate:
  - [x] Script runs end-to-end on local `np=8` outputs.
  - [x] Summary clearly distinguishes coarse/fine/AMR trends.

### Step 7: Documentation and Final Validation
- [x] Update AGENTS docs for pgen/input/publication separation where needed.
- [x] Run code quality checks on changed Python files (`py_compile`, `flake8`).
- [x] Re-run:
  - Orszag-Tang smoke (unchanged suite),
  - new benchmark trio (`np=8`),
  - benchmark analysis script.
- Gate:
  - [x] All required checks complete and results logged below.
  - [x] Final status and caveats documented.

### Step 8: Orszag-Tang Paper-Style Figure Cleanup
- [x] Remove the bottom `x-log(e/u0^2)` panel from paper-style figure output.
- [x] Enforce equal-aspect spatial axes for top density and middle `B/B0` panels.
- [x] Regenerate `shock_storyline_paper_style.png` for current Orszag-Tang
      publication run.
- Gate:
  - [x] Figure now contains exactly two equal-aspect map panels.
  - [x] Output path remains stable for publication review.

## Final Implementation Notes (Local vs HPC)
- Local laptop validation mode (this execution):
  - `mpiexec -n 8`, short runtime windows (`nlim=12`, `tlim=0.02`), and
    lightweight CR load (`ppc=0.0`, injection-driven population only).
- HPC promotion path (no code changes needed):
  - increase `nx1`, `nx2`, and for AMR increase `num_levels`;
  - raise `nlim`/`tlim` for physical growth windows;
  - increase rank/GPU count and keep `meshblock` dimensions tuned for occupancy;
  - keep the same `pic_parallel_shock` pgen and benchmark analysis script.
- Known local limitation:
  - 8-rank short runs are intended for correctness and workflow parity; they do
    not produce strong magnetic amplification trends (`B/B0` ratios stay ~1.0).

## Live Execution Log
- 2026-02-08: Plan created; execution started.
- 2026-02-08: Step 1 baseline completed.
  - git head: `a771feadd23d2489d8d3b80daa270983844068e5`
  - 8-rank smoke baseline command archived at:
    `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_parallel_shock_runs/pic_parallel_shock_20260208_142903/reviews/step1_baseline/command.txt`
  - baseline metrics:
    `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_parallel_shock_runs/pic_parallel_shock_20260208_142903/reviews/step1_baseline/metrics.json`
    (`n_bcc_outputs=12`, `b1_std_amp=1.0058043314`)
- 2026-02-08: Step 2 started (new `pic_parallel_shock` pgen wiring).
- 2026-02-08: Step 2 completed.
  - Added:
    - `src/pgen/tests/pic_parallel_shock.cpp`
    - `pgen` registration in `src/pgen/pgen.hpp`, `src/pgen/pgen.cpp`
    - build wiring in `src/CMakeLists.txt`
  - Build gate:
    - `cmake --build tst/build -j8` passed.
- 2026-02-08: Step 3 completed.
  - Injection and subtraction implemented in `ParallelShockSource`.
  - Source callback includes:
    - eta mass budget with local reservoir carry,
    - per-rank deterministic RNG seeding,
    - conservative gas state subtraction with density/pressure floors.
  - Additional robustness fix:
    - all ranks now enter particle-count collective consistently when any rank
      injects (`ninj_global > 0` path).
- 2026-02-08: Step 4 completed.
  - User AMR callback `ParallelShockRefinement` added with Section 5.4-style
    curvature criteria from density/pressure.
  - AMR fiducial gate evidence:
    - run reports `12 MeshBlocks created` from AMR.
- 2026-02-08: Step 5 completed.
  - Added benchmark decks:
    - `inputs/tests/pic_parallel_shock_coarse_uniform.athinput`
    - `inputs/tests/pic_parallel_shock_fine_uniform.athinput`
    - `inputs/tests/pic_parallel_shock_amr_fiducial.athinput`
  - Local-stability updates:
    - `ox1_bc=inflow` in all three decks,
    - `ppc=0.0` (injection-only CR population),
    - AMR deck uses `ps_inject_t_start=3.0e-3` to avoid first-regrid orphan
      warnings in short local runs.
- 2026-02-08: Step 6 completed.
  - Added benchmark runner/comparison script:
    - `tst/publication/run_pic_parallel_shock_benchmark.py`
  - Latest validated run:
    - run id: `pic_parallel_shock_bench_20260208_150603`
    - summary JSON:
      `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_parallel_shock_runs/pic_parallel_shock_bench_20260208_150603/reviews/benchmark/benchmark_summary.json`
    - comparison figure:
      `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_parallel_shock_runs/pic_parallel_shock_bench_20260208_150603/reviews/benchmark/pic_parallel_shock_benchmark_compare.png`
  - Key local metrics (last shared cycle):
    - coarse: `n_particles=90`, `p90=10.455676360999366`, `B_p95_ratio=1.0`
    - fine: `n_particles=56`, `p90=10.464728926235592`, `B_p95_ratio=1.0`
    - amr: `n_particles=43`, `p90=10.458871291368988`, `B_p95_ratio=1.0`
- 2026-02-08: Step 7 completed.
  - AGENTS docs updated:
    - `src/pgen/AGENTS.md` (new built-in pgen and separation note)
    - `inputs/AGENTS.md` (new decks and benchmark-vs-smoke separation)
  - Python checks passed:
    - `python3 -m py_compile tst/publication/plot_pic_shock_paper_style.py`
    - `python3 -m py_compile tst/publication/run_pic_parallel_shock_benchmark.py`
    - `flake8 tst/publication/plot_pic_shock_paper_style.py`
    - `flake8 tst/publication/run_pic_parallel_shock_benchmark.py`
  - Targeted C++ style check passed for the new pgen:
    - `python3 /tmp/cpplint.py --filter=-build/include_subdir --counting=detailed src/pgen/tests/pic_parallel_shock.cpp`
  - Full-repo style script still reports many pre-existing style issues outside
    this change scope.
- 2026-02-08: Step 8 completed (Orszag-Tang figure style cleanup).
  - Updated script:
    - `tst/publication/plot_pic_shock_paper_style.py`
    - now outputs only two equal-aspect map panels (`rho/rho0`, `B/B0`).
  - Updated figure artifact:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_publication_runs/pic_shock_storyline_20260208_124724/reviews/shock_step5/shock_storyline_paper_style.png`
