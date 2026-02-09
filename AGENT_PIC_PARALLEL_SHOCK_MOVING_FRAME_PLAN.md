# PIC Parallel Shock Moving-Frame/Recentering Plan

## Context and Goal
- Branch snapshot committed before this plan:
  - `codex/pic-shock-storyline`
  - commit `1dd3dd0e91082a782776765c1cbc8b6f036b8f8d`
- Current issue:
  - Coupled run is now numerically stable to `t=3` on 8 ranks, but the shock
    advects to the right edge (`x_shock ~ 3.18` in a `[0,4]` box), leaving
    insufficient residence time in-domain for Section 5.4-style amplification.
- Goal:
  - Add a clean, minimal moving-frame/recentering capability for
    `pic_parallel_shock` so we can run much longer without shock loss, then
    reassess magnetic amplification trends.

## What TRML Does (and What We Reuse)

### Observed in `origin/TRML`
- Hook infrastructure:
  - `src/pgen/pgen.hpp`: `UserWorkInLoopFnPtr`, `bool user_work_in_loop`
  - `src/pgen/pgen.cpp`: parse `problem/user_work_in_loop`, pointer checks
  - `src/hydro/hydro_tasks.cpp` and `src/mhd/mhd_tasks.cpp`: add `WorkInLoop`
    task in `after_timeintegrator`.
- TRML frame tracking implementation:
  - `src/pgen/TRML.cpp`: `UserWorkInLoop -> FrameTracking`.
  - Computes tracked-gas mean velocity via temperature-window reduction.
  - Applies a global Galilean velocity shift to primitives, then `PrimToCons`.

### TRML limitations for this shock benchmark
- Tracking criterion (temperature windows) is problem-specific and noisy.
- It does not exploit known analytic shock speed.
- It is intrusive across Hydro and MHD task classes.
- It does not target PIC coupling details (inflow moment/injection consistency).

## Design Choice for This Branch

### Chosen approach (minimal and clean)
- Add a generic user work-in-loop hook once per cycle.
- Implement a `pic_parallel_shock`-specific analytic frame controller that uses
  known shock speed rather than shock finding.
- First implementation scope:
  - uniform-grid runs only (`adaptive == false` required when enabled),
  - MHD-PIC path only,
  - deterministic, rank-independent control law.

### Not in first pass
- AMR + frame tracking simultaneously.
- Shock-finder-based tracking.
- Complex moving-window remap/regridding of full state.

## Proposed Runtime Model

### New controls under `<problem>`
- `ps_enable_frame_tracking` (bool, default `false`)
- `ps_frame_t_start` (Real): start time for frame control
- `ps_frame_t_ramp` (Real): smooth ramp duration
- `ps_frame_vfrac` (Real): fraction of analytic shock speed to follow
- `ps_frame_dv_max` (Real): max per-step velocity increment
- `ps_frame_apply_to_particles` (bool, default `true`)
- `ps_frame_apply_to_inflow` (bool, default `true`)
- `ps_frame_require_uniform` (bool, default `true`)
- `ps_frame_diag_dcycle` (int): diagnostic cadence for printed frame metrics

### Analytic controller
- Define analytic shock speed from current pgen:
  - `v_shock = 0.5*(gamma-1)*u0`.
- Define desired frame velocity schedule:
  - `V_des(t) = 0` for `t < t_start`.
  - Ramps to `ps_frame_vfrac * v_shock` over `ps_frame_t_ramp`.
- At each cycle (work-in-loop call):
  - compute `dV = V_des(t+dt) - V_des(t)`,
  - clip `dV` by `ps_frame_dv_max`,
  - apply this incremental global shift.
- Determinism benefit:
  - no shock-finder reduction,
  - no accumulated hidden state required for restart parity.

## Code Integration Plan (Step-by-Step)

### Step 0: Hook plumbing
- Files:
  - `src/pgen/pgen.hpp`
  - `src/pgen/pgen.cpp`
  - `src/driver/driver.cpp`
- Changes:
  - Add `UserWorkInLoopFnPtr` and `user_work_in_loop` flag.
  - Parse `problem/user_work_in_loop` in pgen constructors.
  - Validate pointer enrollment if flag is true.
  - Invoke callback once per cycle in `Driver::Execute()` immediately after
    `ExecuteTaskList(..., "after_timeintegrator", 1)` and before time increment.
- Test gate (8 CPUs):
  - build and run baseline deck with hook disabled; results must be unchanged.

### Step 1: Shock pgen enrollment and no-op work function
- File:
  - `src/pgen/tests/pic_parallel_shock.cpp`
- Changes:
  - Add `ParallelShockWorkInLoop(Mesh*)`.
  - Enroll via `user_work_in_loop = true` and `user_work_in_loop_func`.
  - Parse new frame knobs but keep effect disabled by default.
- Test gate (8 CPUs):
  - `ps_enable_frame_tracking=false` run must match pre-change behavior.

### Step 2: Fluid-frame incremental shift
- File:
  - `src/pgen/tests/pic_parallel_shock.cpp`
- Changes:
  - Implement analytic `dV` update each cycle.
  - Shift fluid x-velocity in primitives (`w0(...,IVX,...) += dV`).
  - Rebuild conservative state (`PrimToCons`) over active zones.
  - Leave magnetic field unchanged.
- Test gate (8 CPUs):
  - short run with nonzero `dV` confirms expected momentum/energy change sign.
  - no NaN/negative-state failures.

### Step 3: PIC consistency updates
- File:
  - `src/pgen/tests/pic_parallel_shock.cpp`
- Changes:
  - If enabled, shift particle velocities (`IPVX += dV`) for all local particles.
  - Update inflow conservative boundary state (`u_in(IM1)`, `u_in(IEN)`) when
    `ps_frame_apply_to_inflow=true`.
  - Update injection velocity construction to include frame velocity consistently.
- Test gate (8 CPUs):
  - coupled short run completes.
  - particle counts and injection remain finite.
  - no MPI communication regressions.

### Step 4: Uniform-grid guardrails and diagnostics
- File:
  - `src/pgen/tests/pic_parallel_shock.cpp`
- Changes:
  - If frame tracking enabled and `adaptive==true` with
    `ps_frame_require_uniform=true`, fail fast with clear reason string.
  - Add concise per-`ps_frame_diag_dcycle` diagnostics:
    - current `V_des`, applied `dV`, predicted `x_shock_model`.
- Test gate (8 CPUs):
  - explicit AMR+frame run fails deterministically with correct message.
  - uniform run prints stable diagnostics.

### Step 5: Long-run validation to target times
- Runs:
  - baseline coupled `t=3` (already reference),
  - frame-enabled coupled `t=3`, `t=6`, and if feasible `t=10` on 8 CPUs.
- Outputs:
  - regenerate Section-5.4-style figure using:
    - `tst/publication/plot_pic_parallel_shock_section54.py`
- Test gate:
  - shock remains interior longer than baseline.
  - pipeline remains stable (no rank crashes).
  - B-field metric trend can be tracked over longer times.

### Step 6: Branch parity and HPC readiness
- Add/adjust deck(s):
  - `inputs/tests/pic_parallel_shock_fine_uniform.athinput` with optional
    frame-tracking knobs commented for laptop vs HPC.
- Ensure no hard-coded local assumptions:
  - rank count and mesh scale remain runtime-configurable.
- Test gate:
  - 8-rank parity for `np=1/2/4/8` on short runs.
  - no code changes required to scale to larger domain/rank/GPU runs.

## Validation Matrix (Required)

### Compile and smoke (every step)
- `cmake --build /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/build -j8`
- `mpiexec -n 8 /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/build/src/athena ... time/tlim=1e-3 time/nlim=1`

### Functional checks
- Baseline equivalence with frame disabled.
- Frame enabled with very small `ps_frame_vfrac` to verify no blow-up.
- Coupled mode at `t=3` and beyond with figure regeneration.

### Physics checks
- Compare `x_shock(t)` drift rate vs baseline.
- Compare `B/B0` quantiles downstream over extended time.
- Compare spectrum proxy continuity (no artificial discontinuities at frame events).

## Risk Register and Mitigation
- Risk: reflect-wall physical inconsistency under frame shifts.
  - Mitigation: delay `ps_frame_t_start` until post-formation; cap `dV`; document caveat.
- Risk: AMR interactions.
  - Mitigation: uniform-only first pass with explicit guard.
- Risk: restart parity with dynamic state.
  - Mitigation: deterministic time-schedule `V_des(t)` avoids hidden mutable state.
- Risk: particle-fluid mismatch after shifts.
  - Mitigation: shift both particle and fluid x-velocities in same callback.

## Exit Criteria
- Numerical:
  - stable coupled runs on 8 CPUs through at least `t=6` without MPI/task failures.
- Kinematic:
  - modeled/diagnosed shock remains comfortably away from right boundary.
- Scientific:
  - longer-time runs provide a clear read on whether amplification develops.
  - if still absent, we escalate to a true moving-window remap design (Phase 2).

## Phase 2 (only if needed): true remap recentering
- If velocity-frame approach is insufficient, implement integer-cell recenter remap
  for uniform grids:
  - shift full MHD state left by `nshift` cells when analytic `x_shock_model`
    crosses threshold,
  - refill right-side cells with upstream state,
  - shift particle x positions by `-nshift*dx`,
  - keep this as a separate controlled feature flag.

---

## Execution Status (Updated)

### Completed
- [x] Step 0: Added generic `user_work_in_loop` callback plumbing in:
  - `src/pgen/pgen.hpp`
  - `src/pgen/pgen.cpp`
  - `src/driver/driver.cpp`
- [x] Step 1: Enrolled `pic_parallel_shock` work-in-loop callback and parsed frame knobs.
- [x] Step 2: Added analytic per-cycle velocity-frame update (`dv`) and fluid shift.
- [x] Step 3: Added PIC consistency for velocity-frame mode:
  - particle `IPVX` shift
  - inflow-state update option
  - injection velocity consistency
- [x] Step 4: Added frame guardrails and diagnostics:
  - uniform-grid guard (`no AMR/SMR` when required)
  - periodic frame diagnostics in logs
- [x] Phase-2 implementation path added:
  - deterministic recenter mode (`ps_frame_mode=recenter`)
  - ghost-zone-safe cell shifts (`ps_recenter_shift_cells <= nghost`)
  - particle recenter with rank-safe `CountParticles()` collective handling.

### Validation runs completed (`np=8`)
- Velocity frame mode:
  - `t=3`: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/moving_frame_validation/t3_final`
  - `t=6`: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/moving_frame_validation/t6_final`
  - `t=10`: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/moving_frame_validation/t10_final`
- Recenter mode:
  - `t=3.2`: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/moving_frame_recenter/t32`
  - `t=10`: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/moving_frame_recenter/t10`
  - `t=100`: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/moving_frame_recenter/t100`
- Figure generation:
  - `tst/publication/plot_pic_parallel_shock_section54.py` used for all runs.

### Current findings
- Velocity-only frame mode remains stable but does not keep the shock interior long
  enough for strong amplification in this setup.
- Recenter mode now runs stably to `t=100` on 8 ranks and preserves particle updates.
- In current local-resolution runs, magnetic amplification remains weak
  (`B/B0` high-percentile levels stay close to unity), so Section-5.4-like strong growth
  is not yet reproduced.

### Decomposition parity smoke checks
- Completed short coupled frame-enabled smokes at `np=1/2/4/8` without MPI failures.

## Post-Plan Corrective Pass (Detailed)

### Issue A: recenter validation rejected documented default
- Problem:
  - `ps_recenter_vshock_model=-1` (documented default) was rejected by a faulty
    validation condition.
- Root cause:
  - The condition treated all negative values as invalid, including `-1`.
- Fix:
  - Validation now accepts either:
    - exact default sentinel `-1` (within tolerance), or
    - explicit positive override (`>0`).
- File:
  - `src/pgen/tests/pic_parallel_shock.cpp`

### Issue B: analytic recenter model speed drifted at long times
- Problem:
  - With recenter enabled, long runs drifted because model speed used a strong-shock
    approximation (`0.5*(gamma-1)*u0`) that under-tracked finite-Mach runs in this deck.
- Root cause:
  - Local deck parameters are not in the asymptotic strong-shock limit.
- Fix:
  - Replaced shock-speed estimate with finite-Mach Rankine-Hugoniot form:
    - compute `M_s^2 = u0^2 / (gamma*p0/rho0)`
    - `r = ((gamma+1) M_s^2) / ((gamma-1) M_s^2 + 2)`
    - `v_shock = u0 / (r - 1)`
  - Added fatal guard for non-supersonic/invalid inferred speed.
- File:
  - `src/pgen/tests/pic_parallel_shock.cpp`

### Issue C: recenter robustness and visibility
- Problem:
  - For large `dev`, one-shot `nshift = dev*shift_cells` could exceed safe per-call
    shift assumptions; recenter diagnostics were also sparse.
- Fix:
  - Apply recenter as repeated per-event shifts (`dev` times) with fixed
    `ps_recenter_shift_cells`.
  - Emit recenter diagnostic lines on cadence and on event transitions, including
    `ev_now/ev_next/dev`.
- File:
  - `src/pgen/tests/pic_parallel_shock.cpp`

## New Validation Evidence (After Corrective Pass)

### Build gate
- `cmake --build /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/build -j8`
- Result: pass.

### Coupled recenter short gate (`np=8`, `t=10`)
- Run artifact:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/moving_frame_recenter/check_coupled_t10_v2`
- Key checks:
  - `shock_x ~ 2.164` (held near recenter trigger).
  - no fatal/MPI errors.
  - Section-5.4 plot generated.

### Coupled recenter long gate (`np=8`, `t=100`)
- Run artifact:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/moving_frame_recenter/check_coupled_t100_v2`
- Key checks:
  - completed to `t=100`, cycle `103258`.
  - recenter events continue through end of run (log shows `ev_now/ev_next` increments).
  - final `xshock_model ~ 2.169` from diagnostics; measured `shock_x ~ 2.148`.
  - Section-5.4 figure:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/moving_frame_recenter/check_coupled_t100_v2/fig/section54_check_coupled_t100_v2.png`

### Decomposition parity re-check (`np=1/2/4/8`)
- Command profile:
  - short coupled recenter smokes (`nlim=20`) on all listed rank counts.
- Result:
  - all pass without MPI errors.

## Current Scientific Readout (Local Laptop Scale)
- Moving-frame/recenter implementation is now functionally correct and stable for
  long runs (`t=100`) on 8 CPUs.
- For the current local deck (`u0=3`, `256x128`), magnetic amplification remains
  weak by `t=100`:
  - `B p95 ratio ~ 1.0058`, `p99 ~ 1.024`, `max ~ 1.046` (vs initial).
- Controlled stress probes (`eta` up to `10x`, `u0=10`) showed only modest change
  at this resolution/time window, indicating the remaining limitation is now mostly
  physics-signal scale (not frame-tracking mechanics).
