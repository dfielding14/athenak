> **HISTORICAL ONLY.** Do not execute commands or infer qualification,
> publication, branch, or Frontier-submission authority from this file. Use
> `tst/publication/PIC_PRODUCTION_READINESS_PLAN.md` as the sole controlling plan.

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

## Post-Plan Extension: Explicit Seed Noise Control

### Motivation
- Add a user-facing knob to seed larger initial transverse perturbations and
  accelerate instability onset studies.

### Implemented controls
- New `<problem>` parameters in `pic_parallel_shock`:
  - `ps_seed_noise_amp` (default `0.0`): relative amplitude of seeded transverse
    magnetic perturbations as a fraction of `B0`.
  - `ps_seed_noise_seed` (default `1234`): deterministic phase seed.

### Implementation details
- File:
  - `src/pgen/tests/pic_parallel_shock.cpp`
- The initial-state builder now adds deterministic multi-mode perturbations to
  transverse fields (`By`, `Bz`) as functions of `x` only.
- This preserves discrete divergence consistency in the 2D benchmark geometry
  (`Bx` remains uniform, `By/Bz = By/Bz(x)`).
- Perturbations are included before total-energy initialization, so magnetic
  energy bookkeeping remains consistent.

### Deck updates
- Added documented defaults to:
  - `inputs/tests/pic_parallel_shock_fine_uniform.athinput`
  - `inputs/tests/pic_parallel_shock_coarse_uniform.athinput`
  - `inputs/tests/pic_parallel_shock_amr_fiducial.athinput`

### Validation (8 CPUs)
- Build gate:
  - `cmake --build /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/build -j8`
- Smoke comparisons:
  - `ps_seed_noise_amp=0.0` vs `0.10` with `nlim=1`.
  - Output artifact root:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/seed_noise_validation`
  - Observed:
    - seeded run shows higher initial `|B|` variance and broader high-percentile tail.
- Coupled+recenter compatibility smoke:
  - `nlim=200`, `np=8`, coupled recenter + `ps_seed_noise_amp=0.10`.
  - Run artifact:
    - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/seed_noise_validation/coupled_recenter_noise`
  - Completed without MPI/runtime errors.

## Post-Plan Audit (No-Amplification Root Cause Pass)

### Finding 1: 2D PIC path was effectively 2D2V (major implementation gap)
- Root cause:
  - In `src/particles/particles_pushers.cpp`, the Boris path previously forced
    `vz=0`, `Bz=0`, `Uz=0`, and `cEz=0` whenever `nx3==1`.
  - `src/pgen/tests/pic_parallel_shock.cpp` also forced injected `dirz=0` in 2D.
- Impact:
  - Upstream deposited `Jz` channel was suppressed in 2D runs, removing a key
    transverse coupling path expected for Bell-like amplification behavior.
- Fix implemented:
  - Added runtime flag `<particles>/pic_enable_2d3v` (default `false`) in:
    - `src/particles/particles.hpp`
    - `src/particles/particles.cpp`
  - When enabled in 2D:
    - Boris path keeps `vz/Bz/Uz/cEz` active.
    - interpolation no longer hard-zeros `Bz/Uz`.
    - `pic_parallel_shock` injection samples nonzero `dirz`.
  - Enabled for parallel-shock benchmark decks:
    - `inputs/tests/pic_parallel_shock_fine_uniform.athinput`
    - `inputs/tests/pic_parallel_shock_coarse_uniform.athinput`
    - `inputs/tests/pic_parallel_shock_amr_fiducial.athinput`

### Validation of Finding 1 (`np=8`)
- A/B short gate (`t=1`, coupled, recenter):
  - off:
    `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/audit_2d3v/ab_t1/off`
  - on:
    `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/audit_2d3v/ab_t1/on`
  - metrics:
    `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/audit_2d3v/ab_t1/ab_metrics.json`
  - result:
    - `pic_enable_2d3v=false`: `Jz_up_rms = 0.0`
    - `pic_enable_2d3v=true`: `Jz_up_rms = 2.735e-06`

### Finding 2: even after 2D3V restoration, growth remains weak in local deck
- Long diagnostic run (`np=8`, coupled recenter, `t=20`, seed noise on):
  - run:
    `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/audit_2d3v/diag_t20`
  - metrics:
    `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/audit_2d3v/diag_t20/diag_currents_metrics.json`
- Comparison against prior no-2D3V audit:
  - old:
    `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/audit_no_amp/diag_t20/diag_currents_metrics.json`
  - `Jz` channel is now active, but `B` growth is still weak over local runtime
    windows.

### Finding 3: setup-physics mismatch vs Section 5.4 target is still large
- Section 5.4 reference (`athena++_MHD-PIC.pdf`, Sec. 5.4):
  - high Alfvén Mach (`M_A=30`), `eta=1e-3`, and a much larger 2D box.
- Current local benchmark deck intentionally trades physical scale for laptop cost:
  - `u0=3`, `B0=1` (`M_A~3`), compact domain, low runtime window.
- Consequence:
  - the local deck remains a correctness/pipe-validation setup, not yet a strong
    amplification reproduction setup.

## Critical Gap: Feedback Off/On Invariance (Essential)

### Current status
- This plan did not yet include a hard acceptance gate that requires measurable
  divergence between:
  - `couple_moments_momentum_to_mhd=false`, `couple_moments_energy_to_mhd=false`
  - `couple_moments_momentum_to_mhd=true`,  `couple_moments_energy_to_mhd=true`
- New explicit A/B checks (same seed, same mesh, `np=8`) still show near machine
  precision off/on deltas for local stress levels:
  - `t=10` A/B metrics:
    `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/audit_feedback_ab_2d3v/feedback_ab_t10_metrics.json`
  - boosted short A/B (`deposit_qscale=1e-4`, `t=1`) metrics:
    `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/audit_feedback_ab_2d3v/boosted/feedback_ab_t1_qscale1e-4_metrics.json`
- Diagnostic source-channel outputs in a coupled run are currently zero at dump
  times:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/audit_feedback_ab_2d3v/diag_fb_fields_metrics.json`
  - `prtcl_dpxdt`, `prtcl_dpydt`, `prtcl_dpzdt`, `prtcl_dedt`, `prtcl_ebdot` all
    report zero in the sampled outputs.

### Interpretation
- The essential requirement ("feedback on/off must produce a detectable
  difference") is not satisfied yet in this benchmark path.
- With current evidence, this is either:
  - a scheduling/plumbing issue in when/where feedback moments are consumed, or
  - a diagnostics-timing issue that hides nonzero sources at output time, while
    source impact remains too weak in current local parameters.

### Required follow-up gates (must pass before physics claims)
1. Add an automated `feedback_sensitivity_ab` regression:
   - fixed seed; identical run pairs; only feedback toggles differ.
   - enforce nonzero threshold on at least one fluid delta metric
     (for example, downstream `B` or momentum norm change above tolerance).
2. Add cycle-local instrumentation around feedback application:
   - pre/post `MHD::MHDSrcTerms` (or `EFieldSrc`, per order) norms for
     `IM1/IM2/IM3/IEN`.
   - per-cycle norms of deposited `IMOM_DPXDT/DPYDT/DPZDT/DEDT`.
3. Confirm output observability:
   - write diagnostics at a point guaranteed after deposition and before zeroing,
     or emit explicit history scalars.
4. Only after (1)-(3), rerun long moving-frame benchmarks and reassess magnetic
   amplification and off/on divergence trends.

## Post-Plan Continuation (Latest)

### Fix D: `boris_tsc` 2D interpolation bug (hard blocker)
- Problem:
  - In 2D runs, `InterpolateTSC()` used `wz={0,1,0}` while iterating only one
    z-sample (`nk=1`), so effective z-weight was always zero.
  - This forced sampled particle fields to near-zero in 2D and collapsed
    per-particle delta channels (`IPDP*`, `IPDE`) to zero even when coupled mode
    was enabled.
- Fix:
  - In `src/particles/particles_pushers.cpp`, set 2D TSC weights to
    `wz={1,0,0}` for the collapsed-z path.
  - Kept 3D branch unchanged.
- Validation (`np=8`):
  - Run: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/diag_2d3v/on/run.log`
  - `pr_b_rms` and `pr_dpdt_rms` become nonzero after injection starts.
  - Deposited `IMOM_DPXDT/DPYDT/DPZDT/DEDT` channels become nonzero.

### Fix E: recenter target was parsed but unused
- Problem:
  - `ps_recenter_x_target` was never used in the event-count law, so shock stayed
    near trigger (`~2.9-3.0`) instead of target (`2.0`).
- Fix:
  - In `src/pgen/tests/pic_parallel_shock.cpp`, recenter event count now uses
    target-based overshoot:
    - `events = ceil((x_unshifted - x_target)/dx_shift)` when `x_unshifted > trigger`.
  - Added domain validation for `ps_recenter_x_target` and trigger upper bound.
- Validation (`np=8`, `t=5`):
  - Run: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/recenter_target_check`
  - Measured shock position stabilizes near `x~1.89-1.95` after formation.

### Fix F: delta-feedback density normalization
- Problem:
  - Boris delta channels were deposited as particle-integrated rates, but consumed
    by MHD as conservative density source rates.
- Fix:
  - In `src/particles/particles_moments.cpp`, normalize deposited
    `IMOM_DPXDT/DPYDT/DPZDT/DEDT` by local cell volume.
- Validation (`np=8`):
  - Run: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/diag_after_volnorm/run.log`
  - Source-channel magnitudes increase from near-zero to finite physically-active
    levels (`step_mom_l1`, `step_eng_l1` nonzero and sustained).

### Updated A/B sensitivity status
- Short-moderate coupled A/B remains clearly nonzero:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/postvol_ab/metrics_t2.json`
  - `b_l2 ~ 2.53e-02`, `rho_l2 ~ 4.63e-02`.
- High-drive target-centered A/B at `t=5` still shows weak downstream-mean
  divergence (small relative delta), indicating coupling is active but not yet in
  a strong-growth regime at local laptop scale.

### Long-run outcome after fixes
- Corrected moderate run to `t=100` (`np=8`) remains stable:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/postfix_long/t100_on`
  - Trend: `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/postfix_long/t100_on/bmag_trend.json`
- Section-5.4 style figures regenerated from corrected runs:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/postfix_long/t100_on/pic_parallel_shock_section54_t100.png`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/scan_growth/case4/pic_parallel_shock_section54_case4_t20.png`

### Current conclusion
- The major implementation blockers are resolved:
  1. 2D TSC field sampling bug fixed.
  2. recenter target logic fixed.
  3. delta feedback normalization fixed.
- The pipeline is now physically and numerically consistent enough to scale to
  larger domains/ranks. On this laptop-scale grid/time budget, amplification is
  detectable but still modest; stronger Section-5.4-like growth will require HPC
  scale-up (larger upstream extent, longer runtime, and finer transverse
  resolution) using the same corrected code path.
