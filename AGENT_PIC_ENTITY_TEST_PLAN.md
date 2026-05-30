> **HISTORICAL ONLY.** Do not execute commands or infer qualification,
> publication, branch, or Frontier-submission authority from this file. Use
> `tst/publication/PIC_PRODUCTION_READINESS_PLAN.md` as the sole controlling plan.

# AthenaK MHD-PIC Entity-Mirroring Test Implementation Plan

This document defines the full sequential plan to implement an Entity-mirroring
PIC test suite in AthenaK, including the required runtime isolation modes so
PIC behavior can be validated without significant MHD back-reaction.

## Execution Status

- Step 0 (`Baseline Lock and Test Harness Setup`): complete at commit
  `e997aa2a` with baseline evidence in
  `tst/.codex/pic_entity_suite_baseline/STEP0_BASELINE_REPORT.md`.
- Step 1 (`Add PIC Runtime Isolation Mode Knob`): complete at current working
  tree state with evidence in
  `tst/.codex/pic_entity_suite/step1_runtime_knobs/STEP1_RUNTIME_KNOBS_REPORT.md`.
- Step 2 (`Implement passive_mhd mode`): complete at current working tree state
  with evidence in
  `tst/.codex/pic_entity_suite/step2_passive_mhd/STEP2_PASSIVE_MHD_REPORT.md`.
- Step 3 (`Implement no_mhd mode`): complete at current working tree state with
  evidence in `tst/.codex/pic_entity_suite/step3_no_mhd/STEP3_NO_MHD_REPORT.md`.
- Step 4 (`Upgrade Particle Pusher to Section 2 E+B Midpoint-Boris`): complete
  at current working tree state with evidence in
  `tst/.codex/pic_entity_suite/step4_midpoint_eb/STEP4_MIDPOINT_EB_REPORT.md`.
- Step 5 (`Add PIC Diagnostics Needed for Analytic Validation`): complete at
  current working tree state with evidence in
  `tst/.codex/pic_entity_suite/step5_diagnostics/STEP5_DIAGNOSTICS_REPORT.md`.
- Step 6 (`Finalize Deposit Parity Layer`): complete at current working tree
  state with evidence in
  `tst/.codex/pic_entity_suite/step6_entity_parity/STEP6_ENTITY_PARITY_REPORT.md`.
- Step 7 (`Implement EM Vacuum Wave Analytic Test`): complete at current
  working tree state with evidence in
  `tst/.codex/pic_entity_suite/step7_em_vacuum/STEP7_EM_VACUUM_REPORT.md`.
- Step 8 (`Implement Langmuir Test (Frequency Accuracy)`): complete at current
  working tree state with evidence in
  `tst/.codex/pic_entity_suite/step8_langmuir_proxy/STEP8_LANGMUIR_PROXY_REPORT.md`.
- Step 9 (`Implement Two-Stream Instability Test (Growth Rate)`): complete at
  current working tree state with evidence in
  `tst/.codex/pic_entity_suite/step9_two_stream/STEP9_TWO_STREAM_REPORT.md`.
- Step 10 (`Implement Weibel Instability Test (Growth Rate)`): complete at
  current working tree state with evidence in
  `tst/.codex/pic_entity_suite/step10_weibel/STEP10_WEIBEL_REPORT.md`.
- Step 11 (`Implement Bell Test (Current-Driven Instability)`): complete at
  current working tree state with evidence in
  `tst/.codex/pic_entity_suite/step11_bell/STEP11_BELL_REPORT.md`.
- Step 12 (`Implement Multi-Species Backreaction Oscillation Test`): complete
  at current working tree state with evidence in
  `tst/.codex/pic_entity_suite/step12_multispecies_osc/STEP12_MULTISPECIES_OSC_REPORT.md`.
- Step 13 (`Implement CRSI Linear-Dispersion Test with delta f`): complete at
  current working tree state with evidence in
  `tst/.codex/pic_entity_suite/step13_crsi_deltaf_proxy/STEP13_CRSI_DELTAF_PROXY_REPORT.md`.
- Step 14 (`Implement CRPAI Linear-Dispersion Test`): complete at current
  working tree state with evidence in
  `tst/.codex/pic_entity_suite/step14_crpai_proxy/STEP14_CRPAI_PROXY_REPORT.md`.
- Step 15 (`Implement Driven CRPAI in Expanding/Compressing Box`): complete at
  current working tree state with evidence in
  `tst/.codex/pic_entity_suite/step15_expanding_box_proxy/STEP15_EXPANDING_BOX_PROXY_REPORT.md`.
- Step 16 (`Add Refinement-Boundary Characterization Tests`): complete at
  current working tree state with evidence in
  `tst/.codex/pic_entity_suite/step16_refinement_boundary/STEP16_REFINEMENT_BOUNDARY_REPORT.md`.
- Step 17 (`Add Reduced-Size AMR Shock + Load-Balancing Smoke Test`): complete
  at current working tree state with evidence in
  `tst/.codex/pic_entity_suite/step17_amr_shock_lb/STEP17_AMR_SHOCK_LB_REPORT.md`.
- Step 18 (`Add Full MPI Matrix and Decomposition Robustness Gates`): complete
  at current working tree state with evidence in
  `tst/.codex/pic_entity_suite/step18_mpi_matrix/STEP18_MPI_MATRIX_REPORT.md`.
- Step 19 (`Restart and Safety Guard Coverage`): complete at current working
  tree state with evidence in
  `tst/.codex/pic_entity_suite/step19_restart_safety/STEP19_RESTART_SAFETY_REPORT.md`.
- Step 20 (`Documentation and Integration Closeout`): complete at current
  working tree state (this document, root/subdirectory `AGENTS.md`, and
  PIC handoff/implementation guides updated to reflect as-built behavior).

## 1. Scope and Goal

Goal:
- Add a production-grade PIC regression suite that mirrors Entity test intent as
  closely as AthenaK architecture allows.
- Cover both PIC-only behavior and MHD-PIC interaction behavior.
- Keep all changes gated by deterministic pass/fail criteria in AthenaK's test
  harness.

Primary references:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_IMPLEMENTATION_GUIDE.md`
- `/Users/dbf75/Work/Research/AthenaK/entity`
- `/Users/dbf75/Work/Research/AthenaK/athena++_MHD-PIC.pdf`

## 2. Current Baseline Constraints (Must Be Addressed)

1. CR Boris pusher currently hard-requires MHD fields.
- Evidence:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_pushers.cpp:195`

2. CR Boris path is currently B-only (`E=0`), which is insufficient for
   Langmuir/two-stream/Weibel/Bell physics validation.
- Evidence:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_pushers.cpp:234`

3. Particle CR payload currently stores sampled B components but no sampled E
   components.
- Evidence:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/athena.hpp:62`

4. Existing Entity-mirroring deposit tests are already present and should be
   retained as the first parity layer.
- Evidence:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/inputs/tests/pic_entity_deposit_mink.athinput`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/inputs/tests/pic_entity_deposit_reflect.athinput`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_entity_deposit_mink.py`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_entity_deposit_reflect.py`

### 2.1 Findings from `athena++_MHD-PIC.pdf` Section 2

Reviewed source:
- `/Users/dbf75/Work/Research/AthenaK/athena++_MHD-PIC.pdf`
- Section 2 is on rendered pages 2-4 in:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/pdfs/athena_mhd_pic/page-02.png`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/pdfs/athena_mhd_pic/page-03.png`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/pdfs/athena_mhd_pic/page-04.png`

Key findings to adopt:
1. CR motion uses full Lorentz-force form:
   - `dx/dt = v`
   - `d(p/m)/dt = (q/mc) * (cE + v x B)`
2. Gas feedback is formulated via CR charge/current and momentum/energy source
   terms; not just current-to-electric-field coupling.
3. Electric field for particle push is computed from ideal frozen-in relation:
   - `cE = -u x B`
4. The integration algorithm is a strict two-stage predictor/corrector:
   - Stage 1 predicts midpoint states.
   - Stage 2 performs full update using midpoint-interpolated fields.
5. Midpoint strategy is explicit:
   - move particles by `dt/2` to midpoint,
   - interpolate `u` and `B` at midpoint,
   - compute `E` at midpoint,
   - perform Boris momentum push over full `dt`,
   - move particles by second `dt/2`.
6. Feedback deposit for conservation uses particle momentum/energy changes
   (`delta p`, `delta epsilon`) in stage 2, rather than only raw moments.
7. Interpolation/deposition uses TSC in the basic method.
8. Time step constraints include:
   - MHD CFL,
   - max particle cell crossing per step (paper baseline: <=2),
   - max gyro-rotation angle per step (paper recommendation: <=0.3 rad).
9. When backreaction is negligible, feedback terms can be skipped and CRs can
   be advanced as test particles.

### 2.2 Adopted Design Requirements from PDF Review

The plan below is updated to enforce the Section 2 approach:
1. E+B pusher upgrade must be midpoint-based with a full Lorentz-force Boris
   advance.
2. E-field used by the pusher must be derived from midpoint `u` and `B` using
   the frozen-in condition in the Section 2 baseline path.
3. Coupled backreaction path must support depositing `delta p` and
   `delta epsilon` to gas.
4. Test-particle mode must be explicit for negligible-backreaction runs.
5. TSC interpolation/deposition must remain the default parity path for this
   upgrade.
6. Runtime time-step guards must include max cell-crossing and max gyro-angle
   controls in addition to CFL.

### 2.3 Additional Findings from Full-Paper Review (Sections 3-6)

Reviewed source:
- `/Users/dbf75/Work/Research/AthenaK/athena++_MHD-PIC.pdf`
- Section 3 (expanding/compressing box) on rendered pages:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/pdfs/athena_mhd_pic/page-04.png`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/pdfs/athena_mhd_pic/page-05.png`
- Section 4 (optimization/refinement/`delta f`) on rendered pages:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/pdfs/athena_mhd_pic/page-05.png`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/pdfs/athena_mhd_pic/page-08.png`
- Section 5 benchmark ladder on rendered pages:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/pdfs/athena_mhd_pic/page-08.png`
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tmp/pdfs/athena_mhd_pic/page-16.png`

Key additional findings to adopt:
1. Expanding/compressing-box support is a first-class part of the method with
   explicit source terms in both MHD and particle equations.
2. In expanding coordinates, the particle momentum equation adds `D · (p/m)`
   terms; implementation applies these as half-steps around the Boris push.
3. Performance path relies on intermediate arrays for interpolation/deposition
   locality plus periodic particle sorting, especially at moderate/high ppc.
4. Parallel path relies on non-blocking MPI particle communication in TaskList
   and load balancing that includes particle integration cost.
5. At fine/coarse refinement boundaries, deposition is designed for smoothness
   and homogeneity; strict per-particle feedback conservation can be relaxed.
6. `delta f` is essential for low-anisotropy CRSI/CRPAI regimes: noise is
   reduced strongly, but exact machine-precision conservation is not expected
   with variable particle weights.
7. The benchmark ladder is not only Bell/gyration; it also includes
   multi-species backreaction oscillation, CRSI, CRPAI, and expanding-box CRPAI.

### 2.4 Additional Adopted Requirements from Full-Paper Review

The sequential plan below is further extended with these requirements:
1. Add expanding/compressing-box runtime controls for driven CRPAI workflows.
2. Add an explicit `delta f` runtime path and diagnostics workflow for
   low-anisotropy instability tests.
3. Add refinement-boundary homogeneity tests and explicit conservation
   characterization near coarse/fine interfaces.
4. Add CRSI/CRPAI tests with polarization-resolved spectra and growth-rate fits.
5. Add multi-species CR-backreaction oscillation tests with uniform/SMR/AMR
   parity checks.
6. Add at least one reduced-size AMR shock/load-balancing smoke scenario in
   the AthenaK harness for integration/performance continuity.

## 3. Sequential Execution Plan (Do In Order)

### 3.1 Incremental Execution and Checkpoint Commit Protocol (Required)

All steps in this plan are executed incrementally with hard checkpoints.
Do not start a new step until the current step is fully validated, reviewed,
and committed.

Required loop for every step:
1. Implement the smallest coherent change set for that step only.
2. Run targeted tests for that step plus required style checks for changed files.
3. Perform a thorough review pass before commit:
   - correctness against the step objective and gates,
   - safety/guard behavior (fatal paths and unsupported-mode handling),
   - regression risk against prior passing tests,
   - decomposition/MPI behavior where relevant,
   - diff sanity (no unrelated edits, line-wrap/style compliance).
4. Record validation evidence and key metrics in
   `tst/.codex/pic_entity_suite/<step_name>/`.
5. Stage only the files relevant to that step.
6. Commit with a step-scoped message (for example:
   `PIC plan step N: <short description>`).
7. Record commit SHA in step notes/logs so rollback points are explicit.

Checkpoint rules:
1. No mixed-step commits.
2. If review finds issues, fix and re-run validation before staging.
3. Keep each commit bisectable (buildable and testable for that checkpoint).

### Step 0: Baseline Lock and Test Harness Setup

Objective:
- Freeze current behavior and establish a reproducible baseline before feature
  changes.

Actions:
1. Record HEAD, branch, dirty state, and baseline test outcomes.
2. Create a dedicated log artifact directory:
   - `tst/.codex/pic_entity_suite_baseline/`
3. Re-run current passing PIC tests (serial + MPI) and style gates.

Completion gate:
- Baseline report exists with command lines, pass/fail, and key metrics.

### Step 1: Add PIC Runtime Isolation Mode Knob

Objective:
- Introduce explicit runtime mode selection for PIC coupling context.

Actions:
1. Add runtime option in `<particles>` block:
   - `pic_background_mode = coupled | passive_mhd | no_mhd`
2. Parse and validate in:
   - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp`
   - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp`
3. Preserve current behavior as default:
   - `pic_background_mode = coupled`
4. Add Section 2 control knobs (with conservative defaults):
   - `pic_feedback_mode = coupled | test_particle`
   - `pic_interp_scheme = tsc` (default)
   - `pic_cr_light_speed = <Real>` (artificial particle light speed, C)
   - `pic_max_cell_cross = 2`
   - `pic_theta_max = 0.3`
   - `pic_deltaf_mode = off | on`
   - `pic_sort_interval = <int>` (0 disables runtime re-sorting)
   - `pic_intermediate_arrays = auto | off`
   - `pic_expanding_box_mode = off | on`
   - `pic_expansion_rate_x1/x2/x3 = <Real>`
5. Add validation guardrails for incompatible combinations:
   - `pic_deltaf_mode=on` requires a defined background distribution `f0`
   - expansion controls must be rejected unless expanding-box mode is enabled

Completion gate:
- Existing tests pass unchanged under default mode.
- Invalid mode values fail with clear fatal diagnostics.

### Step 2: Implement `passive_mhd` Mode (First Isolation Path)

Objective:
- Keep MHD containers/field carriers available to PIC while freezing fluid
  evolution so MHD impact is negligible.

Actions:
1. Add runtime branch in MHD task path to freeze fluid updates under
   `passive_mhd`.
2. Preserve only required field communication/coupling operations needed by PIC
   test workflows.
3. Default `passive_mhd` to Section 2-style negligible-backreaction operation
   (`pic_feedback_mode=test_particle`) unless coupling is explicitly requested.
4. Add guardrails preventing unsupported combinations with this mode.

Primary files:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mhd/mhd_tasks.cpp`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp`

Completion gate:
- In dedicated tests, MHD conserved-variable drift stays below defined tolerance
  while PIC deposition/current behavior remains active.

### Step 3: Implement `no_mhd` Mode (True PIC-Only Path)

Objective:
- Enable PIC runs without requiring `<mhd>` module allocation.

Actions:
1. Abstract particle field access away from direct `pmhd` dependency.
2. Remove hard failure in Boris path for `no_mhd` mode.
3. Provide a particle-owned/solver-owned EM field container path for PIC-only
   tests.
4. Preserve the same Section 2 midpoint particle update ordering in `no_mhd`
   mode, with gas feedback operations disabled.

Primary files:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_pushers.cpp`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.hpp`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/mesh/meshblock_pack.cpp`

Completion gate:
- AthenaK can run PIC test decks with `<particles>` and no `<mhd>` block when
  `pic_background_mode=no_mhd`.

### Step 4: Upgrade Particle Pusher to Section 2 E+B Midpoint-Boris

Objective:
- Make PIC instability/wave tests physically meaningful by adopting the
  Section 2 midpoint E+B Boris integration sequence.

Actions:
1. Extend CR particle payload to include sampled E components.
2. Implement explicit midpoint workflow in the pusher:
   - first `dt/2` position advance to midpoint,
   - midpoint interpolation of `u` and `B`,
   - midpoint `E` from frozen-in condition (`cE = -u x B`),
   - full-`dt` Boris momentum update,
   - second `dt/2` position advance to final state.
3. Record per-particle `delta p` and `delta epsilon` needed for Stage 2
   conservation-style gas feedback deposition.
4. Add/update deposit wrappers to consume `delta p`/`delta epsilon` in coupled
   mode.
5. Keep deterministic fallback behavior for existing tests.

Primary files:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/athena.hpp`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles_pushers.cpp`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/particles.cpp`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/particles/AGENTS.md`

Completion gate:
- Dedicated single-particle E/B verification tests pass against analytic
  trajectories or known limits.
- Midpoint electric field orthogonality check (`E dot B ~= 0`) passes for
  frozen-in path.
- Coupled mode conservation checks using deposited `delta p`/`delta epsilon`
  pass to tolerance.

### Step 5: Add PIC Diagnostics Needed for Analytic Validation

Objective:
- Ensure test scripts can measure physically relevant quantities robustly.

Actions:
1. Add/verify outputs for `E`, `B`, `rho`, `J`, and any required mode-analysis
   fields.
2. Add robust helper routines for FFT/mode amplitude/time-series extraction in
   regression scripts.
3. Ensure outputs are deterministic and decomposition-safe.
4. Add diagnostics for Section 2 pusher/backreaction validation:
   - midpoint sample checks (`E dot B`, interpolation sanity),
   - aggregate `delta p`/`delta epsilon` bookkeeping checks.
5. Add full-paper diagnostics required by instability/expansion tests:
   - polarization-resolved spectra (`forward/backward`, `left/right`),
   - mode-by-mode linear growth fits and fit-window quality scores,
   - anisotropy metrics (`xi`, `f(p,mu)/<f>_mu`) for CRPAI workflows,
   - particle phase and energy-drift metrics for gyration checks.

Primary files:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/outputs/outputs.hpp`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/outputs/basetype_output.cpp`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/src/outputs/derived_variables.cpp`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/*.py`

Completion gate:
- Test scripts can compute all target scalar metrics from outputs without
  manual post-processing.
- Section 2-specific diagnostics are available for both `coupled` and
  `test_particle` modes.
- Full-paper benchmark diagnostics are available and reproducible in serial/MPI.

### Step 6: Finalize Deposit Parity Layer (Entity-Mink + Reflect)

Objective:
- Keep existing Entity-style deposit tests as foundational parity checks.

Actions:
1. Re-baseline and tighten pass/fail for:
   - `pic_entity_deposit_mink`
   - `pic_entity_deposit_reflect`
2. Add explicit decomposition-parity checks by field, not only global totals.
3. Add coarse/fine boundary homogeneity checks for uniform particle
   distributions in SMR/AMR layouts.
4. Add explicit conservation-characterization metrics near refinement
   boundaries, with documented expected tolerances.

Completion gate:
- Serial/MPI decomposition parity passes with deterministic thresholds.
- TSC parity path remains stable under the new Section 2 midpoint pusher.
- Refinement-boundary smoothness and tolerance-bounded conservation behavior are
  documented and passing.

### Step 7: Implement EM Vacuum Wave Analytic Test

Objective:
- Add a clean EM-field propagation benchmark with analytic convergence checks.

Actions:
1. Add deck and script mirroring Entity EM-vacuum intent.
2. Compute L1/L2 error norms against analytic wave solution at multiple
   resolutions.
3. Add convergence-rate gate.

Completion gate:
- Error magnitude and convergence thresholds pass in serial and MPI.

### Step 8: Implement Langmuir Test (Frequency Accuracy)

Objective:
- Validate electrostatic oscillation frequency behavior.

Actions:
1. Add Langmuir deck + analysis script.
2. Extract time series of E-field mode amplitude.
3. Fit dominant frequency and compare with expected plasma frequency envelope.

Completion gate:
- Frequency error is within fixed tolerance across serial + MPI.

### Step 9: Implement Two-Stream Instability Test (Growth Rate)

Objective:
- Validate unstable electrostatic growth behavior.

Actions:
1. Add two-stream deck + analysis script.
2. Compute mode amplitude growth in linear phase.
3. Fit exponential growth rate and compare to accepted threshold band.

Completion gate:
- Fitted growth rate falls within tolerance band over defined fit window.

### Step 10: Implement Weibel Instability Test (Growth Rate)

Objective:
- Validate transverse electromagnetic instability growth.

Actions:
1. Add Weibel deck + analysis script.
2. Track magnetic mode energy growth.
3. Fit linear-phase exponential growth rate.

Completion gate:
- Growth rate and phase-window quality checks pass.

### Step 11: Implement Bell Test (Current-Driven Instability)

Objective:
- Validate current-driven unstable mode behavior.

Actions:
1. Add Bell deck + script.
2. Track selected unstable mode metrics and growth envelope.
3. Add 1D baseline with optional 2D/3D parity checks for dispersion trends.
4. Add suppression-relation checks where applicable to linear-phase metrics.
5. Validate against predetermined threshold band and decomposition stability.

Completion gate:
- Growth behavior passes in serial and MPI, with stable decomposition parity.

### Step 12: Implement Multi-Species Backreaction Oscillation Test (Paper 5.3)

Objective:
- Validate CR backreaction with multiple species and mesh-refinement parity.

Actions:
1. Add an electron/positron oscillation deck and harness script mirroring the
   paper's uniform-field setup.
2. Validate oscillation frequency and conserved-energy surrogate against the
   analytic target in the linear regime.
3. Run uniform-grid, SMR, and AMR variants and compare all primary metrics.

Completion gate:
- Frequency/phase/energy metrics pass tolerance.
- Uniform/SMR/AMR runs agree within deterministic parity limits.

### Step 13: Implement CRSI Linear-Dispersion Test with `delta f` (Paper 5.5.1)

Objective:
- Validate low-anisotropy CR streaming instability growth with reduced noise.

Actions:
1. Add isotropic `kappa` CR setup and `delta f` configuration path.
2. Extract polarization-resolved wave spectra and fit linear growth by `k`.
3. Compare measured growth against analytical dispersion bands.

Completion gate:
- Growth-rate fits match expected bands over a documented fit window.
- `delta f` path is stable in serial/MPI and clearly superior to full-`f` noise.

### Step 14: Implement CRPAI Linear-Dispersion Test (Paper 5.5.2)

Objective:
- Validate anisotropy-driven instability and polarization-branch selectivity.

Actions:
1. Add prolate/oblate initialization variants (`xi < 1`, `xi > 1`).
2. Measure polarization-resolved spectra and growth rates in the linear phase.
3. Compare mode growth against analytical expectations and branch selection.

Completion gate:
- Dominant branch selection and growth-rate fits are correct for both
  anisotropy flavors.

### Step 15: Implement Driven CRPAI in Expanding/Compressing Box (Paper 5.5.3)

Objective:
- Validate continuously driven anisotropy workflows in expanding coordinates.

Actions:
1. Implement test setups for expanding and compressing runs with controlled
   anisotropic scale factors.
2. Use an isotropic initial distribution plus evolving `delta f` background as
   specified by the expanding-coordinate formulation.
3. Verify expected anisotropy trend and wave-branch behavior over time.

Completion gate:
- Expanding/compressing runs reproduce expected qualitative/quantitative trends
  in anisotropy and polarized wave growth.

### Step 16: Add Refinement-Boundary Characterization Tests (Paper 4.3)

Objective:
- Verify intended SMR/AMR deposition behavior near coarse/fine interfaces.

Actions:
1. Build dedicated tests for uniform particle distributions crossing
   refinement boundaries.
2. Validate homogeneity/smoothness of deposited moments across boundaries.
3. Record and gate bounded momentum/energy mismatch where non-unit deposition
   weights are expected by design.

Completion gate:
- Boundary smoothness checks pass and conservation-mismatch metrics remain
  within documented thresholds.

### Step 17: Add Reduced-Size AMR Shock + Load-Balancing Smoke Test

Objective:
- Keep paper-level AMR/load-balancing behavior under regression coverage.

Actions:
1. Add a reduced-size non-relativistic parallel-shock scenario compatible with
   the AthenaK harness runtime budget.
2. Validate basic physical smoke metrics (field amplification trend, non-thermal
   tail emergence trend).
3. Capture load-balancing telemetry and verify non-degenerate workload
   redistribution.

Completion gate:
- Physics smoke metrics and load-balancing telemetry pass documented thresholds.

### Step 18: Add Full MPI Matrix and Decomposition Robustness Gates

Objective:
- Ensure all new tests are robust across process layouts.

Actions:
1. Run each new test at minimum `np=1,2,4`.
2. For selected tests, include alternate meshblock decompositions.
3. Compare scalar metrics and field-level max-diff across decompositions.

Completion gate:
- All new tests pass decomposition parity criteria with no unexplained drift.

### Step 19: Restart and Safety Guard Coverage

Objective:
- Ensure no unsafe mode interactions and restart behavior remains explicit.

Actions:
1. Add restart A/B checks for representative new PIC tests.
2. Retain/add negative tests for unsupported runtime combinations.
3. Preserve explicit watch item for per-rank restart restore behavior/messages in
   coupled PIC runs.

Completion gate:
- Restart fidelity tests pass for supported paths.
- Unsupported paths fail with correct, deterministic reason strings.

### Step 20: Documentation and Integration Closeout

Objective:
- Integrate new suite into project docs and test inventory.

Actions:
1. Update:
   - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_HANDOFF.md`
   - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/AGENT_PIC_IMPLEMENTATION_GUIDE.md`
   - subdirectory `AGENTS.md` files touched by implementation
2. Ensure each test has clear purpose, metrics, and thresholds documented.
3. Add final status matrix: implemented/validated/deferred.

Completion gate:
- Documentation reflects as-built behavior and exact test commands.

### Final Status Matrix

| Step | Title | Status | Validation State |
| --- | --- | --- | --- |
| 0 | Baseline lock and harness setup | implemented | validated |
| 1 | Runtime isolation knob | implemented | validated |
| 2 | `passive_mhd` mode | implemented | validated |
| 3 | `no_mhd` mode | implemented | validated |
| 4 | Section-2 E+B midpoint-Boris push | implemented | validated |
| 5 | PIC diagnostics for analytic validation | implemented | validated |
| 6 | Deposit parity layer finalize | implemented | validated |
| 7 | EM vacuum wave anchor | implemented | validated |
| 8 | Langmuir frequency proxy | implemented | validated |
| 9 | Two-stream proxy | implemented | validated |
| 10 | Weibel proxy | implemented | validated |
| 11 | Bell proxy | implemented | validated |
| 12 | Multi-species oscillation parity | implemented | validated |
| 13 | CRSI `delta f` proxy | implemented | validated |
| 14 | CRPAI polarization proxy | implemented | validated |
| 15 | Expanding/compressing-box proxy | implemented | validated |
| 16 | Refinement-boundary characterization | implemented | validated |
| 17 | AMR shock/load-balance smoke | implemented | validated |
| 18 | MPI/decomposition robustness matrix | implemented | validated |
| 19 | Restart/safety guard coverage | implemented | validated |
| 20 | Documentation/integration closeout | implemented | validated |

Deferred:
- none in this execution plan scope.

## 4. Validation Command Matrix (Required)

Run all new tests via AthenaK harness with MPI enabled build when required.

Example pattern:
- `cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst`
- `python run_tests.py particles/<test_name> --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx`

Style gates (changed files):
- `cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/style && bash check_athena_cpp_style_changed.sh`
- `cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/style && bash check_python_style_changed.sh`

Artifact logging convention:
- `tst/.codex/pic_entity_suite/<step_name>/`
- Store per-command logs and extracted numeric summaries.

## 5. Merge Gates for This Plan

No step is considered complete unless all relevant gates are green.

Required final gate before declaring complete:
1. Existing PIC-MHD regressions remain green.
2. New Entity-mirroring PIC suite is green in serial + MPI.
3. CRSI/CRPAI (including driven expanding-box CRPAI) validations are green with
   documented fit windows and thresholds.
4. Refinement-boundary smoothness/conservation-characterization gates are green.
5. Deterministic analytic metrics and thresholds are documented in each test.
6. No blocking parity divergence remains unclassified.
7. Watch item retained and documented:
   - per-rank restart restore behavior/messages for coupled PIC runs.
8. Every completed step has its own reviewed, validated, staged, and committed
   checkpoint with traceable evidence logs.

## 6. Out-of-Scope for This Plan

1. Full architectural rework of AthenaK task scheduler.
2. Non-PIC physics feature additions unrelated to PIC testability.
3. Large refactors not required to enable or validate the test suite.

## 7. Immediate Next Action

Plan complete. Any further work is a follow-on scope (for example, additional
analytic PIC anchors beyond this 20-step ladder).
