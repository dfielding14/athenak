# PIC Entity Plan Step 2 Report: Passive-MHD Isolation Mode

## Scope
- Implemented passive background-mode behavior to keep MHD state containers
  available for PIC while freezing fluid evolution.
- Added runtime guardrails for unsupported passive-mode coupling combinations.
- Added a dedicated regression (`particles/pic_mhd_passive_mode`) validating:
  - particle deposition/current activity remains correct,
  - MHD state drift remains negligible (frozen),
  - guard failures are deterministic for invalid passive-mode combinations.

## Code Changes
- `src/particles/particles.cpp`
  - `pic_feedback_mode` default now depends on background mode:
    - `coupled` for `pic_background_mode=coupled`
    - `test_particle` for `pic_background_mode=passive_mhd`
  - Added guard: `pic_feedback_mode=test_particle` rejects particle-to-MHD
    coupling toggles (`couple_moments_to_mhd`, momentum/energy feedback).
  - Added guard: `pic_background_mode=passive_mhd` requires active `<mhd>`.
  - Added guard: `pic_background_mode=passive_mhd` requires
    `pic_feedback_mode=test_particle` in current implementation stage.
- `src/mhd/mhd_update.cpp`
  - `MHD::RKUpdate` now short-circuits in passive mode.
- `src/mhd/mhd_tasks.cpp`
  - `MHD::MHDSrcTerms` now short-circuits in passive mode.
- `src/mhd/mhd_ct.cpp`
  - `MHD::CT` now short-circuits in passive mode.

## New/Updated Test Inputs and Scripts
- Added `inputs/tests/pic_mhd_passive_mode.athinput`.
- Updated `inputs/tests/pic_mhd_current_coupling.athinput` with staged
  `pic_*` defaults for command-line override coverage.
- Added `tst/scripts/particles/pic_mhd_passive_mode.py`.

## Validation Commands
- `cd tst && python run_tests.py particles/pic_mhd_passive_mode`
- `cd tst && python run_tests.py particles/pic_mhd_passive_mode --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx`
- `cd tst && python run_tests.py particles/pic_decomp_invariance particles/pic_deposit_conservation particles/pic_mhd_current_coupling particles/pic_mhd_coupling_nonperiodic particles/pic_mhd_coupling_multilevel particles/pic_mhd_coupling_decomp particles/pic_mhd_restart_fidelity particles/pic_mhd_passive_mode`
- `cd tst && python run_tests.py particles/pic_decomp_invariance particles/pic_deposit_conservation particles/pic_mhd_current_coupling particles/pic_mhd_coupling_nonperiodic particles/pic_mhd_coupling_multilevel particles/pic_mhd_coupling_decomp particles/pic_mhd_restart_fidelity particles/pic_mhd_passive_mode --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx`
- `cd tst/scripts/style && bash check_athena_cpp_style_changed.sh`
- `python3 -m py_compile tst/scripts/particles/pic_mhd_passive_mode.py`

## Validation Outcomes
- Dedicated passive-mode regression: PASS (serial + MPI).
- Full PIC suite (now 8 tests): PASS (serial + MPI).
- Changed-file C++ style gate: PASS.
- New Python test script syntax check: PASS.

## Passive-Mode Metrics
- Particle integrated totals match deterministic expected values:
  - `Q=1024`, `Jx=512`, `Jy=-256`, `Jz=128`, `npart=1024`.
- MHD drift checks from first to last output in passive mode:
  - `bcc1`, `bcc2`, `bcc3`, `mom1`, `mom2`, `mom3`, `ener` all had
    `L2 diff = 0.0` in serial and MPI runs.

## Guard Evidence
Expected failures observed for:
1. passive mode without active `<mhd>` block,
2. passive mode with `pic_feedback_mode=coupled`,
3. `pic_feedback_mode=test_particle` with coupling toggles enabled.

## Artifacts
- `tst/.codex/pic_entity_suite/step2_passive_mhd/pic_mhd_passive_mode_serial.log`
- `tst/.codex/pic_entity_suite/step2_passive_mhd/pic_mhd_passive_mode_mpi.log`
- `tst/.codex/pic_entity_suite/step2_passive_mhd/pic_suite_serial.log`
- `tst/.codex/pic_entity_suite/step2_passive_mhd/pic_suite_mpi.log`
- `tst/.codex/pic_entity_suite/step2_passive_mhd/style_cpp.log`

## Completion Gate Check
- MHD conserved-variable drift below tolerance while PIC moments remain active: PASS
