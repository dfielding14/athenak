# PIC Entity Plan Step 3 Report: no_mhd PIC-Only Carrier Path

## Scope
- Enabled PIC operation without `<mhd>` by adding a particle-owned field carrier
  path for Boris pushers.
- Removed the hard dependency on `pmhd` for Boris interpolation when
  `pic_background_mode=no_mhd`.
- Added dedicated regression coverage for no-MHD Boris execution and guardrails.

## Code Changes
- `src/particles/particles.hpp`
  - Added no-MHD background field controls:
    - `pic_no_mhd_bx`, `pic_no_mhd_by`, `pic_no_mhd_bz`
  - Added particle-owned field carrier array:
    - `pic_no_mhd_bcc0`
- `src/particles/particles.cpp`
  - `pic_feedback_mode` default now switches to `test_particle` for
    `pic_background_mode=no_mhd`.
  - Parsed no-MHD field controls (`pic_no_mhd_bx/by/bz`).
  - Added guard: `pic_background_mode=no_mhd` requires
    `pic_feedback_mode=test_particle`.
  - Allocated/initialized `pic_no_mhd_bcc0` when `no_mhd` mode is active.
- `src/particles/particles_pushers.cpp`
  - Boris path now accepts either:
    - MHD `bcc0` carrier, or
    - no-MHD particle carrier (`pic_no_mhd_bcc0`) when mode is `no_mhd`.
  - Added explicit fatal when neither carrier path is valid.
  - Updated interpolation routines (`InterpolateLinear`, `InterpolateTSC`) to
    select the active field carrier.

## New Test Input and Script
- Added `inputs/tests/pic_no_mhd_boris.athinput` (no `<mhd>` block).
- Added `tst/scripts/particles/pic_no_mhd_boris.py`.

## Validation Commands
- `cd tst && python run_tests.py particles/pic_no_mhd_boris`
- `cd tst && python run_tests.py particles/pic_no_mhd_boris --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx`
- `cd tst && python run_tests.py particles/pic_decomp_invariance particles/pic_deposit_conservation particles/pic_mhd_current_coupling particles/pic_mhd_coupling_nonperiodic particles/pic_mhd_coupling_multilevel particles/pic_mhd_coupling_decomp particles/pic_mhd_restart_fidelity particles/pic_mhd_passive_mode particles/pic_no_mhd_boris`
- `cd tst && python run_tests.py particles/pic_decomp_invariance particles/pic_deposit_conservation particles/pic_mhd_current_coupling particles/pic_mhd_coupling_nonperiodic particles/pic_mhd_coupling_multilevel particles/pic_mhd_coupling_decomp particles/pic_mhd_restart_fidelity particles/pic_mhd_passive_mode particles/pic_no_mhd_boris --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx`
- `cd tst/scripts/style && bash check_athena_cpp_style_changed.sh`
- `python3 -m flake8 tst/scripts/particles/pic_no_mhd_boris.py`

## Validation Outcomes
- Dedicated no-MHD Boris regression: PASS (serial + MPI).
- Full PIC suite (now 9 tests): PASS (serial + MPI).
- Changed-file C++ style gate: PASS.
- New Python no-MHD test script lint: PASS.

## no_mhd Metrics
- Positive-case integrated totals matched deterministic expectations within
  tolerance:
  - `Q=1024`, `npart=1024`
  - `Jx ~ 409.600006`, `Jy ~ -204.800003`, `Jz ~ 102.400002`
- Serial/MPI parity remained exact within configured tolerances.

## Guard Evidence
Expected failures observed for:
1. Boris with `pic_background_mode=coupled` and no MHD carrier,
2. `pic_background_mode=no_mhd` with `pic_feedback_mode=coupled`,
3. no-MHD test-particle mode with `couple_moments_to_mhd=true`.

## Artifacts
- `tst/.codex/pic_entity_suite/step3_no_mhd/pic_no_mhd_boris_serial.log`
- `tst/.codex/pic_entity_suite/step3_no_mhd/pic_no_mhd_boris_mpi.log`
- `tst/.codex/pic_entity_suite/step3_no_mhd/pic_suite_serial.log`
- `tst/.codex/pic_entity_suite/step3_no_mhd/pic_suite_mpi.log`
- `tst/.codex/pic_entity_suite/step3_no_mhd/style_cpp.log`
- `tst/.codex/pic_entity_suite/step3_no_mhd/style_py.log`

## Completion Gate Check
- AthenaK runs PIC test decks with `<particles>` and no `<mhd>` block when
  `pic_background_mode=no_mhd`: PASS
