# PIC Entity Plan Step 1 Report: Runtime Isolation Knob Parse/Guard Layer

## Scope
- Added staged PIC runtime-control enums + member storage in:
  - `src/particles/particles.hpp`
- Added parse/validation guards in:
  - `src/particles/particles.cpp`
- Extended negative guard coverage in:
  - `tst/scripts/particles/pic_deposit_conservation.py`
- Added default definitions for new override keys in:
  - `inputs/tests/pic_deposit_conservation.athinput`
- Updated touched AGENTS docs:
  - `AGENTS.md`
  - `inputs/AGENTS.md`
  - `src/particles/AGENTS.md`

## Runtime Controls Added (Step 1 parse/guard only)
- `pic_background_mode = coupled|passive_mhd|no_mhd`
- `pic_feedback_mode = coupled|test_particle`
- `pic_interp_scheme = tsc`
- `pic_cr_light_speed > 0`
- `pic_max_cell_cross > 0`
- `pic_theta_max > 0`
- `pic_deltaf_mode = off|on`
- `pic_deltaf_f0` required if `pic_deltaf_mode=on`
- `pic_sort_interval >= 0`
- `pic_intermediate_arrays = auto|off`
- `pic_expanding_box_mode = off|on`
- `pic_expansion_rate_x1/x2/x3` require `pic_expanding_box_mode=on`

## Validation Commands
- `cd tst && python run_tests.py particles/pic_deposit_conservation`
- `cd tst && python run_tests.py particles/pic_deposit_conservation --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx`
- `cd tst && python run_tests.py particles/pic_decomp_invariance particles/pic_deposit_conservation particles/pic_mhd_current_coupling particles/pic_mhd_coupling_nonperiodic particles/pic_mhd_coupling_multilevel particles/pic_mhd_coupling_decomp particles/pic_mhd_restart_fidelity`
- `cd tst && python run_tests.py particles/pic_decomp_invariance particles/pic_deposit_conservation particles/pic_mhd_current_coupling particles/pic_mhd_coupling_nonperiodic particles/pic_mhd_coupling_multilevel particles/pic_mhd_coupling_decomp particles/pic_mhd_restart_fidelity --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx`
- `cd tst/scripts/style && bash check_athena_cpp_style_changed.sh`
- `cd tst/scripts/style && bash check_python_style_changed.sh`

## Validation Outcomes
- `particles/pic_deposit_conservation`: pass (serial)
- `particles/pic_deposit_conservation`: pass (MPI build; includes `mpi2` path)
- Existing baseline PIC regression suite: pass (7/7 serial)
- Existing baseline PIC regression suite: pass (7/7 MPI)
- Changed-file C++ style gate: pass
- Changed-file Python style gate: pass

## Guard Evidence Summary
All expected invalid override paths failed with deterministic diagnostics:
1. invalid `pic_background_mode`
2. invalid `pic_feedback_mode`
3. invalid `pic_interp_scheme`
4. invalid `pic_deltaf_mode`
5. invalid `pic_intermediate_arrays`
6. invalid `pic_expanding_box_mode`
7. `pic_deltaf_mode=on` without `pic_deltaf_f0`
8. non-zero expansion rate with expanding-box mode off

Negative-check total in test harness: `11/11` observed.

## Artifacts
- `tst/.codex/pic_entity_suite/step1_runtime_knobs/pic_deposit_conservation_serial.log`
- `tst/.codex/pic_entity_suite/step1_runtime_knobs/pic_deposit_conservation_mpi.log`
- `tst/.codex/pic_entity_suite/step1_runtime_knobs/pic_suite_serial.log`
- `tst/.codex/pic_entity_suite/step1_runtime_knobs/pic_suite_mpi.log`
- `tst/.codex/pic_entity_suite/step1_runtime_knobs/style_cpp.log`
- `tst/.codex/pic_entity_suite/step1_runtime_knobs/style_py.log`

## Completion Gate Check
- Existing tests pass unchanged under default mode: PASS
- Invalid mode values fail with clear diagnostics: PASS
