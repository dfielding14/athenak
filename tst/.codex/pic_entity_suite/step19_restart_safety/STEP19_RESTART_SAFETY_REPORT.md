# Step 19 Restart and Safety Guard Report

## Scope

Plan step:
- `Step 19: Restart and Safety Guard Coverage`

Head at validation time:
- `7184ef56` (`PIC plan step 18: extend MPI/decomposition robustness matrix`)

New artifacts:
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/inputs/tests/pic_restart_safety_guards.athinput`
- `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/particles/pic_restart_safety_guards.py`

## Commands and Logs

1. MPI restart and safety matrix (new + baseline restart suite)
- Command:
  - `cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst`
  - `python run_tests.py particles/pic_restart_safety_guards particles/pic_mhd_restart_fidelity --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx`
- Log:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_entity_suite/step19_restart_safety/pic_step19_restart_safety.log`

2. Style and syntax checks
- C++ changed-file gate log:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_entity_suite/step19_restart_safety/style_cpp.log`
- Python changed-file gate log:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_entity_suite/step19_restart_safety/style_py.log`
- Python compile check log:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_entity_suite/step19_restart_safety/py_compile.log`
- Flake8 log:
  - `/Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_entity_suite/step19_restart_safety/flake8.log`

## Results Summary

Overall:
- `particles/pic_restart_safety_guards`: passed
- `particles/pic_mhd_restart_fidelity`: passed
- Combined summary: `2 out of 2 test passed`

Restart A/B parity (representative modes in new test):
- `no_mhd` (`np=1,2`): full vs restart match to `abs_tol=1e-8` / `rel_tol=1e-8`
- `passive_mhd` (`np=1,2`): full vs restart match to `abs_tol=1e-8` / `rel_tol=1e-8`
- `coupled_edge_direct` (`np=1,2`): full vs restart match to `abs_tol=1e-8` / `rel_tol=1e-8`

Safety guard checks (new deterministic negative cases):
- `guard_no_mhd_requires_test_particle`:
  matched reason substring
  `"<particles>/pic_background_mode=no_mhd requires <particles>/pic_feedback_mode=test_particle"`
- `guard_passive_mhd_requires_test_particle`:
  matched reason substring
  `"<particles>/pic_background_mode=passive_mhd requires <particles>/pic_feedback_mode=test_particle unless coupled feedback is explicitly implemented for this mode"`
- `guard_test_particle_rejects_coupling_toggles`:
  matched reason substring
  `"<particles>/pic_feedback_mode=test_particle does not support particle-to-MHD coupling toggles"`

Per-rank restart watch item probe:
- Probe mode: coupled direct-edge (`np=2`) with
  `output7/single_file_per_rank=true`
- Observed status: `supported`
- Full vs restart parity also matched with
  `abs_tol=1e-8` / `rel_tol=1e-8`
- Watch item retained for broader coupled-PIC workflows and future topology
  permutations, even though this probe is green.

## Completion Gate Assessment

Step 19 gate requirements:
1. Restart fidelity tests pass for supported paths: **pass**
2. Unsupported paths fail with deterministic reason strings: **pass**
3. Per-rank restart behavior/messages remain explicitly tracked: **pass**

Step 19 status:
- **Complete**
