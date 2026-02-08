# Step 0 Baseline Report

Date: 2026-02-08 04:41:13 UTC
Branch: dev/PIC
HEAD: 14e9048205b333677355561f2201ddf960479018

## Dirty State Snapshot
```
 D AGENT_PIC_PR1_PR4A_REVIEW_BRIEF.md
 M inputs/AGENTS.md
?? AGENT_PIC_ENTITY_TEST_PLAN.md
?? AGENT_PIC_PR4_REVIEW_BRIEF.md
?? inputs/tests/pic_entity_deposit_mink.athinput
?? inputs/tests/pic_entity_deposit_reflect.athinput
?? tmp/
?? tst/.codex/
?? tst/scripts/particles/pic_entity_deposit_mink.py
?? tst/scripts/particles/pic_entity_deposit_reflect.py
```

## Commands Run
1. Serial PIC baseline:
```bash
cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst
python3 run_tests.py particles/pic_decomp_invariance particles/pic_deposit_conservation particles/pic_mhd_current_coupling particles/pic_mhd_coupling_nonperiodic particles/pic_mhd_coupling_multilevel particles/pic_mhd_coupling_decomp particles/pic_mhd_restart_fidelity --log_file .codex/pic_entity_suite_baseline/serial.log
```
2. MPI-enabled PIC baseline:
```bash
cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst
python3 run_tests.py --cmake=-DAthena_ENABLE_MPI=ON --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx particles/pic_decomp_invariance particles/pic_deposit_conservation particles/pic_mhd_current_coupling particles/pic_mhd_coupling_nonperiodic particles/pic_mhd_coupling_multilevel particles/pic_mhd_coupling_decomp particles/pic_mhd_restart_fidelity --log_file .codex/pic_entity_suite_baseline/mpi.log
```
3. Style gates:
```bash
cd /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh
```

## Results Summary
- Serial baseline: PASS (7/7 tests)
- MPI baseline: PASS (7/7 tests)
- C++ changed-file style gate: PASS
- Python changed-file style gate: PASS

## Key Numeric Continuity/Parity Anchors
- default_edge_mode_delta_vs_cc_convert max = 8.64877701e-01
- continuity_direct_vs_edge_rel_ratio = 1.00000000e+00
- continuity_direct_o2_vs_direct_rel_ratio = 1.00000000e+00

## Artifacts
- /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_entity_suite_baseline/serial.log
- /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_entity_suite_baseline/mpi.log
- /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_entity_suite_baseline/style_cpp.log
- /Users/dbf75/Work/Research/AthenaK/athenak-DF/tst/.codex/pic_entity_suite_baseline/style_py.log
