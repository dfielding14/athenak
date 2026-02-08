# Step 4 Validation Report: Midpoint E+B Boris

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: 0a0bc63c

## Scope

Implemented and validated Step 4 requirements from
`AGENT_PIC_ENTITY_TEST_PLAN.md`:
- midpoint `dt/2 -> sample (u,B) -> cE=-u x B -> Boris(dt) -> dt/2` push
- CR payload extension for sampled midpoint fields and per-step deltas
- coupled feedback path consuming deposited `dp/dt` and `dE/dt`
- midpoint frozen-in orthogonality diagnostic (`cE dot B`)

## Commands

```bash
cd tst
python run_tests.py \
  particles/pic_deposit_conservation \
  particles/pic_mhd_current_coupling \
  particles/pic_mhd_coupling_decomp \
  particles/pic_mhd_coupling_multilevel \
  particles/pic_mhd_coupling_nonperiodic \
  particles/pic_mhd_restart_fidelity \
  particles/pic_mhd_passive_mode \
  particles/pic_no_mhd_boris \
  particles/pic_boris_midpoint_eb \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx
```

```bash
cd tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh
```

## Results

- PIC regression matrix: `9/9` passed
- Changed-file C++ style gate: passed
- Changed-file Python style gate: passed

Key midpoint-E+B metrics from `particles/pic_boris_midpoint_eb`:
- `passive_zero:vavg |meas-exp| = 1.49335796e-07`
- `passive_flow:vavg |meas-exp| = 3.66774009e-07`
- `passive_speed_delta = 3.35087812e-03`
- `frozen_in_cE_dot_B abs_err = 0.0`
- `coupled_momentum_exchange |meas-exp| = 2.67495559e-06`
- `coupled_energy_exchange abs_err = 3.91224054e-06`
- serial/MPI parity checks for midpoint test: all `0.0` measured difference

## Evidence Files

- `tst/.codex/pic_entity_suite/step4_midpoint_eb/pic_suite_mpi.log`
- `tst/.codex/pic_entity_suite/step4_midpoint_eb/style_cpp.log`
- `tst/.codex/pic_entity_suite/step4_midpoint_eb/style_py.log`

## Gate Verdict

Step 4 completion gate status: PASS.
