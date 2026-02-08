# Step 5 Validation Report: PIC Diagnostics Layer

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: b53b1662

## Scope

Implemented Step 5 diagnostics additions for PIC analysis workflows:
- new particle-moment derived outputs:
  - `prtcl_dpxdt`, `prtcl_dpydt`, `prtcl_dpzdt`
  - `prtcl_dedt`
  - `prtcl_ebdot`
- output registration/validation plumbing for the new variables
- derived-variable kernels mapping each variable to particle moments storage
- midpoint Boris regression updated to consume and validate diagnostics
- shared Python analysis utilities added for spectrum/growth/frequency fitting:
  - `tst/scripts/particles/pic_analysis_utils.py`

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

Key Step 5 diagnostic metrics from `particles/pic_boris_midpoint_eb`:
- `frozen_in_ebdot_output abs_err = 6.39277365e-16`
- `coupled_momentum_exchange_diag |meas-exp| = 2.99883474e-06`
- `coupled_energy_exchange_diag abs_err = 3.88927720e-06`

## Evidence Files

- `tst/.codex/pic_entity_suite/step5_diagnostics/pic_suite_mpi.log`
- `tst/.codex/pic_entity_suite/step5_diagnostics/style_cpp.log`
- `tst/.codex/pic_entity_suite/step5_diagnostics/style_py.log`

## Gate Verdict

Step 5 completion gate status: PASS.
