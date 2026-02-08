# Step 17 Validation Report: Reduced-Size AMR Shock + Load-Balancing Smoke

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: c0ecb7bd

## Scope

Implemented Step 17 reduced-size AMR shock/load-balancing smoke coverage:
- added deck:
  - `inputs/tests/pic_amr_shock_lb_smoke.athinput`
- added regression:
  - `tst/scripts/particles/pic_amr_shock_lb_smoke.py`

Design note:
- The smoke case uses an AthenaK-compatible shock-rich Orszag-Tang setup with
  deterministic runtime overrides (`ddens_max`, `nlim`, `tlim`) to trigger AMR
  quickly and exercise load balancing in bounded runtime.
- Field-evolution gate uses growth of `std(bcc1)` (a shock-structure
  amplification proxy) instead of total-field magnitude, which is not a stable
  discriminator in this compact configuration.
- Non-thermal behavior is gated via particle-current tail growth proxies.

## Commands

```bash
cd tst
python run_tests.py particles/pic_amr_shock_lb_smoke \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

python run_tests.py particles/pic_refinement_boundary_characterization \
  particles/pic_amr_shock_lb_smoke \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

cd tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh

python -m py_compile \
  tst/scripts/particles/pic_amr_shock_lb_smoke.py
```

## Results

- `particles/pic_amr_shock_lb_smoke`: passed
  - Serial (`np=1`) metrics:
    - `b1_std_final_amp = 1.00580433e+00`
    - `tail_peak = 3.84911427e+00`
    - `n_created = 1.96000000e+02`
  - MPI (`np=4`) metrics:
    - `b1_std_final_amp = 1.00580433e+00`
    - `tail_peak = 3.87764362e+00`
    - `n_created = 1.96000000e+02`
    - `n_communicated = 1.40000000e+01`
    - `lb_efficiency = 1.00000000e+00`
  - serial/MPI parity gates passed for `b1_std_peak_amp` and `tail_peak`.
- Focused regression matrix:
  - `particles/pic_refinement_boundary_characterization`: passed
  - `particles/pic_amr_shock_lb_smoke`: passed
- Changed-file C++ style check: passed (no changed `src/` files to lint).
- Changed-file Python style check: passed (no changed tracked Python files
  under `tst/` or `vis/` to lint).
- Python bytecode check for new script: passed.

## Evidence Files

- `tst/.codex/pic_entity_suite/step17_amr_shock_lb/pic_amr_shock_lb_smoke.log`
- `tst/.codex/pic_entity_suite/step17_amr_shock_lb/pic_step17_regression_matrix.log`
- `tst/.codex/pic_entity_suite/step17_amr_shock_lb/style_cpp.log`
- `tst/.codex/pic_entity_suite/step17_amr_shock_lb/style_py.log`

## Gate Verdict

Step 17 completion gate status: PASS.
