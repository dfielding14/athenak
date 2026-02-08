# Step 9 Validation Report: Two-Stream Growth Proxy

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: 28d191e7

## Scope

Implemented the Step 9 checkpoint with AthenaK-adapted constraints:
- added optional per-species drift initialization controls in CR setup:
  - `speciesN/vx0`, `speciesN/vy0`, `speciesN/vz0`
- added deck: `inputs/tests/pic_two_stream_growth_proxy.athinput`
- added regression: `tst/scripts/particles/pic_two_stream_growth_proxy.py`

Design note:
- Current AthenaK PIC no-MHD path does not include a self-consistent
  electrostatic field solver, so this checkpoint is implemented as a
  counter-streaming growth-control proxy. The regression fits a density-mode
  growth rate and enforces a non-positive growth bound with serial/MPI parity.

## Commands

```bash
cd tst
python run_tests.py particles/pic_two_stream_growth_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

python run_tests.py particles/pic_no_mhd_boris \
  particles/pic_boris_midpoint_eb particles/pic_two_stream_growth_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

cd tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh

python -m py_compile tst/scripts/particles/pic_two_stream_growth_proxy.py
```

## Results

- `particles/pic_two_stream_growth_proxy`: passed
  - `np1` fitted gamma: `-9.59676264e-02`
  - `np2` fitted gamma: `-4.81984005e-02`
  - `serial_vs_mpi2` gamma difference: `4.77692259e-02`
- Focused regression matrix:
  - `particles/pic_no_mhd_boris`: passed
  - `particles/pic_boris_midpoint_eb`: passed
  - `particles/pic_two_stream_growth_proxy`: passed
- Changed-file C++ style check: passed
- Changed-file Python style check: no tracked changed Python files to lint
- Python bytecode check for new script: passed

## Evidence Files

- `tst/.codex/pic_entity_suite/step9_two_stream/pic_two_stream.log`
- `tst/.codex/pic_entity_suite/step9_two_stream/pic_step9_regression_matrix.log`
- `tst/.codex/pic_entity_suite/step9_two_stream/style_cpp.log`
- `tst/.codex/pic_entity_suite/step9_two_stream/style_py.log`

## Gate Verdict

Step 9 completion gate status: PASS (AthenaK-adapted growth proxy).
