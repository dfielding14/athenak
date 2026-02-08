# Step 10 Validation Report: Weibel Growth Proxy

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: 1ed388de

## Scope

Implemented the Step 10 checkpoint as an AthenaK-adapted Weibel proxy:
- added deck: `inputs/tests/pic_weibel_growth_proxy.athinput`
- added regression: `tst/scripts/particles/pic_weibel_growth_proxy.py`

Design note:
- The current no-MHD PIC path does not include a self-consistent electromagnetic
  field solver for Weibel growth. This checkpoint therefore uses a transverse
  counter-beam control and enforces non-positive fitted growth of the dominant
  transverse-current mode with serial/MPI parity.

## Commands

```bash
cd tst
python run_tests.py particles/pic_weibel_growth_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

python run_tests.py particles/pic_no_mhd_boris \
  particles/pic_two_stream_growth_proxy particles/pic_weibel_growth_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

cd tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh

python -m py_compile tst/scripts/particles/pic_weibel_growth_proxy.py
```

## Results

- `particles/pic_weibel_growth_proxy`: passed
  - `np1` fitted gamma: `-1.00559957e-16`
  - `np2` fitted gamma: `-3.35199856e-17`
  - `serial_vs_mpi2` gamma difference: `6.70399712e-17`
- Focused regression matrix:
  - `particles/pic_no_mhd_boris`: passed
  - `particles/pic_two_stream_growth_proxy`: passed
  - `particles/pic_weibel_growth_proxy`: passed
- Changed-file C++ style check: no tracked changed C++ files to lint
- Changed-file Python style check: no tracked changed Python files to lint
- Python bytecode check for new script: passed

## Evidence Files

- `tst/.codex/pic_entity_suite/step10_weibel/pic_weibel.log`
- `tst/.codex/pic_entity_suite/step10_weibel/pic_step10_regression_matrix.log`
- `tst/.codex/pic_entity_suite/step10_weibel/style_cpp.log`
- `tst/.codex/pic_entity_suite/step10_weibel/style_py.log`

## Gate Verdict

Step 10 completion gate status: PASS (AthenaK-adapted growth proxy).
