# Step 11 Validation Report: Bell Growth Anchor

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: 254dd577

## Scope

Implemented a current-driven Bell-style growth anchor in coupled MHD-PIC mode:
- added deck: `inputs/tests/pic_bell_growth_proxy.athinput`
- added regression: `tst/scripts/particles/pic_bell_growth_proxy.py`

Design note:
- This checkpoint uses a coupled-vs-uncoupled comparison on the same setup.
  The gate requires positive transverse magnetic growth only when
  `couple_moments_to_mhd=true`, mirroring Bell-type current-driven behavior in
  the AthenaK coupling model.

## Commands

```bash
cd tst
python run_tests.py particles/pic_bell_growth_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

python run_tests.py particles/pic_mhd_current_coupling \
  particles/pic_two_stream_growth_proxy particles/pic_weibel_growth_proxy \
  particles/pic_bell_growth_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

cd tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh

python -m py_compile tst/scripts/particles/pic_bell_growth_proxy.py
```

## Results

- `particles/pic_bell_growth_proxy`: passed
  - serial uncoupled gamma: `2.64953581e-17`
  - serial coupled gamma: `9.81431661e-03`
  - serial coupled ratio (`B_perp_last/B_perp_first`): `1.10405993e+00`
  - mpi2 uncoupled gamma: `2.64953581e-17`
  - mpi2 coupled gamma: `9.81431661e-03`
  - serial vs mpi2 coupled ratio abs error: `0`
- Focused regression matrix:
  - `particles/pic_mhd_current_coupling`: passed
  - `particles/pic_two_stream_growth_proxy`: passed
  - `particles/pic_weibel_growth_proxy`: passed
  - `particles/pic_bell_growth_proxy`: passed
- Changed-file C++ style check: no tracked changed C++ files to lint
- Changed-file Python style check: no tracked changed Python files to lint
- Python bytecode check for new script: passed

## Evidence Files

- `tst/.codex/pic_entity_suite/step11_bell/pic_bell.log`
- `tst/.codex/pic_entity_suite/step11_bell/pic_step11_regression_matrix.log`
- `tst/.codex/pic_entity_suite/step11_bell/style_cpp.log`
- `tst/.codex/pic_entity_suite/step11_bell/style_py.log`

## Gate Verdict

Step 11 completion gate status: PASS.
