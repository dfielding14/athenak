# Step 16 Validation Report: Refinement-Boundary Characterization Proxies

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: 7a493ab8

## Scope

Implemented Step 16 refinement-boundary characterization proxies and regression:
- added decks:
  - `inputs/tests/pic_refinement_boundary_smr_proxy.athinput`
  - `inputs/tests/pic_refinement_boundary_amr_proxy.athinput`
- added regression:
  - `tst/scripts/particles/pic_refinement_boundary_characterization.py`

Design note:
- This step intentionally uses characterization-style gates for coarse/fine
  interfaces. It validates bounded smoothness and conservation-drift metrics,
  plus serial/MPI parity, instead of forcing strict conservation that is not
  expected for all refinement-boundary deposition cases.

## Commands

```bash
cd tst
python run_tests.py particles/pic_refinement_boundary_characterization \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

python run_tests.py particles/pic_expanding_box_anisotropy_proxy \
  particles/pic_refinement_boundary_characterization \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

cd tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh

python -m py_compile \
  tst/scripts/particles/pic_refinement_boundary_characterization.py
```

## Results

- `particles/pic_refinement_boundary_characterization`: passed
  - `smr:q_drift_np1 = 2.34100000e+04`
  - `smr:cv_np1 = 1.14172080e+00`
  - `smr:jump_np1 = 4.56215293e+00`
  - `amr_proxy:q_drift_np1 = 1.82031250e+04`
  - `amr_proxy:cv_np1 = 8.01303876e-01`
  - `amr_proxy:jump_np1 = 1.46700086e+00`
  - serial/MPI parity checks passed for both SMR and AMR-style runs.
- Focused regression matrix:
  - `particles/pic_expanding_box_anisotropy_proxy`: passed
  - `particles/pic_refinement_boundary_characterization`: passed
- Changed-file C++ style check: passed (no changed `src/` files to lint).
- Changed-file Python style check: passed (no changed tracked Python files
  under `tst/` or `vis/` to lint).
- Python bytecode check for the new script: passed.

## Evidence Files

- `tst/.codex/pic_entity_suite/step16_refinement_boundary/pic_refinement_boundary.log`
- `tst/.codex/pic_entity_suite/step16_refinement_boundary/pic_step16_regression_matrix.log`
- `tst/.codex/pic_entity_suite/step16_refinement_boundary/style_cpp.log`
- `tst/.codex/pic_entity_suite/step16_refinement_boundary/style_py.log`

## Gate Verdict

Step 16 completion gate status: PASS.
