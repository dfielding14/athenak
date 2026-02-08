# Step 15 Validation Report: Expanding/Compressing Box Proxy

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: fda596ed

## Scope

Implemented a driven expanding/compressing-box proxy checkpoint and regression:
- added pusher behavior:
  - `src/particles/particles_pushers.cpp`
    - when `pic_expanding_box_mode=on`, applies half-step expansion/
      compression source updates around the Boris advance
      (`exp(-0.5*dt*rate_xi)` per component).
- added decks:
  - `inputs/tests/pic_expanding_box_proxy.athinput`
  - `inputs/tests/pic_compressing_box_proxy.athinput`
- added regression:
  - `tst/scripts/particles/pic_expanding_box_anisotropy_proxy.py`

Design note:
- This is a compact proxy for anisotropy-trend behavior in expanding
  coordinates. It validates expected sign/trend behavior and MPI stability in
  current AthenaK staged controls.

## Commands

```bash
cd tst
python run_tests.py particles/pic_expanding_box_anisotropy_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

python run_tests.py particles/pic_crpai_polarization_proxy \
  particles/pic_expanding_box_anisotropy_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

cd tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh

python -m py_compile \
  tst/scripts/particles/pic_expanding_box_anisotropy_proxy.py
```

## Results

- `particles/pic_expanding_box_anisotropy_proxy`: passed
  - expanding run slope (`np=1`):
    - `slope = -1.66611619e-01`
  - compressing run slope (`np=1`):
    - `slope = 5.15362742e-01`
  - slope separation (`np=1`):
    - `compressing - expanding = 6.81974361e-01`
  - serial/MPI slope parity:
    - expanding `|np2-np1| = 1.83183512e-08`
    - compressing `|np2-np1| = 5.78987747e-09`
- Focused regression matrix:
  - `particles/pic_crpai_polarization_proxy`: passed
  - `particles/pic_expanding_box_anisotropy_proxy`: passed
- Changed-file C++ style check: passed (`src/particles/particles_pushers.cpp`)
- Changed-file Python style check: no tracked changed Python files to lint
- Python bytecode check for new script: passed

## Evidence Files

- `tst/.codex/pic_entity_suite/step15_expanding_box_proxy/pic_expanding_box_proxy.log`
- `tst/.codex/pic_entity_suite/step15_expanding_box_proxy/pic_step15_regression_matrix.log`
- `tst/.codex/pic_entity_suite/step15_expanding_box_proxy/style_cpp.log`
- `tst/.codex/pic_entity_suite/step15_expanding_box_proxy/style_py.log`

## Gate Verdict

Step 15 completion gate status: PASS.
