# Step 14 Validation Report: CRPAI Polarization Proxy

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: 41b53ab9

## Scope

Implemented a CRPAI-style anisotropy-branch proxy with prolate/oblate flavors
and polarization-branch gating:
- added decks:
  - `inputs/tests/pic_crpai_prolate_proxy.athinput`
  - `inputs/tests/pic_crpai_oblate_proxy.athinput`
- added regression:
  - `tst/scripts/particles/pic_crpai_polarization_proxy.py`

Design note:
- This checkpoint validates branch-selective behavior and growth stability in
  AthenaK's current proxy ladder; it is not yet a full analytical CRPAI
  dispersion benchmark.

## Commands

```bash
cd tst
python run_tests.py particles/pic_crpai_polarization_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

python run_tests.py particles/pic_crsi_deltaf_proxy \
  particles/pic_crpai_polarization_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

cd tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh

python -m py_compile \
  tst/scripts/particles/pic_crpai_polarization_proxy.py
```

## Results

- `particles/pic_crpai_polarization_proxy`: passed
  - prolate (`np=1`) dominant branch growth:
    - `dom_gamma = 2.32579119e+00`
  - oblate (`np=1`) dominant branch growth:
    - `dom_gamma = 2.08539573e+00`
  - branch-selectivity gate:
    - prolate and oblate choose opposite dominant polarization branches
  - serial/MPI growth consistency (log-ratio metric):
    - prolate: `|log(gamma_np2/gamma_np1)| = 3.54888129e-01`
    - oblate: `|log(gamma_np2/gamma_np1)| = 2.86423063e-01`
- Focused regression matrix:
  - `particles/pic_crsi_deltaf_proxy`: passed
  - `particles/pic_crpai_polarization_proxy`: passed
- Changed-file C++ style check: no tracked changed C++ files to lint
- Changed-file Python style check: no tracked changed Python files to lint
- Python bytecode check for new script: passed

## Evidence Files

- `tst/.codex/pic_entity_suite/step14_crpai_proxy/pic_crpai_proxy.log`
- `tst/.codex/pic_entity_suite/step14_crpai_proxy/pic_step14_regression_matrix.log`
- `tst/.codex/pic_entity_suite/step14_crpai_proxy/style_cpp.log`
- `tst/.codex/pic_entity_suite/step14_crpai_proxy/style_py.log`

## Gate Verdict

Step 14 completion gate status: PASS.
