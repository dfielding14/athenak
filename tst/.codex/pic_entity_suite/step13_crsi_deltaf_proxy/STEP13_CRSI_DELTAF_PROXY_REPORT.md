# Step 13 Validation Report: CRSI Delta-f Proxy

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: 3df4b867

## Scope

Implemented a CRSI-style polarization/growth proxy with explicit delta-f runtime
coverage and staged noise-control behavior:
- added runtime behavior:
  - `src/particles/particles.cpp`
    - when `pic_deltaf_mode=on` and `cr_distribution=random`, particle
      initialization switches to a deterministic low-discrepancy quiet start.
- added deck:
  - `inputs/tests/pic_crsi_deltaf_proxy.athinput`
- added regression:
  - `tst/scripts/particles/pic_crsi_deltaf_proxy.py`

Design note:
- This checkpoint validates a stable delta-f execution path and polarization-
  resolved growth metrics in AthenaK's current proxy ladder. It does not claim
  full kinetic CRSI dispersion parity yet.

## Commands

```bash
cd tst
python run_tests.py particles/pic_crsi_deltaf_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

python run_tests.py particles/pic_multispecies_backreaction_oscillation \
  particles/pic_crsi_deltaf_proxy \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

cd tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh

python -m py_compile \
  tst/scripts/particles/pic_crsi_deltaf_proxy.py
```

## Results

- `particles/pic_crsi_deltaf_proxy`: passed
  - serial `deltaf=on` dominant-branch growth:
    - `dom_gamma = 1.37872346e+00`
    - `dom_ratio = 8.07134187e+02`
  - serial `deltaf=on` noise-vs-off metric:
    - `on_noise - off_noise = -1.52127109e-01` (lower is better; pass)
  - MPI (`np=2`) `deltaf=on` dominant-branch growth:
    - `dom_gamma = 9.24290718e-01`
    - serial/MPI log-ratio metric:
      `|log(gamma_np2/gamma_np1)| = 3.99886669e-01`
- Focused regression matrix:
  - `particles/pic_multispecies_backreaction_oscillation`: passed
  - `particles/pic_crsi_deltaf_proxy`: passed
- Changed-file C++ style check: passed (`src/particles/particles.cpp`)
- Changed-file Python style check: no tracked changed Python files to lint
- Python bytecode check for new script: passed

## Evidence Files

- `tst/.codex/pic_entity_suite/step13_crsi_deltaf_proxy/pic_crsi_deltaf_proxy.log`
- `tst/.codex/pic_entity_suite/step13_crsi_deltaf_proxy/pic_step13_regression_matrix.log`
- `tst/.codex/pic_entity_suite/step13_crsi_deltaf_proxy/style_cpp.log`
- `tst/.codex/pic_entity_suite/step13_crsi_deltaf_proxy/style_py.log`

## Gate Verdict

Step 13 completion gate status: PASS.
