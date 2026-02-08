# Step 12 Validation Report: Multi-Species Backreaction Oscillation

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: fba5d1fe

## Scope

Implemented a multi-species coupled oscillation checkpoint with mesh-hierarchy
parity coverage:
- added decks:
  - `inputs/tests/pic_multispecies_osc_uniform.athinput`
  - `inputs/tests/pic_multispecies_osc_smr.athinput`
  - `inputs/tests/pic_multispecies_osc_amr_proxy.athinput`
- added regression:
  - `tst/scripts/particles/pic_multispecies_backreaction_oscillation.py`

Design note:
- The third variant is an AMR-style proxy using nested static refinement
  regions, used here to enforce uniform/SMR/multilevel parity metrics in the
  current AthenaK harness.

## Commands

```bash
cd tst
python run_tests.py particles/pic_multispecies_backreaction_oscillation \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

python run_tests.py particles/pic_bell_growth_proxy \
  particles/pic_multispecies_backreaction_oscillation \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

cd tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh

python -m py_compile \
  tst/scripts/particles/pic_multispecies_backreaction_oscillation.py
```

## Results

- `particles/pic_multispecies_backreaction_oscillation`: passed
  - dominant frequency (`np1`):
    - uniform: `5.00000000e-02`
    - smr: `5.00000000e-02`
    - amr_proxy: `5.00000000e-02`
  - relative MHD energy drift (`np1`):
    - uniform: `1.27988754e-02`
    - smr: `2.21595189e-02`
    - amr_proxy: `2.97997746e-02`
  - `np1` vs `np2` frequency parity: exact match for all three variants
- Focused regression matrix:
  - `particles/pic_bell_growth_proxy`: passed
  - `particles/pic_multispecies_backreaction_oscillation`: passed
- Changed-file C++ style check: no tracked changed C++ files to lint
- Changed-file Python style check: no tracked changed Python files to lint
- Python bytecode check for new script: passed

## Evidence Files

- `tst/.codex/pic_entity_suite/step12_multispecies_osc/pic_multispecies_osc.log`
- `tst/.codex/pic_entity_suite/step12_multispecies_osc/pic_step12_regression_matrix.log`
- `tst/.codex/pic_entity_suite/step12_multispecies_osc/style_cpp.log`
- `tst/.codex/pic_entity_suite/step12_multispecies_osc/style_py.log`

## Gate Verdict

Step 12 completion gate status: PASS.
