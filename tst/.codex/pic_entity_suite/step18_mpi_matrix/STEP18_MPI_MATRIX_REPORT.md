# Step 18 Validation Report: MPI Matrix + Decomposition Robustness

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: 720f5716

## Scope

Implemented Step 18 matrix/robustness extensions for late-stage Entity-mirroring
proxies:
- updated regressions:
  - `tst/scripts/particles/pic_crsi_deltaf_proxy.py`
  - `tst/scripts/particles/pic_crpai_polarization_proxy.py`
  - `tst/scripts/particles/pic_expanding_box_anisotropy_proxy.py`
  - `tst/scripts/particles/pic_refinement_boundary_characterization.py`
  - `tst/scripts/particles/pic_amr_shock_lb_smoke.py`

Coverage added:
- explicit `np=4` checks (in addition to existing `np=1/2`) for CRSI/CRPAI and
  expanding/compressing-box proxies.
- explicit `np=1/2/4` parity and alternate-meshblock decomposition checks for
  refinement-boundary characterization and AMR shock/load-balance smoke tests.

## Commands

```bash
cd tst
python run_tests.py \
  particles/pic_crsi_deltaf_proxy \
  particles/pic_crpai_polarization_proxy \
  particles/pic_expanding_box_anisotropy_proxy \
  particles/pic_refinement_boundary_characterization \
  particles/pic_amr_shock_lb_smoke \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx

cd tst/scripts/style
bash check_athena_cpp_style_changed.sh
bash check_python_style_changed.sh

python -m py_compile \
  tst/scripts/particles/pic_crsi_deltaf_proxy.py \
  tst/scripts/particles/pic_crpai_polarization_proxy.py \
  tst/scripts/particles/pic_expanding_box_anisotropy_proxy.py \
  tst/scripts/particles/pic_refinement_boundary_characterization.py \
  tst/scripts/particles/pic_amr_shock_lb_smoke.py
```

## Results

- Full Step-18 matrix run: `5/5` passed.
- CRSI proxy:
  - `mpi4_on:dom_gamma = 5.38855380e-01`
  - `mpi4_on:dom_r2 = 1.12168889e-01`
  - `|log(serial_vs_mpi4 dom_gamma)| = 9.39466098e-01` (within gate)
- CRPAI proxy:
  - both prolate and oblate branches passed `np=1/2/4` dominant-growth gates
    and cross-rank log-ratio parity bounds.
- Expanding/compressing-box proxy:
  - slope-sign and separation gates passed at `np=1/2/4` with tiny cross-rank
    slope differences.
- Refinement-boundary characterization:
  - SMR/AMR baseline bounded metrics passed at `np=1/2/4`.
  - alternate meshblock decomposition runs passed bounded-characterization gates.
- AMR shock/load-balance smoke:
  - baseline and alternate decomposition passed in serial.
  - `np=2/4` parity passed.
  - retained explicit load-balance telemetry at `np=4`:
    - `n_communicated = 14`
    - `lb_efficiency = 1.0`
- Style gates:
  - changed-file C++ style: passed (no changed `src/` files)
  - changed-file Python style: passed for all updated scripts.
- Python bytecode checks: passed.

## Evidence Files

- `tst/.codex/pic_entity_suite/step18_mpi_matrix/pic_step18_mpi_matrix.log`
- `tst/.codex/pic_entity_suite/step18_mpi_matrix/style_cpp.log`
- `tst/.codex/pic_entity_suite/step18_mpi_matrix/style_py.log`

## Gate Verdict

Step 18 completion gate status: PASS.
