# Step 6 Validation Report: Entity Deposit Parity Layer

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: 8708ee6d

## Scope

Implemented and validated Step 6 parity foundations:
- Added/validated Entity-mirroring regression decks:
  - `inputs/tests/pic_entity_deposit_mink.athinput`
  - `inputs/tests/pic_entity_deposit_reflect.athinput`
- Added/validated regression scripts:
  - `tst/scripts/particles/pic_entity_deposit_mink.py`
  - `tst/scripts/particles/pic_entity_deposit_reflect.py`
- Tightened decomposition checks to include field-level parity (`rho/jx/jy/jz`)
  in both periodic and reflect-boundary tests.

Note:
- Dedicated coarse/fine refinement-boundary characterization is tracked as
  Step 16 and remains there to avoid duplicate gate logic.

## Commands

```bash
cd tst
python run_tests.py \
  particles/pic_entity_deposit_mink \
  particles/pic_entity_deposit_reflect \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx
```

## Results

- Entity parity matrix: `2/2` passed

Key metrics from evidence log:
- `pic_entity_deposit_mink`: exact `Q/J/npart` agreement with expected totals and
  exact field-level parity (`max_abs_diff = 0`) across `np=1,2,3,4`.
- `pic_entity_deposit_reflect`: expected non-conservative particle retention in
  reflect setup (`retention = 0.75`), with exact decomposition parity for
  integrated moments and field-level differences (`max_abs_diff = 0`) across
  `np=1,2,4`.

## Evidence Files

- `tst/.codex/pic_entity_suite/step6_entity_parity/entity_deposit_mpi.log`

## Gate Verdict

Step 6 completion gate status: PASS (parity-layer component).
