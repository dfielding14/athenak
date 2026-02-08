# Step 7 Validation Report: EM-Vacuum Analytic Anchor

Date: 2026-02-08
Branch: dev/PIC
Baseline HEAD before step commit: d6cc79d6

## Scope

Implemented an EM-vacuum-style analytic convergence anchor adapted to AthenaK's
current MHD-PIC architecture:
- added deck: `inputs/tests/pic_em_vacuum_wave.athinput`
- added regression: `tst/scripts/particles/pic_em_vacuum_wave.py`

Design note:
- AthenaK does not yet expose a standalone SRPIC Maxwell vacuum solver path in
  this branch. The anchor therefore uses the `linear_wave` MHD analytic setup
  with particles active in strict test-particle mode (`pic_feedback_mode=test_particle`,
  no coupling), preserving the intended convergence/parity gate while keeping
  PIC runtime paths exercised.

## Commands

```bash
cd tst
python run_tests.py particles/pic_em_vacuum_wave \
  --cmake=-DAthena_ENABLE_MPI=ON \
  --cmake=-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx
```

## Results

- `particles/pic_em_vacuum_wave`: passed
- serial convergence ratio (`n32/n16`): `3.13717349e-01`
- mpi2 convergence ratio (`n32/n16`): `3.13717349e-01`
- serial vs mpi2 parity (`n32` L1): exact match (`abs_err=0`)

## Evidence Files

- `tst/.codex/pic_entity_suite/step7_em_vacuum/pic_em_vacuum.log`

## Gate Verdict

Step 7 completion gate status: PASS.
