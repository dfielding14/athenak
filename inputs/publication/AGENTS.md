<!-- BEGIN build-memory-table -->
# PIC Campaign Deck Navigation

## Scope

This directory contains publication-oriented and preproduction Bell/shock deck
freezes plus generated parameter matrices. Names such as `publication`,
`paper`, `production`, and `pilot` do not grant launch authority, physical
qualification, or scientific-claim closure. Read each deck's comments and
manifest together with `tst/publication/AGENTS.md`.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change compact nonlinear Bell engineering setup | `pic_bell_nonlinear_saturation_engineering_v1.athinput` | `src/pgen/tests/q019_nonlinear_bell_saturation_engineering.cpp` and its analyzer/tests |
| Change the physics-first nonlinear Bell matrix | `q019_physics_first_nonlinear_bell_successor_v2/` | Same-named materializer/validator under `tst/publication/` and compiled Q019 fingerprint/initializer |
| Change Q019 runtime-controller overlays | `q019_nonlinear_bell_runtime_controller_v1/` | Same-named generator/tests and its source matrix manifest |
| Change the Q023-carrier redesign | `q019_q023_carrier_nonlinear_bell_redesign_v1/` | Same-named materializer/tests; preserve its external-current-surrogate classification |
| Change carrier resource calibration | `q019_q023_carrier_resource_calibration_runtime_controller_v1/` | Carrier and controller manifests plus the calibration validator |
| Change Section 5.4 parallel-shock decks | `pic_parallel_shock_section54_*` | `src/pgen/tests/pic_parallel_shock.cpp` and Q011 parsers, contracts, readiness bindings, and analyzers |
| Change compact full-Hall nonlinear Bell pilot | `pic_cr_hall_bell_nonlinear_pilot.athinput` | Q019 initializer, `MHD_PIC_CR_HALL_CODE_MAP.md`, and its focused reducer |
| Change controlled delayed-injection Hall shock | `pic_parallel_shock_cr_hall_controlled_pilot.athinput` | Shock generator, restart ledger, and controlled-shock reducer |
| Change long full-Hall shock successor | `pic_parallel_shock_cr_hall_long_full.athinput` | Compact pilot, shock generator, restart ledger, and controlled-shock reducer |
| Change 3D full-Hall shock check | `pic_parallel_shock_cr_hall_3d_full.athinput` | Qualified long 2D shock, 3D source-scaling probe, shock generator, and 3D reducer |

## Change flow

Generated Bell directories are manifest-bound outputs: change the owning
materializer, regenerate the complete directory and manifest, update compiled
fingerprints where required, then run the focused validators. Do not hand-edit
one member of a generated matrix.

Parallel-shock deck changes must remain synchronized with the runtime generator,
the selected MHD-PIC physical mode, injection/gas-subtraction contracts, output
schema, and Q011 source/deck tests. A byte change can invalidate SHA-256 bindings
in manifests, readiness records, or control-plane inventories.

## Focused validation

From the repository root:

- `python3 -m tst.publication.q019_physics_first_nonlinear_bell_successor_v2 --validate-checked-in-decks`
- `python3 -m tst.publication.q019_q023_carrier_nonlinear_bell_redesign_v1 --validate-checked-in-decks`
- `python3 -m tst.publication.q019_q023_carrier_resource_calibration_runtime_controller_v1 --validate`
- `python3 -m pytest -q tst/publication/test_bell_saturation_engineering_analysis_v1.py`
- `python3 -m pytest -q tst/publication/test_q011_section54_production_science_successor_v1.py`

Use the matching Q019 controller tests or Q011 preparation/exact-conservation
tests when those families change. A validator is part of the contract even when
it reports that checked-in artifacts need regeneration.

## Local constraints

- Do not silently relax fail-closed role, authority, or qualification fields.
- Keep historical `paper_mhd_pic`, active `paper_mhd_pic_vl2_tsc`, corrected
  current-normalization, and experimental extension decks distinct.
- A source-local deck or passing static contract does not by itself authorize
  cluster execution or establish a scientific result.
<!-- END build-memory-table -->
