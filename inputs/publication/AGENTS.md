<!-- BEGIN build-memory-table -->
# PIC Campaign Deck Navigation

## Scope

This directory contains publication-oriented and preproduction Bell/shock deck
freezes. Names such as `publication`,
`paper`, `production`, and `pilot` do not grant launch authority, physical
qualification, or scientific-claim closure. Read each deck's comments and
manifest together with `tst/publication/AGENTS.md`. Active priorities come from
`MHD_PIC_NEXT_STEPS_GUIDE.md`; durable model and setup documents live under
`docs/source/engineering/`.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change compact nonlinear Bell engineering setup | `pic_bell_nonlinear_saturation_engineering_v1.athinput` | `src/pgen/tests/q019_nonlinear_bell_saturation_engineering.cpp` and its analyzer/tests |
| Change Section 5.4 parallel-shock decks | `pic_parallel_shock_section54_*` | `src/pgen/tests/pic_parallel_shock.cpp` and Q011 parsers, contracts, readiness bindings, and analyzers |
| Change compact full-Hall nonlinear Bell pilot | `pic_cr_hall_bell_nonlinear_pilot.athinput` | Q019 initializer, `docs/source/engineering/pic_cr_hall_code_map.md`, and its focused reducer |
| Change controlled delayed-injection Hall shock | `pic_parallel_shock_cr_hall_controlled_pilot.athinput` | Shock generator, restart ledger, and controlled-shock reducer |
| Change long full-Hall shock successor | `pic_parallel_shock_cr_hall_long_full.athinput` | Compact pilot, shock generator, restart ledger, and controlled-shock reducer |
| Change wide full-Hall shock successor | `pic_parallel_shock_cr_hall_wide_full.athinput` | Long full-Hall deck, shock generator, output schema, and shock plotting tools |
| Change 3D full-Hall shock check | `pic_parallel_shock_cr_hall_3d_full.athinput` | Qualified long 2D shock, 3D source-scaling probe, shock generator, and 3D reducer |
| Change the frozen Mignone-R2 full-Hall shock | `pic_parallel_shock_mignone_r2_full_hall.athinput` | `docs/source/engineering/pic_mignone_r2_shock_setup.md`, generator, tracer smoke test, and retained provenance |

## Change flow

Parallel-shock deck changes must remain synchronized with the runtime generator,
explicit VL2/TSC and coupling controls, injection/gas-subtraction contracts,
output schema, and Q011 source/deck tests. A byte change can invalidate SHA-256
bindings in manifests, readiness records, or control-plane inventories.

## Focused validation

From the repository root:

- `python3 -m pytest -q tst/publication/test_bell_saturation_engineering_analysis_v1.py`
- `python3 -m pytest -q tst/publication/test_q011_section54_production_science_successor_v1.py`
- `cd tst && python3 run_tests.py particles/pic_parallel_shock_mignone_tracer_smoke`

Use the matching Q011 preparation/exact-conservation tests when that family
changes.

## Local constraints

- Do not silently relax fail-closed role, authority, or qualification fields.
- Keep explicit VL2/TSC, corrected current-normalization, full-Hall, and
  Hall-off decks distinct.
- A source-local deck or passing static contract does not by itself authorize
  cluster execution or establish a scientific result.
<!-- END build-memory-table -->
