<!-- BEGIN build-memory-table -->
# Particle Test Navigation

## Scope

This subtree validates particle pushers, deposition, MHD-PIC coupling,
boundaries, MPI/mesh decomposition, restart behavior, and bounded physics
proxies. It does not own implementation in `src/` or authorize publication and
production campaigns governed under `tst/publication/`. Active priorities are
in `MHD_PIC_NEXT_STEPS_GUIDE.md`; durable MHD-PIC contracts live under
`docs/source/engineering/`.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change relativistic pushing or field sampling | `pic_relativistic_gyro_paper.py`, `pic_boris_midpoint_eb.py`, `pic_no_mhd_boris.py` | `src/particles/particles_pushers.cpp` and referenced decks |
| Change deposition or TSC behavior | `pic_deposit_conservation.py`, `pic_entity_deposit_*.py`, `pic_paper_smooth_tsc_interface.py` | `src/particles/particles_moments.cpp`, boundary exchange, and the standalone TSC oracle |
| Change MHD backreaction or stage ordering | `pic_paper_coupling_conservation_vl2_tsc.py`, `pic_parallel_shock_rk_stage_budget_vl2_tsc.py`, `pic_mhd_coupling_decomp.py` | `src/mhd/mhd_tasks.cpp`, particle task wiring, and paired VL2 decks |
| Change MPI, AMR, boundary, or restart behavior | `pic_*decomp*.py`, `pic_mhd_coupling_multilevel.py`, `pic_mhd_restart_fidelity.py`, `pic_restart_safety_guards.py` | `src/bvals/`, `src/mesh/`, and restart/output code |
| Change particle migration/destruction compaction | `pic_migration_destruction_compaction.py` | `tst/test_particle_migration_compaction.py` and `src/bvals/particle_compaction.hpp` |
| Change turbulent forcing or its PIC coupling | `pic_turbulent_dynamo_smoke.py` | `tst/test_turb_driver.py`, `src/srcterms/turb_driver.*`, and the MHD/PIC smoke decks |
| Change shock injection or feedback | `pic_parallel_shock_*.py` | `src/pgen/tests/pic_parallel_shock.cpp` and every referenced shock deck |
| Change physical proxy diagnostics | Bell, CRSI, CRPAI, multispecies, expanding-box, and turbulent-dynamo modules | `pic_analysis_utils.py` and the evidence labels in the paired deck/test |

## Test and input flow

A runtime test normally names one or more `inputs/tests/*.athinput` decks,
clears case-specific outputs, launches serial and optional MPI variants with
command-line overrides, then parses logs, `bin`, `pvtk`, or restart output in
`analyze()`. Pairing is not necessarily one-to-one: tests may share decks, use
multiple decks, or delegate to sibling modules. Trace `_INPUT*` constants and
runtime overrides before changing a fixture. Preparation contracts can instead
reference `inputs/publication/`.

Runnable harness modules expose `run()` and `analyze()`. Shared helpers use the
`*_utils.py` or `*_oracle.py` suffix so broad discovery skips them. MPI-aware tests commonly use
`MPIEXEC` (default `mpiexec`) and query `athena -c` before enabling multirank
cases.

MHD-PIC fixtures use explicit controls: Boris cosmic rays store `p/m`,
`pic_cr_initial_state` affects initializer interpretation only,
`time/integrator=vl2` selects staged TSC coupling, and
`pic_cr_hall_mode=off|full` selects induction. Deposition and gas-coupling
toggles default to `false` with `deposit_order=1`; each VL2 fixture must state
the complete coupled order-2 contract.

## Local validation

- `cd tst && python3 run_tests.py particles/pic_deposit_conservation`
- `cd tst && python3 run_tests.py particles/pic_paper_coupling_conservation_vl2_tsc --cmake=-DCMAKE_BUILD_TYPE=Debug`
- `cd tst && python3 run_tests.py particles/pic_mhd_coupling_decomp --cmake=-DAthena_ENABLE_MPI=ON`
- `python3 tst/scripts/particles/pic_paper_smooth_tsc_oracle.py --indent 0`
- `python3 -m flake8 tst/scripts/particles`

## Local constraints

- `*_proxy.py` files are regression or engineering anchors, not physical
  qualification. Legacy `*_publication.py` modules are excluded from broad
  discovery by default and do not establish publication evidence.
- `pic_paper_smooth_tsc_oracle.py` is a standalone oracle without harness
  `run()`/`analyze()` functions; the `*_oracle.py` suffix excludes it from broad
  suite discovery while the compiled interface regression imports it directly.
- PIC or shock preparation contracts may validate frozen source/deck
  structure without executing or qualifying a production campaign.
<!-- END build-memory-table -->
