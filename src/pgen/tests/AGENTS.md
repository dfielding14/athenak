<!-- BEGIN build-memory-table -->
# Built-in Problem Generator Navigation

## Scope

This directory contains problem generators linked into AthenaK's built-in
problem executable, plus a small number of compile-time test fixtures. Built-ins
are selected by `<problem>/pgen_name`; declarations, fresh/restart dispatch, and
build registration live in `src/pgen/pgen.hpp`, `src/pgen/pgen.cpp`, and
`src/CMakeLists.txt`. `driver_user_stop.cpp` is a compile-time `UserProblem`
fixture selected with `-DPROBLEM=tests/driver_user_stop` instead. Active MHD-PIC
priorities are in `MHD_PIC_NEXT_STEPS_GUIDE.md`; durable model and shock setup
documents live under `docs/source/engineering/`.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Add or rename a built-in generator | Its implementation file | Declaration in `src/pgen/pgen.hpp`, both dispatch chains in `src/pgen/pgen.cpp`, `src/CMakeLists.txt`, a matching input, and its regression |
| Change a classical hydro/MHD test | `advection.cpp`, `cpaw.cpp`, `divb_amr.cpp`, `linear_wave.cpp`, `orszag_tang.cpp`, or `shock_tube.cpp` | Matching suite under `tst/scripts/` and `inputs/tests/` |
| Change PIC interface or distribution setup | `pic_paper_smooth_tsc_interface.cpp`, `q006_*`, `q007_*`, `q032_*`, or `q033_*` | `src/particles/`, paired input decks, and `tst/scripts/particles/` |
| Change parallel-shock injection, split removal, tracer sampling, or accounting | `pic_parallel_shock.cpp` | `inputs/tests/pic_parallel_shock_*`, publication shock decks including the Mignone-R2 freeze, and matching particle/Q011 tests |
| Change linear Bell seeds or current normalization | `q023_paper_bell_linear*`, `q043_bell_current_volume_aware.cpp` | Shared headers and corresponding host harnesses under `tst/publication/` |
| Change nonlinear Bell setup | `q019_nonlinear_bell_saturation_engineering.cpp`, `q019_physics_first_nonlinear_bell_successor_v2.hpp` | `inputs/publication/pic_bell_nonlinear_saturation_engineering_v1.athinput`, the full-Hall nonlinear pilot, and the compact Bell analyzer/tests |
| Change orderly driver-stop behavior | `driver_user_stop.cpp` | `inputs/tests/driver_user_stop.athinput` and `tst/publication/test_driver_user_stop.py` |

## Important flow

- A fresh built-in run dispatches on the exact `pgen_name`. A restart restores
  mesh, physics, and particle state before calling the same generator with
  `restart=true`. Re-enroll callbacks and rebuild only derived reservoirs before
  the restart early return; never overwrite restored state afterward.
- MHD Bell initializers build staggered magnetic state through a vector potential
  and discrete curl, derive cell-centered fields, initialize primitives, and
  convert to conserved variables. Preserve high-side face coverage and the
  divergence-free construction when changing those generators.
- `pic_parallel_shock.cpp` prepares one particle-injection transaction per
  physical cycle, replays its gas subtraction with Runge-Kutta stage weights,
  then validates and commits conservation, provenance, and escape ledgers.
  Reordering its callbacks can double-inject particles or corrupt accounting.
  Before preparing that transaction, the callback limits the current timestep
  using the configured injection-velocity envelope; the injected mass budget
  and first particle push must therefore consume the same bounded `dt`.
- The Q019 successor restores its runtime monitor state and re-enrolls its
  work-in-loop, history, and final callbacks before returning on restart.
- Shared Bell headers feed several generators. Trace all consumers before
  changing carrier math or normalization.
- MHD-PIC generators select staged coupling explicitly with
  `<time>/integrator=vl2` and the complete TSC/coupling controls. Boris CR state
  is always `p/m`; `pic_cr_initial_state` controls only initializer input
  interpretation. Hall choices are only `pic_cr_hall_mode=off|full`.

## Focused validation

From the repository root, lightweight Bell/source contracts include:

- `python3 -m pytest -q tst/publication/test_q023_paper_bell_linear_host_harness.py`
- `python3 -m pytest -q tst/publication/test_q043_bell_current_volume_aware_host_harness.py`
- `python3 -m pytest -q tst/publication/test_bell_saturation_engineering_analysis_v1.py`

Runtime shock changes normally require a focused compiled regression, for example:

- `cd tst && python3 run_tests.py particles/pic_parallel_shock_rk_stage_budget_vl2_tsc`
- `cd tst && python3 run_tests.py particles/pic_parallel_shock_restart_controls`
- `cd tst && python3 run_tests.py particles/pic_parallel_shock_outer_x1_escape_restart`
- `cd tst && python3 run_tests.py particles/pic_parallel_shock_split_removal_restart`
- `cd tst && python3 run_tests.py particles/pic_parallel_shock_mignone_tracer_smoke`

## Local constraints

- Treat a PIC/Bell generator, its deck, host math, analyzer, and registration as
  one contract. Many paths intentionally fail closed on geometry, normalization,
  runtime mode, and campaign metadata.
- Keep corrected current-normalization, Hall-off, and full-Hall Bell contracts
  distinct. Sharing seed-carrier math does not make their physical claims
  interchangeable.
- Preparation, proxy, candidate, or engineering labels do not establish
  qualification or authorize a production campaign.
<!-- END build-memory-table -->
