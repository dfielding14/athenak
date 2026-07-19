<!-- BEGIN build-memory-table -->
# Source Terms and Forcing Navigation

## Scope

This subtree owns per-fluid source terms, general cooling/heating, one-time spectral
initial perturbations, the stateful turbulence driver, and the post-timestep frame
tracker. `SourceTerms` objects belong to hydro/MHD/radiation modules; turbulence and frame
tracking are pack-level task participants, while initial perturbations run once after a
new-run pgen.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Add a local fluid or radiation source | `src/srcterms/srcterms.*` | Fluid `*_tasks.cpp`, timestep impact in `srcterms_newdt.cpp` |
| Change general cooling/heating | `src/srcterms/cooling.*` | `docs/source/modules/cooling.md`, units, history output, cooling tests |
| Change turbulence forcing | `src/srcterms/turb_driver.*` | Pack task registration, restart serialization, turbulence docs/tests |
| Change frame tracking | `src/srcterms/frame_tracker.*` | `docs/source/modules/frame_tracking.md`, boundaries, history/restart, AMR rebuild |
| Change initial perturbations | `src/srcterms/initial_perturbations.*` | `src/main.cpp`, MHD discrete curl, example plots/tests |
| Add a user cooling model | Hooks in `src/pgen/pgen.hpp`, `cooling_hooks.hpp` | Host launcher, timestep hook, pgen restart enrollment |

## Important flow

- Hydro/MHD stage tasks pass primitive state, EOS data, explicit stage weight, and source
  history weight to `SourceTerms::ApplySrcTerms`, which updates conserved state. Source
  timestep restrictions feed the module's `NewTimeStep` reduction.
- General cooling is enabled by the standalone `<cooling>` block but executes through the
  relevant fluid's `SourceTerms` object and contributes optional history diagnostics.
- `TurbulenceDriver` evolves stochastic modal amplitudes before the integrator and applies
  forcing during stages; its modal/RNG state is serialized in restarts.
- `FrameTracker` samples and boosts nonrelativistic fluid state after a completed
  integrator, refreshes affected boundaries/timestep data, and stores controller state.
- `ApplyInitialPerturbations` runs only for fresh initial data. Magnetic perturbations are
  constructed as a discrete curl of an edge-centered vector potential.

## Local validation

- `cd tst && python3 run_test_suite.py --test test_suite/cooling/test_cooling_cpu.py`
- `cd tst && python3 run_test_suite.py --test test_suite/turb/test_turb_driving_cpu.py`
- `cd tst && python3 run_test_suite.py --test test_suite/nr/test_nr_frame_tracking_cpu.py`
- `cd tst && python3 run_test_suite.py --test test_suite/mhd/test_mhd_initial_perturbations_cpu.py`

Each command rebuilds a focused CPU test target through the standard harness.

## Local constraints

- Compute explicit fluid sources from primitive state and apply increments to conserved
  state with the stage weight supplied by `Driver`.
- A source that can be more restrictive than the hyperbolic CFL limit must update
  `dtnew` and be included in the owning fluid's timestep reduction.
- Stateful pack-level features require coordinated task registration, restart state,
  AMR pack-pointer refresh, history output, documentation, and restart tests.
- Use the existing module documentation as the parameter contract; update it and the
  shipped inputs when accepted keys or semantics change.
<!-- END build-memory-table -->
