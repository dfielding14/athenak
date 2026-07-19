<!-- BEGIN build-memory-table -->
# Problem Generator Navigation

## Scope

This subtree owns initial conditions, built-in regression problems, custom problem files,
and optional callbacks for user boundaries, source terms, timestep limits, refinement,
history, cooling, and final analysis. It also restores physics payloads on restart before
reenrolling problem-specific behavior.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change dispatch or the pgen callback interface | `src/pgen/pgen.hpp`, `src/pgen/pgen.cpp` | Fresh-run and restart constructor branches |
| Add a built-in regression problem | `src/pgen/tests/` | Declarations and both dispatch tables in `pgen.cpp`, `src/CMakeLists.txt`, `tst/` |
| Add or change a custom compiled problem | `src/pgen/<problem>.cpp` | Configure with `-DPROBLEM=<problem>`, matching `inputs/` example |
| Set hydro/MHD/radiation/Z4c initial state | Closest existing pgen for that physics | EOS conversion, face-field staggering, boundary initialization |
| Enroll a user boundary/source/timestep/refinement/history hook | Function pointers in `pgen.hpp` | Call sites in bvals, fluid tasks, mesh refinement, outputs |
| Change restart payload restoration | Restart constructor in `src/pgen/pgen.cpp` | `src/outputs/restart.cpp`, mesh restart metadata, stateful modules |
| Add final error or analysis output | `pgen_final_func`, `OutputErrors` | Driver finalization and regression assertions |

## Important flow

- With the default build, `<problem>/pgen_name` selects a built-in generator compiled from
  `src/pgen/tests/`. With `-DPROBLEM=name`, CMake compiles
  `src/pgen/name.cpp`, defines `USER_PROBLEM_ENABLED`, and calls that file's
  `ProblemGenerator::UserProblem` instead.
- New runs allocate mesh and physics first, run the selected generator, then apply shared
  initial perturbations. Restart runs rebuild topology, restore each active module's
  arrays, and call the selected generator again with `restart=true` to reenroll callbacks
  and nonserialized state.
- Inputs define the physics and numerical methods independently of the pgen; the generator
  must validate that the modules and variable layouts it needs are present.

## Local validation

- `cd tst && python3 run_test_suite.py --test test_suite/sr/test_sr_lwave1d_cpu.py`:
  representative built-in linear-wave generator and final error analysis.
- `cd tst && python3 run_test_suite.py --test test_suite/nr/test_nr_frame_tracking_examples_cpu.py`:
  custom pgens with user boundaries and moving-frame state.
- For a new custom pgen, configure it explicitly, for example
  `cmake -S . -B build -DPROBLEM=turb`, then build and run its matching input.

## Local constraints

- Every top-level custom pgen defines the same `UserProblem` symbol, so only the file
  selected by CMake may be compiled. Built-in test pgens instead have distinct member names.
- Add a built-in in both fresh-run and restart dispatch branches; omitting the restart
  branch makes its restart files unusable.
- Enroll required callbacks before an early `if (restart) return`; the restart constructor
  verifies that requested user callbacks were restored.
- For MHD initialization, preserve face-centered magnetic staggering and perform the
  matching primitive/conserved initialization used by a nearby pgen.
<!-- END build-memory-table -->
