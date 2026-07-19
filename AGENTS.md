<!-- BEGIN build-memory-table -->
# AthenaK Repository Navigation

## Start here

AthenaK is a C++17/Kokkos block-AMR application for performance-portable astrophysical
fluid, radiation, particle, and numerical-relativity simulations. `README.md` describes
the supported physics, the root `CMakeLists.txt` is the build authority, `src/main.cpp`
owns the executable lifecycle, and `tst/` contains the regression harness.

`kokkos/` is a Git submodule and a separate repository. Treat it as an external
dependency: do not place parent-repository guidance there or mix Kokkos edits into an
AthenaK change.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Understand runtime architecture or change physics code | `src/AGENTS.md` | Nearest subsystem guide and focused regression |
| Change mesh topology, packs, load balancing, or AMR | `src/mesh/AGENTS.md` | Boundary, output/restart, and physics migration paths |
| Change ghost exchange or physical boundaries | `src/bvals/AGENTS.md` | Mesh neighbor metadata and module task ordering |
| Change EOS or primitive recovery | `src/eos/AGENTS.md` | Hydro/MHD constructors and SR/GR tests |
| Change output variables, formats, I/O, or restarts | `src/outputs/AGENTS.md` | `vis/python/`, consuming scripts, restart tests |
| Add initial data or user callbacks | `src/pgen/AGENTS.md` | Matching `inputs/` example and regression test |
| Change cooling, forcing, perturbations, or frame tracking | `src/srcterms/AGENTS.md` | Feature documentation, history/restart, focused tests |
| Change Z4c spacetime or NR diagnostics | `src/z4c/AGENTS.md` | Dynamical matter, NR task assembly, NR/Z4c tests |
| Add, run, or debug regression tests | `tst/AGENTS.md` | CI device class and generated output readers |
| Change build options or compiled sources | `CMakeLists.txt`, `src/CMakeLists.txt` | `config.hpp.in`, CI configurations |
| Change a shipped simulation setup | `inputs/` | Corresponding pgen, module parameter parser, documentation/test input |
| Change user or feature documentation | `docs/source/` | Implemented input keys and shipped examples |
| Change post-processing or visualization | `scripts/`, `vis/python/` | Producing output writer and a representative data file |
| Change cooling-table validation tooling | `tools/` | Cooling module and its documented table contract |

## Architecture and execution flow

1. CMake configures Kokkos and optional MPI/external dependencies, then
   `src/CMakeLists.txt` assembles the `athena` executable.
2. `src/main.cpp` initializes MPI and Kokkos, reads an input or restart, constructs the
   mesh and rank-local block pack, attaches input-selected physics, and runs the problem
   generator.
3. The driver executes dependency-ordered task lists for time integration, boundaries,
   sources, coupled physics, AMR, and diagnostics. Module presence and input blocks shape
   the task graph.
4. Outputs stage device data to host writers. Restart persistence also crosses mesh-tree
   construction and problem-generator restoration, so it is not isolated to one writer.

Follow state across module boundaries. A new evolved quantity can require coordinated
changes to allocation and indices, kernels, boundary buffers, AMR migration, outputs,
restart serialization, and tests.

## Build and run

From the repository root:

```sh
git submodule update --init --recursive
cmake -S . -B build
cmake --build build --parallel
build/src/athena -i inputs/hydro/sod.athinput
```

Enable MPI with `-DAthena_ENABLE_MPI=ON`. Custom problem generators are selected at
configure time with `-DPROBLEM=<name>`; see `src/pgen/AGENTS.md`. Keep builds out of the
source tree.

## Repository-wide validation

Run the current harness from `tst/`:

```sh
cd tst
python3 run_test_suite.py --style
python3 run_test_suite.py --test test_suite/<area>/<test>_cpu.py
```

Use `--cpu`, `--mpicpu`, or `--gpu` for the relevant full suite. Read `tst/AGENTS.md`
before selecting MPI/GPU flags or maintaining a legacy test. Prefer a focused numerical
regression that checks errors, convergence, conservation, or another physical invariant
over a run-only smoke check.

## Repository conventions

- Follow the Google C++ style used by the project and retain C++17/Kokkos portability.
- Add every ordinary compiled source to `src/CMakeLists.txt`; custom pgens are the noted
  configure-time exception.
- Keep input keys, examples, feature documentation, and regression tests synchronized.
- Preserve explicit host/device synchronization and keep device-kernel captures
  device-compatible.
- Add a regression for behavior changes and use the nearest existing test as the harness
  and tolerance model.

## Local guides

- `src/AGENTS.md`: runtime architecture and physics routing.
- `src/mesh/AGENTS.md`: topology, packs, AMR, and load balancing.
- `src/bvals/AGENTS.md`: boundary exchange, flux correction, and physical BCs.
- `src/eos/AGENTS.md`: EOS selection, floors, and primitive recovery.
- `src/outputs/AGENTS.md`: output formats, variables, I/O, and restarts.
- `src/pgen/AGENTS.md`: initial data, built-in/custom problems, and callbacks.
- `src/srcterms/AGENTS.md`: cooling, forcing, perturbations, and frame tracking.
- `src/z4c/AGENTS.md`: spacetime evolution and numerical-relativity diagnostics.
- `tst/AGENTS.md`: current and legacy regression workflows.
<!-- END build-memory-table -->
