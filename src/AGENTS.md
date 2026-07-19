<!-- BEGIN build-memory-table -->
# Runtime Source Navigation

## Scope

`src/` builds the `athena` executable and contains the runtime, mesh, physics, numerical
methods, problem setup, and output implementation. `src/main.cpp` is the lifecycle entry
point; `src/mesh/meshblock_pack.cpp` is the module-registration choke point; and
`src/driver/driver.cpp` executes the dependency-ordered evolution tasks.

## Task lookup

| Task | Start in | Also inspect or validate |
| --- | --- | --- |
| Change global types, variable indices, or Kokkos loop helpers | `src/athena.hpp`, `src/athena_tensor.hpp` | Every kernel/state layout consumer |
| Change CLI, input parsing, startup, or restart entry | `src/main.cpp`, `src/parameter_input.*` | Mesh build, pgen restart, outputs |
| Change time integration or generic task execution | `src/driver/`, `src/tasklist/task_list.hpp` | Module `*_tasks.cpp`, timestep reductions |
| Change mesh, AMR, or load balancing | `src/mesh/AGENTS.md` | Boundary and restart propagation |
| Change ghost exchange or physical boundaries | `src/bvals/AGENTS.md` | Mesh neighbors and physics task ordering |
| Change hydro flux/update behavior | `src/hydro/hydro.hpp`, `hydro_fluxes.cpp`, `hydro_update.cpp` | `hydro_tasks.cpp`, `rsolvers/`, `src/reconstruct/`, EOS |
| Change MHD or constrained transport | `src/mhd/mhd.hpp`, `mhd_fluxes.cpp`, `mhd_corner_e.cpp`, `mhd_ct.cpp` | `mhd_tasks.cpp`, face-boundary/AMR paths, EOS |
| Change radiation transport/coupling | `src/radiation/radiation.hpp`, `radiation_tasks.cpp` | Flux, source, tetrad, fluid coupling, radiation tests |
| Change EOS or primitive recovery | `src/eos/AGENTS.md` | Fluid constructors, FOFC, timestep and source consumers |
| Change dynamical-spacetime matter | `src/dyn_grmhd/` | `src/z4c/AGENTS.md`, `src/tasklist/numerical_relativity.*`, ADM/Tmunu |
| Change Z4c spacetime evolution | `src/z4c/AGENTS.md` | Numerical-relativity tasks, pgens and NR tests |
| Change particles or tracer fields | `src/particles/` | Particle bvals, pack registration, particle outputs/restarts |
| Change cooling, forcing, or moving-frame behavior | `src/srcterms/AGENTS.md` | Fluid tasks, restart/history, feature docs/tests |
| Change initial conditions or user callbacks | `src/pgen/AGENTS.md` | Matching input and regression test |
| Change outputs, derived data, or restart format | `src/outputs/AGENTS.md` | Readers, restart consumers, sharding tests |
| Change diffusion, shearing-box, or two-fluid coupling | `src/diffusion/`, `src/shearing_box/`, `src/ion-neutral/` | Owning Hydro/MHD task chain and focused tests |
| Change coordinate/geometry or utility interpolation | `src/coordinates/`, `src/geodesic-grid/`, `src/utils/`, `src/units/` | Physics kernels and outputs that consume the data |

## Architecture and flow

1. `main.cpp` initializes MPI before Kokkos, parses input/restart data, builds `Mesh` and
   its rank-local `MeshBlockPack`, and attaches coordinates and physics.
2. `MeshBlockPack::AddPhysics` uses input-block presence to construct modules and enroll
   their tasks. A pack holds the shared mesh/coordinate pointers plus active Hydro, MHD,
   radiation, particle, source-term, and numerical-relativity objects.
3. `ProblemGenerator` initializes a fresh state or restores restart payloads and reenrolls
   callbacks. `Driver` then initializes boundaries/primitives, loops over Runge-Kutta or
   ImEx stages, performs outputs/AMR, and runs final analysis.
4. Normal fluid stage flow is reconstruction/Riemann fluxes, AMR flux correction, RK
   update and sources, conserved-state exchange, physical BCs/prolongation, primitive
   conversion, and new-timestep calculation. MHD inserts edge-electric-field exchange and
   constrained transport; coupled modules assemble a combined task graph.

Device field arrays conventionally use `(meshblock, variable, k, j, i)` with face- and
edge-centered wrappers for staggered fields. Host decisions that feed kernels generally
use `DualView`; preserve explicit host/device synchronization.

## Local validation

- `cd tst && python3 run_test_suite.py --style`: source style and Python lint.
- `cd tst && python3 run_test_suite.py --test test_suite/<area>/<test>_cpu.py`: focused
  CPU regression; `_mpicpu.py` and `_gpu.py` select those build modes.
- Read `tst/AGENTS.md` before adding or changing tests. Use a physics-specific convergence
  or invariant test, not only a successful executable exit.

## Local constraints

- Add every new compiled `.cpp` to `src/CMakeLists.txt`, except the single custom pgen
  selected dynamically by `-DPROBLEM=<name>`.
- Input keys are part of the runtime interface. Update shipped inputs, tests, and feature
  documentation when adding or changing accepted parameters.
- A new evolved field commonly propagates through allocation, tasks, boundary buffers,
  AMR migration, output selection, and restart serialization; audit all of those paths.
- Keep code inside Kokkos kernels device-compatible and avoid host-only state or
  unsynchronized views in captures.

## Local guides

- `src/mesh/AGENTS.md`: topology, packs, AMR, and load balancing.
- `src/bvals/AGENTS.md`: cell/face boundary exchange, flux correction, and physical BCs.
- `src/eos/AGENTS.md`: EOS selection, floors, and primitive recovery.
- `src/outputs/AGENTS.md`: output formats, variables, I/O, and restarts.
- `src/pgen/AGENTS.md`: initial data, built-in/custom problems, and callbacks.
- `src/srcterms/AGENTS.md`: cooling, forcing, perturbations, and frame tracking.
- `src/z4c/AGENTS.md`: spacetime evolution and numerical-relativity diagnostics.
<!-- END build-memory-table -->
