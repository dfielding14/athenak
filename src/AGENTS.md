# AGENTS.md

## Purpose
This directory contains the core C++ source for AthenaK: the main entry point,
shared type definitions, input parsing, and all physics/mesh modules. Most
subdirectories have their own `AGENTS.md` with detailed behavior; use those for
module-specific edits. See `../AGENTS.md` for repository-wide conventions and
workflow.

---

## Root-Level Files

### Build and entry point
- `src/CMakeLists.txt`: CMake rules for building the core library/executable.
- `src/main.cpp`: program entry point (MPI/Kokkos init, input parsing, `Mesh`
  construction, driver execution).

### Core types and globals
- `src/athena.hpp`: core typedefs (`Real`), enums for variable indices, and
  shared macros.
- `src/athena_tensor.hpp`: tensor helpers, symmetry tags, and Kokkos subview
  aliases for tensor-like fields.
- `src/parameter_input.hpp` / `src/parameter_input.cpp`: `ParameterInput`
  implementation for `.athinput` files.
- `src/globals.hpp` / `src/globals.cpp`: `global_variable` namespace containing
  the MPI rank and size.

---

## Subdirectory Map
Each subdirectory below contains its own `AGENTS.md` with authoritative
details. Use those docs before making changes.

- `src/bvals/`: boundary conditions and MPI exchanges.
  See `src/bvals/AGENTS.md`.
- `src/coordinates/`: coordinate systems and metric helpers.
  See `src/coordinates/AGENTS.md`.
- `src/diffusion/`: diffusion and conduction operators.
  See `src/diffusion/AGENTS.md`.
- `src/driver/`: time integration drivers and scheduling helpers.
  See `src/driver/AGENTS.md`.
- `src/dyn_grmhd/`: dynamical GRMHD coupling for evolving spacetime with MHD.
  See `src/dyn_grmhd/AGENTS.md`.
- `src/eos/`: equations of state.
  See `src/eos/AGENTS.md`.
- `src/geodesic-grid/`: geodesic angular meshes and spherical interpolation.
  See `src/geodesic-grid/AGENTS.md`.
- `src/hydro/`: hydrodynamics module.
  See `src/hydro/AGENTS.md`.
- `src/ion-neutral/`: ion-neutral (two-fluid) coupling.
  See `src/ion-neutral/AGENTS.md`.
- `src/mesh/`: mesh, MeshBlock, AMR/SMR, and MeshBlockPack.
  See `src/mesh/AGENTS.md`.
- `src/mhd/`: magnetohydrodynamics module.
  See `src/mhd/AGENTS.md`.
- `src/outputs/`: output scheduling and file writers.
  See `src/outputs/AGENTS.md`.
- `src/particles/`: particle infrastructure and pushers.
  See `src/particles/AGENTS.md`.
- `src/pgen/`: problem generators and setup utilities.
  See `src/pgen/AGENTS.md`.
- `src/radiation/`: radiation transport module.
  See `src/radiation/AGENTS.md`.
- `src/reconstruct/`: reconstruction algorithms and limiters.
  See `src/reconstruct/AGENTS.md`.
- `src/shearing_box/`: shearing-box orbital advection and boundaries.
  See `src/shearing_box/AGENTS.md`.
- `src/srcterms/`: source term infrastructure.
  See `src/srcterms/AGENTS.md`.
- `src/tasklist/`: task list definitions and orchestration.
  See `src/tasklist/AGENTS.md`.
- `src/units/`: code-unit and CGS conversions.
  See `src/units/AGENTS.md`.
- `src/utils/`: shared utilities (finite differences, interpolation, runtime
  helpers, etc.).
  See `src/utils/AGENTS.md`.
- `src/z4c/`: Z4c numerical relativity module.
  See `src/z4c/AGENTS.md`.
