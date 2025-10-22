# AthenaK Codebase Overview

## Purpose and High-Level Architecture
AthenaK is a performance-portable astrophysical simulation framework that targets multi-physics magnetohydrodynamics, radiation transport, and numerical relativity on block-structured adaptive meshes. The code is written in modern C++ and is designed to run efficiently across CPUs and GPUs through the [Kokkos](https://kokkos.org) programming model. The simulation lifecycle follows a consistent pattern:

1. **Configuration and Input Parsing** – `src/main.cpp` bootstraps MPI/Kokkos, parses command-line arguments, and loads runtime parameters through `ParameterInput`.【F:src/main.cpp†L23-L140】【F:src/parameter_input.hpp†L14-L117】
2. **Mesh Construction** – A `Mesh` object owns the domain decomposition, adaptive mesh refinement (AMR) hierarchy, boundary conditions, and mesh block metadata. It instantiates `MeshBlockPack` containers that aggregate `MeshBlock` patches for GPU-friendly execution and register the active physics modules.【F:src/mesh/mesh.hpp†L61-L170】【F:src/mesh/meshblock_pack.hpp†L1-L99】
3. **Driver Execution** – The `Driver` orchestrates time integration, choosing explicit or implicit-explicit Runge–Kutta schemes, coordinating task lists, and invoking physics packages for updates. It also tracks diagnostics and wall-clock limits.【F:src/driver/driver.cpp†L1-L120】【F:src/driver/driver.hpp†L1-L73】
4. **Output and Finalization** – `Outputs` manages periodic diagnostics, restart dumps, and mesh/particle visualization products, while the driver handles teardown of physics modules and MPI/Kokkos resources.【F:src/outputs/outputs.hpp†L1-L104】

Key compile-time types, enums, and memory abstractions live in `athena.hpp`, providing aliases for Kokkos array layouts, physics variable indices, and numerical settings shared across the codebase.【F:src/athena.hpp†L1-L129】 Global MPI rank metadata is centralized in `globals.cpp`/`globals.hpp` for lightweight access by subsystems.【F:src/globals.cpp†L1-L18】

## Source Tree Organization

### Core Infrastructure (`src` root)
- **`main.cpp`** – Entry point that initializes MPI/Kokkos, reads input files, configures logging/output paths, and dispatches the main simulation driver.【F:src/main.cpp†L1-L140】
- **`athena.hpp` / `athena_tensor.hpp`** – Foundational typedefs, array aliases, enums, and tensor utilities used across the physics modules.【F:src/athena.hpp†L12-L118】
- **`globals.*`** – Stores process-wide MPI rank counts, populated during startup.【F:src/globals.cpp†L1-L18】
- **`parameter_input.*`** – Parses text input decks into searchable blocks/parameters with thread-safe accessors for integers, reals, booleans, and strings.【F:src/parameter_input.hpp†L30-L117】
- **`outputs/`** – Implements the flexible output system, supporting history tables, binary/coarsened datasets, restarts, PDF diagnostics, particle tracks, and VTK writers. The module exposes `OutputParameters` and `Outputs` classes for scheduling and writing datasets.【F:src/outputs/outputs.hpp†L1-L109】
- **`utils/`** – Miscellaneous utilities: runtime configuration dumps, run directory management, random-number helpers, Lagrange interpolation tables, and finite-difference kernels.【F:src/utils/utils.hpp†L1-L15】

### Mesh and Execution Management (`src/mesh`)
- **`Mesh`** – Reads mesh sizing, boundary, and refinement settings; builds AMR trees; assigns mesh block packs; and manages global counters for mesh blocks and particles.【F:src/mesh/mesh.cpp†L1-L120】
- **`MeshBlock`** – Lightweight holder for per-block geometry, neighbor metadata, and ghost zone layout shared by all physics packages.【F:src/mesh/meshblock.hpp†L1-L60】
- **`MeshBlockPack`** – Aggregates multiple `MeshBlock`s, instantiates coordinate systems, and lazily constructs physics objects (hydro, MHD, relativity, radiation, particles, turbulence, etc.) so that task lists can operate collectively on GPU-friendly batches.【F:src/mesh/meshblock_pack.hpp†L22-L94】
- **`MeshRefinement` / `MeshBlockTree`** – Support static and adaptive refinement by maintaining logical tree structures, neighbor discovery, and prolongation/restriction operators (files under `mesh_refinement*` and `meshblock_tree*`).

### Physics Modules
AthenaK enables multiple physics disciplines through modular namespaces that attach to mesh block packs:
- **Hydrodynamics (`src/hydro`)** – Conservative finite-volume solvers, Riemann solvers (`rsolvers/`), flux computations, source-term handling, and timestep control for Newtonian flows.【F:src/hydro/hydro.hpp†L1-L20】【F:src/hydro/rsolvers/hlle_hyd.hpp†L1-L20】
- **Magnetohydrodynamics (`src/mhd`)** – Extends hydro with constrained transport, corner electric fields, divergence control, and MHD-specific Riemann solvers.【F:src/mhd/mhd.hpp†L1-L18】
- **Dynamical GRMHD (`src/dyn_grmhd`)** – Specializes relativistic MHD evolution in dynamical spacetimes, including dedicated flux, force, and utility routines.【F:src/dyn_grmhd/dyn_grmhd.hpp†L1-L19】
- **Numerical Relativity (`src/z4c`, `src/tasklist/numerical_relativity.*`)** – Implements the Z4c formulation, ADM metric evolution, gravitational waveform extraction, and coupling to matter via stress-energy tensors.【F:src/z4c/z4c.hpp†L1-L30】【F:src/tasklist/numerical_relativity.hpp†L1-L16】
- **Radiation Transport (`src/radiation`)** – Handles moment-based radiation evolution, source coupling, tetrad transforms, opacity models, and dedicated task lists.【F:src/radiation/radiation.hpp†L1-L22】
- **Source Terms (`src/srcterms`)** – Centralized interface for explicit and implicit body forces, tabulated cooling, and turbulence drivers used by multiple physics packages.【F:src/srcterms/srcterms.hpp†L1-L23】
- **Diffusion (`src/diffusion`)** – Optional viscosity, resistivity, and conduction operators with separate coefficient calculations.【F:src/diffusion/conduction.hpp†L1-L15】
- **Ion-Neutral Coupling (`src/ion-neutral`)** – Evolves ion-neutral drift physics and registers related task stages.【F:src/ion-neutral/ion-neutral.hpp†L1-L18】
- **Particles (`src/particles`)** – Lagrangian particle infrastructure, including data layouts, pushers, and mesh-based interaction tasks for tracers, cosmic rays, and star particles.【F:src/particles/particles.hpp†L1-L26】
- **Units (`src/units`)** – Converts between simulation code units and physical (cgs) units to support cooling tables and diagnostics.【F:src/units/units.hpp†L1-L17】

Physics activation is dynamic: `MeshBlockPack::AddPhysics` inspects the input deck and constructs only requested modules, enabling hybrid simulations that mix hydrodynamics, relativity, radiation, and particle dynamics.【F:src/mesh/meshblock_pack.hpp†L57-L94】

### Coordinate Systems and Geometry (`src/coordinates`)
Coordinate classes provide metric tensors, Christoffel symbols, face areas, and Jacobians for different geometries (Cartesian, spherical, Kerr–Schild, etc.). The ADM helper handles background metric data for relativity runs.【F:src/coordinates/coordinates.hpp†L1-L24】 Excision utilities manage interior boundary conditions near compact objects.

### Reconstruction and Riemann Solvers (`src/reconstruct`, `src/*/rsolvers`)
High-order reconstructions (donor-cell, PLM, PPM, WENO-Z) are implemented in `reconstruct/`. Each fluid module ships tailored Riemann solvers for classical, relativistic, and degenerate regimes within its `rsolvers/` subdirectory.【F:src/reconstruct/plm.hpp†L1-L20】【F:src/hydro/rsolvers/roe_hyd.hpp†L1-L18】

### Task Lists (`src/tasklist`)
Task lists break complex updates into dependency-driven stages that can be scheduled flexibly across Kokkos execution spaces. The numerical relativity task list is a representative example that coordinates ADM/Z4c operations per timestep.【F:src/tasklist/task_list.hpp†L1-L68】【F:src/tasklist/numerical_relativity.cpp†L1-L28】

### Source Term Infrastructure (`src/srcterms`)
The source term framework supplies base classes and utilities for turbulence forcing (`turb_driver.*`), tabulated cooling (`cooling_tables.hpp`), and specialized ISM cooling models. It integrates with hydro and MHD updates via the driver’s implicit-explicit staging.【F:src/srcterms/srcterms.hpp†L1-L23】

### Problem Generators (`src/pgen`)
A rich catalogue of problem setups (blast waves, Kelvin–Helmholtz, turbulence boxes, binary neutron stars, radiation tests, etc.) lives under `pgen/`. Each generator populates initial conditions, registers user-defined source terms, and optionally supplies custom analysis hooks executed at initialization/finalization through the `Driver`. The directory also contains regression-focused scenarios within `pgen/tests/` for CI runs.【F:src/pgen/pgen.hpp†L1-L25】

### Additional Subsystems
- **`shearing_box/`** – Implements shearing-box boundary management for differentially rotating systems (coupled to mesh boundary flags).【F:src/shearing_box/shearing_box.hpp†L1-L18】
- **`outputs/pdf.cpp`, `outputs/track_prtcl.cpp`** – Provide specialized diagnostics like probability density functions and particle trajectory dumps.
- **`utils/sn_scheduler.hpp`** – Houses supernova event scheduling utilities used by star formation problems.
- **`srcterms/ismcooling.hpp`** – Adds ISM cooling microphysics interfacing with tabulated rates.

## Testing and Tooling
- **Regression Tests (`tst/`)** – `run_tests.py` dynamically discovers regression suites in `tst/scripts`, builds AthenaK via CMake, executes each problem generator, and analyzes results. Failures raise a `TestError` that CI can capture.【F:tst/run_tests.py†L1-L119】
- **Build Helpers (`configure.sh`, `scripts/`)** – Shell and Python helpers configure Kokkos backends, set compiler flags, and manage test automation.
- **Visualization (`vis/`)** – Contains plotting and analysis utilities (not inspected per instruction but part of delivered tree) for post-processing outputs.

## Execution Flow Summary
1. Users compile via CMake or helper scripts, linking against Kokkos and optional MPI/OpenMP backends.
2. At runtime, `main.cpp` initializes MPI/Kokkos, reads the input deck, sets up restart directories, and constructs a `Mesh` with corresponding `MeshBlockPack`s and physics modules.【F:src/main.cpp†L1-L140】【F:src/mesh/mesh.cpp†L1-L120】
3. The `Driver` selects the appropriate task list(s) (hydro, MHD, GR, radiation, particles) based on enabled modules, iterating through time integration stages, invoking reconstruction, flux computations, source term updates, and mesh refinement as needed.【F:src/driver/driver.cpp†L1-L120】【F:src/mesh/meshblock_pack.hpp†L57-L94】
4. Outputs and diagnostics trigger according to `OutputParameters`, writing history tables, checkpoints, mesh dumps, and particle data for downstream analysis.【F:src/outputs/outputs.hpp†L1-L109】
5. Upon completion or wall-clock expiry, the driver finalizes modules, writes final outputs, and `main.cpp` shuts down Kokkos and MPI, ensuring clean teardown even in parallel environments.【F:src/main.cpp†L129-L188】

## Key Design Considerations for Contributors
- **Performance Portability** – All heavy compute kernels operate on Kokkos `View`/`DualView` types; new kernels should follow the established memory layout conventions defined in `athena.hpp`.【F:src/athena.hpp†L88-L125】
- **Task-Oriented Execution** – Updates are staged via task lists to expose concurrency and maintain explicit dependencies. When adding new physics, supply task entries that plug into the driver’s execution loops rather than hard-coding update sequences.【F:src/tasklist/task_list.hpp†L1-L68】
- **Modular Physics Activation** – New physics packages should mirror existing modules by providing constructor hooks in `MeshBlockPack::AddPhysics` and using `ParameterInput` flags for configuration toggles.【F:src/mesh/meshblock_pack.hpp†L57-L94】
- **Robust I/O** – Extend the output system by adding new identifiers to `var_choice` and wiring data extraction in the appropriate output class, maintaining compatibility with restart and analysis workflows.【F:src/outputs/outputs.hpp†L17-L86】
- **Testing** – Integrate new regression tests under `tst/scripts` to ensure compatibility across supported physics combinations, leveraging the shared harness in `run_tests.py`.【F:tst/run_tests.py†L1-L119】

This overview should equip new contributors with the mental model needed to navigate AthenaK’s core architecture, understand the relationships between mesh management, physics modules, and drivers, and extend the codebase while honoring its performance-portable design.
