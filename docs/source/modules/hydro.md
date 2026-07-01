# Module: Hydrodynamics

## Role in AthenaK
The Hydrodynamics module advances the Euler equations for Newtonian, special-relativistic, and general-relativistic fluids. It owns the mesh-pack level state (`Hydro` in `src/hydro/hydro.{hpp,cpp}`), assembles the task lists that integrate the equations each Runge–Kutta stage (`src/hydro/hydro_tasks.cpp`), and coordinates optional diffusion, source-term, and shearing-box couplings.

## File Layout

| File | Purpose | Key Routines |
|------|---------|--------------|
| `hydro.hpp/cpp` | `Hydro` class, construction, data ownership | constructor, `AssembleHydroTasks`, accessors |
| `hydro_tasks.cpp` | Task graph used by the integrator | `AssembleHydroTasks`, task wrappers (`Fluxes`, `RKUpdate`, …) |
| `hydro_fluxes.cpp` | Template flux kernel per Riemann solver | `CalculateFluxes<Hydro_RSolver::*>` |
| `hydro_newdt.cpp` | CFL and source-term timestep limiter | `Hydro::NewTimeStep` |

## Runtime Overview
`MeshBlockPack::AddPhysics` builds a single `Hydro` object per pack. During construction the class:

- Selects the equation of state and instantiates the appropriate EOS helper (`src/hydro/hydro.cpp:40`).
- Allocates conservative/primitive arrays (and coarse-grid buffers when AMR is enabled).
- Optionally enables viscosity (`diffusion/viscosity.hpp`), hyperviscosity (`diffusion/hyperviscosity.hpp`), and thermal conduction (`diffusion/conduction.hpp`) if the corresponding coefficients appear in the `<hydro>` block.
- Always creates a `SourceTerms` instance for hydro sources.
- Registers all per-stage tasks via `Hydro::AssembleHydroTasks`, wiring MPI exchanges, flux calculations, updates, prolongation/restriction, and the cooling timestep evaluation (`src/hydro/hydro_tasks.cpp:36`).

### Stage Task Order
For each integrator stage (`tl["stagen"]`) the tasks execute in the following sequence:

1. `CopyCons` – prepare stage registers (`u1`) for multi-stage schemes.
2. `Fluxes` – reconstruct interface states, solve Riemann problems, and add diffusive fluxes.
3. `SendFlux`/`RecvFlux` – exchange flux buffers when AMR is active.
4. `RKUpdate` – apply Runge–Kutta increments to `u0`.
5. `HydroSrcTerms` – inject source terms (links to `SourceTerms`).
6. Communication and prolongation/restriction tasks (`SendU_*`, `RecvU_*`, `Prolongate`, `RestrictU`).
7. `ConToPrim` – convert the updated conserved variables back to primitives.
8. `NewTimeStep` – compute the hydro CFL limit and propagate diffusion/source limits (`src/hydro/hydro_newdt.cpp:21`).

Complementary task lists manage boundary posting/clearing before and after each stage.

## State Arrays

| Array | Shape | Contents |
|-------|-------|----------|
| `u0` | `[nmb, nhydro+nscalars, nk, nj, ni]` | Cell-centered conserved variables. |
| `w0` | `[nmb, nhydro+nscalars, nk, nj, ni]` | Cell-centered primitive variables. |
| `u1` | Same as `u0` | Stage accumulation buffer (RK ≥ 2). |
| `uflx.{x1f,x2f,x3f}` | `[nmb, nhydro+nscalars, nk, nj, ni]` | Face-centered fluxes. |
| `coarse_{u0,w0}` | Coarse-grid analogues (allocated when `mesh/multilevel=true`). |
| `fofc`, `utest` | Allocated only when FOFC is active to mark cells requiring first-order fluxes. |

Passive scalars (`nscalars`) are stored immediately after the hydro primitives/conserved fields, guaranteeing consistent reconstruction and updates across all arrays.

## Input Parameters (`<hydro>` block)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `eos` | – (required) | `ideal` or `isothermal`; SR/GR runs require `ideal`. |
| `gamma` | – | Adiabatic index for the ideal EOS (read by the EOS helper). |
| `iso_sound_speed` | – | Constant sound speed for isothermal runs. |
| `nscalars` | `0` | Number of passively advected scalars. |
| `reconstruct` | `plm` | Reconstruction scheme: `dc`, `plm`, `ppm4`, `ppmx`, or `wenoz`. |
| `rsolver` | – (required) | Riemann solver (see compatibility table below). |
| `fofc` | `false` | Enables first-order flux correction; requires additional ghost zones. |
| `viscosity` | absent | Presence activates isotropic viscosity via `diffusion/viscosity.hpp`. |
| `hyperviscosity` | absent | Positive value activates constant fourth-derivative velocity damping via `diffusion/hyperviscosity.hpp`. |
| `conductivity` / `tdep_conductivity` | absent | Presence activates thermal conduction (`diffusion/conduction.hpp`). |

> The global `<time>` block controls `time/evolution`. Non-relativistic kinematic runs support only the linear advection solver.

### Ghost-Zone Requirements

| Reconstruction | Nominal Ghosts | With FOFC |
|----------------|----------------|-----------|
| `dc` | 2 | 2 |
| `plm` | 2 | 3 |
| `ppm4`, `ppmx`, `wenoz` | 3 | 4 |

Violating these limits triggers explicit fatal errors at construction (`src/hydro/hydro.cpp:147`).

### Riemann Solver Compatibility

| Regime | Supported Keywords | Notes |
|--------|-------------------|-------|
| Non-relativistic dynamic | `llf`, `hlle`, `hllc` (ideal EOS only), `roe` | |
| Non-relativistic kinematic | `advect` | Requires `<time>/evolution = kinematic`. |
| Special relativistic dynamic | `llf`, `hlle`, `hllc` | Mapped internally to SR variants (`Hydro_RSolver::*_sr`). |
| General relativistic dynamic | `llf`, `hlle` | Mapped to GR variants (`Hydro_RSolver::*_gr`). |

Invalid combinations (e.g., `hllc` with isothermal EOS or relativistic kinematic runs) raise fatal errors during construction (`src/hydro/hydro.cpp:189`).

## Flux Calculation Pipeline
`Hydro::Fluxes` dispatches to `CalculateFluxes<Hydro_RSolver::*>`, which reconstructs face states, solves the chosen Riemann problem, and fills `uflx` for all faces on every MeshBlock in the pack (`src/hydro/hydro_fluxes.cpp`). Diffusive contributions are added immediately afterwards:

- `Viscosity::IsotropicViscousFlux` when `viscosity` is present (`src/hydro/hydro_tasks.cpp:115`).
- `HyperViscosity::AddHyperViscousFlux` when positive `hyperviscosity` is present.
- `Conduction::AddHeatFlux` for conductive runs.

FOFC is applied either when explicitly enabled or implicitly when a GR run excises cells (black-hole masks), clamping troubled cells to first-order updates (`src/hydro/hydro_tasks.cpp:121`).

## Source-Term Integration
`Hydro::HydroSrcTerms` invokes the shared `SourceTerms` object to apply body forces, turbulence driving, cooling, etc. The same object also supplies the cooling timestep limit inside `Hydro::NewTimeStep`. Couplings to orbital advection and shearing-box boundaries are coordinated through additional tasks (`SendU_OA`, `RecvU_Shr`, …) and rely on the flags set by the `<shearing_box>` block (`src/hydro/hydro_tasks.cpp:52`).

## Timestep Control
`Hydro::NewTimeStep` executes only on the final stage. It evaluates:

1. The hydrodynamic CFL (or pure advection limit for kinematic runs) using local wave speeds, including relativistic characteristic speeds where appropriate (`src/hydro/hydro_newdt.cpp:39`).
2. Diffusion limits for active conduction, viscosity, scalar diffusion, or hyperviscosity.
3. Source-term limits supplied by `SourceTerms::NewTimeStep`.

The minimum value is stored in `dtnew` and reduced into the driver’s global timestep.

## Performance Considerations
- All heavy kernels (`CopyCons`, `CalculateFluxes`, `RKUpdate`, `ConToPrim`) are implemented as `Kokkos::parallel_for`, operating over the entire mesh-pack to maximise GPU occupancy.
- Arrays are dimensioned to the maximum number of mesh blocks that can reside on the rank, avoiding reallocations during AMR operations.
- When FOFC is enabled, additional storage (`fofc`, `utest`) is allocated only once and reused every step.

## Common Pitfalls
- `hllc` cannot be combined with the isothermal EOS; use `hlle` instead.
- Ensure the mesh ghost-zone count satisfies the reconstruction plus FOFC requirements before enabling the option.
- For kinematic runs, only `rsolver=advect` is valid—dynamic solvers expect conservative energy updates.
- Diffusion coefficients must be declared in `<hydro>`; diffusion blocks elsewhere are not consulted by the hydro module.
- Hyperviscosity requires uniform Newtonian meshes with at least two ghost cells; active AMR/SMR and relativistic modes are rejected.
