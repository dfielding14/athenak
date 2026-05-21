# Module: Magnetohydrodynamics

## Role in AthenaK
The Magnetohydrodynamics (MHD) module extends the hydrodynamics integrator to evolve ideal MHD, special relativistic MHD (SRMHD), or general relativistic MHD (GRMHD). It owns face-centred magnetic fields, applies constrained transport, and coordinates optional diffusion modules. The implementation lives in `src/mhd/`.

## File Layout

| File | Purpose | Key Routines |
|------|---------|--------------|
| `mhd.hpp/cpp` | `MHD` class, construction, data ownership | constructor, `AssembleMHDTasks`, accessors |
| `mhd_tasks.cpp` | Task graph for RK stages | wrapper tasks (`Fluxes`, `CornerE`, `CT`, …) |
| `mhd_fluxes.cpp` | Template flux kernel per MHD Riemann solver | `CalculateFluxes<MHD_RSolver::*>` |
| `mhd_corner_e.cpp` | Corner electric field computation | `CornerE` |
| `mhd_ct.cpp` | Constrained transport update | `CT` |
| `mhd_newdt.cpp` | CFL and diffusion/source timestep limit | `MHD::NewTimeStep` |
| `mhd_fofc.cpp` | First-order flux correction driver | `FOFC` |

## Task Flow
`MeshBlockPack::AddPhysics` instantiates the `MHD` object when `<mhd>` exists. If no other fluid modules (Hydro, Radiation, ADM/Z4c) are present, `AssembleMHDTasks` inserts the module’s tasks into the integrator pipeline. For two-fluid runs (`<ion-neutral>`), the ion-neutral driver assembles its own task list that embeds the MHD operations.

Per-stage task ordering (`tl["stagen"]`):

1. `CopyCons` – stage register setup for higher-order integrators.
2. `Fluxes` – reconstruct interface states, run the chosen MHD Riemann solver, add viscosity, resistivity, and conduction fluxes.
3. `SendFlux` / `RecvFlux` – AMR flux exchange.
4. `RKUpdate` – conserved-variable update.
5. `HydroSrcTerms` equivalent – applies shared `SourceTerms` to the MHD state.
6. Communication tasks for face-centred fields (`SendU`, `RecvU`, `SendB`, `RecvB`, etc.).
7. `CornerE` – compute edge-centred electric fields.
8. `CT` – constrained transport update of face-centred `B`.
9. `ConToPrim` – update primitive variables (`w0`/`bcc0`).
10. `NewTimeStep` – evaluate the CFL limit and diffusion/source contributions.

Boundary posting and clearing live in the `before_stagen` / `after_stagen` lists.

## State Arrays

| Array | Shape | Contents |
|-------|-------|----------|
| `u0`, `u1` | `[nmb, nmhd+nscalars, nk, nj, ni]` | Conserved variables + passive scalars. |
| `w0` | Same as `u0` | Primitive variables (density, velocity, pressure). |
| `b0.{x1f,x2f,x3f}` | `[nmb, nk, nj, ni±1]` | Face-centred magnetic fluxes. |
| `bcc0` | `[nmb, 3, nk, nj, ni]` | Cell-centred magnetic fields for primitive reconstruction. |
| `efld.{x1f,x2f,x3f}` / `e*` | Edge-centred electric fields used by CT. |
| `coarse_*` | Coarse-grid analogues allocated when AMR is enabled. |
| `fofc`, `utest`, `bcctest`, `wsaved`, `bccsaved` | Optional buffers used by FOFC and CT safeguard logic. |

Passive scalars are appended to hydro variables in both `u0` and `w0`, ensuring consistent reconstruction and updates.

## Input Parameters (`<mhd>` block)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `eos` | – (required) | `ideal` or `isothermal` (SR/GR require `ideal`). |
| `gamma` | – | Adiabatic index for ideal EOS. |
| `iso_sound_speed` | – | Constant sound speed for isothermal runs. |
| `nscalars` | `0` | Number of passively advected scalars. |
| `reconstruct` | `plm` | Reconstruction method (`dc`, `plm`, `ppm4`, `ppmx`, `wenoz`). |
| `rsolver` | – (required) | Riemann solver keyword (see below). |
| `fofc` | `false` | Enable first-order flux correction. |
| `viscosity` | absent | Enables isotropic viscosity (diffusion module). |
| `ohmic_resistivity` | absent | Enables resistive diffusion. |
| `conductivity` / `tdep_conductivity` | absent | Enables thermal conduction. |

The constructor enforces EOS/solver compatibility and throws fatal errors for unsupported combinations (`src/mhd/mhd.cpp:36`).

### Ghost-Zone Requirements

| Reconstruction | Nominal | With FOFC |
|----------------|---------|-----------|
| `dc` | 2 | 2 |
| `plm` | 2 | 3 |
| `ppm4`, `ppmx`, `wenoz` | 3 | 4 |

(Identical to the Hydro module; see checks in `src/mhd/mhd.cpp:104`.)

### Riemann Solver Compatibility

| Regime | Supported keywords | Notes |
|--------|-------------------|-------|
| Non-relativistic dynamic | `llf`, `hlle`, `hlld` | HLLD available only with `ideal` EOS. |
| Special relativistic dynamic | `llf`, `hlle` | Mapped to `MHD_RSolver::*_sr`. |
| General relativistic dynamic | `llf`, `hlle` | Mapped to `MHD_RSolver::*_gr`. |

No Roe solver is implemented for MHD; specifying it raises a fatal error.

## Diffusion & Source Coupling
- Viscosity, resistivity, and conduction modules are instantiated only when the corresponding coefficients appear in `<mhd>`. Their fluxes are added immediately after the Riemann solver inside `MHD::Fluxes` (`src/mhd/mhd_tasks.cpp:121`).
- `SourceTerms("mhd", …)` applies body forces, turbulence driving, etc., and contributes to the timestep limiter.
- Orbital advection and shearing-box boundary handlers (`OrbitalAdvection{CC,FC}`, `ShearingBoxBoundary{CC,FC}`) are constructed when `<shearing_box>` exists.

## Constrained Transport (CT)
`CornerE` computes edge-centred electric fields using the current fluxes and velocities. `CT` updates each face-centred `B` component via the curl of these EMFs, preserving $\nabla \cdot \mathbf{B} = 0` to machine precision (`src/mhd/mhd_ct.cpp`). Cell-centred fields (`bcc0`) are then refreshed by averaging the face values.

## Timestep Control
`MHD::NewTimeStep` executes on the final stage and evaluates:

1. Fast-magnetosonic CFL speed (or relativistic characteristic speeds if SR/GR is active).
2. Diffusion limits from resistivity and conduction.
3. Source-term limits from the shared `SourceTerms` instance.

The tightest value is stored in `dtnew` (`src/mhd/mhd_newdt.cpp`).

## Compatibility & Couplings
- If both `<hydro>` and `<mhd>` blocks are present, `<ion-neutral>` must also be specified; otherwise construction aborts (`src/mesh/meshblock_pack.cpp:138`).
- Radiation and relativistic metric modules can co-exist; task assembly is handled conditionally in `MeshBlockPack::AddPhysics`.
- Turbulence driving applies to the MHD pack when `<turb_driving>` is present.

## Operational Notes
- For SR/GR runs, only dynamic evolution is supported (`<time>/evolution = dynamic`).
- HLLD requires an ideal EOS; use HLLE otherwise.
- FOFC is automatically enabled in excision regions when running GRMHD with black hole masks (`src/mhd/mhd_tasks.cpp:127`).
- Allocate sufficient ghost zones before enabling higher-order reconstruction plus FOFC.


### $\nabla \cdot \mathbf{B}$ Growth
- Ensure using CT (always enabled in AthenaK)
- Check boundary conditions
- Verify AMR prolongation/restriction if using refinement

### Carbuncle Instability
- Use HLLD solver (most accurate)
- Add small resistivity
- Increase resolution

## See Also
- [Hydro Module](hydro.md) - Hydrodynamics solver
- [Reconstruction Module](reconstruction.md) - Reconstruction methods
- [Riemann Solvers Module](riemann_solvers.md) - Solver details
- Source: `src/mhd/`
