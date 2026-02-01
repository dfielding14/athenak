# AGENTS.md

## Purpose
This directory implements the hydrodynamics module: state storage, reconstruction,
Riemann solves, flux updates, source terms, and boundary communication for hydro-only
runs (and for hydro components of multi-physics when enabled).
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Core class
- `hydro.hpp`: `hydro::Hydro` class, task IDs, data members, and public task functions.
- `hydro.cpp`: constructor/destructor; parses `<hydro>` options, allocates arrays,
  selects EOS, reconstruction, and Riemann solver.

### Time integration pieces
- `hydro_tasks.cpp`: task list assembly and task wrappers (ghost zones, fluxes, sources,
  boundary exchange, prolongation, ConToPrim, timestep).
- `hydro_update.cpp`: explicit RK update of conserved variables (`u0`) from fluxes.
- `hydro_newdt.cpp`: per-pack timestep calculation (CFL + diffusion/source constraints).

### Flux and stability helpers
- `hydro_fluxes.cpp`: reconstruction + Riemann solver pipelines for x1/x2/x3 fluxes.
- `hydro_fofc.cpp`: first-order flux correction (FOFC) and BH excision fallback fluxes.

### Riemann solvers
`rsolvers/` contains solver implementations used by `Hydro::CalculateFluxes`:
- Non-relativistic: `advect`, `llf`, `hlle`, `hllc`, `roe`
- Special relativistic: `llf_sr`, `hlle_sr`, `hllc_sr`
- General relativistic: `llf_gr`, `hlle_gr`
- Single-state LLF for FOFC: `llf_hyd_singlestate.hpp`

---

## Data Layout (major fields)
- `u0`, `w0`: conserved / primitive variables on active meshblocks.
- `u1`: intermediate conserved state for RK stages.
- `uflx`: face-centered fluxes for each conserved variable.
- `coarse_u0`, `coarse_w0`: coarse-grid buffers for SMR/AMR.
- `fofc`, `utest`: FOFC flag + scratch arrays.
- `pbval_u`: cell-centered boundary communication (MPI + buffers).
- `porb_u`, `psbox_u`: orbital advection + shearing box boundary helpers when enabled.
- `pvisc`, `pcond`, `psrc`: viscosity, conduction, and source term objects.

Arrays are dimensioned `[nmb][nvar][k][j][i]` with `nvar = nhydro + nscalars`.

---

## Configuration Inputs (from `<hydro>` block)
Parsed in `hydro.cpp`:
- `eos` (required): `ideal` or `isothermal`
  - `isothermal` is **not** allowed with SR or GR; the constructor exits if
    `<coord>` enables relativistic coordinates.
- `rsolver` (required): depends on relativity & EOS (see below)
- `reconstruct` (optional, default `plm`): `dc`, `plm`, `ppm4`, `ppmx`, `wenoz`
- `nscalars` (optional, default `0`)
- `fofc` (optional, default `false`)
- `viscosity` (optional) -> enables `Viscosity`
- `conductivity` / `tdep_conductivity` (optional) -> enables `Conduction`

Note: `rsolver`, `reconstruct`, and `fofc` are only parsed/validated when
`time/evolution != stationary`.

### Solver compatibility (enforced in `hydro.cpp`)
- **Non-relativistic dynamic**: `llf`, `hlle`, `hllc`, `roe`
  - `hllc` requires ideal EOS (not isothermal).
- **Non-relativistic kinematic**: `advect`
- **Special relativistic**: `llf`, `hlle`, `hllc` (mapped to SR variants)
- **General relativistic**: `llf`, `hlle` (mapped to GR variants)
  - SR/GR kinematic evolution is not implemented; SR/GR require
    `time/evolution = dynamic`.

### Reconstruction + FOFC ghost-zone requirements
- `plm` + FOFC requires `nghost >= 3`
- `ppm4`, `ppmx`, `wenoz` require `nghost >= 3`
- `ppm4`/`ppmx`/`wenoz` + FOFC requires `nghost >= 4`

---

## Task List Flow (Hydro)
`Hydro::AssembleHydroTasks` is only called for **single-fluid hydro** (no `mhd`,
`radiation`, `adm`, or `z4c` blocks). Otherwise, task lists are assembled by the
active coupled module (e.g., Radiation or Ion-Neutral) or by NumericalRelativity
when GR spacetime evolution is enabled.

When assembled, it wires tasks into `MeshBlockPack` task lists:

**before_stagen**
- `InitRecv` (post MPI receives and initialize boundary status)

**stagen**
- `CopyCons` -> `Fluxes` -> `SendFlux` -> `RecvFlux` -> `RKUpdate`
- `HydroSrcTerms`
- `SendU_OA` / `RecvU_OA` (orbital advection, shearing box)
- `RestrictU` -> `SendU` -> `RecvU`
- `SendU_Shr` / `RecvU_Shr` (shearing box BCs)
- `ApplyPhysicalBCs` -> `Prolongate` -> `ConToPrim` -> `NewTimeStep`

**after_stagen**
- `ClearSend` -> `ClearRecv`

---

## Notes on Key Behaviors
- **Flux computation**: Reconstruction (dc/plm/ppm/wenoz) feeds Riemann solvers in
  `Hydro::CalculateFluxes`; scalar fluxes are upwinded using density flux sign.
- **FOFC**: estimates an updated state, flags cells requiring floors, and replaces
  fluxes with first-order LLF fluxes in those regions. Also used for GR excision.
- **Source terms**: `SourceTerms` hooks add gravity/cooling/shearing-box terms and
  user sources; GR coordinate sources are added when appropriate.
- **Timestep**: `NewTimeStep` computes dt from wave speeds, with diffusion and source
  term constraints forwarded to `pcond` and `psrc`.

---

## Extension Points
- **New Riemann solver**: add implementation in `rsolvers/`, extend `Hydro_RSolver`,
  map in `hydro.cpp`, and wire into `CalculateFluxes`.
- **New reconstruction**: implement in `reconstruct/` and add selection in `hydro.cpp`.
- **New source term**: extend `SourceTerms` and call in `Hydro::HydroSrcTerms`.
- **New diagnostics**: add outputs in `src/outputs` and map to hydro variables.
