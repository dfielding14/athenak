# AGENTS.md

## Purpose
This directory implements the MHD module: conserved/primitive state, magnetic fields,
Riemann solvers, constrained transport, source terms, and boundary communication for
magnetohydrodynamics (including SR/GR variants when enabled).
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Core class
- `mhd.hpp`: `mhd::MHD` class, task IDs, data members, and public task functions.
- `mhd.cpp`: constructor/destructor; parses `<mhd>` options, allocates arrays, selects
  EOS, reconstruction, and Riemann solver.

### Time integration pieces
- `mhd_tasks.cpp`: task list assembly and wrappers (ghost zones, fluxes, CT, sources,
  boundary exchange, prolongation, ConToPrim, timestep).
- `mhd_update.cpp`: explicit RK update of conserved variables (`u0`) from fluxes.
- `mhd_newdt.cpp`: per-pack timestep calculation (CFL + diffusion/source constraints).

### Flux/field helpers
- `mhd_fluxes.cpp`: reconstruction + Riemann solver pipelines and face EMFs.
- `mhd_corner_e.cpp`: corner electric field construction (SG07) and GR/SR variants.
- `mhd_ct.cpp`: constrained transport update for face-centered magnetic fields.
- `mhd_fofc.cpp`: first-order flux correction (FOFC) and BH excision fallback fluxes.

### Riemann solvers
`rsolvers/` contains solver implementations used by `MHD::CalculateFluxes`:
- Non-relativistic: `advect`, `llf`, `hlle`, `hlld`
- Special relativistic: `llf_sr`, `hlle_sr`
- General relativistic: `llf_gr`, `hlle_gr`
- Single-state LLF for FOFC: `llf_mhd_singlestate.hpp`

---

## Data Layout (major fields)
- `u0`, `w0`: conserved / primitive variables on active meshblocks.
- `b0`: face-centered magnetic field (x1f/x2f/x3f).
- `bcc0`: cell-centered magnetic field.
- `u1`, `b1`: intermediate RK registers for U and B.
- `uflx`: face-centered fluxes for conserved variables.
- `efld`: edge-centered electric fields for CT.
- `e?x?` scratch arrays: face-centered EMFs returned by RS (for `CornerE`).
- `coarse_u0`, `coarse_w0`, `coarse_b0`: coarse-grid buffers for SMR/AMR.
- `fofc`, `utest`, `bcctest`: FOFC flags and scratch states.
- `pbval_u`, `pbval_b`: boundary communication for cell- and face-centered fields.
- `porb_u/b`, `psbox_u/b`: orbital advection + shearing box boundary helpers.
- `pvisc`, `presist`, `pcond`, `psrc`: viscosity, resistivity, conduction, sources.

Arrays are dimensioned `[nmb][nvar][k][j][i]` with `nvar = nmhd + nscalars`.

---

## Configuration Inputs (from `<mhd>` block)
Parsed in `mhd.cpp`:
- `eos` (required): `ideal` or `isothermal`
  - `isothermal` is **not** allowed with SR or GR; the constructor exits if
    relativistic coordinates are enabled.
- `rsolver` (required): depends on relativity (see below)
- `reconstruct` (optional, default `plm`): `dc`, `plm`, `ppm4`, `ppmx`, `wenoz`
- `nscalars` (optional, default `0`)
- `fofc` (optional, default `false`)
- `viscosity` (optional) -> enables `Viscosity`
- `ohmic_resistivity` (optional) -> enables `Resistivity`
- `conductivity` / `tdep_conductivity` (optional) -> enables `Conduction`

Note: `rsolver`, `reconstruct`, and `fofc` are only parsed/validated when
`time/evolution != stationary`.

### Solver compatibility (enforced in `mhd.cpp`)
- **Non-relativistic dynamic**: `llf`, `hlle`, `hlld`
- **Non-relativistic kinematic**: `advect`
- **Special relativistic**: `llf`, `hlle` (mapped to SR variants)
- **General relativistic**: `llf`, `hlle` (mapped to GR variants)
  - SR/GR kinematic evolution is not implemented; SR/GR require
    `time/evolution = dynamic`.

### Reconstruction + FOFC ghost-zone requirements
- `plm` + FOFC requires `nghost >= 3`
- `ppm4`, `ppmx`, `wenoz` require `nghost >= 3`
- `ppm4`/`ppmx`/`wenoz` + FOFC requires `nghost >= 4`

---

## Task List Flow (MHD)
`MHD::AssembleMHDTasks` is only called for **single-fluid MHD** (no `hydro`,
`radiation`, `adm`, or `z4c` blocks). Otherwise, task lists are assembled by the
active coupled module (e.g., Radiation or Ion-Neutral) or by NumericalRelativity
when GR spacetime evolution is enabled.

When assembled, it wires tasks into `MeshBlockPack` task lists:

**before_timeintegrator**
- `SaveMHDState` (optional snapshot of primitives/B for time derivatives)

**before_stagen**
- `InitRecv` (post MPI receives for U and B, plus fluxes when needed)

**stagen**
- `CopyCons` -> `Fluxes` -> `SendFlux` -> `RecvFlux` -> `RKUpdate`
- `MHDSrcTerms`
- `SendU_OA` / `RecvU_OA` (orbital advection)
- `RestrictU` -> `SendU` -> `RecvU`
- `SendU_Shr` / `RecvU_Shr` (shearing box BCs)
- `CornerE` -> `EFieldSrc` -> `SendE` -> `RecvE` -> `CT`
- `SendB_OA` / `RecvB_OA` (orbital advection for B)
- `RestrictB` -> `SendB` -> `RecvB`
- `SendB_Shr` / `RecvB_Shr` (shearing box BCs for B)
- `ApplyPhysicalBCs` -> `Prolongate` -> `ConToPrim` -> `NewTimeStep`

**after_stagen**
- `ClearSend` -> `ClearRecv`

---

## Notes on Key Behaviors
- **Flux computation**: Reconstruction (dc/plm/ppm/wenoz) feeds MHD RS; EMFs are
  produced on faces and assembled to corner electric fields (`CornerE`) for CT.
- **Constrained transport**: `CT` updates face-centered B using edge EMFs to preserve
  divergence-free fields.
- **FOFC**: estimates updated U and Bcc, flags floor/excision cells, and replaces
  fluxes with first-order LLF fluxes in those regions (also used for GR excision).
- **Source terms**: `SourceTerms` adds gravity/cooling/shearing-box terms; GR/ADM
  coordinate sources are injected when appropriate.
- **Timestep**: `NewTimeStep` uses fast magnetosonic speeds in Newtonian/SR, but
  clamps characteristic speeds to 1.0 in GR/dynamical GR.

### PR2 Particle Coupling Hooks
- `MHD::EFieldSrc` keeps shearing-box behavior and adds an optional particle term:
  - enabled only when `<particles>/couple_moments_to_mhd=true`
  - additive update to edge-centered `efld` using
    `<particles>/couple_j_to_efield_coeff`
  - representation branch:
    - `cell_centered`: reads deposited moments directly
    - `edge_staggered`: reads converted edge-current arrays from particles.
- Optional fluid momentum/energy feedback is split from E-coupling:
  - `couple_fluid_feedback_order=mhd_src_terms`: feedback applied in
    `MHD::MHDSrcTerms`
  - `couple_fluid_feedback_order=efield_src`: feedback applied in `MHD::EFieldSrc`
  - both paths are stage-1-only and gated by per-target toggles
    (`couple_moments_momentum_to_mhd`,
    `couple_moments_energy_to_mhd`).

---

## Extension Points
- **New Riemann solver**: add implementation in `rsolvers/`, extend `MHD_RSolver`,
  map in `mhd.cpp`, and wire into `CalculateFluxes`.
- **New reconstruction**: implement in `reconstruct/` and add selection in `mhd.cpp`.
- **New EMF / CT scheme**: modify `mhd_corner_e.cpp` and/or `mhd_ct.cpp`.
- **New source term**: extend `SourceTerms` and call in `MHD::MHDSrcTerms`.
