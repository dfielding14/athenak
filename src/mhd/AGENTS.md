# AGENTS.md

## Purpose
This directory implements the MHD module: conserved/primitive state, magnetic fields,
Riemann solvers, constrained transport, source terms, and boundary communication for
magnetohydrodynamics (including SR/GR variants when enabled).
See `../../AGENTS.md` for repository-wide conventions and workflow.
For MHD-PIC equations and signs, use
`docs/source/engineering/pic_mhd_model_contract.md` and
`docs/source/engineering/pic_cr_hall_code_map.md`; active priorities are in
`MHD_PIC_NEXT_STEPS_GUIDE.md`.

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
- `MHDSrcTerms` -> `ApplyPICWaveDamping`
- `CornerE` -> `EFieldSrc` -> `SendE` -> `RecvE` -> `CT`
- PIC expanding-box field/state/feedback hooks
- `SendU_OA` / `RecvU_OA` (orbital advection) -> `RestrictU` -> `SendU` ->
  `RecvU` -> shearing-box U exchange
- `SendB_OA` / `RecvB_OA` (orbital advection for B) -> `RestrictB` -> `SendB` ->
  `RecvB` -> shearing-box B exchange
- `ApplyPhysicalBCs` -> `Prolongate` -> `ConToPrim` -> `NewTimeStep`

`Fluxes` adds the full CR-Hall face induction and energy flux before its internal
FOFC call. With `time/integrator=vl2`, full-f TSC particle tasks are inserted
after `CopyCons` and before `Fluxes`; see `src/particles/AGENTS.md` for their
stage ordering.

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
  For full-f VL2/TSC, the trial includes both the Hall-enriched face fluxes and
  the matching stage particle source. On replaced faces, the full-Hall path then
  restores donor-cell Hall induction and energy terms before clearing the flags.
- **Source terms**: `SourceTerms` adds gravity/cooling/shearing-box terms; GR/ADM
  coordinate sources are injected when appropriate.
- **Timestep**: `NewTimeStep` uses fast magnetosonic speeds in Newtonian/SR, but
  clamps characteristic speeds to 1.0 in GR/dynamical GR.

### Particle Coupling Hooks
- `MHD::EFieldSrc` keeps shearing-box behavior and adds an optional particle term:
  - enabled only when `<particles>/couple_moments_to_mhd=true`
  - additive update to edge-centered `efld` using
    `<particles>/couple_j_to_efield_coeff`
  - representation branch:
    - `cell_centered`: reads deposited moments directly
    - `edge_staggered`: reads particle-owned edge-current arrays; those arrays
      come from either CC-to-edge conversion or direct staggered deposition.
- Optional fluid momentum/energy feedback is split from E-coupling:
  - `couple_fluid_feedback_order=mhd_src_terms`: feedback applied in
    `MHD::MHDSrcTerms`
  - `couple_fluid_feedback_order=efield_src`: feedback applied in `MHD::EFieldSrc`
  - non-VL2 runs build deposited rates once per cycle and reuse them across RK
    stages with the stage beta coefficient; targets are gated by
    `couple_moments_momentum_to_mhd` and `couple_moments_energy_to_mhd`.
  - full-f VL2/TSC instead applies the analytic deposited charge/current force
    in stage 1 and the opposite of the deposited, realized particle
    `dp/dt`/`dE/dt` impulse in stage 2. This same helper updates the live state
    and the FOFC trial state.
  - other coupled Boris paths consume deposited per-step particle deltas with
    opposite sign for conservative exchange; non-Boris paths retain legacy
    `J x B` / `J dot B` source handling.
  - every ideal-MHD feedback path validates the post-source trial state against
    density, pressure/internal-energy, temperature, and entropy floors before
    C2P; a would-be floor repair aborts with global event counts so a coupled
    conservation run cannot silently acquire floor energy.

For the uniform-grid full-f VL2/TSC path, `pic_cr_hall_mode=off|full` selects
ideal-MHD or full CR-Hall induction. The full closure's predictor uses deposited
charge/current to build `v_H`; the stage-2 particle push uses the
predicted midpoint field, then deposits the realized particle impulse. The grid
corrector uses that impulse directly,
`cE_H = -(dp_CR/dt)/(alpha_i rho_g)`, rather than reusing the predictor `v_H`.
`Fluxes` adds the limited Hall face EMFs and matching `(cE_H x B)` energy flux
before FOFC, flux communication, and RK update. `CornerE` assembles the total
face/cell EMF for GS07, CT updates the staggered field, and admissibility is
checked against the updated face field. Start with
`docs/source/engineering/pic_cr_hall_code_map.md` before changing this route.

### Passive-MHD Isolation Hook
- When `<particles>/pic_background_mode=passive_mhd` is active, MHD fluid
  evolution is frozen by short-circuiting:
  - `MHD::RKUpdate`
  - `MHD::MHDSrcTerms`
  - `MHD::CT`
- Task-list structure and MHD state containers remain present so PIC push and
  deposited-moment diagnostics can run against stable background fields.

---

## Extension Points
- **New Riemann solver**: add implementation in `rsolvers/`, extend `MHD_RSolver`,
  map in `mhd.cpp`, and wire into `CalculateFluxes`.
- **New reconstruction**: implement in `reconstruct/` and add selection in `mhd.cpp`.
- **New EMF / CT scheme**: modify `mhd_corner_e.cpp` and/or `mhd_ct.cpp`.
- **New source term**: extend `SourceTerms` and call in `MHD::MHDSrcTerms`.
