# AGENTS.md

## Purpose
This directory implements dynamical general-relativistic MHD (DynGRMHD) coupling
for runs that evolve spacetime (`<adm>` or `<z4c>`) with MHD fluids. It replaces
standard MHD fluxes with GR-aware fluxes, supplies GR coordinate source terms,
and provides stress-energy (`Tmunu`) for Z4c.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Core classes
- `dyn_grmhd.hpp`: `dyngr::DynGRMHD` base class and `DynGRMHDPS` template with EOS
  and error-policy selection, plus task-queue interface.
- `dyn_grmhd.cpp`: factory (`BuildDynGRMHD`), configuration parsing, task queuing,
  primitive/conserved conversion, GR source terms, and `Tmunu` construction.

### Fluxes and FOFC
- `dyn_grmhd_fluxes.cpp`: reconstruction + Riemann solver pipelines and scalar
  flux upwinding; calls FOFC when enabled.
- `dyn_grmhd_fofc.cpp`: first-order flux correction (FOFC) and GR excision fluxes.
- `dyn_grmhd_util.hpp`: helper routines for primitives, B-fields, and flux packing.

### Riemann solvers
`rsolvers/` contains GRMHD solvers used by DynGRMHD:
- `llf_dyn_grmhd.hpp`
- `hlle_dyn_grmhd.hpp`
- `flux_dyn_grmhd.hpp` (shared helpers)

---

## How DynGRMHD Is Activated
- In `MeshBlockPack::AddPhysics`, DynGRMHD is constructed only when:
  - `<mhd>` is present, and
  - either `<adm>` or `<z4c>` is present.
- If `<adm>`/`<z4c>` appears with `<hydro>` (without `<mhd>`), the code exits.
- When DynGRMHD is active, `NumericalRelativity` builds task lists and calls
  `DynGRMHD::QueueDynGRMHDTasks`.

---

## Configuration Inputs (from `<mhd>` block)
Parsed in `dyn_grmhd.cpp`:
- `rsolver` (required): `llf` or `hlle` for dynamical GR.
- `fofc_method` (optional, default `llf`): stored in `fofc_method`.
- `dyn_eos` (required): `ideal`, `piecewise_poly`, or `compose`.
- `dyn_error` (required): `reset_floor`.
- `dyn_scratch` (optional, default `0`): Kokkos scratch level for flux kernels.
- `enforce_maximum` (optional, default `true`): enable DMP checks in FOFC.
- `dmp_M` (optional, default `1.2`): DMP threshold multiplier.
- `fixed` (optional, default `false`): disables MHD evolution when true.

---

## Data Ownership and Layout
DynGRMHD operates on the existing MHD state in `MeshBlockPack`:
- Conserved: `pmhd->u0`, `pmhd->u1`
- Primitive: `pmhd->w0`
- Magnetic fields: face-centered `pmhd->b0`, cell-centered `pmhd->bcc0`
- Fluxes and EMFs: `pmhd->uflx`, `pmhd->e?x?`

Additional dependencies:
- Metric and gauge data from `padm->adm`.
- Stress-energy tensor storage in `ptmunu->tmunu`.
- DynGRMHD EOS/primitive solver: `PrimitiveSolverHydro<EOSPolicy, ErrorPolicy>`
  stored as `DynGRMHDPS::eos`.

---

## Task List Flow (NumericalRelativity)
DynGRMHD queues tasks into `NumericalRelativity`, which then assembles the
`before_stagen`, `stagen`, and `after_stagen` lists.

**before_stagen**
- `MHD::InitRecv` (MHD_Recv)

**stagen**
- `MHD::CopyCons` (MHD_CopyU)
- `DynGRMHDPS::CalcFluxes` (MHD_Flux) with GR Riemann solvers
- `DynGRMHD::SetTmunu` (MHD_SetTmunu) only when Z4c is active
- `MHD::SendFlux` -> `MHD::RecvFlux`
- `MHD::RKUpdate` (depends on `MHD_SetTmunu` when present)
- `MHD::MHDSrcTerms` (includes `AddCoordTerms` for dynamical GR)
- `MHD::RestrictU` -> `MHD::SendU` -> `MHD::RecvU`
- CT chain: `MHD::CornerE` -> `MHD::SendE` -> `MHD::RecvE` -> `MHD::CT`
  -> `MHD::RestrictB` -> `MHD::SendB` -> `MHD::RecvB`
- `MHD::ApplyPhysicalBCs` (DynGRMHD has its own BC method, but it is not queued)
- `MHD::Prolongate`
- `DynGRMHDPS::ConToPrim` (optional dependency on `Z4c_Excise`)
- `MHD::NewTimeStep`

**after_stagen**
- `MHD::ClearSend`
- `MHD::ClearRecv`

---

## Algorithm Notes
- **Fluxes**: `CalcFluxes` reconstructs primitives and computes GRMHD fluxes using
  `LLF_DYNGR` or `HLLE_DYNGR`. Scalar fluxes are upwinded using the density flux.
- **FOFC and excision**: `CalcFluxes` calls `FOFC` when `pmhd->use_fofc` or when
  `coord_data.bh_excise` is set. FOFC estimates updated states, flags cells via
  DMP checks, and replaces fluxes with single-state LLF/HLLE fluxes at flagged
  faces and excision regions.
- **Coordinate sources**: `MHD::MHDSrcTerms` calls `DynGRMHD::AddCoordTerms` when
  `pcoord->is_dynamical_relativistic` is true; this adds GR source terms to the
  conserved variables using metric derivatives and the dynamical EOS.
- **Stress-energy**: `SetTmunu` builds the perfect-fluid stress-energy tensor
  from the MHD state and the ADM metric; it is queued only when Z4c is active.
- **Primitive/conserved conversions**:
  - `ConToPrim`/`ConToPrimBC` use the dynamical EOS via `PrimitiveSolverHydro`.
  - `PrimToConInit` is used during BC handling and some problem generators.
  - `ConvertInternalEnergyToPressure` is used by `pgen/lorene_bns.cpp` to update
    pressure from internal energy for dynamical EOS runs.
- **Fixed evolution**: when `<mhd>/fixed = true`, fluxes, `ConToPrim`,
  `AddCoordTerms`, and `SetTmunu` return early without modifying state.

---

## Cautions
- `fofc_method` is parsed and stored but is not used to select the FOFC solver;
  FOFC currently uses the same solver template as `CalcFluxes`.
- DynGRMHD implements its own `ApplyPhysicalBCs`, but the task list uses
  `MHD::ApplyPhysicalBCs` (the DynGRMHD version is commented out in the queue).
- Orbital-advection/shearing-box send/recv tasks are not included in the
  NumericalRelativity queue for DynGRMHD.

---

## Extension Points
- **New EOS policy**: add an enum entry in `DynGRMHD_EOS`, extend `BuildDynGRMHD`,
  instantiate `DynGRMHDPS` templates, and ensure EOS bindings are available in
  `PrimitiveSolverHydro`.
- **New error policy**: add an enum entry in `DynGRMHD_Error`, extend the factory,
  and provide a matching `PrimitiveSolverHydro` error policy.
- **New Riemann solver**: add implementation in `rsolvers/`, extend
  `DynGRMHD_RSolver`, update the solver selection in `DynGRMHD::DynGRMHD`, and
  instantiate `CalcFluxes`/`FOFC` for the new solver.
- **Task sequencing**: modify `QueueDynGRMHDTasks` to insert additional GR-specific
  steps or to switch to `DynGRMHD::ApplyPhysicalBCs`.
