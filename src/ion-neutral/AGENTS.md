# AGENTS.md

## Purpose
This directory implements the ion-neutral (two-fluid) coupling used when both Hydro
(neutrals) and MHD (ions) are enabled. It adds implicit drag/ionization/recombination
updates and wires a custom task list that advances both fluids together.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files and Responsibilities

### Core class
- `ion-neutral.hpp`: `ion_neutral::IonNeutral` class, task IDs, coupling coefficients,
  and implicit RK update interfaces.
- `ion-neutral.cpp`: constructor and `<ion-neutral>` parameter parsing.
- `ion-neutral_tasks.cpp`: task-list assembly and the implicit coupling algorithms
  (`FirstTwoImpRK`, `ImpRKUpdate`).

---

## How Ion-Neutral Is Activated
- Constructed in `MeshBlockPack::AddPhysics` when an `<ion-neutral>` block exists.
- Requires both `<hydro>` and `<mhd>` blocks; otherwise the code exits with a fatal
  error. Ion-neutral is also rejected if `<adm>` or `<z4c>` are present.
- The task list is built by `IonNeutral::AssembleIonNeutralTasks` in this directory.

---

## Configuration Inputs (from `<ion-neutral>` block)
Parsed in `ion-neutral.cpp`:
- `drag_coeff` (required): ion-neutral coupling coefficient, gamma.
- `ionization_coeff` (optional, default `0.0`): ionization rate, xi.
- `recombination_coeff` (optional, default `0.0`): recombination rate, alpha.

---

## Data Ownership and Layout
- Ion-neutral stores no separate state arrays; it operates directly on the existing
  MHD and Hydro conserved arrays (`pmhd->u0`, `phydro->u0`) and the MHD magnetic
  fields (`pmhd->b0`).
- It uses `Driver::impl_src` as a scratch array for stiff source terms. This array is
  allocated in `Driver::Initialize` when `pionn != nullptr`.
- The implicit-source index mapping (from `ImpRKUpdate`) is:
  - `ru_(s,m,0..2,...)` -> ion momentum `IM1..IM3`
  - `ru_(s,m,3..5,...)` -> neutral momentum `IM1..IM3`
  - `ru_(s,m,6,...)` -> ion density `IDN`
  - `ru_(s,m,7,...)` -> neutral density `IDN`

---

## Task List Flow (Ion-Neutral)
Assembled in `ion-neutral_tasks.cpp`:

**before_stagen**
- `MHD::InitRecv` (ions)
- `Hydro::InitRecv` (neutrals)

**stagen**
- `IonNeutral::FirstTwoImpRK` (copies states, performs two implicit stages)
- MHD explicit chain: `Fluxes` -> `SendFlux` -> `RecvFlux` -> `RKUpdate`
- Hydro explicit chain: `Fluxes` -> `SendFlux` -> `RecvFlux` -> `RKUpdate`
- `IonNeutral::ImpRKUpdate` (implicit drag/ionization/recombination)
- `MHD::RestrictU` -> `Hydro::RestrictU`
- `MHD::SendU` + `Hydro::SendU` -> `MHD::RecvU` + `Hydro::RecvU`
- MHD CT chain: `CornerE` -> `SendE` -> `RecvE` -> `CT` -> `RestrictB`
  -> `SendB` -> `RecvB`
- `MHD::ApplyPhysicalBCs` and `Hydro::ApplyPhysicalBCs`
- `MHD::Prolongate` and `Hydro::Prolongate`
- `MHD::ConToPrim` and `Hydro::ConToPrim`
- `MHD::NewTimeStep` and `Hydro::NewTimeStep`

**after_stagen**
- `MHD::ClearSend`
- `Hydro::ClearSend`

---

## Algorithm Notes
- `FirstTwoImpRK` runs only on stage 1. It copies `u0 -> u1` (and `b0 -> b1` for MHD),
  applies two implicit updates (estage = -1, 0), and then converts conserved to
  primitive variables for both fluids over the full domain including ghost zones.
- `ImpRKUpdate`:
  - Adds previous implicit contributions using `Driver::a_twid` (for istage > 1).
  - Applies analytic implicit updates for drag/ionization/recombination using
    `Driver::a_impl * dt`.
  - Computes and stores the stiff source term `R(U^n)` in `Driver::impl_src` for
    use in later stages.
- Ions are evolved with the MHD solver; neutrals are evolved with the Hydro solver.

---

## Constraints and Cautions
- Ion-neutral runs require an ImEx integrator (`time/integrator` must be `imex2`,
  `imex3`, or `imex+`). The driver exits if `pionn != nullptr` and `nimp_stages == 0`.
- The ion-neutral task list does **not** include Hydro/MHD source-term tasks,
  orbital-advection tasks, shearing-box exchanges, or `ClearRecv`. If you need those
  in two-fluid runs, extend `IonNeutral::AssembleIonNeutralTasks` explicitly.
- `pionn` is allocated in `MeshBlockPack::AddPhysics` and stored on the pack. If you
  add resources inside `IonNeutral`, ensure cleanup is handled (the pack destructor
  currently does not delete `pionn`).

---

## Related Areas
- `src/mesh/meshblock_pack.cpp`: constructs `IonNeutral` and enforces block
  compatibility.
- `src/driver/driver.cpp`: ImEx integrator setup and `impl_src` allocation.
- `src/srcterms/turb_driver.cpp`: turbulence forcing inserts tasks using ion-neutral
  task IDs when `pionn` is present.

---

## Extension Points
- Add new coupling terms in `ImpRKUpdate` and update the `impl_src` mapping/size.
- Add or reorder tasks in `AssembleIonNeutralTasks` if new physics needs to be
  staged with the explicit/implicit updates.
