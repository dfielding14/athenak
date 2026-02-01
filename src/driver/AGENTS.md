# AGENTS.md

## Purpose
This directory defines the simulation driver. It wires time integration, task list
execution, output scheduling, AMR coordination, and run initialization/finalization.
See `../../AGENTS.md` for repository-wide conventions and workflow.

---

## Key Files

- `driver.hpp`: `Driver` class interface and runtime parameters (integrator settings,
  time/step limits, wall-clock tracking).
- `driver.cpp`: integrator setup, main execution loop, output scheduling, AMR calls, and
  initialization/finalization logic.

---

## Driver Flow (high level)

1. **Construction**: Reads `<time>` parameters (`evolution`, `integrator`, `tlim`, `nlim`,
   `ndiag`) and sets integrator coefficients.
2. **Initialize**:
   - Sets ghost zones / primitives via `InitBoundaryValuesAndPrimitives`.
   - Computes the initial timestep (if time-evolving).
   - Writes initial outputs for new runs.
   - Allocates ImEx buffers if needed (`impl_src`).
3. **Execute** (time-evolving runs):
   - Executes task lists in order:
     `before_timeintegrator` → per-stage `before_stagen`/`stagen`/`after_stagen`
     → `after_timeintegrator`.
   - Advances time/cycle counters, schedules outputs, and triggers AMR.
   - Recomputes timestep after AMR.
4. **Finalize**: Writes final outputs, calls `pgen_final_func`, and prints diagnostics.

---

## Task Lists and Stage Semantics
- Task lists live in `MeshBlockPack::tl_map` and are executed by name.
- Stages are integer IDs passed to task lists:
  - stage `0` is used for pre-integrator work.
  - stages `1..nexp_stages` map to explicit RK substages.
  - stage `1` is used for `after_timeintegrator` work.

---

## Integrators (as implemented)
The constructor handles these integrator strings:
- `rk1`, `rk2`, `rk3`, `rk4`
- `imex2`, `imex3`, `imex+`

ImEx integrators allocate `impl_src` for stiff sources (currently used for ion-neutral
two-fluid runs). The driver enforces ImEx when ion-neutral physics is enabled.

---

## Output Scheduling
Outputs are triggered when either:
- `time >= last_time + dt` (checked in float precision), or
- `ncycle % dcycle == 0` (if `dcycle` is set).

All outputs are written once in `Finalize` regardless of scheduling.

---

## Initialization Details
`InitBoundaryValuesAndPrimitives` handles:
- Z4c, Hydro, MHD, and Radiation ghost-zone fills and physical BCs.
- Shearing-box boundary communications for Hydro/MHD.
- Prolongation and primitive conversion where appropriate.

---

## Extension Points
- **New integrator**: add coefficients and update the selection logic in `Driver::Driver`.
- **New task list**: register in `MeshBlockPack` and call `ExecuteTaskList` here.
- **New output scheduling policy**: extend the output trigger logic in `Execute`.
