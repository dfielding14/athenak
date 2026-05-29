# Module: Hydrodynamics

The Hydro module owns cell-centered fluid conserved and primitive state and
advances it through flux, update, boundary, source, and timestep tasks. The
same `<hydro>` block is used for non-relativistic, special-relativistic, and
fixed-background general-relativistic configurations; coordinates determine
which implemented EOS and Riemann-solver variants are selected. Dynamical
metric evolution (`<adm>` or `<z4c>`) cannot be combined with `<hydro>` in
the public physics-assembly path.

## Implementation Map

| Source | Responsibility |
| --- | --- |
| `src/hydro/hydro.hpp`, `hydro.cpp` | State, constructor parsing, solver selection |
| `src/hydro/hydro_tasks.cpp` | Stage task wiring and conditional couplings |
| `src/hydro/hydro_fluxes.cpp` | Reconstruct and solve interface fluxes |
| `src/hydro/hydro_update.cpp` | Runge-Kutta state update and primitive conversion |
| `src/hydro/hydro_newdt.cpp` | CFL, diffusion, and source-term timestep limits |

The public state includes `u0` and `w0`, optional coarse-grid copies on
multilevel meshes, stage storage `u1`, face fluxes `uflx`, and FOFC scratch
storage only when FOFC is enabled.

## `<hydro>` Configuration

| Parameter | Requirement/default | Verified behavior |
| --- | --- | --- |
| `eos` | Required | `ideal` or `isothermal`; isothermal is rejected for SR/GR |
| `gamma` | Required by ideal EOS where used | Adiabatic index |
| `iso_sound_speed` | Required by isothermal EOS | Isothermal sound speed |
| `nscalars` | Default `0` | Passive scalar count |
| `reconstruct` | Default `plm` | `dc`, `plm`, `ppm4`, `ppmx`, or `wenoz` |
| `rsolver` | Required for current dynamic, kinematic, and static construction | See table below |
| `fofc` | Default `false` | First-order flux correction |
| `viscosity` | Optional | Constructs `Viscosity("hydro", ...)` when present |
| `conductivity`, `tdep_conductivity` | Optional | Either parameter constructs conduction support |

Source terms are not automatically constructed. A `<hydro_srcterms>` block
constructs the shared `SourceTerms` object; a `<shearing_box>` block
separately constructs orbital-advection and shearing-box helpers.

## Solver And Ghost-Zone Constraints

| Regime | `rsolver` values |
| --- | --- |
| Non-relativistic dynamic | `llf`, `hlle`, `hllc`, `roe`; `hllc` requires ideal EOS |
| Non-relativistic kinematic | `advect` |
| Special relativistic dynamic | `llf`, `hlle`, `hllc` |
| General relativistic dynamic | `llf`, `hlle` |

SR/GR non-dynamic evolution, including `kinematic` and `static`, is rejected
by the constructor. Although the driver accepts `time/evolution = static`,
the current Newtonian Hydro constructor treats it through the non-dynamic
selection branch and requires `rsolver = advect`.

| Reconstruction | Without FOFC | With FOFC |
| --- | --- | --- |
| `dc` | `nghost >= 2` | `nghost >= 2` |
| `plm` | `nghost >= 2` | `nghost >= 3` |
| `ppm4`, `ppmx`, `wenoz` | `nghost >= 3` | `nghost >= 4` |

SMR/AMR additionally requires an even number of ghost zones, so a refined run
that otherwise needs three must set `nghost >= 4`.

## Task And Coupling Behavior

`Hydro::AssembleHydroTasks()` wires flux calculation, Runge-Kutta update,
conditional source application, boundary exchange, prolongation/restriction,
primitive recovery, and timestep evaluation. In particular:

- `psrc->ApplySrcTerms(...)` is called only when `<hydro_srcterms>` created
  `psrc`.
- Shearing-box source application is separate and occurs only when a
  `<shearing_box>` block constructed `psbox_u`.
- A problem generator can enroll `user_srcs_func` and `user_bcs_func`.
- `Hydro::NewTimeStep()` adds shared source-term limits only when `psrc`
  exists.

## First Verified Use

The public Sod input is a simple Hydro validation path:

```bash
./build/src/athena -i inputs/hydro/sod.athinput -d run-sod
```

See [Quickstart](../quickstart.md) for the build command and verified output
files.

## See Also

- [Reconstruction](reconstruction.md)
- [Riemann Solvers](riemann_solvers.md)
- [Source Terms](srcterms.md)
- [Outputs](outputs.md)
