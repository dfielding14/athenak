# Module: Magnetohydrodynamics

The MHD module evolves fluid variables together with face-centered magnetic
fields and uses constrained transport (CT) to update those fields. Ordinary
MHD handles non-relativistic, SRMHD, and fixed-background GRMHD paths.
Supplying `<adm>` or `<z4c>` with `<mhd>` selects the separate
[DynGRMHD](dyn_grmhd.md) construction path and its required `dyn_eos` and
`dyn_error` parameters.

## Implementation Map

| Source | Responsibility |
| --- | --- |
| `src/mhd/mhd.hpp`, `mhd.cpp` | State, parameter parsing, solver selection |
| `src/mhd/mhd_tasks.cpp` | Stage tasks and optional couplings |
| `src/mhd/mhd_fluxes.cpp` | Flux kernels |
| `src/mhd/mhd_corner_e.cpp`, `mhd_ct.cpp` | Electric fields and CT update |
| `src/mhd/mhd_fofc.cpp` | First-order flux correction |
| `src/mhd/mhd_newdt.cpp` | Timestep limits |

The public state contains conserved/primitive arrays `u0` and `w0`,
face-centered `b0`, cell-centered `bcc0`, stage registers, face fluxes, and
edge-centered electric fields. Coarse-grid arrays are allocated on
multilevel meshes.

## `<mhd>` Configuration

| Parameter | Requirement/default | Verified behavior |
| --- | --- | --- |
| `eos` | Required | `ideal` or `isothermal`; isothermal is rejected for SR/GR |
| `gamma` | Required by ideal EOS where used | Adiabatic index |
| `iso_sound_speed` | Required by isothermal EOS | Isothermal sound speed |
| `nscalars` | Default `0` | Passive scalar count |
| `reconstruct` | Default `plm` | `dc`, `plm`, `ppm4`, `ppmx`, or `wenoz` |
| `rsolver` | Required for current dynamic, kinematic, and static construction | See table below |
| `fofc` | Default `false` | First-order flux correction |
| `viscosity` | Optional | Constructs MHD viscosity support when present |
| `ohmic_resistivity` | Optional | Constructs resistivity support when present |
| `conductivity`, `tdep_conductivity` | Optional | Either parameter constructs conduction support |

Source terms are conditional: `<mhd_srcterms>` constructs `psrc`. A
`<shearing_box>` block separately constructs orbital-advection and
shearing-box objects for both cell- and face-centered state.

## Solver And Ghost-Zone Constraints

| Regime | `rsolver` values |
| --- | --- |
| Non-relativistic dynamic | `llf`, `hlle`, `hlld` |
| Non-relativistic kinematic | `advect` |
| Special relativistic dynamic | `llf`, `hlle` |
| General relativistic dynamic | `llf`, `hlle` |

The public non-relativistic HLLD implementation contains ideal-gas and
isothermal paths; it is not limited to ideal EOS in the constructor.
SR/GR non-dynamic evolution, including `kinematic` and `static`, is rejected.
Although the driver accepts `time/evolution = static`, the current Newtonian
MHD constructor treats it through the non-dynamic selection branch and
requires `rsolver = advect`.

| Reconstruction | Without FOFC | With FOFC |
| --- | --- | --- |
| `dc` | `nghost >= 2` | `nghost >= 2` |
| `plm` | `nghost >= 2` | `nghost >= 3` |
| `ppm4`, `ppmx`, `wenoz` | `nghost >= 3` | `nghost >= 4` |

SMR/AMR additionally requires an even `nghost`.

## Constrained Transport And Couplings

`MHD::AssembleMHDTasks()` wires flux calculations, edge electric-field
construction, CT, boundary communication, optional refinement transfer, and
timestep evaluation.

- CT updates `b0` through edge-centered electric fields and refreshes
  cell-centered magnetic state for reconstruction.
- Shared source application occurs only if `<mhd_srcterms>` constructed
  `psrc`; shearing-box source terms are a separate coupling.
- Viscosity, resistivity, and conduction fluxes are included only when their
  corresponding parameters constructed those helpers.
- Supplying both `<hydro>` and `<mhd>` uses the ion-neutral coupled route only
  when `<ion-neutral>` is also supplied; unsupported combinations are rejected
  during physics assembly.

## First Verified Use

The built-in Orszag-Tang input is a simple public MHD route:

```bash
./build/src/athena -i inputs/mhd/orszag_tang.athinput \
  -d run-ot time/nlim=1
```

A verified one-cycle run writes `vtk/OrszagTang.mhd_w.00000.vtk` and
`vtk/OrszagTang.mhd_bcc.00000.vtk` in the run directory, followed by the
cycle-1 files.

## See Also

- [Reconstruction](reconstruction.md)
- [Riemann Solvers](riemann_solvers.md)
- [Source Terms](srcterms.md)
- [Shearing Box](shearing_box.md)
