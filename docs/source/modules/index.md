# AthenaK Modules Documentation

Complete documentation for AthenaK modules, organized by category.

## Core Infrastructure

| Module | Directory | Description | Documentation |
|--------|-----------|-------------|---------------|
| **Mesh** | `src/mesh/` | AMR infrastructure, domain decomposition | `mesh.md` |
| **Driver** | `src/driver/` | Time integration and evolution control | `driver.md` |
| **TaskList** | `src/tasklist/` | Task-based execution system | `tasklist.md` |
| **Coordinates** | `src/coordinates/` | Coordinate systems and metrics | `coordinates.md` |

## Physics Modules

| Module | Directory | Description | Documentation |
|--------|-----------|-------------|---------------|
| **Hydrodynamics** | `src/hydro/` | Euler equations solver | `hydro.md` |
| **MHD** | `src/mhd/` | Magnetohydrodynamics with constrained transport | `mhd.md` |
| **Self-Gravity** | `src/gravity/`, `src/multigrid/` | Poisson self-gravity with a multigrid solver | [self_gravity.md](self_gravity.md) |
| **Radiation** | `src/radiation/` | Radiation transport | `radiation.md` |
| **Z4c** | `src/z4c/` | Numerical relativity | `z4c.md` |
| **DynGRMHD** | `src/dyn_grmhd/` | GR MHD in dynamical spacetimes | `dyn_grmhd.md` |
| **Ion-Neutral** | `src/ion-neutral/` | Two-fluid ion-neutral MHD | `ion_neutral.md` |
| **Particles** | `src/particles/` | Lagrangian particle tracking | `particles.md` |

## Numerical Methods

| Module | Directory | Description | Documentation |
|--------|-----------|-------------|---------------|
| **Reconstruction** | `src/reconstruct/` | Spatial reconstruction | `reconstruction.md` |
| **Riemann Solvers** | `src/*/rsolvers/` | Flux computation algorithms | `riemann_solvers.md` |
| **EOS** | `src/eos/` | Equations of state | `eos.md` |
| **Diffusion** | `src/diffusion/` | Physical diffusion processes | `diffusion.md` |

## Support Systems

| Module | Directory | Description | Documentation |
|--------|-----------|-------------|---------------|
| **Outputs** | `src/outputs/` | Data I/O | `outputs.md` |
| **Boundary Values** | `src/bvals/` | Boundary conditions and MPI | `boundaries.md` |
| **Source Terms** | `src/srcterms/` | External sources and turbulence | `srcterms.md` |
| **Shearing Box** | `src/shearing_box/` | Orbital advection | `shearing_box.md` |
| **Problem Generators** | `src/pgen/` | Initial conditions | `pgen.md` |

## Input Parameters

Module parameters are catalogued in
[Input Parameters Reference](../reference/input_parameters.md).

```{toctree}
:hidden:
:maxdepth: 1

self_gravity
```
