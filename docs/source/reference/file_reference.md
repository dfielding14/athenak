# Public File Reference

This page maps public documentation topics to source directories and shipped
inputs present on the public implementation baseline. It is a navigation aid,
not a promise that every file in a development branch is publicly supported.

## Runtime Core

| Path | Responsibility |
| --- | --- |
| `src/main.cpp` | CLI processing, input/restart loading, run-directory switch, execution setup |
| `src/parameter_input.cpp`, `src/parameter_input.hpp` | Input parsing and command-line parameter overrides |
| `src/driver/driver.cpp`, `src/driver/driver.hpp` | Time integrator and output/finalization loop |
| `src/mesh/` | Mesh, MeshBlock packs, refinement, load balancing, tree construction |
| `src/bvals/` | Boundary exchange, physical boundaries, prolongation and flux correction |
| `src/tasklist/` | Shared task-list facilities |

## Physics And Numerics

| Directory | Public contents | Related page |
| --- | --- | --- |
| `src/hydro/` | Hydro state, tasks, updates, timestep, FOFC and hydro solvers | [Hydrodynamics](../modules/hydro.md) |
| `src/mhd/` | MHD state, constrained transport, EMFs, tasks, FOFC and MHD solvers | [MHD](../modules/mhd.md) |
| `src/eos/` | Ideal/isothermal and relativistic equation-of-state implementations | [EOS](../modules/eos.md) |
| `src/reconstruct/` | `dc.hpp`, `plm.hpp`, `ppm.hpp`, `wenoz.hpp` | [Reconstruction](../modules/reconstruction.md) |
| `src/diffusion/` | Viscosity, conduction and resistivity | [Diffusion](../modules/diffusion.md) |
| `src/radiation/` | Radiation state, sources, fluxes, timestep and tetrads | [Radiation](../modules/radiation.md) |
| `src/ion-neutral/` | Ion-neutral coupling and tasks | [Ion-Neutral](../modules/ion_neutral.md) |
| `src/particles/` | Particle state, pushers and tasks | [Particles](../modules/particles.md) |
| `src/shearing_box/` | Shearing-box and orbital-advection implementation | [Shearing Box](../modules/shearing_box.md) |
| `src/coordinates/`, `src/units/` | Coordinates, ADM helpers, excision and units | [Coordinates](../modules/coordinates.md) |
| `src/z4c/`, `src/dyn_grmhd/` | Numerical-relativity and dynamic GRMHD implementation | [Z4c](../modules/z4c.md), [Dyn GRMHD](../modules/dyn_grmhd.md) |

## Sources, Outputs, And Problem Setup

| Path | Public role |
| --- | --- |
| `src/srcterms/srcterms.cpp`, `src/srcterms/srcterms_newdt.cpp`, `src/srcterms/ismcooling.hpp` | Public fluid/radiation source terms and associated timestep limits |
| `src/srcterms/turb_driver.cpp`, `src/srcterms/turb_driver.hpp` | Public turbulence driver |
| `src/outputs/` | Registered output types including table, history, VTK, binary, PDF and restart writers |
| `src/pgen/pgen.cpp`, `src/pgen/pgen.hpp` | Built-in dispatch and custom problem-generator interface |
| `src/pgen/tests/` | Built-in verification problem implementations |
| `src/pgen/*.cpp` | Source-selected custom problem generators configured with `-DPROBLEM=<stem>` |

The public source-terms directory does **not** contain a CGM cooling-table
implementation. Development records about such work are kept under
[Developer Notes](../engineering/index.md), outside supported public
configuration guidance.

## Shipped Input Families

| Directory | Examples of covered workflows |
| --- | --- |
| `inputs/hydro/` | Sod, implosion, blast, Kelvin-Helmholtz, Rayleigh-Taylor, turbulence |
| `inputs/mhd/` | Brio-Wu, Orszag-Tang, field loop, MHD blast and resistivity |
| `inputs/tests/` | Verification and convergence inputs for built-in generators |
| `inputs/radiation/` | Beam, diffusion, relaxation, shadow and hohlraum cases |
| `inputs/ion-neutral/` | Two-fluid blast and C-shock cases |
| `inputs/particles/` | Particle drift input |
| `inputs/shearing_box/` | Shearing-wave, epicycle and MRI inputs |
| `inputs/grhydro/`, `inputs/grmhd/`, `inputs/dyngr/`, `inputs/z4c/` | Relativistic and numerical-relativity workflows |

An input file may still require a custom build-selected problem generator.
Check [Problem Generators](../modules/pgen.md) and the example page before
assuming that an input runs with the default executable.

## Build, Analysis, And Documentation

| Path | Responsibility |
| --- | --- |
| `CMakeLists.txt`, `src/CMakeLists.txt` | Top-level options and executable source selection |
| `kokkos/` | Bundled Kokkos source configured by CMake |
| `vis/python/` | Public Python output readers, conversion and plotting helpers |
| `docs/source/` | Sphinx/MyST documentation source |
| `.github/workflows/docs.yml` | Documentation publishing workflow |

## Where To Start

| Goal | First files to inspect |
| --- | --- |
| Add or debug an initial condition | `src/pgen/pgen.hpp`, selected `src/pgen/*.cpp`, matching `inputs/` deck |
| Diagnose a shared configuration option | [Input Parameters](input_parameters.md), then the named owning source |
| Add an output format | `src/outputs/outputs.cpp`, `src/outputs/outputs.hpp` |
| Work on mesh refinement | `src/mesh/mesh.cpp`, `src/mesh/mesh_refinement.cpp`, `src/mesh/refinement_criteria.cpp` |
| Analyze public output | `vis/python/` and [Visualization Utilities](../tools/visualization.md) |
