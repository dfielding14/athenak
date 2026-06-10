# File Reference

Complete reference of all source files in AthenaK organized by directory.

## Quick Navigation

| Directory | Purpose | Key Files |
|-----------|---------|-----------|
| `src/` | Main source code | `main.cpp`, `athena.hpp` |
| `src/mesh/` | Mesh infrastructure | `mesh.cpp`, `meshblock.cpp` |
| `src/hydro/` | Hydrodynamics | `hydro.cpp`, `hydro_fluxes.cpp` |
| `src/mhd/` | Magnetohydrodynamics | `mhd.cpp`, `mhd_ct.cpp` |
| `src/pgen/` | Problem generators | 50+ problem setups |
| `inputs/` | Input files | Example configurations |

## Core Files (`src/`)

### Main Entry Point
- **`main.cpp`** - Program entry, initialization
- **`athena.hpp`** - Global definitions, constants
- **`globals.cpp/hpp`** - Global variables

### Parameter Management
- **`parameter_input.cpp/hpp`** - Input file parsing

## Mesh Infrastructure (`src/mesh/`)

### Core Mesh
- **`mesh.cpp/hpp`** - Main mesh class [Details](../modules/mesh.md)
- **`meshblock.cpp/hpp`** - MeshBlock implementation
- **`meshblock_pack.cpp/hpp`** - GPU optimization layer
- **`meshblock_tree.cpp/hpp`** - AMR tree structure

### AMR
- **`mesh_refinement.cpp/hpp`** - Refinement algorithms
- **`build_tree.cpp`** - Tree construction
- **`load_balance.cpp`** - MPI load balancing

## Physics Modules

### Hydrodynamics (`src/hydro/`)
- **`hydro.cpp/hpp`** - Main hydro class [Details](../modules/hydro.md)
- **`hydro_fluxes.cpp`** - Flux calculation
- **`hydro_update.cpp`** - Conservative update
- **`hydro_tasks.cpp`** - Task registration
- **`hydro_newdt.cpp`** - Timestep calculation

### MHD (`src/mhd/`)
- **`mhd.cpp/hpp`** - Main MHD class [Details](../modules/mhd.md)
- **`mhd_ct.cpp`** - Constrained transport
- **`mhd_corner_e.cpp`** - Corner EMF
- **`mhd_fluxes.cpp`** - MHD flux calculation

### Radiation (`src/radiation/`)
- **`radiation.cpp/hpp`** - Main radiation class [Details](../modules/radiation.md)
- **`radiation_fluxes.cpp`** - Radiation flux
- **`radiation_source.cpp`** - Source terms
- **`radiation_tetrad.cpp/hpp`** - GR tetrad formalism

### General Relativity (`src/z4c/`)
- **`z4c.cpp/hpp`** - Z4c evolution [Details](../modules/z4c.md)
- **`z4c_calcrhs.cpp`** - Right-hand side
- **`z4c_gauge.cpp`** - Gauge conditions
- **`z4c_amr.cpp`** - AMR criteria
- **`tmunu.cpp/hpp`** - Stress-energy tensor

### GRMHD (`src/dyn_grmhd/`)
- **`dyn_grmhd.cpp/hpp`** - Dynamical GRMHD [Details](../modules/dyn_grmhd.md)
- **`dyn_grmhd_fluxes.cpp`** - GRMHD fluxes
- **`dyn_grmhd_fofc.cpp`** - First-order flux correction

### Particles (`src/particles/`)
- **`particles.cpp/hpp`** - Particle management [Details](../modules/particles.md)
- **`particles_pushers.cpp`** - Integration algorithms
- **`particles_tasks.cpp`** - Particle tasks

## Numerical Methods

### Reconstruction (`src/reconstruct/`)
- **`dc.hpp`** - Donor cell (1st order)
- **`plm.hpp`** - Piecewise linear (2nd order)
- **`ppm.hpp`** - Piecewise parabolic (3rd order)
- **`wenoz.hpp`** - WENO-Z (5th order)

### Riemann Solvers

#### Hydro (`src/hydro/rsolvers/`)
- **`llf_hyd.hpp`** - Local Lax-Friedrichs
- **`hlle_hyd.hpp`** - HLLE solver
- **`hllc_hyd.hpp`** - HLLC solver
- **`roe_hyd.hpp`** - Roe solver

#### MHD (`src/mhd/rsolvers/`)
- **`llf_mhd.hpp`** - LLF for MHD
- **`hlle_mhd.hpp`** - HLLE for MHD
- **`hlld_mhd.hpp`** - HLLD solver

### Equations of State (`src/eos/`)
- **`eos.cpp/hpp`** - Base EOS class [Details](../modules/eos.md)
- **`ideal_hyd.cpp`** - Ideal gas hydro
- **`ideal_mhd.cpp`** - Ideal gas MHD
- **`isothermal_*.cpp`** - Isothermal EOS

## Support Systems

### Boundaries (`src/bvals/`)
- **`bvals.cpp/hpp`** - Boundary conditions [Details](../modules/boundaries.md)
- **`bvals_cc.cpp`** - Cell-centered BCs
- **`bvals_fc.cpp`** - Face-centered BCs
- **`physics/*.cpp`** - Physics-specific BCs

### Outputs (`src/outputs/`)
- **`outputs.cpp/hpp`** - Output manager [Details](../modules/outputs.md)
- **`vtk_mesh.cpp`** - VTK output
- **`binary.cpp`** - Binary output
- **`restart.cpp`** - Restart files
- **`history.cpp`** - Time series

### Source Terms (`src/srcterms/`)
- **`srcterms.cpp/hpp`** - Source term manager [Details](../modules/srcterms.md)
- **`turb_driver.cpp/hpp`** - Turbulence driving
- **`cooling_tables.hpp`** - Cooling functions

### Coordinates (`src/coordinates/`)
- **`coordinates.cpp/hpp`** - Coordinate systems [Details](../modules/coordinates.md)
- **`cartesian_ks.hpp`** - Kerr-Schild
- **`adm.cpp/hpp`** - ADM formalism

### Diffusion (`src/diffusion/`)
- **`viscosity.cpp/hpp`** - Viscous diffusion
- **`conduction.cpp/hpp`** - Thermal conduction
- **`resistivity.cpp/hpp`** - Ohmic resistivity

### Task System (`src/tasklist/`)
- **`task_list.hpp`** - Task management [Details](../modules/tasklist.md)
- **`numerical_relativity.cpp/hpp`** - NR tasks

### Driver (`src/driver/`)
- **`driver.cpp/hpp`** - Main evolution loop [Details](../modules/driver.md)

### Utilities (`src/utils/`)
- **`utils.hpp`** - Utility functions
- **`random.hpp`** - Random number generation
- **`lagrange_interpolator.cpp/hpp`** - Interpolation

## Problem Generators (`src/pgen/`)

### Test Problems (`src/pgen/tests/`)
- **`linear_wave.cpp`** - Linear wave convergence
- **`shock_tube.cpp`** - Sod shock tube
- **`orszag_tang.cpp`** - Orszag-Tang vortex
- **`cpaw.cpp`** - Circularly polarized Alfven wave

### Astrophysical Problems
- **`blast.cpp`** - Blast waves
- **`turb.cpp`** - Driven turbulence
- **`gr_torus.cpp`** - Relativistic torus
- **`mri3d.cpp`** - 3D MRI
- **`dyngr_tov.cpp`** - TOV neutron star

[See all 50+ problem generators](../modules/pgen.md)

## Input Files (`inputs/`)

### Directory Structure
```
inputs/
├── tests/          # Test problems
├── hydro/          # Hydro examples
├── mhd/            # MHD examples
├── gr/             # General relativity
├── radiation/      # Radiation problems
└── particles/      # Particle examples
```

### Key Input Files
- **`tests/linear_wave_hydro.athinput`** - Hydro convergence test
- **`mhd/mri.athinput`** - MRI setup
- **`gr/torus.athinput`** - Black hole accretion

## Build System

### CMake Files
- **`CMakeLists.txt`** - Main build configuration
- **`src/CMakeLists.txt`** - Source compilation

### Configuration
- **`cmake/`** - CMake modules
- **`kokkos/`** - Kokkos submodule

## Documentation Files

### This Documentation
- **`docs/source/`** - Sphinx source files
- **`docs/source/modules/`** - Module documentation
- **`docs/source/reference/`** - Reference guides

### Key Documentation
- **`README.md`** - Project overview
- **`CLAUDE.md`** - AI assistant guide
- **`LICENSE`** - License information

## File Naming Conventions

### Source Files
- **`.cpp`** - Implementation files
- **`.hpp`** - Header files
- **`_tasks.cpp`** - Task registration
- **`_fluxes.cpp`** - Flux computation
- **`_newdt.cpp`** - Timestep calculation

### Problem Generators
- Must define `ProblemGenerator()` function
- File name becomes build flag: `-DPROBLEM=filename`

### Input Files
- **`.athinput`** - Configuration files
- Organized by physics type

## Quick File Lookup

### By Functionality

| Need to modify... | Look in... |
|-------------------|------------|
| Initial conditions | `src/pgen/` |
| Boundary conditions | `src/bvals/physics/` |
| Output format | `src/outputs/` |
| Riemann solver | `src/*/rsolvers/` |
| Time integration | `src/driver/` |
| Mesh refinement | `src/mesh/mesh_refinement.cpp` |

### By Physics

| Physics | Primary Files |
|---------|---------------|
| Hydro | `src/hydro/hydro*.cpp` |
| MHD | `src/mhd/mhd*.cpp` |
| Radiation | `src/radiation/radiation*.cpp` |
| Relativity | `src/z4c/z4c*.cpp` |
| Particles | `src/particles/particles*.cpp` |

## See Also
- [Module Index](../modules/index.md) - Detailed module documentation
- [Input Parameters](input_parameters.md) - All configuration options
- [API Reference](api_reference.md) - Function documentation