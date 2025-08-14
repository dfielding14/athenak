# AthenaK Modules Documentation

Complete documentation for all AthenaK modules, organized by category.

## Core Infrastructure (4 modules)

| Module | Directory | Description | Documentation |
|--------|-----------|-------------|---------------|
| **Mesh** | `src/mesh/` | AMR infrastructure, domain decomposition | [mesh.md](mesh.md) |
| **Driver** | `src/driver/` | Time integration and evolution control | [driver.md](driver.md) |
| **TaskList** | `src/tasklist/` | Task-based execution system | [tasklist.md](tasklist.md) |
| **Coordinates** | `src/coordinates/` | Coordinate systems and metrics | [coordinates.md](coordinates.md) |

## Physics Modules (7 modules)

| Module | Directory | Description | Documentation |
|--------|-----------|-------------|---------------|
| **Hydrodynamics** | `src/hydro/` | Euler equations solver | [hydro.md](hydro.md) |
| **MHD** | `src/mhd/` | Magnetohydrodynamics with CT | [mhd.md](mhd.md) |
| **Radiation** | `src/radiation/` | Radiation transport (M1 closure) | [radiation.md](radiation.md) |
| **Z4c** | `src/z4c/` | Numerical relativity (Einstein equations) | [z4c.md](z4c.md) |
| **DynGRMHD** | `src/dyn_grmhd/` | GR MHD in dynamical spacetimes | [dyn_grmhd.md](dyn_grmhd.md) |
| **Ion-Neutral** | `src/ion-neutral/` | Two-fluid ion-neutral MHD | [ion_neutral.md](ion_neutral.md) |
| **Particles** | `src/particles/` | Lagrangian particle tracking | [particles.md](particles.md) |

## Numerical Methods (4 modules)

| Module | Directory | Description | Documentation |
|--------|-----------|-------------|---------------|
| **Reconstruction** | `src/reconstruct/` | Spatial reconstruction (PLM, PPM, WENOZ) | [reconstruction.md](reconstruction.md) |
| **Riemann Solvers** | `src/*/rsolvers/` | Flux computation algorithms | [riemann_solvers.md](riemann_solvers.md) |
| **EOS** | `src/eos/` | Equations of state | [eos.md](eos.md) |
| **Diffusion** | `src/diffusion/` | Physical diffusion processes | [diffusion.md](diffusion.md) |

## Support Systems (5+ modules)

| Module | Directory | Description | Documentation |
|--------|-----------|-------------|---------------|
| **Outputs** | `src/outputs/` | Data I/O (14 formats) | [outputs.md](outputs.md) |
| **Boundary Values** | `src/bvals/` | Boundary conditions and MPI | [boundaries.md](boundaries.md) |
| **Source Terms** | `src/srcterms/` | External sources and turbulence | [srcterms.md](srcterms.md) |
| **Shearing Box** | `src/shearing_box/` | Orbital advection | [shearing_box.md](shearing_box.md) |
| **Problem Generators** | `src/pgen/` | Initial conditions | [pgen.md](pgen.md) |

## Additional Components

| Component | Directory | Description |
|-----------|-----------|-------------|
| **Units** | `src/units/` | Physical unit system |
| **Utils** | `src/utils/` | Utility functions |
| **Geodesic Grid** | `src/geodesic-grid/` | Spherical harmonic analysis |

## Module Interaction Diagram

```{mermaid}
flowchart TD
    Main[main.cpp] --> Mesh[Mesh Module]
    Mesh --> MB[MeshBlocks]
    MB --> Coords[Coordinates]
    
    MB --> Physics{Physics Modules}
    Physics --> Hydro[Hydro]
    Physics --> MHD[MHD]
    Physics --> Rad[Radiation]
    Physics --> Particles[Particles]
    
    Hydro --> Recon[Reconstruction]
    MHD --> Recon
    Hydro --> Riemann[Riemann Solvers]
    MHD --> Riemann
    
    Hydro --> EOS[EOS]
    MHD --> EOS
    
    MB --> Driver[Driver]
    Driver --> Tasks[TaskList]
    Tasks --> Outputs[Outputs]
    
    MB --> BVals[Boundary Values]
    BVals --> MPI[MPI Comm]
    
    style Main fill:#f9f,stroke:#333,stroke-width:4px
    style Mesh fill:#bbf,stroke:#333,stroke-width:2px
    style Physics fill:#fbf,stroke:#333,stroke-width:2px
```

## Quick Reference

### Adding a New Physics Module
1. Create directory in `src/`
2. Inherit from base physics class
3. Register with MeshBlock
4. Add tasks to TaskList
5. Document in `docs/source/modules/`

### Module Requirements
Each module must:
- Have clear interfaces with other modules
- Support device (GPU) execution via Kokkos
- Handle AMR refinement properly
- Provide restart capability
- Include unit tests

## Input Parameters

All module parameters are documented in [Input Parameters Reference](../reference/input_parameters.md).

## See Also
- [Architecture Overview](../flowcharts/runtime.md)
- [File Reference](../reference/file_reference.md)
- [Migration Guide](../migration/from_athena_plus_plus.md)