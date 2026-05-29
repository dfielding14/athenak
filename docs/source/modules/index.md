# AthenaK Modules Documentation

Stable module entry points for the public AthenaK implementation, organized by
category. Developer records for work that is not part of this baseline are
listed separately under [Developer Notes](../engineering/index.md).

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
| **Radiation** | `src/radiation/` | Angle-resolved intensity transport on a geodesic angular grid | [radiation.md](radiation.md) |
| **Z4c** | `src/z4c/` | Numerical relativity (Einstein equations) | [z4c.md](z4c.md) |
| **DynGRMHD** | `src/dyn_grmhd/` | GR MHD in dynamical spacetimes | [dyn_grmhd.md](dyn_grmhd.md) |
| **Ion-Neutral** | `src/ion-neutral/` | Two-fluid ion-neutral MHD | [ion_neutral.md](ion_neutral.md) |
| **Particles** | `src/particles/` | Current cosmic-ray drift particle path | [particles.md](particles.md) |

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
| **Outputs** | `src/outputs/` | Data I/O (12 registered formats) | [outputs.md](outputs.md) |
| **Boundary Values** | `src/bvals/` | Boundary conditions and MPI | [boundaries.md](boundaries.md) |
| **Source Terms** | `src/srcterms/` | External sources and turbulence | [srcterms.md](srcterms.md) |
| **Shearing Box** | `src/shearing_box/` | Shearing source terms and boundaries | [shearing_box.md](shearing_box.md) |
| **Problem Generators** | `src/pgen/` | Initial conditions | [pgen.md](pgen.md) |

## Development Records

Problem-specific records whose implementation is not present in the public
baseline, including the CGM cooling-flow note, are intentionally routed through
[Developer Notes](../engineering/index.md) rather than the stable module or
worked-example paths.

## Additional Components

| Component | Directory | Description |
|-----------|-----------|-------------|
| **Units** | `src/units/` | Physical unit system |
| **Utils** | `src/utils/` | Utility functions |
| **Geodesic Grid** | `src/geodesic-grid/` | Angular grids and spherical-surface interpolation |

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
    
    Main --> Driver[Driver]
    Driver --> Tasks[TaskList]
    Driver --> Outputs[Outputs]
    
    MB --> BVals[Boundary Values]
    BVals --> MPI[MPI Comm]
    
    style Main fill:#f9f,stroke:#333,stroke-width:4px
    style Mesh fill:#bbf,stroke:#333,stroke-width:2px
    style Physics fill:#fbf,stroke:#333,stroke-width:2px
```

## Quick Reference

### Extending Physics

New modules must be wired through the existing `MeshBlockPack` and task-list
patterns and documented only to the level verified by their implementation and
tests. Start with [Problem Generators](pgen.md) for problem-specific initial
conditions and [Migration](../migration/index.md) for porting guidance.

## Input Parameters

Shared public controls and owning-source pointers are documented in
[Input Parameters Reference](../reference/input_parameters.md); consult
module sources and shipped input decks for specialized settings.

## See Also
- [Architecture Overview](../flowcharts/runtime.md)
- [File Reference](../reference/file_reference.md)
- [Migration Guide](../migration/index.md)

```{toctree}
:hidden:
:maxdepth: 1

mesh
driver
tasklist
coordinates
hydro
mhd
radiation
z4c
dyn_grmhd
ion_neutral
particles
reconstruction
riemann_solvers
eos
diffusion
outputs
boundaries
srcterms
shearing_box
pgen
```
