# AthenaK System Overview

AthenaK is a Kokkos-based astrophysical simulation code. The public source
baseline documented here contains hydrodynamics, MHD, radiation intensity
transport, numerical-relativity components, source terms, particles, mesh
refinement, output writers, and problem generators. Not every possible
combination is valid; module pages document verified compatibility constraints.

## Execution Model

```{mermaid}
flowchart LR
    INPUT["Input deck and CLI overrides"] --> MESH["Mesh / MeshBlockPack"]
    MESH --> MODS["Coordinates and enabled physics"]
    MODS --> PGEN["ProblemGenerator"]
    PGEN --> DRIVER["Driver and task lists"]
    DRIVER --> OUTPUTS["Configured outputs"]
    DRIVER --> AMR["Optional mesh refinement"]
    AMR --> DRIVER
```

The startup and timestep details are mapped to implementation entry points in
the [Runtime Flowchart](flowcharts/runtime.md).

## Choose A Starting Calculation

| Goal | Verified starting point | Build selection |
| --- | --- | --- |
| First hydro run | `inputs/hydro/sod.athinput` | default executable |
| First MHD visualization run | `inputs/mhd/orszag_tang.athinput` | default executable |
| Driven hydro turbulence | `inputs/hydro/turb.athinput` | `-DPROBLEM=turb` |
| Hydro blast example | `inputs/hydro/blast_hydro.athinput` | `-DPROBLEM=blast` |
| 3D shearing-box MRI | `inputs/shearing_box/mri3d_unstratified.athinput` | default executable |

Start with the [Quickstart](quickstart.md) for commands and observed artifacts,
then consult [Examples](examples/index.md) for build requirements and known
limitations.

## Main Public Components

| Component | Source area | Documented entry point |
| --- | --- | --- |
| Mesh, MeshBlocks, refinement | `src/mesh/` | [Mesh](modules/mesh.md) |
| Time integration | `src/driver/` | [Driver](modules/driver.md) |
| Hydro / MHD | `src/hydro/`, `src/mhd/` | [Hydro](modules/hydro.md), [MHD](modules/mhd.md) |
| Radiation | `src/radiation/`, `src/geodesic-grid/` | [Radiation](modules/radiation.md) |
| Relativity | `src/z4c/`, `src/dyn_grmhd/` | [Z4c](modules/z4c.md), [DynGRMHD](modules/dyn_grmhd.md) |
| Outputs | `src/outputs/` | [Outputs](modules/outputs.md) |
| Initial conditions | `src/pgen/` | [Problem Generators](modules/pgen.md) |

## Configuration Boundaries

- `<hydro>` and `<mhd>` select fluid modules; using both requires the
  ion-neutral coupling path.
- Coordinate settings select non-relativistic, SR, or GR fluid paths where
  implemented; relativistic hydro/MHD reject `eos = isothermal`.
- Radiation stores angle-resolved intensity on a geodesic angular grid. It is
  not documented as an M1-closure module. Its required and optional parameters
  are listed in [Radiation](modules/radiation.md).
- Refinement settings belong under `<mesh_refinement>`; criteria are configured
  with `<amr_criterion*>` blocks as documented in [Mesh](modules/mesh.md).

## For Developers

The central developer interfaces are the Kokkos array aliases and loop helpers
in `src/athena.hpp`, the pack/module assembly in
`src/mesh/meshblock_pack.cpp`, the driver loop in `src/driver/driver.cpp`, and
problem-generator hooks under `src/pgen/`. Use the [Kokkos Guide](kokkos_guide.md),
[API Reference](reference/api_reference.md), and
[Migration Guide](migration/index.md) before extending code.

Implementation records under [Developer Notes](engineering/index.md) are
deliberately separated from public supported guidance; some describe branch
work not present in this public source baseline.
