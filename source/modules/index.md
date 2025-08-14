# Module Documentation

Comprehensive guides for AthenaK's major components and subsystems.

## Core Modules

::::{grid} 2
:gutter: 3

:::{grid-item-card} 🏗️ Mesh
:link: mesh.md

AMR infrastructure, domain decomposition, MeshBlocks
:::

:::{grid-item-card} ⏱️ Driver
:link: driver.md

Time integration, task execution, main loop
:::

:::{grid-item-card} 📋 Task System
:link: tasklist.md

Task-based parallelism, dependencies, scheduling
:::

:::{grid-item-card} 🌐 Coordinates
:link: coordinates.md

Coordinate systems, metrics, transformations
:::

::::

## Physics Modules

::::{grid} 2
:gutter: 3

:::{grid-item-card} 💧 Hydrodynamics
:link: hydro.md

Euler equations, shock capturing, conservation laws
:::

:::{grid-item-card} 🧲 MHD
:link: mhd.md

Magnetohydrodynamics, constrained transport, divergence cleaning
:::

:::{grid-item-card} ☢️ Radiation
:link: radiation.md

Radiation transport, moment methods, opacity
:::

:::{grid-item-card} ⚫ Relativity
:link: gr.md

GR hydro/MHD, numerical relativity (Z4c)
:::

:::{grid-item-card} 🎯 Particles
:link: particles.md

Lagrangian particles, test particles, cosmic rays
:::

::::

## Numerical Methods

::::{grid} 2
:gutter: 3

:::{grid-item-card} 📐 Reconstruction
:link: reconstruction.md

PLM, PPM, WENO-Z spatial reconstruction
:::

:::{grid-item-card} 🔀 Riemann Solvers
:link: riemann.md

HLLE, HLLC, HLLD, Roe solvers
:::

:::{grid-item-card} 🌡️ Equation of State
:link: eos.md

Ideal, isothermal, tabulated EOS
:::

:::{grid-item-card} 🚧 Boundary Conditions
:link: boundaries.md

Periodic, outflow, reflecting, user-defined
:::

::::

## Support Systems

::::{grid} 2
:gutter: 3

:::{grid-item-card} 💾 I/O System
:link: outputs.md

VTK, HDF5, restart files, history output
:::

:::{grid-item-card} 🌊 Source Terms
:link: sources.md

Turbulence driving, cooling, gravity
:::

:::{grid-item-card} 🔄 Diffusion
:link: diffusion.md

Viscosity, resistivity, thermal conduction
:::

:::{grid-item-card} 📦 Problem Generators
:link: pgen.md

Initial conditions, test problems, setups
:::

::::

## Quick Links

- **Getting Started**: [Quickstart Guide](../quickstart.md)
- **Architecture**: [System Overview](../flowcharts/runtime.md)
- **Migration**: [From Athena++](../migration/index.md)
- **API Reference**: [Complete API](../api/index.rst)