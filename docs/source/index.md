# AthenaK Documentation

<div style="width: 100%; max-width: 900px; margin: 1em auto; border: 1px solid #ddd; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 8px rgba(0,0,0,0.1);">
<iframe src="_static/athenak_fluid_sim.html" style="width: 100%; aspect-ratio: 3/1; border: none;"></iframe>
</div>

<details style="margin-top: 0.5em; font-size: 0.85em;">
<summary style="cursor: pointer; font-style: italic; color: #666;">More on this simulation...</summary>
<p style="max-width: 700px; margin: 0.5em 0; color: #444;">
This is a real-time 2D incompressible Navier-Stokes simulation running in WebGL2. The velocity field is evolved by tracing fluid parcels backward in time to find where they originated (semi-Lagrangian advection), then smoothing the velocity to model viscous drag (implicit diffusion solves for the smoothed state directly, which is numerically stable). Pressure projection via Jacobi iteration enforces the divergence-free constraint, ensuring mass conservation. A passive scalar field (the visible "dye") is carried along by the flow and subject to source/sink terms: it grows where the underlying "AthenaK" text mask is defined and decays elsewhere, causing the letters to continuously reform after being disturbed.
</p>
</details>

<div style="margin-top: 2em;"></div>

Welcome to the AthenaK knowledge base. Use this hub to launch simulations quickly, explore the architecture, or deep-dive into physics modules and engineering notes.

```{admonition} New to AthenaK?
:class: tip
Start with the **Getting Started** playbook below, then follow the quick links for your role (simulation author, performance tuner, or contributor).
```

## Getting Started

| Goal | Resource |
| --- | --- |
| First successful run in minutes | [5-Minute Quickstart](quickstart.md) |
| Build with the right toolchain | [Building Guide](building.md) |
| Configure inputs confidently | [Configuration Guide](configuration.md) |
| Learn runtime CLI + workflows | [Running Simulations](running.md) |
| Unblock common issues fast | [Troubleshooting Guide](troubleshooting.md) |

```{toctree}
:hidden:
:maxdepth: 1
:caption: Getting Started

quickstart
building
configuration
running
troubleshooting
```

## Architecture & Execution

- [System Overview](overview.md) — mental model of the codebase
- [Architecture Flowchart](flowcharts/system_architecture.md) — component relationships
- [Runtime Flowchart](flowcharts/runtime.md) — timestep/task sequencing
- [Kokkos Guide](kokkos_guide.md) — performance portability patterns

```{toctree}
:hidden:
:maxdepth: 2
:caption: System Architecture

overview
flowcharts/system_architecture
flowcharts/runtime
kokkos_guide
```

## Module Reference

Navigate the full stack of physics and infrastructure modules:

- Core infrastructure: [Mesh](modules/mesh.md), [Driver](modules/driver.md), [Task Lists](modules/tasklist.md), [Coordinates](modules/coordinates.md)
- Physics: [Hydro](modules/hydro.md), [MHD](modules/mhd.md), [Radiation](modules/radiation.md), [Z4c](modules/z4c.md), [Dyn GRMHD](modules/dyn_grmhd.md), [Ion-Neutral](modules/ion_neutral.md), [Particles](modules/particles.md)
- Numerical methods: [Reconstruction](modules/reconstruction.md), [Riemann Solvers](modules/riemann_solvers.md), [EOS](modules/eos.md), [Diffusion](modules/diffusion.md), [Outputs](modules/outputs.md), [Boundaries](modules/boundaries.md), [Source Terms](modules/srcterms.md), [Shearing Box](modules/shearing_box.md), [Problem Generators](modules/pgen.md)

```{toctree}
:hidden:
:maxdepth: 2
:caption: Modules

modules/index
modules/mesh
modules/driver
modules/tasklist
modules/coordinates
modules/hydro
modules/mhd
modules/radiation
modules/z4c
modules/dyn_grmhd
modules/ion_neutral
modules/particles
modules/reconstruction
modules/riemann_solvers
modules/eos
modules/diffusion
modules/outputs
modules/boundaries
modules/srcterms
modules/shearing_box
modules/pgen
```

## Reference Library

- [Input Parameters](reference/input_parameters.md) — authoritative block/parameter catalogue
- [File Reference](reference/file_reference.md) — generated API & file index
- [API Reference](reference/api_reference.md) — autodoc coverage

```{toctree}
:hidden:
:maxdepth: 2
:caption: Reference

reference/input_parameters
reference/file_reference
reference/api_reference
```

## Migration Guides

Planning a port from Athena++? Start here:

- [Migration Index](migration/index.md)
- [Common Gotchas](migration/common_gotchas.md)

```{toctree}
:hidden:
:maxdepth: 2
:caption: Migration

migration/index
migration/common_gotchas
```

## Worked Examples

Compare against shipped problem generators and diagnostics:

- [Shock Tube](examples/shock_tube.md)
- [Blast Wave](examples/blast_wave.md)
- [Driven Turbulence](examples/turbulence.md)
- [MRI Turbulence](examples/mri_turbulence.md)
- [Binary Merger](examples/binary_merger.md)
- [CGM Cooling Flow (Metals)](cgm_cooling_flow_metals.md)

```{toctree}
:hidden:
:maxdepth: 1
:caption: Examples

examples/shock_tube
examples/blast_wave
examples/turbulence
examples/mri_turbulence
examples/binary_merger
cgm_cooling_flow_metals
```

## Tools & Utilities

- [Visualization Toolkit](tools/visualization.md) — binary readers, converters, and plotting helpers for Athena++ outputs

```{toctree}
:hidden:
:maxdepth: 1
:caption: Tools & Utilities

tools/visualization
```

## Contributor Toolkit

- [Contributing Docs](contributing_docs.md) — style and workflow expectations

```{toctree}
:hidden:
:maxdepth: 1
:caption: Contributing

contributing_docs
```

## Engineering Notes

Internal design dossiers and ongoing investigations live here:

- AMR turbulence evolution, divergence control, particle merging, and more.

```{toctree}
:hidden:
:maxdepth: 1
:caption: Engineering Notes

engineering/amr_turbulence_implementation
engineering/divb_amr_completion_plan
engineering/divb_amr_final_status
engineering/divb_amr_implementation_plan
engineering/divb_amr_implementation_summary
engineering/face_field_correction_implementation
engineering/particle_merge_cr_pushers
engineering/particle_merge_implementation
engineering/particle_merge_unified_data
```

Need to fix something or found a discrepancy? Please log it in `documentation_audit_log.md` and submit a pull request!
