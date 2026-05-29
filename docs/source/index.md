# AthenaK Documentation

Welcome to the AthenaK knowledge base. Choose a path below, or explore the
documentation by topic.

## Find Your Path

| Goal | Start here | Continue with |
| --- | --- | --- |
| Run AthenaK for the first time | [5-Minute Quickstart](quickstart.md) | [Running Simulations](running.md) |
| Configure or adapt a simulation | [Worked Examples](examples/index.md) | [Configuration Guide](configuration.md) |
| Understand or extend the code | [System Overview](overview.md) | [Modules](modules/index.md) and [Kokkos Guide](kokkos_guide.md) |
| Port code from Athena++ | [Migration Guide](migration/index.md) | [Common Gotchas](migration/common_gotchas.md) |
| Review ongoing implementation work | [Developer Notes](engineering/index.md) | Stable guidance in [Modules](modules/index.md) |

<div style="width: 100%; max-width: 900px; margin: 1em auto; border: 1px solid #ddd; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 8px rgba(0,0,0,0.1);">
<iframe src="_static/athenak_fluid_sim.html" title="Interactive AthenaK fluid-simulation banner" style="width: 100%; aspect-ratio: 3/1; border: none;"></iframe>
</div>

<details style="margin-top: 0.5em; font-size: 0.85em;">
<summary style="cursor: pointer; font-style: italic; color: #666;">More on this simulation...</summary>
<p style="max-width: 700px; margin: 0.5em 0; color: #444;">
This is a real-time 2D incompressible Navier-Stokes simulation running in WebGL2. The velocity field is evolved by tracing fluid parcels backward in time to find where they originated (semi-Lagrangian advection), then smoothing the velocity to model viscous drag (implicit diffusion solves for the smoothed state directly, which is numerically stable). Pressure projection via Jacobi iteration enforces the divergence-free constraint, ensuring mass conservation. A passive scalar field (the visible "dye") is carried along by the flow and subject to source/sink terms: it grows where the underlying "AthenaK" text mask is defined and decays elsewhere, causing the letters to continuously reform after being disturbed.
</p>
</details>

<div style="margin-top: 2em;"></div>

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

## Modules

Use the [Modules index](modules/index.md) to browse infrastructure, physics,
numerical methods, and support systems without expanding the complete
catalogue in the global navigation.

```{toctree}
:hidden:
:maxdepth: 1
:caption: Modules

modules/index
```

## Reference Library

- [Input Parameters](reference/input_parameters.md) - verified shared controls and source pointers
- [File Reference](reference/file_reference.md) - public source and input-deck map
- [API Reference](reference/api_reference.md) - implemented developer-facing interfaces

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

Start from the [Worked Examples index](examples/index.md) to find a test
problem by physics goal and connect it to configuration and output guidance.

```{toctree}
:hidden:
:maxdepth: 1
:caption: Examples

examples/index
```

## Tools & Utilities

- [Visualization Toolkit](tools/visualization.md) - readers and converters for current AthenaK outputs

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

## Developer Notes

Implementation records and ongoing investigations are collected separately
from stable module and user guidance in [Developer Notes](engineering/index.md).

```{toctree}
:hidden:
:maxdepth: 1
:caption: Developer Notes

engineering/index
```

Need to fix something or found a discrepancy? Please log it in `documentation_audit_log.md` and submit a pull request!
