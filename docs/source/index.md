# AthenaK Documentation

Welcome to the AthenaK knowledge base. This site is organised to help you get productive quickly, understand the architecture, and dive into the full physics and infrastructure modules that power the code.

```{toctree}
:maxdepth: 1
:caption: Getting Started

quickstart
building
configuration
running
troubleshooting
```

```{toctree}
:maxdepth: 2
:caption: System Architecture

overview
flowcharts/system_architecture
flowcharts/runtime
kokkos_guide
```

```{toctree}
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

```{toctree}
:maxdepth: 2
:caption: Reference

reference/input_parameters
reference/file_reference
reference/api_reference
glossary
```

```{toctree}
:maxdepth: 2
:caption: Migration

migration/index
migration/common_gotchas
```

```{toctree}
:maxdepth: 1
:caption: Examples

examples/shock_tube
examples/blast_wave
examples/turbulence
examples/mri_turbulence
examples/binary_merger
cgm_cooling_flow_metals
```

```{toctree}
:maxdepth: 1
:caption: Contributing

contributing_docs
```

```{toctree}
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
