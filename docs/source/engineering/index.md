# Developer Notes

These pages record implementation work, design investigations, and validation
status. They may describe branches or proposed work that is not present in the
public `main` implementation used by this documentation site. They are useful
when developing or reviewing AthenaK, but the
[Modules](../modules/index.md) and [Reference](../reference/api_reference.md)
sections should be used for stable interface guidance.

```{warning}
Do not treat an engineering record as a supported run workflow unless the same
capability is documented in the stable module/reference sections and is
present in the public implementation.
```

## Adaptive Mesh Refinement And Divergence Control

- [AMR Turbulence Implementation](amr_turbulence_implementation.md)
- [Divergence-Control AMR Completion Plan](divb_amr_completion_plan.md)
- [Divergence-Control AMR Final Status](divb_amr_final_status.md)
- [Divergence-Control AMR Implementation Plan](divb_amr_implementation_plan.md)
- [Divergence-Control AMR Implementation Summary](divb_amr_implementation_summary.md)
- [Face Field Correction Implementation](face_field_correction_implementation.md)

## Particle Implementation Records

- [Particle Merge And Cosmic-Ray Pushers](particle_merge_cr_pushers.md)
- [Particle Merge Implementation](particle_merge_implementation.md)
- [Particle Merge Unified Data](particle_merge_unified_data.md)

## Development-Only Problem Records

- [CGM Cooling Flow With Metals And Supernovae](../cgm_cooling_flow_metals.md)

```{toctree}
:hidden:
:maxdepth: 1

amr_turbulence_implementation
divb_amr_completion_plan
divb_amr_final_status
divb_amr_implementation_plan
divb_amr_implementation_summary
face_field_correction_implementation
particle_merge_cr_pushers
particle_merge_implementation
particle_merge_unified_data
../cgm_cooling_flow_metals
```
