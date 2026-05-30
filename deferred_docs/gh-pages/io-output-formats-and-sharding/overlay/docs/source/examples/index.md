# Worked Examples

These examples connect complete simulation setups to modules and analysis
workflows. The IO example documents output formats and readback utilities with
input decks exercised by the source-tree test suite.

| Example | Purpose | Related page |
| --- | --- | --- |
| [Shock Tube](shock_tube.md) | Hydro verification | [Hydrodynamics](../modules/hydro.md) |
| [Blast Wave](blast_wave.md) | Multidimensional shocks | [Outputs](../modules/outputs.md) |
| [Driven Turbulence](turbulence.md) | Driven flows | [Source Terms](../modules/srcterms.md) |
| [MRI Turbulence](mri_turbulence.md) | MHD instability | [MHD](../modules/mhd.md) |
| [Binary Merger](binary_merger.md) | Relativistic application | [DynGRMHD](../modules/dyn_grmhd.md) |
| [IO Outputs And Sharding](io_outputs_and_sharding.md) | PDFs, slices, node shards, restarts, and Python readback | [Outputs](../modules/outputs.md) |

Development records that are not runnable public examples are kept under
[Developer Notes](../engineering/index.md), separately from this selection path.

```{toctree}
:hidden:
:maxdepth: 1

shock_tube
blast_wave
turbulence
mri_turbulence
binary_merger
io_outputs_and_sharding
```
