# Example: Binary Neutron Star Initial Data

The public tree includes `inputs/dyngr/sgrid_bns_amr.athinput` and the custom
generator `src/pgen/sgrid_bns.cpp` for evolution initialized from an external
SGRID binary-neutron-star data set.

```{warning}
This is not a first-run example. It requires the external SGRID library and a
compatible initial-data directory, neither of which is shipped in this
repository. The documentation audit therefore verifies its configuration path
against source and input files, not a complete local simulation run.
```

## External Prerequisites

- Build or obtain a compatible `libsgrid.a` and place it at
  `sgrid/lib/libsgrid.a`; `src/pgen/sgrid_bns.cpp` requires SGRID to be built
  without OpenMP support.
- Provide the initial-data directory referenced by the input deck:
  `sgrid_26_26_26`.
- Review the substantial 192-cubed mesh and seven-level refinement setup
  before allocating production resources.

## Build Selection

The generator is selected at configure time. CMake links the external archive
only for this selection:

```bash
cmake -S . -B build-sgrid -DPROBLEM=sgrid_bns
cmake --build build-sgrid
./build-sgrid/src/athena -i inputs/dyngr/sgrid_bns_amr.athinput
```

The executable cannot be built from a fresh repository checkout until the
required SGRID archive has been supplied. A run additionally requires the
configured initial-data directory.

## Input And Source Contract

The shipped input enables MHD, ADM and Z4c evolution. It configures adaptive
refinement with a user criterion and sets `<problem>/user_hist = true`.
`src/pgen/sgrid_bns.cpp` enrolls the corresponding history and refinement
callbacks, loads `datadir`, and stages SGRID work beneath the default `SGRID`
output directory.

The input deck configures:

| Stream | Content |
| --- | --- |
| `output1` | History diagnostics |
| `output2` | Binary `mhd_w` slice |
| `output3` | Binary `mhd_u` slice |
| `output4` | Binary `z4c` slice |
| `output5` | Binary `adm` slice |
| `output6` | Restart data |

It also configures Z4c waveform-extraction radii in the `<z4c>` block.

## See Also

- [Z4c](../modules/z4c.md)
- [Dyn GRMHD](../modules/dyn_grmhd.md)
- [Outputs](../modules/outputs.md)
