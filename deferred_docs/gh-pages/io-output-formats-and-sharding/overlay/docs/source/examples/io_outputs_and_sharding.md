# IO Outputs And Sharding

This example suite demonstrates the analysis-oriented IO formats and
distribution choices implemented in AthenaK. The input decks live in
`inputs/io/` in the code repository and are exercised by regression tests.

## Build And Run A Shared-Output Example

Build AthenaK using your normal configuration, then run:

```bash
./build/src/athena \
  -i inputs/io/output_formats.athinput \
  -d run_io_formats
```

This deck demonstrates:

- modern PDF output;
- binary `sphslice` output; and
- coexistence with the existing `file_type = sph` output.

Read the generated products with the public example helper:

```bash
python vis/python/examples/read_io_outputs.py \
  pdf run_io_formats/pdf_nd3_coord_abscostheta_vel_sph_r/io_formats.00000.pdf

python vis/python/examples/read_io_outputs.py \
  sphslice run_io_formats/bin/io_formats.density.r_0.25.00000.sph.bin
```

Use the exact emitted paths from the run directory if the example basename or
stream ID differs after local overrides.

## Four-Dimensional And Scalar-Weighted PDFs

```bash
./build/src/athena \
  -i inputs/io/output_pdf_scalar_weight.athinput \
  -d run_io_weighted_pdf

python vis/python/examples/read_io_outputs.py \
  pdf run_io_weighted_pdf/pdf_nd4_coord_r_hydro_w_s_0_vel_sph_r/io_pdf_scalar.00000.pdf
```

This input uses the modern PDF syntax:

- `variable_1` through `variable_4` for contiguous scalar axes;
- `scaleN = linear`, `log`, or `symlog`;
- `linthreshN` for a symlog axis; and
- `weight = variable` plus `weight_variable` for an integrated scalar weight.

The test suite also checks rejected invalid scale combinations.

## Output Timing And Final Policy

```bash
./build/src/athena \
  -i inputs/io/runtime_policy.athinput \
  -d run_io_policy \
  time/output_timing=true \
  time/final_output_policy=restart_only
```

The input enables:

```ini
<time>
output_timing = true
final_output_policy = restart_only
```

It produces timing lines shaped like:

```text
[output-io] event=final block=output2 type=rst distribution=shared elapsed_max_s=...
```

`restart_only` allows a final restart checkpoint while suppressing additional
terminal analysis outputs. The default policy, `all`, retains normal
final-output behavior.

## MPI Node-Sharded Output

Use an MPI build:

```bash
mpirun -np 2 ./build-mpi/src/athena \
  -i inputs/io/node_sharded_outputs.athinput \
  -d run_io_node
```

The input enables `single_file_per_node = true` for node-sharded binary,
full-volume coarsened-binary, PDF, spherical-slice, and restart streams.
Files are organized below `node_########/` directories. For example:

```text
run_io_node/bin/node_00000000/
run_io_node/pdf_radius_density_hydro_w_d/node_00000000/
run_io_node/rst/node_00000000/
```

Read a node-sharded analysis product by passing one shard and requesting
assembly:

```bash
python vis/python/examples/read_io_outputs.py \
  bin run_io_node/bin/node_00000000/io_node_example.density.00000.bin \
  --assemble-shards

python vis/python/examples/read_io_outputs.py \
  cbin run_io_node/cbin_coarse_density_2/node_00000000/io_node_example.coarse_density.00000.cbin \
  --assemble-shards

python vis/python/examples/read_io_outputs.py \
  pdf run_io_node/pdf_radius_density_hydro_w_d/node_00000000/io_node_example.00000.pdf

python vis/python/examples/read_io_outputs.py \
  sphslice run_io_node/bin/node_00000000/io_node_example.density.r_0.25.00000.sph.bin
```

If local overrides change IDs, basenames, or output numbering, use the
corresponding emitted paths in these readback commands.

## Restart From The Public Manifest

Node-sharded restart output creates:

```text
run_io_node/rst/<basename>.<number>.rst
run_io_node/rst/node_00000000/<basename>.<number>.g<generation>.payload.rst
```

Restart from the `.rst` manifest at the root of `rst/`:

```bash
mpirun -np 2 ./build-mpi/src/athena \
  -r run_io_node/rst/<basename>.<number>.rst \
  -d run_io_resumed
```

Do not point `-r` at a node payload. The runtime validates manifest
completeness and payload metadata before reconstructing the checkpoint.

## Compatibility Boundaries

- Existing unsharded old-style PDF inputs preserve legacy output bytes; new
  PDF configurations should use the modern `variable_N`, `scaleN`, and
  `weight` interface.
- `sphslice` is not a replacement for the existing `file_type = sph` output.
- Node-sharded full-volume `cbin` is covered by this example; sliced
  node-sharded `cbin` is outside the promoted workflow pending a separate
  legacy slice/readback repair.
- The automated MPI example can run multiple ranks on one physical node.
  Before production use, qualify node sharding and node restart on a real
  multi-node allocation and the target parallel filesystem.
