# Running Simulations

Commands below assume AthenaK is built under `build/`; substitute the path to
an MPI or accelerator build where needed.

## Basic Execution

```bash
./build/src/athena -i inputs/hydro/sod.athinput
./build/src/athena -i inputs/hydro/sod.athinput mesh/nx1=512 time/tlim=1.0
```

Validate an input without running:

```bash
./build/src/athena -i inputs/io/output_formats.athinput -n
```

## MPI Execution

```bash
mpirun -np 16 ./build-mpi/src/athena -i inputs/hydro/sod.athinput
```

The MPI build must enable the AthenaK and Kokkos MPI options. For node-sharded
IO, ranks on each shared-memory node cooperate to write that node's shard.

```bash
mpirun -np 16 ./build-mpi/src/athena \
  -i inputs/io/node_sharded_outputs.athinput \
  -d run_node_io
```

Before depending on node-sharded checkpoints in production, run this workflow
across more than one physical node and verify readback and restart.

## Command-Line Options

| Flag | Meaning |
| --- | --- |
| `-i <file>` | Input file for a new run. |
| `-r <file>` | Restart from a checkpoint or public node manifest. |
| `-d <dir>` | Runtime/output root directory. |
| `-n` | Parse input and exit. |
| `-m` | Write mesh structure information and exit. |
| `-c` | Show compiled configuration and exit. |
| `-t hh:mm:ss` | Wall-clock time limit; final behavior follows `time/final_output_policy`. |
| `block/name=value` | Override an input parameter after parsing. |

## Output Distribution And Readback

Use `single_file_per_rank = true` for rank shards or
`single_file_per_node = true` for node shards. Do not set both in one stream.

| Product | Shared path form | Node-sharded path form | Readback |
| --- | --- | --- | --- |
| Binary | `bin/<name>.bin` | `bin/node_########/<name>.bin` | `bin_convert.py` |
| Coarsened binary | `cbin_<id>_<factor>/<name>.cbin` | `cbin_<id>_<factor>/node_########/<name>.cbin` | `bin_convert.py` |
| PDF V2 | `pdf_<id>*/<name>.pdf` | `pdf_<id>*/node_########/<name>.pdf` | `read_pdf.py` |
| Spherical slice | `bin/<name>.sph.bin` | `bin/node_########/<name>.sph.bin` | `read_sphslice.py` |
| Restart | `rst/<name>.rst` | public manifest plus `rst/node_########/*.payload.rst` | `athena -r` |

Example Python readback:

```bash
python vis/python/examples/read_io_outputs.py \
  bin run_node_io/bin/node_00000000/io_node_example.density.00000.bin \
  --assemble-shards

python vis/python/examples/read_io_outputs.py \
  pdf run_node_io/pdf_radius_density_hydro_w_d/node_00000000/io_node_example.00000.pdf
```

For detailed public APIs, see the Visualization Utilities page.

## Restarting A Node-Sharded Job

With `file_type = rst` and `single_file_per_node = true`, AthenaK writes a
public manifest plus node payload shards. Use the manifest in `-r`:

```bash
mpirun -np 16 ./build-mpi/src/athena \
  -r run_node_io/rst/simulation.00001.rst \
  -d continued_run
```

Do not restart from `rst/node_########/*.payload.rst`. The public manifest is
validated for complete ordered node coverage, matching generations, payload
sizes, and completion before its state is consumed.

## Timing And Final Writes

Enable per-output timing and choose final writes in the input:

```ini
<time>
output_timing = true
final_output_policy = restart_only
```

Timing output is printed once for each written stream:

```text
[output-io] event=scheduled block=output_pdf type=pdf distribution=node elapsed_max_s=...
```

For MPI jobs, `elapsed_max_s` reports the slowest rank. Values:

| Policy | Final-write use case |
| --- | --- |
| `all` | Preserve existing behavior; analysis and checkpoints write at completion. |
| `restart_only` | Production runs that need a recovery point without extra final diagnostics. |
| `none` | Runs where no terminal output is wanted. |

A resumed terminal checkpoint is not overwritten merely because it is loaded.

## Promoted IO Examples

| Input | Purpose |
| --- | --- |
| `inputs/io/output_formats.athinput` | Modern PDF and spherical-slice readback, with existing spherical output retained. |
| `inputs/io/output_pdf_scalar_weight.athinput` | Four-dimensional PDF and scalar weighting. |
| `inputs/io/runtime_policy.athinput` | Timing and final-output policy. |
| `inputs/io/node_sharded_outputs.athinput` | MPI node sharding and restart manifest behavior. |

## Performance And Recovery Notes

- Keep checkpoint cadence appropriate to queue-walltime risk.
- Use `output_timing` on representative jobs to measure the slowest-rank cost
  of each stream rather than inferring output cost from total runtime.
- Node-sharded output reduces file counts but introduces node-level
  aggregation; qualify it on the target filesystem and node topology.
- Full-volume coarsened binary node output is supported. Avoid documenting or
  depending on sliced node-sharded `cbin` until its existing slice/readback
  extent issue is repaired separately.
