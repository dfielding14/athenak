## Output Files

| File type | Description | Typical files |
| --- | --- | --- |
| `vtk` | Legacy VTK mesh data for ParaView/VisIt | `vtk/basename.id.00000.vtk` |
| `tab` | ASCII sliced tables | `tab/basename.id.00000.tab` |
| `hst` | History diagnostics | `basename.hydro.hst`, `basename.mhd.hst` |
| `bin` | AthenaK binary mesh dumps | `bin/basename.id.00000.bin` |
| `cbin` | Uniform 3D active-zone full-volume coarsened binary output | `cbin_<id>_<factor>/basename.id.00000.cbin` |
| `pdf` | Legacy or modern V2 histograms | `pdf_<id>*/basename.00000.pdf` |
| `sphslice` | Fixed-radius binary angular slices | `bin/basename.id.r_<radius>.00000.sph.bin` |
| `rst` | Restart checkpoint or public node manifest | `rst/basename.00000.rst` |

Additional registered output types are `log`, `pvtk`, `trk`, `cart`, and
`sph`; see [Outputs](modules/outputs.md).

Supported `cbin` production is restricted to uniform 3D active-zone
full-volume output. Every emitted axis extent must be divisible by
`coarsen_factor`.

Use `single_file_per_rank = true` for rank shards or
`single_file_per_node = true` for node shards. Do not set both in one stream.
Node-sharded products are written beneath `node_########/` directories.
Modern PDF V2 and spherical-slice binary payload scalars currently use
host-native byte order.

Read analysis products with the public Python tools:

```bash
python vis/python/examples/read_io_outputs.py \
  bin run/bin/node_00000000/simulation.density.00000.bin \
  --assemble-shards

python vis/python/examples/read_io_outputs.py \
  pdf run/pdf_density/node_00000000/simulation.00000.pdf
```

### Restarting A Node-Sharded Job

With `file_type = rst` and `single_file_per_node = true`, use the public
manifest in `-r`:

```bash
mpirun -np 16 ./build-mpi/src/athena \
  -r run/rst/simulation.00001.rst \
  -d continued_run
```

Do not restart from `rst/node_########/*.payload.rst`. Native resume validates
the manifest and reads routed MeshBlock spans directly from node payloads. It
does not create a shared `.assembled` staging file. Public-manifest rename is
the commit point: postcommit cleanup failures preserve the resumable manifest
and declared payloads.

### Timing And Final Writes

Configure terminal writes and optional timing under `<time>`:

```ini
<time>
output_timing = true
final_output_policy = restart_only
```

`final_output_policy` accepts `all`, `restart_only`, or `none`. Under MPI,
`elapsed_max_s` reports the slowest rank for each timed output stream.

### Promoted IO Examples

The shipped `inputs/io/` directory includes `output_formats.athinput`,
`output_pdf_scalar_weight.athinput`, `runtime_policy.athinput`, and
`node_sharded_outputs.athinput`. See
[IO Outputs And Sharding](examples/io_outputs_and_sharding.md) for runnable
commands and compatibility limits. Before production use, qualify node
sharding and restart on a real multi-node allocation and the target parallel
filesystem.
