# Output File Sharding

AthenaK supports three output layouts for the formats that can write large mesh
or diagnostic payloads:

| Mode | Input parameter | Files per dump | Coordination scope | Typical use |
| --- | --- | ---: | --- | --- |
| Shared | default | 1 | all MPI ranks | small and moderate runs, simplest file layout |
| Per rank | `single_file_per_rank=true` | `nranks` | each rank independently | debugging or runs where file count is still manageable |
| Per node | `single_file_per_node=true` | about `nnodes` | ranks on the same physical node | large MPI runs where shared files contend and per-rank output creates too many files |

`single_file_per_node` groups MPI ranks by physical node using a node-local
communicator. Simulation ownership does not change: MeshBlocks and physics data
remain distributed by MPI rank. Only the disk layout and the communication used
to assemble each file shard change.

## Configuration

Use one of the sharding flags in each supported output block:

```ini
<output1>
file_type = bin
dt = 0.1
variable = hydro_w
single_file_per_node = true
```

Only one sharding flag may be true in an output block. Setting both
`single_file_per_rank=true` and `single_file_per_node=true` is an error.

Supported sharded output types:

| `file_type` | Shared | Per rank | Per node | Notes |
| --- | --- | --- | --- | --- |
| `rst` | yes | yes | yes | per-node mode writes a shared manifest plus node payload shards |
| `bin` | yes | yes | yes | mesh binary output |
| `cbin` | yes | yes | yes | coarsened mesh binary output |
| `sphslice` | yes | yes | yes | partitioned files store sparse angle ownership |
| `pdf` | yes | yes | yes | partitioned files store sparse COO histograms |

Other output formats continue to use their existing layout.

## Directory Layout

Shared output writes the historical single file at the format-specific root.
Per-rank and per-node output write under deterministic shard directories:

```text
bin/rank_00000000/<basename>.00000.bin
bin/node_00000000/<basename>.00000.bin
rst/rank_00000000/<basename>.00000.rst
rst/node_00000000/<basename>.00000.rst
```

The directory number is zero-padded and based on the MPI rank for per-rank mode
or the dense node id for per-node mode.

## Restart Files

Restart output has one additional rule in per-node mode. The global metadata is
written once in a shared manifest, and the heavy field payload is written in one
file per node:

```text
rst/<basename>.00000.rst
rst/node_00000000/<basename>.00000.rst
rst/node_00000001/<basename>.00000.rst
```

The manifest records enough layout metadata for the reader to locate payload
shards, preserve MeshBlock ordering, and validate the payload stride. The reader
accepts either the manifest path or a path inside a `node_XXXXXXXX/` or
`rank_XXXXXXXX/` directory. Examples:

```bash
./athena -r rst/run.00000.rst
./athena -r rst/node_00000000/run.00000.rst
./athena -r rst/rank_00000000/run.00000.rst
```

Shared and per-rank restart files retain the historical on-disk layout and
remain readable through their established paths. Only per-node restart output
uses the manifest extension.

For per-node restarts, each node opens its payload shard collectively over the
node communicator. Large reads are coalesced where practical and are issued
through chunked byte MPI-IO helpers so the transfer does not depend on MPI
datatype counts fitting in `int`.

## Spherical-Slice Files

Shared spherical-slice output writes the full dense surface. Partitioned
rank/node output writes only the angular points owned by each shard:

1. ASCII preheader and copied input block.
2. `int32` angle indices, where `index = itheta*nphi + iphi`.
3. `float32` values in variable-major order.

The Python reader `vis/python/read_sphslice.py` can be pointed at any shard. It
discovers matching sibling `rank_*` or `node_*` files, validates the shard
metadata, checks that angle coverage is complete and non-overlapping, and
returns a dense array with shape `(ntheta, nphi, nvars)`.

## PDF Files

Shared PDF output writes a dense histogram:

```text
double time
double values[total_bins]
```

Partitioned rank/node PDF output writes sparse COO records:

```text
double time
uint32 nnz
uint32 idx[nnz]
double val[nnz]
```

The header file records `format = dense` or `format = sparse_coo`,
`distribution = shared|rank|node`, dimensional metadata, and `total_bins`. The
reader `vis/python/read_pdf.py` can be pointed at any shard and sums all matching
rank/node sparse records into a dense histogram.

PDF inputs accept the historical one- or two-variable form, including
`mass_weighted=true`, and the N-dimensional form:

```ini
<output2>
file_type = pdf
variable_1 = hydro_w_d
bin1_min = 1.0e-3
bin1_max = 1.0e2
nbin1 = 128
scale1 = log
variable_2 = hydro_w_vx
bin2_min = -1.0
bin2_max = 1.0
nbin2 = 128
scale2 = symlog
linthresh2 = 1.0e-4
weight = variable
weight_variable = hydro_w_d
single_file_per_node = true
```

The N-dimensional form supports up to four `variable_N` axes with
`scaleN=linear|log|symlog`. For `symlog`, specify a positive `linthreshN`.
The `weight` option accepts `volume`, `mass`, or `variable`; when
`weight=variable`, also set `weight_variable`. Do not set `mass_weighted` and
`weight` inconsistently.

## Choosing A Mode

Use shared output when the run is small enough that a single file is not a
contention point and you want the simplest post-processing path.

Use per-rank output when debugging rank-local output or when the number of ranks
is small enough that the filesystem can handle one file per rank per dump.

Use per-node output for large MPI jobs where shared files are slow or unstable,
but per-rank output creates an impractical number of files. This is the intended
default scaling path for large restart dumps.

## Robustness Notes

- Per-node collectives include ranks with zero local payload, so a sparse shard
  cannot leave part of the communicator outside a collective MPI-IO call.
- Large MPI-IO reads and writes are chunked as bytes to avoid MPI count
  overflow.
- Restart path normalization is path-component based, so relative paths,
  absolute paths, and paths inside shard directories are handled consistently.
- A root restart path is treated as a per-node manifest only when its stored
  restart configuration selects per-node output and a matching node payload is present.
- Partitioned Python readers validate shard metadata before reassembly.

## Tests

The branch includes focused regression coverage:

| Test | Coverage |
| --- | --- |
| `tst/scripts/restart/per_node_restart.py` | per-node manifest/payload creation, restart from manifest and shard paths, shared/per-rank compatibility, stale-node collision handling, forced small MPI-IO chunks |
| `tst/scripts/hydro/binary_partitioned.py` | shared, per-rank, and per-node binary and coarsened-binary equivalence, including an empty first rank shard, through `vis/python/bin_convert.py` |
| `tst/scripts/hydro/sphslice_partitioned.py` | shared, per-rank, and per-node spherical-slice equivalence plus sparse output stats |
| `tst/scripts/hydro/pdf_partitioned.py` | shared, per-rank, and per-node PDF equivalence, legacy mass-weighted compatibility, three-dimensional `symlog` input, and variable weighting |

These tests require an MPI-enabled build and `mpiexec`.
