# Visualization And IO Reading Utilities (`vis/python`)

## Supported Public Tools

AthenaK now exposes one canonical converter for native mesh binary output and
dedicated readers for new diagnostic formats:

| Script | Public responsibility |
| --- | --- |
| `bin_convert.py` | Read and convert `.bin` and `.cbin`, including shard assembly and HDF5/XDMF conversion. |
| `read_pdf.py` | Read legacy PDF text products and modern V2 shared/rank/node PDF data. |
| `read_sphslice.py` | Read shared/rank/node `.sph.bin` fixed-radius slices. |
| `examples/read_io_outputs.py` | Minimal command-line readback for promoted IO examples. |
| `athena_read.py` | Existing history/table utilities. |

`bin_convert.py` is the sole supported native binary conversion API. Do not
import or document a separate `bin_convert_new.py` module.

## Binary And Coarsened Binary

### Python API

```python
from vis.python import bin_convert

raw = bin_convert.read_binary("run/bin/simulation.prim.00000.bin")
raw_node = bin_convert.read_binary(
    "run/bin/node_00000000/simulation.prim.00000.bin",
    assemble_shards=True,
)
cbin_node = bin_convert.read_coarsened_binary(
    "run/cbin_coarse_2/node_00000000/simulation.coarse.00000.cbin",
    assemble_shards=True,
)
```

When `assemble_shards=True`, the reader discovers sibling `rank_########` or
`node_########` files and assembles their meshblock records into one logical
output. It validates file structure and consistent metadata across discovered
shards. New node-sharded binary and uniform 3D active-zone full-volume
coarsened-binary files add
`distribution`, `node`, `number of nodes`, and `number of meshblocks`
preheader fields. The canonical reader uses these optional fields to reject
missing or duplicate node IDs while accepting explicit empty shards. Shared
output needs no assembly flag. Metadata records, parameter dumps, and
accumulated MeshBlock payloads are bounded while each shard is read. Aggregate
payload and metadata totals are checked incrementally after each bounded shard
load and before combined aggregate materialization. The canonical readers
require uniform emitted MeshBlock extents within each file. Athdf-like
conversion helpers validate grid and logical-location metadata plus each
MeshBlock's exact logical physical interval before reconstruction, using zero relative
  tolerance and a storage-aware absolute tolerance capped at one eighth of the
  logical block width. They preflight cumulative coordinate,
  field, level-map, and restriction-map allocations plus NumPy
  coordinate/prolongation/restriction-generation temporaries derived from parsed
  metadata. Dense assembly helpers require a matching `num_ghost` argument for
  ghost-bearing products, and malformed singleton-axis ghost widths are
  rejected. Requested ghost zones are placed by interior MeshBlock width with
  coordinates extended outside root bounds; lower crop bounds retain cells they
  intersect; uncovered partial-shard level-map cells are initialized to `-1`.
  The preserved legacy single-MeshBlock ATHDF helper has no ghost argument and
  rejects ghost-bearing or sliced emitted extents rather than truncating them.

### Conversion

`convert_file` writes an Athena HDF5 file and XDMF companion beside the source
binary while preserving existing conversion helpers:

```python
from vis.python import bin_convert

bin_convert.convert_file("run/bin/simulation.prim.00000.bin")
```

The canonical module retains `write_athdf`, `write_xdmf_for`,
`convert_file`, and the read/assembly interfaces needed by existing workflows.

For node-sharded products, read/assemble the logical product before downstream
analysis; the public example script demonstrates the corresponding command
line:

```bash
python vis/python/examples/read_io_outputs.py \
  bin run/bin/node_00000000/simulation.prim.00000.bin \
  --assemble-shards
```

Uniform 3D active-zone full-volume node-sharded `.cbin` is supported.
Lower-dimensional, ghost-zone-expanded, static-refinement, AMR, and sliced
`cbin` are rejected explicitly before publication.
The canonical reader intentionally remains able to read validated historical
`.cbin` products outside the current writer matrix. That read-only
compatibility does not promote those layouts as new producer workflows.

## Reading PDFs

`read_pdf.py` automatically handles:

- legacy unsharded text PDFs written by unchanged old-style inputs;
- modern dense shared V2 PDF files; and
- modern sparse per-rank or per-node V2 shards.

Modern writers emit `AKPDFV2`. The reader also retains historical compatibility for
transitional unversioned dense and sparse binary payloads.

```python
from vis.python import read_pdf

legacy = read_pdf.read_pdf("run/pdf_legacy/simulation.00000.pdf")
modern = read_pdf.read_pdf("run/pdf_rho/simulation.00000.pdf")
node = read_pdf.read_pdf(
    "run/pdf_rho/node_00000000/simulation.00000.pdf",
)
```

Modern V2 validation covers the `AKPDFV2` preamble, dimensionality and edge
metadata, time/cycle agreement between shards, sparse index range, duplicate
records within an individual shard, truncation, and unexpected trailing
records. Sparse reconstruction also requires canonical sibling directory
names, matching path/header shard identifiers, and rank-sharded payload writer
ranks that match the rank directory. New V2 sparse shards declare a complete
rank or node inventory; historical transitional artifacts may omit those
additive fields. Node payloads still carry writer ranks rather than node IDs.
A declared V2 header requires a V2 payload preamble with the same cycle, sparse
V2 families require complete inventory declarations, and sibling headers must
agree on V2 declaration. Metadata must be finite where required, and a modern
dense file must use the shared distribution. Metadata and payload sizes are
  bounded before each shard load. Retained reference state is included while
  each replacement shard header is parsed. Explicit, generated, and legacy bin
  edges are preflighted before token-list or NumPy materialization. Legacy
  numeric rows are ASCII-only, parsed strictly, retained cumulatively, and bounded again
  before stacking. Payload-copy and duplicate-validation temporaries are included
  in the cumulative retained peak before materialization and aggregate
  accumulation, and incorporated local headers and sparse arrays are released
  before the next sibling read. Reconstructed aggregates remove
  shard-local identifiers. The reader sums normal contributions to one bin from
  distinct shards.

The public example script can be used for a quick summary:

```bash
python vis/python/examples/read_io_outputs.py \
  pdf run/pdf_rho/simulation.00000.pdf
```

## Reading Spherical Slices

Use `read_sphslice.py` for `file_type = sphslice` output:

```python
from vis.python import read_sphslice

shared = read_sphslice.read_sphslice(
    "run/bin/simulation.rho_shell.r_5.0000000000000000e-01.00000.sph.bin"
)
node = read_sphslice.read_sphslice(
    "run/bin/node_00000000/simulation.rho_shell.r_5.0000000000000000e-01.00000.sph.bin",
)
```

The assembled object contains the angular grid and sampled values. Sharded
readback rejects inconsistent headers, layout/distribution mismatches,
truncated payloads, incomplete shard inventory, and missing or duplicate
angular ownership. Shared files declare `layout=dense`; sharded files declare
`layout=sparse_angles`. Empty shards are valid inventory members. Whole-file
  reads, cumulative header bytes, individual metadata lines, variable-token
  expansion before splitting, retained reference variable-metadata summary through final
  coordinate generation, coordinate allocations, embedded input dumps, payload copies,
  duplicate-validation and ownership-diagnostic temporaries, and embedded-header
  offsets are bounded before loading or materialization. Incorporated sparse
  arrays are released before the next sibling read, and reconstructed aggregates
  remove shard-local identifiers.
The filename radius component uses a deterministic round-trip scientific
token; for example, `slice_r = 0.5` emits
`r_5.0000000000000000e-01`.

```bash
python vis/python/examples/read_io_outputs.py \
  sphslice run/bin/node_00000000/simulation.rho_shell.r_5.0000000000000000e-01.00000.sph.bin
```

`sphslice` is a new binary analysis product; it does not replace the existing
`file_type = sph` VTK-oriented output.

## Matching Readers To Distribution

| Output type | Shared | Per-rank | Per-node |
| --- | --- | --- | --- |
| `.bin` | `bin_convert.read_binary` | `read_binary(..., assemble_shards=True)` | `read_binary(..., assemble_shards=True)` |
| `.cbin` full volume | `bin_convert.read_coarsened_binary` | `read_coarsened_binary(..., assemble_shards=True)` | `read_coarsened_binary(..., assemble_shards=True)` |
| PDF | `read_pdf.read_pdf` | `read_pdf` on one shard, automatic assembly | `read_pdf` on one shard, automatic assembly |
| `.sph.bin` | `read_sphslice.read_sphslice` | `read_sphslice` on one shard, automatic assembly | `read_sphslice` on one shard, automatic assembly |
| `.rst` | `athena -r <file>` | `athena -r <file>` | `athena -r <public-manifest.rst>` |

## Reproducible Example Inputs

The source tree includes:

```text
inputs/io/output_formats.athinput
inputs/io/output_pdf_scalar_weight.athinput
inputs/io/runtime_policy.athinput
inputs/io/node_sharded_outputs.athinput
```

These examples are exercised by the IO test suite, so the filenames, parser
keys, and public readback tools described here are kept aligned with runnable
inputs.
