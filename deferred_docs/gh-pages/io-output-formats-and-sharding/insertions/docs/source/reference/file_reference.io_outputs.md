## IO Output Artifact Reference

Merge this section into `docs/source/reference/file_reference.md` in the
output-file reference area after the IO feature branch has merged.

### Analysis Output Families

| Product | Shared file form | Rank/node shard form | Public reader |
| --- | --- | --- | --- |
| Native mesh binary (`bin`) | `bin/<basename>.<id>.<number>.bin` | `bin/rank_########/...` or `bin/node_########/...` | `vis/python/bin_convert.py` |
| Full-volume coarsened binary (`cbin`) | `cbin_<id>_<factor>/<basename>.<id>.<number>.cbin` | `<directory>/rank_########/...` or `<directory>/node_########/...` | `vis/python/bin_convert.py` |
| Modern PDF V2 (`pdf`) | `pdf_<id>[_<axes>]/<basename>.<number>.pdf` | `<directory>/rank_########/...` or `<directory>/node_########/...` | `vis/python/read_pdf.py` |
| Spherical slice (`sphslice`) | `bin/<basename>.<id>.r_<radius>.<number>.sph.bin` | `bin/rank_########/...` or `bin/node_########/...` | `vis/python/read_sphslice.py` |
| Node restart (`rst`) | `rst/<basename>.<number>.rst` public manifest | `rst/node_########/<basename>.<number>.g<generation>.payload.rst` | `athena -r <manifest>` |

Rank and node directory components are zero-padded numeric identifiers.
Readers accept a shard path and find the sibling family needed to reconstruct
one logical output. For `.bin` and `.cbin`, request assembly through
`bin_convert.py`; PDF and spherical-slice readers discover sharded families
when given a shard path.

Node-sharded `.bin` and full-volume `.cbin` files add `distribution`, `node`,
`number of nodes`, and `number of meshblocks` preheader fields. These are
additive: legacy shared and rank files retain their existing schema. Node
writers publish explicit empty shards, and the canonical reader requires the
complete dense node inventory when the additive metadata is present. Each
nonempty `.bin` or `.cbin` file must also use one uniform emitted MeshBlock
extent; malformed mixed-extent files are rejected before dense reconstruction.

### PDF Compatibility And V2 Payloads

Pure legacy unsharded PDF blocks retain the existing text artifact pair:

```text
pdf_<id>[_<axes>]/<basename>.bins.pdf
pdf_<id>[_<axes>]/<basename>.<number>.pdf
```

Modern shared PDF output writes:

```text
pdf_<id>[_<axes>]/<basename>.header.pdf
pdf_<id>[_<axes>]/<basename>.<number>.pdf
```

Modern rank/node output writes each header beside its payload shard:

```text
pdf_<id>[_<axes>]/rank_########/<basename>.header.pdf
pdf_<id>[_<axes>]/rank_########/<basename>.<number>.pdf
pdf_<id>[_<axes>]/node_########/<basename>.header.pdf
pdf_<id>[_<axes>]/node_########/<basename>.<number>.pdf
```

The ASCII header identifies AthenaK PDF format version 2, dimensions, edges,
weighting, distribution, and dense or sparse layout. Payload files begin with
the `AKPDFV2` magic and include version, layout, dimensionality, writer rank,
count, time, and cycle. Shared V2 products store a dense flattened array.
Rank/node V2 products store sparse indexed contributions that are summed by
`read_pdf.py`. Sparse headers record `rank`/`number_of_ranks` or
`node`/`number_of_nodes`; path and header shard IDs must match, and
rank-sharded payload writer ranks must also match the rank directory. Node
payloads still carry the writer rank rather than the node ID. The reader
requires canonical sibling directories. New V2 sparse shards declare a
complete dense inventory; historical transitional unversioned dense and sparse
binary payloads remain readable and may omit those additive inventory fields.
Sibling headers must agree on V2 declaration. When a header declares V2, the
reader requires a V2 payload preamble with the same cycle, bounds metadata and
payload reads before construction, and removes shard-local identifiers from
the reconstructed aggregate. Each modern header or payload file publishes
independently through `<file>.tmp` followed by atomic rename.

### Spherical Slice Payloads

`file_type = sphslice` writes `.sph.bin` files containing the requested
radius, angular dimensions, output metadata, variable names, and payload.
Shared files declare `layout=dense` and contain the full angular surface.
Rank/node files declare `layout=sparse_angles`, carry angular ownership and
shard-inventory records, and preserve valid explicit empty shards.
`read_sphslice.py` reassembles them while checking for layout mismatches,
inventory gaps, and duplicate or missing angular ownership. Writers publish
through `<file>.tmp` followed by rename. Reader whole-file, coordinate-array,
and embedded-header-offset bounds apply before loading, and reconstructed
aggregates remove shard-local identifiers.

The writer accepts native state-backed scalar variables for spherical slices.
It rejects derived-array variables because trilinear angular sampling can
require derived ghost-zone values that are not populated by the current
derived-variable kernels.

Do not conflate `.sph.bin` with the existing `file_type = sph` output; both
are retained and have separate workflows.

### Node Restart Manifests And Payloads

For `single_file_per_node = true`, the root `.rst` is a public manifest. It
records the ordered node payloads written for one checkpoint generation.
Payloads are published before the manifest is committed, and the runtime
rejects:

- absolute, traversing, or malformed payload paths;
- payload paths that escape the manifest directory through symlinks;
- missing, duplicate, or out-of-order node identifiers;
- inconsistent generations or incomplete manifests;
- mismatching payload byte counts or coverage; and
- replicated payload headers that differ from canonical payload zero;
- oversized declared payload or segment inventories and non-positive segments; and
- absent or corrupt node-payload content markers.

Always restart through the manifest:

```bash
mpirun -np <ranks> ./build-mpi/src/athena -r run/rst/<basename>.<number>.rst
```

Do not use a `*.payload.rst` file directly. Generated payloads include a content
marker after the replicated parameter dump, so hard links and byte copies are
also rejected as ordinary restart entry points.

After validation, AthenaK routes each rank's local MeshBlock spans directly
from the declared node payloads using chunked positioned reads. Native resume
does not create a shared `<manifest>.assembled` staging file.

### Qualification Boundary

Node-sharded binary, full-volume coarsened binary, modern PDF, spherical
slice, and restart behavior are covered by automated multi-rank tests on one
physical node. A multi-node MPI run on the deployment filesystem remains a
production qualification step. Sliced `cbin` is deliberately not advertised:
construction rejects emitted extents that are incompatible with supported
coarsening.
