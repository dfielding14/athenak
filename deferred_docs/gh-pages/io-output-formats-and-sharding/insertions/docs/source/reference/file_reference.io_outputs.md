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
the `AKPDFV2` magic and include version, layout, dimensionality, shard
identifier, count, time, and cycle. Shared V2 products store a dense flattened
array. Rank/node V2 products store sparse indexed contributions that are
summed by `read_pdf.py`.

### Spherical Slice Payloads

`file_type = sphslice` writes `.sph.bin` files containing the requested
radius, angular dimensions, output metadata, variable names, and payload.
Shared files contain the full angular surface. Rank/node files carry angular
ownership records that `read_sphslice.py` reassembles while checking for
duplicates or gaps.

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
- missing, duplicate, or out-of-order node identifiers;
- inconsistent generations or incomplete manifests;
- mismatching payload byte counts or coverage.

Always restart through the manifest:

```bash
mpirun -np <ranks> ./build-mpi/src/athena -r run/rst/<basename>.<number>.rst
```

Do not use a `*.payload.rst` file directly.

### Qualification Boundary

Node-sharded binary, full-volume coarsened binary, modern PDF, spherical
slice, and restart behavior are covered by automated multi-rank tests on one
physical node. A multi-node MPI run on the deployment filesystem remains a
production qualification step. Sliced node-sharded `cbin` is deliberately not
advertised pending repair of a pre-existing sliced coarsened-output readback
issue.
