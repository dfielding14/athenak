## IO Output Artifact Reference

| Product | Shared file form | Rank/node shard form | Public reader |
| --- | --- | --- | --- |
| Native mesh binary (`bin`) | `bin/<basename>.<id>.<number>.bin` | `bin/rank_########/...` or `bin/node_########/...` | `vis/python/bin_convert.py` |
| Uniform 3D active-zone full-volume coarsened binary (`cbin`) | `cbin_<id>_<factor>/<basename>.<id>.<number>.cbin` | `<directory>/rank_########/...` or `<directory>/node_########/...` | `vis/python/bin_convert.py` |
| Modern PDF V2 (`pdf`) | `pdf_<id>[_<axes>]/<basename>.<number>.pdf` | `<directory>/rank_########/...` or `<directory>/node_########/...` | `vis/python/read_pdf.py` |
| Spherical slice (`sphslice`) | `bin/<basename>.<id>.r_<radius>.<number>.sph.bin` | `bin/rank_########/...` or `bin/node_########/...` | `vis/python/read_sphslice.py` |
| Node restart (`rst`) | `rst/<basename>.<number>.rst` public manifest | `rst/node_########/<basename>.<number>.g<generation>.payload.rst` | `athena -r <manifest>` |

Rank and node directory components are zero-padded numeric identifiers.
Readers accept a shard path and discover the sibling family needed to
reconstruct one logical output.

### PDF Compatibility And V2 Payloads

Pure legacy unsharded PDF blocks retain their existing text artifacts. Modern
PDF output writes an ASCII header plus an `AKPDFV2` payload. Shared V2 products
store a dense flattened array. Rank/node V2 products store sparse indexed
contributions assembled by `read_pdf.py`. Each node-sharded V2 header records
the writer leader as `payload_rank`; the reader requires that value to match the
rank embedded in the binary preamble so misplaced payloads fail closed.
PDF V2 binary scalars currently use host-native byte order.

### Spherical Slice Payloads

`file_type = sphslice` writes `.sph.bin` files containing the requested
radius, angular dimensions, output metadata, variable names, and payload.
Shared files contain the full angular surface. Rank/node files carry angular
ownership records and are assembled by `read_sphslice.py`.
Spherical-slice binary scalars currently use host-native byte order. The
header-only reader validates intrinsic dimensions, positive radius, bounded
point count, layout, and the selected shard's declared identity; the full
reader discovers siblings and rejects incomplete or inconsistent families.
Neither path can prove the writer-side domain-interior condition without the
original mesh domain.

### Node Restart Manifests And Payloads

For `single_file_per_node = true`, the root `.rst` is the public manifest.
Always restart through that manifest. Do not use a `*.payload.rst` file
directly. Native resume validates the manifest and reads routed local
MeshBlock spans directly from the declared node payloads. It does not create a
shared `<manifest>.assembled` staging file. Manifest rename is the publication
commit point. Postcommit cleanup failures preserve the resumable manifest and
payloads.

### Qualification Boundary

Automated multi-rank tests exercise node-sharded output and restart on one
physical node. A real multi-node MPI run on the deployment filesystem remains
a production qualification step. Lower-dimensional, ghost-zone-expanded,
static-refinement, AMR, and sliced `cbin` producer workflows are deliberately
not advertised. For every supported `cbin` producer workflow, each emitted
axis extent must be divisible by `coarsen_factor`.
