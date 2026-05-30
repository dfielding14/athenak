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
shards. Tests cover emitted empty/non-owning binary-slice cases. Shared output
needs no assembly flag.

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

Full-volume node-sharded `.cbin` is supported. Sliced `cbin` is not a promoted
workflow in this feature because the pre-existing shared sliced-output path
still requires a separate extent/readback repair.

## Reading PDFs

`read_pdf.py` automatically handles:

- legacy unsharded text PDFs written by unchanged old-style inputs;
- modern dense shared V2 PDF files; and
- modern sparse per-rank or per-node V2 shards.

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
records. The reader sums normal contributions to one bin from distinct
shards.

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
    "run/bin/simulation.rho_shell.r_0.5.00000.sph.bin"
)
node = read_sphslice.read_sphslice(
    "run/bin/node_00000000/simulation.rho_shell.r_0.5.00000.sph.bin",
)
```

The assembled object contains the angular grid and sampled values. Sharded
readback rejects inconsistent headers, truncated payloads, and missing or
duplicate angular ownership.

```bash
python vis/python/examples/read_io_outputs.py \
  sphslice run/bin/node_00000000/simulation.rho_shell.r_0.5.00000.sph.bin
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
