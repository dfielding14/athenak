## Shipped Utilities

| Script | Public behavior |
| --- | --- |
| `athena_read.py` | Reads `.tab`, `.hst`, and `.athdf` data into Python objects. |
| `plot_tab.py`, `plot_hst.py` | Plot selected table and history variables. |
| `bin_convert.py` | Reads and converts `.bin` and `.cbin`, including shard assembly and HDF5/XDMF conversion. |
| `read_pdf.py` | Reads legacy PDF text products and modern V2 shared/rank/node PDF data; provides `read_pdf_header()` for validated header-only inspection. |
| `read_sphslice.py` | Reads shared/rank/node `.sph.bin` fixed-radius slices and provides `read_sphslice_header()`. New explicit shards require complete inventories; historical version-1 shared/rank layouts remain readable. |
| `examples/read_io_outputs.py` | Provides minimal command-line readback for promoted IO examples. |
| `make_athdf.py` | Batch `.bin`/`.cbin` wrapper around `bin_convert.py` for a filename stem, with shard assembly and reader-budget flags; excludes the separate `.sph.bin` family. |
| `plot_mesh.py`, `plot_slice.py`, `cartgrid.py` | Additional mesh, slice, and Cartesian-grid helpers. |

`bin_convert.py` is the sole supported native binary conversion API.
The public readers accept `io_reader_common.ReaderLimits` overrides for live
memory, header-read, and payload-read budgets.
