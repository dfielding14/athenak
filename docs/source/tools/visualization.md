# Visualization Utilities (`vis/python`)

## Overview
The `vis/python/` toolbox ships with AthenaK to streamline post-processing of binary
outputs, generate HDF5/XDMF companions, and provide quick-look visualisations.  The
helpers depend only on NumPy and operate on the raw Athena++ binary format, so they can be
used on systems where the full Fortran post-processing stack is unavailable.

## Directory Highlights
| Script | Purpose | Notes |
|--------|---------|-------|
| `athena_read.py` | Legacy readers for tabular (`.tab`) and history files. | Maintained for backward compatibility with regression tooling (`vis/python/athena_read.py:1`). |
| `bin_convert_new.py` | High-performance binary reader/converter with AMR, coarsening, and multi-rank support. | Preferred interface for modern workflows (see below). |
| `bin_convert.py` | Original converter and HDF5 writer. | Still provides the `write_athdf` helper used by `bin_convert_new.py` when writing to disk (`vis/python/bin_convert.py:1484`). |
| `plot_mesh.py`, `plot_slice.py`, `plot_tab.py` | Quick plotting hooks targeting mesh topology, slices, and ASCII tables. | Rely on Matplotlib; useful for exploratory analysis. |

## Binary Conversion Toolkit (`bin_convert_new.py`)
`bin_convert_new.py` consolidates Athena++ binary parsing into a single module.  Its
private `_read_binary_impl` routine handles header parsing, AMR mesh-block metadata, and
field extraction, and is shared by all public entry points (`vis/python/bin_convert_new.py:1378`).

### Supported Data Sets
- Single-rank binary dumps via `read_binary` and `read_coarsened_binary`
  (`vis/python/bin_convert_new.py:74`).
- Multi-rank outputs written with `single_file_per_rank=false` by scanning every `rank_*`
  file through `read_all_ranks_binary` (`vis/python/bin_convert_new.py:87`).
- Coarsened equivalents for on-the-fly averaged dumps through
  `read_all_ranks_coarsened_binary` (`vis/python/bin_convert_new.py:168`).

All readers return dictionaries mirroring the structure of an `athdf` file: mesh-block
index arrays, geometry metadata, and a `mb_data` mapping from variable name to NumPy
arrays.

### athdf-Compatible Views
The module can synthesise an in-memory athdf representation without writing to disk:

| Function | Description |
|----------|-------------|
| `read_binary_as_athdf` | Prolongates/restricts mesh blocks into a uniform grid, supports sub-domain selection, and records optional refinement levels (`vis/python/bin_convert_new.py:242`). |
| `read_all_ranks_binary_as_athdf` | Multi-rank variant that stitches blocks across ranks before the athdf reshaping logic runs (`vis/python/bin_convert_new.py:386`). |
| `read_rank_binary_as_athdf` | Flattens a single-rank dump into a contiguous grid using the full domain extents stored in the header (`vis/python/bin_convert_new.py:1248`). |

#### `read_binary_as_athdf` Arguments
| Option | Default | Behavior |
|--------|---------|----------|
| `raw` | `False` | Return the raw mesh-block dictionary without athdf-style assembly. |
| `data` | `None` | Existing output buffer; pass a dict to reuse arrays instead of allocating new ones. |
| `quantities` | `None` | Iterable of variable names to extract; defaults to all variables in the dump. |
| `dtype` | `np.float32` | Data type for coordinate and field arrays. |
| `level` | Finest level | Target AMR level to project onto; coarser levels are prolongated, finer levels restricted. |
| `return_levels` | `False` | When `True`, populate a `Levels` array with the originating AMR level for every cell. |
| `subsample` | `False` | If enabled, skip restriction and simply copy blocks that match the requested level. |
| `fast_restrict` | `False` | Use a simplified (volume-blind) restriction path for coarse projections. |
| `x{1,2,3}_min`, `x{1,2,3}_max` | `None` | Spatial window for trimming the domain before interpolation. |
| `vol_func` | Cartesian volume | Callback returning a cell volume for restriction; override for curvilinear grids. |
| `vol_params` | `None` | Optional extra arguments forwarded to `vol_func` (reserved for future use). |
| `face_func_{1,2,3}` | `None` | Custom face-coordinate generators; default is a uniform linspace between bounds. |
| `center_func_{1,2,3}` | Midpoint | Custom cell-center evaluators given face positions. |
| `num_ghost` | `0` | Number of ghost zones to pad around the projected grid. |

All variants accept `quantities` filters, optional dtype overrides, and geometric helper
callbacks so downstream consumers can provide custom volume/centroid definitions.

### File Conversion and Visualisation Hooks
- `write_xdmf_for` emits ParaView/VisIt companion files while preserving mesh-block
  topology and separating cell-centred and face-centred fields
  (`vis/python/bin_convert_new.py:1148`).
- `convert_file` orchestrates the common “`.bin` → `.athdf` + `.xdmf`” pipeline.  It
  delegates the actual HDF5 write to the legacy `write_athdf` helper, so ensure
  `bin_convert.py` is importable before calling this wrapper
  (`vis/python/bin_convert_new.py:1222`).
- `athinput` reads Athena input decks into nested dictionaries, useful for embedding run
  metadata alongside converted outputs (`vis/python/bin_convert_new.py:1331`).

### Exported API
The `__all__` whitelist keeps the public surface area focused on conversion and metadata
helpers (`vis/python/bin_convert_new.py:1513`).  Downstream code should import only from
that list so future refactors can change internals without breaking user scripts.

### Example Workflow
```python
from pathlib import Path
from vis.python import bin_convert_new
from vis.python import bin_convert  # provides write_athdf

dump = Path("run/output/data.bin")

# Read into memory for analysis
athview = bin_convert_new.read_binary_as_athdf(dump, quantities=["rho", "vel1"])

# Write an HDF5/XDMF pair for ParaView
bin_convert_new.convert_file(str(dump))
```

`convert_file` writes `data.athdf` and `data.athdf.xdmf` beside the original dump.  When
working with multi-rank outputs, point the helper at the rank zero file and the module
will automatically ingest the sibling ranks before writing the consolidated product.

## Interop with Plotting Scripts
The conversion helpers integrate cleanly with the plotting utilities.  For example, the
`.athdf` generated by `convert_file` can be loaded with the routines in
`plot_slice.py`, while the neutral `athinput` reader simplifies annotating plots with the
exact solver settings extracted from the run deck.
