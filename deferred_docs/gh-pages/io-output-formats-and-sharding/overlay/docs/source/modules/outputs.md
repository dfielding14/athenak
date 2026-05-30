# Module: Outputs

## Overview

AthenaK output blocks write solution state, reduced diagnostics, checkpoints,
and analysis-oriented products at user-defined cadences. The IO extensions
described on this page add:

- shared, per-rank, and per-node output distributions for supported binary
  products and restarts;
- one- through four-dimensional PDFs with versioned modern files;
- fixed-radius spherical slices (`file_type = sphslice`);
- output timing diagnostics; and
- an explicit end-of-run output policy.

Existing formats remain available. In particular, `file_type = sph` is the
existing spherical VTK-style output and is distinct from the new binary
`sphslice` product.

## Output Blocks

An output stream is configured in an `<output#>` block:

```ini
<output1>
file_type = bin
id = density
variable = hydro_w_d
dt = 0.1
```

Cadence is controlled by `dt` or `dcycle`. Common diagnostic products include:

| `file_type` | Product | Reader or consumer |
| --- | --- | --- |
| `bin` | Native mesh binary dump | `vis/python/bin_convert.py` |
| `cbin` | Coarsened mesh binary dump | `vis/python/bin_convert.py` |
| `pdf` | Histogram/PDF output | `vis/python/read_pdf.py` |
| `sphslice` | Binary fixed-radius angular slice | `vis/python/read_sphslice.py` |
| `sph` | Existing spherical output | Existing VTK workflow |
| `rst` | Exact restart checkpoint | `athena -r` |
| `vtk`, `tab`, `hst`, `log` | Existing visualization or diagnostic streams | Existing tools |

## File Distribution

Supported binary output streams choose one file distribution. Omitting both
sharding flags preserves shared-file behavior.

| Setting | Meaning | Output layout |
| --- | --- | --- |
| Neither flag set | Shared output file | One product per output event |
| `single_file_per_rank = true` | Independent file per MPI rank | `rank_########/` shards |
| `single_file_per_node = true` | One aggregate file per shared-memory MPI node | `node_########/` shards |

The two flags are mutually exclusive. An input that enables both is rejected.
The node mode constructs a shared-memory node communicator only when a
node-sharded stream or node restart manifest requires it.

### Supported Per-Node Products

| Product | Per-node support | Notes |
| --- | --- | --- |
| `bin` | Yes | Full output and Cartesian slicing are reader-tested. |
| `cbin` | Yes for full-volume output | Do not promote sliced `cbin`; existing shared sliced readback has a separate extent defect. |
| Modern `pdf` | Yes | Shards use the sparse-coordinate V2 layout. |
| `sphslice` | Yes | Angular ownership is assembled by the reader. |
| `rst` | Yes | Public manifest references transactional node payloads. |

The feature tests exercise MPI ranks sharing a physical node. Production
adoption of `single_file_per_node` should also include a real multi-node run
that checks empty/non-owning node cases and restart recovery.

## N-D PDF Output

### Modern Configuration

Modern PDFs use one to four named axes:

```ini
<output_pdf>
file_type = pdf
id = rho_r
dt = 0.1

variable_1 = coord_r
bin1_min = 1.0e-2
bin1_max = 1.0e1
nbin1 = 128
scale1 = log

variable_2 = hydro_w_d
bin2_min = 1.0e-5
bin2_max = 1.0e1
nbin2 = 96
scale2 = log

weight = mass
```

| Parameter | Accepted value or default | Description |
| --- | --- | --- |
| `variable_1` ... `variable_4` | Scalar output variable | Axes are contiguous from axis 1; one axis is required. |
| `nbinN` | Positive integer | Interior bin count for axis `N`. |
| `binN_min`, `binN_max` | Real values with min less than max | Interior bin range for axis `N`. |
| `scaleN` | `linear` (default), `log`, `symlog` | Axis-edge transform. |
| `linthreshN` | Positive real | Required for `scaleN = symlog`. |
| `weight` | `volume` (default), `mass`, `variable` | Contribution accumulated in each bin. |
| `weight_variable` | Scalar output variable | Required only for `weight = variable`. |

`scaleN = log` requires positive bin limits. Modern syntax is selected by
`variable_1` and by other modern options such as sharding, `scaleN`, or
`weight`; use these spellings in new inputs.

PDF contributions use active zones. Weighting incorporates cell volume, so
`mass` integrates density over cell volume and `variable` integrates the
selected scalar over cell volume.

### Variables For Diagnostics

The PDF axes and variable weights accept normal scalar output variables and
these derived scalar families:

| Family | Available names |
| --- | --- |
| Cartesian position | `coord_x`, `coord_y`, `coord_z` |
| Spherical position | `coord_r`, `coord_theta`, `coord_phi`, `coord_costheta`, `coord_abscostheta` |
| Cylindrical position | `coord_cyl_R`, `coord_cyl_phi`, `coord_cyl_z` |
| Velocity projections | `vel_sph_r`, `vel_sph_theta`, `vel_sph_phi`, `vel_cyl_R`, `vel_cyl_phi` |
| Spherical mass flux | `mdot_sph`, `mdot_sph_out`, `mdot_sph_in` |
| Vertical mass flux | `mdot_vert`, `mdot_vert_out`, `mdot_vert_in` |
| Spherical energy flux | `edot_sph`, `edot_sph_out`, `edot_sph_in`, `edot_sph_kin`, `edot_sph_th`, `edot_sph_mag` |
| Vertical energy flux | `edot_vert`, `edot_vert_out`, `edot_vert_in` |
| Passive scalars | `hydro_u_s_N`, `hydro_w_s_N`, `mhd_u_s_N`, `mhd_w_s_N` |

The suffix `N` in passive scalar variables is a zero-based scalar index, for
example `hydro_w_s_0`.

Unqualified `mdot_*`, `edot_*`, and `vel_*` diagnostics are defined for
single-fluid Hydro or MHD output only. Inputs enabling `<ion-neutral>` contain
both fluids and reject these generic names until module-qualified two-fluid
diagnostics are defined.

### Legacy Compatibility

Existing unsharded PDF inputs using the older keys remain supported:

```ini
<output_legacy_pdf>
file_type = pdf
variable = hydro_w_d
bin_min = 0.0
bin_max = 2.0
nbin = 32
mass_weighted = true
```

Pure legacy one- and two-axis unsharded configurations continue to emit the
legacy text files byte-for-byte relative to the frozen `origin/main`
fixtures. `mass_weighted` remains valid in that legacy path. New inputs
should use the modern axis and `weight` keys.

### Modern V2 PDF Files

Modern output writes an ASCII description file and a versioned binary payload:

```text
pdf_<id>[_<axis-labels>]/
  <basename>.header.pdf
  <basename>.<output-number>.pdf
```

For rank/node sharding, each shard writes its V2 header beside its payload
below the corresponding shard directory:

```text
pdf_<id>[_<axis-labels>]/rank_00000000/<basename>.header.pdf
pdf_<id>[_<axis-labels>]/rank_00000000/<basename>.<output-number>.pdf
pdf_<id>[_<axis-labels>]/node_00000000/<basename>.header.pdf
pdf_<id>[_<axis-labels>]/node_00000000/<basename>.<output-number>.pdf
```

The header identifies `AthenaK PDF format version=2`, distribution, axes,
edges, layout, and weight. Binary files begin with the `AKPDFV2` magic and
carry version, layout, dimensionality, rank/node identifier, count, time, and
cycle before their values.

| Distribution | V2 layout | Contents |
| --- | --- | --- |
| Shared | `dense` | Complete flattened histogram of doubles. |
| Per-rank or per-node | `sparse_coo` | `(flattened-bin-index, value)` contributions assembled by the reader. |

Use `vis/python/read_pdf.py` for both legacy PDFs and modern shared or sharded
V2 files. It validates header/preamble consistency, dimensions, edge arrays,
truncated records, duplicate sparse indices within a shard, and shard
metadata. Contributions to the same bin from different shards are summed.

## Spherical Slice Output

`sphslice` samples a scalar field on a fixed spherical radius and writes an
angular array:

```ini
<output_slice>
file_type = sphslice
id = rho_shell
variable = hydro_w_d
slice_r = 0.5
ntheta = 64
nphi = 128
dt = 0.1
single_file_per_node = true
```

| Parameter | Default | Meaning |
| --- | --- | --- |
| `variable` | Required | Native scalar field backed by simulation state, for example `hydro_w_d` or `mhd_w_d`. |
| `slice_r` | Required | Radius inside the domain. |
| `ntheta` | `64` | Number of polar angular cells; must be at least 2. |
| `nphi` | `128` | Number of azimuthal angular cells; must be at least 2. |

Files use a `.sph.bin` suffix and include the requested radius:

```text
bin/<basename>.<id>.r_<radius>.<output-number>.sph.bin
bin/rank_00000000/<basename>.<id>.r_<radius>.<output-number>.sph.bin
bin/node_00000000/<basename>.<id>.r_<radius>.<output-number>.sph.bin
```

Derived-array fields, including `coord_*`, `mdot_*`, `edot_*`, and `vel_*`,
are rejected for `sphslice` because its trilinear sampling may require
ghost-zone values that those derived arrays do not yet populate.

Use `vis/python/read_sphslice.py` to read and automatically assemble sharded
angular data when a shard path is supplied. The
reader rejects truncated files, inconsistent metadata, and missing or
duplicate angular ownership in sharded products.

## Node-Sharded Restart Files

For a restart stream configured with `single_file_per_node = true`, the public
restart filename is a manifest and each node writes a transactionally
published payload:

```ini
<output_restart>
file_type = rst
id = checkpoint
dcycle = 1000
single_file_per_node = true
```

```text
rst/<basename>.<output-number>.rst
rst/node_00000000/<basename>.<output-number>.g<generation>.payload.rst
rst/node_00000001/<basename>.<output-number>.g<generation>.payload.rst
...
```

Restart from the public manifest, never from an individual payload:

```bash
./build-mpi/src/athena -r run/rst/simulation.00001.rst
```

The runtime validates relative payload paths, ordered contiguous node IDs,
consistent generations, file byte counts, segment coverage, and the manifest
completion marker before reconstructing restart input. Malformed or incomplete
manifests are rejected.

## End-Of-Run Output Policy

Configure final writes under `<time>`:

```ini
<time>
final_output_policy = restart_only
output_timing = true
```

| `final_output_policy` value | Behavior at a normal terminal condition or wall-clock stop |
| --- | --- |
| `all` | Write every output stream eligible for a final write. This is the default. |
| `restart_only` | Write final restart streams only. |
| `none` | Write no additional final output. |

Final writes increment output numbering normally. A resumed run already at a
terminal state does not overwrite its existing final restart checkpoint.

When `output_timing = true`, rank zero emits one timing line per written
output:

```text
[output-io] event=scheduled block=output_pdf type=pdf distribution=node elapsed_max_s=...
```

`event` is `initial`, `scheduled`, or `final`. Under MPI,
`elapsed_max_s` is the maximum write time over ranks and is therefore the
useful quantity for identifying synchronization-limited IO.

## Executable Examples And Readers

The code branch supplies input decks in `inputs/io/`:

| Input deck | Demonstrates |
| --- | --- |
| `output_formats.athinput` | Modern PDF, `sphslice`, and existing `sph` coexistence. |
| `output_pdf_scalar_weight.athinput` | Four-axis PDF and scalar-variable weighting. |
| `runtime_policy.athinput` | Timing and final-output behavior. |
| `node_sharded_outputs.athinput` | Node binary/coarsened/PDF/slice/restart products. |

The associated readback entry point is:

```bash
python vis/python/examples/read_io_outputs.py pdf run/pdf_nd3_coord_abscostheta_vel_sph_r/io_formats.00000.pdf
python vis/python/examples/read_io_outputs.py sphslice run/bin/node_00000000/io_node_example.density.r_0.25.00000.sph.bin
python vis/python/examples/read_io_outputs.py bin run/bin/node_00000000/io_node_example.density.00000.bin --assemble-shards
```

See the Visualization Utilities page and the worked IO example for complete
commands and compatibility limits.
