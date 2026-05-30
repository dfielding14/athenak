### Output Blocks

Each `<output#>` block defines one stream:

```ini
<output1>
file_type = bin
id = density
variable = hydro_w_d
dt = 0.1
single_file_per_node = true
```

| Field | Purpose |
| --- | --- |
| `file_type` | Product type, including `bin`, `cbin`, `pdf`, `sphslice`, `sph`, `rst`, `vtk`, `tab`, or `hst`. |
| `id` | Stable stream label used in paths and filenames. |
| `variable` | Field/group for products that require it. |
| `dt`, `dcycle` | Output cadence controls; supply one positive cadence. |
| `single_file_per_rank` | Write independent rank shards when supported. |
| `single_file_per_node` | Aggregate each supported stream into one file per shared-memory MPI node. |

`single_file_per_rank` and `single_file_per_node` cannot both be true in the
same output block. Supported node-sharded products are `bin`, uniform 3D
active-zone full-volume `cbin`, modern `pdf`, `sphslice`, and `rst`.

For `file_type = cbin`, set `coarsen_factor` to a power of two between `2`
and the shortest MeshBlock dimension, inclusive. Every emitted axis extent
must be divisible by that factor. Lower-dimensional, ghost-zone-expanded,
static-refinement, AMR, and sliced `cbin` configurations are rejected before
publication.

### Final Output And Timing

The `<time>` block controls terminal writes and optional timing diagnostics:

```ini
<time>
final_output_policy = restart_only
output_timing = true
```

| Key | Values | Default | Meaning |
| --- | --- | --- | --- |
| `final_output_policy` | `all`, `restart_only`, `none` | `all` | Select streams eligible for the final write. |
| `output_timing` | boolean | `false` | Emit one timing diagnostic for each output write. |

For MPI runs, reported `elapsed_max_s` values are the maximum across ranks.

### Modern PDF Configuration

Modern PDFs support one to four contiguous scalar axes and volume, mass, or
scalar-variable weighting:

```ini
<output_pdf>
file_type = pdf
id = flux_pdf
dt = 0.1
variable_1 = coord_r
nbin1 = 128
bin1_min = 1.0e-2
bin1_max = 10.0
scale1 = log
weight = volume
```

Use `scaleN = linear`, `log`, or `symlog`; a symlog axis also requires a
positive `linthreshN`. Legacy unsharded PDF keys, including `mass_weighted`,
remain accepted for compatibility. Use the modern form for new inputs.
`weight = mass` aborts on non-finite or non-positive conserved density.
`weight = variable` accepts finite, possibly signed scalar values. See
[Outputs](modules/outputs.md) for diagnostic semantics.
Modern PDF V2 binary payload scalars currently use host-native byte order.

### Spherical Slice Configuration

```ini
<output_shell>
file_type = sphslice
id = shell_density
variable = hydro_w_d
slice_r = 0.5
ntheta = 64
nphi = 128
dt = 0.1
```

`sphslice` writes a binary fixed-radius data product. It is separate from the
existing `file_type = sph` output. The origin-centered spherical surface must
fit inside a 3D domain: `slice_r` must be positive and strictly interior to
every domain face. It accepts
native state-backed scalar fields and native multi-field groups, and rejects
derived-array fields until angular interpolation has ghost-zone-safe derived
sampling. Writers and readers reject non-finite spherical-slice samples,
including finite in-memory values that overflow when narrowed to the
serialized float payload. Spherical-slice binary payload scalars currently use
host-native byte order.

### Restart Configuration

For node-sharded restart output, restart from the visible `.rst` manifest, not
from a node payload:

```ini
<output_restart>
file_type = rst
dcycle = 500
single_file_per_node = true
```

Native node resume validates the manifest and reads routed local MeshBlock
spans directly from node payloads. It does not create a shared `.assembled`
staging file. Public-manifest rename is the commit point: failures before it
roll back publication, while failures after it preserve the resumable manifest
and declared payloads. Qualify node-sharded runs on the target multi-node
system before relying on them for checkpoint recovery.
