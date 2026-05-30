# Configuring AthenaK

AthenaK reads `.athinput` files containing labelled blocks (`<blockname>`)
and `name = value` entries. Later assignments override earlier ones, including
runtime overrides passed as `block/name=value`.

## Minimal Layout

```ini
<job>
basename = simulation_name

<mesh>
nx1 = 256
x1min = 0.0
x1max = 1.0
ix1_bc = periodic
ox1_bc = periodic

<time>
cfl_number = 0.4
tlim = 10.0
integrator = rk2
final_output_policy = all
output_timing = false

<hydro>
gamma = 1.4
reconstruct = plm
rsolver = hllc
```

Validate or override an input at launch:

```bash
./build/src/athena -i inputs/hydro/sod.athinput -n
./build/src/athena -i inputs/hydro/sod.athinput mesh/nx1=512
```

## Essential Blocks

| Block | Purpose | Common keys |
| --- | --- | --- |
| `<job>` | Naming and metadata | `basename` |
| `<mesh>`, `<meshblock>` | Domain, boundaries, and block decomposition | `nx1`, `x1min`, `ix1_bc`, block `nx1` |
| `<time>` | Integration and end-of-run behavior | `cfl_number`, `tlim`, `nlim`, `integrator`, `final_output_policy`, `output_timing` |
| `<output#>` | Individual output streams | `file_type`, cadence, `id`, `variable`, sharding keys |
| `<hydro>` / `<mhd>` | Fluid model | equation of state, reconstruction, solver |

## Final Output And Timing

The `<time>` block controls what is written when the run reaches its terminal
condition or a wall-clock stop:

```ini
<time>
final_output_policy = restart_only
output_timing = true
```

| Key | Values | Default | Meaning |
| --- | --- | --- | --- |
| `final_output_policy` | `all`, `restart_only`, `none` | `all` | Select output streams eligible for the final write. |
| `output_timing` | boolean | `false` | Emit one timing diagnostic for each output write. |

`restart_only` is useful for production jobs that need a final recovery point
without writing every analysis product again. Timing lines use:

```text
[output-io] event=final block=output_restart type=rst distribution=node elapsed_max_s=...
```

For MPI runs, the reported elapsed time is the maximum across ranks.

## Output Blocks

A stream defines a file product and a cadence:

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
same output block.

For `file_type = cbin`, set `coarsen_factor` to a power of two between `2`
and the shortest MeshBlock dimension, inclusive. The supported contract is
uniform three-dimensional active-zone full-volume output. The writer rejects
lower-dimensional, ghost-zone-expanded, static-refinement, AMR, and sliced
configurations before publication. Full-volume node-sharded `cbin` is
supported within that contract.

### Distribution Examples

Shared output:

```ini
<output_shared>
file_type = bin
id = shared
variable = hydro_w_d
dt = 0.1
```

Rank-sharded output:

```ini
<output_rank>
file_type = bin
id = rank
variable = hydro_w_d
dt = 0.1
single_file_per_rank = true
```

Node-sharded output:

```ini
<output_node>
file_type = bin
id = node
variable = hydro_w_d
dt = 0.1
single_file_per_node = true
```

Supported node-sharded products are `bin`, uniform 3D active-zone full-volume
`cbin`, modern `pdf`, `sphslice`, and `rst`. Node-sharded binary and uniform
3D active-zone full-volume coarsened-binary
writers emit additive inventory metadata and valid explicit empty shards.
Lower-dimensional, ghost-zone-expanded, static-refinement, AMR, and sliced
`cbin` are rejected explicitly and are not promoted workflows.

## Modern PDF Configuration

Modern PDFs support one to four scalar axes and volume, mass, or
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
variable_2 = mdot_sph
nbin2 = 128
bin2_min = -1.0
bin2_max = 1.0
scale2 = symlog
linthresh2 = 1.0e-5
weight = volume
single_file_per_node = true
```

| PDF key | Accepted values |
| --- | --- |
| `variable_1` ... `variable_4` | One to four contiguous scalar output-variable names |
| `nbinN` | Positive bin count |
| `binN_min`, `binN_max` | Range with minimum less than maximum |
| `scaleN` | `linear`, `log`, `symlog` |
| `linthreshN` | Positive threshold required for `symlog` |
| `weight` | `volume`, `mass`, `variable` |
| `weight_variable` | Required scalar name for `weight = variable` |
| `max_writer_allocation_bytes` | Positive per-writer allocation cap for histogram and edge arrays, host mirrors, copied and derived fields, metadata, and serialized staging; default `536870912` |

Legacy unsharded PDF keys remain accepted for compatibility, including
`mass_weighted`; use the modern form for new configurations.

PDF histogram updates use backend-portable atomics. MPI reductions stage
through host memory rather than requiring GPU-aware MPI.

The generic `mdot_*`, `edot_*`, and `vel_*` derived names are supported only
for single-fluid Hydro or MHD configurations. They are rejected for
`<ion-neutral>` two-fluid inputs because no unqualified quantity can identify
which fluid is meant.
These names are Newtonian-only. Total and thermal energy fluxes require an
ideal-gas total-energy module, `edot_sph_mag` requires MHD, and derived-array
outputs reject `ghost_zones = true`. Modern PDFs sample active zones only and
reject two-fluid `weight = mass`.

## Spherical Slice Configuration

```ini
<output_shell>
file_type = sphslice
id = shell_density
variable = hydro_w_d
slice_r = 0.5
ntheta = 64
nphi = 128
dt = 0.1
single_file_per_node = true
```

`slice_r` must be positive and strictly inside an origin-centered 3D domain;
`ntheta` and `nphi` must each be at least 2. The optional positive
`max_writer_allocation_bytes` cap defaults to `536870912`. The `variable`
must be a native scalar field such as `hydro_w_d` or
`mhd_w_d`; derived-array fields are rejected until spherical interpolation has
ghost-zone-safe derived sampling. `sphslice` is a binary fixed-radius data
product; the existing `file_type = sph` remains a separate output format.
Radius filename components use deterministic round-trip scientific tokens, so
`slice_r = 0.5` emits `r_5.0000000000000000e-01`.

## Restart Configuration

```ini
<output_restart>
file_type = rst
id = checkpoint
dcycle = 500
single_file_per_node = true
```

In node mode the visible `.rst` file is a manifest; restart using that
manifest, not its payload files. The manifest is the only public node restart
entry point. Native resume validates paths, symlink containment, inventory,
bounded positive segment records, payload sizes, replicated headers, and a
node-payload content marker before reading routed local MeshBlock spans
directly from node payloads without `.assembled` staging. Marked hard-link and
byte-copy aliases are rejected outside validated manifest loading.

## Physics Blocks And Tips

Presence of `<hydro>` or `<mhd>` activates the corresponding fluid module.
Problem generators expose `<problem>` and other module-specific blocks as
needed. General operational tips:

- content following `#` is ignored;
- use `./build/src/athena -i <input> -n` to validate before a long run;
- test output layouts and public Python readers with a small run before
  enabling node-sharded production output; and
- qualify node-sharded runs on the target multi-node system before relying on
  them for checkpoint recovery.
