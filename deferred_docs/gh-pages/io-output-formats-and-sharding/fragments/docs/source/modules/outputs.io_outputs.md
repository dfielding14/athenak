## Registered Formats

| `file_type` | Public purpose |
| --- | --- |
| `tab` | Formatted one-dimensional table slice |
| `hst` | Integrated time-series diagnostics |
| `log` | Event-counter log |
| `vtk` | Mesh fields in legacy VTK |
| `pvtk` | Particle VTK output |
| `trk` | Tracked-particle stream |
| `cbin` | Uniform 3D active-zone full-volume coarsened binary mesh data |
| `pdf` | Legacy or modern V2 histograms |
| `bin` | Native binary mesh data |
| `cart` | Cartesian interpolated-grid output |
| `sph` | Spherical-surface VTK output |
| `sphslice` | Fixed-radius binary angular slice |
| `rst` | Restart checkpoint |

No other output format is registered by the public `Outputs` constructor.

## Shared Configuration

Each `<output#>` block defines a stream with `file_type`, `dt` or `dcycle`,
and format-specific options. Field-bearing formats also require `variable`
unless their modern syntax defines axis variables instead. Common optional
fields include `id`, `ghost_zones`, `gid`, and Cartesian slice positions.

## File Distribution

Supported writers accept shared output, per-rank shards through
`single_file_per_rank = true`, or per-node shards through
`single_file_per_node = true`. The two shard selectors are mutually exclusive.
Node-sharded products are `bin`, uniform 3D active-zone full-volume `cbin`,
modern `pdf`, `sphslice`, and `rst`.

Uniform 3D active-zone full-volume node-sharded `cbin` requires a
`coarsen_factor` that is a power of two between `2` and the shortest MeshBlock
dimension. Every emitted axis extent must be divisible by that factor.
Lower-dimensional, ghost-zone-expanded, static-refinement, AMR, and sliced
`cbin` configurations are rejected before publication.

## Modern PDF Output

Modern PDFs use one to four contiguous scalar axes:

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
weight = mass
```

`scaleN` accepts `linear`, `log`, or `symlog`; symlog axes require a positive
`linthreshN`. `weight` accepts `volume`, `mass`, or `variable`.
`weight = mass` accumulates positive conserved density times cell volume and
aborts on non-finite or non-positive conserved density. `weight = variable`
requires `weight_variable` and accumulates that finite, possibly signed value
times cell volume. Legacy unsharded PDF inputs remain supported. Modern output
uses an ASCII header plus a versioned `AKPDFV2` payload. Use
`vis/python/read_pdf.py` for legacy, shared V2, and sharded V2 products.

The `edot_sph_out`, `edot_sph_in`, `edot_vert_out`, and `edot_vert_in`
diagnostics partition the signed total energy flux. For MHD output, that total
includes the Poynting contribution, so its sign can differ from the gas
velocity. Total-energy channels require an ideal-gas total-energy fluid module.
`edot_sph_mag` requires MHD because it reports magnetic energy transport.

## Spherical Slice Output

`sphslice` samples native state-backed scalar fields or groups on a fixed
spherical radius and writes a `.sph.bin` angular array:

```ini
<output_slice>
file_type = sphslice
id = rho_shell
variable = hydro_w_d
slice_r = 0.5
ntheta = 64
nphi = 128
dt = 0.1
```

The origin-centered spherical surface must fit inside a 3D domain: `slice_r`
must be positive and strictly interior to every domain face. Derived-array
fields are rejected until spherical
interpolation has ghost-zone-safe derived sampling. `sphslice` is distinct
from the existing `file_type = sph` output. Use `vis/python/read_sphslice.py`
for shared, per-rank, and per-node products. Writers and readers reject
non-finite spherical-slice samples, including finite in-memory values that
overflow when narrowed to the serialized float payload.

## Node-Sharded Restart Files

For a restart stream configured with `single_file_per_node = true`, the public
restart filename is a manifest and each node writes a payload below a
`node_########/` directory. Restart from the manifest, never from an individual
payload. Native resume validates the manifest and reads each rank's routed
MeshBlock spans directly from node payloads. It does not create a shared
`.assembled` staging file. Public-manifest rename is the commit point: failures
before it roll back publication, while failures after it must preserve the
resumable manifest and declared payloads.

## New Binary Payload Portability

Modern PDF V2 and `sphslice` payload scalars currently use host-native byte
order. Same-platform reader/writer round trips are supported. Cross-endian
portability is not claimed for these new binary payloads.

## End-Of-Run Output Policy

Configure final writes under `<time>`:

```ini
<time>
final_output_policy = restart_only
output_timing = true
```

`final_output_policy` accepts `all`, `restart_only`, or `none`; the default is
`all`. Under MPI, timed output records report the maximum elapsed write time
across ranks.

## Executable Examples And Readers

The code branch supplies focused input decks in `inputs/io/` and a matching
readback helper at `vis/python/examples/read_io_outputs.py`. See
[IO Outputs And Sharding](../examples/io_outputs_and_sharding.md) for runnable
commands. Automated MPI coverage exercises multiple ranks on one physical
node. Real multi-node IO and restart qualification remains required before
production rollout.
