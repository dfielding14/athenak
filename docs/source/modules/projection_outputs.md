# Projection Outputs

## Purpose

`file_type = proj` reduces any supported cell-centered output variable, including
derived variables, along one or two Cartesian axes. A three-dimensional calculation can
therefore produce either a two-dimensional map or a one-dimensional profile. Projections
support physical bounds on every reduced axis and can write integrals, volume-weighted
means, mass-weighted means, or weighted line-of-sight standard deviations.

The production layout is `projection_layout = native_amr`: the simulation writes
additive projected patches at each active leaf MeshBlock's native transverse resolution.
Post-processing composes those patches to the desired display or analysis resolution.
This retains AMR structure and avoids allocating a finest-level image on every MPI rank.

## Production Configuration

An AMR slab diagnostic with a mass-weighted mean and dispersion is:

```ini
<output1>
file_type = proj
projection_layout = native_amr
id = narrow_slab_mass_energy
variable = hydro_w_e
projection_axes = x1
projection_x1_min = 0.0
projection_x1_max = 0.1
weighting = mass
single_file_per_node = true
dt = 0.05

<output2>
file_type = proj
projection_layout = native_amr
id = narrow_slab_energy_stddev
variable = hydro_w_e
projection_axes = x1
projection_x1_min = 0.0
projection_x1_max = 0.1
weighting = mass
statistic = stddev
single_file_per_node = true
dt = 0.05
```

A three-dimensional calculation can instead be reduced to a one-dimensional bounded
profile:

```ini
<output3>
file_type = proj
projection_layout = native_amr
id = bounded_density_profile
variable = hydro_w_d
projection_axes = x1,x2
projection_x1_min = 0.0
projection_x1_max = 0.1
projection_x2_min = -0.1
projection_x2_max = 0.1
weighting = integral
single_file_per_node = true
dt = 0.05
```

| Parameter | Values | Default | Meaning |
| --- | --- | --- | --- |
| `variable` | any supported output variable | required | Stored or derived field to reduce. Groups produce one projected component per field. |
| `projection_layout` | `native_amr`, `uniform` | `native_amr` | Native additive patches for production, or a dense guarded reference image. |
| `projection_axes` | one or two of `x1,x2,x3` | unset | Coordinates to reduce. One axis returns a 2D map; two axes return a 1D profile. |
| `projection_axis` | `x1`, `x2`, `x3` | `x3` | Single-axis alias used when `projection_axes` is absent. |
| `projection_xN_min/max` | real | full selected coordinate | Physical reduction interval on a selected axis. |
| `projection_min/max` | real | full selected coordinate | Single-axis compatibility aliases. |
| `weighting` | `integral`, `volume`, `mass` | `integral` | Reduction measure; `mass` uses Hydro or MHD density. |
| `statistic` | `value`, `stddev` | `value` | Reduced value or weighted standard deviation; `stddev` requires `volume` or `mass`. |
| `single_file_per_node` | boolean | `false` | In native layout, gather only within each node and write additive node shards. Recommended at scale. |
| `single_file_per_rank` | boolean | `false` | In native layout, write rank-local patches without communication; useful when file count is acceptable. |

`ghost_zones`, `gid`, and `slice_x*` are rejected because projection bounds define the
spatial selection.

## Moment Definition

For every native patch pixel and requested variable \(q\), the writer stores additive
moments. Let \(V_{c,p,D}\) be the overlap of active leaf cell \(c\) with the retained
patch pixel \(p\) and the bounded reduced axes \(D\). For an integral, the reader uses

```{math}
S_1(p) = \frac{1}{A_p}\sum_c q_c V_{c,p,D},
```

where \(A_p\) is the retained pixel measure implicit in patch composition. For normalized
volume or mass statistics the stored records are

```{math}
W(p)  = \sum_c w_c V_{c,p,D},\qquad
S_1(p)= \sum_c w_c q_c V_{c,p,D},\qquad
S_2(p)= \sum_c w_c q_c^2 V_{c,p,D},
```

with \(w_c=1\) for volume weighting and \(w_c=\rho_c\) for mass weighting. After every
overlapping patch and shard has been added, the reader evaluates

```{math}
\langle q\rangle_W = S_1/W,\qquad
\sigma_W(q)=\left[S_2/W-(S_1/W)^2\right]^{1/2}.
```

Storing moments rather than finalized means is essential: independent line-of-sight
MeshBlocks, MPI ranks, and per-node shards may contribute to the same projected location.

## AMR And Bounds

Native output writes one patch for each active leaf MeshBlock intersecting the selected
projection interval. Patch resolution follows the leaf resolution in the retained axes:
coarse outer regions remain coarse, while a refined central region remains fine.
Patches at different levels may overlap in the projected plane because they correspond
to different depths along the removed coordinates; their moments are added during
composition.

Projection bounds are applied before host staging. A slab such as
`projection_x1_min = 0.0`, `projection_x1_max = 0.1` excludes non-intersecting
MeshBlocks and stages only intersecting `x1` cell ranges within retained blocks.
Derived fields are still calculated through the normal device derived-variable path
before the selected cells are copied into projection patches.

In native layout, each completed local patch is serialized immediately rather than
retaining a second in-memory copy of the entire local patch collection. Per-node mode
therefore communicates only serialized records within a node; it never constructs a
global dense image during output.

## Guarded Uniform Reference Layout

Dense uniform projection is retained for bounded, small quicklook products and
numerical reference comparisons only. It must be explicitly requested and explicitly
sized:

```ini
<output9>
file_type = proj
projection_layout = uniform
id = reference_column
variable = hydro_w_d
projection_axis = x3
weighting = integral
projection_level = 1
projection_max_megabytes = 256
dt = 0.05
```

`projection_level` has no default in this layout. `projection_max_megabytes` limits the
dense arrays together with selected-cell staging and local patch composition storage per
rank, and defaults to 256 MiB. Uniform output uses a shared rank-zero write and rejects
shard flags; at least one reduced axis must specify bounds strictly inside the domain.
It is not the production path for large two-dimensional maps.

## Files, Sharding, And Composition

Native files are binary `.proj` streams containing an ASCII metadata preheader, the
input parameter dump, and patch records. Each patch stores its leaf level, retained
physical bounds, dimensions, `W`, `S1`, and `S2` when standard deviation is requested.

| Layout or mode | Path example | Communication |
| --- | --- | --- |
| Native shared | `proj/run.map.00000.proj` | Patch payload gathered globally; small jobs only. |
| Native per node | `proj/node_00000000/run.map.00000.proj` | Patch payload gathered inside each physical node only. |
| Native per rank | `proj/rank_00000000/run.map.00000.proj` | No projection communication; more files. |
| Uniform | `proj/run.map.00000.proj` | Dense global reductions; quicklook/reference only. |

`vis/python/projection.py` discovers sibling rank or node shards automatically and
composes native patches at a requested physical AMR level:

```python
from projection import read_projection, plot_projection

projection = read_projection(
    "proj/node_00000000/blast_projection_amr.narrow_slab_mass_energy.00001.proj",
    level=1,
)
mean_energy = projection["fields"]["eint"]
plot_projection(
    "proj/node_00000000/blast_projection_amr.narrow_slab_mass_energy.00001.proj",
    variable="eint",
    level=1,
    output="narrow_slab_energy.png",
)
```

If `level` is omitted, the reader composes at the finest level present among the
selected patch records. A lower level performs conservative coarsening during
composition.

## Shipped Examples

### Fixed-Grid Hydro Blast: Map And 3D-To-1D Profile

`inputs/hydro/projection/blast_column_density.athinput` uses a three-dimensional
spherical blast. Both plotted products are composed from native projection records: a
full-depth `x3` column-density map and a one-dimensional profile reduced over
\(0 \le x_1 \le 0.1\), \(-0.1 \le x_2 \le 0.1\).

```{image} ../_static/projection_blast_column_density.png
:alt: Fixed-grid blast-wave column density composed from native projection patches
:width: 72%
:align: center
```

```{image} ../_static/projection_blast_bounded_profile.png
:alt: One-dimensional bounded blast projection composed from native projection patches
:width: 72%
:align: center
```

### Four-Level AMR Blast: Fractional-Domain Statistics And Native Binning

`inputs/hydro/projection/blast_amr_weighted.athinput` exercises adaptive refinement and
allows four adaptive levels above the root mesh (`num_levels = 5`, because the root
counts as level zero). It projects only the slab \(0 \le x_1 \le 0.1\) in a
\([-0.5,0.5]^3\) domain. The figure composes its additive moments at the finest native
level present: the first two panels show a mass-weighted internal-energy mean and
line-of-sight deviation, while the third reports the finest native leaf patch
contributing to each projected bin. White outlines mark finest-level MeshBlock
footprints, making the adaptive output representation visible independently of the
projected field. At the displayed evolved output, the slab contains 776 additive patch
records spanning native leaf levels 1 through 4; composition at level 4 gives the
`384 x 384` display grid without requiring that dense grid during simulation output.

```{image} ../_static/projection_blast_amr_slab_statistics.png
:alt: Four-level native-AMR blast slab statistics and contributing leaf-patch levels
:width: 100%
:align: center
```

### Derived Variable In A True 3D MHD Problem

Orszag-Tang is fundamentally two-dimensional. The MHD example instead uses
`inputs/mhd/projection/field_loop3d_current.athinput`, a rotated three-dimensional
field loop (`iprob = 4`). It projects derived `mhd_j2` through
\(-0.2 \le x_3 \le 0.2\) and also writes a volume-weighted deviation of derived
`mhd_bmag`.

```{image} ../_static/projection_field_loop3d_current.png
:alt: Native projection of bounded derived current density in a 3D field loop
:width: 72%
:align: center
```

## Efficiency Validation

Setting `ATHENAK_OUTPUT_IO_STATS=1` prints one telemetry record per projection shard
writer. It reports patch count, staged selected cells, shard payload bytes,
within-shard communicated bytes, estimated peak projection-owned bytes including local
staging/serialization buffers, and elapsed output time.

The regression suite includes:

| Test | Coverage |
| --- | --- |
| `tst/test_suite/nr/test_nr_projection_cpu.py` | Native fixed-grid moments, bounded multi-axis profiles, derived variables, standard deviations, native-to-uniform equivalence, uniform bounds and memory guards, SMR and AMR composition. |
| `tst/test_suite/nr/test_nr_projection_mpicpu.py` | Per-node shard reconstruction, bounded AMR profile correctness, and central-versus-broad refinement scaling metrics. |
| `tst/inputs/projection_scale_central.athinput` | Mixed-level central refinement with full and bounded line-of-sight statistics. |
| `tst/inputs/projection_scale_broad.athinput` | Broad refinement workload for patch/payload growth comparisons. |
