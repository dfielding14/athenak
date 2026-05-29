# Module: Outputs

The output manager scans each `<output*>` input block, applies its cadence and
format-specific parameters, and constructs the corresponding writer from
`src/outputs/outputs.cpp`.

## Registered Formats

| `file_type` | Writer | Public purpose |
| --- | --- | --- |
| `tab` | `FormattedTableOutput` | Formatted one-dimensional table slice |
| `hst` | `HistoryOutput` | Integrated time-series diagnostics |
| `log` | `EventLogOutput` | Event-counter log |
| `vtk` | `MeshVTKOutput` | Mesh fields in legacy VTK |
| `pvtk` | `ParticleVTKOutput` | Particle VTK output |
| `trk` | `TrackedParticleOutput` | Tracked-particle stream |
| `cbin` | `CoarsenedBinaryOutput` | Coarsened binary mesh data |
| `pdf` | `PDFOutput` | One- or two-dimensional histogram |
| `bin` | `MeshBinaryOutput` | Native binary mesh data |
| `cart` | `CartesianGridOutput` | Cartesian interpolated-grid output |
| `sph` | `SphericalSurfaceOutput` | Spherical-surface VTK output |
| `rst` | `RestartOutput` | Restart checkpoint |

No other output format is registered by the public `Outputs` constructor.

## Shared Configuration

```ini
<output1>
file_type = vtk
variable  = hydro_w
dt        = 0.1
```

| Parameter | Requirement/default | Notes |
| --- | --- | --- |
| `file_type` | Required | One registered value above |
| `dt` or `dcycle` | Required for an active stream | Time or cycle cadence |
| `variable` | Required except `hst`, `rst`, and `log` | Must name a supported output variable/group for the enabled physics |
| `id` | Defaults to `variable` | Used in filenames for field streams |
| `ghost_zones` | Default `false` | Include ghost cells where applicable |
| `gid` | Default `-1` | Restrict to one MeshBlock; invalid for a one-block mesh |
| `slice_x1`, `slice_x2`, `slice_x3` | Optional | Coordinate cuts inside the domain |
| `data_format` | Default `%12.5e` | Formatted table/history/PDF bin output |

The constructor permits no more than one history stream, one event-log stream,
or one restart stream in a run.

## Observed File Layout

The public quickstart and focused audit runs establish these paths:

| Configuration | Generated path pattern |
| --- | --- |
| Sod `tab`, `variable = hydro_w` | `tab/Sod.hydro_w.00000.tab` |
| Sod `hst` | `Sod.hydro.hst` |
| Orszag-Tang `vtk`, `variable = mhd_w` | `vtk/OrszagTang.mhd_w.00000.vtk` |
| `bin`, `variable = hydro_w` | `bin/Sod.hydro_w.00000.bin` |
| `rst` | `rst/Sod.00000.rst` |

For `bin`, `cbin`, and `rst`, `single_file_per_rank = true` changes the
parallel file layout as implemented by their writers.

## Format-Specific Controls

### Table Output

`tab` output is a one-dimensional slice. A multidimensional problem must
provide enough `slice_x*` coordinates to leave one unsliced direction.

### History And Restart

`hst` supports `user_hist_only = true` only when
`<problem>/user_hist = true` and the selected problem generator enrolled a
history callback. Configure restart state with:

```ini
<output3>
file_type = rst
dt        = 0.1
```

The initialization checkpoint is numbered `00000`; restart from a later file
when demonstrating continuation.

### Binary And Coarsened Binary

`bin` accepts `single_file_per_rank`. `cbin` additionally requires
`coarsen_factor` and accepts `compute_moments`:

```ini
<output2>
file_type      = cbin
variable       = hydro_w
dt             = 0.2
coarsen_factor = 2
```

Use [Visualization Utilities](../tools/visualization.md) to convert public
binary output to `.athdf`/XDMF.

### PDF Histograms: Registered But Currently Failing In A Smoke Test

The public `PDFOutput` implements one- or two-dimensional histograms, not a
general N-dimensional histogram interface. Its parser and writer source expose
the following interface:

```ini
<output4>
file_type     = pdf
variable      = hydro_w_d
dt            = 0.1
bin_min       = 0.01
bin_max       = 1.1
nbin          = 64
logscale      = false
mass_weighted = false
```

| Parameter | Requirement/default |
| --- | --- |
| `variable`, `bin_min`, `bin_max`, `nbin` | Required for the first dimension |
| `logscale` | Default `true`; if true, `bin_min` must be positive |
| `mass_weighted` | Default `false`; otherwise cell volume is weighted by density |
| `variable_2` | Optional second quantity; alone does not activate 2-D binning |
| `nbin2` | Default `0`; set greater than `1` to activate a second histogram dimension |
| `bin2_min`, `bin2_max`, `logscale2` | Second-bin controls; `logscale2` defaults to `true`, so a 2-D logarithmic setup needs positive `bin2_min` or must set `logscale2 = false` |

The implementation intends to write `pdf_<id>/basename.bins.pdf` and numbered
`.pdf` data files; the directory includes `variable_2` for two-dimensional
histograms. However, a focused run using the public Sod input plus the
one-dimensional configuration above currently aborts in `src/outputs/pdf.cpp`
when a `Kokkos::deep_copy` tries to copy mismatched view extents. Treat `pdf`
as a registered but presently unvalidated output path, and use `tab`, `bin`,
or `vtk` for working field-output examples until this failure is fixed.

## Implementation Entry Points

| File | Responsibility |
| --- | --- |
| `src/outputs/outputs.cpp` | Block parsing, writer dispatch, compatibility checks |
| `src/outputs/basetype_output.cpp` | Resolves selected variables and gathers output data |
| `src/outputs/derived_variables.cpp` | Computed variables requested through output names |
| `src/outputs/history.cpp` | Integrated diagnostics |
| `src/outputs/restart.cpp` | Checkpoint writing and read-compatible state |
| `src/outputs/pdf.cpp` | One-/two-dimensional histogram implementation |

## See Also

- [Running Simulations](../running.md)
- [Public Input Parameters](../reference/input_parameters.md)
- [Visualization Utilities](../tools/visualization.md)
