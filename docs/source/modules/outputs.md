# Module: Outputs

## Overview
The Outputs module builds one output handler per `<output[n]>` block, stores
configuration in `OutputParameters`, and writes mesh/particle/diagnostic files
through `BaseTypeOutput` subclasses constructed by `Outputs::Outputs`.
([src/outputs/outputs.cpp:58], [src/outputs/outputs.hpp:132])

## Source Location
`src/outputs/`

## Output Formats (from `Outputs::Outputs`)

### Mesh data

| file_type | Class | Extension | Source |
|-----------|-------|-----------|--------|
| `vtk` | `MeshVTKOutput` | `.vtk` | `vtk_mesh.cpp` ([src/outputs/outputs.cpp:215], [src/outputs/vtk_mesh.cpp:55]) |
| `bin` | `MeshBinaryOutput` | `.bin` | `binary.cpp` ([src/outputs/outputs.cpp:345], [src/outputs/binary.cpp:51]) |
| `cbin` | `CoarsenedBinaryOutput` | `.cbin` | `coarsened_binary.cpp` ([src/outputs/outputs.cpp:224], [src/outputs/coarsened_binary.cpp:301]) |
| `rst` | `RestartOutput` | `.rst` | `restart.cpp` ([src/outputs/outputs.cpp:350], [src/outputs/restart.cpp:169]) |
| `tab` | `FormattedTableOutput` | `.tab` | `formatted_table.cpp` ([src/outputs/outputs.cpp:204], [src/outputs/formatted_table.cpp:59]) |

### Particle data

| file_type | Class | Extension | Source |
|-----------|-------|-----------|--------|
| `pvtk` | `ParticleVTKOutput` | `.part.vtk` | `vtk_prtcl.cpp` ([src/outputs/outputs.cpp:218], [src/outputs/vtk_prtcl.cpp:81]) |
| `trk` | `TrackedParticleOutput` | `.trk` | `track_prtcl.cpp` ([src/outputs/outputs.cpp:221], [src/outputs/track_prtcl.cpp:96]) |

### Diagnostics

| file_type | Class | Extension | Source |
|-----------|-------|-----------|--------|
| `hst` | `HistoryOutput` | `.hst` | `history.cpp` ([src/outputs/outputs.cpp:207], [src/outputs/history.cpp:343]) |
| `log` | `EventLogOutput` | `.log` | `eventlog.cpp` ([src/outputs/outputs.cpp:211], [src/outputs/eventlog.cpp:71]) |
| `pdf` | `PDFOutput` | `.pdf` + `.header.pdf` | `pdf.cpp` ([src/outputs/outputs.cpp:232], [src/outputs/pdf.cpp:272]) |

## Output Block Parameters

### Required keys
- `file_type`: selects one of the formats above. ([src/outputs/outputs.cpp:88])
- `dt` or `dcycle`: outputs are only registered if `dt > 0` or `dcycle > 0`.
  ([src/outputs/outputs.cpp:72])
- `variable`: required for `vtk`, `bin`, `cbin`, `tab`, `pvtk`, and `pdf` outputs.
  For `pdf`, `variable` is a legacy alias for `variable_1`. (`hst`, `log`, `rst`,
  and `trk` do not read `variable`.) ([src/outputs/outputs.cpp:90],
  [src/outputs/outputs.cpp:266])

Only one `hst`, `log`, and `rst` block is allowed. ([src/outputs/outputs.cpp:367])

### Common optional keys
| Parameter | Default | Notes |
|-----------|---------|-------|
| `id` | `variable` | Used in filenames for output types that parse `variable`. ([src/outputs/outputs.cpp:90]) |
| `file_number` | `0` | Starting index for file numbering. ([src/outputs/outputs.cpp:86]) |
| `last_time` | `-1.0` | When negative, outputs always trigger; otherwise advanced by `dt`. ([src/outputs/outputs.cpp:72], [src/outputs/restart.cpp:181]) |
| `ghost_zones` | `false` | Include ghost zones for mesh-style outputs that honor `include_gzs`. ([src/outputs/outputs.cpp:106], [src/outputs/vtk_mesh.cpp:82]) |
| `gid` | `-1` | Output only one MeshBlock by ID; bounds are validated. ([src/outputs/outputs.cpp:109]) |
| `slice_x1/x2/x3` | unset | Slice planes are validated against the mesh extent. ([src/outputs/outputs.cpp:124]) |
| `data_format` | `%12.5e` | Numeric format for `tab` and `hst` outputs. ([src/outputs/outputs.cpp:197], [src/outputs/history.cpp:387]) |
| `user_hist_only` | `false` | For `hst`: emit only user history data when `<problem>/user_hist` is enabled. ([src/outputs/outputs.cpp:186]) |

### Coarsened binary (`cbin`)
| Parameter | Default | Notes |
|-----------|---------|-------|
| `coarsen_factor` | none | Required; must divide each MeshBlock dimension or the run aborts. ([src/outputs/outputs.cpp:225], [src/outputs/coarsened_binary.cpp:216]) |
| `compute_moments` | `false` | If true, writes 1st–4th raw moments of each variable in the coarsened cell. ([src/outputs/outputs.cpp:228], [src/outputs/coarsened_binary.cpp:336]) |
| `single_file_per_rank` | `false` | Emit per-rank files in `cbin_<id>_<factor>/rank_XXXXXXXX/`. ([src/outputs/outputs.cpp:225], [src/outputs/coarsened_binary.cpp:319]) |

### Binary + restart (`bin`, `rst`)
| Parameter | Default | Notes |
|-----------|---------|-------|
| `single_file_per_rank` | `false` | Emit per-rank files in `bin/rank_XXXXXXXX/` or `rst/rank_XXXXXXXX/`. ([src/outputs/outputs.cpp:345], [src/outputs/restart.cpp:169]) |

### PDF (N-D)
PDFs support 1–4 dimensions and can be weighted by volume, mass, or a variable.
The parser accepts the N-D form (`variable_1`, `bin1_min`, `nbin1`, ...), and
also keeps a legacy 1D/2D form (`variable`, `bin_min`, `nbin`, `variable_2`, etc.).
([src/outputs/outputs.cpp:266])

| Parameter | Default | Notes |
|-----------|---------|-------|
| `variable_1..variable_4` | required | Names for each dimension. The highest index present defines `pdf_ndim`. ([src/outputs/outputs.cpp:266]) |
| `binN_min/binN_max` | required | Bin limits per dimension. ([src/outputs/outputs.cpp:296]) |
| `nbinN` | required | Number of bins per dimension. ([src/outputs/outputs.cpp:296]) |
| `logscaleN` | `false` | Logarithmic bins per dimension; requires `binN_min > 0`. ([src/outputs/outputs.cpp:303], [src/outputs/outputs.cpp:317]) |
| `weight` | `volume` | `volume`, `mass`, or `variable`. ([src/outputs/outputs.cpp:240]) |
| `weight_variable` | none | Required when `weight=variable`. ([src/outputs/outputs.cpp:254]) |

Legacy 1D/2D keys still work:
- `variable`, `bin_min`, `bin_max`, `nbin`, `logscale` (logscale defaults to `true`).
- `variable_2`, `bin2_min`, `bin2_max`, `nbin2` (default `10`), `logscale2`.
([src/outputs/outputs.cpp:276])

PDF variables must be scalar. Group variables (e.g., `hydro_u`, `mhd_u`) are
rejected, and any variable that expands to multiple components will error at
PDF construction. ([src/outputs/outputs.cpp:170], [src/outputs/pdf.cpp:62])

The histogram weight is always multiplied by the cell volume; `weight=mass`
adds density, and `weight=variable` multiplies by the supplied variable.
([src/outputs/pdf.cpp:238])

The legacy `mass_weighted` flag is rejected; use `weight=mass` instead.
([src/outputs/outputs.cpp:233])

## Variable Names and Derived Quantities
Selectable variables are enumerated in `var_choice` in `outputs.hpp`.
([src/outputs/outputs.hpp:24])

Derived variables are computed in `derived_variables.cpp` and include:
- Cartesian, spherical, and cylindrical coordinate scalars (`coord_*`).
- Spherical and cylindrical velocity components (`vel_sph_*`, `vel_cyl_*`).
- Radial and vertical mass/energy fluxes (`mdot_*`, `edot_*`).

These use `CellCenterX` for `x1/x2/x3` coordinates, so the coordinate-based
quantities assume `x1/x2/x3` map to Cartesian axes. ([src/outputs/derived_variables.cpp:1563])

## File Naming Conventions
- `vtk`: `vtk/<basename>.<id>[.<gid>].<NNNNN>.vtk`. ([src/outputs/vtk_mesh.cpp:57])
- `pvtk`: `pvtk/<basename>.<id>[.<gid>].<NNNNN>.part.vtk`. ([src/outputs/vtk_prtcl.cpp:81])
- `bin`: `bin/<basename>.<id>.<NNNNN>.bin` (or `bin/rank_XXXXXXXX/` when
  `single_file_per_rank=true`). ([src/outputs/binary.cpp:56])
- `cbin`: `cbin_<id>_<factor>/[rank_XXXXXXXX/]<basename>.<id>.<NNNNN>.cbin`.
  ([src/outputs/coarsened_binary.cpp:306])
- `tab`: `tab/<basename>.<id>.<NNNNN>.tab`. ([src/outputs/formatted_table.cpp:60])
- `rst`: `rst/<basename>.<NNNNN>.rst` (or `rst/rank_XXXXXXXX/`).
  ([src/outputs/restart.cpp:169])
- `hst`: `<basename>.<physics>.hst` (physics suffix per module). ([src/outputs/history.cpp:345])
- `log`: `<basename>.log`. ([src/outputs/eventlog.cpp:76])
- `trk`: `trk/<basename>.trk`. ([src/outputs/track_prtcl.cpp:99])
- `pdf`: `pdf_<id>[_<var2>...]/<basename>.header.pdf` and
  `pdf_<id>[_<var2>...]/<basename>.<NNNNN>.pdf`. ([src/outputs/pdf.cpp:281])

## Usage Examples

### Basic VTK Output
```ini
<output1>
file_type = vtk
variable = hydro_u
id = hydro
dt = 0.1
```

### N-D PDF Output (2D)
```ini
<output2>
file_type = pdf
id = phase
variable_1 = temperature
bin1_min = 1.0e-2
bin1_max = 1.0e2
nbin1 = 128
logscale1 = true

variable_2 = hydro_w_d
bin2_min = 1.0e-4
bin2_max = 1.0e1
nbin2 = 128
logscale2 = true

weight = mass
```
