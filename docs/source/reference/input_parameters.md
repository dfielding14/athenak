# Public Input Parameter Reference

This reference covers shared public input blocks verified against the
implementation on the public documentation baseline. Problem-generator
parameters are intentionally documented with their shipped input decks and
source files, because they change with the generator selected at build or run
time.

```{warning}
Parameters from development-only CGM cooling, tiled turbulence, or staged AMR
turbulence work are not public configuration options and are not listed here.
```

## Core Run Blocks

### `<job>`

| Parameter | Requirement | Meaning | Source |
| --- | --- | --- | --- |
| `basename` | Required when output is configured | Prefix for generated files | `src/outputs/outputs.cpp` |

### `<mesh>`

| Parameter | Requirement/default | Meaning |
| --- | --- | --- |
| `x1min`, `x1max`, `x2min`, `x2max`, `x3min`, `x3max` | Required | Physical domain extent |
| `nx1`, `nx2`, `nx3` | Required | Global active-zone counts; unused dimensions use `1` |
| `nghost` | Default `2` | Ghost zones on each side |
| `ix1_bc`, `ox1_bc` | Required | Inner/outer x1 boundary flags |
| `ix2_bc`, `ox2_bc` | Required for multidimensional meshes | Inner/outer x2 boundary flags |
| `ix3_bc`, `ox3_bc` | Required for three-dimensional meshes | Inner/outer x3 boundary flags |

Boundary flags are parsed by `src/bvals/bvals.cpp`. `inflow` additionally
requires a problem generator that initializes the constant inflow state.
`shear_periodic` is limited to appropriate x1 faces in shearing-box setups.

### `<meshblock>`

| Parameter | Default | Meaning | Source |
| --- | --- | --- | --- |
| `nx1` | `<mesh>/nx1` | Per-MeshBlock x1 active zones | `src/mesh/mesh.cpp` |
| `nx2` | `<mesh>/nx2` | Per-MeshBlock x2 active zones | `src/mesh/mesh.cpp` |
| `nx3` | `<mesh>/nx3` | Per-MeshBlock x3 active zones | `src/mesh/mesh.cpp` |

### `<time>`

| Parameter | Requirement/default | Meaning | Source |
| --- | --- | --- | --- |
| `evolution` | Required | `static`, `kinematic`, or `dynamic` | `src/driver/driver.cpp` |
| `integrator` | Default `rk2` for evolving runs | `rk1`, `rk2`, `rk3`, `rk4`, `imex2`, `imex2+`, or `imex3` | `src/driver/driver.cpp` |
| `tlim` | Required for `kinematic`/`dynamic`; not read for `static` | Ending simulation time | `src/driver/driver.cpp` |
| `nlim` | Default `-1` for `kinematic`/`dynamic`; not read for `static` | Cycle limit; negative disables this stop condition | `src/driver/driver.cpp` |
| `ndiag` | Default `1` for `kinematic`/`dynamic`; not read for `static` | Diagnostic print cadence in cycles | `src/driver/driver.cpp` |
| `cfl_number` | Used by active physics modules | CFL multiplier for timestep estimates | Module timestep implementation |

## Fluid Blocks

### `<hydro>` And `<mhd>`

| Parameter | Default/requirement | Notes |
| --- | --- | --- |
| `eos` | Required | Selects the equation-of-state implementation |
| `gamma` | Required for ideal EOS | Ratio of specific heats |
| `iso_sound_speed` | Required for isothermal EOS | Isothermal sound speed |
| `nscalars` | Default `0` | Number of passive scalars |
| `reconstruct` | Default `plm` | `dc`, `plm`, `ppm4`, `ppmx`, or `wenoz` |
| `rsolver` | Required for current dynamic, kinematic, and static fluid construction | Newtonian Hydro dynamic: `llf`, `hlle`, `hllc`, `roe`; Newtonian MHD dynamic: `llf`, `hlle`, `hlld`; Newtonian kinematic and current Newtonian `static`: `advect`; SR/GR dynamic solver sets are narrower, as listed on the module pages |
| `fofc` | Default `false` | First-order flux correction |
| `viscosity`, `conductivity` | Optional | Enables corresponding diffusion when present |
| `ohmic_resistivity` | Optional in `<mhd>` | Enables resistivity when present |

`ppm4`, `ppmx`, and `wenoz` require `<mesh>/nghost >= 3`. FOFC with
`plm` requires `nghost >= 3`; FOFC with any of those higher-order
reconstructions requires `nghost >= 4`.

Supplying both `<hydro>` and `<mhd>` selects the ion-neutral two-fluid route
and requires an `<ion-neutral>` block. Conversely, `<ion-neutral>` requires
both fluid blocks and is rejected with dynamical spacetime inputs
(`<adm>` or `<z4c>`).

The driver accepts `time/evolution = static`, but the current Newtonian
Hydro/MHD constructors still enter their non-dynamic solver-selection path and
therefore require `rsolver = advect`. SR/GR constructors reject non-dynamic
selection, including `static`. Treat this as executable current behavior
rather than assuming a solver-free static initialization mode.

## Mesh Refinement

| Block/parameter | Requirement/default | Meaning | Source |
| --- | --- | --- | --- |
| `<mesh_refinement>/refinement` | Default `none` | `none`, `static`, or `adaptive` mesh refinement mode | `src/mesh/mesh.cpp` |
| `<mesh_refinement>/num_levels` | Adaptive only; default `1` | Number of adaptive refinement levels | `src/mesh/build_tree.cpp` |
| `<mesh_refinement>/max_nmb_per_rank` | Required for adaptive refinement | Per-rank MeshBlock capacity | `src/mesh/build_tree.cpp` |
| `<mesh_refinement>/ncycle_check` | Adaptive only; default `1` | AMR check cadence | `src/mesh/mesh_refinement.cpp` |
| `<mesh_refinement>/refinement_interval` | Adaptive only; default `5` | Minimum cycles between adaptive changes | `src/mesh/mesh_refinement.cpp` |
| `<mesh_refinement>/prolong_primitives` | Multilevel optional | Select primitive-variable prolongation for static or adaptive refinement transfers | `src/mesh/mesh_refinement.cpp` |

Static refinement reads `<refined_region*>` blocks and does not require the
adaptive allocation controls above. Adaptive setups may also specify initial
`<refined_region*>` blocks before subsequent criterion-driven evolution.
Adaptive refinement reads one or more `<amr_criterion*>` blocks. Each
criterion requires `method`, with supported
methods `min_max`, `slope`, `second_deriv`, `location`, and `user`.
Non-location/non-user criteria require `variable`; selecting `method = user`
requires the selected problem generator to enroll `user_ref_func`, because
the adaptive path invokes that callback. The
optional thresholds and locations are `value_min`, `value_max`, `location_x1`,
`location_x2`, `location_x3`, and `location_rad`.

## Output Blocks

Any block whose name starts with `<output` defines an output stream.

| Parameter | Requirement/default | Meaning |
| --- | --- | --- |
| `file_type` | Required | `tab`, `hst`, `log`, `vtk`, `pvtk`, `trk`, `cbin`, `pdf`, `bin`, `cart`, `sph`, or `rst` |
| `dt` or `dcycle` | One cadence required | Output interval by time or cycle |
| `variable` | Required except for `hst`, `rst`, and `log` | Output variable or group |
| `id` | Defaults to `variable` | Output identifier in filenames |
| `ghost_zones` | Default `false` | Include ghost zones |
| `gid` | Default `-1` | Restrict output to one MeshBlock where valid |
| `slice_x1`, `slice_x2`, `slice_x3` | Optional | Coordinate slices within the domain |
| `data_format` | Default `%12.5e` | Formatted-output numeric formatting |
| `file_number`, `last_time` | Maintained/defaulted by output manager | Numbering and cadence state |

Format-specific controls:

| Format | Additional parameters |
| --- | --- |
| `hst` | `user_hist_only` (requires `<problem>/user_hist = true`) |
| `bin`, `rst` | `single_file_per_rank` |
| `cbin` | `single_file_per_rank`, required `coarsen_factor`, optional `compute_moments` |
| `pdf` | Required `bin_min`, `bin_max`, `nbin`; optional `logscale`, `mass_weighted`; a configured second dimension also needs `variable_2`, `nbin2 > 1`, and valid second-bin scaling controls |

Only one `hst`, one `rst`, and one `log` stream may be configured in a run.

## Public Source Terms

Source-term blocks are attached to the active module:
`<hydro_srcterms>`, `<mhd_srcterms>`, or `<rad_srcterms>`.

| Parameter | Default/requirement | Applies when enabled |
| --- | --- | --- |
| `const_accel` | Default `false` | Requires `const_accel_val` and `const_accel_dir` |
| `ism_cooling` | Default `false` | Requires `hrate` and a `<units>` block used by cooling conversion |
| `rel_cooling` | Default `false` | Requires `crate_rel`; `cpower_rel` defaults to `1.0` |
| `rad_beam` | Default `false` | Requires `dii_dt`, `pos_1..3`, `dir_1..3`, `width`, and `spread` |

The separate `<turb_driving>` block creates the turbulence driver:

| Parameter | Default |
| --- | --- |
| `nlow` | `1` |
| `nhigh` | `2` |
| `driving_type` | `0` |
| `expo` | `5.0/3.0` |
| `exp_prp` | `5.0/3.0` |
| `exp_prl` | `0.0` |
| `dedt` | `0.0` |
| `tcorr` | `0.0` |

See [Source Terms](../modules/srcterms.md) and
[Driven Turbulence](../examples/turbulence.md) for verified use.

## Problem Generator Parameters

The default executable dispatches the built-in value of `<problem>/pgen_name`
listed in [Problem Generators](../modules/pgen.md). A custom executable built
with `-DPROBLEM=<file-stem>` obtains its `<problem>` controls from
`src/pgen/<file-stem>.cpp`.

Use the shipped input deck for the selected generator as the parameter
template. For example:

| Workflow | Input template | Generator selection |
| --- | --- | --- |
| Sod shock tube | `inputs/hydro/sod.athinput` | Built-in `shock_tube` |
| Orszag-Tang vortex | `inputs/mhd/orszag_tang.athinput` | Built-in `orszag_tang` |
| Turbulence | `inputs/hydro/turb.athinput` | Build `-DPROBLEM=turb` |
| Blast wave | `inputs/hydro/blast_hydro.athinput` | Build `-DPROBLEM=blast` |

## Specialized Blocks And Owning Sources

These public blocks are specialized; their module pages and shipped decks
should be consulted together with their owning implementation.

| Block | Owning public source | Representative inputs |
| --- | --- | --- |
| `<coord>`, `<units>` | `src/coordinates/`, `src/units/` | `inputs/grhydro/`, `inputs/grmhd/` |
| `<radiation>` | `src/radiation/` | `inputs/radiation/` |
| `<particles>` | `src/particles/` | `inputs/particles/` |
| `<ion-neutral>` | `src/ion-neutral/` | `inputs/ion-neutral/` |
| `<shearing_box>` | `src/shearing_box/` | `inputs/shearing_box/` |
| `<z4c>`, `<z4c_amr>` | `src/z4c/` | `inputs/z4c/` |

This boundary avoids presenting unverified, generator-specific, or
development-only settings as universal public parameters.
