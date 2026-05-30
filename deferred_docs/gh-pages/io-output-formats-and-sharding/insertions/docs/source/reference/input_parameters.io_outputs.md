## IO Output Formats And Runtime Policy

Merge this reviewed section into `docs/source/reference/input_parameters.md`
near the output/time parameter catalogue after the IO code branch has merged.
It records parameters introduced or clarified by the IO feature branch.

### `<time>` Output Finalization Parameters

| Parameter | Type | Default | Accepted values | Description |
| --- | --- | --- | --- | --- |
| `output_timing` | boolean | `false` | `true`, `false` | Report one output-duration line per written stream. Under MPI, `elapsed_max_s` is reduced as the maximum across ranks. |
| `final_output_policy` | string | `all` | `all`, `restart_only`, `none` | Select streams eligible for final writes when the run reaches a terminal or wall-clock condition. |

Invalid `final_output_policy` values are rejected during driver
initialization. Final outputs advance their output number normally. Loading a
checkpoint already at the terminal state does not overwrite the terminal
restart file.

### Shard Selection In Supported `<output#>` Blocks

| Parameter | Type | Default | Formats in this feature | Description |
| --- | --- | --- | --- | --- |
| `single_file_per_rank` | boolean | `false` | Existing supported formats; extended readers for `bin`, `cbin`, modern `pdf`, and `sphslice`; restart retained | Write one output shard per MPI rank. |
| `single_file_per_node` | boolean | `false` | `bin`, full-volume `cbin`, modern `pdf`, `sphslice`, `rst` | Write one shard per MPI shared-memory node. |

`single_file_per_rank = true` and `single_file_per_node = true` in the same
block are mutually exclusive and cause input rejection. Sliced node-sharded
`cbin` is not a promoted workflow because shared sliced `cbin` readback has a
pre-existing meshblock-extent defect requiring a separate fix.

### Modern `<output#>` PDF Parameters

Use `file_type = pdf` with contiguous dimensions from axis 1 through axis 4.

| Parameter | Type | Required/default | Accepted values or validation |
| --- | --- | --- | --- |
| `variable_1` ... `variable_4` | string | Axis 1 required for modern syntax | Scalar output variable; dimensions must not contain gaps. |
| `nbin1` ... `nbin4` | integer | Required for each active axis | Positive count. |
| `bin1_min` ... `bin4_min` | real | Required for each active axis | Less than matching maximum. |
| `bin1_max` ... `bin4_max` | real | Required for each active axis | Greater than matching minimum. |
| `scale1` ... `scale4` | string | `linear` | `linear`, `log`, `symlog`; `log` requires positive limits. |
| `linthresh1` ... `linthresh4` | real | Required only with matching `scaleN = symlog` | Positive threshold; rejected when supplied for non-symlog axis. |
| `weight` | string | `volume` | `volume`, `mass`, `variable`. |
| `weight_variable` | string | Required with `weight = variable` | Scalar output variable. |

Legacy unsharded PDF configurations using `variable`, optional `variable_2`,
legacy bin/log keys, and `mass_weighted` remain accepted. A pure legacy
unsharded configuration preserves the legacy text file layout. New modern or
sharded configurations use the V2 PDF representation.

### `<output#>` Spherical Slice Parameters

Use `file_type = sphslice` for fixed-radius binary angular samples:

| Parameter | Type | Required/default | Validation |
| --- | --- | --- | --- |
| `variable` | string | Required | Native scalar field such as `hydro_w_d` or `mhd_w_d`; derived-array fields are rejected until ghost-zone-safe angular interpolation is implemented. |
| `slice_r` | real | Required | Requested radius must be valid for the configured three-dimensional domain. |
| `ntheta` | integer | `64` | Must be at least 2. |
| `nphi` | integer | `128` | Must be at least 2. |
| `single_file_per_rank` | boolean | `false` | Optional rank-sharded angular ownership. |
| `single_file_per_node` | boolean | `false` | Optional node-sharded angular ownership. |

This product is separate from the pre-existing `file_type = sph` output.

### Diagnostic Scalar Names Available To Modern PDFs

In addition to existing scalar output fields, modern PDF axes and variable
weights can use:

```text
coord_x coord_y coord_z coord_r coord_theta coord_phi
coord_cyl_R coord_cyl_phi coord_cyl_z coord_costheta coord_abscostheta
mdot_sph mdot_sph_out mdot_sph_in mdot_vert mdot_vert_out mdot_vert_in
edot_sph edot_sph_out edot_sph_in edot_vert edot_vert_out edot_vert_in
vel_sph_r vel_sph_theta vel_sph_phi vel_cyl_R vel_cyl_phi
edot_sph_kin edot_sph_th edot_sph_mag
hydro_u_s_N hydro_w_s_N mhd_u_s_N mhd_w_s_N
```

For passive scalar names, `N` is a zero-based scalar index.

The generic `mdot_*`, `edot_*`, and `vel_*` diagnostic names apply to
single-fluid Hydro or MHD only. They are rejected for `<ion-neutral>`
two-fluid configurations because module-qualified flux semantics are not yet
defined.
