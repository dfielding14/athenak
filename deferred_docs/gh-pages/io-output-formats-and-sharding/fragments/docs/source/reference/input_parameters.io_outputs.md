## Output Blocks

Any block whose name starts with `<output` defines an output stream.

| Parameter | Requirement/default | Meaning |
| --- | --- | --- |
| `file_type` | Required | `tab`, `hst`, `log`, `vtk`, `pvtk`, `trk`, `cbin`, `pdf`, `bin`, `cart`, `sph`, `sphslice`, or `rst` |
| `dt` or `dcycle` | One cadence required | Output interval by time or cycle |
| `variable` | Required for ordinary field-bearing formats | Output variable or group |
| `id` | Defaults to `variable` | Output identifier in filenames |
| `ghost_zones` | Default `false` | Include ghost zones where supported |
| `gid` | Default `-1` | Restrict output to one MeshBlock where valid |
| `single_file_per_rank` | Default `false` | Write supported streams as rank shards |
| `single_file_per_node` | Default `false` | Write `bin`, uniform 3D active-zone full-volume `cbin`, modern `pdf`, `sphslice`, or `rst` as node shards |

The two shard selectors are mutually exclusive.

### `<time>` Output Finalization Parameters

| Parameter | Type | Default | Accepted values | Description |
| --- | --- | --- | --- | --- |
| `output_timing` | boolean | `false` | `true`, `false` | Report one output-duration line per written stream. Under MPI, `elapsed_max_s` is reduced as the maximum across ranks. |
| `final_output_policy` | string | `all` | `all`, `restart_only`, `none` | Select streams eligible for final writes when the run reaches a terminal or wall-clock condition. |

### Coarsened-Binary `<output#>` Parameters

For `file_type = cbin`, `coarsen_factor` is required and must be a power of
two between `2` and the shortest MeshBlock dimension, inclusive.
Every emitted axis extent must be divisible by that factor. `compute_moments`
defaults to `false`. The supported producer matrix is uniform
three-dimensional active-zone full-volume output.

### Modern `<output#>` PDF Parameters

Modern PDFs use contiguous axis numbers from 1 through 4:

| Parameter | Required/default | Accepted values or validation |
| --- | --- | --- |
| `variable_1` ... `variable_4` | Axis 1 required for modern syntax | Scalar output variable; dimensions must not contain gaps. |
| `nbin1` ... `nbin4` | Required for each active axis | Positive count. |
| `bin1_min` ... `bin4_min`, `bin1_max` ... `bin4_max` | Required for each active axis | Minimum must be less than maximum. |
| `scale1` ... `scale4` | `linear` | `linear`, `log`, `symlog`; `log` requires positive limits. |
| `linthresh1` ... `linthresh4` | Required for a matching symlog axis | Positive threshold. |
| `weight` | `volume` | `volume`; `mass` uses conserved density times cell volume and aborts on non-finite or non-positive conserved density; `variable` uses a finite, possibly signed `weight_variable` value times cell volume. |
| `weight_variable` | Required with `weight = variable` | Scalar output variable. |
| `max_writer_allocation_bytes` | `536870912` | Positive per-writer allocation cap. |

Legacy unsharded PDF configurations using `variable`, optional `variable_2`,
legacy bin/log keys, and `mass_weighted` remain accepted.

### `<output#>` Spherical Slice Parameters

For `file_type = sphslice`, `variable` is required and must name a native
state-backed scalar or multi-field group. The origin-centered spherical
surface must fit inside a 3D domain. `slice_r` is required, positive, and
strictly interior to every domain face. `ntheta` defaults to `64`, `nphi`
defaults to `128`, and both must be
at least `2`.
`max_writer_allocation_bytes` defaults to `536870912`.

### Generic Diagnostic Scalar Names

Ordinary scalar output streams and modern PDF axes or variable weights can
use coordinate projections plus the `mdot_*`, `edot_*`, and `vel_*` scalar
families. These generic fluid diagnostics apply to Newtonian single-fluid
Hydro or MHD only. They are rejected for `<ion-neutral>` two-fluid
configurations because module-qualified semantics are not defined. Total-energy
channels require an ideal-gas total-energy fluid module. `edot_sph_mag`
requires MHD because it reports magnetic energy transport.

For the formulas below, $\rho$ is the conserved density, $\mathbf{v}$ is the
fluid velocity, $u$ is the internal-energy density, $\gamma$ is the ideal-gas
adiabatic index, and $\mathbf{B}$ is the cell-centered magnetic field. In Hydro,
set $\mathbf{B}=0$.

| Diagnostic | Newtonian single-fluid definition |
| --- | --- |
| `mdot_sph` | $\rho v_r$ |
| `mdot_vert` | $\operatorname{sign}(z)\rho v_z$ |
| `edot_sph_kin` | $\frac{1}{2}\rho |\mathbf{v}|^2 v_r$ |
| `edot_sph_th` | $\gamma u v_r$, the enthalpy contribution |
| `edot_sph_mag` | $|\mathbf{B}|^2 v_r - (\mathbf{v}\cdot\mathbf{B})B_r$ |
| `edot_sph` | $(\frac{1}{2}\rho |\mathbf{v}|^2 + \gamma u + |\mathbf{B}|^2)v_r - (\mathbf{v}\cdot\mathbf{B})B_r$ |
| `edot_vert` | $\operatorname{sign}(z)[(\frac{1}{2}\rho |\mathbf{v}|^2 + \gamma u + |\mathbf{B}|^2)v_z - (\mathbf{v}\cdot\mathbf{B})B_z]$ |

Each `_out` channel is `max(base_value, 0)` and each `_in` channel is
`min(base_value, 0)`, so inflow remains signed. Radial projections are
canonically zero at the origin. Cylindrical radial and azimuthal projections
are canonically zero on the cylindrical axis. Vertical diagnostics are
canonically zero on the midplane. These diagnostics abort on non-positive
density or any non-finite input or derived value.
