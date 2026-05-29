# Module: Z4c

The public Z4c module evolves spacetime fields and constraint-monitoring
state, with optional wave extraction, compact-object tracking, horizon dumps,
and refinement helpers. Use shipped Z4c input decks as the starting point for
this expert-facing functionality.

## State And Sources

`src/z4c/z4c.cpp` allocates the named Z4c fields:

```text
z4c_chi
z4c_gxx z4c_gxy z4c_gxz z4c_gyy z4c_gyz z4c_gzz
z4c_Khat
z4c_Axx z4c_Axy z4c_Axz z4c_Ayy z4c_Ayz z4c_Azz
z4c_Gamx z4c_Gamy z4c_Gamz
z4c_Theta
z4c_alpha
z4c_betax z4c_betay z4c_betaz
```

Constraint output names include `con_C`, `con_H`, `con_M`, `con_Z`, and
`con_Mx`/`con_My`/`con_Mz`.

| Source | Responsibility |
| --- | --- |
| `src/z4c/z4c.*` | Allocation and parameter parsing |
| `z4c_calcrhs.cpp`, `z4c_update.cpp`, `z4c_tasks.cpp` | Evolution and task wiring |
| `z4c_gauge.cpp`, `z4c_Sbc.cpp` | Gauge and boundary work |
| `z4c_adm.cpp` | ADM-derived fields |
| `z4c_calculate_weyl_scalars.cpp`, `z4c_wave_extr.cpp` | Wave extraction |
| `z4c_amr.*` | Z4c refinement strategies |
| `compact_object_tracker.*`, `horizon_dump.*` | Tracking and local dumps |

## `<z4c>` Controls Parsed In The Constructor

| Parameter | Default |
| --- | --- |
| `chi_psi_power` | `-4.0` |
| `chi_div_floor` | `-1000.0` |
| `chi_min_floor` | `1e-12` |
| `diss` | `0.0` |
| `eps_floor` | `1e-12` |
| `damp_kappa1`, `damp_kappa2` | `0.0` |
| `lapse_harmonicf` | `1.0` |
| `lapse_harmonic` | `0.0` |
| `lapse_oplog` | `2.0` |
| `lapse_advect` | `1.0` |
| `slow_start_lapse` | `false` |
| `ssl_damping_amp`, `ssl_damping_time`, `ssl_damping_index` | `0.6`, `20.0`, `1` |
| `shift_Gamma`, `shift_advect`, `shift_alpha2Gamma`, `shift_H`, `shift_eta` | `1.0`, `1.0`, `0.0`, `0.0`, `2.0` |
| `use_z4c` | `true` |
| `user_Sbc` | `false` |
| `excise_chi` | `0.0625` |
| `extrap_order` | Clamped by `nghost` into the implemented range |

## Wave Extraction And Tracking

Wave extraction is disabled by default:

| Parameter | Default/role |
| --- | --- |
| `nrad_wave_extraction` | `0`; number of extraction grids |
| `extraction_nlev` | `10`; geodesic grid level |
| `extraction_radius_<n>` | `10`; radius per extraction grid, with zero-based `n = 0 ... nrad_wave_extraction - 1` |
| `waveform_dt` | `1`; output cadence |

When extraction is active the code writes into `waveforms/`. Compact-object
trackers are enrolled by successive `<z4c>/co_<n>_type` parameters; horizon
dumps are enabled by `dump_horizon_<n>`.

`inputs/z4c/z4c_boosted_puncture.athinput` is not currently an exact runnable
reference: it numbers extraction radii as `_1` through `_4`, while the
constructor reads `_0` through `_3`, and it requests
`pgen_name = z4c_one_puncture` rather than the built-in
`z4c_boosted_puncture` route that consumes the configured puncture velocity.

## Z4c-Specific Refinement

`src/z4c/z4c_amr.cpp` parses `<z4c_amr>/method` values `trivial`, `tracker`,
`chi`, and `dchi`. `chi` reads `chi_min` (default `0.2`), and `dchi` reads
`dchi_max` (default `0.01`). Optional `radius_<n>_rad` and
`radius_<n>_reflevel` pairs add radial refinement floors.

This is distinct from the generic `<amr_criterion*>` interface described on
the [Mesh](mesh.md) page; use only values accepted by the current parser and
do not assume every shipped input deck is synchronized with it.

To invoke Z4c-specific refinement in an adaptive run, a problem generator
must install `user_ref_func`, and an `<amr_criterion*>` block must select
`method = user`; that callback then evaluates the `<z4c_amr>` method above.
The boosted-puncture deck demonstrates this bridge, but remains subject to the
generator and extraction-radius defects noted above.

## Shipped Starting Points

- `inputs/tests/linear_wave_z4c.athinput` selects the built-in
  `z4c_linear_wave` generator and outputs `z4c`.
- `inputs/z4c/onepuncture/z4c_onepuncture.athinput` is paired with the custom
  generator `src/pgen/z4c_one_puncture.cpp`; configure with
  `-DPROBLEM=z4c_one_puncture` before running it. It outputs `adm`, `z4c`,
  and constraints.
- More specialized binary and AMR inputs are under `inputs/z4c/`, but some
  currently use values not accepted by the public implementation:
  `z4c_onepuncture_amr.athinput` uses `dchi_max`, and
  `z4c_twopuncture_amr_criterion.athinput` uses `L2_sphere_in_sphere`.
  The implemented `<z4c_amr>/method` values are `trivial`, `tracker`, `chi`,
  and `dchi`. Both named adaptive decks also omit the required
  `<amr_criterion*>` block, so correcting only the Z4c method value does not
  make either input runnable.

## See Also

- [Dynamical GRMHD](dyn_grmhd.md)
- [Coordinates](coordinates.md)
- [Outputs](outputs.md)
