# Frame Tracking Recipes And Migration

Use these blocks as conservative starting configurations. They assume the
problem generator converts frame-aware lab-frame boundary or source values
using `x_lab = x_grid + FrameDisplacement(axis)` and
`v_grid = v_lab - FrameVelocity(axis)`.

## `cloud_tracking_conservative`

This block matches the cooled-cloud example input
`inputs/hydro/cloud_crushing_snr.athinput`:

```ini
<frame_tracking>
tracked_fluid = hydro
start_time = 0.02
apply_every = 1
axes = x1
x1_target = 3.0
target = density
target_min = 5.0
target_max = 200.0
mode = pd
tau_avg = 0.01
tau_relax = 0.05
tau_vel = 0.05
max_abs_boost = 50.0
max_boost_change_mode = per_time
max_boost_change_rate = 5.0
```

Use it only with a boundary/source implementation that explicitly transforms
lab-frame inflow values. The shipped cloud-crushing problem generator provides
that transformation in Sedov-boundary mode.

## `mixing_layer_tracking_conservative`

This block matches `inputs/hydro/TRML/TRML_frame_tracking.athinput`:

```ini
<frame_tracking>
tracked_fluid = hydro
start_time = 0.0
apply_every = 1
axes = x3
x3_target = 0.0
target = temperature
target_min = 0.015
target_max = 0.08
mode = pd
tau_avg = 0.1
tau_relax = 1.0
tau_vel = 0.25
max_abs_boost = 0.2
max_boost_change_mode = per_time
max_boost_change_rate = 0.02
```

The shipped TRML problem generator provides frame-aware `x3` boundaries.
Do not apply this block unchanged to a problem with lab-frame boundaries or
source terms that have not been transformed.

## Removed Configuration Aliases

Configuration aliases are no longer accepted. Startup fails with the
canonical replacement when a removed key is present. Restart state is
separate: older restart records containing `frame_velocity_x*` and
`frame_displacement_x*` remain readable with a non-exact-continuation warning.

| Removed key(s) | Canonical key |
| --- | --- |
| `enable` | `enabled` |
| `ft_apply_every` | `apply_every` |
| `t_start`, `t_frame_tracking_start` | `start_time` |
| `ndiag` | `diagnostic_every` |
| `target_variable`, `derived_variable` | `target` |
| `min`, `ft_temp_lo` | `target_min` |
| `max`, `ft_temp_hi` | `target_max` |
| `center` | `target_center` |
| `reacquire_scale` | `target_scale` |
| `axis`, `directions` | `axes` |
| `z_target`, `ft_z_target` | `x3_target` |
| `z_lower_limit`, `ft_lowzlim`, `ft_use_lowzlim` | `x3_lower_limit` |
| `z_upper_limit`, `ft_uppzlim`, `ft_use_uppzlim` | `x3_upper_limit` |
| `ft_mode` | `mode` |
| `ft_position_signal` | `position_signal` |
| `ft_position_blend` | `position_blend` |
| `ft_tau_avg`, `ft_tau_relax`, `ft_tau_vel`, `ft_tau_int` | `tau_avg`, `tau_relax`, `tau_vel`, `tau_int` |
| `ft_int_max_abs`, `ft_int_leak_tau`, `ft_int_unsat_only` | `int_max_abs`, `int_leak_tau`, `int_unsat_only` |
| `ft_weight_floor`, `ft_min_global_weight` | `weight_floor`, `min_global_weight` |
| `ft_max_abs_boost`, `ft_max_boost_change`, `ft_max_boost_change_rate` | `max_abs_boost`, `max_boost_change`, `max_boost_change_rate` |
| `ft_max_boost_change_mode` | `max_boost_change_mode` |
| `ft_reacquire_expand_factor`, `ft_reacquire_max_expand`, `ft_reacquire_recover_updates`, `ft_reacquire_leak_on_miss` | Corresponding key without the `ft_` prefix |
| `ft_boundary_guard_enable`, `ft_boundary_guard_cells`, `ft_boundary_guard_min_scale` | Corresponding key without the `ft_` prefix |
| `xN_target_position`, `target_xN` | `xN_target` |
| `xN_min_limit`, `xN_low_limit` | `xN_lower_limit` |
| `xN_max_limit`, `xN_upp_limit` | `xN_upper_limit` |
