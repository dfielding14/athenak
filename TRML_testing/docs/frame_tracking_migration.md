# TRML Frame-Tracking Migration (Modern-Only)

Frame tracking in `src/pgen/TRML.cpp` now accepts only the modern `ft_*` schema plus existing `use_frame_tracking`, `t_frame_tracking_start`, and `z_interface`.

Any legacy frame-tracking key now fails fast at startup with a migration error.

## Key Mapping

| Legacy key | Modern key | Notes |
|---|---|---|
| `boost_every` / `n_frame_track` / `ft_boost_interval` | `ft_apply_every` | Actuation cadence in cycles |
| `max_vz_frame_tracking` | `ft_max_abs_boost` | Absolute cap on command |
| `T_lo` / `T_ft_lo_over_T_cold` | `ft_temp_lo` | Explicit temperature lower bound |
| `T_hi` / `T_ft_hi_over_T_cold` | `ft_temp_hi` | Explicit temperature upper bound |
| `min_tracked_mass` | `ft_weight_floor` | Minimum global tracked weight |
| `ft_use_uppzlin` / `ft_uppzlin` | `ft_use_uppzlim` / `ft_uppzlim` | Typo fix |

## Removed Concepts

| Removed key | Replacement |
|---|---|
| `history_len`, `ft_history_length` | Use EMA with `ft_tau_avg` |
| `shift_smooth_beta` | Use EMA and/or `ft_max_boost_change` |
| `ft_alpha_z`, `ft_alpha_v`, `ft_z_target` | Use timescale-based control with `z_interface`, `ft_tau_relax`, `ft_tau_vel` |
| `ft_weight_mode`, `ft_apply_to` | Not supported in TRML modern controller |

## Minimal Modern Example

```ini
use_frame_tracking = true
t_frame_tracking_start = 1.0
z_interface = 0.0

ft_apply_every = 10
ft_mode = pd
ft_tau_avg = 0.5
ft_tau_relax = 10.0
ft_tau_vel = 2.5

ft_temp_lo = 0.032
ft_temp_hi = 0.056
ft_weight_floor = 1.0e-6

ft_max_abs_boost = 0.01
ft_max_boost_change = 0.01
```
