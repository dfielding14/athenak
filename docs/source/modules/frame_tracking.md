# Module: Frame Tracking

`FrameTracker` is a shared post-timestep helper for simulations that are easier
to run in a moving Galilean frame. It samples user-selected material, estimates
its position and velocity, and applies a uniform velocity boost to the active
fluid so the selected material remains near a chosen grid-coordinate target.

The tracker is enabled by a `<frame_tracking>` input block. It is owned by
`MeshBlockPack`, scheduled on the existing `after_timeintegrator` task list,
updates its pack pointer after AMR rebuilds, refreshes the selected fluid
timestep after an applied boost, and stores versioned controller state in
restart files.

## Runtime Contract

Frame tracking does not move mesh coordinates or remap AMR blocks. The fluid
variables are stored in the grid frame after each boost, while the tracker keeps
the lab-frame grid-origin displacement as bookkeeping.

Problem generators that impose lab-frame physics at boundaries should use:

```text
x_lab = x_grid + FrameDisplacement(axis)
v_grid = v_lab - FrameVelocity(axis)
```

where `FrameVelocity(axis)` is the lab-frame velocity of the moving grid.

## Public API

Problem generators can include `srcterms/frame_tracker.hpp` and query the
current frame from `MeshBlockPack::pframe_tracker`:

| API | Meaning |
| --- | --- |
| `FrameVelocity(axis)` | Lab-frame grid velocity on axis `0`, `1`, or `2`. |
| `FrameDisplacement(axis)` | Lab-frame grid-origin displacement on axis `0`, `1`, or `2`. |
| `FrameVelocity()` | Three-component `std::array<Real, 3>` frame velocity. |
| `FrameDisplacement()` | Three-component `std::array<Real, 3>` frame displacement. |

The pointer is null when `<frame_tracking>` is absent or disabled, so pgen
boundary code should keep a zero-frame fallback.

## Input Parameters

### Core Selection

| Parameter | Default | Description |
| --- | ---: | --- |
| `enabled` | `true` | Master switch. |
| `tracked_fluid` | inferred when unambiguous | Fluid to sample and boost: `hydro` or `mhd`. Required when both exist. |
| `target` | `density` | Selected variable: `density`, `temperature`, `pressure`, `entropy`, `internal_energy`, `scalar`/`scalarN`, `v1`/`v2`/`v3`, or `speed`. |
| `target_min`, `target_max` | unbounded | Inclusive target range for tracked cells. |
| `target_center` | midpoint or geometric midpoint | Center used when expanding the target range for reacquisition. |
| `target_scale` | `auto` | `auto`, `linear`, or `log` expansion around `target_center`. |
| `weight` | `mass` | Moment weighting: `mass`, `volume`, or absolute `target` value. |
| `weight_floor` | `0.0` | Skip actuation when selected material weight is less than or equal to this value. |
| `min_global_weight` | `0.0` | Optional observable-weight threshold for skipping actuation. |

### Axes And Targets

| Parameter | Default | Description |
| --- | ---: | --- |
| `axes` | highest active dimension | Tracking axes: `x1`, `x2`, `x3`, `all`, or a comma-separated list. |
| `track_x1`, `track_x2`, `track_x3` | unset | Per-axis booleans that override `axes`. |
| `x1_target`, `x2_target`, `x3_target` | domain center | Desired tracked-material position on each active axis. |
| `x*_lower_limit`, `x*_upper_limit` | disabled | Optional limits on the next predicted control position. |

### Controller

| Parameter | Default | Description |
| --- | ---: | --- |
| `mode` | `pd` | `velocity`, `position`, or `pd`. |
| `position_signal` | `blend` | Position observable: centroid, selected-band midpoint, or a blend. |
| `position_blend` | `0.7` | Weight of the selected-band midpoint when `position_signal=blend`. |
| `apply_every` | `1` | Number of cycles between tracking updates. |
| `start_time` | `0.0` | Simulation time after which boosts begin. |
| `tau_avg` | `1.0` | Exponential smoothing time for sampled position and velocity. |
| `tau_relax` | `1.0` | Position-feedback timescale. |
| `tau_vel` | `1.0` | Velocity-feedback timescale. |
| `tau_int` | `0.0` | Integral feedback timescale; zero disables integral control. |
| `int_max_abs` | `0.0` | Absolute cap on the integral term; zero disables integral control. |
| `int_leak_tau` | `1.0` | Leak time for the integral term. |
| `int_unsat_only` | `true` | Prevents additional integral windup while the boost is saturated. |

### Limits And Reacquisition

| Parameter | Default | Description |
| --- | ---: | --- |
| `max_abs_boost` | `0.0` | Absolute per-update boost cap; zero disables this cap. |
| `max_boost_change` | `0.0` | Compatibility per-update slew limit when `max_boost_change_mode=per_apply`. |
| `max_boost_change_rate` | `-1.0` | Recommended per-time boost slew limit when `max_boost_change_mode=per_time`. |
| `max_boost_change_mode` | `per_apply` for compatibility | Interpret the slew limit per update or per unit time. New configurations should set `per_time` explicitly. |
| `reacquire_expand_factor` | `1.0` | Factor by which to expand the target range after missed samples. |
| `reacquire_max_expand` | `1.0` | Maximum target-range expansion. |
| `reacquire_recover_updates` | `1` | Good updates required before backing off one missed-sample level. |
| `reacquire_leak_on_miss` | `true` | Leak the integral term while no valid tracked material is found. |
| `boundary_guard_enable` | `false` | Smoothly reduces boosts when tracked material nears a domain edge. |
| `boundary_guard_cells` | `8` | Width of the edge guard in root-grid cells. |
| `boundary_guard_min_scale` | `0.15` | Minimum multiplicative boost scale at the domain edge. |

### Diagnostics And Restart State

| Parameter | Default | Description |
| --- | ---: | --- |
| `diagnostic_every` | `-1` | Compatibility cadence for rank-0 stdout diagnostic blocks; negative disables printing. |
| `verbose` | `false` | Prints parsed setup at startup. |
| `state_version` | restart only | Restart-state schema version; current value is `1`. |
| `frame_velocity_x*` | restart only | Stored lab-frame grid velocity. |
| `frame_displacement_x*` | restart only | Stored lab-frame grid-origin displacement. |

When any `file_type=hst` output is configured, the tracker writes a dedicated
`<basename>.frame_tracker.hst` file. It does not add columns to hydro/MHD
history or require a pgen `user_hist` callback. Each active axis has five
columns:

| Column | Quantity |
| --- | --- |
| `ft_vf_xN` | Lab-frame grid velocity. |
| `ft_dx_xN` | Lab-frame grid-origin displacement. |
| `ft_pos_xN` | Current selected-material control position. |
| `ft_err_xN` | Filtered position error relative to its target. |
| `ft_dv_xN` | Last applied velocity boost. |

Shared fields are `ft_weight`, `ft_misses`, `ft_recov`, `ft_limit`, and
`ft_skip`. The full three-axis output has 20 fields. Controller values are
global: on MPI runs rank 0 contributes them to the normal history reduction and
other ranks contribute zeros.

Restart output with `state_version=1` stores the complete mutable controller
state, including timing, reacquisition counters, filtering, integral control,
slew-limit state, and per-axis history. A restart containing only
`frame_velocity_x*` and `frame_displacement_x*` is accepted for compatibility,
but it emits one rank-0 warning and cannot be an exact controller continuation.

## Conservative Configuration Patterns

New inputs should choose the tracked fluid explicitly and use time-based slew
limiting so controller behavior is insensitive to output cadence or timestep
changes:

```ini
# Cooled cloud tracking
<frame_tracking>
tracked_fluid = hydro
axes = x1
x1_target = 3.0
target = density
target_min = 5.0
target_max = 200.0
mode = pd
max_abs_boost = 50.0
max_boost_change_mode = per_time
max_boost_change_rate = 5.0
```

```ini
# Mixing-layer tracking
<frame_tracking>
tracked_fluid = hydro
axes = x3
x3_target = 0.0
target = temperature
target_min = 0.015
target_max = 0.08
mode = pd
max_abs_boost = 0.05
max_boost_change_mode = per_time
max_boost_change_rate = 0.02
```

## Examples

- `src/pgen/cloud_crushing.cpp` and
  `inputs/hydro/cloud_crushing_snr.athinput` demonstrate a cooled cloud
  impacted by a Sedov-Taylor boundary inflow.
- `src/pgen/TRML_frame_tracking.cpp` and
  `inputs/hydro/TRML/TRML_frame_tracking.athinput` provide a compact
  turbulent-radiative-mixing-layer example.
- `inputs/hydro/frame_tracking_smoke.athinput` exercises the shared tracker
  with a built-in shock-tube pgen.

The cloud-crushing and TRML examples use moving-frame user boundaries that
sample lab-frame initial/boundary states at `x_grid + FrameDisplacement()` and
store grid-frame velocities as `v_lab - FrameVelocity()`.

## Compatibility Notes

- Frame tracking supports non-relativistic hydro and MHD. It rejects special
  relativistic, general relativistic, and dynamical-GR coordinate systems.
- When both hydro and MHD are present, set `tracked_fluid` explicitly; an
  ambiguous configuration fails at startup.
- `target=entropy` requires an ideal-gas EOS.
- `target=scalarN` requires `0 <= N < nscalars`.
- Compatibility aliases are accepted for an interim deprecation period and
  generate rank-0 warnings. New inputs should use the parameter names listed
  above.
- Arbitrary user boundaries and source terms cannot be made frame-aware
  automatically. Only examples that explicitly apply the moving-frame
  coordinate and velocity transformations are validated.
- AMR changes preserve the tracker object and refresh its pack pointer after
  mesh-block redistribution.

## Related Pages

- [Frame Tracking Production Plan](frame_tracking_next_steps.md)
- [Sedov Cloud Crushing With Frame Tracking](../examples/cloud_crushing_snr.md)
- [TRML Frame Tracking Example](../examples/trml_frame_tracking.md)
