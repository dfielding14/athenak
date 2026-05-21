# Module: Frame Tracking

`FrameTracker` is a shared post-timestep helper for simulations that are easier
to run in a moving Galilean frame. It samples user-selected material, estimates
its position and velocity, and applies a uniform velocity boost to the active
fluid so the selected material remains near a chosen grid-coordinate target.

The tracker is enabled by a `<frame_tracking>` input block. It is owned by
`MeshBlockPack`, scheduled on the existing `after_timeintegrator` task list,
updates its pack pointer after AMR rebuilds, and stores its frame state in
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
| `max_boost_change` | `0.0` | Per-update boost slew limit when `max_boost_change_mode=per_apply`. |
| `max_boost_change_rate` | `-1.0` | Per-time boost slew limit when `max_boost_change_mode=per_time`. |
| `max_boost_change_mode` | `per_apply` | Interpret the boost-change limit per update or per unit time. |
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
| `diagnostic_every` | `-1` | Cycle cadence for rank-0 diagnostic blocks; negative disables printing. |
| `verbose` | `false` | Prints parsed setup at startup. |
| `frame_velocity_x*` | restart only | Stored lab-frame grid velocity. |
| `frame_displacement_x*` | restart only | Stored lab-frame grid-origin displacement. |

Restart output writes `frame_velocity_x1`, `frame_velocity_x2`,
`frame_velocity_x3`, `frame_displacement_x1`, `frame_displacement_x2`, and
`frame_displacement_x3` into the `<frame_tracking>` block.

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

- Exactly one active hydro-like fluid is expected. In a hydro+MHD ion-neutral
  setup, the tracker currently samples and boosts the MHD state.
- `target=entropy` requires an ideal-gas EOS.
- `target=scalarN` requires `0 <= N < nscalars`.
- The tracker applies a uniform Galilean boost after the timestep, so it is
  appropriate for non-relativistic hydro/MHD workflows.
- AMR changes preserve the tracker object and refresh its pack pointer after
  mesh-block redistribution.

## Related Pages

- [Sedov Cloud Crushing With Frame Tracking](../examples/cloud_crushing_snr.md)
- [TRML Frame Tracking Example](../examples/trml_frame_tracking.md)
