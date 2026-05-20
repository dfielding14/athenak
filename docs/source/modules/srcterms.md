# Frame Tracking

`FrameTracker` is a shared post-timestep helper enabled with a
`<frame_tracking>` input block. It samples selected material, estimates its
position and velocity, and applies global Galilean boosts so that material can
remain near a chosen target location in grid coordinates.

The tracker is owned by `MeshBlockPack`, scheduled on the existing
`after_timeintegrator` task list, updates its pack pointer across AMR rebuilds,
and stores its state in restart files through the `<frame_tracking>` block.

Problem generators can query the moving-frame state with:

| API | Meaning |
|-----|---------|
| `FrameVelocity(axis)` | Lab-frame velocity of the grid frame on axis `0`, `1`, or `2`. |
| `FrameDisplacement(axis)` | Lab-frame displacement of the grid-frame origin on axis `0`, `1`, or `2`. |
| `FrameVelocity()` / `FrameDisplacement()` | Three-component `std::array<Real, 3>` versions. |

Common parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `enabled` | `true` | Master switch. |
| `target` | `density` | Cell selector: `density`, `temperature`, `pressure`, `entropy`, `internal_energy`, `scalar`/`scalarN`, `v1`/`v2`/`v3`, or `speed`. |
| `target_min`, `target_max` | unbounded | Inclusive target range used to select tracked cells. |
| `axes` | highest active dimension | Tracking axes: `x1`, `x2`, `x3`, `all`, or comma-separated combinations. |
| `x1_target`, `x2_target`, `x3_target` | domain center | Target centroid location on each tracked axis. |
| `mode` | `pd` | Controller mode: `velocity`, `position`, or `pd`. |
| `apply_every` | `1` | Number of cycles between tracking updates. |
| `start_time` | `0.0` | Simulation time after which updates begin. |
| `diagnostic_every` | `-1` | Cycle cadence for tracker diagnostics; negative disables printing. |
| `tau_avg`, `tau_relax`, `tau_vel` | `1.0` | Filter, position-feedback, and velocity-feedback timescales. |
| `max_abs_boost` | `0.0` | Absolute boost limit; `0.0` disables the cap. |
| `max_boost_change` | `0.0` | Per-update slew limit when positive. |

The frame velocity is the lab velocity of the moving grid frame. The velocity
stored in the fluid is grid-frame velocity, so lab-frame boundary data should
be transformed as `v_grid = v_lab - V_frame`.

For cloud-crushing runs with ISM cooling, `hrate_auto = true` can set
`hydro_srcterms/hrate` from `problem/pressure_over_k`, and
`cooling_timestep_factor` scales the cooling timestep estimate.
