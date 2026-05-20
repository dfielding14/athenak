# Sedov Cloud Crushing With Frame Tracking

This example uses `src/pgen/cloud_crushing.cpp` and
`inputs/hydro/cloud_crushing_snr.athinput` to initialize a cold ISM cloud in
pressure balance with a warm ISM, then drives the left boundary with a
Taylor-von Neumann-Sedov blast wave.

The problem is built around three pieces:

- ISM cooling and heating set the cold and warm thermal equilibria at the
  requested `problem/pressure_over_k` and `hydro/hrate`.
- The cold cloud starts as a smoothed density sphere embedded in the warm phase.
- The `inner_x1` user boundary samples the full xi-dependent TVNS interior
  profile, not just a constant post-shock state.

## Units

The shipped input uses cgs-backed code units:

```ini
<units>
length_cgs  = 3.0856775809623245e18  # 1 pc
time_cgs    = 3.15576e13             # 1 Myr
mass_cgs    = 4.91416144649125e31    # 1 amu cm^-3 * (1 pc)^3
mu          = 1.0
```

With the default cooling/heating choices, the equilibrium phases are

| Quantity | Value |
| --- | ---: |
| Cold temperature | `50.2123 K` |
| Warm temperature | `6296.21 K` |
| Cold number density | `62.9781 cm^-3` |
| Warm number density | `0.502251 cm^-3` |
| Ambient pressure | `4.36599e-13 dyn cm^-2` |
| Ambient mass density | `8.40076e-25 g cm^-3` |
| Cold-cloud mass density | `1.05339e-22 g cm^-3` |

## Sedov Boundary Options

The Sedov boundary is controlled from the `<problem>` block:

| Parameter | Meaning |
| --- | --- |
| `sedov_energy_cgs` | Explosion energy in ergs. |
| `sedov_beta` | Dimensionless Sedov similarity constant. |
| `sedov_start_time` | Simulation time at which the blast evolution begins. |
| `sedov_origin_distance` | Default distance from `x1min` to the explosion origin. |
| `sedov_origin_x1`, `sedov_origin_x2`, `sedov_origin_x3` | Optional explicit lab-frame origin. |
| `sedov_radius_at_start` | Shock radius at `sedov_start_time`; together with `E` and ambient density this fixes the Sedov age. |

At each ghost cell, the code computes

```text
x_lab = x_grid + X_frame
r_lab = |x_lab - x_snr|
xi    = r_lab / R_s(t)
v_grid = v_TVNS,lab(r_lab,t) - V_frame
```

where `X_frame` and `V_frame` come from the shared `FrameTracker` when
`<frame_tracking>` is enabled. This is the key moving-frame contract: the grid
coordinates remain grid-frame coordinates, but user boundary conditions can
recover the lab-frame sample point from the exposed frame displacement.

If the sample point lies outside the blast radius, the boundary fills warm ISM
gas moving at `-V_frame` in grid coordinates. If it lies inside the shock, the
boundary fills the TVNS density, pressure, and radial velocity at that `xi`.

## Frame Tracking

The input tracks the cold cloud by selecting dense gas:

```ini
<frame_tracking>
axes       = x1
x1_target  = 3.0
target     = density
target_min = 5.0
target_max = 200.0
mode       = pd
```

The tracker applies post-timestep Galilean velocity boosts. It also exposes
the current frame velocity and displacement through `FrameVelocity(axis)` and
`FrameDisplacement(axis)`. Restart files store these values in the
`<frame_tracking>` block as `frame_velocity_x*` and `frame_displacement_x*`.

The frame velocity is the lab velocity of the grid frame. The boost applied to
the fluid is the opposite cumulative velocity, so problem generators should
transform lab-frame velocities as

```text
v_grid = v_lab - V_frame
```

Frame tracking does not remap mesh coordinates or move AMR blocks. The
displacement is the bookkeeping needed to evaluate lab-frame physics
consistently while the grid-frame fluid variables are Galilean shifted.

## Build And Run

Configure and build the problem:

```bash
cmake -S . -B build_cloud_crushing -DPROBLEM=cloud_crushing
cmake --build build_cloud_crushing -j 4
```

A fast low-resolution validation run is:

```bash
rm -rf run_cloud_crushing_lowres
mkdir run_cloud_crushing_lowres
./build_cloud_crushing/src/athena \
  -i inputs/hydro/cloud_crushing_snr.athinput \
  -d run_cloud_crushing_lowres \
  mesh/nx1=32 mesh/nx2=16 mesh/nx3=16 \
  meshblock/nx1=16 meshblock/nx2=8 meshblock/nx3=8 \
  time/nlim=160 time/tlim=0.004 time/cfl_number=0.01 \
  output1/dt=2.0e-4 output2/dt=5.0e-4 \
  frame_tracking/diagnostic_every=5 \
  frame_tracking/max_abs_boost=50.0 \
  frame_tracking/max_boost_change=2.0 \
  2>&1 | tee run_cloud_crushing_lowres/athena.log
```

This run is intentionally low resolution. It is a wiring and consistency test,
not a production cloud-crushing calculation.

## Validation Plots

Plot the analytic boundary state in cgs units through `t = 10` code units:

```bash
python3 scripts/plot_sedov_boundary.py \
  --output docs/source/_static/sedov_boundary_cgs.png \
  --csv-output docs/source/_static/sedov_boundary_cgs.csv \
  --t-max 10
```

![Sedov boundary state](../_static/sedov_boundary_cgs.png)

Plot the low-resolution run diagnostics:

```bash
python3 scripts/plot_cloud_crushing_validation.py \
  run_cloud_crushing_lowres \
  --output-prefix docs/source/_static/cloud_crushing_lowres
```

The validation script reads `athena.log` for tracker diagnostics and the
binary VTK `hydro_w` snapshots for an independent cold-cloud centroid check.
The run used for this documentation produced 31 frame-tracking diagnostic
blocks and 5 VTK snapshots. The final diagnostic frame state was
`v_frame=-256.259` and `X_frame=-0.142559` in code units at
`t=0.00151389`; the final VTK cold-cloud centroid was `x1=3.1529` at
`t=0.00153755`.

![Frame-tracking validation](../_static/cloud_crushing_lowres_validation.png)

The midplane density snapshots show the dense cloud and the shocked inflow
entering from the left boundary in the same run.

![Low-resolution midplane density](../_static/cloud_crushing_lowres_midplane.png)
