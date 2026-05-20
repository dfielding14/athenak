# TRML Frame Tracking Example

This example uses `src/pgen/TRML_frame_tracking.cpp` and
`inputs/hydro/TRML/TRML_frame_tracking.athinput` to initialize a
pressure-balanced hot/cold mixing layer and keep the interface near a target
`x3` position with the shared frame tracker.

The pgen is intentionally compact. It demonstrates the moving-frame contract
without pulling in the older production TRML branch history:

- `rho_0`, `pgas_0`, and `density_contrast` define a pressure-balanced hot
  phase and cold phase.
- `velocity` sets the shear speed across the interface.
- `drift_velocity_x3` gives the layer a lab-frame drift so the tracker has a
  motion to cancel.
- `cooling_enabled = true` enables a local top-hat cooling source and a
  problem-specific cooling timestep limit.
- Passive scalar 0 stores the cold fraction.
- User `x3` boundaries sample the lab-frame layer at
  `x_grid + FrameDisplacement()` and set grid-frame velocities with
  `v_grid = v_lab - FrameVelocity()`.

## Frame Tracking Block

The shipped input tracks intermediate-temperature gas around the interface:

```ini
<frame_tracking>
axes       = x3
x3_target  = 0.0
target     = temperature
target_min = 0.015
target_max = 0.08
mode       = pd
```

With the default problem parameters, `T_cold = 0.01` and `T_hot = 1.0`, so this
range selects the cooling/mixing layer rather than the cold or hot bulk.

## Build And Run

Configure and build the custom pgen:

```bash
cmake -S . -B build_trml_frame_tracking -DPROBLEM=TRML_frame_tracking
cmake --build build_trml_frame_tracking -j 4
```

A fast wiring test is:

```bash
rm -rf run_trml_frame_tracking_lowres
mkdir run_trml_frame_tracking_lowres
./build_trml_frame_tracking/src/athena \
  -i inputs/hydro/TRML/TRML_frame_tracking.athinput \
  -d run_trml_frame_tracking_lowres \
  mesh/nx1=16 mesh/nx2=16 mesh/nx3=32 \
  meshblock/nx1=16 meshblock/nx2=16 meshblock/nx3=32 \
  time/nlim=8 time/tlim=0.04 \
  output1/dt=0.02 output2/dt=10 output3/dt=10 \
  frame_tracking/diagnostic_every=1
```

This is not a production TRML configuration. It is a small, PR-friendly example
that verifies the frame-tracking hooks, restart-state path, pgen-local cooling
timestep, and moving-frame boundary transformation.
