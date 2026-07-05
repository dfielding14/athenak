# CR Pusher Accuracy Test Plan

This plan validates the updated CR particle pusher, especially the
per-particle gyro subcycling fast path:

```text
particles/subcycle_per_particle_gyro = true
particles/exchange_gyro_only_substeps = false
```

The production speedup comes from avoiding the old behavior where every
particle paid the globally worst gyro substep count. The tests below are meant
to verify that this change preserves the expected Boris-pusher accuracy in
controlled fields.

Production snapshot A/B tests are intentionally excluded here. These are clean
unit/integration-style pusher tests.

## Test 1: Uniform Magnetic Field

Purpose: isolate the Boris integrator and gyro subcycling with an analytic
reference.

Setup:

```text
B = B0 zhat
E = 0
MHD evolution = off/frozen
boundary = periodic
particles = several masses/species and pitch angles
```

Run matrix:

```text
subcycle_per_particle_gyro = false, gyro_fraction = 0.05   # old global reference
subcycle_per_particle_gyro = true,  gyro_fraction = 0.05
subcycle_per_particle_gyro = true,  gyro_fraction = 0.10
subcycle_per_particle_gyro = true,  gyro_fraction = 0.20
subcycle_per_particle_gyro = true,  gyro_fraction = 0.30
```

Diagnostics:

```text
gyro radius error
gyro period / phase error
parallel velocity conservation
total speed / kinetic energy conservation
magnetic moment conservation
long-time phase drift over many gyroperiods
```

Expected result:

```text
gyro_fraction=0.05 per-particle and global-subcycle should agree closely.
Larger gyro_fraction values should show smooth convergence toward the strict reference.
```

## Test 2: Smooth Divergence-Free Magnetic Field

Purpose: test field interpolation and adiabatic motion in a nonuniform but
controlled magnetic field.

Construct the field from a vector potential so `div B = 0`. One simple example:

```text
A_z = A0 sin(kx) sin(ky)
B_x =  dA_z/dy
B_y = -dA_z/dx
B_z = B0
```

Choose `B0` large enough that the field does not reverse.

Setup:

```text
E = 0
MHD evolution = off/frozen
boundary = periodic
particles = several masses, gyroradii, and pitch angles
```

Strict numerical reference:

```text
subcycle_per_particle_gyro = true
gyro_fraction = 0.01
```

Run matrix:

```text
subcycle_per_particle_gyro = true,  gyro_fraction = 0.02
subcycle_per_particle_gyro = true,  gyro_fraction = 0.05
subcycle_per_particle_gyro = true,  gyro_fraction = 0.10
subcycle_per_particle_gyro = true,  gyro_fraction = 0.20
subcycle_per_particle_gyro = true,  gyro_fraction = 0.30
subcycle_per_particle_gyro = false, gyro_fraction = 0.05
```

Diagnostics:

```text
position error relative to strict reference
velocity error relative to strict reference
magnetic moment drift
total speed / kinetic energy conservation
pitch-angle evolution
species-dependent error scaling
per-particle vs global-subcycle agreement at comparable particle timestep
```

Expected result:

```text
Errors should decrease monotonically as gyro_fraction decreases.
No species should show anomalous nonconvergence.
Per-particle subcycling should match the global path when the same particle-level timestep is used.
```

## Test 3: MeshBlock Boundary / Exchange Stress Test

Purpose: verify correctness when particles cross block boundaries frequently.

Setup:

```text
field = uniform B or smooth div-free B
MHD evolution = off/frozen
MeshBlocks = small enough to force frequent crossings
particles = large parallel velocity plus moderate perpendicular velocity
checks = consistency and motion-bound checks enabled
```

Run matrix:

```text
subcycle_per_particle_gyro = false, gyro_fraction = 0.05
subcycle_per_particle_gyro = true,  gyro_fraction = 0.05
subcycle_per_particle_gyro = true,  gyro_fraction = 0.10
```

Diagnostics:

```text
particle count conservation
duplicate tag detection
missing tag detection
invalid PGID detection
position/velocity continuity through MeshBlock boundaries
agreement with single-block or low-crossing reference
```

Expected result:

```text
No particle loss or duplication.
No invalid ownership.
No trajectory discontinuity at MeshBlock boundaries.
The code should fall back to the exchange-safe global path whenever cell or MeshBlock crossing constraints require intermediate ownership exchange.
```

## Implementation Notes

Suggested input files:

```text
inputs/particles/cr_pusher_uniform_b.athinput
inputs/particles/cr_pusher_smooth_divb0.athinput
inputs/particles/cr_pusher_boundary_crossing.athinput
```

Suggested driver/analyzer scripts:

```text
scripts/run_cr_pusher_accuracy_tests.sh
scripts/analyze_cr_pusher_accuracy.py
```

The analyzer should write tables like:

```text
test  mode  gyro_fraction  species  max_dx  rms_dx  max_dv  rms_dv  dmu/mu  dE/E
```

For the uniform-field test, compare against the analytic helix. For the
smooth-field and boundary tests, compare against the strict numerical reference.

## Acceptance Criteria

Before relaxing production subcycling beyond the current conservative setting:

```text
gyro_fraction=0.05 per-particle must match or improve the old global-subcycle result.
gyro_fraction=0.10 is acceptable only if trajectory and moment errors remain small.
gyro_fraction=0.20 or 0.30 remains experimental unless convergence is very clean.
No boundary-crossing test may lose or duplicate particles.
Energy/speed conservation must remain consistent with Boris expectations.
Magnetic moment drift must converge with decreasing gyro_fraction.
```

## Recommended Work Order

1. Implement the uniform-B test and analytic analyzer.
2. Validate old global path vs new per-particle path at `gyro_fraction=0.05`.
3. Add the smooth divergence-free field test with strict-reference comparison.
4. Add the boundary/exchange stress test.
5. Run the gyro-fraction sweep.
6. Summarize error vs runtime and decide whether production can safely use
   anything looser than `gyro_fraction=0.05`.
