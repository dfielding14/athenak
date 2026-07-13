# `simple_TRML` problem generator and inputs

`simple_TRML` is the minimal production turbulent radiative mixing-layer
problem generator. It initializes the original smooth hot/cold shear layer,
applies radiative cooling, optionally ramps the cooling on, initializes a
cold-material scalar, and can supply a fixed hot reservoir at outer x3.

The earlier experimental implementation is preserved as
`src/pgen/simple_TRML_extended.cpp`. It retains adaptive-frame-aware boundaries
and the experimental lower-x3 boundary choices, but those features are not part
of the minimal `simple_TRML` pgen.

## Building

Build the minimal pgen with:

```bash
cmake -S . -B build-trml -DPROBLEM=simple_TRML
cmake --build build-trml --parallel 16
```

Build the preserved extended pgen with:

```bash
cmake -S . -B build-trml-extended -DPROBLEM=simple_TRML_extended
cmake --build build-trml-extended --parallel 16
```

## Physical setup

The pgen initializes a pressure-balanced shear layer:

```text
rho(z) = smooth tanh from rho_cold at low x3 to rho_hot at high x3
P(z)   = pres
vx(z)  = smooth tanh from -velocity/2 to +velocity/2
vy(z)  = 0
vz(z)  = initial_vz + configured pgen-local perturbation
```

The standard normalization in the supplied inputs is:

```text
rho_hot = 1
pres = 1
gamma = 5/3
chi = rho_cold/rho_hot
T_hot = pres/rho_hot
T_cold = pres/rho_cold
```

For the xi comparison inputs:

```text
chi = 10^1.75 = 56.23413251903491
Mach_rel = 0.5
velocity = v_rel = 0.6454972243679028
t_shear = (x1max - x1min)/velocity = 1.5491933384829668
xi = t_shear/t_cool_min
```

## Initial-condition controls

The primary controls are:

```ini
<problem>
rho_hot = 1.0
rho_cold = 56.23413251903491
pres = 1.0
velocity = 0.6454972243679028
phase_sharpness = 96
initial_vz = 0.0
```

`initial_vz` is literal. The supplied value is stored as the initial vertical
velocity everywhere, before the optional perturbation is added. There is no
frame sign conversion:

```text
initial_vz = -0.035  ->  initial vz = -0.035
initial_vz = +0.035  ->  initial vz = +0.035
```

The pgen-local sinusoidal vertical perturbation is controlled by:

```ini
max_init_perturb_log2freq = ...
min_init_perturb_log2freq = ...
init_perturb_sharpness = ...
init_perturb_vel_frac = ...
```

Set `init_perturb_vel_frac = 0.0` when using the shared
`<initial_perturbations>` module so that the velocity is perturbed only once.

When passive scalars are enabled, every scalar is initialized as conserved
cold-material density:

```text
scalar = rho*cold_fraction
```

The cold fraction varies smoothly from one on the cold side to zero on the hot
side using the same tanh profile as the gas layer.

## X3 boundaries

### Inner x3

The minimal pgen does not implement an inner-x3 user boundary. Select an
AthenaK native boundary in the mesh block, for example:

```ini
ix3_bc = outflow
```

or use AthenaK's native reflecting boundary when desired. Setting
`ix3_bc = user` is rejected at startup because the pgen would otherwise leave
the inner ghost zones unfilled.

### Outer x3

With:

```ini
ox3_bc = user
```

the pgen fills the outer ghost zones with the hot reservoir:

```text
rho = rho_hot
P = pres
vx = +velocity/2
vy = vy in the adjacent active cell
vz = vz in the adjacent active cell
cold fraction = 0
```

Thus rho, pressure, and shear velocity are fixed while `vy` and `vz` are
zero-gradient primitive velocities. Total energy is rebuilt from the imposed
pressure and these velocities.

Set `ox3_bc = outflow` to use AthenaK's native outer outflow instead.

## Cooling

The cooling controls are:

```ini
<problem>
t_cool_min = ...
T_cutoff_over_T_cold = ...
T_ci_over_T_cold = ...
T_ih_over_T_cold = ...
beta = ...
```

The exact cooling update is clamped so that cooling cannot reduce the gas below
`T_cold`. Cooling is disabled at or below `T_cold` and above the configured
upper-temperature cutoff.

The optional startup ramp is:

```ini
cooling_ramp_time = -1.0
cooling_ramp_min_factor = 0.0
```

`cooling_ramp_time <= 0` disables the ramp. Otherwise the cooling rate is
multiplied by:

```text
s = clamp(time/cooling_ramp_time, 0, 1)
ramp = cooling_ramp_min_factor
       + (1 - cooling_ramp_min_factor)*(3*s^2 - 2*s^3)
```

The factor starts at `cooling_ramp_min_factor` and reaches exactly one at
`cooling_ramp_time`.

## History diagnostics

The minimal pgen writes 25 history variables:

- cooling removal rate: `cooling_rate`;
- grid-frame boundary estimates: `M_flux_top`, `M_flux_bot`, `E_flux_top`,
  and `E_flux_bot`;
- phase-binned velocity, kinetic-energy, momentum, and volume diagnostics;
- interface/shear vertical-position diagnostics;
- the instantaneous cooling multiplier, `cool_ramp`.

The boundary flux histories are estimates made from the boundary-adjacent
active cells. They are not the Riemann solver's actual face fluxes.

## Supplied minimal inputs

The three xi comparison inputs use 128 x 128 x 256 cells, eight MeshBlocks,
native lower outflow, the user upper hot reservoir, and a one-shear-time cooling
ramp:

```text
TRML_xi1_M0p5_chi10p1p75_outflowz_coolramp1ts_128x128x256.athinput
TRML_xi10_M0p5_chi10p1p75_outflowz_coolramp1ts_128x128x256.athinput
TRML_xi100_M0p5_chi10p1p75_outflowz_coolramp1ts_128x128x256.athinput
```

`TRML_with_Tracers_and_Tracking.athinput` remains the small tracer-enabled
example despite its historical filename. It now uses the minimal pgen, native
lower outflow, the user upper reservoir, and no adaptive frame tracking.

The `TRML_frame_tracking*.athinput` files belong to the separate
`TRML_frame_tracking.cpp` pgen, not to `simple_TRML`.
