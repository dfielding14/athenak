# `simple_TRML` pgen and TRML input guide

This directory contains inputs for the `simple_TRML` AthenaK problem generator.
The current `simple_TRML` pgen is intentionally clean: it initializes the
original smooth turbulent radiative mixing layer and does not contain the
experimental 1D-front/table initialization path.

## Build target

Configure this pgen with:

```bash
./configure.sh
```

or directly with:

```bash
cmake -S . -B build-frontier-hip-cpe25.09-cce20-rocm6.4.2 -DPROBLEM=simple_TRML
cmake --build build-frontier-hip-cpe25.09-cce20-rocm6.4.2 --parallel 16
```

## Physical setup

The pgen initializes a pressure-balanced hot/cold shear layer.

The standard normalization used by the current inputs is:

```text
rho_hot = 1
P0 = 1
gamma = 5/3
chi = rho_cold/rho_hot
T_hot = P0/rho_hot
T_cold = P0/rho_cold
```

For the current comparison matrix:

```text
chi = 10^1.75 = 56.23413251903491
Mach_rel = 0.5
cs_hot = sqrt(gamma*P0/rho_hot) = 1.2909944487358056
v_rel = Mach_rel*cs_hot = 0.6454972243679028
t_shear = (x1max - x1min)/v_rel = 1.5491933384829668
```

The cooling parameter is:

```text
xi = t_shear/t_cool_min
```

so the current xi sweep uses:

```text
xi = 1:   t_cool_min = 1.5491933384829668
xi = 10:  t_cool_min = 0.15491933384829668
xi = 100: t_cool_min = 0.015491933384829668
```

## Initial conditions

The active-domain IC is the original smooth TRML setup:

- density is a tanh transition from `rho_cold` at low x3 to `rho_hot` at high x3;
- pressure is uniform at `pres`;
- `v_x` is a tanh shear transition from `-0.5*velocity` to `+0.5*velocity`;
- the original lab-frame `v_z` is zero;
- `scalar0` stores conserved cold-material density, `rho*cold_fraction`.

The pgen-local sinusoidal vertical perturbation remains controlled by:

```ini
max_init_perturb_log2freq = ...
min_init_perturb_log2freq = ...
init_perturb_sharpness = ...
init_perturb_vel_frac = ...
```

Most current inputs set `init_perturb_vel_frac = 0.0` and use the shared
`<initial_perturbations>` block instead.

The following old 1D-front options have been removed from this pgen:

```ini
front_j = ...
front_step_x3 = ...
initial_profile = front_step
```

`initial_profile = tanh` and `initial_profile = smooth_tanh` are accepted only
as compatibility aliases for the current smooth IC.

## Velocity-frame convention

All frame-related options use:

```text
v_grid = v_lab - V_frame
v_lab  = v_grid + V_frame
```

The stored AthenaK velocities are grid-frame velocities.

### Constant fixed frame

The pgen-level fixed frame is controlled by:

```ini
<problem>
fixed_frame_velocity_x3 = 0.0
```

This does not require a `<frame_tracking>` block. It only changes how the pgen
and user boundary fills write velocities.

For example, with `fixed_frame_velocity_x3 = +0.1`, the original lab-frame
state `v_z_lab = 0` is stored as:

```text
v_z_grid = -0.1
```

If adaptive frame tracking is also enabled, the effective x3 frame velocity is:

```text
V_frame = fixed_frame_velocity_x3 + adaptive_frame_velocity_x3
```

## Inner-x3 user boundary options

The user lower/inner x3 boundary is selected with:

```ini
<problem>
lower_x3_user_bc = reflect
```

Allowed values are:

```text
reflect
outflow
mass_balance
```

These only apply when the mesh has:

```ini
ix3_bc = user
```

If `ix3_bc = outflow`, AthenaK's normal outflow boundary is used instead.

### `reflect`

`reflect` fills the lower ghost zones with cold gas at pressure `pres` and
reflects the normal velocity in the grid frame. The wall is therefore stationary
in the computational frame, not necessarily in the original lab frame.

Tangential behavior:

- `zero_gradient_vx = true` copies `v_x` from the lower active layer;
- `zero_gradient_vx = false` uses the cold-side reservoir value transformed to
  the grid frame.

### `outflow`

`outflow` copies the adjacent lower active conserved state into the lower ghost
zones.

Use this for a user-managed lower outflow while keeping the pgen's user upper
reservoir active. For true outflow on both z faces, set:

```ini
ix3_bc = outflow
ox3_bc = outflow
```

### `mass_balance`

`mass_balance` fills the lower ghost zones with cold gas and imposes a controlled
normal velocity.

Controls:

```ini
lower_x3_mass_balance_vz = 0.0
lower_x3_mass_balance_velocity_frame = grid
lower_x3_mass_balance_pressure_mode = fixed
lower_x3_mass_balance_tangential_mode = zero_gradient
lower_x3_mass_balance_max_abs_vz = -1.0
```

`lower_x3_mass_balance_vz` is the fallback imposed lower velocity. It is used
unless an enabled frame tracker has a valid live top-Mdot average, in which case
the lower velocity is set from:

```text
v_bot = Mdot_top/(rho_cold*A)
```

where `Mdot_top = integral rho*v3*dA` at the upper x3 face and positive means
outward/+x3.

`lower_x3_mass_balance_velocity_frame` controls how the imposed velocity is
interpreted:

- `grid`: the value is already a stored/grid-frame velocity;
- `lab`: the value is a lab-frame velocity and is stored as
  `v_grid = v_lab - V_frame`.

Pressure modes:

- `fixed`: use `P_bot = P0`;
- `total_pressure`: use the ram-pressure-supported value
  `P_bot = P0 + rho_hot*(v_hot_lab - v_cold_lab)^2`.

Tangential modes:

- `zero_gradient`: copy `v_x` and `v_y` from the lower active layer, then apply
  them at `rho_cold`;
- `reservoir`: use the fixed cold-side tangential reservoir velocity transformed
  to the grid frame.

## Outer-x3 user reservoir

When:

```ini
ox3_bc = user
```

the pgen fills the upper ghost zones with hot gas at `rho_hot`, pressure `pres`,
and cold fraction zero. It enforces the hot-side shear velocity and uses
zero-gradient primitive velocities for `v_y` and `v_z`:

```text
rho = rho_hot
P = pres
v_x,grid = +0.5*v_rel - V_frame,x
v_y,grid = v_y,active
v_z,grid = v_z,active
```

Set `ox3_bc = outflow` to use normal AthenaK outflow instead of this hot
reservoir.

## Cooling options

The base cooling controls are:

```ini
t_cool_min = ...
T_cutoff_over_T_cold = ...
T_ci_over_T_cold = ...
T_ih_over_T_cold = ...
beta = ...
```

Cooling is clamped so that it does not cool below `T_cold`. Cooling is also
disabled below `T_cold` and above `T_cutoff`.

### Cooling ramp

The optional startup ramp is:

```ini
cooling_ramp_time = -1.0
cooling_ramp_min_factor = 0.0
```

`cooling_ramp_time <= 0` disables the ramp. Otherwise the cooling rate is
multiplied by:

```text
s = clamp(time/cooling_ramp_time, 0, 1)
smoothstep = 3*s^2 - 2*s^3
cooling_ramp_factor = cooling_ramp_min_factor
                      + (1 - cooling_ramp_min_factor)*smoothstep
```

For the xi comparison matrix, use:

```ini
cooling_ramp_time = 1.5491933384829668
cooling_ramp_min_factor = 0.0
```

## Adaptive frame tracking

Leave frame tracking disabled by omitting the `<frame_tracking>` block or setting:

```ini
<frame_tracking>
enabled = false
```

The advanced top-Mdot finite-window mode remains available:

```ini
mode = top_mdot_window
mdot_velocity_sign = flipped
mdot_velocity_factor = 1.0
mdot_density_mode = geometric_mean
mdot_start_time = 0.5
mdot_window_time = 1.5491933384829668
```

Sign convention:

```text
Mdot_top = integral rho*v3*dA at upper x3
positive Mdot_top = outward/+x3
v_frame_target = sign * factor * <Mdot_top>_window/(rho_frame*A)
```

`mdot_velocity_sign = flipped` means `sign = -1`.

## History diagnostics

The pgen writes 29 history variables. Important columns include:

- `cooling_rate`;
- `M_flux_top`, `M_flux_bot`;
- `E_flux_top`, `E_flux_bot`;
- `Mtop_lab`, `Mbot_lab`;
- `Etop_lab`, `Ebot_lab`;
- `zavg_i`, `zmin_intermediate`, `zmax_intermediate`;
- `cool_ramp`, the cooling-ramp factor.

The grid-frame fluxes use stored velocities. The `*_lab` fluxes add the current
frame velocity back before forming the flux.

## Current input files

`TRML_with_Tracers_and_Tracking.athinput`

: Clean small tracer-enabled example using user x3 reservoirs. Frame tracking is
  present but disabled by default.

`TRML_xi1_M0p5_chi10p1p75_outflowz_coolramp1ts_128x128x256.athinput`

: xi=1 outflow-z production comparison input.

`TRML_xi10_M0p5_chi10p1p75_outflowz_coolramp1ts_128x128x256.athinput`

: xi=10 outflow-z production comparison input.

`TRML_xi100_M0p5_chi10p1p75_outflowz_coolramp1ts_128x128x256.athinput`

: xi=100 outflow-z production comparison input.

The xi inputs use:

```text
128 x 128 x 256 cells
64 x 64 x 128 meshblocks
8 total meshblocks
ix3_bc = outflow
ox3_bc = outflow
cooling_ramp_time = 1 t_shear
tlim = 31 t_shear
```

## Plotting workflow

The current Frontier scratch/project workflow uses:

```text
/lustre/orion/ast207/proj-shared/dfielding/TRML/simple
```

with plotting helpers:

```text
plot_history.py
plot_profiles.py
plot_slice_panels.py
```

After a run finishes, generate:

- history plots;
- top/bottom mass-flux plots;
- cooling-rate plots;
- vertical profiles at selected times;
- x1/x2/x3 slice panels at evenly spaced times.
