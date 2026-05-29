# Example: Driven Turbulence

AthenaK's turbulence driver evolves a modal Ornstein-Uhlenbeck forcing state
and applies its rendered acceleration to hydro, MHD, or ion-neutral fluid
states. The canonical input interface, AMR and restart invariants, and
localization options are documented in
[Turbulence Driving](../modules/turbulence_driving.md).

## Basic Run

The legacy baseline deck remains available with the current parameter
interface:

```bash
cmake -S . -B build
cmake --build build -j
./build/src/athena -i inputs/hydro/turb.athinput
```

For validated examples that exercise normalization, spatial tiling,
localization, AMR, and restart behavior, use the maintained example page:

- [Localized Turbulence Driving](turbulence_driving.md)

## Canonical Configuration

```ini
<turb_driving>
turb_flag     = 2
tcorr         = 0.25
dt_update     = 0.01
normalization = edot
dedt          = 0.1
nlow          = 1
nhigh         = 3
npeak         = 2
spectrum      = parabolic
sol_fraction  = 1.0
rseed         = 2718
```

Choose `normalization = accel_rms` with `accel_rms = <value>` instead of
`dedt` when the applied acceleration amplitude is the controlled quantity.
For MHD AMR calculations, use conserved prolongation
(`prolong_primitives = false`) so magnetic-field divergence control follows
the validated AMR path.

The provisional names `constant_edot`, `dt_turb_update`, `spect_form`,
`tile_driving`, and `x_turb_scale_height` are not accepted by the current
driver.
