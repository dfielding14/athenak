# Initial Fourier Perturbations

`InitialPerturbations` applies one-time Fourier perturbations to the initial
state after the problem generator runs. The hook is skipped on restarts, so a
checkpointed run does not replay the perturbation.

The feature is enabled with an `<initial_perturbations>` block. It can perturb
density, velocity, and face-centered magnetic fields in non-relativistic
single-fluid hydro or MHD runs. Magnetic perturbations are generated through an
edge-centered vector potential and added as a discrete curl, which preserves
the constrained-transport divergence to roundoff.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `variables` | empty | Comma- or space-separated list containing `density`, `velocity`, and/or `magnetic`. |
| `perturb_density`, `perturb_velocity`, `perturb_magnetic` | `false` | Boolean alternatives to `variables`. A positive RMS also enables the corresponding field. |
| `density_rms` / `rho_rms` | `0.0` | RMS density perturbation. |
| `density_fractional` | `true` | Interpret `density_rms` as a fractional perturbation when true. |
| `velocity_rms` / `v_rms` | `0.0` | RMS velocity perturbation. |
| `magnetic_rms` / `b_rms` | `0.0` | RMS magnetic perturbation. Requires MHD. |
| `spectral_slope` / `expo` | `5/3` | Power-law slope for Fourier amplitudes. |
| `nlow`, `nhigh` | `1`, `3` | Inclusive integer wavenumber shell. `kmin` and `kmax` are accepted aliases. |
| `min_k*`, `max_k*` | `0`, `nhigh` | Optional per-axis integer mode bounds. Inactive dimensions are forced to zero. |
| `f_solenoidal` | `1.0` | Velocity-mode blend: `1` is solenoidal, `0` is compressive. |
| `rseed` | `-1` | RNG seed. Non-negative values give reproducible mode amplitudes. |
| `localization` | `none` | `include` applies a Gaussian window; `exclude` applies one minus that window. |
| `x1/x2/x3_center` | domain center | Gaussian center for localization. |
| `x1/x2/x3_scale` | `-1.0` | Gaussian scale length per coordinate. Positive values activate that coordinate. |
| `remove_density_mean`, `remove_velocity_mean` | `false` | Remove the volume-weighted mean before RMS normalization. |

Example:

```ini
<initial_perturbations>
variables       = density, velocity, magnetic
density_rms     = 1.0e-2
velocity_rms    = 2.0e-2
magnetic_rms    = 1.0e-2
spectral_slope  = 1.0
nlow            = 1
nhigh           = 4
f_solenoidal    = 0.75
rseed           = 314159
localization    = include
x1_center       = 0.0
x2_center       = 0.0
x3_center       = 0.0
x1_scale        = 1.0
x2_scale        = 0.5
x3_scale        = 0.5
```
