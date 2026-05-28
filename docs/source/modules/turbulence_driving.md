# Turbulence Driving

## Scope

`TurbulenceDriver` applies a stochastic acceleration field to hydro, MHD, or
ion-neutral states.  The driver constructs a band-limited Fourier field,
evolves its amplitudes with an Ornstein-Uhlenbeck (OU) process, removes the
net accelerated momentum, and normalizes the remaining field before applying
it during the integrator stages.

This implementation adds three controls needed for multiphase and refined
calculations:

- constant volume-averaged energy injection or constant RMS acceleration;
- periodic repetition of one driving realization on spatial tiles;
- Gaussian inclusion or exclusion envelopes around a chosen region.

The same input controls operate on uniform meshes, static mesh refinement, and
adaptive mesh refinement (AMR).

## Forcing Construction

For selected wavevectors \(\boldsymbol{k}\), the raw acceleration is assembled
from random Fourier amplitudes and cached cell-centered basis functions:

```{math}
\boldsymbol{a}_0(\boldsymbol{x},t)
 = \sum_{\boldsymbol{k}} \mathrm{Re}\left[
   \boldsymbol{A}_{\boldsymbol{k}}(t)
   \exp(i\boldsymbol{k}\cdot\boldsymbol{x})\right].
```

The OU update cadence is `dt_turb_update`.  For a non-zero correlation time,

```{math}
\boldsymbol{A}^{n+1}
 = f_{\rm corr}\boldsymbol{A}^{n}
 + g_{\rm corr}\boldsymbol{\xi}^{n}, \qquad
f_{\rm corr}=\exp(-\Delta t_{\rm turb}/t_{\rm corr}),\quad
g_{\rm corr}=\sqrt{1-f_{\rm corr}^{2}},
```

where \(\boldsymbol{\xi}^{n}\) is a new seeded random realization.  When
`tcorr <= 1e-6`, the driver uses white-noise forcing.

`sol_fraction` blends the solenoidal and compressive Fourier projections:
`1.0` is purely solenoidal and `0.0` is purely compressive.

## Tiled Driving

With `tile_driving = true`, the Fourier realization is evaluated using local
tile coordinates rather than global coordinates.  In coordinate \(q\),

```{math}
q_{\rm tile} = q - q_{\min}
 - L_{\rm tile}\left\lfloor\frac{q-q_{\min}}{L_{\rm tile}}\right\rfloor,
\qquad
L_{\rm tile} = \frac{L}{N_{\rm tile}}.
```

This repeats one statistically identical driving pattern on
`tile_nx * tile_ny * tile_nz` subdomains.  Each tile count must evenly divide
the root-grid cell count in that direction.  Counts in absent dimensions are
automatically reduced to one.

## Spatial Localization

The optional envelope is constructed from each positive scale height:

```{math}
G(\boldsymbol{x}) =
\exp\left[-\frac{1}{2}
\sum_{q\in\{x_1,x_2,x_3\}}
\left(\frac{q-q_0}{H_q}\right)^2\right].
```

Only directions with `*_turb_scale_height > 0` contribute to the sum.
`localization = include` applies \(W=G\), concentrating forcing near the
specified centre.  `localization = exclude` applies \(W=1-G\), suppressing
forcing in that region while retaining it outside.  For backward compatibility,
a positive scale height without an explicit `localization` value selects the
include form.

After applying the envelope, the code subtracts the mass-weighted mean field:

```{math}
\boldsymbol{a} \leftarrow W\boldsymbol{a}_0
- \frac{\int \rho W\boldsymbol{a}_0\,dV}{\int \rho\,dV}.
```

## Normalization Modes

All normalization integrals use physical cell volumes.  This is essential on
AMR meshes: refined blocks represent smaller volumes and must not receive
extra statistical weight merely because they contain more cells.

### Constant Energy Injection

With `constant_edot = true`, `dedt` is the target volume-averaged energy
injection rate.  If the final acceleration is \(s\boldsymbol{a}\), the driver
computes

```{math}
m_0 = \frac{1}{V}\int
\frac{1}{2}\rho |\boldsymbol{a}|^2\Delta t\,dV,\qquad
m_1 = \frac{1}{V}\int
\boldsymbol{p}\cdot\boldsymbol{a}\,dV
```

and applies the non-negative solution of

```{math}
m_0s^2+m_1s=\dot{e}_{\rm target}.
```

### Constant RMS Acceleration

With `constant_edot = false`, `accel_rms` is imposed directly:

```{math}
a_{\rm rms}^2 = \frac{1}{V}\int |\boldsymbol{a}|^2\,dV,\qquad
s = \frac{a_{\rm rms,target}}{a_{\rm rms}}.
```

If `accel_rms` is omitted, the code retains donor-input compatibility by using
`sqrt(2*dedt)` and emits a warning.  New inputs should specify `accel_rms`
explicitly.

## AMR Operation

The forcing and basis arrays are allocated with the same MeshBlock capacity
used by the evolved state arrays (`nmb_maxperrank`).  Before each forcing
update, `EnsureBasisSize` observes changes in the active MeshBlocks and AMR
creation/deletion counters.  When refinement changes the layout it rebuilds
cell-centered Fourier bases from the current MeshBlock geometry, including
tile-local coordinates, before applying localization and volume-weighted
normalization.

Consequently:

- newly refined cells receive a force evaluated at their own physical centres;
- Gaussian envelopes and tile boundaries are defined in physical coordinates;
- energy or acceleration normalization is independent of coarse/fine cell
  multiplicity;
- AMR input files must still set a sufficient `max_nmb_per_rank` capacity.

The AMR rebuild reconstructs the instantaneous forcing on the new block
layout; it is not a claim of bitwise identity with an alternative run that
never refined.

## Input Reference

The following keys belong in `<turb_driving>`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `turb_flag` | `2` | `1` drives for `tdriv_duration`; `2` drives continuously. |
| `tdriv_start` | `0.0` | Time at which driving begins. |
| `tdriv_duration` | `tcorr` | Duration when `turb_flag = 1`. |
| `tcorr` | `0.0` | OU correlation time. |
| `dt_turb_update` | `0.01` | OU refresh cadence; must be positive. |
| `constant_edot` | `true` | Select constant `dedt` or constant `accel_rms`. |
| `dedt` | `0.0` | Target energy injection rate for constant-Edot mode. |
| `accel_rms` | `sqrt(2*dedt)` | Target RMS acceleration when `constant_edot = false`. |
| `nlow`, `nhigh` | `1`, `3` | Inclusive driven mode-radius bounds. |
| `npeak` / `kpeak` | `kpeak = 4*pi` | Spectral peak for the parabolic spectrum; `npeak` is in tile-local x1 mode units, while `kpeak` is an explicit wavenumber. |
| `spect_form` | `1` | `1` parabolic spectrum; `2` power law. |
| `expo` | `5/3` | Isotropic power-law exponent. |
| `min_kx/y/z`, `max_kx/y/z` | `0`, `nhigh` | Optional directional mode bounds. |
| `sol_fraction` | `1.0` | Solenoidal fraction; range `[0,1]`. |
| `rseed` | `-1` | Non-negative values select reproducible random sequences. |
| `tile_driving` | `false` | Enable tiled evaluation of the same realization. |
| `tile_factor` | `1` | Common tile count fallback. |
| `tile_nx/y/z` | `tile_factor` | Per-direction tile counts. |
| `localization` | `none` | `none`, `include`, or `exclude`. |
| `x/y/z_turb_scale_height` | `-1.0` | Positive Gaussian scale heights. |
| `x/y/z_turb_center` | `0.0` | Gaussian centre coordinates. |

## Validation Cases

The repository provides three short runnable inputs under
`inputs/hydro/turb_driving/`:

| Input | Purpose |
| --- | --- |
| `constant_edot_fixed_grid.athinput` | Uniform-grid reference using fixed `dedt`. |
| `tiled_localized_include.athinput` | Tiled forcing with fixed `accel_rms` and a central inclusion envelope. |
| `amr_localized_exclude.athinput` | Adaptive refinement with a central exclusion envelope and fixed `dedt`. |

See {doc}`../examples/turbulence_driving` for commands and simulation-derived
projections of these configurations.
