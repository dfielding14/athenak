# Turbulence Driving

## Scope

`TurbulenceDriver` applies a stochastic acceleration source to hydro, MHD, or
ion-neutral fluid states. It supports uniform meshes and adaptive mesh
refinement (AMR), with three controls intended for localized multiphase
calculations:

- normalization by fixed energy injection rate or fixed RMS acceleration;
- periodic repetition of one forcing realization on spatial tiles;
- Gaussian inclusion or exclusion of a specified physical region.

The implementation keeps stochastic state in Fourier coefficient space, not
in cell arrays. Mesh changes therefore alter only the evaluation geometry.

## Modal Evolution

For selected wavevectors \(\boldsymbol{k}\), the unnormalized acceleration is

```{math}
\boldsymbol{a}_0(\boldsymbol{x},t)
 = \sum_{\boldsymbol{k}} \operatorname{Re}\left[
 \boldsymbol{A}_{\boldsymbol{k}}(t)
 \exp(i\boldsymbol{k}\cdot\boldsymbol{x})\right].
```

At cadence `dt_update`, each complex modal coefficient is updated by an
Ornstein-Uhlenbeck (OU) recurrence,

```{math}
\boldsymbol{A}^{n+1}
= f_{\rm corr}\boldsymbol{A}^{n}
+ g_{\rm corr}\boldsymbol{\xi}^{n},\qquad
f_{\rm corr}=\exp(-\Delta t_{\rm update}/t_{\rm corr}),\qquad
g_{\rm corr}=\sqrt{1-f_{\rm corr}^{2}},
```

where \(\boldsymbol{\xi}^{n}\) is a newly generated, projected modal
innovation. With negligible `tcorr`, each update replaces the previous state
with a white-noise realization. `sol_fraction=1` gives solenoidal forcing and
`sol_fraction=0` gives compressive forcing.

`mode_amp_real` and `mode_amp_imag` are the authoritative OU state. The
cell-centered `force` field is rendered from those coefficients and the
current MeshBlock geometry.

## Tiled Evaluation

`tile_nx`, `tile_ny`, and `tile_nz` specify the number of repetitions in each
coordinate direction. A count of one disables tiling in that direction. For
coordinate \(q\),

```{math}
q_{\rm tile} = q - q_{\min}
 - L_{\rm tile}\left\lfloor\frac{q-q_{\min}}{L_{\rm tile}}\right\rfloor,
\qquad L_{\rm tile}=L/N_{\rm tile}.
```

The basis is evaluated at \(q_{\rm tile}\), so every tile receives the same
modal realization. Tile counts must be positive and evenly divide the
root-grid cell count; inactive dimensions require a tile count of one.

## Spatial Localization

When `localization=include` or `localization=exclude`, positive widths among
`sigma_x1`, `sigma_x2`, and `sigma_x3` define

```{math}
G(\boldsymbol{x}) =
\exp\left[-\frac{1}{2}\sum_q
\left(\frac{x_q-\mathrm{center}_{x_q}}{\sigma_{x_q}}\right)^2\right].
```

Directions with non-positive widths are omitted. The applied envelope is
\(W=G\) for `include` and \(W=1-G\) for `exclude`. The driver then removes
the mass-weighted mean acceleration,

```{math}
\boldsymbol{a} \leftarrow W\boldsymbol{a}_0
-\frac{\int \rho W\boldsymbol{a}_0\,dV}{\int \rho\,dV}.
```

`localization=none` rejects positive `sigma_x*` parameters, and localized
modes require at least one positive width.

## Normalization

All integrals include physical cell volumes. On AMR meshes this prevents fine
cells from receiving extra statistical weight solely because they are more
numerous.

### Energy Injection

With `normalization=edot`, `dedt` specifies the target volume-averaged energy
injection rate. For a scale factor \(s\), the driver evaluates

```{math}
m_0=\frac{1}{V}\int \frac{1}{2}\rho|\boldsymbol{a}|^2\Delta t\,dV,
\qquad
m_1=\frac{1}{V}\int \boldsymbol{p}\cdot\boldsymbol{a}\,dV,
```

and uses the non-negative solution to

```{math}
m_0s^2+m_1s=\dot e_{\rm target}.
```

### RMS Acceleration

With `normalization=accel_rms`, `accel_rms` is imposed directly:

```{math}
a_{\rm rms}^2=\frac{1}{V}\int|\boldsymbol{a}|^2\,dV,\qquad
s=\frac{a_{\rm rms,target}}{a_{\rm rms}}.
```

`dedt` and `accel_rms` are mutually exclusive: the selected normalization
mode requires its target and rejects the other target.

## AMR And Restart Invariants

On every refinement topology change, the driver resizes capacity if required
and recomputes only the trigonometric basis on active MeshBlocks. It does not
interpolate, restrict, or regenerate modal OU coefficients. Thus newly
created fine cells evaluate the same instantaneous realization at their
physical centers, and tile/localization boundaries remain physical-coordinate
operations.

Restart output records:

- a versioned configuration record and completed OU update count;
- random-number generator state;
- the real and imaginary coefficient arrays.
- accumulated injected work when `record_injected_work = true`.

It does not checkpoint the rendered cell-centered force array. On restart,
the force is regenerated on the restored mesh from the authoritative modal
state. A checkpoint is rejected if turbulence-driving configuration differs
from the run reading it. Prototype checkpoint layouts are intentionally not
supported.

Real-valued defaults added at runtime are read back from their serialized
parameter representation on the initial run, so an unchanged restart uses
the same exact configuration record as the uninterrupted calculation.

MPI reductions use `MPI_ATHENA_REAL`, so the driver source is valid for both
double and single `Real` builds. Full single-precision executable builds of
the current upstream tree remain blocked by pre-existing non-turbulence
sources; see the example validation notes.

## Input Reference

The following keys belong in `<turb_driving>`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `turb_flag` | `2` | `1` drives for `tdriv_duration`; `2` drives continuously. |
| `tdriv_start` | `0.0` | Time at which driving begins. |
| `tdriv_duration` | `tcorr` | Duration when `turb_flag=1`. |
| `tcorr` | `0.0` | OU correlation time; must be non-negative. |
| `dt_update` | `0.01` | Modal update cadence; must be positive. |
| `normalization` | `edot` | `edot` or `accel_rms`. |
| `dedt` | required for `edot` | Target volume-averaged injection rate. |
| `accel_rms` | required for `accel_rms` | Target volume-weighted RMS acceleration. |
| `nlow`, `nhigh` | `1`, `3` | Inclusive driven mode-radius bounds. |
| `npeak` / `kpeak` | `kpeak=4*pi` | Parabolic spectral peak; `npeak` is tile-local mode number. |
| `spectrum` | `parabolic` | `parabolic` or `power_law`. |
| `expo`, `exp_prp`, `exp_prl` | `5/3`, `5/3`, `0` | Power-law spectrum exponents. |
| `min_kx/y/z`, `max_kx/y/z` | `0`, `nhigh` | Optional directional mode bounds. |
| `driving_type` | `0` | `0` for three-dimensional; `1` for planar driving. |
| `physical_k_shell` | `false` | Apply `nlow`/`nhigh` to `abs(k)/k_shell_unit`. |
| `k_shell_unit` | `0.0` | Positive reference wavenumber required by `physical_k_shell`. |
| `isotropic_power_spectrum` | `false` | Apply `expo` to total `abs(k)` for planar driving. |
| `sol_fraction` | `1.0` | Solenoidal fraction in `[0,1]`. |
| `rseed` | `-1` | Non-negative values select reproducible sequences. |
| `record_injected_work` | `false` | Restart exact accumulated forcing work for a single nonrelativistic ideal/CGL fluid. |
| `tile_nx`, `tile_ny`, `tile_nz` | `1` | Tile repetitions in each direction. |
| `localization` | `none` | `none`, `include`, or `exclude`. |
| `sigma_x1`, `sigma_x2`, `sigma_x3` | `-1.0` | Positive Gaussian widths in active directions. |
| `center_x1`, `center_x2`, `center_x3` | `0.0` | Gaussian center coordinates. |

## Rejected Prototype Names

This interface has no deprecated aliases. In particular, inputs containing
`constant_edot`, `tile_factor`, or `x_turb_scale_height` terminate with an
input error rather than silently selecting a new behavior. The driver also
rejects the provisional names `dt_turb_update`, `spect_form`,
`tile_driving`, `y_turb_scale_height`, `z_turb_scale_height`, and
`x/y/z_turb_center`; use the canonical keys above.

## Included Problems

Three runnable demonstrations are provided under `inputs/hydro/turb_driving/`:

| Input | Purpose |
| --- | --- |
| `constant_edot_fixed_grid.athinput` | Uniform-grid reference with `normalization=edot`. |
| `tiled_localized_include.athinput` | Tiled fixed-RMS forcing with a Gaussian inclusion region. |
| `amr_localized_exclude.athinput` | AMR calculation with a Gaussian exclusion region. |

See {doc}`../examples/turbulence_driving` for commands, measured checks, and
projections generated from these inputs.
