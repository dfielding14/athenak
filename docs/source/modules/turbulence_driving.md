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

where \(\boldsymbol{\xi}^{n}\) is a newly generated Gaussian modal innovation,
weighted by the spectrum and selected projection policy. With negligible
`tcorr`, each update replaces the previous state with a white-noise realization.
For three-dimensional `solenoidal_compressive` driving, `sol_fraction=1` gives
solenoidal forcing and `sol_fraction=0` gives compressive forcing.

Mode selection, radial spectrum, and force projection are independent.
All policies select the same spherical shell and evaluate `spectrum` using
$|\boldsymbol{k}|$. `spectrum=parabolic` is the default; `power_law` uses `expo`.
The projection policy selects the acceleration components:

- `solenoidal_compressive` blends projections parallel and perpendicular to
  $\boldsymbol{k}$ using `sol_fraction`.
- `mks24_random_unprojected` uses independent random amplitudes in all three
  components and ignores `sol_fraction`.
- `mks24_alfvenic_perpendicular` projects the acceleration perpendicular to
  both `mean_field_direction` and $\boldsymbol{k}$. It requires
  `sol_fraction=1` and retains variation along the reference field. When
  $\boldsymbol{k}$ is parallel to that direction, both transverse
  polarizations are allowed. The direction is fixed rather than taken from
  the instantaneous local magnetic field.

The default mode set is the full signed spherical shell. The driver applies
the requested exclusions, then stores each real-field conjugate pair only
once. For example, this excludes all modes with zero parallel wavenumber
relative to a diagonal reference field, together with two individual mode
pairs:

```ini
<turb_driving>
nlow = 1
nhigh = 3
mean_field_direction = 1 1 1
exclude_k_parallel_zero = true
exclude_modes = 1 1 0; 1 -1 1
projection_policy = mks24_alfvenic_perpendicular
sol_fraction = 1.0
spectrum = parabolic
```

`mean_field_direction` is a nonzero physical Cartesian vector, normalized
internally; its default is `0 0 1`. The zero-parallel filter removes modes
with $\boldsymbol{k}\cdot\boldsymbol{B}_0=0$, using
$k_i=2\pi n_i/L_{\mathrm{tile},i}$. Thus a reference field along x removes
$k_x=0$, while a field along $(1,1,1)$ removes $k_x+k_y+k_z=0$.
The latter is equivalent to $n_x+n_y+n_z=0$ only for cubic tiles.
This filter controls spatial variation independently of force polarization.

`exclude_modes` contains semicolon-separated integer triplets in tile-local
Fourier mode indices. Excluding $(n_x,n_y,n_z)$ also excludes its negative.
An empty list, written as `none`, excludes no individual modes.

Optional `min_kx/y/z` and `max_kx/y/z` remain signed integer search bounds.
Defaults cover the full shell: $-n_{\mathrm{high}}$ through
$n_{\mathrm{high}}$ for index shells, or the required component bounds for
physical shells and tile lengths. Inactive dimensions use zero. The
physical mode set includes the conjugate partners of the searched modes;
overlapping pairs are evaluated only once. The full signed default increases
angular coverage and changes fixed-seed realizations from the older
positive-octant default.

The direction configures the forcing; it does not initialize or rotate the
magnetic field. The current CGL paper problem generator initializes its mean
field along z, so its input and z-based history checks should retain
`mean_field_direction = 0 0 1` until that initializer is extended.

`physical_k_shell` changes mode selection only. The retained parabolic
weighting uses `kpeak = npeak * 2*pi/L_tile,x1` (or explicit `kpeak`), with
width endpoints `nlow * min(dkx,dky,dkz)` and
`nhigh * max(dkx,dky,dkz)`, where `dk_i = 2*pi/L_tile,i`. These do not
automatically become `nlow * k_shell_unit` and `nhigh * k_shell_unit`
in a physical shell.
Check the intended peak and weighting separately for an elongated box.
Parabolic power is zero beyond the parabola's roots, even when those roots
lie inside the selected shell. On a continuation from an older executable,
stored amplitudes in those modes decay on `tcorr` and receive no further
random increments. Use the same forcing version for
matched simulations.

`mode_amp_real` and `mode_amp_imag` are the authoritative OU state. The
cell-centered `force` field is rendered from those coefficients and the
current MeshBlock geometry.

For Newtonian fluids without radiation, modal updates run before RK evolution;
the force is normalized against the current conserved state and applied once as
a full-`dt` kick at the final fluid stage. The existing boundary exchange and
primitive recovery follow the kick, so no separate pre-RK refresh is needed.
Radiation and relativistic task graphs retain their previous pre-RK scheduling.
OU coefficients remain fixed between `dt_update` boundaries; due updates are
processed at cycle start without applying extra fluid kicks. Driving start and
stop times also use cycle-start time, without shortening a crossing timestep.
The single kick gives first-order splitting for smooth forcing; moving its
position does not establish higher temporal accuracy of the coupled evolution.

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
injection rate. A zero target disables forcing, including when the field opposes
the flow. For a scale factor \(s\), the driver evaluates

```{math}
m_0=\frac{1}{V}\int \frac{1}{2}\rho|\boldsymbol{a}|^2\Delta t\,dV,
\qquad
m_1=\frac{1}{V}\int \boldsymbol{p}\cdot\boldsymbol{a}\,dV,
```

and uses the non-negative solution to

```{math}
m_0s^2+m_1s=\dot e_{\rm target}.
```

For ion-neutral runs, these integrals use the sum of the two fluids' densities
and conserved momenta.

### Conserved Momentum And Energy

For each nonrelativistic fluid, the normalized acceleration applies

```{math}
\boldsymbol{p}'=\boldsymbol{p}+\rho\boldsymbol{a}\Delta t,\qquad
E'=E+\boldsymbol{p}\cdot\boldsymbol{a}\Delta t
       +\tfrac12\rho|\boldsymbol{a}|^2\Delta t^2.
```

The work uses that fluid's own pre-kick conserved momentum and density, so a
pure kick preserves thermal energy up to roundoff even when primitives are
stale. Active CGL retains its thermal energy and conserved anisotropy. Isothermal
fluids receive only the momentum increment; passive CGL also leaves both J/A
thermodynamic slots unchanged. The kick leaves magnetic fields unchanged.
Net-momentum removal acts on the acceleration before this update.

With `record_injected_work=true`, each full kick adds its measured total-energy
change to the cumulative work counter for an active ideal/CGL fluid, or its
kinetic-energy change for isothermal/passive CGL. This diagnostic supports a
single nonrelativistic fluid and does not add work to passive thermodynamics.

### RMS Acceleration

With `normalization=accel_rms`, `accel_rms` is imposed directly:

```{math}
a_{\rm rms}^2=\frac{1}{V}\int|\boldsymbol{a}|^2\,dV,\qquad
s=\frac{a_{\rm rms,target}}{a_{\rm rms}}.
```

This fixes acceleration amplitude, not injected power or the resulting Mach
number. Work depends on the evolving fluid state and can have either sign.

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
- the exact retained integer mode list;
- random-number generator state;
- the real and imaginary coefficient arrays.
- accumulated injected work when `record_injected_work = true`.

It does not checkpoint the rendered cell-centered force array. On restart,
the force is regenerated on the restored mesh from the authoritative modal
state. A checkpoint is rejected if turbulence-driving configuration differs
from the run reading it. The unified driver writes turbulence restart
version 5; older checkpoint layouts are rejected rather than reinterpreted.
Validation includes the normalized reference direction, zero-parallel
filter, and exact retained mode list, so different exclusions with equal
mode counts cannot silently reuse an incompatible modal state.

The final-stage kick does not preserve trajectories from the older pre-RK
binary; record the executable revision with every continuation.

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
| `spectrum` | `parabolic` | Radial `parabolic` or `power_law` spectrum for every projection policy. |
| `expo` | `5/3` | Radial power-law spectrum exponent. |
| `min_kx/y/z`, `max_kx/y/z` | Full signed shell | Optional inclusive signed search bounds; inactive dimensions use zero. |
| `mean_field_direction` | `0 0 1` | Nonzero physical Cartesian reference direction, normalized internally; does not initialize the magnetic field. |
| `exclude_k_parallel_zero` | `false` | Remove every mode with physical `k` perpendicular to the reference direction. |
| `exclude_modes` | `none` | Semicolon-separated integer triplets; remove each listed mode and its conjugate partner. |
| `projection_policy` | `solenoidal_compressive` | `solenoidal_compressive`, `mks24_random_unprojected`, or `mks24_alfvenic_perpendicular`; independent of shell and spectrum. |
| `physical_k_shell` | `false` | Apply `nlow`/`nhigh` to `abs(k)/k_shell_unit`. |
| `k_shell_unit` | `0.0` | Positive reference wavenumber required by `physical_k_shell`. |
| `sol_fraction` | `1.0` | Solenoidal/compressive amplitude blend in `[0,1]`; unused by unprojected driving and required to be 1 for the Alfvénic policy. |
| `rseed` | `-1` | Positive values select reproducible sequences; zero and negative values use seed 1. |
| `record_injected_work` | `false` | Restart measured forcing work for one nonrelativistic fluid; kinetic work for isothermal/passive CGL. |
| `tile_nx`, `tile_ny`, `tile_nz` | `1` | Tile repetitions in each direction. |
| `localization` | `none` | `none`, `include`, or `exclude`. |
| `sigma_x1`, `sigma_x2`, `sigma_x3` | `-1.0` | Positive Gaussian widths in active directions. |
| `center_x1`, `center_x2`, `center_x3` | `0.0` | Gaussian center coordinates. |

## Rejected Prototype Names

This interface has no deprecated aliases. The removed `driving_type`,
`isotropic_power_spectrum`, `exp_prp`, and `exp_prl` keys are rejected. Select
the desired spectrum and projection explicitly when migrating an input.
The four legacy planar example decks now use 3D perpendicular solenoidal
forcing; this deliberately changes their old planar spatial forcing.
In particular, inputs containing
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
