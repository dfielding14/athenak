# Turbulence Driving

## Scope

`TurbulenceDriver` applies a stochastic acceleration source to hydro, MHD, or
ion-neutral fluid states. It supports uniform meshes and adaptive mesh
refinement (AMR), with three controls intended for localized multiphase
calculations:

- normalization by fixed RMS acceleration or fixed energy injection rate;
- periodic repetition of one forcing realization on spatial tiles;
- Gaussian inclusion or exclusion of a specified physical region.

The implementation keeps stochastic state in Fourier coefficient space, not
in cell arrays. Mesh changes therefore alter only the evaluation geometry.

## Modal Evolution

For selected wavevectors $\boldsymbol{k}$, the unnormalized acceleration is

$$
\boldsymbol{a}_0(\boldsymbol{x},t)
 = \sum_{\boldsymbol{k}} \operatorname{Re}\left[
 \boldsymbol{A}_{\boldsymbol{k}}(t)
 \exp(i\boldsymbol{k}\cdot\boldsymbol{x})\right].
$$

At cadence `dt_update`, each complex modal coefficient is updated by an
Ornstein-Uhlenbeck (OU) recurrence,

$$
\boldsymbol{A}^{n+1}
= f_{\rm corr}\boldsymbol{A}^{n}
+ g_{\rm corr}\boldsymbol{\xi}^{n},\qquad
f_{\rm corr}=\exp(-\Delta t_{\rm update}/t_{\rm corr}),\qquad
g_{\rm corr}=\sqrt{1-f_{\rm corr}^{2}},
$$

where $\boldsymbol{\xi}^{n}$ is a newly generated, projected modal
innovation. With negligible `tcorr`, each update replaces the previous state
with a white-noise realization. `sol_fraction=1` gives solenoidal forcing and
`sol_fraction=0` gives compressive forcing.

`mode_amp_real` and `mode_amp_imag` are the authoritative OU state. The
cell-centered `force` field is rendered from those coefficients and the
current MeshBlock geometry.

The modal update cadence and random sequence are independent of the RK stages.
The actual acceleration is obtained after localization, mass-weighted mean
removal, and normalization. Its RMS is therefore prescribed in `accel_rms`
mode, but its temporal statistics are not those of an unmodified Gaussian OU
field.

## Tiled Evaluation

`tile_nx`, `tile_ny`, and `tile_nz` specify the number of repetitions in each
coordinate direction. A count of one disables tiling in that direction. For
coordinate $q$,

$$
q_{\rm tile} = q - q_{\min}
 - L_{\rm tile}\left\lfloor\frac{q-q_{\min}}{L_{\rm tile}}\right\rfloor,
\qquad L_{\rm tile}=L/N_{\rm tile}.
$$

The basis is evaluated at $q_{\rm tile}$, so every tile receives the same
modal realization. Tile counts must be positive and evenly divide the
root-grid cell count; inactive dimensions require a tile count of one.

## Spatial Localization

When `localization=include` or `localization=exclude`, positive widths among
`sigma_x1`, `sigma_x2`, and `sigma_x3` define

$$
G(\boldsymbol{x}) =
\exp\left[-\frac{1}{2}\sum_q
\left(\frac{x_q-\mathrm{center}_{x_q}}{\sigma_{x_q}}\right)^2\right].
$$

Directions with non-positive widths are omitted. The applied envelope is
$W=G$ for `include` and $W=1-G$ for `exclude`. The driver then removes
the mass-weighted mean acceleration,

$$
\boldsymbol{a} \leftarrow W\boldsymbol{a}_0
-\frac{\int \rho W\boldsymbol{a}_0\,dV}{\int \rho\,dV}.
$$

`localization=none` rejects positive `sigma_x*` parameters, and localized
modes require at least one positive width.

## Normalization

The normalization modes are shared by Newtonian and special-relativistic
driving. The full-timestep kick and exact-work guarantees below apply to
Newtonian hydro and MHD, including ion-neutral fluids. Single-fluid
special-relativistic driving retains its legacy stage updates; ion-neutral
turbulence driving requires Newtonian fluids.
General-relativistic and dynamical-spacetime turbulence driving are unsupported.

All integrals include physical cell volumes. On AMR meshes this prevents fine
cells from receiving extra statistical weight solely because they are more
numerous. For Newtonian fluids, the force is rendered and normalized from the
current conserved state immediately before its full-timestep kick.
Non-finite normalization integrals are rejected before scaling, with a
separate diagnostic from a zero forcing field.

`normalization=accel_rms` is the default. The strength must still be supplied
explicitly as `accel_rms`; there is no default forcing amplitude. Existing
power-controlled inputs must specify `normalization=edot` alongside `dedt`.

### RMS Acceleration

With `normalization=accel_rms`, the volume-weighted RMS of the applied
acceleration is imposed directly:

$$
a_{\rm rms}^2=\frac{1}{V}\int|\boldsymbol{a}|^2\,dV,\qquad
s=\frac{a_{\rm rms,target}}{a_{\rm rms}}.
$$

The resulting work depends on the fluid state and can be positive or
negative. This mode does not prescribe an energy injection rate.

### Energy Injection

With `normalization=edot`, `dedt` specifies the target volume-averaged energy
injection rate. For a scale factor $s$, the driver evaluates

$$
m_0=\frac{1}{V}\int \frac{1}{2}\rho|\boldsymbol{a}|^2\Delta t\,dV,
\qquad
m_1=\frac{1}{V}\int \boldsymbol{p}\cdot\boldsymbol{a}\,dV,
$$

where $\boldsymbol{p}=\rho\boldsymbol{v}$, and uses the non-negative solution to

$$
m_0s^2+m_1s=\dot e_{\rm target}.
$$

Each Newtonian kick then does work
$\Delta E_{\rm box}=\texttt{dedt}\,V\Delta t$.
This is added to total energy for an ideal EOS; an isothermal EOS does not
evolve total energy.
A zero target produces no kick. Fixed-power normalization is a feedback
prescription: its amplitude depends on the fluid state and the timestep.
For positive power, an unscaled force opposing the flow can require a
finite impulse even as the timestep tends to zero. At zero force-flow
correlation, the impulse scales as $\sqrt{\Delta t}$. Thus the exact work
budget does not establish a regular continuous-time limit for this feedback
prescription; it is not an externally prescribed OU acceleration.

`dedt` and `accel_rms` are mutually exclusive: the selected normalization
mode requires its target and rejects the other target.

## Time Integration

Before RK evolution, the Newtonian driver advances its modal OU state through
all updates due at the cycle-start time. The resulting realization is held
for that cycle. It applies one kick at the final RK stage, after the
cell-centered fluid update and source terms, before boundary exchange and
primitive recovery.

A cycle is driven when its start time is at least `tdriv_start` and, for
`turb_flag=1`, less than `tdriv_start + tdriv_duration`. Cycle zero is included
when it satisfies these conditions. The driver does not shorten timesteps
at start, stop, or OU update boundaries. An active cycle receives a full
kick even if it crosses the stopping time; a cycle starting before the
driving window receives none. Pending OU updates are processed together at
the next active cycle, retaining the prescribed recurrence and random
sequence but resolving changes only at cycle boundaries.

The kick uses the current conserved density and momentum, with the full
timestep $\Delta t$:

$$
\boldsymbol{p}'=\boldsymbol{p}+\rho\boldsymbol{a}\Delta t,\qquad
E'=E+\boldsymbol{p}\cdot\boldsymbol{a}\Delta t
+\frac12\rho|\boldsymbol{a}|^2\Delta t^2.
$$

The energy increment is exactly the kinetic-energy change. Density and
magnetic field are unchanged, so the kick preserves thermal energy to
roundoff for an ideal EOS. Isothermal fluids receive only the momentum
update. In two-fluid calculations, each fluid uses its own conserved state
and EOS. Removing the mass-weighted mean immediately before the kick also
gives zero net momentum injection, summed over both fluids when present.

For smooth fixed-RMS forcing, this single operator-split kick gives
first-order coupling to fluid evolution, even when the fluid integrator is
RK2, RK3, or RK4. Exact kick work does not establish higher-order accuracy
of the coupled evolution, or temporal convergence of the singular
fixed-power feedback cases described above.

Single-fluid special-relativistic driving keeps its previous scheduling and arithmetic,
including the first-cycle skip. The Newtonian kick and accuracy statements
above do not describe that path.

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
| `normalization` | `accel_rms` | `accel_rms` or `edot`; supply the corresponding strength explicitly. |
| `dedt` | required for `edot` | Target volume-averaged injection rate. |
| `accel_rms` | required for `accel_rms` | Target volume-weighted RMS acceleration. |
| `nlow`, `nhigh` | `1`, `3` | Inclusive driven mode-radius bounds. |
| `npeak` / `kpeak` | `kpeak=4*pi` | Parabolic spectral peak; `npeak` is tile-local mode number. |
| `spectrum` | `parabolic` | `parabolic` or `power_law`. |
| `expo`, `exp_prp`, `exp_prl` | `5/3`, `5/3`, `0` | Power-law spectrum exponents. |
| `min_kx/y/z`, `max_kx/y/z` | `0`, `nhigh` | Optional directional mode bounds. |
| `driving_type` | `0` | `0` for three-dimensional; `1` for planar driving. |
| `sol_fraction` | `1.0` | Solenoidal fraction in `[0,1]`. |
| `rseed` | `-1` | Non-negative values select reproducible sequences. |
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

See {doc}`../examples/turbulence_driving` for commands, verification checks,
and illustrative projections.
