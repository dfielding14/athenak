# CGL Landau-Fluid Code Guide

This guide maps the CGL/Landau-fluid method onto the AthenaK implementation.
For the physics model, start with
[CGL Method Physics Primer](cgl_mhd_method.md). For validation and campaign
workflows, use [CGL Landau-Fluid Validation](cgl_landau_fluid_validation.md).

## Source Map

| Area | Primary files | Role |
| --- | --- | --- |
| CGL EOS and primitive recovery | `src/eos/cgl_mhd.cpp`, `src/eos/ideal_c2p_mhd.hpp`, `src/eos/cgl_physics.hpp` | CGL pressure recovery, conserved anisotropy, hard-wall projection, limiter predicates, heat-flux ratio helpers. |
| LF parabolic operator | `src/diffusion/cgl_landau_fluid.hpp`, `src/diffusion/cgl_landau_fluid.cpp`, `src/diffusion/cgl_landau_fluid_arithmetic.hpp` | Runtime parsing, face-state construction, LF heat-flux kernels, diagnostics, profiling, safe arithmetic. |
| MHD integration | `src/mhd/mhd.cpp`, `src/mhd/mhd_tasks.cpp`, `src/mhd/mhd_sts.cpp`, `src/mhd/rsolvers/hlle_cgl.hpp` | Construction, task graph, STS sweep lifecycle, passive/active CGL fluxes, LF-only fast paths. |
| Diagnostics and tests | `src/outputs/history.cpp`, `src/pgen/tests/cgl_landau_fluid.cpp`, `src/pgen/tests/cgl_lf_paper.cpp`, `tst/scripts/cgl/`, `tst/test_suite/cgl/` | History columns, quantitative LF pgen, reduced paper pgen, regression workflows. |

The LF closure is registered as the MHD parabolic process
`mhd/cgl_heat_flux`. It is not ordinary isotropic conduction, and ordinary
`<mhd>/conductivity` is rejected when `eos = cgl`.

## Conserved And Primitive State

In ordinary CGL MHD state, AthenaK uses:

| Slot | Meaning |
| --- | --- |
| `IDN` | density |
| `IM1`, `IM2`, `IM3` | momentum |
| `IEN` | total energy |
| `IAN` | conserved CGL anisotropy outside LF sweeps |
| `IPR` | primitive parallel pressure `p_parallel` |
| `IPP` | primitive perpendicular pressure `p_perp` |

The conserved anisotropy is produced by `CGLConservedAnisotropy()` in
`src/eos/ideal_c2p_mhd.hpp` and recovered during CGL primitive conversion in
`src/eos/cgl_mhd.cpp`.

During a Landau-fluid split sweep, the `IAN` slot is temporarily reinterpreted
as magnetic moment:

$$
\mu_B = p_\perp / |B|.
$$

The sweep begins by converting anisotropy to magnetic moment, refreshes
primitive variables from magnetic moment between STS stages, and converts back
to anisotropy at the end. This is why normal hyperbolic fluxes, outputs, and
restart state must never observe a half-completed LF sweep.

## Runtime Setup

The CGL/LF path is selected in the `<mhd>` block:

```ini
<time>
sts_integrator = rkl2

<mhd>
eos = cgl
rsolver = hlle
cgl_heat_flux = landau_fluid
cgl_heat_flux_integrator = sts
lf_k_parallel = 32.0
lf_coefficient_mode = local
```

The constructor for `CGLLandauFluid` parses the LF-specific parameters and
rejects invalid closure scales, coefficient modes, runtime modes, and
incompatible fast-path combinations.

## Split-Sweep Lifecycle

For STS LF transport, MHD owns the split lifecycle:

1. `BeginCGLLandauFluidSTSSweep` converts `IAN` from conserved anisotropy to
   magnetic moment.
2. Each STS stage clears the LF flux slots and calls
   `CGLLandauFluid::AddHeatFluxes`.
3. `AddHeatFluxes` precomputes `T_parallel`, `T_perp`, and `|B|`, constructs
   x1/x2/x3 face states using the supplied CT face field, evaluates capped
   parallel and perpendicular LF heat
   fluxes, and writes only `IEN` and `IAN` face fluxes.
4. `STSUpdateU` applies the RKL2 update to energy and magnetic moment.
5. `CGLRefreshPrimFromMagneticMoment` rebuilds CGL primitives between stages.
6. `RecordAdmissibility` updates LF health counters and optionally aborts in
   strict mode.
7. `EndCGLLandauFluidSTSSweep` converts `IAN` back to conserved anisotropy
   and refreshes both cached pressures from the stored E/A state.

`CGLCollisions` applies walls after RK for every CGL configuration. Without
split LF it also applies rates with the full cycle timestep.
`STSPostSweepCGLCollisions` applies walls after every CGL parabolic sweep.
For non-LF STS processes, such as viscosity or resistivity, both sweep
checkpoints use `walls_only`; rates remain in the post-RK hook. With split LF,
A is restored first, and the post-sweep checkpoint applies full-step rates
before walls. The driver also invokes that LF post checkpoint when the
post-sweep has no stages. Each sweep checkpoint follows the final field and
primitive refresh. Both collision hooks refresh coarse state after their
updates, and the post-sweep timestep estimate uses the corrected state.
The post-RK hook refreshes the local
timestep from the corrected state before the driver's global minimum.
It is the sole post-RK CGL timestep owner; ordinary MHD retains its last-stage
`NewTimeStep` task. Initialization, regridding, and post-parabolic refreshes remain.

The scalar rate map retains both the target anisotropy and its pressure exchange.
`expm1` preserves tiny background kicks; the initially smaller pressure receives
the exchange directly, and the larger pressure is recovered from the target.
This avoids small-pressure cancellation while retaining the monotone stiff limit.
The grid operator still changes only conserved A and refreshes pressures from E/A.

`RecordAdmissibility` counts each unprojected LF stage. Intermediate wall
crossings are diagnostic; invalid pressures and floors remain fatal in strict
mode. Separate entry/post-wall checks enforce walls without incrementing stage
counts. After AMR, walls run after the final field refresh.

## LF Face Closure

`BuildCGLLFFaceState` is the face-level admission gate. It averages density and
pressures to the face, applies floors, and divides the CT-normal/averaged-tangential
field vector by the mean adjacent cell magnitudes. It selects local or background
`c_parallel`, combines background and limiter collision frequencies, and returns
false when the mean magnitude is at or below `bfloor`. Transport, diagnostics,
and profile probes all use this geometry; see the
[face regularization](cgl_landau_fluid.md) for its physical limits.

The closure formulas are shared through `src/eos/cgl_physics.hpp`:

- `ParallelHeatFluxRatio`
- `PerpendicularHeatFluxRatio`
- `LimitedParallelHeatFlux`
- `LimitedPerpendicularHeatFlux`

`CGLLFFlux` is the safe path using `ScaledValue` arithmetic from
`cgl_landau_fluid_arithmetic.hpp`. `CGLLFFluxFast` is the normal-range
production path using direct `Real` arithmetic after the same ratio and cap
logic.

## Runtime Modes

The current implementation exposes three LF performance/safety switches:

| Parameter | Choices | Default | Notes |
| --- | --- | --- | --- |
| `cgl_lf_diagnostics` | `full`, `none` | `full` | `full` collects q-face, cap, and q-work diagnostics with reductions. `none` skips those reductions and leaves heat-flux diagnostic columns inactive while preserving admissibility counters. |
| `cgl_lf_arithmetic` | `safe`, `fast` | `safe` | `safe` uses overflow-protected scaled arithmetic. `fast` uses direct `Real` arithmetic for normal states and is the production performance path. |
| `cgl_lf_sts_flux` | `weighted`, `physical` | `weighted` | `weighted` embeds STS weights in LF face fluxes. `physical` writes physical fluxes and lets `STSUpdateU` apply RKL2 weights. |

`cgl_lf_sts_flux = physical` requires both
`cgl_lf_diagnostics = none` and `cgl_lf_arithmetic = fast`; construction
aborts otherwise.

For production performance timing, use:

```ini
<mhd>
cgl_lf_diagnostics = none
cgl_lf_arithmetic = fast
cgl_lf_sts_flux = physical
cgl_lf_profile = false
cgl_lf_profile_detail = false
```

For conservative debugging or reproducing the original protected path, use the
defaults:

```ini
<mhd>
cgl_lf_diagnostics = full
cgl_lf_arithmetic = safe
cgl_lf_sts_flux = weighted
```

The runtime modes can also be overridden by environment variables:

```bash
ATHENAK_CGL_LF_DIAGNOSTICS=none
ATHENAK_CGL_LF_ARITHMETIC=fast
ATHENAK_CGL_LF_STS_FLUX=physical
ATHENAK_CGL_LF_PROFILE=0
ATHENAK_CGL_LF_PROFILE_DETAIL=0
```

Profiling mode intentionally adds Kokkos fences, and detailed profiling launches
extra probe kernels. It should be used for attribution, not for production wall
time.

## STS Fast Paths

The MHD STS update has a CGL-LF-only path when no other MHD STS process is
active. In that case, flux clearing, state copies, and the update kernel touch
only `IEN` and `IAN` instead of all MHD variables. This preserves the same RKL2
recurrence while avoiding full-array work for variables that LF does not
advance.

With `cgl_lf_sts_flux = physical`, LF face fluxes are ordinary physical fluxes.
`STSUpdateU` then applies

$$
\Delta U = -dt_{\rm sweep}\,\nabla\cdot F
$$

and the normal RKL2 coefficients in the same update expression used by other
cell-centered STS processes. The weighted path remains available for debugging
and robustness comparisons.

## Diagnostics

`CGLLFDiagnostics` stores cumulative LF counters. Normal `.mhd.hst` output
appends LF columns when the closure is active:

- stage and admissibility counters: `lf_nstage`, `lf_dfloor`, `lf_pfloor`,
  `lf_nonfin`, `lf_nonpos`, `lf_mirror`, `lf_firehs`, `lf_hardbd`,
  `lf_hwproj`;
- heat-flux diagnostics in full mode: `lf_qface`, `lf_qprcap`, `lf_qpr10`,
  `lf_qpecap`, `lf_qpe10`, `lf_qprwrk`, `lf_qpewrk`;
- optional retained pressure-work diagnostics: `lf_cpwrk`, `lf_cawrk` when
  `cgl_lf_record_pressure_work = true`.

Strict admissibility is controlled by `cgl_lf_strict_admissibility`. In strict
mode, floors and non-finite or non-positive thermodynamic state are fatal at
every LF stage; hard-bound violations are fatal at split-sweep checkpoints.
Backup walls and scattering require an explicit `backup_limiters = true` in
both strict and relaxed modes.

## Active And Passive CGL

Active CGL runs apply the anisotropic pressure tensor in the dynamic MHD fluxes.
Passive-Delta runs evolve CGL/LF pressures diagnostically while the flow uses
the isothermal-MHD passive path. The passive branch affects Riemann fluxes and
timestep estimates, so active/passive comparisons should be documented as a
model-path comparison rather than a single-line source-term toggle.

## Restrictions

Current restrictions are intentionally conservative:

- CGL is Newtonian MHD only; SR, GR, and dynamical-GR MHD reject it.
- CGL dynamic runs use `rsolver = hlle`.
- Ordinary isotropic conduction is incompatible with `eos = cgl`.
- CGL LF STS cannot be combined with another MHD STS process in the same run.
- Explicit LF is a reference mode and cannot be combined with another active
  MHD parabolic process.
- CGL LF with AMR supports primitive fine/coarse boundary prolongation in both
  the ordinary anisotropy representation and the temporary LF STS magnetic-
  moment representation. Mesh creation and deletion use CGL-specific primitive
  rebuild wrappers rather than treating `IAN` as a passive scalar. Active block
  creation now prolongates `rho`, velocity, CGL internal energy `U`, and
  pressure anisotropy `Delta`, then slope-scales child states before final
  admissibility projection.
- Current-threshold LF/STS AMR refine-up and forced moving-pattern
  refine/derefine churn smokes have passed on GPUs with strict admissibility.
  Dedicated CGL AMR repair counters report whether AMR transfer needed density,
  energy, pressure, low-field, hard-wall, anisotropy, interval, or slope-scaling
  repairs.
- CGL LF pressure-work recording is disabled for AMR primitive prolongation.
  Keep `<mhd>/cgl_lf_record_pressure_work = false` for LF/STS AMR runs until
  the diagnostic flux-communication path is separately audited.

## Validation Hooks

Use the built-in `cgl_landau_fluid` pgen for local LF operator tests and
`cgl_lf_paper` for reduced Majeski/Squire-style turbulence mechanics. Routine
test entry points live under `tst/scripts/cgl/` and `tst/test_suite/cgl/`.

The user-facing validation workflows and interpretation limits are documented
in [CGL Landau-Fluid Validation](cgl_landau_fluid_validation.md). Do not treat
the presence of input decks, source code, or a private qualitative comparison
as proof that a production comparison has been completed.
