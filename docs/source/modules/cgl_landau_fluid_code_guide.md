# CGL Landau-Fluid Code Guide

This guide maps the CGL/Landau-fluid method onto the AthenaK implementation.
For the physics model, start with
[CGL Method Physics Primer](cgl_mhd_method.md). For validation and campaign
workflows, use [CGL Landau-Fluid Validation](cgl_landau_fluid_validation.md).

## Source Map

| Area | Primary files | Role |
| --- | --- | --- |
| CGL EOS and primitive recovery | `src/eos/cgl_mhd.cpp`, `src/eos/ideal_c2p_mhd.hpp`, `src/eos/cgl_physics.hpp` | CGL pressure recovery, conserved anisotropy, finite backup relaxation, limiter predicates, heat-flux ratio helpers. |
| LF parabolic operator | `src/diffusion/cgl_landau_fluid.hpp`, `src/diffusion/cgl_landau_fluid.cpp`, `src/diffusion/cgl_landau_fluid_arithmetic.hpp` | Runtime parsing, face-state construction, LF heat-flux kernels, diagnostics, profiling, safe arithmetic. |
| MHD integration | `src/mhd/mhd.cpp`, `src/mhd/mhd_tasks.cpp`, `src/mhd/mhd_sts.cpp`, `src/mhd/rsolvers/hlle_cgl.hpp` | Construction, task graph, STS sweep lifecycle, passive/active CGL fluxes, LF-only fast paths. |
| Diagnostics and tests | `src/outputs/history.cpp`, `src/pgen/tests/cgl_landau_fluid.cpp`, `src/pgen/cgl_lf_paper.cpp`, `tst/scripts/cgl/`, `tst/test_suite/cgl/` | History columns, quantitative LF pgen, reduced paper pgen, regression workflows. |

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
   magnetic moment. On refined/shearing meshes, `RefreshCGLLFBoundaryState`
   rebuilds halos in that representation before stage receives are posted.
   This makes the first and subsequent RHSs use the same boundary operator.
2. Each STS stage calls `CGLLandauFluid::AddHeatFluxes`; LF-only stages skip
   flux clearing because the face kernels overwrite every consumed `IEN`/`IAN` face.
3. `AddHeatFluxes` precomputes `T_parallel`, `T_perp`, and `|B|`, constructs
   x1/x2/x3 face states using the supplied CT face field, evaluates capped
   parallel and perpendicular LF heat
   fluxes, and writes only `IEN` and `IAN` face fluxes. On refined meshes,
   existing coarse/fine restriction and a symmetric mean of cached same-level
   estimates give adjacent cells one shared LF flux before the update.
4. `STSUpdateU` applies the RKL2 update to energy and magnetic moment.
5. `CGLRefreshPrimFromMagneticMoment` rebuilds CGL primitives between stages.
6. `RecordAdmissibility` updates LF health counters and optionally aborts in
   strict mode.
7. `EndCGLLandauFluidSTSSweep` converts `IAN` back to conserved anisotropy
   and refreshes both cached pressures from the stored E/A state.

`CGLCollisions` applies full-cycle rates after RK only without split LF.
With split LF, `STSPostSweepCGLCollisions` applies them after the final post-LF
chunk, after A is restored. The driver invokes that post hook even when the
sweep has no stages. Other sweep hooks retain state/coarse synchronization and
timestep refresh without applying collisions. Neither regridding nor split
boundaries project anisotropy. The two-tier finite update preserves its backup
residual through ordinary active E/A or passive J/A encoding.
On refined/shearing LF meshes, the driver completes pending communication and
rebuilds anisotropy halos after each split checkpoint; it computes the post-RK
and post-LF bounds only afterward. Other CGL runs refresh the post-RK bound
inside the collision hook. Entry boundary refreshes use stage zero and do not
increment `lf_nstage`. Ordinary MHD retains its last-stage `NewTimeStep` task.
Initialization and regridding retain their existing boundary refreshes.

On uniform periodic strict-LF runs without shear/orbital advection, stage 1
exchanges the full state. Owner and halo density/scalar repairs are identical
with synchronized frozen B and EOS. Later LF stages change only IEN/IAN, so
compact communication remains valid without an extra exchange or reduction.
This includes AMR meshes whose leaves all share one level, determined from the
global leaf list during neighbor setup after initialization, restart or regridding.
Those later LF stages also skip unused coarse restriction and prolongation,
while retaining owned-cell primitive recovery and repairs, and shared-face flux
communication. Mixed-level meshes retain the full boundary path.

The scalar rate map retains both the target anisotropy and its pressure exchange.
`expm1` preserves tiny background kicks; the initially smaller pressure receives
the exchange directly, and the larger pressure is recovered from the target.
This avoids small-pressure cancellation within each selected branch; the
two-tier map has jumps at backup thresholds. The active grid operator changes
only conserved A and refreshes pressures from E/A. The passive operator
re-encodes J/A at fixed thermal energy, up to representation roundoff.

`RecordAdmissibility` counts each unprojected LF stage. Threshold crossings and
finite density/pressure repairs are recorded without stopping. Invalid states
remain fatal in strict mode. Nonfinite raw density/energy/moment and decoded LF
pressures are rejected before flooring can hide them. Entry and sweep-end checks
retain finite/positive-state validation without
incrementing stage counts or rejecting finite backup overshoots.

## LF Face Closure

`BuildCGLLFFaceState` is the face-level admission gate. It averages density and
pressures to the face, applies floors, and divides the CT-normal/averaged-tangential
field vector by the larger of its own magnitude and the mean adjacent cell
magnitudes. This retains cancellation attenuation while keeping its magnitude
at most one. It selects local or background
`c_parallel`, combines background and limiter collision frequencies, and returns
false when the mean magnitude is at or below `bfloor`. Transport, diagnostics,
and the timestep bound all use this geometry; see the
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
| `cgl_lf_arithmetic` | `safe`, `fast` | `safe` | `safe` uses overflow-protected scaled arithmetic. `fast` uses direct `Real` arithmetic for normal states and remains available for historical comparisons. |
| `cgl_lf_sts_flux` | `weighted`, `physical` | `weighted` | `weighted` embeds STS weights in LF face fluxes. `physical` writes physical fluxes and lets `STSUpdateU` apply RKL2 weights. |

`cgl_lf_sts_flux = physical` requires both
`cgl_lf_diagnostics = none` and `cgl_lf_arithmetic = fast`; construction
aborts otherwise.

For active rebuild production timing, use:

```ini
<mhd>
cgl_lf_diagnostics = none
cgl_lf_arithmetic = safe
cgl_lf_sts_flux = weighted
cgl_lf_strict_admissibility = true
cgl_lf_profile = false
```

For conservative debugging or reproducing the original protected path, use the
defaults:

```ini
<mhd>
cgl_lf_diagnostics = full
cgl_lf_arithmetic = safe
cgl_lf_sts_flux = weighted
```

Runtime modes come from the input parameters; LF mode environment overrides have
been removed. Profiling adds Kokkos fences and retains whole-direction timings.
The replay sub-kernels are removed: `cgl_lf_profile_detail=false` remains accepted
for old inputs, while `true` is rejected. Use profiling for attribution and the
ordinary elapsed interval for production timing.

For a matched throughput comparison, restart the same checkpoint on the same
rank layout with identical forcing, arithmetic, and STS settings. Compare
elapsed time over the same physical interval and compute LF evaluations from
the change in `lf_nstage` divided by the owned cell count on a fixed mesh.
This separates the cost per evaluation from a change in the number of stages.

To diagnose a cold cell, replay from before the first drop in the LF timestep.
Track its pressures and neighboring magnetic fields across the pre-LF sweep,
hyperbolic update, cooling, post-LF sweep, and scattering. If the first extreme
ratio appears inside LF, inspect that sweep's individual stages. In a cooled
run, `cool_hits` counts updates clipped by the cooling guard, while passive
`ratio_hits` counts pressure redistributions at recovery. Compare at the same
times with the original executable and forcing realization.

## STS Fast Paths

The MHD STS update has a CGL-LF-only path when no other MHD STS process is
active. In that case, state copies and the update kernel touch only `IEN` and
`IAN` instead of all MHD variables. LF-only stages skip the flux clear; generic
diffusion retains it. This preserves the same RKL2 recurrence while avoiding
full-array work for variables that LF does not advance.

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
mode, nonfinite or nonpositive recovered state is fatal at LF stages and split
checkpoints. Finite density/pressure repairs increment existing restartable
`lf_dfloor`/`lf_pfloor`; the EOS event log records `eos_dfloor`/`eos_efloor` per
completed cycle, including repairs outside LF. These event counts can include
halos and repeated refreshes; they are not unique cells or mass/energy budgets.
Threshold crossings also remain diagnostic. No restart field or diagnostic pass
is added for the nonfatal repair policy.
Backup relaxation and scattering require an explicit `backup_limiters = true` in
both strict and relaxed modes.

## Active And Passive CGL

Active CGL runs apply the anisotropic pressure tensor in the dynamic MHD fluxes.
Passive-Delta runs evolve CGL/LF pressures diagnostically while the flow uses
the native MHD passive path (`passive_eos=isothermal` or `adiabatic`). The passive branch affects Riemann fluxes and
timestep estimates, so active/passive comparisons should be documented as a
model-path comparison rather than a single-line source-term toggle. The
[passive model contract](cgl_passive.md) specifies J/A outside LF, U/moment
inside LF, restart encoding and the absence of irreversible shock heating in
the passive CGL pair. Adiabatic passive flow carries its own MHD energy.

## Restrictions

Current restrictions are intentionally conservative:

- CGL is Newtonian MHD only; SR, GR, and dynamical-GR MHD reject it.
- CGL dynamic runs use `rsolver = hlle`.
- Ordinary isotropic conduction is incompatible with `eos = cgl`.
- CGL LF STS cannot be combined with another MHD STS process in the same run.
- Explicit LF is a reference mode and cannot be combined with another active
  MHD parabolic process.
- Active CGL LF with AMR supports primitive fine/coarse boundary prolongation in both
  the ordinary anisotropy representation and the temporary LF STS magnetic-
  moment representation. Mesh creation and deletion use CGL-specific primitive
  rebuild wrappers rather than treating `IAN` as a passive scalar. Active block
  creation conservatively prolongates density, momentum, total energy, and
  scalar densities, then rebuilds the CGL encoding from prolonged pressure
  anisotropy and the finalized magnetic field. Valid primary conserved fields
  are preserved. AMR does not apply a limiter projection.
- Earlier LF/STS AMR refine-up and forced moving-pattern refine/derefine
  smokes passed on GPUs with the former scheduled-wall operator; they do not
  qualify the restored finite-backup model.
  Dedicated CGL AMR repair counters report whether AMR transfer needed density,
  energy, pressure, low-field, anisotropy, interval, or slope-scaling repairs.
  The legacy hard-wall counters remain reserved and zero.
- Passive CGL and positive LF chunk caps reject AMR/SMR; see the
  [passive contract](cgl_passive.md) and [STS controls](super_time_stepping.md).
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
