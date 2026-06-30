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
   x1/x2/x3 face states, evaluates capped parallel and perpendicular LF heat
   fluxes, and writes only `IEN` and `IAN` face fluxes.
4. `STSUpdateU` applies the RKL2 update to energy and magnetic moment.
5. `CGLRefreshPrimFromMagneticMoment` rebuilds CGL primitives between stages.
6. `RecordAdmissibility` updates LF health counters and optionally aborts in
   strict mode.
7. `EndCGLLandauFluidSTSSweep` converts `IAN` back to conserved anisotropy.

When LF and CGL collisions or limiter scattering are active, the driver applies
collision updates after each LF half-sweep. This gives the chronological split
`L(dt/2) C(dt/2) H(dt) L(dt/2) C(dt/2)`.

## LF Face Closure

`BuildCGLLFFaceState` is the face-level admission gate. It averages density and
pressures to the face, applies floors, constructs the field direction, selects
local or background `c_parallel`, combines background and limiter collision
frequencies, and returns false for faces with unusable magnetic-field strength.

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
mode, floors, non-finite or non-positive thermodynamic state, or emergency hard
bound violations abort the run during LF primitive refresh. In relaxed mode,
the emergency backup limiter is enabled when instability limiters are active.

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
- CGL LF with AMR requires conserved prolongation; primitive prolongation is
  rejected.

## Validation Hooks

Use the built-in `cgl_landau_fluid` pgen for local LF operator tests and
`cgl_lf_paper` for reduced Majeski/Squire-style turbulence mechanics. Routine
test entry points live under `tst/scripts/cgl/` and `tst/test_suite/cgl/`.

The user-facing validation workflows and interpretation limits are documented
in [CGL Landau-Fluid Validation](cgl_landau_fluid_validation.md). Do not treat
the presence of input decks, source code, or a private qualitative comparison
as proof that a production comparison has been completed.
