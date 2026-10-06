# CGL Landau-Fluid Heat Flux

## Scope

AthenaK supports non-relativistic anisotropic MHD with the CGL equation of
state and an optional Landau-fluid (LF) heat-flux closure. The feature is
implemented as a distinct MHD parabolic process, `mhd/cgl_heat_flux`, rather
than as ordinary isotropic thermal conduction.

For a physics-first explanation of the model, see
[CGL Method Physics Primer](cgl_mhd_method.md). For the source-level map,
runtime modes, and performance switches, see
[CGL Landau-Fluid Code Guide](cgl_landau_fluid_code_guide.md). For validation
workflows and evidence boundaries, see
[CGL Landau-Fluid Validation](cgl_landau_fluid_validation.md). For static and
adaptive refinement, see [CGL Adaptive Mesh Refinement](cgl_amr.md).

Use CGL with the HLLE solver:

```ini
<time>
sts_integrator = rkl2
sts_safety = 0.9

<mhd>
eos = cgl
passive = false
rsolver = hlle
gamma = 1.666666666666667
cgl_heat_flux = landau_fluid
cgl_heat_flux_integrator = sts
lf_k_parallel = 32.0
lf_coefficient_mode = local
firehose_threshold = 2.0
mirror_threshold = 1.0
cgl_lf_strict_admissibility = false
```

`lf_k_parallel` is the closure wavenumber magnitude, not a conductivity.
Set `lf_coefficient_mode = background` with a positive `lf_c_parallel0` to
hold the closure speed fixed while gradients and evolved variables remain
local.

## State And Update

CGL evolves density, momentum, total energy, and a sixth MHD state storing
conserved pressure anisotropy outside the LF sweep. Primitive arrays use
`IPR` for parallel pressure and `IPP` for perpendicular pressure.

The LF flux advances total energy and magnetic moment. During each LF split
sweep the MHD task graph performs:

1. Convert conserved anisotropy to magnetic moment.
2. Evaluate field-aligned LF energy and magnetic-moment fluxes.
3. Advance those two conserved quantities through either RKL2 STS stages or
   one explicit reference stage.
4. Refresh CGL primitives between stages in the magnetic-moment
   representation.
5. Convert magnetic moment back to conserved anisotropy at sweep completion.

This lifecycle prevents ordinary hyperbolic fluxes, output, and restart state
from interpreting magnetic moment as pressure anisotropy.

CGL collision rates run once per full timestep. With LF, the pre-sweep and
hyperbolic boundaries apply only the configured backup walls and the unconditional
fluid firehose wall; the post-sweep applies exact background decay and monotone
backward-Euler soft-limiter relaxation over the full timestep, then those walls.
Without LF, rates and walls run after the hyperbolic timestep. Primitive recovery
does not apply soft-threshold scattering.

AthenaK records LF-stage health metrics in normal MHD history output whenever
this closure is active. Set `cgl_lf_strict_admissibility = true` for
verification runs to terminate immediately if an LF refresh produces
non-finite or non-positive thermodynamic state, activates density or pressure
floors. Hard walls are checked strictly at sweep entry and after the scheduled
end-of-sweep projection. Intermediate LF crossings remain counted in `lf_hardbd`
but do not abort: LF can cross a wall before that projection, even from an
admissible initial state. The fluid firehose bound is always checked; configured
backup walls are checked only when `backup_limiters` is enabled. No per-stage
wall projection or soft scattering is added.

At a face, $\nu_{\rm eff}$ is the sum of background collisions, one
`limiter_nu_coll` contribution when an enabled soft threshold is exceeded, and
one `limiter_backup_nu` contribution when backup limiting is enabled and a
configured or fluid wall is exceeded. The two limiter contributions add.

The perpendicular closure uses the BGK moment coefficient from SHD97 eq. 49
and Sharma et al. (2006) eq. 12:

$$
\chi_\perp = \frac{2c_\parallel^2}
{\sqrt{2\pi}c_\parallel k_\parallel+2\nu_{\rm eff}}.
$$

The factor $2\nu_{\rm eff}$ deliberately differs from the $+\nu_{\rm eff}$
printed in Squire et al. (2023) eq. 2.7. It recovers
$\chi_\perp\to c_\parallel^2/\nu_{\rm eff}$ for strong collisions.

## Closure Controls

| Parameter | Default | Meaning |
| --- | --- | --- |
| `eos` | required | Set to `cgl` for this feature. |
| `passive` | `false` | `true` is disabled: the passive thermal energy equation is inconsistent, pending WO2. |
| `iso_sound_speed` | unused in active CGL | Retained for the disabled passive model; it does not bypass the `passive=true` fence. |
| `cgl_heat_flux` | absent | Set to `landau_fluid` to enable LF transport. |
| `cgl_heat_flux_integrator` | `sts` | `sts` for production runs or `explicit` for reference verification. |
| `lf_k_parallel` | required | Finite positive closure wavenumber magnitude. |
| `lf_coefficient_mode` | `local` | `local` or `background`. |
| `lf_c_parallel0` | required in background mode | Finite positive fixed parallel thermal speed. |
| `nu_coll` | `0.0` | Background anisotropy-relaxation frequency. |
| `mirror_limiter` | `false` | Enable mirror-limiter relaxation. |
| `firehose_limiter` | `false` | Enable firehose-limiter relaxation. |
| `firehose_threshold` | `2.0` | Positive $\Lambda_{\rm FH}$; soft firehose threshold $\Delta p=-\Lambda_{\rm FH}B^2/2$. |
| `mirror_threshold` | `1.0` | Positive $\Lambda_{\rm M}$; soft mirror threshold $\Delta p=+\Lambda_{\rm M}B^2/2$. |
| `mirror_backup_factor` | `2.0` | Multiplier of the soft mirror threshold; at least 1. |
| `firehose_backup_factor` | `1.0` | Multiplier of the soft firehose threshold; at least 1. The wall cannot lie below $-B^2$. |
| `limiter_backup_nu` | `1e10` | Nonnegative LF heat-flux suppression frequency in inverse code time. |
| `cgl_firehose_threshold` | absent | Legacy alias: `oblique` = 1.4, `parallel` = 2.0. Conflicting explicit numeric values are rejected. |
| `limiter_nu_coll` | required when either limiter is enabled | Nonnegative finite soft-limiter relaxation frequency. |
| `limiter_hardwall` | `false` | Legacy `true` is rejected. Use finite `limiter_nu_coll` for soft-threshold relaxation. |
| `backup_limiters` | `false` | Project onto both configured backup walls; requires at least one soft limiter to be enabled. LF never enables this flag implicitly. |
| `cgl_lf_strict_admissibility` | `false` | Fail immediately on LF floors or invalid states; enforce hard walls at sweep entry and after the scheduled end projection. |
| `cgl_lf_record_pressure_work` | `false` | Retain RK-integrated applied CGL pressure-traction work diagnostics. |
| `cgl_lf_diagnostics` | `full` | `full` collects heat-flux face/cap/work diagnostics; `none` skips those reductions for production runs. |
| `cgl_lf_arithmetic` | `safe` | `safe` uses overflow-protected scaled arithmetic; `fast` uses direct normal-range `Real` arithmetic. |
| `cgl_lf_sts_flux` | `weighted` | `weighted` embeds STS weights in LF fluxes; `physical` writes physical fluxes and requires `cgl_lf_diagnostics = none` plus `cgl_lf_arithmetic = fast`. |
| `cgl_lf_profile` | `false` | Enable coarse LF timing regions with Kokkos fences; for profiling only. |
| `cgl_lf_profile_detail` | `false` | Add directional probe kernels for detailed profiling when `cgl_lf_profile = true`. |

The performance/safety switches can be overridden by
`ATHENAK_CGL_LF_DIAGNOSTICS`, `ATHENAK_CGL_LF_ARITHMETIC`,
`ATHENAK_CGL_LF_STS_FLUX`, `ATHENAK_CGL_LF_PROFILE`, and
`ATHENAK_CGL_LF_PROFILE_DETAIL`. Production wall-time measurements should keep
profiling disabled. An available performance configuration is:

```ini
<mhd>
cgl_lf_diagnostics = none
cgl_lf_arithmetic = fast
cgl_lf_sts_flux = physical
cgl_lf_profile = false
cgl_lf_profile_detail = false
```

At an operator face with `|B| <= bfloor`, LF does not construct a local field
direction and applies zero heat-flux contribution at that face. The local
accuracy campaign exercises this shutdown behavior with strict monitoring.

The paper inputs explicitly set `firehose_threshold = 2.0` and
`mirror_threshold = 1.0`, corresponding to $\Delta p=-B^2$ and $B^2/2$.
Their stiff soft-limiter rate is `limiter_nu_coll = 1e10`. The finite update
approaches the soft threshold with residual
$(\Delta p_0-\Delta p_{\rm threshold})/(1+\nu_{\rm lim}\,dt)$;
it is applied once per cycle. This replaces the former
`limiter_hardwall = true` soft projection during primitive recovery and AMR
transfer. Remove that legacy setting when migrating an input, and set
`backup_limiters` explicitly. Existing finite limiter rates are never overwritten
by the parser. LF strictness does not change the backup setting.

`lf_mirror` and `lf_firehs` count physical threshold occupancy. `lf_hardbd`
counts violations of the unconditional fluid wall and enabled backup walls.
The historical `lf_hwproj` column remains present and is zero for this closure.

Normal `.mhd.hst` output appends cumulative columns when LF is active:
`lf_nstage`, `lf_dfloor`, `lf_pfloor`, `lf_nonfin`, `lf_nonpos`,
`lf_mirror`, `lf_firehs`, `lf_hardbd`, `lf_qface`, `lf_qprcap`,
`lf_qpr10`, `lf_qpecap`, `lf_qpe10`, `lf_qprwrk`, `lf_qpewrk`, and
`lf_hwproj`. When
`cgl_lf_record_pressure_work = true`, it additionally appends `lf_cpwrk`
and `lf_cawrk`. The
face-count columns record owned LF faces and parallel/perpendicular unlimited
heat-flux ratios exceeding `q_max` or `10*q_max`; a shared face is retained
once, and a coarse/fine interface is owned by its fine-side closure faces.
Differences between successive rows give interval counts; normalize limiter
counts by `lf_nstage` and heat-flux-cap counts by `lf_qface`. All cumulative LF
diagnostic columns are preserved through CGL-LF restart files so interval
analysis remains continuous across segments. `lf_hwproj` is retained for old
history readers and is zero in new runs; old nonzero values describe the
pre-WO1 primitive-recovery soft projection.
`lf_qprwrk` and `lf_qpewrk` are cumulative RKL2-applied owned-face
contractions of the capped heat fluxes with their corresponding temperature
jumps. They characterize the closure-generated face fluxes; shearing-box
radial flux reconciliation occurs afterward, so these two channel diagnostics
are not exact seam-applied contractions. On refined meshes, each coarse/fine
interface uses the same fine-side ownership convention as the face-count
columns because those fluxes are restricted into the coarse update.
These are signed operator contractions; they are not required to be positive,
equal an offline snapshot proxy, or close a total energy budget. The existing
`aam-D` history column remains the conserved anisotropy variable for
compatibility.

When `cgl_lf_diagnostics = none`, LF admissibility and threshold-occupancy
counters remain active, but heat-flux face, cap, and q-work reductions are not
collected. In that mode the corresponding heat-flux diagnostic columns should
be treated as intentionally inactive rather than as measured zero cap
occupancy.

When `cgl_lf_record_pressure_work = true`, `lf_cpwrk` is the cumulative
explicit-RK-applied contraction of velocity with the retained CGL
pressure-traction divergence, and `lf_cawrk` is its `Delta p` anisotropic
component. The retained face traction is corrected through the same AMR flux
exchange used by the momentum update before the contraction is evaluated.
Archived passive-Delta runs retain zeros for both fields because their
diagnostic CGL pressures were not applied to flow momentum. New passive runs
are disabled.

## LF Timestep

LF estimates its parabolic reference step from the actual face stencil. The
estimate includes staggered normal magnetic fields, directional cell widths,
face-to-cell density ratios, transverse VL4 limiter derivatives, and the
coupling between parallel and perpendicular temperatures from magnetic-field
gradients. A normal face field divided by the mean neighboring magnetic
magnitude need not have magnitude at most one; replacing it with a normalized
cell-centered direction can miss strong discrete stiffness.

The RKL2 controller uses this reference step multiplied by
`<time>/sts_safety`, whose default is `0.9` and allowed range is finite
`0 < sts_safety <= 1`. This factor applies to every process assigned to STS.
The advective CFL and explicitly integrated parabolic processes continue to use
`<time>/cfl_number`. Both initial stage selection and the refreshed post-sweep
budget use `sts_safety`; restart files retain the setting. The optional
`sts_max_dt_ratio` cycle cap uses this safety-scaled budget.

The estimate bounds the magnitude of the local frozen temperature Jacobian.
It does not establish a nonlinear RKL2 stability theorem for state-dependent
closure coefficients, scattering switches, projection walls, or the composite
AMR operator. See [LF Discrete Timestep Bound](cgl_lf_timestep.md) for the
assumptions, derivation, exceptional arithmetic, and regression cases.

## Current Restrictions

- CGL is not available for SR, GR, or dynamical-GR MHD.
- `mhd/passive = true` is disabled pending the WO2 thermal-energy redesign.
- LF split integration rejects inflow and user boundary conditions because
  they do not have a magnetic-moment-aware `IAN` contract. Periodic, outflow,
  reflecting, and diode boundaries remain supported.
- CGL dynamic runs use `rsolver = hlle`; LLF and HLLD are rejected.
- Ordinary `<mhd>/conductivity` is rejected with `eos = cgl`.
- CGL LF with `cgl_heat_flux_integrator = sts` cannot be combined with
  another MHD process selecting STS in the same run.
- CGL LF with `cgl_heat_flux_integrator = explicit` runs the same protected
  LF split lifecycle with a one-stage Euler half-sweep. It requires
  `sts_integrator = none` and cannot yet be combined with another active MHD
  parabolic process.
- CGL LF with mesh refinement supports both conserved prolongation and
  CGL-aware primitive prolongation. With
  `<mesh_refinement>/prolong_primitives = true`, fine/coarse primitive
  boundaries and active block creation/deletion rebuild the CGL thermodynamic
  state instead of treating `IAN` as a passive scalar. During LF STS stages the
  AMR path is explicit about the temporary magnetic-moment representation
  `IAN = p_perp/|B|`.
- CGL LF pressure-work recording remains disabled with AMR primitive
  prolongation because that diagnostic flux-communication path has not been
  separately audited. Keep `<mhd>/cgl_lf_record_pressure_work = false` for
  LF/STS AMR primitive-prolongation runs.
- CGL LF STS is compatible with three-dimensional shearing-periodic
  boundaries on uniform grids. Each parabolic stage reconciles radial heat
  fluxes across the shearing seam, updates and remaps the conserved CGL state
  while `IAN` stores magnetic moment, then refreshes CGL primitives before
  the next stage. The pre-sweep uses the displacement at `t` and the
  post-sweep uses `t + dt`.
- Modal `<turb_driving>` forcing is supported with CGL LF. One kick and one
  Ornstein-Uhlenbeck advance occur per cycle before the RK state copy. The
  momentum kick and zero-net-momentum correction each add their exact kinetic
  energy change to total energy, preserving thermal energy. This schedule
  applies to all driven fluids, including two-fluid runs.

## Verification

Focused unit problems in `inputs/unit_tests/` exercise CGL transforms, CGL
FOFC, LF parallel and perpendicular decay, magnetic-field-gradient coupling,
flux limiting, limiter suppression, and a field-aligned wave. Routine
regressions also exercise analytic uniform collisional relaxation and both
firehose threshold policies. The strengthened acceptance suite additionally
checks noninteger-period wave convergence, one-e-fold decay, exact cellwise
limiter relaxation, oblique two-dimensional decay across block/rank layouts,
and two-level SMR decay with total-energy conservation. See the
{ref}`measured acceptance table <wo1-acceptance-checks>`. The LF
quantitative pgen is the built-in `src/pgen/tests/cgl_landau_fluid.cpp`.

A reduced forced-turbulence initializer is registered as
`pgen_name = cgl_lf_paper`. Its active smoke deck is
`inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput`; the passive deck
is retained as a disabled reference pending WO2. These decks initialize `rho0 = 1`, `B0` along `z`, and
`p_parallel0 = p_perp0 = beta0 B0^2/2`, use the explicit MKS24
`firehose_threshold = 2.0`, `mirror_threshold = 1.0` policy, and exercise the shared
turbulence driver. `mhd/passive = true` now fails at construction because its
thermal energy equation is inconsistent. Runtime regressions check this fence;
direct unit checks retain coverage of the isothermal passive signal-speed path.
These are reduced smoke cases, not standard paper-resolution runs. Their
forcing-orientation, seed-continuation, and multi-cycle OU/RK source-work checks
qualify reduced mechanics, not paper-scale statistics or
figure diagnostics.

The paper pgen `.user.hst` output retains volume-integrated mass, kinetic,
magnetic and CGL thermal energies, `b2`, `b4`, `delta_p`, `abs_delta_p`,
local-beta integral, mirror/firehose/hard-bound volumes, effective collision
rate, and instantaneous forcing power (written using the compact history
labels `therm_cgl`, `abs_dp`, `mirror_vol`, `fire_vol`, `hard_vol`, and
`force_pwr`) in addition to forcing-orientation quantities. Paper inputs
also enable `record_injected_work`, adding cumulative exact net forcing-source
work, including the zero-net-momentum projection, as `force_work`. Each once-per-cycle kick adds its measured energy change to this counter
before the RK stages; the analyzer uses it with conserved `tot-E` to report
an active-CGL global energy residual. These histories support reduced global
summaries such as `C_B2`; `.mhd.hst` supplies operator-face heat-flux-cap
activity and retained snapshots supply spatial diagnostics.

Paper-standard input definitions and the limiter-frequency scan live under
`inputs/cgl_lf_paper/`. They encode the standard `192x192x384` domain,
duration, stiff finite-rate baseline, physical forcing shell, binary snapshot
cadence, and analysis window. The nine standard definitions cover all eight
active/passive, Alfvenic/random beta-10/beta-100 series in MKS24 Figure 2(b)
plus the active Alfvenic beta-1 case. Two `paper-heat-flux` definitions supply
the nonnominal active beta-10 Figure 12 variants; its nominal active and
passive comparisons reuse standard definitions. Two `paper-scale-separation`
definitions supply the nonstandard `n_perp = 96` and `384` Figure 11 cases;
its `n_perp = 192` comparison reuses the standard active Alfvenic beta-10
definition. Two `paper-compressive` definitions supply the active random
beta-1 and beta-100 sonic-correlation Figure 3 cases; four other Figure 3
cases reuse standard definitions. The `paper-standard`, `paper-nulim`,
`paper-heat-flux`, `paper-compressive`, and `paper-scale-separation`
workflows require explicit production authorization;
the presence of these decks is not evidence that paper-scale runs have been
executed. Passive definitions are retained for provenance, but executable
workflows omit them and list them in `disabled_cases` until WO2.

For CGL `mhd_w` or `mhd_w_bcc` output, the existing `eint` field retains its
legacy meaning of `p_parallel`; output now also includes `p_perp`. Paper
snapshot analysis must use both fields when constructing `Delta p`.

`python3 scripts/cgl_lf_workflow.py paper-analyze --output-dir <bundle>`
regenerates reduced-history summaries and, for retained binary snapshots,
current PDF, compressive-velocity, normalized-density,
thermal/magnetic-pressure, local-field pressure/velocity-gradient spectral,
pressure-transfer, CGL pressure-work decomposition, and alignment
diagnostics with synthetic numerical checks. The pressure-work product splits
`p_perp div(u)` from `-Delta p (b b : grad u)` and compares the latter with
the direct transfer integral. The pressure-transfer product additionally
reports the MKS24 dimensionless ratio
`T_Delta_p/T_total`, using the stated estimate
`T_total ~= E_K (2 pi u_rms/L_perp)`. Multi-snapshot ensembles include trapezoidal
time-integral estimates for this product and the reconstructed heat-flux
proxy; both remain snapshot-derived. For decks enabling
`cgl_lf_record_pressure_work`, interval history analysis additionally reports
the applied hyperbolic pressure ledger from `lf_cpwrk` and `lf_cawrk`.
The workflow also renders available diagnostic figures under
`figures/paper/`. For
paper-standard bundles it uses each case's declared analysis window to
produce ensemble-average products and interval heat-flux-cap fractions.
With `--eddy-samples <count>`, `paper-analyze` additionally computes
deterministic, local-field-conditioned three-point structure functions and
emits `eddy_anisotropy.velocity_perp` and
`eddy_anisotropy.magnetic_perp`; `--eddy-bins` and `--eddy-seed` are
archived with the analysis products.
If a pinned MKS24 staging manifest is available, the analysis bundle retains
its archive and source-TeX checksums as reference provenance.
The snapshot products also include Figure 2(a)-coordinate joint PDFs of
`delta p_parallel` and `delta p_perp` versus
`<p> delta rho/<rho>`, rendered for each analyzed case. This supplies the
joint diagnostic in AthenaK units. The paper raster is decoded into sampled
surfaces with `scripts/digitize_cgl_lf_mks24_fig2a.py`. For this beta-10
panel, the pinned source and matching input deck establish pressure scale
`s = 0.5`; the extractor applies the corresponding two-dimensional PDF
Jacobian and emits comparison-manifest surfaces while retaining raw donor
samples for audit.
They also include the Figure 3 compressive-flow projection
`E_{khat dot u}(k_perp)`, a normalized density-fluctuation spectrum, and
separate `p_parallel`, `p_perp`, and AthenaK `B^2/2` spectra needed for
Figures 4(b) and 6(a). The products are rendered for inspection, but their
paper curves remain excluded until spectral ordinate transformations are
qualified. The manuscript definitions and cited donor plotted units do not
state the discrete Fourier normalization needed for absolute spectral or
strain ordinate conversion; those comparisons require external numeric data,
donor diagnostic code, or an explicit normalization statement. Figure 3's
sonic-correlation and beta-1 random definitions now
exist under the guarded `paper-compressive` workflow, but have not been run.
With repeated `--reference-curves <manifest.json>` options, `paper-analyze`
also accepts
external numerical or digitized curves and sampled joint-PDF surfaces only
when their manifest records provenance, SHA-256 digests, and positive
per-sample uncertainties. It reports uncertainty-normalized residuals
against supported PDF, spectrum,
raw or MKS24-normalized transfer, selected-shell alignment-distribution, alignment-peak-versus-
`k_perp`, eddy-anisotropy, threshold-volume history, and
`pressure_density_joint.parallel`/`.perpendicular` products and renders
curve or surface comparison figures. Use `--alignment-shells` with
`paper-analyze` when a comparison manifest requires an alignment-peak curve
over additional shells.
The staged-reference tooling includes pinned vector extraction for Figure
2(b), Figure 4(a) normalized-density PDFs, Figure 5(b) normalized eddy
anisotropy, Figure 7 lower-panel and Figure
13(d) normalized transfer curves, the
dimensionless Figure 9, Figure 11 lower-panel, and Figure 12 alignment
curves, and the dimensionless `beta Delta` PDF curves in Figure 13(b).
Figure 8 is additionally admitted as selected-shell `alignment.<shell>` PDFs:
its checksum-pinned RGB heatmaps are decoded using the labeled linear
colorbar, checked against the published per-shell unit normalization, and
emitted with an explicit raster-extraction uncertainty. The
paper states `p0 = 100` in code
units for its beta-100 limiter runs, while the equivalent AthenaK
`v_A = 1` normalization uses `p0 = 50`; until transforms for the listed
spectral and strain observables are qualified, dimensional Figure 11 upper
spectra, Figure 12
spectra, and Figure 13(a),(c) curves are retained only as excluded audit
context. Figure 2(a)'s sampled-surface comparison route is implemented, and
its labeled raster/color mapping is admitted through its source-derived
beta-10 pressure conversion and surface-density Jacobian.
Figure 5(b)'s dimensionless
curves are admitted through the opt-in
conditioned-structure-function analysis; no paper-scale comparison has been
executed.
These products are analysis infrastructure; they do not by themselves
establish statistically converged paper comparisons.

For user-facing validated runs and retained result summaries, follow the
documented workflows in [CGL Landau-Fluid Validation](cgl_landau_fluid_validation.md).
For developer regression execution, run
`python run_tests.py cgl/cgl_landau_fluid` from `tst/`. The routine CPU tests
compare the explicit split against STS capped with
`time/sts_max_dt_ratio=1.0`.
The focused shearing-box regression is:

```bash
cd tst
python run_test_suite.py --mpicpu \
  --test test_suite/cgl/test_cgl_lf_sbox_mpicpu.py
```

It compares serial and two-rank remote shearing partners and checks capped STS
against the one-stage explicit LF split. It also requires both LF heat-flux
channels to be active, verifies admissibility, checks magnetic fluxes and
`div B`, and exercises an MPI restart.
The routine CPU interaction regression uses
`inputs/tests/cgl_lf_turb_driving_amr.athinput` to check strict LF
admissibility, refinement, rendered forcing, and deterministic modal restart
while turbulence driving is active.

See also [CGL Method Physics Primer](cgl_mhd_method.md),
[CGL Landau-Fluid Code Guide](cgl_landau_fluid_code_guide.md),
[Super Time Stepping](super_time_stepping.md),
[Magnetohydrodynamics](mhd.md), and
[Turbulence Driving](turbulence_driving.md).
