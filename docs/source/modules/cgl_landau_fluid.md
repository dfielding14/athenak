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

The perpendicular coefficient uses the collisionally modified 3+1 closure of
[Snyder, Hammett & Dorland (1997), equation 49](https://w3.pppl.gov/~hammett/gyrofluid/papers/snyder_pop97.pdf),
also given by [Sharma, Hammett & Quataert (2003), equation 35](https://w3.pppl.gov/~hammett/gyrofluid/papers/2003/sharma03.pdf):

$$
\chi_\perp=\frac{c_\parallel^2}{\sqrt{\pi/2}\,c_\parallel k_L+\nu_{\rm eff}}
=\frac{2c_\parallel^2}{\sqrt{2\pi}\,c_\parallel k_L+2\nu_{\rm eff}}.
$$

Here $c_\parallel=\sqrt{p_\parallel/\rho}$ in local mode, and the background
collision frequency satisfies $\partial_t(p_\perp-p_\parallel)
=-\nu_{\rm coll}(p_\perp-p_\parallel)$ without other sources. This deliberately
differs from the single $\nu_{\rm eff}$ printed in the perpendicular denominator
of [Squire et al. (2023), equation 2.7](https://arxiv.org/html/2303.00468v2#S2.E7),
which uses the same pressure-relaxation convention. The collisionless coefficient
is unchanged, while $\chi_\perp\nu_{\rm eff}/c_\parallel^2\to1$ at strong
collisionality. The parallel coefficient and heat-flux caps are unchanged;
this choice does not make the complete closure an exact Braginskii transport model.

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

Collision rates are applied once per cycle with the full outer timestep.
Without split LF, rates and walls follow RK. If another STS process is enabled,
walls also follow each parabolic half-sweep, without repeating collision rates.
With split LF, the order is
`L(dt/2) W H(dt) W L(dt/2) C(dt) W`, where `W` applies walls and `C` applies
background and finite-rate limiter relaxation. A skipped LF sweep still receives
its wall/rate checkpoint. After regridding, walls use the final refreshed field
and the coarse state is refreshed. This collision schedule is not a symmetric
split and does not establish second-order accuracy for the complete model.

At fixed mean pressure, background relaxation multiplies
`Delta p = p_perp - p_parallel` by `exp(-nu_coll*dt)`. Beyond an enabled soft
threshold `Delta_*`, limiter relaxation then gives
`Delta_* + (Delta p - Delta_*)/(1 + limiter_nu_coll*dt)`.
The fluid-firehose wall `Delta p >= -B^2` is mandatory at wall checkpoints,
including configurations with no optional limiters. Explicitly enabled backup
walls add the configured firehose/mirror bounds. Energy, momentum and B remain
unchanged; corrected pressures are decoded from the stored E/A state.

LF histories count intermediate wall crossings without projecting each RKL
stage. Strict admissibility rejects floors and nonfinite/nonpositive state at
every stage, and rejects wall violations at entry/post-wall checkpoints.
Checkpoint checks do not add LF RHS counts. Backup scattering is enabled only
by `backup_limiters`; relaxed admissibility does not enable it implicitly.

## Closure Controls

| Parameter | Default | Meaning |
| --- | --- | --- |
| `eos` | required | Set to `cgl` for this feature. |
| `passive` | `false` | Use `true` for [passive J/A thermodynamics](cgl_passive.md) on native isothermal flow; requires positive `iso_sound_speed` and uniform periodic Newtonian MHD. |
| `cgl_heat_flux` | absent | Set to `landau_fluid` to enable LF transport. |
| `cgl_heat_flux_integrator` | `sts` | `sts` for production runs or `explicit` for reference verification. |
| `lf_k_parallel` | required | Finite positive closure wavenumber magnitude. |
| `lf_coefficient_mode` | `local` | `local` or `background`. |
| `lf_c_parallel0` | required in background mode | Finite positive fixed parallel thermal speed. |
| `nu_coll` | `0.0` | Finite nonnegative background anisotropy-relaxation frequency. |
| `mirror_limiter` | `false` | Enable mirror-limiter relaxation. |
| `firehose_limiter` | `false` | Enable firehose-limiter relaxation. |
| `firehose_threshold` | `2.0` | Finite positive coefficient: soft threshold is `-firehose_threshold*B^2/2`. |
| `mirror_threshold` | `1.0` | Finite positive coefficient: soft threshold is `mirror_threshold*B^2/2`. |
| `cgl_firehose_threshold` | absent | Legacy `oblique` maps to 1.4; `parallel` maps to 2.0. A conflicting numeric threshold is rejected. |
| `limiter_nu_coll` | required when a soft limiter is enabled | Finite nonnegative soft relaxation frequency; supplied values are checked even with limiters disabled. |
| `limiter_hardwall` | `false` | `true` is rejected. Use an explicit finite rate for soft relaxation. |
| `backup_limiters` | `false` | Enable both backup walls and their extra LF scattering; requires at least one soft limiter. |
| `firehose_backup_factor` | `1.0` | Finite factor >=1 multiplying the soft firehose threshold, capped by the mandatory fluid wall. |
| `mirror_backup_factor` | `2.0` | Finite factor >=1 multiplying the soft mirror threshold. |
| `limiter_backup_nu` | `1e10` | Finite nonnegative backup contribution to LF scattering. |
| `cgl_lf_strict_admissibility` | `false` | Reject invalid/floor-repaired stages and violations at wall checkpoints; count intermediate wall crossings. |
| `cgl_lf_record_pressure_work` | `false` | Retain RK-integrated applied CGL pressure-traction work diagnostics. |
| `cgl_lf_diagnostics` | `full` | `full` collects heat-flux face/cap/work diagnostics; `none` skips those reductions for production runs. |
| `cgl_lf_arithmetic` | `safe` | `safe` uses overflow-protected scaled arithmetic; `fast` uses direct normal-range `Real` arithmetic. |
| `cgl_lf_sts_flux` | `weighted` | `weighted` embeds STS weights in LF fluxes; `physical` writes physical fluxes and requires `cgl_lf_diagnostics = none` plus `cgl_lf_arithmetic = fast`. |
| `cgl_lf_profile` | `false` | Enable coarse LF timing regions with Kokkos fences; for profiling only. |
| `cgl_lf_profile_detail` | `false` | Compatibility key: `false` is accepted; the removed replay mode (`true`) is rejected. |

CGL requires finite positive `dfloor`, `pfloor`, `bfloor`, and `sigma_max`, and
finite nonnegative `tfloor` and `sfloor`. Zero collision rates remain valid.
These input guards do not qualify single-precision evolution or the numerical
repairs planned in later rebuild pieces.

Select runtime modes through the input parameters; LF mode environment overrides
have been removed. The active rebuild's production timing configuration is:

```ini
<mhd>
cgl_lf_diagnostics = none
cgl_lf_arithmetic = safe
cgl_lf_sts_flux = weighted
cgl_lf_strict_admissibility = true
cgl_lf_profile = false
```

Fast/physical mode remains available for separately labeled historical controls.

LF faces use the CT normal field and arithmetic means of the adjacent
cell-centered tangential components. They divide this vector by
$\bar B=(|\boldsymbol B_L|+|\boldsymbol B_R|)/2$, without renormalizing it
to unit length. The same $\bar B$ supplies inverse-field factors and the
magnetic pressure in the effective scattering rate. Faces with
$\bar B\le B_{\rm floor}$ contribute zero heat flux.

For neighboring fields $(g,+B,0)$ and $(g,-B,0)$, the normal thermal
conduction factor in the small-gradient limit is $g^2/(g^2+B^2)$.
Normalizing the canceled vector instead would incorrectly give unity.
This is a discrete regularization of unresolved reversals; a smooth sheet
of fixed physical width and nonzero guide field recovers the continuum
direction under refinement. The CT normal component is a separate staggered
value, so the resulting vector is not universally bounded by unit length.
This geometry alone does not establish a stability bound for the full stencil.

The default soft thresholds match `beta Delta = -2` and `+1`. Shipped decks
that formerly requested `limiter_hardwall=true` now explicitly use stiff
`limiter_nu_coll=1e10` relaxation where intended; existing explicit rates are
preserved. This is a once-per-cycle finite-rate update, not primitive-recovery
projection to the soft threshold. LF faces add background, applicable soft and
enabled backup frequencies once each; equality adds no scattering. Threshold
occupancy diagnostics may count equality, which is an admissible wall state.

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
analysis remains continuous across segments. `lf_hwproj` is retained for history
and restart compatibility; the scheduled collision/wall operator does not
increment this legacy primitive-refresh counter.
`lf_qprwrk` and `lf_qpewrk` are cumulative RKL2-weighted owned-face
contractions of the capped heat fluxes with their corresponding temperature
jumps. They characterize the closure-generated face fluxes before shearing-seam or
refined-mesh shared-face reconciliation. The cap counts likewise describe the
owned local closure estimate. These diagnostics are not exact contractions of
the reconciled fluxes applied at those faces. On refined meshes, each coarse/fine
interface uses the same fine-side ownership convention as the face-count
columns because those fluxes are restricted into the coarse update.
These are signed operator contractions; they are not required to be positive,
equal an offline snapshot proxy, or close a total energy budget. The existing
`aam-D` history column remains the conserved anisotropy variable for
compatibility.

When `cgl_lf_diagnostics = none`, LF stage-admissibility counters remain active,
but heat-flux face, cap, and q-work reductions are not
collected. In that mode the corresponding heat-flux diagnostic columns should
be treated as intentionally inactive rather than as measured zero cap
occupancy.

When `cgl_lf_record_pressure_work = true`, `lf_cpwrk` is the cumulative
explicit-RK-applied contraction of velocity with the retained CGL
pressure-traction divergence, and `lf_cawrk` is its `Delta p` anisotropic
component. The retained face traction is corrected through the same AMR flux
exchange used by the momentum update before the contraction is evaluated.
Passive-Delta runs retain zeros for both fields because their diagnostic CGL
pressures are not applied to flow momentum.

## Current Restrictions

- CGL is not available for SR, GR, or dynamical-GR MHD.
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
  `IAN = p_perp/|B|`. Halos are rebuilt in magnetic moment before the first LF
  evaluation and in anisotropy after split checkpoints, so interpolation and
  nonlinear slot conversion occur in a consistent order. Physical corners are
  refilled after prolongation. Same-level
  LF face estimates are averaged before either adjacent block updates, while
  coarse/fine flux restriction retains conservation. The timestep bound covers
  the local frozen stencil; it does not compose derivatives through prolongation
  and shared-face reconciliation. Finite convergence checks qualify tested AMR
  states, not a general nonlinear AMR stability envelope.
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
- Modal `<turb_driving>` forcing is supported with CGL LF. It applies one
  full-timestep kick at the final fluid stage, outside the protected LF
  magnetic-moment sweep. Existing boundary exchange and primitive recovery
  follow the kick; CGL primitive restriction retains its CT dependency.
  Coupling remains first-order split. The routine strict AMR/restart regression
  exercises this combined path.

## Verification

Focused unit problems in `inputs/unit_tests/` exercise CGL transforms, CGL
FOFC, LF parallel and perpendicular decay, magnetic-field-gradient coupling,
flux limiting, limiter suppression, and a field-aligned wave. Routine
regressions also exercise analytic uniform collisional relaxation and both
firehose threshold policies. The LF
quantitative pgen is the built-in `src/pgen/tests/cgl_landau_fluid.cpp`.

The reduced forced-turbulence initializer is registered as
`pgen_name = cgl_lf_paper`; the supported active smoke deck is
`inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput`. It
initializes `rho0 = 1`, `B0` along `z`, and
`p_parallel0 = p_perp0 = beta0 B0^2/2`, uses the explicit MKS24
`cgl_firehose_threshold = parallel` policy, and exercises the shared
turbulence driver. The passive smoke deck uses the separate
[J/A thermodynamic model](cgl_passive.md) and its required isothermal sound speed.
These are reduced smoke cases, not standard paper-resolution runs. Their
forcing-orientation, seed-continuation, and multi-cycle OU/full-kick work checks
qualify reduced mechanics, not paper-scale active/passive statistics or
figure diagnostics.

The paper pgen `.user.hst` output retains volume-integrated mass, kinetic,
magnetic and CGL thermal energies, `b2`, `b4`, `delta_p`, `abs_delta_p`,
local-beta integral, mirror/firehose/hard-bound volumes, effective collision
rate, and instantaneous forcing power (written using the compact history
labels `therm_cgl`, `abs_dp`, `mirror_vol`, `fire_vol`, `hard_vol`, and
`force_pwr`) in addition to forcing-orientation quantities. Paper inputs
also enable `record_injected_work`, adding cumulative exact net forcing-source
work, including the zero-net-momentum projection, as `force_work`. Each full
kick records its measured total-energy change at the final fluid stage for
active CGL, or its kinetic-energy change for passive CGL without changing J/A.
The analyzer uses it with conserved `tot-E` to report an active-Delta global
energy residual. These histories support reduced global
summaries such as `C_B2`; `.mhd.hst` supplies operator-face heat-flux-cap
activity and retained snapshots supply spatial diagnostics.

Paper-standard input definitions and the limiter-frequency scan live under
`inputs/cgl_lf_paper/`. They encode the standard `192x192x384` domain,
duration, hard-wall baseline, physical forcing shell, binary snapshot
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
executed.

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
