# CGL Landau-Fluid Heat Flux

## Scope

AthenaK supports non-relativistic anisotropic MHD with the CGL equation of
state and an optional Landau-fluid (LF) heat-flux closure. The feature is
implemented as a distinct MHD parabolic process, `mhd/cgl_heat_flux`, rather
than as ordinary isotropic thermal conduction.

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
cgl_firehose_threshold = oblique
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

When LF and CGL collisions or limiter scattering are active, each pre/post LF
half-sweep is followed by a collision source update over that same half-cycle
duration. Chronologically this is `L(dt/2) C(dt/2) H(dt) L(dt/2) C(dt/2)`,
where `L`, `C`, and `H` denote LF, collision, and hyperbolic updates. The two
collision calls therefore advance one physical `dt`, matching the single
full-step collision interval used when LF is disabled.

AthenaK records LF-stage health metrics in normal MHD history output whenever
this closure is active. Set `cgl_lf_strict_admissibility = true` for
verification runs to terminate immediately if an LF refresh produces
non-finite or non-positive thermodynamic state, activates density or pressure
floors, or crosses an emergency mirror/firehose bound. Emergency-bound
reporting is independent of whether `backup_limiters` is enabled.

## Closure Controls

| Parameter | Default | Meaning |
| --- | --- | --- |
| `eos` | required | Set to `cgl` for this feature. |
| `passive` | `false` | When `true`, evolve CGL/LF pressures diagnostically while mass, momentum, and magnetic-field fluxes use the isothermal-MHD passive-Delta path. |
| `iso_sound_speed` | required when `passive = true` | Isothermal sound speed used by passive-Delta momentum fluxes and signal speeds. |
| `cgl_heat_flux` | absent | Set to `landau_fluid` to enable LF transport. |
| `cgl_heat_flux_integrator` | `sts` | `sts` for production runs or `explicit` for reference verification. |
| `lf_k_parallel` | required | Positive closure wavenumber magnitude. |
| `lf_coefficient_mode` | `local` | `local` or `background`. |
| `lf_c_parallel0` | required in background mode | Positive fixed parallel thermal speed. |
| `nu_coll` | `0.0` | Background anisotropy-relaxation frequency. |
| `mirror_limiter` | `false` | Enable mirror-limiter relaxation. |
| `firehose_limiter` | `false` | Enable firehose-limiter relaxation. |
| `cgl_firehose_threshold` | `oblique` | `oblique` activates at `beta Delta <= -1.4`; `parallel` activates at `beta Delta <= -2`. |
| `limiter_nu_coll` | `0.0` | Limiter relaxation frequency. |
| `limiter_hardwall` | `false` | With an enabled instability limiter, replace its pressure-relaxation update by an energy-preserving projection to the selected mirror/firehose threshold. |
| `backup_limiters` | `false` | Apply rapid correction after an emergency bound is crossed. |
| `cgl_lf_strict_admissibility` | `false` | Fail an LF split stage on unsafe state, LF floors, or hard-bound violations. |
| `cgl_lf_record_pressure_work` | `false` | Retain RK-integrated applied CGL pressure-traction work diagnostics. |

At an operator face with `|B| <= bfloor`, LF does not construct a local field
direction and applies zero heat-flux contribution at that face. The local
accuracy campaign exercises this shutdown behavior with strict monitoring.

MKS24 production simulations use the mirror threshold `beta Delta >= 1` and
the parallel-firehose threshold `beta Delta <= -2`; paper reproduction inputs
must therefore set `cgl_firehose_threshold = parallel` and
`limiter_hardwall = true` explicitly for their hard-wall closure. In this
mode, CGL primitive recovery pins newly unstable pressures to the selected
threshold while preserving `0.5*p_parallel + p_perp`; ordinary `nu_coll`
relaxation remains enabled, but finite-rate limiter pressure relaxation is
replaced by the algebraic constraint. The
oblique-firehose policy remains the default to preserve behavior of existing
feature-branch inputs. In both policies `lf_mirror` and `lf_firehs` count
physical threshold occupancy, while `lf_hardbd` counts emergency numerical
overshoot and is a strict-validation failure.

Normal `.mhd.hst` output appends cumulative columns when LF is active:
`lf_nstage`, `lf_dfloor`, `lf_pfloor`, `lf_nonfin`, `lf_nonpos`,
`lf_mirror`, `lf_firehs`, `lf_hardbd`, `lf_qface`, `lf_qprcap`,
`lf_qpr10`, `lf_qpecap`, `lf_qpe10`, `lf_qprwrk`, `lf_qpewrk`, and
`lf_hwproj`. When
`cgl_lf_record_pressure_work = true`, it additionally appends `lf_cpwrk`
and `lf_cawrk`. The
face-count columns record evaluated LF faces and parallel/perpendicular
unlimited heat-flux ratios exceeding `q_max` or `10*q_max`. Differences
between successive rows give interval counts; normalize limiter counts by
`lf_nstage` and heat-flux-cap counts by `lf_qface`. All cumulative LF
diagnostic columns are preserved through CGL-LF restart files so interval
analysis remains continuous across segments. `lf_hwproj` counts applications
of the algebraic hard-wall constraint during CGL primitive-refresh task
ranges, including refreshed support/ghost states; it is not a normalized
active-cell occupancy statistic. It is expected to be nonzero in
limiter-active hard-wall production intervals and is not a safety violation.
`lf_qprwrk` and `lf_qpewrk`
are cumulative RKL2-applied
owned-face contractions of the capped heat fluxes with their corresponding
temperature jumps. On refined meshes, each coarse/fine interface is owned by
the fine-side closure faces whose flux is restricted into the coarse update.
These are signed operator contractions; they are not required to be positive,
equal an offline snapshot proxy, or close a total energy budget. The existing
`aam-D` history column remains the conserved anisotropy variable for
compatibility.

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
- CGL LF with mesh refinement currently requires conserved prolongation;
  `<mesh_refinement>/prolong_primitives = true` is rejected.
- Modal `<turb_driving>` forcing is supported with CGL LF. It is applied in
  the ordinary source-term task graph, outside the protected LF
  magnetic-moment sweep; the routine strict AMR/restart regression exercises
  this combined path.

## Verification

Focused unit problems in `inputs/unit_tests/` exercise CGL transforms, CGL
FOFC, LF parallel and perpendicular decay, magnetic-field-gradient coupling,
flux limiting, limiter suppression, and a field-aligned wave. Routine
regressions also exercise analytic uniform collisional relaxation and both
firehose threshold policies. The LF
quantitative pgen is the built-in `src/pgen/tests/cgl_landau_fluid.cpp`.

An active/passive-Delta reduced forced-turbulence initializer is registered as
`pgen_name = cgl_lf_paper`; its smoke decks are
`inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput` and
`inputs/cgl_lf_paper/cgl_lf_paper_smoke_passive_beta10.athinput`. They
initialize `rho0 = 1`, `B0` along `z`, and
`p_parallel0 = p_perp0 = beta0 B0^2/2`, use the explicit MKS24
`cgl_firehose_threshold = parallel` policy, and exercise the shared
turbulence driver. Passive mode requires `mhd/passive = true` and
`problem/passive_delta = true`; a routine regression verifies that changing
stable diagnostic initial anisotropy does not change its driven flow fields.
These are reduced smoke cases, not standard paper-resolution runs. Their
forcing-orientation, seed-continuation, and multi-cycle OU/RK source-work checks
qualify reduced mechanics, not paper-scale active/passive statistics or
figure diagnostics.

The paper pgen `.user.hst` output retains volume-integrated mass, kinetic,
magnetic and CGL thermal energies, `b2`, `b4`, `delta_p`, `abs_delta_p`,
local-beta integral, mirror/firehose/hard-bound volumes, effective collision
rate, and instantaneous forcing power (written using the compact history
labels `therm_cgl`, `abs_dp`, `mirror_vol`, `fire_vol`, `hard_vol`, and
`force_pwr`) in addition to forcing-orientation quantities. Paper inputs
also enable `record_injected_work`, adding cumulative exact net forcing-source
work, including the zero-net-momentum projection, as `force_work`. The counter
is advanced as an explicit-RK companion quantity so earlier stage increments
receive the same final-state weighting as the evolved conserved variables; the
analyzer uses it with conserved `tot-E` to report an active-Delta global
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
The routine CPU interaction regression uses
`inputs/tests/cgl_lf_turb_driving_amr.athinput` to check strict LF
admissibility, refinement, rendered forcing, and deterministic modal restart
while turbulence driving is active.

See also [Super Time Stepping](super_time_stepping.md) and
[Magnetohydrodynamics](mhd.md), and
[Turbulence Driving](turbulence_driving.md).
