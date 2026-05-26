# MKS24 Reproduction-First CGL-LF Manuscript and Later Extension Plan

Status: execution plan, revised 2026-05-26. The first paper-production
segment, mapped case `R16/s00_t0_t2`, completed on Frontier as job `4674731`
on 2026-05-25 and was recorded only as a clean partial diagnostic result
(`t = 1.66204848945`, `1.253333` node-hours), not an accepted scientific
result. After that legacy shared-MPI-I/O segment incurred an approximately
twenty-minute pause at its `t = 1` snapshot/checkpoint boundary, commit
`1f52784d` changed retained production snapshots and restarts to rank-local
files with grouped checksum/analysis support. A fresh ranked-output
`R16/s01_rankio_t0_t2` segment began running as Frontier job `4675322` on
2026-05-25; its complete initial eight-rank snapshot/restart groups were
successfully read by the rank-aware analyzer, and its complete `t = 1`
eight-rank snapshot/restart groups were emitted without the legacy-scale
boundary pause. It completed as an inspected clean partial at
`t = 1.87769834013` using `1.253056` node-hours. Same-digest continuation
`R16/s02_rankio_t1p877698_t2` (job `4676557`) reached `t = 2.0`, passed
inspection with zero strict LF safety counters and complete terminal
rank-local products, and used `0.106389` node-hours. It is an accepted prefix
gate, not a completed `t = 10` case. Same-digest continuation
`R16/s03_rankio_t2_t3p5` (job `4676696`) completed on 2026-05-25 at
16:50:22 EDT, reached exactly `t = 3.5`, passed formal inspection with zero
strict LF safety counters and complete terminal rank-local snapshot/restart
products, and used `1.383611` node-hours. It extends the accepted ranked
prefix only. Continuation `R16/s04_rankio_t3p5_t5` (job `4679803`)
completed on 2026-05-25 at 21:53:04 EDT, reached exactly `t = 5.0`, passed
formal inspection with zero strict LF safety counters and complete terminal
rank-local snapshot/restart products, and used `1.567778` node-hours.
Continuation `R16/s05_rankio_t5_t6p5` (job `4681476`) terminated cleanly
on its Athena wall-clock limit on 2026-05-26 at 00:11:57 EDT, reached
`t = 6.43828031751 < 6.5`, passed partial-continuation inspection with zero
strict LF safety counters and complete terminal ranked products, and used
`1.753333` node-hours. Short continuation `R16/s06_rankio_t6p438280_t6p5`
(job `4683291`) completed on 2026-05-26 at 00:34:30 EDT, reached exactly
`t = 6.5`, passed formal inspection with zero strict LF safety counters and
complete terminal ranked products, and used `0.080556` node-hours. It extends
the accepted ranked prefix to `t = 6.5`. Continuation
`R16/s07_rankio_t6p5_t7p5` (job `4683918`) completed on 2026-05-26 at
01:58:33 EDT, reached exactly `t = 7.5`, passed formal inspection with zero
strict LF safety counters and complete terminal ranked products, and used
`1.281944` node-hours. It extends the accepted ranked prefix to `t = 7.5`.
Continuation `R16/s08_rankio_t7p5_t8p5` was submitted as job `4686032`
from the accepted `s07` terminal restart set toward the first prefix
containing samples in the `t = 8`--`10` analysis window.
Utility commits `476f9dbd`, `9c5c1b40`, and `bad8ba05` select only the
explicit accepted restart lineage during bundling and admit normal
sampled-history separation across authenticated restart boundaries while
rejecting gaps inconsistent with the recorded cadence. They assembled the
ranked prefix through `t = 2.0` into diagnostic bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t2_prefix_20260525`, excluding
the abandoned shared-output branch. The default production-window
`paper-analyze` action selected no prefix snapshots because its declared
window is `t = 8`--`10` and wrote one history-only diagnostic figure; a
separately labeled transient-window (`t = 0`--`2`) analyzer check read all ten
retained snapshots, selected 101 history rows, passed its synthetic check,
and wrote nine diagnostic figures. Neither action is reproduction evidence.
The first scientific objective is an independent AthenaK reproduction
of the attainable numerical results in Majeski, Kunz, and Squire (2024,
hereafter MKS24), before any new weak-guide or background-collisionality
experiment is executed.

Primary manuscript source: `docs/cgl_lf_validation.tex`

Governing prose guide: `docs/writing_style_guide.md`

Detailed implementation and evidence record:
`docs/cgl_lf_mks24_reproduction_implementation_plan.md`

Pinned reference source: official arXiv `2405.02418v2`, staged by
`scripts/stage_cgl_lf_mks24_reference.py` with checksums retained outside the
tracked source tree.

Completion boundary: this document defines what must be done. It is not
evidence that MKS24 has been reproduced, that a new scientific result exists,
or that a manuscript is complete.

## 1. Revision Decision

The previous version of this plan did **not** include everything needed to
reproduce MKS24. It defined a new cubic-box experiment with weak and stronger
guide-field anchors, a scan in background collision frequency `nu_coll`, and
new target values for magnetic fluctuation amplitude and sonic Mach number.
Those are extension experiments, not the published MKS24 experiment.

MKS24 instead studies collisionless high-beta CGL-Landau-fluid (CGL-LF)
turbulence in an elongated mean-field box. Its collisionality study is a
variation of the anomalous microinstability scattering limiter
`nu_lim`, applied only after mirror or firehose thresholds are crossed, not a
general uniform background-collision scan. The manuscript must reproduce this
distinction exactly before using the extension to ask new questions.

The project is reorganized into two non-interchangeable stages:

| Stage | Objective | Permission boundary |
| --- | --- | --- |
| I. MKS24 reproduction | Reproduce the published CGL-LF simulation setup, numerical comparisons, and all quantitatively attainable figure claims; account honestly for the external hybrid-kinetic comparator. | This is the only production-science campaign authorized for planning under the current `4000` node-hour ceiling. Execution still requires the gates below. |
| II. New extension | Study weak versus stronger guide fields and background collisionality after Stage I establishes a credible baseline. | Retained as a later candidate protocol; it is not presently authorized and must be recosted and approved after Stage I. |

The manuscript claim discipline is:

| Claim level | Permitted statement | Required evidence |
| --- | --- | --- |
| Implementation | AthenaK implements the MKS24 CGL-LF case definitions and diagnostics used here. | Input/workflow tests, physical-equation audit, forcing and limiter audit, restart/GPU qualification, archived provenance. |
| CGL-LF numerical reproduction | AthenaK independently reproduces a specified MKS24 CGL-LF result within stated comparison uncertainty. | Matched production cases, explicit normalizations, reference data or admissible digitization, convergence/statistical checks, and a per-panel pass/fail table. |
| Full-paper reproduction | Every claimed MKS24 result, including non-AthenaK external comparisons, has been reproduced or independently re-evaluated. | In addition to CGL-LF reproduction, a qualified path for MKS24 Figure 10 and any other externally generated result. |
| Extension inference | A new guide-field or background-collisionality result extends the reproduced baseline. | Stage I acceptance followed by separately approved, funded, executed, and analyzed Stage II cases. |
| Not permitted | The method is perfect or reproduces results that have not passed the stated gates. | No numerical campaign establishes perfection; absent data are not evidence. |

## 2. MKS24 Source Contract

The source to reproduce is:

```text
S. Majeski, M. W. Kunz, and J. Squire (2024)
Self-organization in collisionless, high-beta turbulence
Journal of Plasma Physics 90, 535900601
arXiv:2405.02418v2
```

Before modifying reproduction inputs or interpreting a panel, future agents
must restage or verify the pinned source and reread the governing passages:

| Source portion | Reproduction relevance |
| --- | --- |
| `MKS24.tex:204-239` | CGL-MHD equations, 3+1 Landau-fluid heat fluxes, cold-electron and microinstability-model assumptions. |
| `MKS24.tex:463-469` | Domain, resolution, beta, standard closure parameters, forcing, and duration. |
| `MKS24.tex:478-496` | Fourier spectra, pressure-stress transfer, and scale-dependent alignment diagnostics. |
| `MKS24.tex:503-655` | Main numerical results through heat-flux sensitivity. |
| `MKS24.tex:662-698` | Microinstability limiter-rate scan and alternative firehose-threshold discussion. |
| `MKS24.tex:760-918` | Appendices relating the reduced theory to RKMHD, Braginskii stress, and larger density fluctuations. |

The implementation record is already substantially advanced: it defines
guarded production inputs, diagnostic products, reference extractors, known
normalization gaps, reduced tests, and GPU qualification evidence. This
manuscript plan does not supersede those technical records. It determines
which of them must be completed and accepted before writing scientific claims.

### 2.1 What "completely reproduce" can mean

MKS24 contains theory, CGL-LF simulations, and a hybrid-kinetic comparison
whose data originated outside the CGL-LF solver. An AthenaK CGL-LF campaign
can independently reproduce:

1. The stated CGL-LF closure, numerical regime, active/passive construction,
   forcing family, limiter-rate variation, and CGL-LF diagnostics.
2. Every CGL-LF simulation panel for which matching observables and a
   qualified normalization/reference comparison are available or obtained.
3. The paper's analytic scalings and their connection to the reproduced
   numerical observables, by a documented derivation audit.

It cannot silently claim independent reproduction of MKS24 Figure 10, which
uses a `Pegasus++` hybrid-kinetic simulation reported from prior work. A
full-paper reproduction claim therefore requires one of the following:

1. Obtain the underlying hybrid-kinetic dataset and reproduce the plotted
   alignment analysis with provenance.
2. Obtain and execute a scientifically equivalent hybrid-kinetic workflow,
   including its model assumptions and convergence evidence.
3. Restrict the claim explicitly to complete reproduction of MKS24's CGL-LF
   simulation results, treating Figure 10 as external comparison context.

Until that decision is resolved, the manuscript must use the third wording
and must not say that every result in MKS24 was independently reproduced.

## 3. Stage I Physical and Numerical Specification

### 3.1 CGL-LF model

Use the CGL pressure equations with the MKS24 3+1 Landau-fluid closure:

```math
q_{\perp}=-{v_{\rm th,\parallel}\over\sqrt{\pi}|k_\parallel|}
\left[\rho\nabla_\parallel\left({p_\perp\over\rho}\right)
-p_\perp\left(1-{p_\perp\over p_\parallel}\right)
{\nabla_\parallel B\over B}\right],
```

```math
q_{\parallel}=-{2v_{\rm th,\parallel}\over\sqrt{\pi}|k_\parallel|}
\rho\nabla_\parallel\left({p_\parallel\over\rho}\right),\qquad
v_{\rm th,\parallel}=\sqrt{2p_\parallel/\rho}.
```

Record AthenaK's unit conversion and coefficient mapping in every production
bundle. Preserve the paper assumptions of initially isotropic ion pressure,
cold electrons for the directly reproduced CGL-LF calculation, and
unresolved ion-Larmor-scale instabilities represented through the declared
limiter closure.

### 3.2 Published setup to match

| Quantity | MKS24 reproduction value |
| --- | --- |
| Domain | Fully periodic `[L_x, L_y, L_z] = [1, 1, 2]`, with `L_z = L_parallel` and `L_x = L_y = L_perp`. |
| Guide field | Uniform `B_0` along `z`; AthenaK decks state the selected unit mapping explicitly. |
| Standard mesh | `n_perp = 192`, `n_parallel = 384`. |
| Scale-separation meshes | Active Alfvenic beta-10 cases at `n_perp = 96, 192, 384`, with `n_parallel = 2 n_perp`. |
| Initial beta values | `beta0 = 1, 10, 100`. |
| Standard LF scale | `|k_parallel| = 4*pi/L_parallel = 2*pi` in the stated box. |
| Standard limiter | Hard-wall equivalent of `nu_lim = 10^10 v_A/L_perp`. |
| Instability thresholds | Mirror `beta Delta > 1`; parallel/fluid firehose `beta Delta < -2`. |
| Background collisions | `nu_coll = 0.0` for direct MKS24 CGL-LF reproduction. |
| Duration | At least `t_f = 10 L_perp/v_A`; statistical windows beyond `t v_A/L_perp = 6` and no shorter than `2 L_perp/v_A`. |

The paper reports a beta-100 pressure normalization in its Figure 13
discussion that does not by itself fix all AthenaK dimensional spectral
ordinate transforms. The input normalization must remain explicit, and
dimensional panel comparisons must be blocked until their mapping is
qualified rather than normalized by guesswork.

### 3.3 Published forcing to match

| Quantity | MKS24 reproduction value |
| --- | --- |
| Process | Ornstein-Uhlenbeck correlated velocity forcing. |
| Fiducial injection | `d_t E_K = 0.32 rho0 v_A^2 L_perp^3`. |
| Fiducial correlation time | `t_corr = L_parallel/v_A`. |
| Forced shell | `k in (2*pi/L_parallel) * [1, 3]`. |
| Power distribution | Proportional to `k^-2`. |
| Alfvenic forcing | Perturb only velocity perpendicular to `z`, incompressibly at the outer scale. |
| Random forcing | Velocity forcing without the Alfvenic directional constraint, while retaining declared time correlation. |
| Sonic-correlation case | The beta-100 random-forcing compressive test uses `t_corr` corresponding to `L_parallel/v_th`. |

If the MKS24 source does not publish a forcing random seed or an output needed
for pointwise identity, record that as unavailable source information. The
reproduction target then becomes statistical agreement of declared
observables, not equality of turbulent fields.

### 3.4 Collisionality terminology

For Stage I, distinguish two different inputs:

| Parameter | Meaning | Stage I role |
| --- | --- | --- |
| `nu_coll` | Uniform background collisional relaxation in AthenaK. | Fixed to `0.0`; a nonzero scan is not an MKS24 reproduction case. |
| `nu_lim` / AthenaK limiter-scattering representation | Anomalous scattering activated only in cells beyond mirror/firehose thresholds and driving pressure back toward marginality. | Reproduce MKS24's limiter sensitivity at beta 100 using `20`, `200`, and hard-wall-equivalent `10^10 v_A/L_perp`. |

The later uniform-`nu_coll` experiment remains useful, but it answers a new
question and must not be described as the reproduction of MKS24's
microinstability-limiter scan.

## 4. Stage I Run Matrix

The repository already defines guarded inputs for the required CGL-LF case
families. Reuse one accepted run wherever a plotted role is physically and
numerically identical; do not execute aliases twice merely because a figure
uses them twice.

### 4.1 Sixteen mapped CGL-LF simulations and one unresolved definition

| ID | Existing input role | Configuration | MKS24 purpose |
| --- | --- | --- | --- |
| R02 | `standard_active_alfvenic_beta10` | Active, Alfvenic, beta `10`, standard mesh/closure | Central backbone: Figures 1, 2, 4-9, 11, and 12 nominal role. |
| R03 | `standard_active_alfvenic_beta100` | Active, Alfvenic, beta `100`, hard wall | Figures 3-4 beta results and Figure 13 hard-wall role. |
| R04 | `standard_active_random_beta10` | Active, random, beta `10`, standard closure | Forcing comparison. |
| R05 | `standard_active_random_beta100` | Active, random, beta `100`, Alfvenic correlation time | Forcing/density comparison. |
| R06 | `standard_passive_alfvenic_beta10` | Passive, Alfvenic, beta `10` | Active/passive comparison. |
| R07 | `standard_passive_alfvenic_beta100` | Passive, Alfvenic, beta `100` | Unstable-volume/beta comparison and Figure 13 passive context. |
| R08 | `standard_passive_random_beta10` | Passive, random, beta `10` | Forcing/passive comparison. |
| R09 | `standard_passive_random_beta100` | Passive, random, beta `100` | Forcing/passive and Figure 3 beta comparison. |
| R10 | `compressive_active_random_beta1` | Active, random, beta `1` | Figure 3 beta-dependent compressive spectrum. |
| R11 | `compressive_active_random_beta100_sonic` | Active, random, beta `100`, sonic `t_corr` | Figure 3 sonic-correlation experiment. |
| R12 | `heat_flux_beta10_strong` | Active, Alfvenic, beta `10`, `|k_parallel|/100` | Figure 12 stronger-heat-flux limit. |
| R13 | `heat_flux_beta10_weak` | Active, Alfvenic, beta `10`, `100 |k_parallel|` | Figure 12 weak/double-adiabatic-like limit. |
| R14 | `nulim_beta100_20` | Active, Alfvenic, beta `100`, `nu_lim = 20` | Figure 13 limiter sensitivity. |
| R15 | `nulim_beta100_200` | Active, Alfvenic, beta `100`, `nu_lim = 200` | Figure 13 limiter sensitivity. |
| R16 | `scale_separation_beta10_nperp96` | Active, Alfvenic, beta `10`, `96 x 96 x 192` | Figure 11 low-resolution case. |
| R17 | `scale_separation_beta10_nperp384` | Active, Alfvenic, beta `10`, `384 x 384 x 768` | Figure 11 high-resolution case. |

The repository also defines `standard_active_alfvenic_beta1`. The current
checksum-qualified extractor mappings do not map that case to a displayed
MKS24 result: the source explicitly identifies Figure 3's beta-1 curve as
randomly driven (`R10`), while the admitted Figure 4 and Figure 9 mappings use
`R02`-`R09`. Treat the separate Alfvenic beta-1 deck as unresolved inventory,
not as an authorized production case, unless the Phase A source audit
identifies and records a specific published claim requiring it.

`R03` supplies the otherwise duplicated beta-100 hard-wall `nu_lim` role, and
`R02` supplies the `n_perp = 192` scale-separation and nominal heat-flux roles.
The currently mapped run matrix is therefore sixteen unique calculations
after role and input-alias reuse.

### 4.2 Alias and execution contract

The workflow exposes some published roles through more than one campaign
entry. The production manifest must resolve these aliases before reserving
compute:

| Physical case | Existing workflow role(s) | Required handling |
| --- | --- | --- |
| `R02` | `paper_standard_active_alfvenic_beta10`; reused by the nominal Figure 12 comparison and the `n_perp = 192` Figure 11 comparison. | Execute once and link all three panel roles to the same accepted run bundle. |
| `R03` | `paper_standard_active_alfvenic_beta100` and `paper_nulim_beta100_hardwall`. | Compare parsed physics, numerical, forcing, output, and analysis parameters after ignoring the basename/comments; if they remain identical, execute once and link the hard-wall role to that output. |

If an intended alias has a material parameter difference, it is not an alias:
assign a new `Rxx` identifier, revise the run count and budget, and document
why the additional calculation is necessary before submission. In particular,
the five paper workflow commands must not be executed independently without
deduplicating their overlapping case roles.

### 4.3 Reproduction is not achieved by input files alone

The existing input definitions are necessary but not sufficient. Before a
case contributes to a claim, it must have:

1. A retained submitted input, executable hash, code revision, module/build
   provenance, forcing metadata, and case identity.
2. Strict numerical safety and restart evidence appropriate to its limiter
   model and resolution.
3. A validated late-time analysis interval satisfying the published duration
   and averaging requirements.
4. Retained source snapshots/histories, derived products, analysis
   configuration, checksums, and cost/storage accounting.
5. A comparison outcome for every MKS24 panel to which it contributes.

## 5. Figure-by-Figure Reproduction Gate

The paper shall not state "MKS24 is reproduced" until the following table has
been completed with pass/fail/blocked outcomes and justified uncertainties.

| MKS24 figure or result | Required cases/products | Reproduction gate |
| --- | --- | --- |
| Figure 1: beta-10 active/passive `beta Delta` slices | `R02`, `R06`; common-time or statistically declared slice visual; instability mask. | Recreate the qualitative active/passive spatial contrast with source-provenance caption; do not treat a hand-selected attractive snapshot as quantitative proof. |
| Figure 2(a): pressure-density joint PDFs | `R02`; joint `delta p_parallel`, `delta p_perp`, density surface. | Use the existing source-qualified sampled-surface transform or a superior author dataset; report residuals and uncertainty. |
| Figure 2(b): unstable-volume histories | Active/passive Alfvenic/random beta `10`, `100` roles (`R02`-`R09`). | Compare all eight curves against the retained checksum-qualified reference histories and the reported active/passive volume-fraction conclusion. |
| Figure 3: compressive velocity spectra | Panel (a): `R03`, `R05`, `R11`; panel (b): `R10`, `R04`, `R11`, `R09`. | Qualify absolute spectral normalization or obtain reference data; then test forcing-correlation and beta-dependence conclusions. |
| Figure 4: density PDF and density spectra | The eight Alfvenically correlated active/passive beta-10/beta-100 roles, `R02`-`R09`. | Compare admitted normalized density PDFs; unblock spectral ordinate mapping before claiming reproduction of panel (b). |
| Figure 5: energy spectra and local-field eddy anisotropy | `R02`, `R04`; local-field structure/eddy products. | Compare admitted normalized eddy-scale curves; qualify panel (a) spectral ordinate mapping. |
| Figure 6: pressure/magnetic and strain spectra | `R02`, `R04`. | Establish comparison normalization for dimensional spectra and test pressure-balance/parallel-strain-suppression claims. |
| Figure 7: pressure-anisotropy gradient and transfer | `R02`, `R04`, `R06`; normalized transfer product. | Compare admitted dimensionless transfer curves; unblock upper-panel dimensional spectra before a complete panel claim. |
| Figure 8: scale-dependent alignment PDFs | `R02`, `R06`. | Compare existing calibrated selected-shell alignment PDFs and preserve color-map uncertainty. |
| Figure 9: alignment peaks across cases | The eight standard active/passive beta-10/beta-100 roles, `R02`-`R09`. | Compare the existing admitted peak-alignment curves and test active/passive separation. |
| Figure 10: hybrid-kinetic comparison | External `Pegasus++` beta-16 result, not generated by CGL-LF AthenaK. | Block a full-paper reproduction claim until qualified source data or a compatible kinetic rerun is obtained; otherwise label it external context. |
| Figure 11: resolution/scale separation | `R16`, `R02`, `R17`. | Compare dimensionless alignment peaks; unblock energy-spectrum ordinate mapping; verify the scale at which alignment ends trends with dissipation scale. |
| Figure 12: heat-flux sensitivity | `R12`, `R02`, `R13`, plus nominal passive `R06`. | Compare dimensionless alignment peaks; unblock lower-panel gradient-spectrum mapping; test claimed insensitivity across factor-100 `|k_parallel|` variations. |
| Figure 13: limiter-induced collisionality | `R14`, `R15`, `R03`, with passive hard-wall context `R07` where plotted. | Compare admitted `beta Delta` PDF and normalized transfer curves; unblock energy/strain dimensional panels; test nonmonotonic response to `nu_lim`. |
| Figure 13 associated quantitative text | `R15` postprocessing with parallel and oblique firehose masks; outer-scale rate estimate; external kinetic context. | Test the reported `10.4%` versus `18.5%` unstable-volume reclassification over the same late-time interval and the estimated rate near `235 v_A/L_perp`; treat the reported hybrid-kinetic `17.9%` comparison as external unless its data are obtained. |
| Analytic/reduced-model conclusions and appendices | Equation audit plus reproduced diagnostics above. | Trace each claimed analytic prediction to an equation and a reproduced numerical test; explicitly separate derivation verification from independent simulation evidence. |

### 5.1 Theory and appendix audit

The numerical panel gate alone does not reproduce all of MKS24. The
derivations and stated analytic limits must be checked independently:

| MKS24 theory component | Required reproduction work | Acceptance record |
| --- | --- | --- |
| Main-text high-beta reduced CGL-MHD ordering and invariants | Re-derive the retained ordering, pressure-stress suppression, heat-flux scaling, and conserved quantities from the stated closure and conventions. | Equation-by-equation derivation notebook or manuscript appendix, with assumptions and any discrepancies listed. |
| Appendix: comparison with reduced kinetic MHD | Reproduce the reduction and identify exactly which magneto-immutability signatures survive and which non-local pressure effects do not. | Checked derivation linked to the main-text claim table. |
| Appendix: Braginskii viscous stress | Reproduce the collisional-limit argument and the criterion `beta k_parallel v_A/nu >> 1`, keeping it distinct from the threshold-activated `nu_lim` simulations. | Checked derivation plus terminology review against Stage I Section 3.4. |
| Appendix: larger density fluctuations | Reproduce the amended ordering and identify which conclusions depend on the small-density-fluctuation assumption. | Checked limitation statement connected to reproduced density diagnostics. |

The staged source audit identifies thirteen displayed result figures and
these three theoretical appendices. A complete claim must account for both
the panel gate above and this derivation audit.

### 5.2 Blocking reference-data work

The implementation record currently documents a real comparison boundary:
dimensionless products and a pressure-coordinate surface transform are
admitted, while several absolute spectral and strain ordinates are not yet
qualified because the source does not specify the discrete transform
normalization required for a faithful mapping to AthenaK products.

For a complete CGL-LF reproduction, those panels cannot remain silently
excluded. Before scientific manuscript completion, do at least one of:

1. Obtain machine-readable panel data or analysis code from the authors or an
   archival deposit.
2. Obtain an explicit normalization convention sufficient to implement and
   validate the transform.
3. Reanalyze source simulation dumps if made available with their analysis
   definitions.

If none is available, the manuscript may still report a carefully delimited
partial reproduction, but it must list the blocked panels and may not use
"complete reproduction."

## 6. Diagnostics, Visuals, and Acceptance Rules

### 6.1 Required MKS24 observables

| Observable family | Minimum retained products |
| --- | --- |
| Basic state | `rho`, `u`, `B`, `p_parallel`, `p_perp`, beta, `Delta p`, forcing work, and energy accounting. |
| Instability exposure | Mirror/firehose active fractions, limiter action, `nu_lim` policy, hard-wall projection/counters, and any invalid-state counters. |
| Pressure/density | Joint pressure-density PDFs, normalized density PDFs, and density fluctuation spectra. |
| Compressibility | Fourier-projected compressive velocity spectrum and its fraction of flow energy. |
| Cascade/pressure structure | Kinetic, magnetic, `p_parallel`, `p_perp`, magnetic-pressure, `Delta p`, and local parallel/perpendicular gradient spectra. |
| Energetics | `Delta p`-stress transfer with the MKS24 stated normalization, applied pressure/LF work, forcing work, and global accounting quality. |
| Organization | Scale-dependent rate-of-strain/local-field alignment PDFs and peak curves; local-field structure-function eddy anisotropy. |
| Reproduction comparison | Per-panel reference data, normalization metadata, residuals, uncertainty, and pass/fail/blocked status. |

### 6.2 Additional credibility diagnostics

These diagnostics strengthen the reproduction without changing the target
experiment:

| Addition | Reason |
| --- | --- |
| Low-field fractions `f_low(q)` and field-direction regularization counters | Protect local-field operations from hidden near-zero-field failures. |
| Block/bootstrap uncertainty over the declared late-time interval | Prevent one realization or burst from determining a result. |
| Exact forcing/pressure/LF work ledgers and restart identity | Distinguish physical stress transfer from implementation/accounting error. |
| Resolution and output-cadence sensitivity for selected products | Demonstrate that the plotted conclusion is not created by numerical or sampling choices. |

### 6.3 Visual-selection rule

For reproduction figures, first use the panel definition in MKS24. For any
additional AthenaK visual:

1. Select snapshots by a predeclared statistical rule, such as the saved
   output closest to the midpoint of the accepted averaging window.
2. Use fixed or pooled-percentile color limits documented in provenance.
3. Label supplementary event-focused snapshots as selected by an explicit
   criterion; do not substitute them for statistics.

### 6.4 Acceptance categories

Each panel/product must be assigned one of:

| Status | Meaning |
| --- | --- |
| `passed` | Setup and observable are matched, comparison uncertainty is declared, and AthenaK is quantitatively consistent with the reference criterion. |
| `failed` | A qualified comparison is inconsistent; retain it and revise the claim. |
| `blocked_reference` | A matching AthenaK product exists, but a defensible published-data transform/reference dataset is unavailable. |
| `external_model` | The published panel uses physics/code outside the AthenaK CGL-LF reproduction scope, such as Figure 10. |
| `not_run` | The required production simulation has not been executed. |

## 7. Stage I Implementation Work Before Production

The detailed record should remain the source of truth for completed
implementation evidence. For manuscript execution, the immediate required
work is:

| Work item | Existing starting point | Required completion evidence |
| --- | --- | --- |
| Freeze MKS24 case inventory | `inputs/cgl_lf_paper/mks24_stage_i_manifest.json` and guarded workflow `paper-mks24-stage-i` define sixteen source-mapped executions, reuse aliases, and exclude the unmapped active-Alfvenic beta-1 deck. | Preserve the committed manifest and validate it against every submitted input bundle. |
| Close reference-data boundary | Many dimensionless curves/surfaces are extracted; dimensional panels remain blocked. | Author/archive data, explicit transform proof, or manuscript-scoped blocked-panel decision. |
| Verify numerical model identity | Corrected collisions, limiter threshold selection, forcing work, passive feedback, and reduced/GPU gates are recorded. | Formal pre-production review against Stage I Sections 3-6 and passing required tests. |
| Production accounting path | `scripts/frontier/cgl_lf_stage_i.py` is separate from the debug-only utility, validates mapped cases and aliases, uses `batch` with default production `normal` QOS, archives executable/input and production-utility provenance, enforces sequential Stage I reservations, and requires an inspection record before an output may be recorded as `accepted` or retained as a clean partial prefix. Rank-local output sets are retained and checksum-verified as grouped products. Continuations may use only the inspected terminal restart set from a parent with matching case, input digest, and executable digest. | Record and review each completed segment before preparing its continuation or the next case. |
| Production analysis/plot orchestration | Existing analyzer and MKS24 extractors implement many products; `paper-analyze` composes repeated checksum-qualified split reference manifests while rejecting duplicate product identifiers, and commit `449e297b` adds explicit partial-case comparison recording plus header-only snapshot-window selection; `cgl_lf_stage_i.py bundle-case` follows the explicit accepted restart lineage, merges sampled histories that need not repeat authenticated restart boundary rows while bounding gaps by their retained cadence and boundary timesteps, and links time-deduplicated shared or rank-local snapshot sets only after a case reaches its required final time, while `bundle-campaign` requires all sixteen completed mapped cases for cross-case comparisons. | Workflow command regenerating the per-panel gate table and retained figures from accepted output bundles. |
| Manuscript conversion | Illustrated validation note exists. | TeX manuscript structured around reproduction claims, blocked boundaries, and only later an extension. |

No weak/strong-guide case family or uniform-`nu_coll` extension case should be
implemented as a production priority while these reproduction blockers remain.

## 8. Stage I Compute and Storage Budget

The user-raised ceiling is `4000` Frontier node-hours for the broader project.
Stage I is prioritized within that ceiling. The detailed implementation record
provides the current conservative MKS24 production estimate from corrected
reduced-GPU evidence:

Any `1000` node-hour ceiling retained in that earlier implementation record
predates this revision; the `4000` node-hour project ceiling and Stage I-first
reservation in this document control future submission authorization.

| Stage I allocation | Calculation basis | Node-hours |
| --- | --- | ---: |
| Previously consumed qualification ledger | Retained CGL-LF debug qualification evidence | `0.851670` |
| Conservative mapped-run reservation including unresolved-definition contingency | Earlier seventeen-run estimate: fourteen mapped standard-layout cases plus one held beta-1 contingency slot, one `n_perp = 96`, and one `n_perp = 384` case | `534.444444` |
| Runtime contingency for measurement error, segment extensions, and reruns | `2x` planning ceiling retaining the unresolved-definition contingency until Phase A disposition | `534.444445` |
| Stage I production reservation plus prior consumption | `0.851670 + 1068.888889` | `1069.740559` |
| Remaining project ceiling after Stage I reservation | `4000 - 1069.740559` | `2930.259441` |

The estimate is not a measured standard-run benchmark. Prior to production:

1. Implement or approve the production submission/accounting path.
2. Execute no more than one production case at a time.
3. Segment runs at restart boundaries and review strict counters, work
   accounting, storage, and measured cost after the first segment.
4. Stop and recost if projected Stage I consumption exceeds its reserved
   ceiling.

The existing estimated sequential-retention storage allocation, conservatively
including the unresolved beta-1 definition as a contingency, is approximately
`880.368 GB` with margin; the no-pruning alternative is approximately
`1140.924 GB`. Confirm the storage reservation and retention rule before
executing the mapped matrix, and reduce these envelopes only after Phase A
records the beta-1 disposition and updates the cost model.

### 8.1 Active production record

Stage I execution has begun under the production path; the following rows
are execution provenance, not accepted reproduction evidence:

| Case/segment | Frontier job | Requested allocation | Input/build provenance | Acceptance gate |
| --- | ---: | --- | --- | --- |
| `R16/s00_t0_t2` | `4674731`, `COMPLETED`, `0:0` | `1` node, `01:15:12` charged elapsed; `1.253333` actual node-hours | Legacy shared-MPI-I/O input at committed revision `e129a7ff8be60eb177c8268da0a73e7700d8b543`; executable SHA-256 `de99066f3546976093f6ff656660952032e30f31c939736c0212629cfb015157`; archived submitted input and matrix in the run manifest directory. | `clean_partial`: terminal `t = 1.66204848945 < 2.0`, zero strict LF safety counters and retained terminal restart; retained as diagnostic evidence only and not continued into the changed output-protocol digest. |
| `R16/s01_rankio_t0_t2` | `4675322`, `COMPLETED`, `0:0` | `1` node, `01:15:11` charged elapsed; `1.253056` actual node-hours | Rank-local retained-output input/utility revision `1f52784d291718fff3bf1dc984a682ada52780ab`; executable revision `e129a7ff8be60eb177c8268da0a73e7700d8b543` and SHA-256 `de99066f3546976093f6ff656660952032e30f31c939736c0212629cfb015157`; complete `t = 1` ranked snapshot/restart boundary and terminal ranked products retained, initial snapshot read through the rank-aware analyzer, and no legacy-scale output pause observed. | `clean_partial`: terminal `t = 1.87769834013 < 2.0`, zero strict LF safety counters and complete terminal rank-local snapshot/restart products; eligible for same-digest continuation only. |
| `R16/s02_rankio_t1p877698_t2` | `4676557`, `COMPLETED`, `0:0` | `1` node, `00:06:23` charged elapsed; `0.106389` actual node-hours | Continued from inspected terminal eight-rank restart set of `s01` with unchanged input/executable digests. | `accepted`: reached `t = 2.0`, zero strict LF safety counters, and complete terminal rank-local snapshot/restart products; establishes the ranked prefix gate only. |
| `R16/s03_rankio_t2_t3p5` | `4676696`, `COMPLETED`, `0:0` | `1` node, `01:23:01` charged elapsed; `1.383611` actual node-hours | Continued from inspected accepted terminal eight-rank restart set of `s02` with unchanged input/executable digests; terminal outputs contain five complete eight-rank snapshot groups at `t = 2.50014742786`, `2.75033411668`, `3.00023812685`, `3.25023774201`, and `3.5`, plus the terminal eight-rank restart group. | `accepted`: reached `t = 3.5`, formal inspection records zero strict LF safety counters and complete terminal rank-local products; establishes a longer ranked prefix gate only. |
| `R16/s04_rankio_t3p5_t5` | `4679803`, `COMPLETED`, `0:0` | `1` node, `01:34:04` charged elapsed; `1.567778` actual node-hours | Continued from the inspected accepted terminal eight-rank `.00004.rst` set of `s03` with unchanged input/executable digests. Terminal outputs contain six complete eight-rank snapshot groups at `t = 3.75032903725`, `4.00000477772`, `4.25028947444`, `4.49999984569`, `4.75028920409`, and `5.0`, plus the terminal eight-rank restart group. | `accepted`: reached `t = 5.0`, formal inspection records zero strict LF safety counters and complete terminal rank-local products; establishes a longer ranked prefix gate only. |
| `R16/s05_rankio_t5_t6p5` | `4681476`, `COMPLETED`, `0:0` | `1` node, `01:45:12` charged elapsed; `1.753333` actual node-hours | Continued from the inspected accepted terminal eight-rank `.00005.rst` set of `s04` using unchanged input/executable digests. Terminal outputs contain six complete eight-rank snapshot groups at the `t = 5.25`, `5.5`, `5.75`, `6.0`, `6.25`, and terminal `6.43828031751` states, plus complete scheduled `t = 6.0` and terminal restart groups. | `clean_partial`: clean application-walltime termination at `t = 6.43828031751 < 6.5`; formal inspection records zero strict LF safety counters and complete terminal ranked products, allowing same-digest continuation only. |
| `R16/s06_rankio_t6p438280_t6p5` | `4683291`, `COMPLETED`, `0:0` | `1` node, `00:04:50` charged elapsed; `0.080556` actual node-hours | Continued from the inspected clean-partial terminal eight-rank `.00007.rst` set of `s05` using unchanged input/executable digests and a `time/tlim=6.5` override; terminal outputs include complete eight-rank snapshot and restart groups at `t = 6.5`. | `accepted`: reached `t = 6.5`, formal inspection records zero strict LF safety counters and complete terminal rank-local products; promotes the accepted ranked prefix through the clean-partial parent. |
| `R16/s07_rankio_t6p5_t7p5` | `4683918`, `COMPLETED`, `0:0` | `1` node, `01:16:55` charged elapsed; `1.281944` actual node-hours | Continued from the inspected accepted terminal eight-rank `.00008.rst` set of `s06` using unchanged input/executable digests and a `time/tlim=7.5` override; terminal outputs contain three complete eight-rank snapshot groups at `t = 7.00026203265`, `7.25019568336`, and `7.5`, plus the terminal eight-rank restart group. | `accepted`: reached `t = 7.5`, formal inspection records zero strict LF safety counters and complete terminal rank-local products; extends the accepted ranked prefix but remains before the analysis window. |
| `R16/s08_rankio_t7p5_t8p5` | `4686032`, `PENDING` at submission check on `2026-05-26` | `1` node, `02:00:00`, `batch` / default `normal` QOS; `2.000000` reserved node-hours | Submitted from the inspected accepted terminal eight-rank `.00009.rst` set of `s07` using unchanged input/executable digests and a `time/tlim=8.5` override. | Require exact target completion, formal inspection, zero strict LF safety counters, and complete terminal rank-local products before admitting the first production-window prefix evidence. |

The ranked `t = 0`--`2` prefix was assembled as the diagnostic bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t2_prefix_20260525` after the
bundler was made lineage-aware, then reassembled successfully under the
cadence-bounded implementation as
`runs/bundles/diagnostic-prefixes/R16_rankio_t2_prefix_cadence_20260525`.
Both bundles contain only `s01` and `s02` and link ten retained rank-local
snapshots. The first bundle's default `paper-analyze`
execution produced only a history summary because a `t = 2` prefix has no
samples in the production `t = 8`--`10` analysis window. A distinct,
infrastructure-only transient-window check over `t = 0`--`2` then read all
ten snapshots, selected 101 history rows, passed the synthetic analyzer
check, and rendered nine diagnostic figures. It is not a final-time `R16`
bundle and is not panel reproduction evidence.

After `s03` acceptance, the lineage-checked bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t3p5_prefix_cadence_20260525`
was assembled from only `s01`, `s02`, and `s03`; its manifest records
`accepted_final_time = 3.5`, `status = accepted_for_analysis`, and fifteen
deduplicated rank-local snapshots. It verifies accepted-prefix assembly at
the new gate, but it still has no samples in the required production window
`t = 8`--`10`: a default-window analyzer check found `176` merged history
rows overall but selected zero history rows and zero snapshots in that
window, while passing its synthetic diagnostic test. It is not panel
reproduction evidence.

After `s04` acceptance, the further lineage-checked bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t5_prefix_cadence_20260525`
was assembled from only `s01` through `s04`; its manifest records
`accepted_final_time = 5.0`, `status = accepted_for_analysis`, and twenty-one
deduplicated rank-local snapshots. A default-window analyzer check found
`251` merged history rows overall but selected zero history rows and zero
snapshots in `t = 8`--`10`, while passing its synthetic diagnostic test.
This is accepted continuation evidence, not panel reproduction evidence.

After `s06` acceptance, including the inspected `s05` clean partial in its
continued lineage, the accepted-prefix bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t6p5_prefix_cadence_20260526`
was assembled from `s01` through `s06`; its manifest records
`accepted_final_time = 6.5`, `status = accepted_for_analysis`, and
twenty-eight deduplicated ranked snapshots. A default-window analyzer check
found `326` merged history rows overall but selected zero history rows and
zero snapshots in `t = 8`--`10`, while passing its synthetic diagnostic test.
This is accepted continuation evidence, not panel reproduction evidence.

After `s07` acceptance, the accepted-prefix bundle
`runs/bundles/diagnostic-prefixes/R16_rankio_t7p5_prefix_cadence_20260526`
was assembled from `s01` through `s07`; its manifest records
`accepted_final_time = 7.5`, `status = accepted_for_analysis`, and
thirty-one deduplicated ranked snapshots. A generic default-window analysis
found `376` merged history rows overall but selected zero history rows and
zero snapshots in `t = 8`--`10`, while passing its synthetic diagnostic
test and writing one history-only figure. A requested Figure 11 reference
comparison fails closed because the present case selects an empty analysis;
this prefix is not panel reproduction evidence.

Commit `449e297b` adds an explicitly scoped
`--allow-partial-reference-cases` comparison mode, which records products
omitted because their cases are not present in a single-case bundle, and
selects snapshot windows from binary headers before loading full fields. A
local terminal-snapshot plumbing probe at `t = 2.0` evaluated only the
available Figure 11 `R16` product and recorded the absent `R02` and `R17`
products as omissions. That transient, out-of-window probe is not production
comparison evidence.

Job `4683918` is inspected and recorded as an accepted continuation through
`t = 7.5` after the clean partial `4681476` at `t = 6.43828031751` and
accepted terminal-gate segment `4683291`; job `4686032` is submitted toward
the `t = 8.5` acceptance gate and the first production-window prefix.
Do not describe `R16` as a completed accepted case until ranked-output
segments satisfy the stated gates through the required final time.

### 8.2 Why Stage II is not currently reserved in full

The former full extension plan reserved `3735.627090` node-hours including
the existing ledger. Its incremental requirement after the already consumed
ledger is `3734.775420` node-hours. Combining that extension matrix with the
Stage I MKS24 reservation would require:

```text
0.851670 + 1068.888889 + 3734.775420 = 4804.515979 node-hours,
```

which exceeds the `4000` node-hour ceiling by `804.515979` node-hours. The
full extension matrix is therefore not simultaneously authorized. After
Stage I, use measured costs and results to design a reduced Stage II within
the remaining budget or request a larger allocation.

## 9. Stage I Manuscript Program

`docs/cgl_lf_validation.tex` presently documents method validation. Convert
it into the reproduction manuscript only after the production inputs,
reference-data decision, and execution authorization have been reviewed.

| Manuscript section | Required argument and evidence |
| --- | --- |
| Abstract | State the independent reproduction question, the exact CGL-LF scope achieved, principal quantitative outcomes once available, and any external/normalization limitation. |
| Introduction | Motivate magneto-immutability and the need to reproduce MKS24 before testing new guide-field/collisionality regimes. |
| MKS24 target contract | State the closure, box, forcing, beta/limiter/LF-scale cases, and what direct reproduction excludes unless externally supplied. |
| AthenaK method validation | Condense existing operator, forcing, limiter, restart, GPU, and work-accounting evidence needed to trust the independent reproduction. |
| Reproduction campaign | Present the sixteen mapped cases, any audited beta-1 disposition, output/analysis protocol, cost/storage provenance, and per-panel comparison method. |
| Reproduced equation of state and compressive behavior | Figures 1-4 comparisons, with blocked status if required normalization remains unresolved. |
| Reproduced cascade and self-organization | Figures 5-9 and 11 comparisons: spectra, transfer, alignment, structure, and resolution. |
| Reproduced closure sensitivities | Figures 12-13: heat-flux scale and limiter-induced scattering rate. |
| External comparison boundary | Treat Figure 10 and any unregenerated hybrid-kinetic content explicitly. |
| Discussion and conclusions | State which MKS24 CGL-LF results passed, failed, or remain blocked; introduce Stage II only as subsequent work unless it has later been executed. |
| Reproducibility appendix/data statement | Inputs, hashes, source checksums, extraction/analysis configuration, comparison table, ledger, and retained-data locations. |

The main reproduction visuals should follow the published panel logic rather
than replacing it with a new narrative. Additional numerical-control and
provenance figures may be placed in appendices or supplemental material.

### 9.1 Manuscript acceptance checklist

A manuscript may claim completion of the Stage I reproduction only when:

1. It follows `docs/writing_style_guide.md` and makes no perfection claim.
2. Every Stage I run used for a claim is accepted, archived, and linked to its
   case role, analysis interval, provenance, and budget entry.
3. Every MKS24 simulation figure/result is assigned a status from Section 6.4.
4. All status `passed` statements have quantitative comparison data and
   uncertainty; all blocked or external panels are disclosed.
5. A full-paper reproduction claim is made only if the Figure 10 external
   model gate and any other external-result gates have been resolved.
6. TeX, bibliography, compact plotted source data, rendered figures, and
   provenance are tracked; large production data are retained externally with
   checksums.
7. The PDF builds through the declared TeX/BibTeX workflow without missing
   figures or references.

## 10. Stage II: Later Guide-Field and Background-Collisionality Extension

Stage II preserves the new experiment requested before the MKS24 priority was
identified. It is not part of Stage I reproduction and must not begin until
the Stage I acceptance decision is recorded.

### 10.1 New questions

After reproducing MKS24, ask:

1. How do magneto-immutability, anisotropic pressure work, and structural
   diagnostics change when field wandering is much larger than in MKS24?
2. How does a uniform background collisional relaxation rate `nu_coll`
   change the turbulent state, separately from MKS24's threshold-activated
   anomalous limiter scattering?
3. Are guide-field-dependent responses robust across beta and LF closure
   scale?

### 10.2 Extension anchors

Use a periodic cubic box and isotropic solenoidal driving only as an
extension design, not a reproduction configuration. Define:

| Extension guide regime | Collisionless anchor target |
| --- | --- |
| Weak guide `wg` | Late-time `R_B = B_rms/B_mean = 2.0 +/- 0.2`, `M_s = 0.50 +/- 0.05`, and primary `beta_rms approximately 10`. |
| Stronger guide `sg` | Late-time `R_deltaB = delta_B_rms/B_mean = 0.50 +/- 0.05`, `M_s = 0.50 +/- 0.05`, and primary `beta_rms approximately 10`. |

The retained target definitions are:

```math
\boldsymbol{B}_{\rm mean}=\langle\boldsymbol{B}\rangle_V,\qquad
B_{\rm mean}=|\boldsymbol{B}_{\rm mean}|,\qquad
B_{\rm rms}=\langle|\boldsymbol{B}|^2\rangle_V^{1/2},
```

```math
\delta B_{\rm rms}=
\langle|\boldsymbol{B}-\boldsymbol{B}_{\rm mean}|^2\rangle_V^{1/2},
\qquad
p_{\rm iso}={2p_\perp+p_\parallel\over3},
```

```math
M_s={u_{\rm rms}\over
\sqrt{(5/3)\langle p_{\rm iso}\rangle_V/\langle\rho\rangle_V}},
\qquad
\beta_{\rm rms}={2\langle p_{\rm iso}\rangle_V\over B_{\rm rms}^2}.
```

The weak-guide definition deliberately uses total `B_rms`; the stronger-guide
definition deliberately uses fluctuating `delta_B_rms`. Report both magnetic
ratios for both anchors and do not interchange these targets during tuning.

Calibrate anchors with a calibration-only seed and hold science seeds out of
tuning. Measure each collisionless anchor turnover time
`tau_eddy,0,G` and define the background-collision scan:

```math
C_G = \nu_{\rm coll}\tau_{{\rm eddy},0,G}
    \in \{0,\ 0.1,\ 1,\ 10,\ 100\}.
```

For nonzero `C_G`, keep each accepted anchor's imposed field, forcing,
thermodynamics, LF scale, and seed fixed while changing only `nu_coll` and
case metadata. Drift in magnetic amplitude, Mach number, or beta is then a
physical response rather than a calibration failure. Separately tuned
state-matched runs may be secondary comparisons only.

### 10.3 Candidate Stage II matrix

This is the retained scientific design to be recosted after Stage I, not an
execution authorization. A post-reproduction approval may reduce it only by
recording which resulting inference is relinquished.

| Tier | Guide regimes and `C_G` | Resolutions | Models and seeds | Purpose |
| --- | --- | --- | --- | --- |
| Calibration | `wg`, `sg`; `0` | `32^3`, then `64^3` | Active finalists and matched passive confirmation; calibration seed `161803` only. | Freeze anchor inputs and `tau_eddy,0,G`; never enter scientific averages. |
| Collisionless anchor ladder | `wg`, `sg`; `0` | `64^3`, `96^3`, `192^3`, `384^3` | Active/passive; seed `271828`. | Resolve guide-field and feedback baselines. |
| Collisionless seed repeat | `wg`, `sg`; `0` | `96^3`, `192^3` | Active/passive; seed `314159`. | Estimate realization uncertainty. |
| Background-collision sweep | `wg`, `sg`; `0.1`, `1`, `10`, `100` | `96^3`, `192^3` | Active/passive; seed `271828`. | Measure fixed-input `nu_coll` response. |
| Strongly collisional endpoint | `wg`, `sg`; `100` | `384^3` | Active/passive; seed `271828`. | Test resolution of the largest collision contrast. |
| Collisional seed repeat | `wg`, `sg`; `1`, `100` | `192^3` | Active/passive; seed `314159`. | Check sampling at the transition and endpoint. |
| Selected sensitivities | Both guide regimes; collisionless and `C_G = 1` roles selected before execution. | `96^3`, with `192^3` beta checks only if budgeted. | Active/passive beta checks at `beta_rms approximately 1,100`; active-only LF scale factors `1/2,2`; seed `271828`. | Bound beta and closure-scale dependence of principal inferences. |

Stable production identities should encode guide regime, collision coordinate,
feedback choice, resolution, and seed, for example
`g_wg_c1_active_N192_s271828`. For every active/passive pair, reuse the
identical archived forcing realization. For every fixed-input collision
series, change only `nu_coll` and identifying metadata after its collisionless
anchor has been frozen.

### 10.4 Extension safeguards and acceptance contract

| Safeguard or secondary check | Requirement |
| --- | --- |
| Held-out calibration | Use tuning-only seed `161803`; science seeds `271828` and `314159` are not opened before anchor freeze. |
| Low-field validity | Record `min(|B|/B_mean)`, `f_low(q)` for `q = 0.01, 0.05, 0.10`, and field-direction regularization counters; no unqualified claim if local-field validity fails. |
| Quantitative structure | Use isotropic shell spectra as primary weak-guide scale-space measures and local-field structure functions for anisotropy; guide-relative spectra are labeled secondary. |
| Visual selection | Pre-register midpoint snapshots, central planes, and pooled percentile color limits before reviewing visuals. |
| Beta sensitivity | After primary `beta_rms approximately 10`, assess selected cases at `beta_rms approximately 1` and `100`. |
| LF-scale sensitivity | At minimum vary `lf_k_parallel` by factors of two around the baseline in selected anchor/transitional cases. |
| Estimands | Separate active/passive feedback effects, fixed-input background-collisionality effects, and weak/strong-guide contrasts. |
| Late-time sampling | Require a stationary accepted window spanning at least four measured turnover times, with block/bootstrap uncertainty retained for reported statistics. |
| Safety and accounting | Require zero invalid-state/undefined-direction failures and retain forcing, applied pressure/LF work, limiter, and energy-residual evidence. |
| Resolution claims | Describe a guide/collision effect as resolution-robust only after its sign survives the declared seed checks and the applicable `384^3` endpoint comparison. |

The minimum Stage II analysis bundle must contain target/stationarity
decisions, low-field and limiter exposure, active/passive and collision
contrasts with uncertainty, isotropic spectra and local-field structure
products, selected common-scale slices, exact submitted-input/build hashes,
and compute/storage accounting. No Stage II result enters the manuscript
unless that bundle records its accepted case and comparison role.

### 10.5 Stage II recosting rule

Do not reuse the former `3735.627090` node-hour extension reservation after
Stage I is inserted ahead of it. After Stage I:

1. Use measured MKS24 performance and remaining ledger rather than the prior
   speculative scaling.
2. Decide whether a reduced extension can answer the main new question within
   the remaining `2930.259441` node-hours or whatever remainder is actually
   measured.
3. If the full two-guide, five-collisionality, beta/LF-sensitivity, paired
   `384^3` program remains scientifically necessary, request an amended total
   budget before submission.

## 11. Execution Order and Stop Conditions

| Phase | Work | Deliverable | Stop condition |
| --- | --- | --- | --- |
| A. Reproduction specification audit | Verify pinned MKS24 source, sixteen-run mapped alias map, disposition of the unmapped active-Alfvenic beta-1 definition, closure/forcing/limiter identity, and figure/status table. | Reviewed Stage I protocol. | A published simulation role or observable remains unidentified. |
| B. Reference-data closure | Resolve dimensional panel normalization/data boundary and the Figure 10 external-model decision. | Qualified reference manifests or explicit scoped limitation. | A "complete" claim would depend on unqualified data. |
| C. Production readiness | Validate required local/GPU/restart/work-accounting tests and implement reviewed production accounting. | Passing readiness review and fail-closed ledger. | Numerical gates or submission controls fail. |
| D. MKS24 production | Execute the sixteen mapped Stage I runs conditionally, plus the beta-1 case only if Phase A assigns it a published role; run one at a time with segment reviews. | Accepted run bundles and measured ledger/storage. | Safety, stationarity, accounting, storage, or budget gate fails. |
| E. MKS24 analysis | Generate each required product and complete the figure-by-figure status table. | Quantitative reproduction report. | Any claimed result lacks qualified comparison evidence. |
| F. Reproduction manuscript | Transform the TeX note into a buildable Stage I manuscript. | Manuscript and reproducibility package. | Claims exceed the accepted status table. |
| G. Extension decision | Recalculate the Stage II design using measured Stage I results/costs. | Approved extension matrix or deferred-work record. | Insufficient budget or unresolved baseline reproduction. |
| H. Extension execution/manuscript integration | Execute only an approved Stage II matrix and extend the manuscript if scientifically supported. | New-result supplement or expanded manuscript. | Any extension gate fails or funding is absent. |

## 12. Immediate Handoff

As of 2026-05-26:

1. The repository contains the TeX validation note, the writing style guide,
   a detailed MKS24 implementation/evidence plan, substantial MKS24 analysis
   infrastructure, guarded inputs for sixteen source-mapped CGL-LF production
   roles plus one unresolved active-Alfvenic beta-1 definition, a committed
   sixteen-run mapped execution manifest, a production-QOS accounting utility
   with terminal-time inspection and accepted-bundle assembly gates, and
   reduced/GPU qualification evidence.
2. Frontier job `4674731` completed as the legacy shared-output
   `R16/s00_t0_t2` diagnostic segment and was recorded `clean_partial` at
   `t = 1.66204848945`. The ranked-output replacement is inspected and
   accepted through its `t = 6.5` prefix via jobs `4675322`, `4676557`,
   `4676696`, `4679803`, `4681476`, and `4683291`;
   its lineage-only diagnostic prefix bundle passed a separately labeled
   transient-window analyzer pipeline check without promoting that prefix to
   a reproduced result. Segment `R16/s05_rankio_t5_t6p5`
   terminated cleanly on application walltime at `t = 6.43828031751 < 6.5`,
   was recorded `clean_partial` with complete terminal ranked products and
   zero strict LF safety counters, and used `1.753333` node-hours. Short
   continuation `R16/s06_rankio_t6p438280_t6p5` reached exact `t = 6.5`,
   was recorded `accepted` after formal inspection, and used `0.080556`
   node-hours. The accepted-prefix `t = 6.5` bundle contains twenty-eight
   deduplicated ranked snapshots but remains outside the production
   comparison window. Continuation `R16/s07_rankio_t6p5_t7p5` reached
   exact `t = 7.5`, was recorded `accepted`, and used `1.281944`
   node-hours; its generic prefix analysis still selects zero samples in the
   `t = 8`--`10` window. Continuation `R16/s08_rankio_t7p5_t8p5` is
   submitted as job `4686032` toward the first production-window prefix.
   No complete `t = 10` `R16` run bundle or reproduced panel exists yet.
3. It does **not** contain completed paper-scale MKS24 production
   simulations, completed quantitative comparisons for the full figure
   program, resolved dimensional reference mappings for all panels,
   independent reproduction of the Figure 10 hybrid-kinetic comparator, or a
   completed paper manuscript.
4. The next task is not to launch the weak-guide extension. It is to monitor,
   inspect and account `R16/s08_rankio_t7p5_t8p5` after job `4686032`
   completes, continue only inspected Stage I segments sequentially to their
   required durations, assemble accepted per-case analysis bundles, and then
   complete the panel/status and manuscript program.

## 13. Reproducibility Record

For every accepted Stage I or later Stage II calculation, archive:

1. Git revision, executable/build/platform metadata, exact submitted input,
   workflow role, model/forcing/seed metadata, parent restart, and case alias
   mapping.
2. Model controls: `beta0`, domain/resolution, LF scale, thresholds,
   `nu_coll`, limiter scattering policy/rate, forcing shell/injection/time
   correlation, and output cadence.
3. Consumed/reserved node-hours, walltime, storage footprint, retention and
   pruning decisions, and allocation identifier.
4. Analysis intervals, raw histories, required snapshots, products, scripts,
   comparison reference checksums, normalization rules, uncertainties, figure
   files, and pass/fail/blocked decisions.
5. All failures and excluded results; an unsuccessful case may not be silently
   replaced by a selected successful realization.

Before drafting a result statement, classify its evidence:

| Verb | Use only when |
| --- | --- |
| `reproduces` | A matched published result passes the declared quantitative comparison and provenance gates. |
| `shows` | A result is directly measured in accepted simulations and survives its stated uncertainty/gates. |
| `suggests` | An inference is supported but limited by model, sampling, reference-data, or resolution bounds. |
| `does not establish` | The desired claim, including perfection or complete reproduction across an unresolved external-model boundary, exceeds evidence. |

## 14. References and Requirement Traceability

### 14.1 References and local sources

1. Majeski, S., Kunz, M. W., and Squire, J. (2024),
   *Self-organization in collisionless, high-beta turbulence*,
   Journal of Plasma Physics 90, 535900601, arXiv:2405.02418.
2. Squire, J., Schekochihin, A. A., Quataert, E., and Kunz, M. W. (2019),
   *Magneto-immutable turbulence in weakly collisional plasmas*,
   Journal of Plasma Physics, arXiv:1811.12421.
3. Squire, J., Kunz, M. W., Arzamasskiy, L., Johnston, Z., Quataert, E.,
   and Schekochihin, A. A. (2023), *Pressure anisotropy and viscous heating
   in weakly collisional plasma turbulence*, Journal of Plasma Physics,
   arXiv:2303.00468.
4. `docs/cgl_lf_mks24_reproduction_implementation_plan.md`: implementation,
   reference-extraction, qualification, and current evidence record.
5. `docs/writing_style_guide.md`: manuscript structure and claim-calibration
   requirements.
6. `docs/cgl_lf_validation.tex`: current manuscript source to convert.

### 14.2 Traceability

| Requirement | Where specified | Evidence required before completion |
| --- | --- | --- |
| Reproduce MKS24 before extension | Sections 1-2 and 11 | Stage I status table and accepted manuscript claim boundary. |
| Match MKS24 CGL-LF physical/numerical setup | Section 3 | Audited manifests, accepted production cases, and provenance. |
| Include MKS24 collisionality result correctly | Sections 3.4, 4, and 5 | Beta-100 `nu_lim` comparisons for `20`, `200`, and hard wall; no conflation with `nu_coll`. |
| Cover all MKS24 numerical figures/results and theory appendices | Section 5 | Per-panel pass/fail/blocked/external-model table and checked derivation record. |
| Address non-CGL Figure 10 honestly | Sections 2.1 and 5 | External dataset/kinetic rerun or explicitly limited reproduction claim. |
| Follow `writing_style_guide.md` in a complete manuscript | Section 9 | Converted TeX, figure/provenance package, and successful PDF build. |
| Use up to `4000` node-hours responsibly | Section 8 | Fail-closed ledger, measured updates, and Stage I-first reservation. |
| Retain requested weak/strong-guide and uniform-collisionality study | Section 10 | Post-reproduction, recosted and approved Stage II protocol. |
