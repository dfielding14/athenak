# Weak-Guide-Field CGL-LF Turbulence Manuscript and Validation Plan

Status: execution plan, written 2026-05-25. No weak-guide production run is
authorized or claimed by this document.

Primary manuscript source: `docs/cgl_lf_validation.tex`

Governing prose guide: `docs/writing_style_guide.md`

Related implementation record:
`docs/cgl_lf_mks24_reproduction_implementation_plan.md`

Completion boundary: this document defines the work required to produce the
paper; it is not evidence that the weak-guide runs or manuscript results
already exist. Future agents should change a phase from planned to complete
only when the named artifacts and gates in this document have been retained.

## 1. Purpose and claim discipline

The central physical uncertainty is whether the AthenaK CGL-Landau-fluid
(CGL-LF) implementation remains numerically controlled and physically
informative in turbulent states where the local magnetic field is not a small
perturbation of a strong imposed guide field. A weak guide field makes this a
particularly severe test: the field direction wanders, the magnitude of the
field can vary strongly unless the dynamics suppress that variation, and
local pressure-anisotropy and heat-flux operations must remain robust.

The requested conclusion that the method is "perfect" is not a valid
scientific objective. No finite numerical campaign can establish perfection,
and `docs/writing_style_guide.md` specifically requires claims no stronger
than their assumptions and evidence. The paper objective is therefore:

> Determine whether the CGL-LF method passes a demanding weak-guide-field,
> subsonic driven-turbulence validation campaign, identify any failure modes,
> and quantify which physical conclusions are supported after numerical and
> statistical convergence tests.

The manuscript may make the following claims only after the indicated gates
pass.

| Claim level | Permitted claim | Required evidence |
| --- | --- | --- |
| Implementation | The weak-guide case is implemented and reproducible. | Archived inputs, executable revision, forcing seed, manifest, restart test, and diagnostic definitions. |
| Numerical validation | The run is numerically controlled in the tested regime. | Strict safety gates, energy/work accounting, convergence of reported statistics, and LF/limiter diagnostics. |
| Target attainment | The suite samples the requested state. | Stationary late-time windows with defined `R_B` and `M_s` inside their tolerances. |
| Physical inference | Pressure anisotropy changes weak-guide turbulence in measured ways consistent with magneto-immutability. | Active/passive matched comparisons, uncertainty estimates, and robust signatures in strain, heating, residual energy, and structure. |
| Not permitted | The method is perfect, universally valid, or a kinetic-plasma replacement. | These claims are not established by fluid-model simulations. |

## 2. Source motivation and the new regime

Squire et al. (2019) proposed magneto-immutability: in weakly collisional,
high-beta plasma turbulence, pressure-anisotropy stresses organize motions to
resist changes in magnetic-field strength. Their Braginskii-MHD simulations
used periodic driven boxes, including cubic, isotropically forced,
trans-Alfvenic turbulence; they highlighted reduced
`\hat{b}\hat{b}:\nabla u`, modified residual energy, spectra, local-field
structure functions, and the need for higher-resolution structural tests.

Squire et al. (2023) studied mean-field Alfvenic turbulence with a
CGL-Landau-fluid closure. That work motivates using active versus passive
pressure feedback, instability-limited pressure anisotropy, pressure-work and
heat-flux energetics, spectra, PDFs, and structure functions to determine
whether a CGL-LF state retains an MHD-like cascade while dissipating energy
through pressure anisotropy near the driving scale.

The proposed study is not a reproduction of either paper. It combines the
weak-guide, large-fluctuation stress of the 2019 motivation with the CGL-LF
closure and diagnostics motivated by the 2023 work. Its distinguishing
late-time requirements are

```text
weak guide field:  R_B = B_rms / B_mean approximately 2
subsonic flow:     M_s approximately 0.5
```

The study tests whether the same physical organization survives when a global
mean-field decomposition is no longer the dominant geometrical guide.

### 2.1 What is inherited and what is new

The motivating papers constrain the question but do not already answer it.
Squire et al. (2019) report trans-Alfvenic cubic runs with
`\delta B_\perp/B_0 approximately 1`, isotropic forcing, spectra, and
local-field structure functions; their reported highest trans-Alfvenic
resolution is `192^3`, and they explicitly identify higher-resolution
structural measurements as future work. Under the definitions adopted below,
the requested `R_B = 2` corresponds approximately to
`\delta B_{\rm rms}/B_{\rm mean} = \sqrt{3}`. The proposed calculation is
therefore a weaker-guide, larger-field-wandering extension, not a direct
replication of their trans-Alfvenic state.

Squire et al. (2023) supply the closer closure precedent: their CGL-LF
calculations motivate active/passive comparisons, pressure-anisotropy
heating, heat-flux accounting, spectra, PDFs, and scale-dependent structure.
Their mean-field setup does not establish the weak-guide limit. The present
paper must consequently make two comparisons separately: whether AthenaK
recovers closure-consistent numerical behavior, and whether its new
weak-guide statistics support a magneto-immutability interpretation.

## 3. Definitions that must be frozen before running

The phrase `B_rms ~ 2 * B_mean` is ambiguous unless `B_rms` is defined.
Production input manifests and the paper must record both measures:

```math
\boldsymbol{B}_{\rm mean} = \langle \boldsymbol{B} \rangle_V,\qquad
B_{\rm mean} = |\boldsymbol{B}_{\rm mean}|,
```

```math
B_{\rm rms} = \langle |\boldsymbol{B}|^2 \rangle_V^{1/2},\qquad
\delta B_{\rm rms}
= \langle |\boldsymbol{B}-\boldsymbol{B}_{\rm mean}|^2\rangle_V^{1/2}.
```

The primary interpretation of the requested target is

```math
R_B = B_{\rm rms}/B_{\rm mean} = 2.0 \mathbin{+/-} 0.2.
```

For a conserved uniform mean field, this corresponds to
`\delta B_{\rm rms}/B_{\rm mean} \simeq \sqrt{3}` when the volume-average
cross term vanishes. The paper must also show
`R_{\delta B} = \delta B_{\rm rms}/B_{\rm mean}`. If the scientific intent is
instead `R_{\delta B} approximately 2`, that target must be amended before
production execution rather than silently substituted during analysis.

For an anisotropic-pressure calculation, "sonic Mach number" also needs an
explicit scalar convention. Use

```math
p_{\rm iso} = (2p_{\perp}+p_{\parallel})/3,\qquad
c_s = \left(\gamma\langle p_{\rm iso}\rangle_V/
                 \langle\rho\rangle_V\right)^{1/2},\quad \gamma=5/3,
```

```math
M_s = u_{\rm rms}/c_s,\qquad
u_{\rm rms}=\langle|\boldsymbol{u}-\langle\boldsymbol{u}\rangle_V|^2
\rangle_V^{1/2}.
```

The target is the late-time window average

```math
\overline{M_s} = 0.50 \mathbin{+/-} 0.05.
```

Directional CGL acoustic/characteristic-speed diagnostics should be archived
as secondary quantities, but not substituted for the declared manuscript
metric.

Use code-unit beta definitions consistent with the current AthenaK input
convention `p0 = 0.5 * beta0 * B0^2`:

```math
\beta_{\rm mean} = {2\langle p_{\rm iso}\rangle_V\over B_{\rm mean}^2},
\qquad
\beta_{\rm rms} = {2\langle p_{\rm iso}\rangle_V\over B_{\rm rms}^2}.
```

At `R_B = 2`, `beta_mean = 4 beta_rms`. The baseline campaign should target
late-time `beta_rms approximately 10`, requiring an initial calibration near
`beta0 approximately 40` when `beta0` is defined against the imposed mean
field. This keeps the turbulent field in a high-beta CGL-LF regime without
confusing mean-field beta with the dynamically sampled rms-field beta.

## 4. Required observables

The following dimensionless controls and responses must be present in every
production manifest or analysis bundle.

| Quantity | Definition or purpose |
| --- | --- |
| `R_B`, `R_deltaB` | Requested weak-guide target and its fluctuation-resolved interpretation. |
| `M_s` | Requested subsonic target using the frozen scalar convention above. |
| `beta_mean`, `beta_rms` | Separate imposed-field and sampled-field thermodynamic regimes. |
| `M_A,mean`, `M_A,rms` | Alfvenic amplitudes against mean and rms fields; prevent an imprecise use of "trans-Alfvenic." |
| `Delta p = p_perp - p_parallel` and `beta Delta` | Pressure-anisotropy response and instability proximity. |
| `\hat{b}\hat{b}:\nabla u` | Direct magneto-immutability diagnostic motivated by the 2019 paper. |
| `-Delta p \hat{b}\hat{b}:\nabla u` | Anisotropic pressure-work/heating density. |
| LF heat-flux work/contraction and counters | Closure behavior and numerical safety. |
| Mirror/firehose/hard-wall fractions | Extent to which inferred physics depends on modeled microinstability bounds. |
| Forcing work and energy residual | Global accounting gate. |
| Magnetic/kinetic residual energy | Direct comparison to a major 2019 signature. |
| Local-field structure functions and spectra | Scale-dependent structure and convergence. |

The interruption-number language should be used carefully. The 2019
Braginskii interruption number and the 2023 CGL-LF interpretation motivate
regime plots, but the manuscript must state which closure-specific quantity is
computed and may not identify a CGL-LF measurement with a Braginskii
coefficient without a derivation.

## 5. Existing implementation state and required additions

### 5.1 What exists

The current paper-standard deck
`inputs/cgl_lf_paper/cgl_lf_paper_standard_active_alfvenic_beta10.athinput`
provides an active CGL-LF turbulence baseline with a strong imposed
`b0 = 1`, `beta0 = 10`, LF heat flux, forcing-work accounting, strict
limiter handling, and an analysis interval at `t = 8` to `10`.
`src/pgen/tests/cgl_lf_paper.cpp` currently initializes a uniform guide field
and exposes a turbulence paper mode. `scripts/analyze_cgl_lf_paper.py`
already contains substantial pressure, energy, spectrum, and paper-comparison
infrastructure.

The retained Frontier record establishes only a reduced corrected nonlinear
qualification case, not a paper-scale production result. It records
`0.851670` node-hours already used, no active reservation, and a corrected
one-node/eight-GPU `96 x 96 x 192` run reaching `t = 2.0` in `2080` allocated
seconds. All new reservations must be charged against that same ledger.

### 5.2 Code and input work before node-hour spending

Implement and test the following before submitting calibration jobs:

1. Add weak-guide input templates under `inputs/cgl_lf_paper/` for calibration,
   active production, and matched passive-feedback control. Use a periodic
   cubic domain and a verified isotropic solenoidal Ornstein-Uhlenbeck forcing
   configuration. Do not use forcing tied to a fixed perpendicular direction
   unless the paper explicitly treats it as a separate mean-field bridge case.
2. Add history outputs or deterministic snapshot reductions for
   `B_mean`, `B_rms`, `delta_B_rms`, `u_rms`, `M_s`, both beta definitions,
   both Alfvenic Mach definitions, residual energy, and the selected
   magneto-immutability/heating diagnostics.
3. Extend `scripts/analyze_cgl_lf_paper.py` with weak-guide target attainment,
   stationarity, active/passive comparison, uncertainty, slice, and
   convergence products. Each product must archive its formula and analysis
   window.
4. Add workflow case names and manifests without weakening the existing
   Frontier authorization and budget checks. A new production workflow must
   require explicit reservation and must refuse execution when its projected
   ceiling exceeds the remaining ledger.
5. Segment long output runs into a sparse spin-up stage and a dense analysis
   stage restarted from the spin-up checkpoint. This preserves high-quality
   visuals without writing full-resolution snapshots unnecessarily throughout
   the transient.

### 5.3 Numerical qualification gates

No statistical-production case is admissible until all applicable gates pass:

| Gate | Passing evidence |
| --- | --- |
| Unit and reduced regressions | Existing CGL-LF, hard-wall, restart, forcing-work, and analyzer tests remain passing after new diagnostics. |
| Weak-field robustness | A reduced run traverses the calibrated field-amplitude range without unaccounted floors, invalid field-direction operations, NaNs, or strict-bound violations. |
| Active/passive identity | Matching initial state and forcing manifest, differing only in declared pressure-feedback physics. |
| Restart identity | A split weak-guide run reproduces an uninterrupted reduced reference within established tolerances. |
| GPU sizing | A short representative cubic run records memory, output size, runtime, and LF timestep behavior before a `192^3` or larger production submission. |
| Energy accounting | Active production intervals close the forcing/energy/pressure-work accounting within a declared tolerance validated on reduced tests. |

### 5.4 Implementation and artifact map

Do not reimplement analysis products that are already present. The current
analyzer supplies snapshot pressure-anisotropy PDFs, local parallel strain,
pressure-work decomposition, pressure-transfer products, kinetic/magnetic and
pressure spectra, a heat-flux transport proxy, and opt-in local-field eddy
anisotropy. The weak-guide campaign needs to expose those products through a
new case family and add only the target, stationarity, visualization, and
comparison layers required by this experiment.

| File or new artifact | Required change | Acceptance evidence |
| --- | --- | --- |
| `inputs/cgl_lf_paper/cgl_lf_paper_weak_guide_active.athinput` | Add a cubic active base deck with isotropic forcing and placeholders fixed by calibration. | Parsed input metadata matches the contract in Section 6.1. |
| `inputs/cgl_lf_paper/cgl_lf_paper_weak_guide_passive.athinput` | Add the matched passive-feedback base deck. | A test confirms it differs from active only in feedback-declaring controls. |
| `src/pgen/tests/cgl_lf_paper.cpp` and history support | Retain vector mean field, `B_rms`, `delta_B_rms`, volume-weighted `u_rms`, scalar `M_s`, both beta values, and residual energy at history cadence. | Synthetic or reduced-run regression evaluates each formula and preserves it over restart. |
| `scripts/analyze_cgl_lf_paper.py` | Reuse existing turbulence products; add weak-guide target windows, block statistics, seed/resolution comparisons, and an acceptance JSON product. | Reduced synthetic and snapshot tests exercise pass, fail, off-target, and nonstationary cases. |
| `scripts/plot_cgl_lf_paper.py` | Add target histories/convergence and two-dimensional field/derived-field slice rendering with common active/passive limits. | Figure-generation regression lists all required output files. |
| `scripts/cgl_lf_workflow.py` | Add `paper-weak-guide-calibration`, `paper-weak-guide-production`, and `paper-weak-guide-analyze` case families and manifest metadata. | Workflow tests enumerate the fixed case matrix and refuse expensive execution without authorization. |
| Production allocation/accounting path | Provide reviewed production-QOS reservation/submission/accounting outside `scripts/frontier/cgl_lf_frontier.py`. That existing utility is debug-only and deliberately rejects paper production. | Offline budget self-test plus retained production reservation records; the debug-only rejection remains passing. |
| `docs/cgl_lf_validation.tex` | Convert the note into the manuscript structure and figures specified below. | Warning-free PDF build plus figure/data provenance table. |

## 6. Simulation campaign

### 6.1 Physical setup

The primary weak-guide suite should use:

| Choice | Baseline |
| --- | --- |
| Domain | Periodic cubic box, `L_x = L_y = L_z = L = 1`. |
| Mean field | Uniform `B_mean` along `z`; its amplitude is calibrated, not presumed. |
| Mesh | `N x N x N`; choose the meshblock decomposition after the GPU sizing run without changing `N` or physical parameters. |
| Forcing mode | `<turb_driving>/driving_type = 0`, the implemented isotropic solenoidal Ornstein-Uhlenbeck velocity-forcing mode. |
| Forcing shell | `physical_k_shell = true`, `k_shell_unit = 2*pi/L = 6.283185307179586`, `nlow = 1`, `nhigh = 3`, `isotropic_power_spectrum = true`, `expo = 2.0`. |
| Seeds | Primary ladder seed `271828`; second-seed cases `314159`; calibration may use the primary seed only. |
| Forcing amplitude/time | Calibrate `dedt`; start the bracket from the currently exercised `dedt = 0.32`. Freeze `tcorr` after requiring it to be within a factor of two of the measured outer-scale turnover time; start from `tcorr = 2.0`. |
| Primary target | Late-time `R_B = 2.0 +/- 0.2`, `M_s = 0.50 +/- 0.05`, and `beta_rms approximately 10`. |
| Thermodynamics | `rho0 = 1`; begin calibration with `beta0 = 40` (the initial `beta_mean`), `p_parallel0 = p_perp0 = 0.5 * beta0 * b0^2`, then freeze the accepted calibrated value. |
| Closure | Active CGL-LF with currently qualified hard-wall instability treatment and strict diagnostics. |
| Matched control | Passive-pressure-feedback case with identical forcing realization and thermodynamic initialization. |
| LF scale | Baseline `lf_k_parallel = 2*pi/L = 6.283185307179586`, interpreted as a closure assumption in a wandering-field box; active sensitivity cases use `pi/L` and `4*pi/L`. |
| Limiter/numerics | `cgl_firehose_threshold = parallel`, `limiter_hardwall = true`, strict LF admissibility, applied pressure-work histories, forcing-work history, and explicit metadata for integrator/CFL/STS controls. |
| Time interval | Provisional `t_end = 12`; use `t = 0` to `8` for spin-up and `t = 8` to `12` for analysis only if the stationarity and turnover-span criteria below pass. |

The cubic, isotropically forced baseline follows the most relevant structural
test in Squire et al. (2019), where local-field anisotropy could be measured
without anisotropy being imposed by elongated outer-scale forcing. It differs
deliberately from the existing elongated, strong-guide MKS24-oriented deck.
An optional strong-guide bridge case may be added later, but it must not
replace the requested weak-guide campaign.

### 6.2 Calibration procedure

The targets cannot be selected reliably from initial conditions alone because
the saturated fluctuation amplitude depends on forcing, pressure feedback,
limiter activity, and heat flux. Calibrate before increasing resolution:

1. At `32^3`, run the primary-seed active bracket shown below. These are
   tuning runs and must never enter scientific averages.
2. Interpolate or extend the bracket only after plotting late-time
   `(R_B, M_s, beta_rms)`. Select the two nearest candidates to the target
   point; do not select by visual field appearance.
3. At `64^3`, rerun those two candidates through `t = 8`, including a matched
   passive case for the best active candidate. Adjust `p0`/`beta0` only to
   place the chosen `beta_rms` regime after the magnetic and velocity targets
   are plausible.
4. If the best candidate is not inside both target tolerances at `64^3`,
   carry out a second, locally bracketed calibration round and charge it to
   the calibration reservation. Do not advance the resolution ladder.
5. Select a single fixed parameter set before the resolution ladder. Freeze
   the forcing algorithm, shell, correlation time, target definitions, output
   cadence, closure settings, and seeds in a versioned manifest.
6. At each production resolution, require the late-time target tolerance and
   stationarity checks. If either target drifts with resolution, report the
   drift and recalibrate or alter the scientific question; do not relabel an
   off-target run as a successful weak-guide test.

The first active bracket is concrete but revisable based on its measured
response:

| Calibration variable | Values | Reason |
| --- | --- | --- |
| Resolution/duration/seed | `32^3`, `t_end = 4`, `271828` | Cheap response map, excluded from science ensemble. |
| `b0` | `0.25`, `0.50`, `1.00` | Bracket field wandering relative to the imposed field. |
| `dedt` | `0.08`, `0.32`, `1.28` | Bracket driven velocity and magnetic fluctuation amplitudes around the exercised forcing power. |
| `beta0` | `40` for all first-round cases | Starting hypothesis corresponding to `beta_rms approximately 10` at `R_B approximately 2`. |
| Number/cost | Nine active cases; estimated `9 x 0.042798 = 0.385185` node-hours under the conservative scaling | Leaves room for confirmation, sizing, and local refinements inside the `20` node-hour calibration reserve. |

The accepted parameter manifest must record the selected `b0`, `dedt`,
`tcorr`, `beta0`, derived initial pressure, `lf_k_parallel`, seed list, and
the calibration data used to select them.

### 6.3 Late-time acceptance contract

For each reportable case, define the outer-scale turnover time from the
accepted analysis window,

```math
\tau_{\rm eddy} = L_{\rm force}/\overline{u_{\rm rms}},
\qquad L_{\rm force} = 2\pi/k_{\rm peak},
```

where `k_peak` is the energy-injection-spectrum peak archived from the forcing
manifest. The provisional `t = 8` to `12` window is admissible only if it
spans at least four measured `tau_eddy`; otherwise extend the window before
submitting a higher resolution.

Split an admissible window into four equal-duration blocks and emit
`analysis/weak_guide_acceptance.json` with the following decisions:

| Test | Passing definition |
| --- | --- |
| Magnetic target | Full-window mean `R_B` lies in `[1.8, 2.2]`; report `R_deltaB` without replacing this metric. |
| Sonic target | Full-window mean `M_s` lies in `[0.45, 0.55]`. |
| Target stationarity | The difference between first and last block means is below `0.10` for `R_B` and below `0.025` for `M_s`. |
| Energetic stationarity | Conserved total-energy change versus recorded forcing work is consistent with the reduced-test accounting tolerance or a tighter production tolerance frozen before the ladder; applied closure-work histories are reported separately to interpret internal transfers. |
| Safety | Zero nonfinite, nonpositive, floor, emergency-bound, or invalid-direction failures; report hard-wall and cap activity as physics/model exposure rather than silently excluding it. |
| Sample sufficiency | Four turnover times and at least four block summaries; bootstrap intervals are emitted for paper statistics. |

If `t = 8` to `12` is too short to establish those conditions at `96^3`,
stop the high-resolution ladder and spend only the reserved contingency on a
reviewed time-extension plan. An off-target or nonstationary case is retained
as a failed calibration/result, not omitted from the record.

### 6.4 Resolution ladder and matched cases

Use one fixed forcing seed at every baseline resolution and a second seed at
intermediate resolutions to quantify realization variance economically.

| Tier | Resolution | Active CGL-LF | Passive control | Seeds | Purpose |
| --- | ---: | ---: | ---: | ---: | --- |
| Calibration | `32^3` | tuning grid | optional finalist only | common | Find `b0` and `dedt`; no manuscript statistics. |
| Confirmation | `64^3` | yes | yes | 1 | Confirm target and workflow after tuning. |
| Low science | `96^3` | yes | yes | 2 | First reportable statistics and LF sensitivity anchor. |
| Fiducial | `192^3` | yes | yes | 2 | Statistical and resolution comparison. |
| High | `384^3` | yes | yes | 1 | Highest approved resolution for convergence and visuals. |
| Sensitivity | `96^3` | two extra active LF-scale variants | no | baseline seed | Dependence on LF closure scale. |

The high-resolution case is paired with a passive control because a beautiful
single active run cannot isolate pressure-feedback physics. The second seed is
spent at `96^3` and `192^3`, where it can test uncertainty without consuming
the high-resolution reserve.

The production manifest should use stable case identities:

| Case pattern | Included values |
| --- | --- |
| `wg_active_N{N}_s271828` and `wg_passive_N{N}_s271828` | `N = 64, 96, 192, 384`. |
| `wg_active_N{N}_s314159` and `wg_passive_N{N}_s314159` | `N = 96, 192`. |
| `wg_active_N96_s271828_kpar_half` | Active-only `lf_k_parallel = pi/L`. |
| `wg_active_N96_s271828_kpar_double` | Active-only `lf_k_parallel = 4*pi/L`. |
| `wg_passive_state_N{N}_s271828` | Conditional state-matched passive controls at `N = 96, 192` only if triggered below. |

Each active/passive pair must reuse an identical archived OU seed and forcing
configuration. It is an implementation failure if the passive run silently
uses a separately generated random realization.

The target gate is mandatory for the primary active ladder. The passive cases
in the fixed matrix are forcing-matched controls: they answer how removing
active pressure feedback changes a calculation with the same input forcing
and mean field. They must report their attained `(R_B, M_s)` values, but they
must not be called state-matched controls if those values fall outside the
active target tolerances.

If a forcing-matched passive control falls outside either tolerance, use the
contingency allocation for passive state-matching runs at `96^3` and `192^3`
with the primary seed, adjusting only the declared calibration controls and
archiving the adjustment. These two additional `t = 12` cases cost
`31.200001` node-hours on the conservative basis. Use forcing-matched
comparisons to discuss feedback at fixed input and state-matched comparisons
to discuss conditional turbulence statistics. Do not add a state-matched
`384^3` case without recomputing and approving the budget.

## 7. Node-hour budget and maximum admissible resolution

The hard ceiling is `1000` Frontier node-hours, including the existing
`0.851670` node-hours in the retained CGL-LF ledger. This plan uses the
corrected G013 nonlinear timing only as a conservative planning basis:
`96 x 96 x 192` through `t = 10` is estimated at `2.888889` node-hours on
one node/eight GPUs. Because the new primary box is cubic and has fewer cells
than that elongated layout at the same nominal resolution, retaining the
elongated scaling below is conservative until a measured cubic sizing run
replaces it.

For provisional `t_end = 12`, cubic-case reservations are conservatively
estimated as

```math
C_N = 2.888889\ {\rm node\ hours}
      \left({N\over96}\right)^3 {12\over10}.
```

| Resolution | Conservative node-hours per `t = 12` case |
| ---: | ---: |
| `32^3` | `0.128395` |
| `64^3` | `1.027161` |
| `96^3` | `3.466667` |
| `192^3` | `27.733334` |
| `384^3` | `221.866675` |

The proposed allocation is:

| Allocation | Cases represented | Node-hours |
| --- | --- | ---: |
| Previously consumed ledger | Retained CGL-LF qualification evidence | `0.851670` |
| Implementation/sizing/calibration reservation | Short GPU sizing plus `32^3`/`64^3` tuning and reruns | `20.000000` |
| Baseline paired ladder | Active and passive at `64^3`, `96^3`, `192^3`, `384^3`, one seed | `508.187674` |
| Second-seed statistical pairs | Active and passive at `96^3` and `192^3` | `62.400002` |
| LF-scale sensitivity | Two additional active `96^3` cases | `6.933334` |
| Failure/time/statistics contingency | Conditional state-matched passive controls, time extensions, or reruns only after a documented gate review | `200.000000` |
| Total committed ceiling | Includes consumed ledger and contingency | `798.372680` |
| Uncommitted margin below `1000` | Preserved for accounting error or a separately reviewed amendment | `201.627320` |

This allocation assumes that an accepted statistical case ends at `t = 12`.
If the turnover-span or stationarity gate requires longer runs, agents must
recompute every affected case reservation before submitting the extension.
The `200` node-hour contingency is not permission to overrun the ledger
silently; it is the only preallocated source for an approved extension or
rerun.

The largest presently admissible paired production resolution is `384^3`.
Under the same conservative evidence, one `512^3` `t = 12` case costs
`525.906193` node-hours; the required active/passive pair alone costs
`1051.812386` node-hours, before calibration or contingency. A measured cubic
benchmark may later justify an amendment, but agents must not schedule a
`512^3` pair on an assumed cell-count speedup.

The existing `scripts/frontier/cgl_lf_frontier.py` utility is not the
submission path for this suite: it is an intentionally debug-only
qualification tool and rejects paper-production decks. Before any production
submission, either implement a separately reviewed production accounting
driver with the same fail-closed `1000` node-hour ledger or retain equivalent
manual reservation records approved for the production queue. Do not remove
the debug rejection to make the new cases executable.

Before every submitted job:

1. Regenerate the ledger summary from the workflow tooling.
2. Reserve the job's conservative upper bound, including queue/runtime limit.
3. Refuse submission if consumed plus active reservations plus the new ceiling
   exceeds `1000.000000` node-hours.
4. Replace estimates with measured costs after each accepted sizing or
   production case, without retrospectively weakening the ceiling.

## 8. Output and storage strategy

Full-cadence high-resolution output through spin-up is unnecessary and may be
more restrictive than compute time. Run each production case as two
restart-connected segments:

| Segment | Time range | Retained output |
| --- | --- | --- |
| Spin-up | `0 <= t <= 8` | History `dt = 0.02`; binary snapshots `dt = 2.0`; restart at least at `t = 4` and `8`, retaining the `t = 8` parent of the analysis segment. |
| Analysis | `8 <= t <= 12` | History `dt = 0.02`; binary snapshots `dt = 0.25`; restart `dt = 2.0`, retaining the final two checkpoints. |

Implement the cadence change by two input manifests joined by a recorded
restart, rather than by manually deleting transient output after a
single high-cadence run. If an accepted analysis interval is extended, retain
the same `dt = 0.25` snapshot cadence unless a storage amendment records a
different choice before execution.

The G011 storage measurements in the existing plan imply that a conservative
`384` elongated-layout binary snapshot is about `4.08 GB` and a restart is
about `14.7 GB`. Retaining approximately 22 snapshots and two final restarts
would therefore reserve about `120 GB` per high-resolution case before the
expected cubic reduction. Confirm actual cubic file sizes in the sizing job,
record the storage reservation alongside node-hours, and prune only according
to an archived retention rule.

Never prune the manifest, submitted input, history data, final restart,
selected figure snapshots, analysis JSON, log, executable revision, or
checksums needed to reproduce a plotted result.

## 9. Analysis and figure program

Every figure should serve an explicit argument, with captions written as
mini-arguments in the manner required by `docs/writing_style_guide.md`.
Figures must distinguish measurements from interpretation and carry the
resolution, model, seed treatment, and time window needed to interpret them.

### 9.1 Target and numerical-control figures

| Product | Argument it must support |
| --- | --- |
| Time histories of `R_B`, `R_deltaB`, `M_s`, `beta_rms`, and injected/total energy | The calculation actually reaches a stationary weak-guide, subsonic state in the declared window. |
| Target-plane plot of late-time `(R_B, M_s)` by resolution/model/seed | Target attainment is not an accident of one resolution or realization. |
| Energy-work residual, hard-wall fraction, LF counters, and extrema | The state is numerically admissible rather than sustained by hidden failure handling. |

### 9.2 Field-structure visuals

At a common statistically selected time, render matched active/passive slices
with identical limits and a documented plane-selection rule:

| Slice field | Reason to show it |
| --- | --- |
| `rho/<rho>` and `|B|/B_mean` | Establish compressibility level and magnetic-strength structure. |
| `Delta p` and `beta Delta` | Show anisotropy and proximity to firehose/mirror bounds. |
| `\hat{b}\hat{b}:\nabla u` | Display the strain that magneto-immutability is predicted to reduce. |
| `-Delta p \hat{b}\hat{b}:\nabla u` | Locate anisotropic pressure work/heating. |
| Limiter mask and LF activity/cap measure | Reveal where modeled regulation affects the apparent physical structures. |

Add magnetic-field-line or streamline overlays only if they are generated from
the same archived three-dimensional snapshot and remain readable. A visual is
supporting evidence, not a substitute for the quantitative distributions and
convergence plots below.

### 9.3 Quantitative turbulence figures

| Figure family | Required comparisons and interpretation |
| --- | --- |
| PDFs and joint PDFs | Active versus passive distributions of `beta Delta`, `|B|/B_mean`, `\hat{b}\hat{b}:\nabla u`, pressure work, and the joint relation between strain and pressure work. |
| Energy and residual-energy statistics | Determine whether active pressure feedback produces the magnetic dominance emphasized in the 2019 motivation. |
| Spectra | Kinetic, magnetic, density, `Delta p`, parallel/perpendicular pressure, strain, and pressure-work spectra; compare active/passive and resolved-scale convergence. |
| Scale-dependent heating/transfer | Determine where injected energy is removed by anisotropic pressure work and LF transport, motivated by the 2023 results. |
| Local-field structure functions | Measure scale-dependent anisotropy without relying on the weak global guide field; follow local-field conditioning in the 2019 analysis. |
| Resolution and seed convergence | Plot statistic differences from the `384^3` case and realization uncertainty; identify which conclusions are stable. |
| LF-scale sensitivity | Show which principal conclusions persist when `lf_k_parallel` changes by factors of two. |

For external-paper comparisons, dimensionless definitions with demonstrably
matching conventions may be discussed directly. Absolute spectral ordinates
must remain excluded from quantitative cross-paper overlays unless the
normalization/provenance gap recorded in the MKS24 implementation plan is
resolved.

### 9.4 Statistical protocol

1. Establish a single late-time window from declared stationarity tests before
   examining preferred physical comparisons.
2. Use turnover-time block averages and block bootstrap intervals for scalar
   measurements and binned diagnostics.
3. Use the two independent seeds at `96^3` and `192^3` to distinguish
   realization variance from resolution trends.
4. Compare spectra only over overlapping resolved ranges and record the
   de-aliasing/binning/FFT normalization metadata.
5. Treat failure of convergence, target attainment, energy closure, or LF
   sensitivity as a scientific result that limits the claim.

### 9.5 Quantitative decision and artifact contract

Before inspecting the `384^3` case, freeze the principal tests in the
production manifest:

| Principal statistic | Comparison |
| --- | --- |
| Parallel-strain suppression | Active/passive ratio of RMS `\hat{b}\hat{b}:\nabla u`. |
| Magnetic-strength regulation | Active/passive comparison of the distribution width of `|B|/B_mean`, while still requiring target attainment separately. |
| Pressure-feedback energetics | Active anisotropic pressure work and LF work normalized by recorded injected work. |
| Residual energy | Active/passive magnetic-minus-kinetic residual-energy comparison. |
| Limiter/closure exposure | Mirror, firehose, hard-wall, and LF-cap fractions and their LF-scale sensitivity. |

An active/passive effect claimed from the first four statistics may be
described as resolution-robust only if its sign is the same for both seeds at
`96^3` and `192^3`, is retained by the `384^3` primary-seed comparison, and
the normalized `192^3` to `384^3` change is no greater than the larger of
`10%` or the measured two-seed uncertainty at `192^3`. Limiter/closure
exposure instead qualifies those claims: a strong or unresolved sensitivity
prevents an unqualified physical conclusion. If a gate fails, the manuscript
may show the result but must label its resolution or model dependence. Apply
the equivalent criterion to spectral or structure-function statements only on
the common resolved range.

Every accepted production bundle must retain these machine-readable and
rendered products:

| Artifact | Contents |
| --- | --- |
| `analysis/weak_guide_acceptance.json` | Target values, block trends, turnover-span decision, safety counters, energy gate, and pass/fail reasons. |
| `analysis/weak_guide_statistics.json` | Block/bootstrap estimates, active/passive effects, seed scatter, resolution comparisons, and LF-scale sensitivity. |
| `analysis/weak_guide_provenance.json` | Executable/input hashes, parent restart, case identity, seed, model choices, output cadence, analysis configuration, and ledger entry. |
| `figures/weak_guide_targets.pdf` | Target/stationarity histories and `(R_B, M_s)` target plane. |
| `figures/weak_guide_slices.pdf` | Common-scale active/passive slices of fields, anisotropy, strain, pressure work, and limiter/LF exposure. |
| `figures/weak_guide_pdfs.pdf` | PDFs and joint PDFs used for magneto-immutability inferences. |
| `figures/weak_guide_energetics.pdf` | Injection, applied pressure/LF work, residual energy, and accounting quality. |
| `figures/weak_guide_spectra.pdf` | Spectra and resolved-range convergence markers. |
| `figures/weak_guide_structure.pdf` | Local-field structure functions/anisotropy and resolution comparison. |
| `figures/weak_guide_sensitivity.pdf` | Seed, resolution, LF-scale, and safety/limiter exposure summary. |

## 10. Manuscript conversion plan

The current `docs/cgl_lf_validation.tex` is an illustrated method-validation
note. Convert it into a complete paper only after the weak-guide diagnostics
and at least the calibration/qualification evidence exist. Preserve useful
method material, but restructure the argument around the physical uncertainty
and the new test.

| Manuscript section | Required content and evidence |
| --- | --- |
| Abstract | Begin with the uncertainty in CGL-LF turbulence under weak guide fields; state the tested regime and only the results supported by completed runs; end with the principal limitation. |
| Introduction | Explain magneto-immutability and why weak guide fields challenge closure/numerics; position the study relative to the 2019 Braginskii and 2023 CGL-LF works. |
| Model and regime parameters | CGL-LF equations, limiter and heat-flux model, `R_B`, `M_s`, beta definitions, Alfvenic measures, and any interruption-related measure with closure-specific qualification. |
| Numerical method and validation | Condense existing operator tests, hard-wall/energy-work corrections, restart/GPU qualification, and new weak-field tests. |
| Turbulence experiment | Forcing, cubic domain, calibration procedure, active/passive pairing, seeds, resolution ladder, analysis window, and budget/provenance. |
| Results I: attained state | Demonstrate stationary `R_B`, `M_s`, beta and safety/energy gates before presenting physical claims. |
| Results II: structures and magneto-immutability | Slices, PDFs, strain reduction, residual energy, limiter occupancy, and active/passive differences. |
| Results III: energetics and cascade | Pressure/LF work, spectra, transfer/heating scale dependence, structure functions, convergence and sensitivity. |
| Discussion | State what is robust, how it agrees or differs from the source motivations, and what remains model-, limiter-, or closure-scale-dependent. |
| Conclusions | Make the strongest defensible claim; do not replace validation with a perfection claim. |
| Reproducibility appendix/data statement | Case table, hashes, ledger, definitions, analysis commands, retained-data location, and figure-to-product mapping. |

Planned main-text figure sequence:

| Figure | Purpose |
| ---: | --- |
| 1 | Regime diagram and campaign design: `R_B`, `M_s`, beta, active/passive cases, and resolution ladder. |
| 2 | Target attainment and numerical-control histories. |
| 3 | Matched active/passive slices of `|B|`, anisotropy, strain, and pressure work. |
| 4 | PDFs/joint PDFs showing suppression or absence of suppression of field-strength-changing strain. |
| 5 | Energy histories and residual energy with convergence and seed uncertainty. |
| 6 | Kinetic, magnetic, pressure, and strain spectra with resolution comparison. |
| 7 | Heating/transfer spectra or cumulative scale-dependent energy pathways. |
| 8 | Local-field structure functions and anisotropy. |
| 9 | Closure-scale and limiter/safety sensitivity summary. |

### 10.1 Mapping the existing TeX note into the paper

Retain `docs/cgl_lf_validation.tex` as the editable manuscript source unless
a reviewed journal-template conversion requires a new filename. The existing
content should be transformed as follows:

| Current note content | Manuscript treatment |
| --- | --- |
| Purpose and implementation narrative | Replace with a styled abstract and introduction beginning from the weak-guide uncertainty. |
| Conserved variables, LF closure, heat-flux cap, and threshold policy | Retain in a compact methods section; move detailed operator exposition to an appendix if it interrupts the physical argument. |
| Method validation problems and five method figures | Retain as numerical-method evidence, condensed into one validation section plus appendix material as needed. |
| Reduced MKS24-oriented product-path figures | Do not present as physical results. Replace in the main narrative with accepted weak-guide figures; retain only where explicitly labeled as diagnostic validation. |
| MKS24 reproduction boundary | Recast as a data/provenance and limitation statement; it is not the central result of the new paper. |
| Current bibliography | Add Squire et al. (2019) and the weak-guide interpretation, retain Squire et al. (2023), and add any directly used method/diagnostic references. |

### 10.2 Complete-manuscript acceptance checklist

The manuscript is complete only when all items below are true:

1. `docs/cgl_lf_validation.tex` has the section order in Section 10, follows
   `docs/writing_style_guide.md`, and nowhere asserts perfection or kinetic
   generality.
2. All claims in the abstract, results, discussion, and conclusion are
   traceable to accepted cases and to the claim verbs defined in Section 12.
3. The case table includes every successful and failed calibration or
   production case used to make a decision, its role, resolution, seed,
   analysis window, and target/gate outcome.
4. The main figures use the accepted analysis artifacts from Section 9.5;
   each caption states the comparison, parameters, principal pattern, and
   conclusion or limitation.
5. The reproducibility appendix records executable/input hashes, analysis
   commands, data-product paths/checksums, output retention, and final
   node-hour accounting.
6. A manuscript PDF builds without warning-producing missing references or
   figures from a clean checkout plus the declared retained-data products.

After the workflow and analysis additions in this plan exist, the expected
paper-product sequence is:

```bash
python3 scripts/cgl_lf_workflow.py paper-weak-guide-analyze \
  --output-dir /path/to/accepted-weak-guide-bundle
cd docs
pdflatex -halt-on-error -interaction=nonstopmode \
  -output-directory /tmp/cgl_lf_manuscript_build cgl_lf_validation.tex
pdflatex -halt-on-error -interaction=nonstopmode \
  -output-directory /tmp/cgl_lf_manuscript_build cgl_lf_validation.tex
```

The workflow command is a required future interface, not a currently
available command; Phase B must implement and test it before this sequence is
used.

## 11. Execution phases and stop conditions

| Phase | Work | Deliverable | Stop condition |
| --- | --- | --- | --- |
| A. Specification | Freeze target definitions, model comparison, forcing choice, manuscript claims, and budget ledger. | Reviewed protocol and input schema. | Definitions remain ambiguous or budget tooling cannot enforce ceiling. |
| B. Implementation | Add weak-guide histories, analyzer products, input/workflow support, and tests. | Passing local/reduced validation bundle. | Existing numerical gates regress or new observables are not reproducible. |
| C. Calibration | Tune `b0`, `dedt`, and thermodynamic regime at `32^3`/`64^3`. | Frozen production manifest meeting target tolerances provisionally. | Targets cannot be reached without unsafe behavior or unstable transients. |
| D. GPU qualification | Run representative short case and restart/output sizing. | Updated measured runtime/storage reservation. | Runtime, memory, LF behavior, energy closure, or strict counters are unacceptable. |
| E. Resolution ladder | Run paired `64^3`, `96^3`, `192^3`, then `384^3` only as gates remain satisfied. | Archived production cases and ledger. | Off-target state, nonstationarity, failed accounting, or budget margin breach. |
| F. Robustness | Run second seeds and LF-scale variants. | Uncertainty and closure-sensitivity products. | Principal inference changes without a qualified explanation. |
| G. Manuscript | Revise TeX, generate figures/tables, write conclusions and caveats. | Buildable complete manuscript and reproducibility package. | Results cannot support stated claims; rewrite conclusions accordingly. |

### 11.1 Current handoff status and first execution sequence

As of 2026-05-25, this repository contains the restored style guide, the
illustrated TeX validation note, existing strong-guide CGL-LF paper-analysis
infrastructure, and reduced GPU qualification/budget evidence. It does not
contain weak-guide base decks, weak-guide target/stationarity products, a
weak-guide slice renderer, a production submission/accounting path, accepted
weak-guide production runs, or a completed weak-guide manuscript.

A future implementation agent should therefore proceed in this order:

1. Implement the target histories, weak-guide base inputs, case-family
   metadata, synthetic target/stationarity tests, and active/passive identity
   test without submitting an HPC job.
2. Implement slice and quantitative analysis products and exercise them on
   reduced local/GPU qualification data; keep all figures explicitly labeled
   non-statistical until production cases are accepted.
3. Establish reviewed production allocation/accounting that preserves the
   debug-only protection of `scripts/frontier/cgl_lf_frontier.py`.
4. Execute calibration and short GPU sizing only after the local/reduced
   gates pass; freeze accepted physical parameters and updated measured cost.
5. Execute the resolution ladder conditionally, analyze accepted cases, and
   then convert the TeX note using the manuscript checklist above.

The high-resolution submission is conditional: do not execute `384^3` merely
because it appears in this plan. It is authorized only after lower-resolution
cases meet targets, demonstrate stationarity, pass numerical gates, and update
the ledger with measured costs.

## 12. Reproducibility record future agents must maintain

For each accepted calculation, archive:

1. Git revision, executable build metadata, platform/compiler/GPU information,
   exact submitted input, case role, forcing seed, and parent restart if any.
2. Consumed and reserved node-hours before and after execution, walltime,
   allocation identifier, storage footprint, and pruning record.
3. Target definitions and late-time intervals, with the raw histories from
   which target decisions were made.
4. Analysis configuration, JSON products, plotted source data, plot scripts,
   figure files, and checksums.
5. Gate results, including failures. Failed cases must remain identifiable and
   may not be quietly replaced by a selected successful realization.

Before drafting a result statement in the manuscript, classify it:

| Verb | Use only when |
| --- | --- |
| `shows` | A defined result is directly visible in accepted simulations and survives its stated uncertainty/gates. |
| `suggests` | A physically motivated inference is supported but limited by model, sampling, or resolution. |
| `may` | The result motivates a test or implication not established here. |
| `does not establish` | A tempting claim, including model perfection or kinetic generality, exceeds the evidence. |

## 13. References and local sources

1. Squire, J., Schekochihin, A. A., Quataert, E., and Kunz, M. W. (2019),
   *Magneto-immutable turbulence in weakly collisional plasmas*,
   Journal of Plasma Physics, arXiv:1811.12421.
   <https://arxiv.org/abs/1811.12421>
2. Squire, J., Kunz, M. W., Arzamasskiy, L., Johnston, Z., Quataert, E.,
   and Schekochihin, A. A. (2023), *Pressure anisotropy and viscous heating
   in weakly collisional plasma turbulence*, Journal of Plasma Physics,
   arXiv:2303.00468. <https://arxiv.org/abs/2303.00468>
3. `docs/writing_style_guide.md`: governing manuscript organization, claim
   calibration, figure-caption, and revision guidance.
4. `docs/cgl_lf_validation.tex`: current CGL-LF validation manuscript source.
5. `docs/cgl_lf_mks24_reproduction_implementation_plan.md`: current
   implementation, evidence, runtime, storage, and budget ledger context.

## 14. Requirement traceability

| Requested element | Where this plan specifies it | Completion evidence for the eventual paper project |
| --- | --- | --- |
| Follow `writing_style_guide.md` | Sections 1, 9, 10, and 12 | Manuscript claim/caption review and buildable TeX/PDF. |
| Turn `cgl_lf_validation.tex` into a complete manuscript | Sections 10 and 11 | Converted TeX source, generated PDF, data/figure provenance, and checked claims. |
| Weak guide field with `B_rms approximately 2 B_mean` | Sections 3, 6.1, 6.2, and 6.3 | Acceptance JSON with declared `R_B` definition and passed late-time target. |
| Sonic Mach number `approximately 0.5` | Sections 3, 6.1, and 6.3 | Acceptance JSON with scalar `M_s` definition and passed late-time target. |
| Resolution range as large as possible within `1000` node-hours | Sections 6.4 and 7 | Fail-closed ledger, measured updates, conditional paired ladder through the largest admissible resolution. |
| Robust turbulence analysis and strong visuals | Sections 4 and 9 | Retained JSON products, common-scale slices, quantitative figures, convergence and uncertainty tests. |
| Motivation from arXiv:1811.12421 and arXiv:2303.00468 | Sections 2 and 13 | Literature-positioned introduction/discussion and diagnostics tied to those motivations without claiming reproduction. |
| Lead future agents | Sections 5, 6, 7, 8, 9.5, 10.2, 11, and 12 | File-level work map, cases, gates, budget, retained products, and immediate work order. |
