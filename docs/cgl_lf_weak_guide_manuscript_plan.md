# Guide-Field and Collisionality Dependence of CGL-LF Turbulence: Manuscript and Validation Plan

Status: execution plan, written 2026-05-25. No simulation in the proposed
campaign has been run or is claimed by this document. Revised to include
collisionless-to-collisional and weak-to-strong-guide-field comparisons under
a `4000` node-hour ceiling.

Primary manuscript source: `docs/cgl_lf_validation.tex`

Governing prose guide: `docs/writing_style_guide.md`

Related implementation record:
`docs/cgl_lf_mks24_reproduction_implementation_plan.md`

Completion boundary: this document defines the work required to produce the
paper; it is not evidence that the proposed runs or manuscript results already
exist. Future agents should change a phase from planned to complete only when
the named artifacts and gates in this document have been retained.

## 1. Purpose and claim discipline

The central physical uncertainty is how collisionality changes turbulent
self-organization in the AthenaK CGL-Landau-fluid (CGL-LF) model, and whether
that change depends on the strength of the imposed guide field. A weak guide
field makes this a particularly severe numerical and physical test: the field
direction wanders, the magnitude of the field can vary strongly unless the
dynamics suppress that variation, and local pressure-anisotropy and heat-flux
operations must remain robust. A matched stronger-guide complement determines
which conclusions are specific to that severe weak-guide geometry.

The requested conclusion that the method is "perfect" is not a valid
scientific objective. No finite numerical campaign can establish perfection,
and `docs/writing_style_guide.md` specifically requires claims no stronger
than their assumptions and evidence. The paper objective is therefore:

> Determine how pressure-anisotropy feedback and heating in subsonic driven
> turbulence change from collisionless to increasingly collisional CGL-LF
> regimes, compare weak- and stronger-guide states, identify any failure
> modes, and quantify which conclusions survive numerical and statistical
> convergence tests.

The manuscript may make the following claims only after the indicated gates
pass.

| Claim level | Permitted claim | Required evidence |
| --- | --- | --- |
| Implementation | Both guide-regime anchors and their collisionality cases are implemented and reproducible. | Archived inputs, executable revision, forcing seeds, collisional manifests, restart tests, and diagnostic definitions. |
| Numerical validation | The run is numerically controlled in the tested regime. | Strict safety gates, energy/work accounting, convergence of reported statistics, and LF/limiter diagnostics. |
| Target attainment | The anchor suite samples the requested weak- and strong-guide states. | Stationary late-time windows with defined magnetic-amplitude targets and `M_s` inside their tolerances before the collisionality scan is frozen. |
| Physical inference | Pressure anisotropy and collisionality change guide-field-dependent turbulence in measured ways consistent with or contrary to magneto-immutability. | Active/passive comparisons, fixed-input collisionality trends, uncertainty estimates, and robust signatures in strain, heating, residual energy, and structure. |
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
through pressure anisotropy near the driving scale. It also directly
motivates the collisionality axis: the reported pressure-anisotropy heating
and its interpretation differ between collisionless and collisional regimes.

The proposed study is not a reproduction of either paper. It combines the
weak-guide, large-fluctuation stress of the 2019 motivation with the CGL-LF
closure, heating diagnostics, and collisionality dependence motivated by the
2023 work. Its two collisionless anchor requirements are

```text
weak guide anchor:      R_B = B_rms / B_mean approximately 2
strong guide anchor:    R_deltaB = delta_B_rms / B_mean approximately 0.5
both anchor flows:      M_s approximately 0.5
```

Holding each accepted anchor's forcing and initialization fixed while
increasing collisionality then tests how the turbulent state itself changes.
The study tests whether the same physical organization survives when a global
mean-field decomposition is no longer dominant and whether collisions restore,
erase, or reorganize the measured signatures.

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
heating, heat-flux accounting, spectra, PDFs, scale-dependent structure, and
comparison of collisionless with collisional behavior. Their mean-field setup
does not establish the weak-guide limit. The present paper must consequently
make three comparisons separately: whether AthenaK recovers closure-consistent
numerical behavior, whether weak- and strong-guide anchors differ, and how
increasing collisionality alters each anchor's turbulence.

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

The stronger-guide complement is defined by fluctuation amplitude rather than
total RMS field:

```math
R_{\delta B} = \delta B_{\rm rms}/B_{\rm mean}
             = 0.50 \mathbin{+/-} 0.05.
```

For that anchor, report `R_B` as well; in the ideal conserved-mean-field
relation its expected value is `sqrt(1 + R_deltaB^2) approximately 1.118`.
The two anchors therefore isolate the consequences of substantial field
wandering versus an ordered-guide regime while using the same declared
magnetic-amplitude conventions.

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

At weak-guide `R_B = 2`, `beta_mean = 4 beta_rms`. The primary anchors should
target late-time `beta_rms approximately 10`: the weak-guide initialization
therefore begins calibration near `beta0 approximately 40`, whereas the
stronger-guide initialization begins near
`beta0 approximately 10(1 + 0.5^2) = 12.5`. These starting values are
calibration hypotheses, not assumed outcomes. This keeps both anchors in the
same sampled-field beta regime without confusing mean-field beta with the
dynamically sampled rms-field beta.

### 3.1 Collisionality axis

The code input `nu_coll` is a frequency. Define the comparison coordinate
using the measured collisionless turnover time of each accepted guide-field
anchor `G`:

```math
\mathcal{C}_G = \nu_{\rm coll}\tau_{{\rm eddy},0,G},\qquad
\nu_{\rm coll}(\mathcal{C}_G,G)
= {\mathcal{C}_G\over\tau_{{\rm eddy},0,G}}.
```

Here `tau_eddy,0,G` is measured after the collisionless anchor for guide
regime `G` passes its target and stationarity gates. Freeze that value before
generating the collisional cases; also report the a posteriori live value
`nu_coll * tau_eddy` in every collisional analysis window. Use the sequence

```text
C_G = 0, 0.1, 1, 10, 100.
```

The `C_G = 0` anchors set `nu_coll = 0.0`. The nonzero sequence spans weak
collisional modification, turnover-scale competition, and strongly
collisional response. The principal collisionality comparison holds each
anchor's imposed field, forcing, thermodynamics, LF scale, and seed fixed:
drift in `R_B`, `R_deltaB`, `M_s`, or beta is therefore a result rather than a
calibration failure. Conditional state-matched follow-ups may be added only as
separately labeled secondary comparisons.

## 4. Required observables

The following dimensionless controls and responses must be present in every
production manifest or analysis bundle.

| Quantity | Definition or purpose |
| --- | --- |
| Guide regime `G`, `R_B`, `R_deltaB` | Distinguish the weak-guide `R_B` target from the stronger-guide fluctuation target and report both measures for both anchors. |
| `M_s` | Requested subsonic target using the frozen scalar convention above. |
| `beta_mean`, `beta_rms` | Separate imposed-field and sampled-field thermodynamic regimes. |
| `nu_coll`, `C_G`, and live `nu_coll tau_eddy` | Record the prescribed collision coordinate and the realized collision/turnover competition. |
| `M_A,mean`, `M_A,rms` | Alfvenic amplitudes against mean and rms fields; prevent an imprecise use of "trans-Alfvenic." |
| `Delta p = p_perp - p_parallel` and `beta Delta` | Pressure-anisotropy response and instability proximity. |
| `\hat{b}\hat{b}:\nabla u` | Direct magneto-immutability diagnostic motivated by the 2019 paper. |
| `-Delta p \hat{b}\hat{b}:\nabla u` | Anisotropic pressure-work/heating density. |
| LF heat-flux work/contraction and counters | Closure behavior and numerical safety. |
| Mirror/firehose/hard-wall fractions | Extent to which inferred physics depends on modeled microinstability bounds. |
| `min(|B|/B_mean)` and `f_low(q)` for `q = 0.01, 0.05, 0.10` | Quantify exposure of local-field operations to near-zero field magnitude. |
| Field-direction floor/regularization counters | Detect any analysis or evolution result that depends on undefined or regularized `\hat{b}`. |
| Forcing work and energy residual | Global accounting gate. |
| Magnetic/kinetic residual energy | Direct comparison to a major 2019 signature. |
| Isotropic shell spectra | Primary scale-space measure in a weak-guide cubic box without presupposing a global preferred direction. |
| Local-field structure functions; global-mean-field spectra as secondary | Scale-dependent anisotropy and a stronger-guide bridge, with the geometry of each measure kept explicit. |

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
limiter handling, `nu_coll = 0.0`, and an analysis interval at `t = 8` to
`10`. It is not an input suite for either new anchor or for a collisionality
scan.
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

1. Add weak- and stronger-guide input templates under `inputs/cgl_lf_paper/`
   for calibration, active production, and matched passive-feedback control.
   Generate collisional variants from a recorded manifest containing
   `nu_coll`, `C_G`, and the frozen collisionless turnover time used to compute
   it. Use a periodic cubic domain and a verified isotropic solenoidal
   Ornstein-Uhlenbeck forcing configuration. Do not use forcing tied to a
   fixed perpendicular direction unless the paper explicitly treats it as a
   separate mean-field bridge case.
2. Add history outputs or deterministic snapshot reductions for
   `B_mean`, `B_rms`, `delta_B_rms`, `u_rms`, `M_s`, both beta definitions,
   both Alfvenic Mach definitions, residual energy, prescribed and realized
   collisionality, low-field fractions and regularization counters, and the
   selected magneto-immutability/heating diagnostics.
3. Extend `scripts/analyze_cgl_lf_paper.py` with guide-regime target
   attainment, stationarity, fixed-input collisionality trends,
   active/passive comparison, uncertainty, slice, isotropic-shell-spectrum,
   and convergence products. Each product must archive its formula and
   analysis window.
4. Add workflow case names, calibration-only versus held-out science seeds,
   and manifests without weakening the existing
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
| Low-field robustness | A reduced run traverses both calibrated field-amplitude ranges without unaccounted floors, invalid field-direction operations, NaNs, or strict-bound violations; its `f_low(q)` products are reproducible. |
| Active/passive identity | Matching initial state and forcing manifest, differing only in declared pressure-feedback physics. |
| Fixed-input collision identity | Cases along a `C_G` scan differ from their anchor only by declared `nu_coll` and case metadata; a state-matched follow-up is separately named. |
| Calibration separation | The inputs selected with calibration seed `161803` are frozen before any held-out science seed is analyzed. |
| Restart identity | A split guide-field/collisionality case reproduces an uninterrupted reduced reference within established tolerances. |
| GPU sizing | A short representative cubic run records memory, output size, runtime, and LF timestep behavior before a `192^3` or larger production submission. |
| Energy accounting | Active production intervals close the forcing/energy/pressure-work accounting within a declared tolerance validated on reduced tests. |

### 5.4 Implementation and artifact map

Do not reimplement analysis products that are already present. The current
analyzer supplies snapshot pressure-anisotropy PDFs, local parallel strain,
pressure-work decomposition, pressure-transfer products, kinetic/magnetic and
pressure spectra, a heat-flux transport proxy, and opt-in local-field eddy
anisotropy. The expanded campaign needs to expose those products through new
guide-field and collisionality case families and add only the target,
low-field, stationarity, visualization, and comparison layers required by
this experiment.

| File or new artifact | Required change | Acceptance evidence |
| --- | --- | --- |
| `inputs/cgl_lf_paper/cgl_lf_paper_weak_guide_active.athinput` | Add a cubic active base deck with isotropic forcing and placeholders fixed by calibration. | Parsed input metadata matches the contract in Section 6.1. |
| `inputs/cgl_lf_paper/cgl_lf_paper_weak_guide_passive.athinput` | Add the matched passive-feedback base deck. | A test confirms it differs from active only in feedback-declaring controls. |
| `inputs/cgl_lf_paper/cgl_lf_paper_strong_guide_active.athinput` and `..._passive.athinput` | Add matched stronger-guide decks with `R_deltaB` calibrated to `0.5`. | Parsed inputs differ from weak-guide cases only in calibrated anchor parameters and declared metadata. |
| Collisional case manifests/overrides | Instantiate each anchor at frozen `C_G` values through explicit `nu_coll` overrides. | Automated comparison proves the fixed-input scan changes only collisionality metadata and `nu_coll`. |
| `src/pgen/tests/cgl_lf_paper.cpp` and history support | Retain vector mean field, `B_rms`, `delta_B_rms`, volume-weighted `u_rms`, scalar `M_s`, both beta values, residual energy, live collisionality, and low-field metrics/counters at history cadence. | Synthetic or reduced-run regression evaluates each formula and preserves it over restart. |
| `scripts/analyze_cgl_lf_paper.py` | Reuse existing turbulence products; add two-anchor target windows, collisional trends, low-field checks, block statistics, seed/resolution comparisons, isotropic shell spectra, and an acceptance JSON product. | Reduced synthetic and snapshot tests exercise pass, fail, off-target, collisional, and nonstationary cases. |
| `scripts/plot_cgl_lf_paper.py` | Add target/collisionality histories, convergence plots, and two-dimensional field/derived-field slice rendering with predeclared comparable limits. | Figure-generation regression lists all required output files. |
| `scripts/cgl_lf_workflow.py` | Add `paper-guide-calibration`, `paper-guide-collisionality-production`, and `paper-guide-collisionality-analyze` case families and manifest metadata. | Workflow tests enumerate the fixed matrix, enforce held-out seeds, and refuse expensive execution without authorization. |
| Production allocation/accounting path | Provide reviewed production-QOS reservation/submission/accounting outside `scripts/frontier/cgl_lf_frontier.py`. That existing utility is debug-only and deliberately rejects paper production. | Offline budget self-test plus retained production reservation records; the debug-only rejection remains passing. |
| `docs/cgl_lf_validation.tex` | Convert the note into the manuscript structure and figures specified below. | Warning-free PDF build plus figure/data provenance table. |

## 6. Simulation campaign

### 6.1 Physical setup

Both primary guide-field suites should use:

| Choice | Baseline |
| --- | --- |
| Domain | Periodic cubic box, `L_x = L_y = L_z = L = 1`. |
| Mean field | Uniform `B_mean` along `z`; its amplitude is calibrated separately for each guide regime, not presumed. |
| Mesh | `N x N x N`; choose the meshblock decomposition after the GPU sizing run without changing `N` or physical parameters. |
| Forcing mode | `<turb_driving>/driving_type = 0`, the implemented isotropic solenoidal Ornstein-Uhlenbeck velocity-forcing mode. |
| Forcing shell | `physical_k_shell = true`, `k_shell_unit = 2*pi/L = 6.283185307179586`, `nlow = 1`, `nhigh = 3`, `isotropic_power_spectrum = true`, `expo = 2.0`. |
| Seeds | Calibration-only seed `161803`; primary science seed `271828`; second science seed `314159`. No science seed may be used to tune anchor parameters. |
| Forcing amplitude/time | Calibrate `dedt`; start the bracket from the currently exercised `dedt = 0.32`. Freeze `tcorr` after requiring it to be within a factor of two of the measured outer-scale turnover time; start from `tcorr = 2.0`. |
| Anchor targets | Late-time `M_s = 0.50 +/- 0.05` and `beta_rms approximately 10` for both; weak guide uses `R_B = 2.0 +/- 0.2`, stronger guide uses `R_deltaB = 0.50 +/- 0.05`. |
| Thermodynamics | `rho0 = 1`; begin weak-guide calibration near `beta0 = 40` and stronger-guide calibration near `beta0 = 12.5`, using `p_parallel0 = p_perp0 = 0.5 * beta0 * b0^2`; freeze each accepted value. |
| Closure | Active CGL-LF with currently qualified hard-wall instability treatment and strict diagnostics. |
| Matched control | Passive-pressure-feedback case with identical forcing realization and thermodynamic initialization. |
| LF scale | Baseline `lf_k_parallel = 2*pi/L = 6.283185307179586`, interpreted as a closure assumption in a wandering-field box; active sensitivity cases use `pi/L` and `4*pi/L`. |
| Collisionality | Calibrate anchors at `nu_coll = 0.0`; compute nonzero `nu_coll` from the frozen `C_G` grid in Section 3.1 without retuning the fixed-input scan. |
| Limiter/numerics | `cgl_firehose_threshold = parallel`, `limiter_hardwall = true`, strict LF admissibility, applied pressure-work histories, forcing-work history, and explicit metadata for integrator/CFL/STS controls. |
| Time interval | Provisional `t_end = 12`; use `t = 0` to `8` for spin-up and `t = 8` to `12` for analysis only if the stationarity and turnover-span criteria below pass. |

The anchor labels and target contract are:

| Label | Magnetic target | Scientific role |
| --- | --- | --- |
| `wg` | `R_B = 2.0 +/- 0.2`; report `R_deltaB` | Large field wandering and the most demanding local-field geometry. |
| `sg` | `R_deltaB = 0.50 +/- 0.05`; report `R_B` | Ordered-guide complement with the same forcing geometry and sampled-field beta goal. |

The cubic, isotropically forced baseline follows the structural motivation in
Squire et al. (2019), where local-field anisotropy can be measured without
imposing it through elongated outer-scale forcing. It differs deliberately
from the existing elongated, strong-guide MKS24-oriented deck. An elongated
mean-field bridge may be added only as a separately budgeted supplementary
study and must not substitute for either primary anchor.

### 6.2 Calibration procedure

The targets cannot be selected reliably from initial conditions alone because
the saturated fluctuation amplitude depends on forcing, pressure feedback,
limiter activity, and heat flux. Calibrate the two collisionless anchors
before increasing resolution or assigning collisional frequencies:

1. At `32^3`, run active `wg` and `sg` brackets with calibration-only seed
   `161803` and `nu_coll = 0.0`. These tuning runs must never enter science
   averages or uncertainty estimates.
2. Interpolate or extend each bracket only after plotting late-time magnetic
   target, `M_s`, and `beta_rms`. Select finalists by these declared metrics,
   not by visual field appearance.
3. At `64^3`, rerun the two nearest candidates for each anchor through at
   least `t = 8`, including a matched passive case for each best active
   candidate. Adjust `p0`/`beta0` only with calibration-seed data to place
   the selected sampled-field beta regime after the magnetic and velocity
   targets are plausible.
4. If either best active candidate is outside its magnetic or sonic tolerance
   at `64^3`, perform a locally bracketed calibration round and charge it to
   the reservation. Do not expose either held-out science seed or advance
   that anchor's ladder.
5. Freeze one collisionless parameter set per guide regime, including its
   measured `tau_eddy,0,G`, before generating `nu_coll` values. Freeze the
   forcing algorithm, shell, correlation time, target definitions, output
   cadence, closure settings, and science seeds in a versioned manifest.
6. Use science seed `271828` for the first reportable collisionless cases.
   If a frozen anchor drifts outside target tolerance with resolution, retain
   and report that result; any recalibration is a new, explicitly versioned
   campaign and cannot overwrite the held-out observation.
7. For every nonzero `C_G`, reuse the frozen anchor input and change only
   `nu_coll` plus identifying metadata. Its altered `R_B`, `R_deltaB`, `M_s`,
   and beta are measured collisional responses, not failures to match the
   collisionless anchor.

The first collisionless active bracket is concrete but revisable based on its
measured response:

| Calibration variable | Weak guide `wg` | Stronger guide `sg` | Reason |
| --- | --- | --- | --- |
| Resolution/duration/seed | `32^3`, `t_end = 4`, `161803` | `32^3`, `t_end = 4`, `161803` | Cheap response map excluded from science inference. |
| `b0` | `0.25`, `0.50`, `1.00` | `1.00`, `2.00`, `4.00` | Bracket fluctuation amplitude relative to the imposed field. |
| `dedt` | `0.08`, `0.32`, `1.28` | `0.08`, `0.32`, `1.28` | Bracket velocity and magnetic fluctuations around exercised forcing power. |
| `beta0` | Begin at `40` | Begin at `12.5` | Starting hypotheses for `beta_rms approximately 10`. |
| `nu_coll` | `0.0` | `0.0` | Collisionality is varied only after collisionless anchors are frozen. |
| Number/cost | Nine active cases | Nine active cases | Eighteen cases cost approximately `0.770370` node-hours under conservative scaling, leaving confirmation, sizing, and local refinements inside the `80` node-hour reservation. |

Each accepted anchor manifest must record the selected `b0`, `dedt`, `tcorr`,
`beta0`, derived initial pressure, `lf_k_parallel`, calibration seed and
excluded data, frozen science seeds, measured `tau_eddy,0,G`, generated
`nu_coll(C_G,G)`, and the calibration data used to select the anchor.

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
`analysis/guide_collisionality_acceptance.json` with the following decisions:

| Test | Passing definition |
| --- | --- |
| Collisionless anchor magnetic target | For `wg`, full-window mean `R_B` lies in `[1.8, 2.2]`; for `sg`, full-window mean `R_deltaB` lies in `[0.45, 0.55]`; report both quantities for both anchors. |
| Collisionless anchor sonic target | Full-window mean `M_s` lies in `[0.45, 0.55]` for both `G` at `C_G = 0`. |
| Collisional response reporting | For `C_G > 0`, report departure from its fixed-input anchor in magnetic amplitude, `M_s`, beta, and live `nu_coll tau_eddy`; do not reject a case merely because collisions alter those values. |
| Stationarity | For the anchor magnetic metric, `M_s`, and principal heating measures, the difference between first and last block means is below the pre-registered tolerance (`0.10` for `wg R_B`, `0.025` for `sg R_deltaB` and `M_s` initially). Expand this table before execution if a collisional-response metric needs its own tolerance. |
| Energetic stationarity | Conserved total-energy change versus recorded forcing work is consistent with the reduced-test accounting tolerance or a tighter production tolerance frozen before the ladder; applied closure-work histories are reported separately to interpret internal transfers. |
| Safety | Zero nonfinite, nonpositive, floor, emergency-bound, or invalid-direction failures; report hard-wall and cap activity as physics/model exposure rather than silently excluding it. |
| Low-field interpretability | Report `min(|B|/B_mean)` and `f_low(q)` for `q = 0.01, 0.05, 0.10`; no `\hat{b}` floor or regularization event is permitted. If the accepted-window `f_low(0.05) > 10^-3`, retain the data but withhold an unqualified local-closure physical claim pending a reviewed sensitivity study. |
| Sample sufficiency | Four turnover times and at least four block summaries; bootstrap intervals are emitted for paper statistics. |

If `t = 8` to `12` is too short to establish those conditions at `96^3`,
stop the high-resolution ladder and spend only the reserved contingency on a
reviewed time-extension plan. An off-target or nonstationary case is retained
as a failed calibration/result or measured collisional response, not omitted
from the record. The low-field threshold above is a conservative pre-analysis
interpretability rule, not a claim that a larger fraction is necessarily
incorrect.

### 6.4 Resolution ladder and matched cases

The pre-registered core matrix uses forcing-matched active/passive pairs. Its
collisionless anchors establish resolution behavior in both guide regimes;
its collisional sweep establishes trends at useful intermediate resolutions;
and its strongly collisional high-resolution endpoints test whether the
limiting contrast is resolved.

| Tier | Guide regimes | `C_G` | Resolution(s) | Models/seeds | Purpose |
| --- | --- | --- | --- | --- | --- |
| Calibration | `wg`, `sg` | `0` | `32^3`, then `64^3` | active finalists plus confirmation passive; `161803` only | Select collisionless anchors; never enter science averages. |
| Collisionless anchor ladder | `wg`, `sg` | `0` | `64^3`, `96^3`, `192^3`, `384^3` | active/passive; `271828` | Anchor convergence, principal slices, and maximum-resolution comparison. |
| Collisionless seed repeat | `wg`, `sg` | `0` | `96^3`, `192^3` | active/passive; `314159` | Realization uncertainty without high-resolution duplication. |
| Collisionality sweep | `wg`, `sg` | `0.1`, `1`, `10`, `100` | `96^3`, `192^3` | active/passive; `271828` | Fixed-input collisional response curves. |
| Strongly collisional endpoint | `wg`, `sg` | `100` | `384^3` | active/passive; `271828` | High-resolution endpoint for the largest collision contrast. |
| Collisional seed repeat | `wg`, `sg` | `1`, `100` | `192^3` | active/passive; `314159` | Sampling check at turnover-scale and strongly collisional endpoints. |

The high-resolution cases remain paired with passive controls because a
single active result cannot isolate pressure-feedback physics. Intermediate
resolutions carry the denser `C_G` and seed matrix because they quantify trend
and variance without exhausting the high-resolution reserve.

The production manifest should use stable case identities:

| Case pattern | Included values |
| --- | --- |
| `g{wg|sg}_c0_{active|passive}_N{N}_s271828` | `N = 64, 96, 192, 384`; collisionless anchor ladder. |
| `g{wg|sg}_c0_{active|passive}_N{N}_s314159` | `N = 96, 192`; collisionless anchor repeat. |
| `g{wg|sg}_c{0p1|1|10|100}_{active|passive}_N{N}_s271828` | `N = 96, 192`; nonzero-collisionality sweep. |
| `g{wg|sg}_c100_{active|passive}_N384_s271828` | High-resolution strongly collisional endpoint. |
| `g{wg|sg}_c{1|100}_{active|passive}_N192_s314159` | Collisional seed repeats. |
| `g{wg|sg}_c{0|1}_active_N96_s271828_kpar_{half|double}` | Active-only LF-scale sensitivity cases. |
| `g{wg|sg}_c{C}_passive_state_N{N}_s271828` | Conditional state-matched passive controls only if triggered and budget-reviewed. |

Each active/passive pair must reuse an identical archived OU seed and forcing
configuration. It is an implementation failure if the passive run silently
uses a separately generated random realization.

The target gate is mandatory for the collisionless active anchors. All passive
cases and all nonzero-`C_G` cases in the fixed matrix are forcing-matched or
fixed-input comparisons: they answer how changing feedback or collisionality
changes a calculation initialized from the same declared anchor. They must
report their attained state, but they must not be called state-matched
controls.

If attained states differ enough to make a fixed-input contrast difficult to
interpret, use the contingency allocation for explicitly secondary
state-matched controls at `96^3` or `192^3`, adjusting only declared
calibration controls and archiving the adjustment. Use fixed-input
comparisons for causal response to the changed model parameter and
state-matched comparisons for conditional turbulence statistics. Do not add a
state-matched `384^3` case without recomputing and approving the budget.

### 6.5 Secondary beta and LF-scale checks

The primary sampled-field regime is `beta_rms approximately 10`. Because the
strength of anisotropy feedback and collisional relaxation may depend on beta,
reserve targeted secondary tests at `beta_rms approximately 1` and `100`:

| Check | Cases | Interpretation |
| --- | --- | --- |
| Collisionless beta span | Both guide regimes, both secondary beta values, active/passive, `96^3` and `192^3`, primary science seed. | Determine whether guide-field conclusions persist across weak and strong anisotropy-feedback regimes. |
| Turnover-collisional beta span | Both guide regimes, both secondary beta values, `C_G = 1`, active/passive, `96^3`, primary science seed. | Check whether the central collisional transition is qualitatively beta-dependent. |
| LF-scale sensitivity | Both guide regimes at primary beta, `C_G = 0` and `1`, active-only `96^3` cases at `lf_k_parallel = pi/L` and `4*pi/L`. | Bound closure-scale dependence for anchor and transitional collisional states. |

Tune secondary-beta collisionless starting thermodynamics with
calibration-only data and freeze them before opening science outputs. For
each `C_G = 1` secondary-beta case, hold its corresponding collisionless
secondary-beta input fixed except for `nu_coll`. These cases are secondary and
cannot rescue a failed primary-beta claim.

## 7. Node-hour budget and maximum admissible resolution

The hard ceiling is `4000` Frontier node-hours, including the existing
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
| Implementation/sizing/calibration reservation | Short GPU sizing plus tuning, confirmation, and secondary-beta refinements using calibration-only data | `80.000000` |
| Collisionless anchor paired ladders | Two guide regimes; active/passive at `64^3`, `96^3`, `192^3`, and `384^3`; primary seed | `1016.375348` |
| Collisionless second-seed pairs | Two guide regimes; active/passive at `96^3` and `192^3` | `124.800004` |
| Collisionality sweep | Four nonzero `C_G` values, two guide regimes; active/passive at `96^3` and `192^3`; primary seed | `499.200016` |
| Strongly collisional high-resolution endpoint | `C_G = 100`, two guide regimes; active/passive at `384^3`; primary seed | `887.466700` |
| Collisional second-seed pairs | `C_G = 1` and `100`, two guide regimes; active/passive at `192^3` | `221.866672` |
| Beta sensitivity | Secondary `beta_rms` values `1` and `100` at `C_G = 0` (`96^3`, `192^3`) and `C_G = 1` (`96^3`), both guide regimes and active/passive models | `277.333344` |
| LF-scale sensitivity | `C_G = 0` and `1`, both guide regimes; two active-only `96^3` variants | `27.733336` |
| Failure/time/state-matching contingency | Conditional state-matched controls, time extensions, or reruns only after a documented gate review | `600.000000` |
| Total committed ceiling | Includes consumed ledger and contingency | `3735.627090` |
| Uncommitted margin below `4000` | Preserved for accounting error or a separately reviewed amendment | `264.372910` |

This allocation assumes that an accepted statistical case ends at `t = 12`.
If the turnover-span or stationarity gate requires longer runs, agents must
recompute every affected case reservation before submitting the extension.
The `600` node-hour contingency is not permission to overrun the ledger
silently; it is the only preallocated source for an approved extension or
rerun.

The largest fully paired, pre-registered production resolution is `384^3`.
Under the same conservative evidence, one `512^3` `t = 12` case costs
`525.906193` node-hours; the required active/passive pair alone costs
`1051.812386` node-hours. Replacing any required `384^3` pair with that
`512^3` pair adds more than the present margin and would sacrifice either the
pre-registered comparison matrix or contingency protection. A measured cubic
benchmark may later justify a reviewed amendment, but agents must not
schedule a `512^3` pair on an assumed cell-count speedup.

The existing `scripts/frontier/cgl_lf_frontier.py` utility is not the
submission path for this suite: it is an intentionally debug-only
qualification tool and rejects paper-production decks. Before any production
submission, either implement a separately reviewed production accounting
driver with the same fail-closed `4000` node-hour ledger or retain equivalent
manual reservation records approved for the production queue. Do not remove
the debug rejection to make the new cases executable.

Before every submitted job:

1. Regenerate the ledger summary from the workflow tooling.
2. Reserve the job's conservative upper bound, including queue/runtime limit.
3. Refuse submission if consumed plus active reservations plus the new ceiling
   exceeds `4000.000000` node-hours.
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
expected cubic reduction. The fixed matrix contains eight `384^3` cases:
active/passive for both guide regimes at `C_G = 0` and `100`. Reserve at
least `1.0 TB` for those conservative high-resolution products alone and
obtain a `1.5 TB` campaign reservation before production to cover
lower-resolution products, derived arrays, and metadata unless measured cubic
sizing supports a documented reduction. Confirm actual cubic file sizes in
the sizing job, record the storage reservation alongside node-hours, and
prune only according to an archived retention rule.

Pre-register figure snapshot selection before reviewing visual candidates:
use the saved snapshot nearest the midpoint of each accepted analysis window,
use central coordinate planes for the principal slice grid, and set shared
color limits from pooled, declared percentiles of the cases being compared.
Additional event-focused views may be shown only as labeled supplementary
figures selected by a recorded quantitative criterion.

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
| Time histories of `R_B`, `R_deltaB`, `M_s`, `beta_rms`, and injected/total energy | Each collisionless anchor reaches its declared subsonic state and collisional cases expose any fixed-input state drift. |
| Anchor target planes and collisional response curves by guide/resolution/model/seed | Anchor attainment and subsequent `C_G` trends are not artifacts of one resolution or realization. |
| Energy-work residual, hard-wall fraction, LF counters, low-field fractions, and extrema | The state is numerically admissible rather than sustained by hidden failure handling or undefined field directions. |

### 9.2 Field-structure visuals

At the pre-registered midpoint snapshot and central planes from Section 8,
render matched active/passive and selected collisionless/collisional slices
with pooled percentile limits recorded in the figure provenance:

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
| PDFs and joint PDFs | Within each guide regime, compare active/passive and `C_G` trends in `beta Delta`, `|B|/B_mean`, strain, pressure work, and their joint relation. |
| Energy and residual-energy statistics | Determine how active feedback and increasing collisions alter magnetic dominance and heating for `wg` versus `sg`. |
| Primary isotropic shell spectra | Kinetic, magnetic, density, `Delta p`, pressure, strain, and pressure-work spectra in the common cubic geometry; compare model, `C_G`, and resolved-scale convergence without presuming a global guide direction. |
| Secondary guide-relative spectra | Show global-mean-field `k_parallel/k_perp` products only as a geometry-labeled complement, most interpretable for `sg`; do not use them as the sole weak-guide structural diagnostic. |
| Scale-dependent heating/transfer | Determine where injected energy is removed by anisotropic pressure work and LF transport as collisionality increases, motivated by the 2023 results. |
| Local-field structure functions | Measure scale-dependent anisotropy without relying on the weak global guide field; follow local-field conditioning in the 2019 analysis. |
| Resolution and seed convergence | Plot anchor and `C_G = 100` statistic differences from the `384^3` cases and intermediate-`C_G` realization uncertainty. |
| Beta and LF-scale sensitivity | Show whether primary conclusions persist across sampled-field beta and factors-of-two changes in `lf_k_parallel`. |

For external-paper comparisons, dimensionless definitions with demonstrably
matching conventions may be discussed directly. Absolute spectral ordinates
must remain excluded from quantitative cross-paper overlays unless the
normalization/provenance gap recorded in the MKS24 implementation plan is
resolved.

### 9.4 Statistical protocol

1. Establish each case's late-time window from declared stationarity tests
   before examining preferred physical comparisons; record any window
   differences when making a pairwise contrast.
2. Use turnover-time block averages and block bootstrap intervals for scalar
   measurements and binned diagnostics.
3. Use the held-out science seeds at collisionless anchors and at
   `C_G = 1, 100` to distinguish realization variance from resolution and
   collisionality trends; never add a posteriori seed selection.
4. Compare spectra only over overlapping resolved ranges; make isotropic
   shells the primary weak-guide spectral result and record the
   de-aliasing/binning/FFT normalization metadata.
5. Report fixed-input contrasts as primary; label any separately tuned
   state-matched control as secondary and do not pool the estimands.
6. Treat failure of convergence, anchor target attainment, energy closure,
   low-field validity, beta sensitivity, or LF sensitivity as a scientific
   result that limits the claim.

### 9.5 Quantitative decision and artifact contract

Before inspecting the `384^3` case, freeze the principal tests in the
production manifest:

| Principal statistic | Pre-registered comparison |
| --- | --- |
| Parallel-strain suppression | `D_fb(G,C)`: active/passive ratio of RMS `\hat{b}\hat{b}:\nabla u` for matched `G` and `C_G`. |
| Magnetic-strength regulation | `D_fb(G,C)`: active/passive comparison of distribution width of `|B|/B_mean`, while testing anchors separately. |
| Pressure-feedback energetics | Active anisotropic pressure work and LF work normalized by recorded injected work as functions of `G` and `C_G`. |
| Residual energy | `D_fb(G,C)`: active/passive magnetic-minus-kinetic residual-energy comparison. |
| Collisional response | `D_C(G,C) = X_active(G,C) - X_active(G,0)` at fixed anchor inputs, for each principal statistic. |
| Guide-field dependence | `D_G(C) = X_active(wg,C) - X_active(sg,C)` with matched declared beta regime and collision coordinate. |
| Limiter/closure/low-field exposure | Mirror, firehose, hard-wall, LF-cap, field-direction regularization, and `f_low(q)` fractions plus their beta/LF-scale sensitivity. |

An active/passive effect at a collisionless anchor may be described as
resolution-robust only if its sign is the same for both seeds at `96^3` and
`192^3`, is retained by the `384^3` primary-seed comparison, and the
normalized `192^3` to `384^3` change is no greater than the larger of `10%`
or the measured two-seed uncertainty at `192^3`. Apply the same endpoint
criterion to `C_G = 100`, using its second-seed `192^3` check and primary-seed
`384^3` comparison. Intermediate collisional trends may be described as
resolved at `192^3` only when their sign/order survives overlapping
uncertainties and safety gates; they are not high-resolution-converged by
assumption. Limiter, low-field, beta, and LF-scale exposure qualify all
physical claims. Apply equivalent criteria to spectral or structure-function
statements only on the common resolved range.

Every accepted production bundle must retain these machine-readable and
rendered products:

| Artifact | Contents |
| --- | --- |
| `analysis/guide_collisionality_acceptance.json` | Anchor targets, collisional state drift, block trends, turnover-span decision, low-field/safety counters, energy gate, and pass/fail reasons. |
| `analysis/guide_collisionality_statistics.json` | Block/bootstrap estimates, active/passive effects, `D_C`/`D_G` contrasts, seed scatter, resolution comparisons, and beta/LF-scale sensitivity. |
| `analysis/guide_collisionality_provenance.json` | Executable/input hashes, parent restart, case identity, calibration exclusion, seed, `C_G`/`nu_coll`, model choices, output cadence, analysis configuration, and ledger entry. |
| `figures/guide_collisionality_targets.pdf` | Anchor target/stationarity histories and fixed-input collisional response curves. |
| `figures/guide_collisionality_slices.pdf` | Preselected common-scale guide/feedback/collision slices of fields, anisotropy, strain, pressure work, and limiter/LF exposure. |
| `figures/guide_collisionality_pdfs.pdf` | PDFs and joint PDFs used for magneto-immutability and collisionality inferences. |
| `figures/guide_collisionality_energetics.pdf` | Injection, applied pressure/LF work, residual energy, collision trend, and accounting quality. |
| `figures/guide_collisionality_spectra.pdf` | Isotropic spectra, explicitly secondary guide-relative spectra, and resolved-range convergence markers. |
| `figures/guide_collisionality_structure.pdf` | Local-field structure functions/anisotropy and resolution/collisionality comparison. |
| `figures/guide_collisionality_sensitivity.pdf` | Seed, resolution, beta, LF-scale, and safety/limiter/low-field exposure summary. |

## 10. Manuscript conversion plan

The current `docs/cgl_lf_validation.tex` is an illustrated method-validation
note. Convert it into a complete paper only after the guide-field and
collisionality diagnostics and at least the calibration/qualification
evidence exist. Preserve useful method material, but restructure the argument
around the physical uncertainty and the new tests.

| Manuscript section | Required content and evidence |
| --- | --- |
| Abstract | Begin with uncertainty in collisional CGL-LF turbulence across guide-field strength; state only tested regimes and supported results; end with the principal limitation. |
| Introduction | Explain magneto-immutability, why weak guide fields challenge closure/numerics, and why a collisionless-to-collisional comparison is needed; position the study relative to the 2019 Braginskii and 2023 CGL-LF works. |
| Model and regime parameters | CGL-LF equations, limiter, heat-flux and collisional-relaxation model, `R_B`, `R_deltaB`, `M_s`, beta definitions, `C_G`, Alfvenic measures, and any interruption-related measure with closure-specific qualification. |
| Numerical method and validation | Condense existing operator tests, hard-wall/energy-work corrections, restart/GPU qualification, and new low-field and fixed-input-collision tests. |
| Turbulence experiment | Cubic forcing, two anchor calibrations, held-out seeds, fixed-input `C_G` matrix, active/passive pairing, resolution ladder, secondary beta/LF-scale checks, analysis window, and budget/provenance. |
| Results I: attained states | Demonstrate stationary collisionless anchors, low-field/safety/energy gates, and measured collisional state drift before presenting physical claims. |
| Results II: guide and collisional structures | Slices, PDFs, strain response, residual energy, limiter occupancy, and active/passive/guide/collision contrasts. |
| Results III: energetics and cascade | Pressure/LF work, collisional heating trends, isotropic spectra, transfer/heating scale dependence, local-field structure functions, convergence and sensitivity. |
| Discussion | State what is robust, how it agrees or differs from the source motivations, and what remains model-, limiter-, or closure-scale-dependent. |
| Conclusions | Make the strongest defensible claim; do not replace validation with a perfection claim. |
| Reproducibility appendix/data statement | Case table, hashes, ledger, definitions, analysis commands, retained-data location, and figure-to-product mapping. |

Planned main-text figure sequence:

| Figure | Purpose |
| ---: | --- |
| 1 | Campaign design and attained anchor states: guide targets, `M_s`, beta, `C_G`, active/passive matrix, and resolution ladder. |
| 2 | Anchor and collisional numerical-control histories, including low-field and limiter exposure. |
| 3 | Preselected matched slices across guide regime and collision endpoints for `|B|`, anisotropy, strain, and pressure work. |
| 4 | PDFs/joint PDFs showing whether field-strength-changing strain suppression changes with guide field and collisionality. |
| 5 | Principal active/passive, collisional, and guide-field estimands with seed and resolution uncertainty. |
| 6 | Injection, pressure/LF work, residual energy, and collisional heating trends. |
| 7 | Isotropic kinetic, magnetic, pressure, and strain spectra with resolved-range convergence. |
| 8 | Heating/transfer scale dependence and local-field structure functions/anisotropy. |
| 9 | Beta, LF-scale, limiter, and low-field validity sensitivity summary. |

### 10.1 Mapping the existing TeX note into the paper

Retain `docs/cgl_lf_validation.tex` as the editable manuscript source unless
a reviewed journal-template conversion requires a new filename. The existing
content should be transformed as follows:

| Current note content | Manuscript treatment |
| --- | --- |
| Purpose and implementation narrative | Replace with a styled abstract and introduction beginning from guide-dependent collisional-response uncertainty. |
| Conserved variables, LF closure, heat-flux cap, and threshold policy | Retain in a compact methods section; move detailed operator exposition to an appendix if it interrupts the physical argument. |
| Method validation problems and five method figures | Retain as numerical-method evidence, condensed into one validation section plus appendix material as needed. |
| Reduced MKS24-oriented product-path figures | Do not present as physical results. Replace in the main narrative with accepted guide/collisionality figures; retain only where explicitly labeled as diagnostic validation. |
| MKS24 reproduction boundary | Recast as a data/provenance and limitation statement; it is not the central result of the new paper. |
| Current bibliography | Add Squire et al. (2019) and the weak-guide interpretation, retain Squire et al. (2023) and its collisional comparison motivation, and add directly used method/diagnostic references. |

### 10.2 Publication and data mechanics

Keep the manuscript in the current article-format TeX source while the result
structure and figures change rapidly. Introduce a tracked
`docs/cgl_lf_validation.bib` when the manuscript conversion begins and move
its citations through a tested BibTeX build; migrate to a journal class only
after the main text, figures, and references build reliably in the current
format. Track manuscript TeX, bibliography, compact plotted source tables,
figure PDFs, and provenance summaries in the repository. Archive raw
high-volume snapshot bundles externally with checksums and a retained-data
manifest referenced by the paper; do not place raw production outputs in Git.

### 10.3 Complete-manuscript acceptance checklist

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
python3 scripts/cgl_lf_workflow.py paper-guide-collisionality-analyze \
  --output-dir /path/to/accepted-guide-collisionality-bundle
cd docs
pdflatex -halt-on-error -interaction=nonstopmode \
  -output-directory /tmp/cgl_lf_manuscript_build cgl_lf_validation.tex
bibtex /tmp/cgl_lf_manuscript_build/cgl_lf_validation
pdflatex -halt-on-error -interaction=nonstopmode \
  -output-directory /tmp/cgl_lf_manuscript_build cgl_lf_validation.tex
pdflatex -halt-on-error -interaction=nonstopmode \
  -output-directory /tmp/cgl_lf_manuscript_build cgl_lf_validation.tex
```

The workflow command and BibTeX conversion are required future interfaces, not
currently available steps; Phase B and the manuscript conversion must
implement and test them before this sequence is used.

## 11. Execution phases and stop conditions

| Phase | Work | Deliverable | Stop condition |
| --- | --- | --- | --- |
| A. Specification | Freeze two anchor definitions, collision coordinate, estimands, forcing choice, manuscript claims, and budget ledger. | Reviewed protocol and input schema. | Definitions remain ambiguous or budget tooling cannot enforce ceiling. |
| B. Implementation | Add guide/collision/low-field histories, analyzer and plot products, input/workflow support, and tests. | Passing local/reduced validation bundle. | Existing numerical gates regress or new observables are not reproducible. |
| C. Calibration | With seed `161803` only, tune collisionless `wg` and `sg` anchors at `32^3`/`64^3`, measure `tau_eddy,0,G`, and generate `nu_coll` values. | Frozen anchor and collisional-scan manifests meeting target tolerances provisionally. | Either anchor cannot be reached without unsafe behavior or unstable transients. |
| D. GPU qualification | Run representative short case and restart/output sizing. | Updated measured runtime/storage reservation. | Runtime, memory, LF behavior, energy closure, or strict counters are unacceptable. |
| E. Anchor ladders | Run paired collisionless `64^3`, `96^3`, `192^3`, then `384^3` anchors only as gates remain satisfied. | Archived collisionless production cases and ledger. | Anchor target failure, nonstationarity, failed accounting, low-field failure, or budget margin breach. |
| F. Collisionality matrix | Run fixed-input nonzero-`C_G` pairs through `192^3`, and `C_G = 100` through `384^3` only after intermediate gates pass. | Guide-dependent collisional-response products. | Safety, stationarity, accounting, low-field, or budget gate fails. |
| G. Robustness | Run declared second seeds, beta variants, LF-scale variants, and justified secondary state matching. | Uncertainty and model-sensitivity products. | Principal inference changes without a qualified explanation. |
| H. Manuscript | Revise TeX, generate figures/tables, write conclusions and caveats. | Buildable complete manuscript and reproducibility package. | Results cannot support stated claims; rewrite conclusions accordingly. |

### 11.1 Current handoff status and first execution sequence

As of 2026-05-25, this repository contains the restored style guide, the
illustrated TeX validation note, existing strong-guide CGL-LF paper-analysis
infrastructure, and reduced GPU qualification/budget evidence. It does not
contain the new weak/strong-guide anchor decks, collisional scan manifests,
guide/collisionality target or low-field products, a campaign slice renderer,
a production submission/accounting path, accepted campaign production runs,
or a completed manuscript for this study.

A future implementation agent should therefore proceed in this order:

1. Implement the target, collisionality, and low-field histories; both anchor
   base inputs; case-family metadata; synthetic target/stationarity tests; and
   active/passive plus fixed-input-collision identity tests without submitting
   an HPC job.
2. Implement slice and quantitative analysis products and exercise them on
   reduced local/GPU qualification data; keep all figures explicitly labeled
   non-statistical until production cases are accepted.
3. Establish reviewed production allocation/accounting that preserves the
   debug-only protection of `scripts/frontier/cgl_lf_frontier.py`.
4. Execute calibration with seed `161803` and short GPU sizing only after the
   local/reduced gates pass; freeze anchors, collision frequencies, physical
   parameters, science seeds, and updated measured cost.
5. Execute collisionless ladders, the collisionality matrix, and declared
   sensitivities conditionally; analyze accepted cases; then convert the TeX
   note using the manuscript checklist above.

Every high-resolution submission is conditional: do not execute `384^3`
merely because it appears in this plan. Collisionless `384^3` cases require
accepted lower-resolution anchors; collisional `384^3` endpoint cases require
accepted lower-resolution fixed-input trends. Both require stationarity,
numerical gates, low-field checks, and an updated measured-cost ledger.

## 12. Reproducibility record future agents must maintain

For each accepted calculation, archive:

1. Git revision, executable build metadata, platform/compiler/GPU information,
   exact submitted input, case role, guide regime, forcing seed, calibration
   versus science designation, and parent restart if any.
2. Consumed and reserved node-hours before and after execution, walltime,
   allocation identifier, storage footprint, and pruning record.
3. Target definitions, `C_G`, `nu_coll`, frozen `tau_eddy,0,G`, realized
   `nu_coll tau_eddy`, and late-time intervals, with the raw histories from
   which target and collisional-response decisions were made.
4. Analysis configuration, JSON products, plotted source data, plot scripts,
   figure files, and checksums.
5. Gate results, including failures, field-direction regularization counters,
   and `f_low(q)` values. Failed cases must remain identifiable and may not be
   quietly replaced by a selected successful realization.

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
| Stronger-guide complement with `delta_B_rms/B_mean approximately 0.5` | Sections 3, 6.1, 6.2, and 6.3 | Acceptance JSON with declared `R_deltaB` definition and passed collisionless anchor target. |
| Sonic Mach number `approximately 0.5` | Sections 3, 6.1, and 6.3 | Acceptance JSON with scalar `M_s` definition and passed late-time target. |
| Collisionless and increasingly collisional turbulent boxes | Sections 3.1, 6.2, 6.4, 9.5, and 12 | Frozen `C_G` scan, fixed-input manifests, collisional-response estimands, and archived live collision measures. |
| Resolution range as large as useful within `4000` node-hours | Sections 6.4 and 7 | Fail-closed ledger, measured updates, and conditional paired `384^3` anchor and endpoint tests. |
| Robust turbulence analysis and strong visuals | Sections 4, 8, and 9 | Preselected common-scale slices, retained JSON products, isotropic spectra, low-field checks, quantitative figures, convergence, and uncertainty tests. |
| Motivation from arXiv:1811.12421 and arXiv:2303.00468 | Sections 2 and 13 | Literature-positioned introduction/discussion and diagnostics tied to those motivations without claiming reproduction. |
| Held-out calibration, beta/LF sensitivity, low-field validity, and publication mechanics | Sections 5, 6.2, 6.3, 6.5, 8, 9, and 10.2 | Frozen calibration provenance, sensitivity artifacts, safety reports, and manuscript/data build record. |
| Lead future agents | Sections 5, 6, 7, 8, 9.5, 10.3, 11, and 12 | File-level work map, cases, gates, budget, retained products, and immediate work order. |
