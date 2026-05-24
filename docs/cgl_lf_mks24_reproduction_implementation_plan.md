# CGL Landau-Fluid Physical-Correctness and MKS24 Reproduction Implementation Plan

## Status and Purpose

This document is an implementation handoff for future agents working on the
AthenaK CGL Landau-fluid (CGL-LF) feature branch. It supersedes neither the
user documentation in `docs/source/` nor the historical planning note
`next_steps_forCGL.md`. The earlier note describes work that has already been
substantially delivered on `feature/cgl-landau-fluid`: hardening, routine CPU
coverage, explicit reference integration, a validation workflow, AMR smoke
testing, restart coverage, and MPI reproducibility coverage.

The remaining work is different in character. It must:

1. Correct confirmed physical-correctness problems in the current LF/collision
   and instability-limiter implementation.
2. Turn the existing operator-level validation framework into a faithful
   reproduction program for Majeski, Kunz, and Squire 2024 (MKS24), whose
   source is available locally in `docs/arXiv-2405.02418v2/MKS24.tex`.
3. Establish a reproducible workflow for paper-grade driven-turbulence runs,
   their analysis products, and their limitations.

This is a plan, not a claim that the work has already been performed.

## Read This First

Future agents should begin every implementation session by confirming the live
branch and worktree status:

```bash
git status --short --branch
git log --oneline --decorate -12
git fetch origin main CGL-STS-LF gh-pages --prune
```

At the time this plan was written, the relevant baseline was:

```text
branch: feature/cgl-landau-fluid
HEAD:   8f7d13c6 Add CGL-LF MPI reproducibility coverage
```

The local paper source directory and the older planning note were untracked:

```text
?? docs/arXiv-2405.02418v2/
?? next_steps_forCGL.md
```

Do not accidentally stage, delete, or rewrite either untracked path unless the
user explicitly asks for it. The paper source directory is a reference input
for this plan.

## Executive Summary

The existing branch is a good foundation for a CGL-LF numerical feature:

- LF transport is routed through a dedicated parabolic split update.
- The conserved `IAN` slot has a protected representation lifecycle:
  anisotropy during ordinary MHD evolution, temporary magnetic moment only
  inside LF split sweeps.
- STS and explicit reference integration both exist.
- LF admissibility counters, strict validation mode, AMR restrictions,
  restart testing, MPI reproducibility testing, and a Python workflow driver
  already exist.

However, the branch is not ready to support paper-reproduction claims.
Three findings dominate all next work:

1. **Collision splitting is physically wrong or, at minimum, unjustified.**
   When LF and collisions are active, each full cycle currently applies two
   collision relaxations, each using the full mesh timestep. The quantitative
   oracle reproduces that schedule, so tests validate the bug rather than
   detect it.
2. **The firehose activation threshold is not the MKS24 production
   threshold.** The current code activates at the oblique threshold
   `beta Delta = -1.4`; MKS24 production simulations use the parallel-fluid
   threshold `beta Delta = -2`.
3. **The current validation workflow does not implement the MKS24 experiment.**
   It validates linear/operator cases, not forced 3D turbulence, active versus
   passive anisotropy, published diagnostic products, or the limiter scan.

The work must therefore proceed in this order:

1. Correct and validate physical update semantics.
2. Make threshold conventions explicit and reproducible.
3. Strengthen diagnostic truthfulness and numerical accuracy tests.
4. Add a current-API paper problem generator and driven-turbulence workflows.
5. Add paper observables and figure-generation analysis.
6. Execute tiered smoke, convergence, standard, and HPC reproduction runs.

Do not begin performance tuning, new operator composition support, Frontier
GPU qualification, or publication-level claims before phases 1 through 3 are
complete. Frontier GPU correctness and operational validation are required
after the physical-core phases and before any production-readiness claim.

## Meaning of Production Readiness

Completing only the original operator, paper-pgen, and analysis tasks would
not by itself justify describing this feature as fully production ready. It
would justify a much narrower statement: that the corrected CGL-LF model has a
credible, reproducible implementation of the scoped MKS24 CGL-LF experiment.

For this project, **production ready** means all of the following:

1. The implemented equations, closure coefficients, limiter conventions, and
   source-term splitting are physically specified, independently tested, and
   consistent with the claimed application.
2. The CPU and Frontier GPU implementations pass accuracy, robustness,
   restart, MPI-decomposition, and reproducibility gates at representative
   configurations.
3. Production workflows fail loudly on invalid numerical state, record enough
   provenance to recreate runs, and can regenerate scientific analysis from
   retained outputs.
4. MKS24 reproduction cases have been executed at the required scientific
   fidelity or are clearly labelled as smoke/reduced validation rather than
   reproduction.
5. Known limitations, unsupported configurations, performance evidence,
   resource requirements, and operational risks are documented for users.
6. Every serious new finding discovered during implementation or campaign
   execution has either been resolved and regression-tested or is an explicit
   release blocker.

This definition is scoped to the supported CGL-LF configurations established
by this plan. It does **not** imply that unimplemented LF composition with
other parabolic operators, primitive-prolongation AMR, kinetic-scale
microinstability physics, or arbitrary boundary conditions become supported
without their own validation programs.

### Production Readiness Gate Table

Future agents must maintain this gate table as implementation proceeds. A gate
cannot be marked complete solely because code exists; it requires retained
test or run evidence.

| Gate | Required evidence | Release-blocking? | Initial status |
| --- | --- | --- | --- |
| Physical equations and normalization | Equation-to-code review, closure/cap tests, documented code-unit mapping | Yes | Partial; collisionless mapping reviewed |
| Collision source-term timing | Independent analytic regression and operator-ordering documentation | Yes | Failing/current behavior suspect |
| Instability-threshold policy | Explicit MKS24 policy, alternative-policy tests, archived input choice | Yes | Failing/current default mismatched |
| Safety diagnostic truthfulness | Strict tests independent of backup correction | Yes | Failing/current hard-bound gap |
| CPU numerical accuracy | Convergence, asymptotic, cap, directional, restart, MPI/AMR checks | Yes | Partial |
| Frontier GPU numerical equivalence | CPU/GPU comparison suite, MPI/GPU restart checks, strict diagnostics | Yes for Frontier use | Not started |
| Paper pgen and forcing fidelity | Input/pgen review, forcing metadata, restartable reduced smoke | Yes for MKS24 claims | Not started |
| Paper observables and analysis | Synthetic analysis tests and archived reduced-run products | Yes for MKS24 claims | Not started |
| Standard MKS24 results | Required cases, durations, manifests, figure comparisons | Yes for reproduction claim | Not started |
| Operational workflow | Frontier scripts, budget ledger, failure recovery, storage plan | Yes for Frontier use | Not started |
| User documentation | Sphinx build and accurately scoped runbook | Yes | Partial |
| Performance suitability | Representative timing/memory/I/O evidence; no uninvestigated prohibitive bottleneck | Yes for production use | Not started |

The final readiness report must distinguish:

- **production ready for supported CGL-LF use**;
- **MKS24 reproduction complete**;
- **validated on Frontier GPU hardware**;
- **future extensions not yet supported**.

None of these phrases should be used as a substitute for the others.

## Living Plan and Revision Protocol

This plan must be treated as a controlled living document, not as a checklist
that remains fixed when contrary evidence appears. The scientific and
operational work is expected to generate new findings. Revising the plan in
response is required engineering work, not scope drift.

### Required Iteration Loop

For each implementation phase, test campaign, or Frontier execution campaign,
future agents must follow this loop:

1. **Restate the hypothesis.** Identify the physical, numerical, operational,
   or reproducibility claim being tested.
2. **Define acceptance criteria before running.** State the observable,
   tolerance or qualitative decision rule, required diagnostics, anticipated
   compute cost, and failure response.
3. **Execute the smallest discriminating test first.** Do not spend large
   amounts of compute confirming behavior that a reduced case could falsify.
4. **Inspect results rather than only process exit status.** Examine histories,
   strict diagnostics, conservation residuals, convergence behavior, plots,
   resource use, and logs.
5. **Classify every unexpected finding.** Decide whether it is a code bug, an
   analysis bug, a physical-model limitation, an operational/HPC issue, an
   unsupported configuration, or unresolved evidence.
6. **Revise this plan before broadening execution when assumptions change.**
   Add or alter tasks, gates, tolerances, run matrices, or limitations in a
   reviewable documentation commit alongside or before the implementation
   fix.
7. **Add a regression or audit check.** A resolved defect must leave behind a
   test, analysis validation, manifest check, or operational guard that would
   expose recurrence.
8. **Rerun affected gates.** Do not proceed to more expensive tiers until
   corrected lower-tier evidence passes.
9. **Update evidence and budget records.** Archive run products and record
   consumed node-hours, including failed or inconclusive jobs.
10. **Make the next decision explicit.** Proceed, revise, block, or request
   user input; do not silently weaken the objective.

### Stop-the-Line Findings

The following findings immediately block progression to higher-cost or
production-claim tiers:

- a mismatch between implemented equations or threshold policy and the
  documented claimed model;
- a failed independent analytic oracle;
- unexplained pressure floors, nonfinite state, nonpositive pressure, or
  emergency-bound failure in a case expected to be admissible;
- failure to converge with increasing resolution or decreasing timestep where
  convergence is expected;
- unexplained CPU versus Frontier GPU disagreement beyond declared tolerances;
- restart or MPI-decomposition sensitivity beyond declared tolerances;
- analysis products that fail synthetic tests or violate conservation/transfer
  identities;
- forcing that cannot be reconstructed from archived parameters and seed;
- storage, runtime, or I/O behavior that makes the planned campaign unsafe or
  non-reproducible;
- projected or consumed Frontier test usage exceeding the compute budget.

### Required Findings Register

Maintain a tracked findings register in this document or in a linked
follow-up Markdown file once implementation begins. It must use a stable
table such as:

| ID | Date | Phase/run | Finding | Evidence path | Severity | Plan change | Regression added | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| F-001 | TBD | Phase A | Collision updates currently cover two full timesteps per LF cycle | Source audit and corrected analytic test | Blocker | Correct collision split before paper work | Required in Commit 1 | Open |
| F-002 | TBD | Phase B | Current ordinary firehose threshold differs from MKS24 production convention | Source and paper audit | Blocker | Add threshold policy and paper deck setting | Required in Commit 2 | Open |
| F-003 | TBD | Phase B | Hard-bound diagnostic currently depends on backup correction enablement | Source audit and negative test | High | Separate safety from correction policy | Required in Commit 3 | Open |

Each plan-revision commit must update this register or explain why no new
finding was introduced. Each implementation commit resolving a finding must
refer to its finding ID in its commit or PR description.

### Tolerance Revision Rule

Never loosen a numerical tolerance merely to make a failing run pass. A
tolerance change is acceptable only when:

- the old tolerance is shown to be invalid for a stated reason, such as
  reduction-order behavior, output precision, statistical sampling, or a
  mathematically justified discretization error;
- measured evidence from multiple runs supports the new value;
- the change preserves the physical purpose of the test;
- the rationale is recorded in the findings register and relevant docs.

### Readiness Review Cadence

Conduct a formal readiness review:

- after completion of each core physics phase;
- before any Frontier GPU test campaign;
- before any reduced nonlinear science campaign;
- before any standard MKS24 campaign;
- after any stop-the-line finding;
- before stating that the code is production ready.

The readiness review must state which gate-table rows are complete, what
evidence supports them, how much Frontier test budget has been consumed, and
what unresolved issues remain.

## Non-Negotiable Constraints

### Preserve the Current CGL Storage Architecture

Do not redesign the hyperbolic CGL variable layout or Riemann interfaces as part
of this reproduction program.

The required invariant is:

- Outside a dedicated LF split sweep, `IAN` stores the conserved anisotropy
  quantity used by normal CGL evolution, output, restart, FOFC, collision
  handling, and timestep evaluation.
- Inside a dedicated LF split sweep only, `IAN` stores the magnetic-moment
  representation required for the LF parabolic update.
- The begin/end sweep conversions and existing representation assertions are
  the correct architectural boundary.

Relevant current files:

- `src/eos/ideal_c2p_mhd.hpp`
- `src/mhd/mhd.hpp`
- `src/mhd/mhd.cpp`
- `src/mhd/mhd_sts.cpp`
- `src/diffusion/cgl_landau_fluid.cpp`
- `src/outputs/restart.cpp`

### Preserve the Existing LF Split Ownership

LF transport must remain a dedicated split parabolic operator owned through the
current `pcgl_lf` path. Do not fold LF heat flux into ordinary hyperbolic fluxes
or reintroduce older conduction ownership assumptions from donor branches.

### Preserve Existing Composition Restrictions Initially

The present prohibitions on mixing CGL LF with other MHD parabolic processes
remain appropriate until the corrected LF implementation and MKS24 baseline are
stable. Resistivity is especially high risk because the LF magnetic-moment
representation depends directly on the evolving field magnitude and direction.

### Preserve Existing Useful Public Inputs

Keep the current public keys under `<mhd>` unless a narrowly scoped extension is
required:

```text
cgl_heat_flux
cgl_heat_flux_integrator
lf_k_parallel
lf_coefficient_mode
lf_c_parallel0
nu_coll
mirror_limiter
firehose_limiter
limiter_nu_coll
backup_limiters
cgl_lf_strict_admissibility
```

Any new threshold-selection or paper-workflow input must be documented,
recorded in manifests, covered by negative tests where applicable, and set
explicitly in paper reproduction decks.

### Do Not Merge `origin/CGL-STS-LF` Wholesale

The donor branch contains useful paper-validation assets:

```text
docs/cgl_lf_paper_reproduction_plan.md
docs/cgl_lf_paper_validation_runbook.md
inputs/cgl_lf_paper/*.athinput
scripts/analyze_cgl_lf_paper.py
scripts/run_cgl_lf_paper_smoke.sh
src/pgen/cgl_lf_paper.cpp
```

It also has incompatible ancestry, older ownership assumptions, duplicated
physics formulae, and the same incorrect ordinary firehose convention for the
MKS24 target. Port assets selectively and adapt them to the current
architecture; do not merge the branch wholesale.

## Current Baseline Inventory

### Implemented Feature and Safety Infrastructure

The current branch already contains:

| Area | Current asset |
| --- | --- |
| LF operator | `src/diffusion/cgl_landau_fluid.cpp`, `.hpp` |
| Shared CGL constants/helpers | `src/eos/cgl_physics.hpp` |
| CGL EOS and collisions | `src/eos/cgl_mhd.cpp`, `src/eos/ideal_c2p_mhd.hpp` |
| Protected LF representation lifecycle | `src/mhd/mhd_sts.cpp`, `src/mhd/mhd.cpp` |
| Built-in quantitative pgen | `src/pgen/tests/cgl_landau_fluid.cpp` |
| Built-in FOFC pgen | `src/pgen/tests/cgl_fofc.cpp` |
| AMR diagnostic pgen support | `src/pgen/tests/divb_amr.cpp` |
| CPU tests | `tst/test_suite/cgl/test_cgl_landau_fluid_cpu.py` |
| MPI CPU tests | `tst/test_suite/cgl/test_cgl_landau_fluid_mpicpu.py` |
| Workflow driver | `scripts/cgl_lf_workflow.py` |
| Workflow wrapper | `scripts/run_cgl_lf_validation.sh` |
| Existing plots | `scripts/plot_cgl_lf_validation.py` |
| Feature docs | `docs/source/modules/cgl_landau_fluid.md` |
| Validation docs | `docs/source/modules/cgl_landau_fluid_validation.md` |

### Existing Workflow Coverage

The current `scripts/cgl_lf_workflow.py` supports:

| Workflow | Meaning today |
| --- | --- |
| `quick` | One strict one-dimensional LF quantitative run |
| `compare` | STS versus explicit comparison, collisionless and finite-collision |
| `amr` | Strict two-dimensional CGL-LF AMR smoke |
| `full` | Eighteen operator/linear-wave/eigenmode cases plus plots |
| `plot` | Regenerate existing operator-validation plots |
| `summarize` | Regenerate a result summary from retained outputs |

This workflow layer must remain usable throughout the reproduction work.
It should not silently change meaning. In particular, `full` is presently a
closure/operator validation tier, not an MKS24 turbulence reproduction tier.

### Existing Routine Tests

The current CPU suite covers:

- quantitative LF decay and clean diagnostics;
- limiter occupancy while remaining admissible;
- explicit versus capped-STS comparisons;
- finite-collision STS versus explicit consistency;
- CGL FOFC live-flux behavior;
- CGL-LF AMR with conserved prolongation;
- restart reproducibility;
- invalid integrator rejection;
- explicit-mode invalid-composition rejection;
- primitive-AMR-prolongation rejection.

The MPI CPU suite covers one-rank versus four-rank AMR reproducibility.

These tests should remain passing, but several must be revised once the
collision-splitting bug and firehose-threshold policy are corrected.

## Confirmed Defects and Gaps

### Defect 1: Collision Relaxation Advances Two Full Timesteps per LF Cycle

#### Evidence

The driver executes a pre and a post parabolic sweep around the normal
time-integrator work:

```text
src/driver/driver.cpp
  pre STS/parabolic sweep
  main time integrator
  post STS/parabolic sweep
```

The CGL-LF task path invokes collisions after each completed split LF sweep:

```text
src/mhd/mhd_sts.cpp
  STSPostSweepCGLCollisions(...)
    peos->Collisions(...)
```

The collision function always uses:

```cpp
auto &dtc = pmy_pack->pmesh->dt;
SingleColl_CGLMHD(..., dtc, ...);
```

The quantitative pgen oracle explicitly encodes:

```cpp
s = apply_heat_flux(s, 0.5*pm->time);
s = apply_collision(s, pm->time);
s = apply_heat_flux(s, 0.5*pm->time);
s = apply_collision(s, pm->time);
```

#### Physical Consequence

For pure background collision relaxation with decay proportional to
`exp(-nu_coll dt)`, the LF-enabled path applies approximately
`exp(-2 nu_coll dt)` per full cycle. Limiter-driven collision effects are also
over-applied whenever their activation persists across the two sweeps.

#### Why Existing Tests Do Not Protect Against It

Current finite-collision tests compare STS and explicit LF paths that share the
same collision scheduling. They establish integrator consistency under the
current schedule, not physical correctness of that schedule.

#### Required Resolution

Future agents must first determine and document the intended operator split.
Two defensible forms are:

1. Apply collision relaxation once per full physical timestep, preserving the
   normal non-LF collision-update meaning.
2. Retain symmetric placement around the split update but use two
   half-collision updates, each with `dt/2`.

The recommended implementation is the second form only if the task graph
continues to treat collision relaxation as part of the symmetric LF/source
split. Otherwise, remove the LF-specific duplicate call and use the normal
once-per-cycle collision application.

Do not choose by convenience. Derive the operator ordering, document it in the
code and Sphinx docs, and verify it analytically.

### Defect 2: Current Firehose Activation Does Not Match MKS24

#### Current Behavior

`src/eos/cgl_physics.hpp` currently defines:

```cpp
constexpr Real kFirehoseThreshold = -0.7;
constexpr Real kFirehoseHardBound = -1.0;
constexpr Real kMirrorThreshold = 0.5;
constexpr Real kMirrorHardBound = 1.0;
```

In AthenaK normalized units, with `beta = 2 p/B^2` and
`Delta = (p_perp - p_parallel)/p`, these correspond to:

| Code threshold | Physical convention |
| --- | --- |
| `p_perp - p_parallel <= -0.7 B^2` | `beta Delta <= -1.4`, oblique firehose |
| `p_perp - p_parallel <= -1.0 B^2` | `beta Delta <= -2`, parallel firehose |
| `p_perp - p_parallel >= 0.5 B^2` | `beta Delta >= 1`, mirror |

#### MKS24 Requirement

MKS24 states that its production CGL simulations activate collisional
microinstability regulation at:

```text
mirror:   beta Delta > 1
firehose: beta Delta < -2
```

The paper later discusses `beta Delta = -1.4` as an alternative oblique
threshold and notes that it changes unstable-volume fractions. Therefore the
current hard-coded ordinary activation convention cannot reproduce the
published production simulations.

#### Required Resolution

Introduce an explicit threshold policy. Preferred user-facing design:

```text
<mhd>
cgl_firehose_threshold = parallel   # MKS24 production convention: beta Delta = -2
```

with an alternative:

```text
cgl_firehose_threshold = oblique    # beta Delta = -1.4 comparison convention
```

Alternative designs, such as a numeric `firehose_beta_delta_threshold`, are
acceptable only if they are less error-prone in input decks and clearly
documented. In every design:

- paper input decks must specify the threshold explicitly;
- manifests and summaries must record the threshold choice;
- diagnostics must distinguish threshold occupancy from safety/floor failures;
- existing limiter tests must be parameterized for both conventions;
- the default must be intentionally selected and documented.

Because the current feature branch may already have users, do not silently
change defaults without deciding whether backward compatibility or paper
correctness governs the public default. Paper workflows must not depend on the
default either way.

### Defect 3: Hard-Bound Monitoring Depends on Corrective Backup Mode

#### Current Behavior

`CGLLandauFluid::RecordAdmissibility` increments `lf_hardbd` only when
`eos.backup_lim` is enabled. Thus the diagnostic answers:

```text
did the state exceed a hard bound while backup correction was enabled?
```

rather than:

```text
did the state exceed a hard bound?
```

#### Consequence

With strict admissibility enabled but backup limiters disabled, a hard-bound
crossing can be invisible to `lf_hardbd` and therefore invisible to strict
termination.

#### Required Resolution

Separate three concepts:

1. **Physical threshold occupancy:** expected in limiter-regulated paper runs
   and required for reproducing unstable-volume fractions.
2. **Limiter action or effective scattering:** how much regulation is applied.
3. **Safety violation or unacceptable overshoot:** a strict validation failure,
   independent of whether automatic backup correction was enabled.

`lf_hardbd` must be a truthful safety diagnostic or must be renamed and
redefined. The preferred narrow fix is to count safety-bound crossings
regardless of `backup_lim`, while using separate paper diagnostics for physical
threshold occupancy.

### Gap 4: Current Workflows Are Operator Validation, Not Paper Reproduction

The existing `full` matrix covers linear and controlled validation cases:

- parallel/perpendicular heat-flux decay;
- finite collision decay;
- `grad |B|` response;
- heat-flux limiting and limiter suppression;
- mirror/firehose limiter exercise;
- field waves;
- oblique waves;
- pure-CGL versus LF eigenmode comparisons.

MKS24 is instead a forced-turbulence study requiring:

- fully periodic driven three-dimensional turbulence;
- active and passive pressure-anisotropy dynamics;
- Alfvenic and random forcing;
- beta scans;
- limiter-scattering-frequency scans;
- nonlinear steady-state statistics;
- spectra, transfer functions, PDFs, and alignment diagnostics.

No paper-grade result should be claimed from the existing `full` workflow.

### Gap 5: Current Diagnostics Cannot Produce the Published Figures

Existing `lf_*` history counters are useful safety counters, but they do not
provide:

- physical unstable volume fractions at the paper threshold;
- active versus passive comparisons;
- `C_B2 = <B^4>/<B^2>^2 - 1`;
- density and pressure-anisotropy PDFs;
- field-parallel and perpendicular gradient spectra;
- anisotropic pressure-stress transfer functions;
- local-field/rate-of-strain alignment distributions;
- heat-flux power and cap-activity breakdowns;
- forcing injection and closed energy budgets.

These products must be implemented before figure-level reproduction is
possible.

## MKS24 Reproduction Contract

### What Paper Is Being Reproduced

The target source bundle is:

```text
docs/arXiv-2405.02418v2/MKS24.tex
docs/arXiv-2405.02418v2/fig*.pdf
```

The paper is:

```text
Self-organization in collisionless, high-beta turbulence
S. Majeski, M. W. Kunz, and J. Squire
```

Future agents should read the source directly before modifying reproduction
inputs or interpreting plots. At minimum, reread:

```text
MKS24.tex:204-239   theoretical CGL-LF equations and heat-flux closure
MKS24.tex:463-469   simulation setup, limiter convention, forcing, run length
MKS24.tex:478-496   spectra, transfer, and alignment diagnostics
MKS24.tex:503-515   unstable-volume and pressure-density results
MKS24.tex:662-698   limiter-frequency scan and threshold discussion
```

### Governing CGL-LF Closure to Match

The paper uses the CGL pressure equations closed by the 3+1 Landau-fluid heat
fluxes:

```text
q_perp = - v_th_parallel/(sqrt(pi) |k_parallel|)
         [ rho grad_parallel(p_perp/rho)
           - p_perp (1 - p_perp/p_parallel) grad_parallel(B)/B ]

q_parallel = - 2 v_th_parallel/(sqrt(pi) |k_parallel|)
             rho grad_parallel(p_parallel/rho)

v_th_parallel = sqrt(2 p_parallel/rho)
```

The current collisionless code coefficients are consistent with these formulae
under the implemented `c_parallel = sqrt(p_parallel/rho)` normalization. That
mapping must be retained and documented in any paper runbook.

### Required Numerical Setup

The paper production configuration includes:

| Quantity | Paper requirement |
| --- | --- |
| Domain | Fully periodic `[Lx, Ly, Lz] = [1, 1, 2]` |
| Mean magnetic field | Along `z`, `B0 = e_z` in paper units |
| Standard resolution | `n_perp = 192`, `n_parallel = 384` |
| Beta values | `beta0 = 1, 10, 100` |
| LF closure scale | `|k_parallel| = 4 pi/L_parallel` for standard runs |
| Production limiter rate | `nu_lim = 1e10 v_A/L_perp` |
| Run duration | At least `10 L_perp/v_A` |
| Default threshold model | Mirror `beta Delta > 1`; firehose `beta Delta < -2` |

AthenaK uses magnetic fields with `v_A = B/sqrt(rho)`. For initially isotropic
pressure in code units:

```text
beta = 2 p0/B0^2
p0 = 0.5 beta B0^2
```

Paper inputs must state this normalization explicitly.

### Required Forcing

MKS24 forces velocity perturbations through an Ornstein-Uhlenbeck correlated
process. Paper-standard runs require:

| Quantity | Value or mode |
| --- | --- |
| Fiducial energy injection rate | `d_t E_K = 0.32 rho0 v_A^2 L_perp^3` |
| Correlation time | `t_corr = L_parallel/v_A` |
| Forced wavenumber range | `k in (2 pi/L_parallel) * [1, 3]` |
| Power distribution | proportional to `k^-2` |
| Alfvenic forcing | Perpendicular to `e_z`, incompressible at outer scale |
| Random forcing | Velocity forced without Alfvenic directional constraint |

The reproduction implementation must either use existing AthenaK turbulence
forcing infrastructure after verifying it meets these requirements, or extend
it narrowly with deterministic controls. Do not create an unverified forcing
surrogate and label it a paper reproduction.

### Active and Passive Pressure-Anisotropy Modes

The paper compares:

- **active-Delta:** CGL anisotropic pressure feeds back on the turbulent
  momentum evolution;
- **passive-Delta:** MHD-like dynamics evolve while pressure anisotropy is
  advanced diagnostically without anisotropic-pressure feedback on momentum.

The current CGL solver has a `passive` EOS capability in its Riemann path. The
paper pgen and inputs must verify that this mode reproduces the intended
active/passive distinction in the current architecture and that LF evolution
and collision policy remain consistent in both modes.

### Limits of Reproduction Claims

MKS24 uses a CGL-LF fluid surrogate for microinstability scattering; that part
is reproducible with AthenaK once corrected. Other comparisons in the broader
literature use hybrid-kinetic physics or finite electron pressure. Therefore:

- do claim reproduction of the CGL-LF turbulence calculations once the setup,
  closure, threshold convention, forcing, diagnostics, and statistical outputs
  are matched;
- do not claim reproduction of kinetic-scale mirror/firehose fluctuations;
- label cold-electron analogues clearly if comparing against configurations
  requiring nonzero isothermal electron pressure;
- treat later evolving-`nu_lim` or kinetic comparisons as extensions, not
  prerequisites for reproducing the published CGL-LF simulations.

## Paper Figure and Observable Map

The implementation is complete only when it can produce the following
categories of results from archived, reproducible runs. Exact visual agreement
is not the first acceptance gate; first require matched setup, correct
diagnostics, converged statistics, and documented quantitative comparisons.

| Paper result | Required simulation modes | Required products | Current status |
| --- | --- | --- | --- |
| Figure 2: pressure-density distributions and unstable-volume history | Active and passive, Alfvenic, representative beta values | PDFs of `delta p_parallel`, `delta p_perp`, density; volume fraction beyond thresholds versus time | Missing |
| Figures 3-4: compressive/density behavior | Active/passive; beta and forcing variants | Density PDFs and density/compressive spectra | Missing |
| Figures 5-6: turbulent spectra and rate-of-strain suppression | Active/passive, Alfvenic/random | Velocity/magnetic/pressure spectra; parallel/perpendicular velocity-gradient spectra | Missing |
| Figure 7: pressure-anisotropy organization and transfer | Active/passive, Alfvenic/random, beta 10 | `grad_parallel Delta p`, `grad_perp Delta p` spectra; `T_Delta p(k_perp)/T_total` | Missing |
| Figures 8-9: alignment diagnostic | Active/passive; beta/forcing variants | Scale-dependent PDFs and peak locations of local-field/rate-of-strain eigenvector alignment | Missing |
| Figure 10: kinetic comparison context | Not reproduced directly by CGL-LF | Documentation of non-reproducible kinetic comparator and any qualitative CGL-only comparison | Out of direct scope |
| Figure 11: scale-separation behavior | Selected resolutions | Resolution/convergence or scale-separation comparison products | Missing |
| Figure 12: heat-flux sensitivity | Heat-flux parameter variants | Spectra/transfer/diagnostic products under heat-flux changes | Missing |
| Figure 13: limiter-frequency sensitivity | Beta 100, Alfvenic, several `nu_lim` | Energy spectra, `beta Delta` PDFs, strain PDFs, pressure-stress transfer | Missing and threshold-blocked |

## Delivery Strategy

### Guiding Principle

Fix physical correctness before expanding scientific scope. A larger test
matrix built on a doubled collision rate or a mismatched instability convention
would produce reproducible but physically misleading results.

### Recommended Pull Request Stack

Use small reviewable pull requests or a strictly ordered commit stack. Do not
combine all paper infrastructure and all physics fixes into a single
unreviewable change.

| Stage | Purpose | Must land before |
| --- | --- | --- |
| A | Collision-split correction and analytic regression | Any finite-collision paper claim |
| B | Configurable threshold convention and truthful hard-bound diagnostics | Any unstable-volume or limiter paper claim |
| C | Accuracy/convergence and diagnostic semantics reinforcement | Paper pgen review |
| D | Paper pgen, forcing setup, and reduced smoke workflow | Standard turbulence runs |
| E | Analysis/figure pipeline and archived result contract | Paper-grade execution |
| F | Standard and HPC reproduction runs plus interpretation docs | Publication/review claims |

## Phase A: Correct Collision and Source-Term Time Integration

### Objective

Make every collision and limiter-scattering update correspond to the intended
physical interval exactly once per full evolution cycle.

### Required Investigation

Before editing, trace all collision invocations in configurations:

1. CGL without LF.
2. CGL with LF STS.
3. CGL with LF explicit reference.
4. CGL with limiter scattering but no background collisions.
5. CGL with both background and limiter scattering.

Search locations:

```bash
rg -n "Collisions|STSPostSweepCGLCollisions|nu_coll|lim_coll|SingleColl" src
rg -n "parabolic|after_timeintegrator|before_timeintegrator|sts" src/mhd src/driver
```

Produce a short operator-ordering note in the PR description or an associated
developer document. It should show the intended evolution composition, for
example:

```text
H(dt)                     hyperbolic CGL operator
L(dt)                     LF heat-flux operator
C(dt)                     background and limiter collision operator

Candidate symmetric form:
C(dt/2) L(dt/2) H(dt) L(dt/2) C(dt/2)
```

Do not assume this form until checking how non-LF collisions are intended to
compose with the existing driver.

### Preferred Implementation Scope

Likely touched files:

```text
src/eos/cgl_mhd.cpp
src/eos/cgl_mhd.hpp or relevant EOS declaration header
src/mhd/mhd_sts.cpp
src/mhd/mhd_tasks.cpp
src/driver/driver.cpp, only if a timestep fraction must be provided cleanly
src/pgen/tests/cgl_landau_fluid.cpp
tst/test_suite/cgl/test_cgl_landau_fluid_cpu.py
inputs/tests/cgl_lf_decay.athinput or a new focused collisional input
docs/source/modules/cgl_landau_fluid.md
docs/source/modules/cgl_landau_fluid_validation.md
```

Pass an explicit collision substep duration into the collision update if the
function is used in multiple split contexts. Do not hide half-step behavior by
reading a mutable global or by changing `pmesh->dt`.

### Required Tests

Add or revise tests that do not compare two implementations with the same
incorrect operator scheduling:

1. **Pure background collisional relaxation with negligible LF transport.**
   Initialize anisotropy with spatially uniform thermodynamics and nonzero
   `nu_coll`; measure analytic decay over one or multiple timesteps:

   ```text
   Delta p(t) = Delta p(0) exp(-nu_coll t)
   ```

   Use this as the primary guard against doubled collision advancement.

2. **Collision plus weak LF analytic/semi-analytic case.**
   Retain the existing decay test only after changing its oracle to the
   physically justified operator order and timestep fractions.

3. **Explicit-versus-STS finite-collision comparison.**
   Keep this as an integration-consistency check, not the primary correctness
   check.

4. **Limiter-triggered collision duration check.**
   Construct a controlled threshold-active case and verify relaxation depends
   on the intended duration, without unexpected floors or hard-bound events.

5. **Restart regression with finite collision enabled.**
   Once the correction is stable, ensure restart continuation matches an
   uninterrupted path under the corrected split.

### Acceptance Criteria

- Analytic background-collision decay matches the intended `exp(-nu t)` law,
  not `exp(-2 nu t)`, at a tolerance justified by the chosen timestep and
  integrator.
- STS and explicit paths agree under the corrected source-term scheduling.
- Current collisionless tests remain passing.
- Strict diagnostics remain clean for well-posed test inputs.
- The operator-ordering choice is documented and represented in test names or
  comments; no regression test encodes unexplained driver behavior.

### Recommended Commit

```text
Correct CGL-LF collision source-term splitting and add analytic regression
```

## Phase B: Implement Paper-Compatible Threshold Policies and Diagnostics

### Objective

Make the instability-limiter convention an explicit, auditable physical model
choice and make safety counters truthful independently of corrective actions.

### Public Interface Design

Preferred public key:

```text
<mhd>
cgl_firehose_threshold = parallel
```

Accepted values:

| Value | Activation condition | Use |
| --- | --- | --- |
| `parallel` | `beta Delta <= -2` | MKS24 production reproduction |
| `oblique` | `beta Delta <= -1.4` | Alternative comparison experiment |

The mirror threshold remains:

```text
beta Delta >= 1
```

If implementation experience favors a numeric parameter, use a clear name,
validate sign/range, and still provide paper inputs that state the value
explicitly. Do not make paper behavior depend on an undocumented default.

### Threshold Versus Safety Semantics

The present notion of ordinary threshold and backup hard bound requires careful
redesign when the paper activation threshold is itself `beta Delta = -2`.
Future agents must decide and document:

- which threshold triggers physical anomalous scattering;
- whether overshoot beyond that threshold is allowed during finite-rate
  regulation;
- what state constitutes an unacceptable numerical safety violation;
- whether `backup_limiters` is a correction mechanism, a different physical
  limiter model, or both.

Recommended diagnostic separation:

| Diagnostic category | Meaning | Strict failure? |
| --- | --- | --- |
| `mirror_active` / `firehose_active` | State is physically over the chosen instability activation threshold | No; expected in paper runs |
| `nu_eff` or limiter-action aggregate | Scattering applied due to limiter policy | No |
| `hard_bound` or `overshoot_failure` | State crosses a defined emergency numerical bound | Yes in strict validation |
| Floor/nonfinite/nonpositive counters | Numerical state invalid or required repair | Yes in strict validation |

Do not equate paper unstable occupancy with numerical failure.

### Physics Helper Refactor

Keep threshold computations centralized in `src/eos/cgl_physics.hpp`. It should
provide, as needed:

- a small threshold-policy enum or structured parameters;
- conversion between input policy and code-unit coefficients;
- mirror/firehose activation predicates;
- emergency-bound predicates;
- effective collision-rate calculation;
- heat-flux cap constants and helper;
- diagnostics-friendly threshold values.

Use the helper in:

```text
src/eos/cgl_mhd.cpp
src/diffusion/cgl_landau_fluid.cpp
src/pgen/tests/cgl_landau_fluid.cpp
future paper pgen and history diagnostics
```

Never duplicate firehose constants in the paper pgen or analysis layer.

### Required Tests

1. Input parsing tests for supported and unsupported threshold policies.
2. A controlled state scan crossing mirror and both firehose conventions.
3. A paper-policy limiter case that activates at `beta Delta = -2`.
4. An oblique-policy case that activates between `-1.4` and `-2`.
5. A strict diagnostic case showing emergency-bound reporting occurs even
   with backup correction disabled.
6. Manifest/summary checks proving the selected policy is archived.

### Documentation Changes

Update:

```text
docs/source/modules/cgl_landau_fluid.md
docs/source/modules/cgl_landau_fluid_validation.md
docs/source/reference/input_parameters.md
```

Documentation must explicitly state:

- which threshold policy is selected by each reproduction input;
- why MKS24 decks use the parallel-fluid convention;
- why the oblique convention remains available;
- that limiter occupancy is physical activity, while strict safety failures
  are separate.

### Acceptance Criteria

- MKS24 input decks explicitly activate the parallel-firehose policy.
- Alternative oblique-firehose experiments are possible and identifiable.
- Safety counters do not depend on whether automatic correction is enabled.
- No diagnostic called a safety failure is actually a normal paper observable.

### Recommended Commits

```text
Add selectable CGL firehose threshold policies for paper reproduction
Make CGL-LF hard-bound diagnostics independent of backup correction
```

## Phase C: Upgrade Numerical Accuracy and Physics Validation

### Objective

Establish a credible numerical-accuracy basis before running costly nonlinear
turbulence simulations.

### Existing Tests to Preserve

Preserve the role of:

- LF decay;
- heat-flux cap;
- limiter activity;
- oblique linear-wave and eigenmode checks;
- FOFC;
- explicit reference comparisons;
- AMR smoke;
- restart;
- MPI AMR reproducibility.

Revise finite-collision tests after Phase A and threshold-sensitive tests after
Phase B.

### New Required Accuracy Tests

#### C1. Collisionless Damping Convergence

For parallel and perpendicular LF decay:

- use at least three resolutions;
- use at least three timestep caps or explicit small-step references;
- fit the error trend rather than checking one final tolerance;
- report both absolute and relative errors.

#### C2. Finite-Collision Asymptotics

Verify the LF coefficient limits for:

```text
nu_eff << c_parallel |k_parallel|
nu_eff ~  c_parallel |k_parallel|
nu_eff >> c_parallel |k_parallel|
```

Check both `q_parallel` and `q_perp`, and both local and background coefficient
modes.

#### C3. Rotated-Field Invariance

Run the same field-aligned thermal perturbation with the mean field directed:

```text
x
y
z
obliquely
```

Project onto the local/mean field direction and verify equivalent damping to
the measured discretization order.

#### C4. Heat-Flux Cap Extremes

Construct large gradients such that:

```text
|q_unlimited|/q_max > 10^4
```

Require:

- bounded limited flux;
- correct flux sign;
- finite positive pressures;
- no strict-mode failure in admissible configured cases;
- no STS overshoot compared with explicit reference.

#### C5. Low-Field Robustness

Exercise regions approaching the configured magnetic-field floor while LF is
active. Verify:

- no nonfinite field-direction operation;
- controlled cap and coefficient behavior;
- documented treatment of the low-field limit;
- finite pressure and energy states.

#### C6. Conservation and Energy Accounting

For closed non-forced cases, record energy residuals attributable to
discretization and split transport. For forced cases later, close the budget
using injected work, anisotropic-pressure work, heat-flux transport, and
limiter/collision effects where represented.

#### C7. Timestep/STS Sensitivity

For representative smooth and limiter-active problems, scan:

```text
cgl_heat_flux_integrator = explicit
cgl_heat_flux_integrator = sts
sts_max_dt_ratio
resolution
```

Do not state that `sts_max_dt_ratio = 1` makes STS identical to explicit. It
sets a parabolic timestep scale; it does not change method identity.

### Test Tier Placement

| Test group | Routine CPU CI | Manual local validation | Nightly/HPC |
| --- | --- | --- | --- |
| Corrected analytic collision decay | Yes | Yes | Yes |
| Threshold input/strict diagnostic negatives | Yes | Yes | Yes |
| One collisionless damping case | Yes | Yes | Yes |
| One explicit-versus-STS case | Yes | Yes | Yes |
| Rotated-field suite | No | Yes | Yes |
| Resolution/timestep convergence suite | No | Yes | Yes |
| Extreme cap and low-field stress suite | Small case optional | Yes | Yes |
| Nonlinear turbulence | No | Reduced smoke only | Yes |

### Acceptance Criteria

- All physical-model conventions are configurable and archived.
- The corrected finite-collision behavior has an analytic independent oracle.
- Accuracy plots or tables exist for convergence-sensitive tests.
- No paper pgen work begins from a branch with failing corrected core tests.

### Recommended Commit

```text
Extend CGL-LF accuracy validation after source-term and limiter corrections
```

## Phase D: Implement the MKS24 Paper Problem Generator

### Objective

Add a current-architecture problem generator capable of initializing the
MKS24 turbulence runs and associated controlled comparison modes.

### Porting Strategy

Use `origin/CGL-STS-LF:src/pgen/cgl_lf_paper.cpp` as an intent/reference
source only. Its valuable components include:

- paper-oriented mode enumeration;
- initialization ideas for turbulence and selected wave analogues;
- active/passive intent;
- user-history concepts such as `C_B2`, limiter occupancy, effective
  collisions, and heat-flux cap fractions.

Do not copy uncritically because it:

- assumes older LF ownership;
- duplicates threshold and heat-flux physics;
- uses the wrong ordinary firehose convention for MKS24 production;
- does not implement the complete paper figure-analysis pipeline.

### Recommended Location and Registration

Preferred implementation:

```text
src/pgen/tests/cgl_lf_paper.cpp
```

registered in the default test/problem binary analogously to:

```text
src/pgen/tests/cgl_landau_fluid.cpp
src/pgen/tests/divb_amr.cpp
```

with:

```text
<problem>
pgen_name = cgl_lf_paper
```

If production-scale forcing infrastructure cannot cleanly operate through the
built-in test binary, a dedicated build-time pgen may be justified, but this
must be documented and the workflow driver must hide the operational
complexity from users.

### Required Pgen Modes

The pgen should support at least:

| Mode | Purpose | Required for paper baseline? |
| --- | --- | --- |
| `turbulence` | MKS24 driven periodic turbulence | Yes |
| `turbulence_passive` or `passive=true` with `turbulence` | Passive-Delta control | Yes |
| `np_mode` | Controlled non-propagating/pressure-balance analogue | Useful validation, not primary MKS24 figure requirement |
| `fast_wave` | Compressive controlled analogue | Useful validation |
| `oblique_iaw` | Oblique compressive controlled analogue | Useful validation |
| `linear_wave_scan` | Wave/frequency diagnostics | Useful validation |

Avoid adding modes not attached to a defined test or paper diagnostic.

### Initial Conditions for Turbulence

Required paper-standard initialization:

```text
rho0 = 1
B0 = (0, 0, 1)
p_parallel0 = p_perp0 = 0.5 * beta0 * B0^2
domain = [1, 1, 2]
boundary conditions = periodic in all directions
```

Pgen initialization must correctly construct:

- primitive density and velocity;
- face-centered magnetic field;
- cell-centered magnetic field as expected by initialization;
- `p_parallel` and `p_perp`;
- total conserved energy;
- conserved CGL anisotropy representation in `IAN`;
- user-history registration;
- restart-safe behavior.

Never initialize `IAN` in the temporary LF magnetic-moment representation.

### Forcing Integration

Before writing new forcing code, inspect the existing turbulence pgen and
forcing infrastructure:

```text
src/pgen/turb.cpp
related turbulence driver/source modules and input examples
```

The paper pgen must support:

- deterministic random seed recording;
- Alfvenic forcing;
- random forcing;
- the published correlation time and injection rate;
- the published forced-mode range and power distribution;
- restart continuation without changing the statistical forcing process unless
  this limitation is explicitly documented.

If existing forcing cannot provide a reproducible restart stream or cannot
match the published modes, repair or extend it as a separate, narrowly scoped
subtask with its own tests.

### Required Input Decks

Add versioned paper input directories rather than mixing production runs with
unit inputs:

```text
inputs/cgl_lf_paper/
  cgl_lf_paper_smoke_active_beta10.athinput
  cgl_lf_paper_smoke_passive_beta10.athinput
  cgl_lf_paper_standard_active_alfvenic_beta1.athinput
  cgl_lf_paper_standard_active_alfvenic_beta10.athinput
  cgl_lf_paper_standard_active_alfvenic_beta100.athinput
  cgl_lf_paper_standard_passive_alfvenic_beta10.athinput
  cgl_lf_paper_standard_active_random_beta10.athinput
  cgl_lf_paper_nulim_beta100_20.athinput
  cgl_lf_paper_nulim_beta100_200.athinput
  cgl_lf_paper_nulim_beta100_hardwall.athinput
```

Names may be refined, but each deck must be self-identifying and explicitly
state:

- threshold policy;
- LF integrator;
- `lf_k_parallel`;
- coefficient mode;
- `nu_coll` and/or limiter collision rate;
- forcing mode, seed, correlation time, and injection rate;
- output cadence;
- analysis/statistical averaging window;
- whether it is smoke, convergence, or paper-grade.

### Pgen-Level Acceptance Criteria

- A small active turbulence case initializes and evolves with clean strict
  safety diagnostics.
- A small passive case evolves with the documented passive feedback behavior.
- Restart does not corrupt the CGL representation or forcing setup.
- Threshold and LF closure choices are visible from input and manifest output.
- No pgen-local duplicated threshold or LF formula can drift from shared
  physics helpers.

### Recommended Commit

```text
Add current-API CGL-LF paper problem generator and reduced turbulence inputs
```

## Phase E: Add Paper Diagnostics and Analysis Products

### Objective

Produce the observables needed to make quantitative comparisons with MKS24,
while keeping raw simulation output manageable and archived analyses
reproducible.

### Diagnostic Architecture

Use two layers:

1. **In-simulation reduced histories** for inexpensive, frequent global
   diagnostics and immediate failure checks.
2. **Offline snapshot analysis** for spatial PDFs, spectra, transfer
   functions, and alignment statistics.

Do not attempt to encode every spectral product as a history scalar. Do not
require enormous raw dumps at every small timestep.

### Reduced History Diagnostics

Add paper-pgen user history or a reusable CGL diagnostic facility containing
well-defined labels. At minimum include:

| Quantity | Definition or purpose |
| --- | --- |
| Volume and mass | Normalization and sanity checking |
| Kinetic, magnetic, CGL thermal, and total energies | Energy evolution |
| `<B^2>` and `<B^4>` | Compute `C_B2 = <B^4>/<B^2>^2 - 1` |
| Mean beta | State tracking |
| Mean `Delta p` and mean `|Delta p|` | Pressure-anisotropy amplitude |
| Mirror-threshold volume | Paper threshold occupancy |
| Firehose-threshold volume | Paper threshold occupancy under selected policy |
| Emergency-bound volume | Safety interpretation only |
| Mean or integrated `nu_eff` | Limiter-scattering activity |
| Heat-flux cap occupancy | Fractions with `|q_L|/q_max > 1` and `> 10` |
| Heat-flux power proxy | Diagnostic of LF transport influence |
| Injected forcing work | Required for budget closure |
| Anisotropic-pressure work if available robustly | Required for transfer/budget interpretation |

All threshold-volume diagnostics must evaluate the user-selected threshold
policy, and the manifest must identify that policy.

### Snapshot Output Requirements

At scientifically useful cadence, retain enough fields for offline analysis:

```text
rho
velocity components
magnetic-field components
p_parallel
p_perp
possibly total/internal energy and forcing diagnostics
```

The output design should support:

- local beta and `beta Delta` calculations;
- density and pressure PDFs;
- scalar/vector Fourier spectra;
- local field-direction projected gradients;
- transfer-function analysis;
- rate-of-strain eigenvector alignment analysis.

Use sparse high-volume snapshots plus frequent histories for paper-grade runs.
Record exact cadence in manifests.

### Analysis Package

Add a dedicated standard-library-first analysis interface, allowing NumPy and
Matplotlib for scientific post-processing:

```text
scripts/analyze_cgl_lf_paper.py
scripts/plot_cgl_lf_paper.py            # optional split if analysis grows
scripts/cgl_lf_workflow.py              # orchestration and manifest integration
```

The donor analyzer is a useful beginning for history parsing and summary JSON,
but it does not compute spectra, transfer functions, or alignment diagnostics.

### Required Analysis Capabilities

#### E1. History Summaries

For every paper run:

- parse reduced history;
- report safety-counter cleanliness;
- report threshold-volume statistics over both transient and steady windows;
- report energy/injection summaries;
- report `C_B2`;
- report cap and limiter activity;
- write machine-readable JSON and readable Markdown.

#### E2. PDF Products

Compute:

- pressure-density fluctuation PDFs for Figure 2-style comparisons;
- `beta Delta` PDFs for limiter-frequency comparisons;
- selected strain-distribution PDFs for Figure 13-style interpretation.

Record binning, normalization, selected snapshot times, and averaging windows.

#### E3. Spectra

Compute bin-averaged spectra in the convention described by MKS24:

- velocity and magnetic fluctuations;
- density/compressive quantities;
- `Delta p`;
- field-parallel and perpendicular gradients of `Delta p`;
- field-parallel and perpendicular components of velocity gradients.

Implement and test isotropic and `k_perp` shell binning as required. Specify
normalization conventions in documentation and output metadata.

#### E4. Pressure-Anisotropy Transfer Function

Implement the paper-defined `T_Delta p(k_perp)` diagnostic. This is one of the
highest-risk analysis tasks because implementation choices can alter sign and
normalization. Requirements:

- reproduce the formula from the paper in code comments or analysis docs;
- write synthetic-field tests with known zero or sign behavior;
- compare integral transfer against a real-space anisotropic work measure where
  possible;
- archive normalization conventions in summary metadata.

#### E5. Alignment Diagnostic

Implement the local-field/rate-of-strain alignment diagnostic:

- construct the symmetrized velocity-gradient tensor;
- compute real eigenvectors/eigenvalues at each valid analysis point;
- identify stretching and compressing eigenvectors consistently;
- compute alignment with the local magnetic-field direction;
- bin probability distributions by `k_perp` after the required filtering;
- generate peak-alignment summaries.

This requires careful numerical validation on synthetic flow fields before it
is interpreted physically.

### Figure Reproduction Outputs

Each paper analysis workflow should write:

```text
manifest.json
summary.md
diagnostics.json
logs/
history/
snapshots/ or links-to-snapshots/
analysis/
  pdfs/
  spectra/
  transfer/
  alignment/
figures/
```

Generated large run products belong in ignored output/build areas or external
archival storage, not committed directly to git. Commit analysis code, input
decks, small reference metadata, and curated documentation figures only when
needed.

### Acceptance Criteria

- Reduced histories and analysis scripts operate on archived output without
  rerunning simulations.
- Synthetic tests protect spectra/transfer/alignment implementations.
- Figure-generation scripts state every simulation case and analysis window
  required for each panel.
- Paper diagnostics never reuse LF stage-sampling counters as substitutes for
  physical volume/time statistics without clearly stating that difference.

### Recommended Commits

```text
Add reduced CGL-LF paper diagnostics and archived summary products
Add CGL-LF paper PDF, spectra, transfer, and alignment analysis
```

## Phase F: Extend the Reproducible Workflow Driver

### Objective

Make operator validation and paper reproduction easy to launch, archive,
inspect, and rerun without confusing the two tiers.

### Preserve Existing Commands

Keep existing user commands functional:

```bash
python3 scripts/cgl_lf_workflow.py quick
python3 scripts/cgl_lf_workflow.py compare
python3 scripts/cgl_lf_workflow.py amr
python3 scripts/cgl_lf_workflow.py full
python3 scripts/cgl_lf_workflow.py plot --output-dir <existing-bundle>
python3 scripts/cgl_lf_workflow.py summarize --output-dir <existing-bundle>
scripts/run_cgl_lf_validation.sh
```

Document that `full` means complete **operator validation** unless explicitly
renamed with a compatibility alias.

### New Paper Workflow Modes

Add clear modes whose computational cost and scientific meaning are evident:

| Workflow | Purpose | Expected cost |
| --- | --- | --- |
| `paper-smoke` | Reduced active/passive and limiter-policy startup checks | Workstation/local |
| `paper-convergence` | Selected resolution/timestep/threshold studies | Local cluster or modest HPC |
| `paper-standard` | Published-resolution selected MKS24 matrix | HPC |
| `paper-nulim` | Figure 13 limiter-frequency suite | HPC |
| `paper-analyze` | Regenerate paper diagnostics and figures from archived runs | Analysis host |
| `paper-summary` | Regenerate human/machine summaries without plotting | Any host |

The exact CLI naming may change, but avoid overloading `full` with runs that
are orders of magnitude more expensive.

### Manifest Schema Evolution

The current bundle schema should be extended or version-bumped. Every paper
bundle must record:

```text
schema_version
workflow
status
created_utc
git_revision
executable
build_configuration
build_dir
host and MPI metadata where available
cases
diagnostics
analysis_products
```

Each paper case must record:

```text
input_file
exact_command
runtime_overrides
resolution
domain
beta0
lf_k_parallel
lf_coefficient_mode
cgl_heat_flux_integrator
cgl_firehose_threshold or numeric equivalent
mirror/firehose limiter settings
nu_coll
limiter_nu_coll
forcing_mode
forcing_seed
forcing_parameters
output_cadence
analysis_window
log_path
history_paths
snapshot_paths or archive references
case_status
measured_checks
```

### Result Naming and Storage

Retain the current ignored-bundle principle, with a separate paper namespace if
useful:

```text
build-cgl-implementation/cgl_lf_runs/<timestamp>-<workflow>/
```

or:

```text
build-cgl-implementation/cgl_lf_paper_runs/<timestamp>-<workflow>/
```

Paper-grade raw snapshots may be too large for the worktree disk. Support a
user-supplied output/archive directory and record absolute or portable
references to the result products in the manifest.

### Failure Policy

Paper workflows must fail early on:

- nonzero process return code;
- missing required output;
- nonfinite/nonpositive states;
- floors or emergency-bound failures in strict smoke/verification runs;
- manifest inconsistency;
- failed analysis validation.

Paper production runs may legitimately have nonzero physical
mirror/firehose-threshold occupancy. That activity must not be treated as a
workflow failure.

### Recommended Commit

```text
Add tiered CGL-LF MKS24 execution and analysis workflows
```

## Phase G: Documentation and User Workflow

### Objective

Prevent users and reviewers from confusing kernel correctness, robustness
testing, and paper reproduction.

### Documents to Update

Update normal Sphinx sources:

```text
docs/source/modules/cgl_landau_fluid.md
docs/source/modules/cgl_landau_fluid_validation.md
docs/source/reference/input_parameters.md
docs/source/index.md, only if adding a distinct paper-reproduction page
```

Consider adding:

```text
docs/source/modules/cgl_landau_fluid_paper_reproduction.md
```

if the paper workflow becomes substantial enough that it would overwhelm the
validation page.

### Required Documentation Sections

1. **Physical closure and normalization**
   - map paper variables to AthenaK code units;
   - define heat-flux coefficients and caps;
   - identify cold-electron scope.

2. **Threshold-policy conventions**
   - describe parallel versus oblique firehose choices;
   - state exactly which convention MKS24 uses;
   - distinguish limiter occupancy from safety failure.

3. **Collision source-term ordering**
   - document the corrected split and why its timestep fractions are correct.

4. **Operator-validation workflows**
   - retain `quick`, `compare`, `amr`, and `full`;
   - explain what they do and do not establish.

5. **Paper-reproduction workflows**
   - provide copy-paste reduced smoke invocation;
   - provide HPC-standard invocation patterns;
   - show analysis regeneration from archived output.

6. **Result interpretation**
   - explain manifest and summary products;
   - explain physical occupancy, safety counters, spectra, transfer, and
     alignment products;
   - document run averaging windows and randomness.

7. **Known limitations**
   - uniform-grid paper baseline versus AMR robustness;
   - fluid microinstability surrogate versus kinetic physics;
   - cold-electron analogue limits;
   - any missing exact comparator data or digitization uncertainty.

### Documentation Build Gate

Build documentation with warnings treated as errors after each documentation
commit. If the repository uses a separate `gh-pages` publication branch,
verify the integration path there before claiming published documentation is
ready.

### Recommended Commit

```text
Document CGL-LF physical conventions and MKS24 reproducible workflows
```

## Frontier Execution, Budget, and Operational Reproducibility

### Purpose and Authority

The supercomputer execution target for this project is OLCF Frontier. The
authoritative system reference is the
[Frontier User Guide](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html).
Future agents must reread the relevant guide sections immediately before the
first Frontier build or submission in a campaign because installed module
stacks, known issues, and scheduling policies can change.

The following user-specified constraints override convenience:

```text
Project/account: AST207
Project run root: /lustre/orion/ast207/proj-shared/dfielding/CGL
Allowed QOS: debug only
Maximum testing budget: 1000 Frontier node-hours total across this effort
Concurrency: one debug-QOS job in any state at a time
```

All Frontier simulations for this CGL-LF project must place their output,
logs, manifests, restart products, and usage accounting beneath:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL
```

Do not write CGL-LF simulation products into the older `divb` directory or a
personal scratch path merely because a prior script used those locations.

### Debug-QOS Restriction and Its Consequence

The Frontier guide states that the `debug` QOS is for short non-production
debug tasks, that production work and job chaining using that QOS are
prohibited, that each user may have only one `debug` job in any state at a
time, and that requested walltime cannot exceed two hours.

Therefore:

- every Frontier job submitted by agents under this plan must include:

  ```bash
  #SBATCH -p batch
  #SBATCH -q debug
  #SBATCH -t HH:MM:SS   # never greater than 02:00:00
  ```

- agents must submit jobs sequentially and manually inspect completed results
  before submitting the next job;
- agents must not use Slurm dependencies, job arrays, automated job chaining,
  or background waves of `sbatch` submissions in `debug`;
- debug-QOS jobs may be used for build checks, GPU/MPI correctness,
  performance reconnaissance, reduced turbulence tests, restart validation,
  output/analysis verification, and campaign sizing;
- full paper-standard or paper-HPC production simulations must not be launched
  in `debug` while this restriction applies, because OLCF explicitly describes
  that as non-production access.

This distinction matters: the code may be brought to **production quality**
using rigorous debug-queue qualification plus documented evidence, while the
actual expensive published MKS24 production campaign requires a later explicit
user decision and an OLCF-compliant production-queue execution plan.

If the user later authorizes another QOS or partition for paper production,
update this document in a reviewable commit before submitting those runs.

### Frontier Hardware and Rank Mapping Baseline

According to the Frontier guide, a compute node contains four AMD MI250X
accelerators that appear to Slurm and ROCr as eight GPU-visible Graphics
Compute Dies (GCDs). Default core specialization leaves 56 CPU cores available
to user job steps. OLCF presents eight MPI ranks per node, one rank per GPU and
seven CPU cores per rank, as the common mapping case.

Use this as the initial multi-GPU validation mapping:

```text
ranks per node:       8
GPUs per rank:        1
CPUs per rank:        7
GPU binding:          closest
CPU binding:          threads
hardware threads:     1 per physical core
```

Always make the task count explicit in `srun`. Do not rely on
`--ntasks-per-node` alone. The guide documents cases where omitted explicit
task counts can produce an unintended single rank per node.

For the Kokkos HIP executable expected here, do not enable CPU OpenMP
threading by default. Use:

```bash
export OMP_NUM_THREADS=1
```

while still reserving `-c 7` CPU cores per GPU rank for runtime and host-side
work. If a later performance study intentionally enables OpenMP, it must be a
separate build and validation branch in the evidence, not an undocumented
launch change.

Small single-block quantitative pgens may not support eight MPI ranks. Use
them for CPU/local validation or one-rank GPU bring-up only. Use multiblock
paper-smoke or AMR configurations for eight-rank and multinode GPU/MPI
qualification.

### Frontier Directory Layout

Create and preserve the following structure:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL/
  repo/
    athenak-DF/                       # checked-out exact tested revision
  build/
    frontier-hip-<gitsha>-<stack>/    # immutable build per revision/toolchain
  inputs/
    archived/                         # copied submitted input decks if needed
  runs/
    <campaign>/
      <run-name>/
        manifest/
        output/
        restart/
        analysis/
  logs/
    build/
    slurm/
  accounting/
    frontier_node_hours.csv
    frontier_budget_summary.md
  scripts/
    build/
    submit/
    accounting/
```

Rules:

- Build from an exact git revision. Never run an executable whose source
  revision cannot be recovered from a committed SHA.
- Do not reuse a build directory after changing source revision, toolchain,
  CMake flags, or physical-model compile options.
- Keep submitted input decks and exact command overrides with each result
  bundle.
- Retain restart files until the run is accepted and archived or explicitly
  rejected.
- Keep large generated output under `runs/`, not inside the git checkout.

### Node-Hour Budget Policy

The total authorized Frontier **testing** budget for this CGL-LF effort is:

```text
1000 node-hours
```

This is a hard cumulative ceiling unless the user explicitly changes it.
According to OLCF accounting, a job is charged according to requested nodes
multiplied by actual time from entering the running state until exit; unused
cores within requested nodes do not reduce the charge, and failed jobs still
consume allocation.

#### Before Every Submission

An agent must:

1. Inspect the accounting ledger and compute remaining authorized test budget.
2. Calculate the conservative reservation for the proposed job:

   ```text
   reserved_node_hours = requested_nodes * requested_walltime_hours
   ```

3. Refuse to submit if:

   ```text
   consumed_actual_node_hours + reserved_node_hours > 1000
   ```

4. Confirm there is no existing `debug` job by this user in `PD`, `R`, or any
   other state before issuing `sbatch`.
5. Record the proposed job, scientific purpose, acceptance criterion, node
   count, time limit, and reservation in the campaign notes before it runs.

Recommended pre-submission commands:

```bash
squeue -u "$USER" -o "%.18i %.12P %.30j %.8u %.2t %.10M %.6D %R"
tail -n 20 /lustre/orion/ast207/proj-shared/dfielding/CGL/accounting/frontier_node_hours.csv
```

#### After Every Completed or Failed Job

Record actual allocation use before submitting another job. Use `sacct` to
retrieve the allocation record. Where supported, use `-X` to suppress step
records; otherwise select only the top-level allocation and do not sum
`.batch`, `.extern`, or `srun` child steps.

Example:

```bash
JOBID=<completed-job-id>
sacct -X -j "$JOBID" \
  --format=JobIDRaw,JobName,State,ExitCode,AllocNodes,ElapsedRaw,Start,End \
  -P
```

Compute:

```text
actual_node_hours = AllocNodes * ElapsedRaw_seconds / 3600
```

Append one row per Slurm allocation to:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL/accounting/frontier_node_hours.csv
```

Required columns:

```text
job_id,submitted_utc,completed_utc,campaign,run_name,purpose,state,exit_code,
nodes,requested_walltime,elapsed_seconds,reserved_node_hours,actual_node_hours,
cumulative_actual_node_hours,git_revision,executable_sha256,input_file,
output_dir,result,notes
```

Also update:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL/accounting/frontier_budget_summary.md
```

with:

- authorized budget;
- consumed actual node-hours;
- currently reserved but unfinished node-hours, if any;
- remaining budget;
- jobs failed or rejected and what was learned;
- next proposed job and why it is worth its cost.

Use `showusage` as an external allocation cross-check when appropriate, while
retaining the project-specific ledger as the authoritative record for the
1000 node-hour testing constraint.

#### Initial Budget Envelope

Agents must proceed incrementally rather than treating 1000 node-hours as a
target to spend. Reserve budget by decision stage:

| Activity | Initial upper envelope | Expansion rule |
| --- | ---: | --- |
| GPU build/configuration and one-node launch checks | 10 node-hours | Stop when environment is verified |
| Corrected core physics GPU/CPU and restart comparisons | 30 node-hours | Expand only after passing local CPU gates |
| Paper-pgen/analysis smoke and forcing checks | 40 node-hours | Expand only after outputs are interpretable |
| Reduced nonlinear science validation | 120 node-hours | Expand only after formal readiness review |
| Scaling, memory, I/O, and checkpoint reconnaissance | 100 node-hours | Expand only to answer a measured blocker |
| Contingency for defects and reruns | 200 node-hours | Every use requires finding ID |
| Unallocated reserve | 500 node-hours | User review required before consuming |

These envelopes total the 1000 node-hour maximum. They are not entitlements:
agents should use substantially less when reduced tests settle the question.

#### Recommended Accounting Helper

As part of the first Frontier workflow implementation, add a small helper
under:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL/scripts/accounting/
```

or a tracked repository script that is copied there with each campaign. It
should accept a completed Slurm job ID, query the top-level allocation record,
calculate actual node-hours, reject duplicate recording of a job ID, append
the CSV row, and regenerate the Markdown budget summary.

Minimal shell/Python pattern for future agents to turn into a reviewed helper:

```bash
#!/bin/bash
set -euo pipefail

CGL_ROOT=/lustre/orion/ast207/proj-shared/dfielding/CGL
LEDGER="${CGL_ROOT}/accounting/frontier_node_hours.csv"
JOBID="${1:?usage: record_frontier_usage.sh JOBID CAMPAIGN RUN_NAME PURPOSE RESULT}"
CAMPAIGN="${2:?missing campaign}"
RUN_NAME="${3:?missing run name}"
PURPOSE="${4:?missing purpose}"
RESULT="${5:?missing result}"

mkdir -p "${CGL_ROOT}/accounting"
[[ -f "${LEDGER}" ]] || printf '%s\n' \
  'job_id,submitted_utc,completed_utc,campaign,run_name,purpose,state,exit_code,nodes,requested_walltime,elapsed_seconds,reserved_node_hours,actual_node_hours,cumulative_actual_node_hours,git_revision,executable_sha256,input_file,output_dir,result,notes' \
  > "${LEDGER}"

if awk -F, -v id="${JOBID}" 'NR > 1 && $1 == id { found=1 } END { exit !found }' "${LEDGER}"; then
  echo "Job ${JOBID} already appears in ${LEDGER}" >&2
  exit 1
fi

sacct -X -j "${JOBID}" \
  --format=JobIDRaw,State,ExitCode,AllocNodes,TimelimitRaw,ElapsedRaw,Submit,End \
  -n -P > "${CGL_ROOT}/accounting/${JOBID}.sacct.txt"

python3 - "${LEDGER}" "${CGL_ROOT}/accounting/${JOBID}.sacct.txt" \
  "${JOBID}" "${CAMPAIGN}" "${RUN_NAME}" "${PURPOSE}" "${RESULT}" <<'PY'
import csv
import pathlib
import sys

ledger = pathlib.Path(sys.argv[1])
sacct_path = pathlib.Path(sys.argv[2])
job_id, campaign, run_name, purpose, result = sys.argv[3:8]

rows = [
    row for row in csv.reader(sacct_path.read_text().splitlines(), delimiter="|")
    if row and row[0] == job_id
]
if len(rows) != 1:
    raise SystemExit(f"Expected one top-level allocation row for {job_id}, found {len(rows)}")

record = rows[0]
state, exit_code = record[1], record[2]
nodes = int(record[3])
requested_minutes = int(record[4])
elapsed_seconds = int(record[5])
submitted_utc, completed_utc = record[6], record[7]
reserved = nodes * requested_minutes / 60.0
actual = nodes * elapsed_seconds / 3600.0

existing = list(csv.DictReader(ledger.open()))
cumulative = actual + sum(float(row["actual_node_hours"]) for row in existing)
if cumulative > 1000.0 + 1.0e-12:
    raise SystemExit(f"Recorded usage exceeds authorized budget: {cumulative:.6f} node-hours")

with ledger.open("a", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow([
        job_id, submitted_utc, completed_utc, campaign, run_name, purpose,
        state, exit_code, nodes, requested_minutes, elapsed_seconds,
        f"{reserved:.6f}", f"{actual:.6f}", f"{cumulative:.6f}",
        "", "", "", "", result, "",
    ])
print(f"Recorded {actual:.6f} node-hours; cumulative={cumulative:.6f}; remaining={1000.0-cumulative:.6f}")
PY
```

This is a starting pattern, not a substitute for review. Before relying on
it, verify Frontier's actual `sacct` output field format and make the helper
populate git revision, executable checksum, input, and output path from the
run manifest rather than leaving them blank. The helper must fail closed if it
cannot identify exactly one top-level allocation record.

### Frontier Module and Runtime Environment Policy

The supplied prior-working module sequence is the starting baseline:

```bash
module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime
```

This module stack is not permanent truth. Before its first use in a new
campaign:

1. Check that every requested module still exists.
2. Read current Frontier system-update and known-issue notices.
3. Record `module -t list`, `CC --version`, `cmake --version`, and the git SHA
   in the build manifest.
4. If the module versions must change, create a new build directory and run
   the GPU qualification tests again.

The prior run script also supplied MPI/GPU environment tuning:

```bash
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_OFI_NUM_CQ_ENTRIES=131072
export FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software
export HSA_XNACK=1
```

Treat these settings as an explicit run configuration to validate, not as an
unexamined optimization guarantee:

- retain `MPICH_GPU_SUPPORT_ENABLED=1` for GPU-aware MPI testing;
- record all environment settings in each manifest;
- validate correctness and memory behavior under the selected
  `HSA_XNACK`/managed-memory setting;
- do not change communication or memory-tuning variables between comparison
  runs without treating that as a new run configuration;
- do not load performance instrumentation packages in the baseline correctness
  build unless profiling is the stated purpose of the run.

### Improved Frontier Build Script Template

The following template uses the repository's existing `built_in_pgens`
configuration so the CGL-LF test and future registered paper pgens are
available in one executable. Save an adapted version beneath
`$CGL_ROOT/scripts/build/` and archive its output.

```bash
#!/bin/bash
set -euo pipefail

CGL_ROOT=/lustre/orion/ast207/proj-shared/dfielding/CGL
SRC_DIR="${SRC_DIR:-${CGL_ROOT}/repo/athenak-DF}"
BUILD_JOBS="${BUILD_JOBS:-16}"
STACK_TAG="cpe25.09-cce20-rocm6.4.2"

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime

test -d "${SRC_DIR}/.git"
GIT_SHA="$(git -C "${SRC_DIR}" rev-parse HEAD)"
GIT_SHORT="$(git -C "${SRC_DIR}" rev-parse --short=12 HEAD)"
BUILD_DIR="${BUILD_DIR:-${CGL_ROOT}/build/frontier-hip-${GIT_SHORT}-${STACK_TAG}}"
BUILD_LOG_DIR="${CGL_ROOT}/logs/build"
BUILD_MANIFEST_DIR="${CGL_ROOT}/runs/build-manifests/${GIT_SHORT}-${STACK_TAG}"
mkdir -p "${BUILD_DIR}" "${BUILD_LOG_DIR}" "${BUILD_MANIFEST_DIR}"

if ! git -C "${SRC_DIR}" diff --quiet --ignore-submodules -- ||
   ! git -C "${SRC_DIR}" diff --cached --quiet --ignore-submodules --; then
  echo "Refusing to build from tracked source modifications without a committed SHA." >&2
  exit 1
fi

{
  date -u +"created_utc=%Y-%m-%dT%H:%M:%SZ"
  hostname
  echo "source_dir=${SRC_DIR}"
  echo "git_revision=${GIT_SHA}"
  echo "build_dir=${BUILD_DIR}"
  module -t list 2>&1
  CC --version 2>&1 | head -n 4
  cmake --version | head -n 1
} | tee "${BUILD_MANIFEST_DIR}/environment.txt"

cmake -S "${SRC_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DPROBLEM=built_in_pgens \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_VEGA90A=ON \
  -DCMAKE_CXX_COMPILER=CC \
  2>&1 | tee "${BUILD_LOG_DIR}/configure-${GIT_SHORT}-${STACK_TAG}.log"

cmake --build "${BUILD_DIR}" --parallel "${BUILD_JOBS}" \
  2>&1 | tee "${BUILD_LOG_DIR}/build-${GIT_SHORT}-${STACK_TAG}.log"

ATHENA="${BUILD_DIR}/src/athena"
test -x "${ATHENA}"
sha256sum "${ATHENA}" | tee "${BUILD_MANIFEST_DIR}/athena.sha256"
"${ATHENA}" -c > "${BUILD_MANIFEST_DIR}/athena-config.txt"
cp "${BUILD_DIR}/CMakeCache.txt" "${BUILD_MANIFEST_DIR}/CMakeCache.txt"

printf 'Built executable: %s\nGit revision: %s\n' "${ATHENA}" "${GIT_SHA}"
```

Notes:

- The original example used `-DPROBLEM=cgm_cooling_flow_amr_metals`; that is
  not appropriate for CGL-LF testing. The CGL-LF quantitative pgen is already
  registered in `built_in_pgens`, and the future paper pgen should preferably
  be registered there as well.
- The original example loaded `perftools-base perftools` for every build.
  Keep instrumentation out of the correctness baseline and create a separate
  profiling build only when profiling is explicitly required.
- The repository's AMD/HIP build documentation requires the Kokkos HIP and
  VEGA90A options. The explicit Cray compiler wrapper and ZEN3 target above
  follow the supplied Frontier environment; if configuration or linking
  requires additional ROCm include/link flags, record that as a toolchain
  finding and rebuild in a new build directory rather than silently editing a
  proven build.
- Compilation itself must not launch high-concurrency work on shared login
  nodes without respecting OLCF guidance. Reduce `BUILD_JOBS` or compile in an
  allocated debug node when the build load warrants it, accounting for any
  allocated node time.

### Improved Frontier Debug-Run Slurm Template

This is a correctness and reduced-validation template. It intentionally uses
`debug` and therefore is not a paper-production submission template.

```bash
#!/bin/bash
#SBATCH -J cgl_lf_debug
#SBATCH -A AST207
#SBATCH -o /lustre/orion/ast207/proj-shared/dfielding/CGL/logs/slurm/%x.%j.log
#SBATCH -p batch
#SBATCH -q debug
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH --gpus-per-node=8
#SBATCH --threads-per-core=1

set -euo pipefail

CGL_ROOT=/lustre/orion/ast207/proj-shared/dfielding/CGL
CAMPAIGN="${CAMPAIGN:-qualification}"
RUN_NAME="${RUN_NAME:-cgl_lf_debug_${SLURM_JOB_ID}}"
ATHENA="${ATHENA:?Set ATHENA to an archived Frontier HIP executable}"
INPUT="${INPUT:?Set INPUT to the input deck for this run}"
ATHENA_WALLTIME="${ATHENA_WALLTIME:-00:25:00}"
RANKS_PER_NODE="${RANKS_PER_NODE:-8}"
CPUS_PER_TASK="${CPUS_PER_TASK:-7}"
NNODES="${SLURM_NNODES:?Missing SLURM_NNODES}"
NRANKS="$((NNODES * RANKS_PER_NODE))"
RUN_DIR="${CGL_ROOT}/runs/${CAMPAIGN}/${RUN_NAME}"
OUT_DIR="${RUN_DIR}/output"
MANIFEST_DIR="${RUN_DIR}/manifest"
mkdir -p "${OUT_DIR}" "${MANIFEST_DIR}" "${CGL_ROOT}/logs/slurm"

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_OFI_NUM_CQ_ENTRIES=131072
export FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software
export HSA_XNACK=1
export OMP_NUM_THREADS=1

test -x "${ATHENA}"
test -f "${INPUT}"
cp "${INPUT}" "${MANIFEST_DIR}/submitted_input.athinput"
sha256sum "${ATHENA}" > "${MANIFEST_DIR}/athena.sha256"

{
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
  echo "slurm_job_id=${SLURM_JOB_ID}"
  echo "campaign=${CAMPAIGN}"
  echo "run_name=${RUN_NAME}"
  echo "athena=${ATHENA}"
  echo "input=${INPUT}"
  echo "nodes=${NNODES}"
  echo "ranks=${NRANKS}"
  echo "ranks_per_node=${RANKS_PER_NODE}"
  echo "cpus_per_task=${CPUS_PER_TASK}"
  echo "athena_walltime=${ATHENA_WALLTIME}"
  module -t list 2>&1
  env | LC_ALL=C sort | grep -E '^(MPICH|FI_|HSA_|OMP_|ROCR|HIP|CRAY)=' || true
} > "${MANIFEST_DIR}/run_environment.txt"

RUN_ARGS=(-i "${INPUT}")
if [[ -n "${RESTART:-}" ]]; then
  test -f "${RESTART}"
  cp "${RESTART}" "${MANIFEST_DIR}/submitted_restart.rst"
  RUN_ARGS=(-r "${RESTART}")
fi

srun -N "${NNODES}" -n "${NRANKS}" --ntasks-per-node="${RANKS_PER_NODE}" \
  -c "${CPUS_PER_TASK}" --threads-per-core=1 --cpu-bind=threads \
  --gpus-per-task=1 --gpu-bind=closest \
  "${ATHENA}" "${RUN_ARGS[@]}" \
  -d "${OUT_DIR}" \
  -t "${ATHENA_WALLTIME}" \
  job/basename="${RUN_NAME}"

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" >> "${MANIFEST_DIR}/run_environment.txt"
```

Submission example for a reduced multiblock validation case:

```bash
CGL_ROOT=/lustre/orion/ast207/proj-shared/dfielding/CGL
ATHENA="${CGL_ROOT}/build/frontier-hip-<gitsha>-<stack>/src/athena"
INPUT="${CGL_ROOT}/repo/athenak-DF/inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput"
sbatch --export=ALL,ATHENA="${ATHENA}",INPUT="${INPUT}",CAMPAIGN=paper-smoke,RUN_NAME=paper_smoke_active_beta10 \
  "${CGL_ROOT}/scripts/submit/cgl_lf_debug.sbatch"
```

Before using this submission example, confirm the `paper` input exists and
that its meshblock count supports the selected rank count. Until the paper
pgen is implemented, use a current multiblock CGL-LF AMR input only for
hardware/runtime qualification and label the result accordingly.

### Restart and Walltime Handling in Debug

The Slurm time limit must not exceed:

```text
02:00:00
```

Set AthenaK's internal `-t` walltime shorter than the Slurm request so it has
time to emit final output or restart products. Use a margin of at least five
minutes for short jobs and ten minutes for near-two-hour jobs. For example:

| Slurm request | Maximum recommended AthenaK `-t` |
| --- | --- |
| `00:30:00` | `00:25:00` |
| `01:00:00` | `00:52:00` |
| `02:00:00` | `01:50:00` |

Do not automatically resubmit restart segments with dependencies or scripted
job chaining in `debug`. After a job completes:

1. record its node-hours;
2. inspect safety counters, state validity, log errors, and restart presence;
3. decide whether continuation is justified as a further validation segment;
4. submit the next job manually only if it remains within policy and budget.

### Frontier GPU Qualification Ladder

Frontier qualification must occur only after Phases A through C pass on the
ordinary CPU test path. It is a gate for production readiness on Frontier,
not a substitute for local physical correctness.

Run in the following order, stopping on the first unexplained failure:

| Step | Job size/time envelope | Purpose | Required comparison |
| --- | --- | --- | --- |
| G-001 | 1 node, at most 10 minutes | Confirm HIP executable launches and reports configuration | Build manifest, clean startup, device/MPI environment |
| G-002 | 1 node, at most 20 minutes | Run multiblock strict LF smoke on GPUs | Clean LF diagnostics and valid outputs |
| G-003 | 1 node, at most 30 minutes | Compare one-rank and eight-rank same-node result | Tolerance-based physical/history agreement |
| G-004 | 1 node, at most 30 minutes | Corrected finite-collision GPU regression | Agreement with analytic reference, not only CPU |
| G-005 | 1 node, at most 45 minutes | GPU restart continuation | Final-state agreement and clean resumed diagnostics |
| G-006 | 2 to 4 nodes, at most 60 minutes | MPI/GPU decomposition behavior on multiblock problem | Clean strict diagnostics, tolerance-based reduced histories |
| G-007 | 1 to 4 nodes, at most 90 minutes | Reduced paper-pgen active/passive/forcing validation after implementation | Complete manifests and finite analysis products |
| G-008 | Carefully approved small set | Timing, memory, I/O, and checkpoint sizing | Evidence for later production plan |

Required Frontier GPU comparison principles:

- Compare physical fields or reduced conserved/history quantities using
  tolerance-based metrics; do not expect bitwise CPU/GPU identity.
- Establish tolerances from measured convergence and reduction behavior, not
  convenience.
- Use identical model inputs, threshold policy, forcing seed, and output
  cadence when comparing decompositions or platforms.
- Keep strict diagnostics enabled for qualification runs.
- Treat any nonfinite state, floor activation, emergency-bound failure,
  unexplained energy discrepancy, or inability to restart as a stop-the-line
  finding.

### Frontier Performance Evidence Required for Production Quality

Performance optimization remains deferred until profiling demonstrates a real
problem. Performance **measurement** is not deferred: a production-quality
workflow must establish that intended runs are feasible and that output and
checkpoint strategies are safe.

After correctness qualification, record for selected reduced cases:

- elapsed walltime and simulated time advanced;
- zones or cell-updates per second where meaningful;
- requested and used nodes/GPUs;
- peak or observed memory evidence available from runtime/tools;
- output and restart size and write time;
- LF/STAGE metrics and limiter/cap activity relevant to workload;
- projected node-hours and storage for the standard MKS24 matrix.

Only then decide whether to optimize scratch allocation, repeated derived-state
construction, I/O cadence, meshblock layout, or other measured bottlenecks.
Any optimization requires rerunning the corrected physics and Frontier GPU
qualification gates.

### No Authorized Production Submission Template Yet

This plan intentionally provides only a `debug`-QOS validation template,
because the current instruction authorizes agents to use only `debug` and the
OLCF policy prohibits using it for production work. After qualification
evidence exists, future agents must prepare a proposed standard-production
campaign containing:

- required paper-standard cases;
- predicted runtime, node count, node-hours, and storage;
- selected executable SHA and validated module stack;
- checkpoint/restart strategy;
- analysis and archive workflow;
- requested non-debug queue/QOS consistent with current OLCF policy.

Submit that proposal to the user for explicit approval and revise this
document before launching any paper-production jobs.

## Phase H: Execution Tiers and Scientific Acceptance

### Tier 0: Routine CI

Purpose: protect core correctness and cheap regression behavior.

Include:

- analytic collision-time regression;
- threshold-policy parsing and one controlled limiter test;
- strict hard-bound negative test;
- existing collisionless LF case;
- one explicit-versus-STS case;
- FOFC;
- restart;
- AMR CPU smoke;
- MPI AMR decomposition test where the current CI already supports it.

No forced three-dimensional paper simulation belongs in routine CI.

### Tier 1: Local Paper Smoke

Purpose: verify pgen, forcing, diagnostics, workflow packaging, and analysis
execution.

Recommended setup:

```text
resolution: 16x16x32 or 32x32x64 for initial debugging
runtime:    short, enough to invoke forcing and outputs
cases:      active beta10, passive beta10, random beta10,
            limiter-frequency variation, threshold-policy comparison
```

Acceptance:

- successful startup and evolution;
- clean numerical safety diagnostics for strict smoke configurations;
- actual forcing activity;
- finite reduced histories;
- analysis scripts produce all requested smoke summaries and sample plots;
- active/passive modes produce distinguishable, interpretable metadata and
  outputs, without requiring statistically converged physical differences.

This tier may be run on Frontier only as a short `debug`-QOS qualification
job following the budget, ledger, and sequential-submission rules above.

### Tier 2: Reduced Scientific Validation

Purpose: establish qualitative physical behavior and catch major workflow or
diagnostic mistakes before expensive production.

Recommended setup:

```text
resolution: 48x48x96 or 96x96x192
runtime:    several outer-scale times, stated explicitly
cases:      active/passive beta10 Alfvenic,
            active beta10 random,
            active beta100 selected nu_lim scan,
            selected heat-flux sensitivity case
```

Acceptance:

- stable runs and reproducible analysis;
- threshold-policy choice is reflected in occupancy measurements;
- active/passive diagnostics show scientifically plausible differences;
- spectra and transfer code pass synthetic checks and generate finite products;
- energy/injection accounting is understood well enough to proceed.

Do not compare shortened reduced runs directly to steady-state paper statistics
without prominent qualification.

This tier may be executed through carefully limited Frontier `debug` runs only
when its purpose is validation and sizing, not production science. Every run
must remain within the node-hour ledger and formal readiness review process.

### Tier 3: Paper-Standard Production

Purpose: reproduce the principal CGL-LF results in MKS24.

This tier is scientifically required for a completed reproduction claim, but
it is **not currently authorized for execution by future agents** under the
debug-only Frontier instruction. OLCF prohibits production work in the
`debug` QOS. Agents must complete the lower-tier qualification, prepare a
costed production campaign proposal, obtain explicit user authorization for
an appropriate production-compliant submission mode, and revise this document
before submitting Tier 3 jobs.

Required baseline:

```text
resolution: 192x192x384
runtime:    t_final >= 10 L_perp/v_A
averaging:  steady-state window selected consistently with the paper,
            initially after t v_A/L_perp > 6 where applicable
```

Minimum case matrix:

| Group | Cases |
| --- | --- |
| Active Alfvenic beta scan | `beta0 = 1, 10, 100` |
| Passive comparison | At least matching Alfvenic cases needed for published comparisons |
| Random forcing comparison | Required beta values used in targeted figures |
| Heat-flux sensitivity | Runs needed for Figure 12 interpretation |
| Limiter-frequency scan | Beta 100 Alfvenic cases corresponding to Figure 13 |

Acceptance:

- matched numerical setup, threshold convention, forcing configuration, and
  run durations;
- clean numerical safety outcomes or fully explained deviations;
- complete manifests and archived analysis;
- quantitative figure-panel comparisons with uncertainty estimates from time
  variability and, where feasible, resolution checks;
- explicit note where raw reference data are unavailable and comparison relies
  on digitized published figures.

### Tier 4: Convergence and HPC Confirmation

Purpose: establish robustness of scientific conclusions.

Like Tier 3, this is a required later scientific campaign rather than an
authorized `debug`-QOS workload. Only reduced tests needed to determine
feasibility, correctness, or resource requirements may be run under the
current debug-only instruction.

Selected cases should be repeated at:

```text
96x96x192
192x192x384
384x384x768, resources permitting
```

Assess:

- inertial-range spectral trends;
- threshold-volume statistics;
- pressure-stress transfer;
- alignment PDFs;
- sensitivity to output cadence and averaging interval;
- sensitivity to threshold convention as an explicitly separate experiment.

### Data Management

Paper-grade outputs require a storage plan before execution. Each HPC campaign
must define:

- use of `/lustre/orion/ast207/proj-shared/dfielding/CGL` as its run root;
- requested output directory;
- expected raw snapshot volume;
- history and checkpoint cadence;
- restart policy;
- retention policy for raw versus reduced products;
- manifest and environment archival;
- how outputs are transferred back for analysis.

Do not launch long high-resolution runs without this plan.
Do not launch any Frontier test run without updating its projected and actual
node-hour accounting in the project ledger.

## Detailed Test and Validation Matrix

| ID | Test or run | Main question | Tier | Required before paper standard? |
| --- | --- | --- | --- | --- |
| P-001 | Uniform anisotropy, background collision only | Does collision evolve for exactly physical time `dt`? | CI | Yes |
| P-002 | Uniform threshold-active limiter relaxation | Does limiter scattering use intended time fraction? | CI/manual | Yes |
| P-003 | Threshold-policy state scan | Are parallel and oblique conventions distinguished correctly? | CI | Yes |
| P-004 | Strict hard-bound without backup correction | Does strict safety report real violations? | CI | Yes |
| N-001 | Parallel LF decay convergence | Is collisionless `q_parallel` damping accurate? | CI/manual | Yes |
| N-002 | Perpendicular LF decay convergence | Is collisionless `q_perp` damping accurate? | Manual | Yes |
| N-003 | Finite-collision LF asymptotics | Are collision-modified coefficients and splitting accurate? | Manual | Yes |
| N-004 | Rotated-field invariance | Is field projection directionally consistent? | Manual | Yes |
| N-005 | Extreme cap stress | Does free-streaming limitation remain stable and signed? | Manual | Yes |
| N-006 | Low-field stress | Is behavior near `bfloor` safe and documented? | Manual | Recommended |
| N-007 | STS timestep sweep | Does STS converge toward explicit reference? | Manual | Yes |
| R-001 | Restart with corrected finite collisions | Is representation and forcing/collision continuation reproducible? | CI/manual | Yes |
| R-002 | MPI AMR regression | Does decomposition preserve robustness metrics? | CI | Preserve |
| T-001 | Paper pgen initialization smoke | Are CGL fields and forcing initialized correctly? | Local | Yes |
| T-002 | Active/passive reduced turbulence | Do intended model branches execute and archive correctly? | Local/reduced | Yes |
| T-003 | Alfvenic/random forcing reduced run | Are forcing modes correct and observable? | Reduced | Yes |
| G-001 | Frontier HIP launch and device binding | Does the archived executable launch with expected GPU/MPI mapping? | Frontier debug | Yes for Frontier readiness |
| G-002 | Frontier GPU strict LF smoke | Does the GPU path maintain safety and output contracts? | Frontier debug | Yes for Frontier readiness |
| G-003 | Frontier one-rank/eight-rank and multinode comparison | Is GPU/MPI behavior decomposition-stable? | Frontier debug | Yes for Frontier readiness |
| G-004 | Frontier restart and corrected-collision verification | Does GPU execution preserve the corrected physical and restart contract? | Frontier debug | Yes for Frontier readiness |
| G-005 | Frontier runtime/I/O sizing | Is a future paper campaign operationally feasible? | Frontier debug | Yes before proposing production |
| A-001 | History analysis regression | Can archived histories regenerate summaries? | CI/manual | Yes |
| A-002 | Synthetic spectral analysis | Are binning/normalization conventions correct? | Manual/CI if cheap | Yes |
| A-003 | Synthetic transfer analysis | Is transfer sign/normalization defensible? | Manual/CI if cheap | Yes |
| A-004 | Synthetic alignment analysis | Is eigenvector/alignment logic correct? | Manual/CI if cheap | Yes |
| S-001 | MKS24 active/passive beta10 standard | Reproduce Figure 2/7/8 backbone | HPC | Core result |
| S-002 | MKS24 beta scan standard | Reproduce beta dependence | HPC | Core result |
| S-003 | MKS24 random forcing cases | Reproduce forcing dependence | HPC | Core result |
| S-004 | MKS24 heat-flux sensitivity | Reproduce Figure 12-type claim | HPC | Core result |
| S-005 | MKS24 limiter-frequency scan | Reproduce Figure 13-type claim | HPC | Core result |
| H-001 | Resolution confirmation selected cases | Are key conclusions robust to resolution? | HPC | Before final claim |

## Suggested Commit and Review Sequence

The following series is intentionally more granular than the prior feature
series. It prevents large paper infrastructure changes from hiding corrections
to physical behavior.

### Commit 1: Correct Collision Splitting

```text
Correct CGL-LF collision source-term splitting and add analytic regression
```

Deliverables:

- corrected collision timestep handling;
- independent analytic background-collision regression;
- corrected pgen reference oracle;
- updated explicit-versus-STS finite-collision test;
- concise documentation of update ordering.

Gate:

- routine CGL CPU tests;
- existing STS diffusion CPU tests;
- `git diff --check`.

### Commit 2: Add Paper Threshold Policy

```text
Add selectable CGL firehose threshold policies for MKS24 reproduction
```

Deliverables:

- shared threshold policy;
- input parsing and validation;
- parallel and oblique threshold tests;
- documentation and input-parameter reference.

Gate:

- CPU limiter tests;
- input fatal tests;
- operator full workflow as applicable.

### Commit 3: Repair Safety Diagnostic Semantics

```text
Make CGL-LF hard-bound diagnostics independent of corrective limiters
```

Deliverables:

- truthful safety counting;
- strict-mode expected-failure cases;
- summary/manifest semantics updated;
- documentation clarifying physical occupancy versus safety failure.

Gate:

- CPU and AMR CGL tests;
- workflow summary regeneration tests.

### Commit 4: Expand Accuracy Validation

```text
Add CGL-LF convergence, asymptotic, and directional accuracy validation
```

Deliverables:

- collisionless/finite-collision accuracy studies;
- rotated-field tests;
- cap stress and timestep comparisons;
- analysis tables/plots where suitable.

Gate:

- local validation suite;
- docs build if results are documented.

### Commit 5: Add Paper Pgen and Reduced Inputs

```text
Add MKS24 CGL-LF turbulence pgen and reduced smoke inputs
```

Deliverables:

- current-API paper pgen;
- active/passive setup;
- forcing integration;
- reduced smoke inputs;
- initial histories.

Gate:

- build default binary;
- run all paper-smoke cases;
- restart smoke if forcing state is restart-sensitive.

### Commit 6: Add Paper Reduced Diagnostics

```text
Add CGL-LF paper history diagnostics and archived summaries
```

Deliverables:

- threshold-volume histories;
- `C_B2`;
- energy/injection summaries;
- cap and collision activity;
- manifest extensions.

Gate:

- smoke workflows;
- finite history validation;
- summary regeneration from existing bundle.

### Commit 7: Add Paper Spatial Analysis

```text
Add CGL-LF paper PDF, spectra, transfer, and alignment analysis
```

Deliverables:

- offline analysis package;
- synthetic validation;
- figure product structure;
- documentation of definitions.

Gate:

- synthetic tests;
- reduced science bundle reanalysis;
- inspect generated plots.

### Commit 8: Add Paper Workflows and Documentation

```text
Document and automate tiered MKS24 CGL-LF reproduction workflows
```

Deliverables:

- `paper-smoke`, `paper-standard`, and analysis modes;
- runbook;
- result interpretation;
- storage/HPC guidance.

Gate:

- paper-smoke end-to-end;
- operator workflows unchanged;
- Sphinx build with warnings as errors.

### Campaign Work: Execute Standard and HPC Runs

Do not treat expensive generated result bundles as ordinary source commits.
Instead:

- execute approved run matrices;
- archive manifests, reduced diagnostics, and plot products;
- commit only curated small artifacts or result summaries if requested;
- document the exact result archive used for any paper claim.

Under the present user instruction, this campaign section defines required
future scientific work but does not authorize submission of production runs
through Frontier `debug`. A separate approved execution update is required
after the Frontier debug qualification ladder and resource estimate are
complete.

## Repository Validation Gates

### After Every Code Commit

Run at minimum:

```bash
git diff --check
```

Then run the narrowest relevant unit/regression tests for the changed
behavior.

### After Core Physics Corrections

Required:

```bash
cd tst
python run_tests.py cgl/cgl_landau_fluid
```

Also run retained STS tests and any corrected finite-collision local workflow
defined by the implementation.

### After Workflow and Documentation Changes

Required:

```bash
python3 scripts/cgl_lf_workflow.py quick
python3 scripts/cgl_lf_workflow.py compare
python3 scripts/cgl_lf_workflow.py amr
python3 scripts/cgl_lf_workflow.py full
```

Run `plot` and `summarize` on retained bundles. Build Sphinx documentation
using the repository-supported warnings-as-errors command.

### After MPI or AMR-Relevant Changes

Required:

```bash
cd tst
python run_test_suite.py --mpicpu
```

or the live repository equivalent if the test runner interface has changed.
Also rerun the existing divB AMR validation needed by the CGL AMR dependency.

### Before Review or Push

Run all relevant project gates, including:

```bash
git diff --check
python run_test_suite.py --cpu
python run_test_suite.py --mpicpu
python run_test_suite.py --style
```

Use the actual repository commands present at that time. The prior branch work
identified an existing style-wrapper problem involving
`docs/source/conf.py`; confirm whether it has been resolved rather than
silently accepting blocked CI.

## Files Future Agents Should Inspect First

### Physics and Evolution

```text
src/eos/cgl_physics.hpp
src/eos/cgl_mhd.cpp
src/eos/ideal_c2p_mhd.hpp
src/diffusion/cgl_landau_fluid.cpp
src/diffusion/cgl_landau_fluid.hpp
src/mhd/mhd.cpp
src/mhd/mhd.hpp
src/mhd/mhd_sts.cpp
src/mhd/mhd_tasks.cpp
src/driver/driver.cpp
src/mhd/rsolvers/hlle_cgl.hpp
```

### Existing Test and Workflow Foundation

```text
src/pgen/tests/cgl_landau_fluid.cpp
src/pgen/tests/cgl_fofc.cpp
src/pgen/tests/divb_amr.cpp
inputs/tests/cgl_lf_decay.athinput
inputs/tests/cgl_lf_limiter.athinput
inputs/tests/cgl_lf_amr_2d.athinput
inputs/tests/cgl_lf_restart.athinput
inputs/unit_tests/cgl_lf_*.athinput
tst/test_suite/cgl/test_cgl_landau_fluid_cpu.py
tst/test_suite/cgl/test_cgl_landau_fluid_mpicpu.py
scripts/cgl_lf_workflow.py
scripts/plot_cgl_lf_validation.py
scripts/generate_cgl_lf_eigenmode_inputs.py
```

### Documentation and Paper Reference

```text
docs/source/modules/cgl_landau_fluid.md
docs/source/modules/cgl_landau_fluid_validation.md
docs/source/reference/input_parameters.md
docs/arXiv-2405.02418v2/MKS24.tex
docs/arXiv-2405.02418v2/fig*.pdf
```

### Selective Donor Assets

Read through `git show` rather than merging:

```bash
git show origin/CGL-STS-LF:src/pgen/cgl_lf_paper.cpp
git show origin/CGL-STS-LF:scripts/analyze_cgl_lf_paper.py
git show origin/CGL-STS-LF:scripts/run_cgl_lf_paper_smoke.sh
git show origin/CGL-STS-LF:docs/cgl_lf_paper_validation_runbook.md
git show origin/CGL-STS-LF:docs/cgl_lf_paper_reproduction_plan.md
```

## Decision Log Required During Implementation

Future agents should maintain a short decision log in PR descriptions or in a
dedicated tracked implementation note if the work spans multiple PRs. Record:

1. Collision operator ordering and why it represents the intended physics.
2. Threshold-policy API, default choice, and paper input convention.
3. Definition of safety bounds versus physical instability occupancy.
4. Forcing implementation source and its reproducibility properties.
5. Output cadence and statistical averaging intervals.
6. Analysis definitions and normalization choices.
7. Which paper panels are exactly reproduced, compared qualitatively, or out
   of scope.
8. Storage/archive location for long-run result bundles.

Without this log, later agents will have to rediscover physical and numerical
assumptions from code and generated files.

## Open Technical Decisions

These issues must be resolved explicitly; do not bury them in implementation
details.

### Collision Split Form

Decision required:

- one full collision update per evolution cycle; or
- two symmetric half-collision updates.

Criterion:

- analytic physical-time correctness first;
- symmetry and numerical accuracy second;
- minimal architectural disruption third.

### Threshold Default

Decision required:

- preserve current oblique default for backward compatibility while paper decks
  explicitly choose parallel; or
- change the public default to parallel because the feature is not yet
  established and the principal reproduction target uses parallel.

Criterion:

- user compatibility and explicit documentation;
- no hidden behavioral change in archived runs.

### Emergency-Bound Definition Under Parallel Threshold

Decision required:

- define a stricter numerical overshoot bound beyond `beta Delta = -2`;
- use finite overshoot tolerance tied to regulator model;
- or use strict safety only for nonfinite/floor failures while reporting all
  physical threshold excursions separately.

Criterion:

- the diagnostic must not mark normal finite-rate paper physics as a numerical
  failure;
- it must still reveal unstable numerical excursions.

### Pgen Integration With Forcing

Decision required:

- reuse and extend the current turbulence forcing implementation; or
- implement a narrowly scoped paper pgen forcing hook.

Criterion:

- matches MKS24 forcing definition;
- restart and random-seed reproducibility;
- minimal duplication.

### Reference Data for Quantitative Paper Comparison

Decision required:

- obtain machine-readable original data from authors/archive if available; or
- digitize plotted reference curves with documented uncertainty.

Criterion:

- comparisons must state reference-data provenance and error bars.

### Paper Source and Figure Provenance on Frontier

Decision required:

- add the permitted MKS24 source/figure bundle to an appropriate tracked or
  archived project-data location; or
- add a documented retrieval/staging procedure that produces an immutable
  reference copy beneath the Frontier CGL run root.

Context:

- the local `docs/arXiv-2405.02418v2/` directory was untracked when this plan
  was authored;
- a pushed code branch alone will therefore not deliver that reference bundle
  to Frontier;
- simulations need not read the paper source, but figure comparison,
  provenance, and agent interpretation do require stable access to the exact
  reference version.

Criterion:

- every paper comparison campaign must identify the exact paper source/figure
  reference and retain either its checksum or its immutable archive path.

## Explicitly Deferred Work

Do not include the following in the initial physical-correction and MKS24
reproduction series unless necessary to fix a demonstrated blocker:

- performance optimization or kernel fusion before representative measurement;
- GPU performance optimization or performance claims before Frontier GPU
  correctness qualification and sizing measurements;
- LF composition with viscosity, resistivity, or scalar diffusion;
- primitive-prolongation AMR support;
- kinetic-scale mirror/firehose physics;
- evolving microinstability scattering models beyond the paper's specified
  fixed `nu_lim` runs;
- finite-electron-pressure EOS extensions unless a separately approved
  comparison requires them;
- broad redesign of CGL conserved-variable storage or Riemann solver layout.

Frontier GPU correctness qualification, MPI/restart behavior, runtime sizing,
and budget accounting are **not** deferred: they are required for a
production-readiness claim on the intended platform after core physics
corrections pass. The deferred item is optimization or broader performance
claim-making before measured evidence identifies a need.

Once MKS24 reproduction is established, the remaining feature extensions can
be separately planned and measured.

## Future-Agent Start Checklist

Use this sequence when beginning an implementation slice:

1. Confirm live branch, dirty files, and remote refs.
2. Read this plan section corresponding to the assigned phase.
3. Reopen the relevant MKS24 source lines and current code paths; do not rely
   solely on summaries.
4. Identify whether existing user or generated changes overlap your write
   scope; do not revert unrelated work.
5. Make the smallest coherent code change for the assigned physical or
   workflow objective.
6. Add or revise tests that would fail before your change and pass afterward.
7. Run narrow validation immediately, then broader gates appropriate to the
   phase.
8. Update documentation for new public behavior in the same reviewable change.
9. Record decisions, test commands, generated bundle paths, and any remaining
   uncertainty.
10. If a result invalidates an assumption in this document, register the
    finding and revise the plan before increasing scope or compute cost.
11. Before any Frontier submission, confirm `debug`-QOS compliance, the
    single-job rule, run-root placement, and remaining node-hour budget.
12. After every Frontier job, record actual node-hours and inspect the result
    before submitting another job.
13. Do not claim paper reproduction until the required standard runs and
    analyses have been executed and archived.

## Definition of Completion

The CGL-LF implementation can be described as accurate, reliable, usable, and
capable of reproducing the MKS24 CGL-LF numerical study only when all of the
following are true:

### Physical Correctness

- Collision evolution uses a documented, analytically validated physical
  timestep.
- Paper runs use the MKS24 firehose and mirror threshold conventions
  explicitly.
- Heat-flux closure and caps are verified against the stated paper equations.
- Safety diagnostics distinguish numerical failure from expected limiter
  physics.

### Numerical Reliability

- Corrected analytic, convergence, cap, directional, restart, AMR, and MPI
  tests pass.
- Strict smoke runs remain free of floors, nonfinite states, nonpositive
  pressures, and defined emergency-bound failures.
- Operator validation remains reproducible through archived workflow bundles.

### Frontier Operational Qualification

- A committed revision is built for Frontier with an archived module/CMake
  manifest and executable checksum.
- GPU/MPI strict smoke, corrected finite-collision behavior, decomposition
  comparison, and restart qualification pass using the documented launch
  mapping.
- Timing, memory, I/O, checkpoint, and storage measurements support a costed
  production-campaign proposal.
- All authorized Frontier test jobs run only under `debug`, remain beneath the
  cumulative 1000 node-hour limit, and are recorded in the project ledger.
- Paper-production jobs are not submitted under `debug`; an authorized,
  policy-compliant execution plan is approved and documented before such work.

### Usability and Reproducibility

- Users can run operator validation separately from paper reproduction with
  unambiguous commands.
- Paper workflow manifests record every model and forcing choice necessary for
  reruns.
- Analysis can be regenerated from retained outputs without rerunning
  simulation jobs.
- Sphinx documentation explains both direct use and the paper reproduction
  procedure.

### Paper Reproduction

- The driven active/passive turbulence infrastructure matches the paper setup.
- The required beta, forcing, heat-flux, and limiter-frequency case matrices
  have been run at appropriate resolution and duration.
- Diagnostics required for Figures 2 through 13 that are within CGL-LF scope
  are produced and compared quantitatively.
- Non-reproducible kinetic-only comparisons and remaining modeling limitations
  are stated plainly.
- Result bundles, analysis products, reference-data provenance, and
  interpretation are archived for review.

### Adaptive Review and Release Decision

- The findings register is current, and every blocking/high-severity finding
  is resolved or explicitly identified as a release blocker.
- Every tolerance revision, model-policy change, workflow change, and
  Frontier-environment change is justified with retained evidence.
- A final readiness report maps retained evidence to each production gate and
  states exactly which configurations are supported.
- Any remaining scientific or operational limitation is documented before
  users are directed to run the feature.

## Immediate Next Three Implementation Tasks

No future agent should begin by porting the large paper pgen. The immediate
critical path is:

1. **Correct collision-time integration.**
   Trace the source-term task graph, implement a physically correct update
   duration, and replace the self-referential finite-collision oracle with an
   independent analytic regression.

2. **Implement explicit instability-threshold policy and truthful safety
   diagnostics.**
   Add the MKS24 parallel-firehose option, retain an explicit oblique option,
   separate physical occupancy from emergency failure, and cover both in
   tests and manifests.

3. **Strengthen accuracy validation before nonlinear paper work.**
   Add convergence, finite-collision asymptotic, and directional tests so that
   the eventual turbulence program rests on a defensible corrected numerical
   operator.

Only after these three tasks pass review should work begin on the MKS24 driven
turbulence pgen, diagnostics, and paper-grade execution workflows.
