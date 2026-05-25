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
   reproduction program for Majeski, Kunz, and Squire 2024 (MKS24), using a
   checksum-recorded staging copy of the official arXiv `2405.02418v2`
   source instead of assuming paper contents are vendored in the checkout.
3. Establish a reproducible workflow for paper-grade driven-turbulence runs,
   their analysis products, and their limitations.

Implementation began on 2026-05-24. Phases A and B have working-tree
implementations and focused CPU evidence recorded below. Phase D now has an
active-Delta reduced turbulence initializer, repaired deterministic forcing,
RK-companion-integrated reduced forcing ledger with focused CPU oracles, and a passing
active/passive `paper-smoke` workflow. Reduced passive-Delta flow
decoupling is now regression-tested and pinned paper-source staging has been
exercised. Paper forcing inputs now state their physical mode shell and
`k^-2` power contract; standard and limiter-scan input definitions now exist
behind an explicit execution guard. Snapshot analysis foundations are
implemented, including manifest-selected time-window aggregation and
operator-face heat-flux-cap activity, local-field velocity-gradient/strain
products, manifest-qualified cell-centered heat-flux smoothing proxies, and
snapshot pressure-work/heat-flux time-quadrature estimates, and generic
diagnostic figure rendering. A local `accuracy`
workflow now retains collisionless resolution/timestep, finite-collision,
extreme-cap, rotated-field, and low-field evidence. Paper forcing now retains
restartable cumulative RK-integrated applied source work for active global
energy-residual analysis. LF histories now also retain restartable, RKL2-applied capped-face
heat-flux contractions with fine-side ownership on AMR interfaces. Optional
pressure-work histories now retain RK-applied total and anisotropic CGL
traction work after the same AMR correction as momentum. Optional
reference-curve manifests now provide checksum- and uncertainty-qualified
comparison plumbing. Long-time energy qualification, steady nonlinear
convergence, populated reference comparison,
production interpretation of the applied local-work ledgers, Frontier
production execution, and production-run phases remain
open. A tracked Frontier debug-campaign
utility and HIP build script now
implement offline-verifiable root, QOS, walltime, sequential-submission, and
budget-ledger policy enforcement. Live Frontier debug qualification now
retains corrected GPU restart, cross-node decomposition, reduced
paper-smoke evidence through G-007, post-F-026 pressure-work GPU restart and
AMR/MPI comparisons, scheduler-launched MPI CPU AMR comparison, an explicit
AMD/HIP managed-memory-policy confirmation, a partial G-008 sizing baseline,
an opt-in reduced shared-MPI-I/O timing follow-up, and a standard-layout
startup sizing measurement with a preliminary storage envelope; it does not
constitute paper-scale production evidence. The
current Sphinx runbook builds with warnings
treated as errors, and the repository style gate passes in an isolated
environment containing the declared documentation and test dependencies.

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

The earlier local paper source directory and the older planning note were
reported as untracked when the plan was authored:

```text
?? docs/arXiv-2405.02418v2/
?? next_steps_forCGL.md
```

Do not accidentally stage, delete, or rewrite either untracked path unless the
user explicitly asks for it. The implemented reference path is now the
checksum-recorded staging procedure described below.

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
Three initial findings define the corrected core and the remaining work:

1. **Collision splitting was physically wrong.** This working tree now gives
   each of the two LF post-sweep collision calls its half-sweep duration and
   checks uniform collisional relaxation against an independent analytic
   oracle.
2. **The baseline firehose activation threshold did not match MKS24
   production.** This working tree now exposes explicit `oblique` and
   `parallel` policies; reduced paper smoke inputs select the MKS24
   `beta Delta = -2` convention.
3. **The validation workflow still does not implement the full MKS24
   experiment.** It now includes reduced active and passive forcing smoke
   cases and snapshot-diagnostic foundations, but not long-time paper forcing
   calibration, figure/reference comparisons, execution of the limiter scan,
   or standard-resolution execution.

The work must therefore proceed in this order:

1. Correct and validate physical update semantics.
2. Make threshold conventions explicit and reproducible.
3. Strengthen diagnostic truthfulness and numerical accuracy tests.
4. Add a current-API paper problem generator and driven-turbulence workflows.
5. Add paper observables and figure-generation analysis.
6. Execute tiered smoke, convergence, standard, and HPC reproduction runs.

Do not begin performance tuning, new operator composition support, Frontier
GPU qualification, or publication-level claims until remaining numerical
accuracy gates and Phase D blockers are resolved. Frontier GPU correctness and
operational validation are required before any production-readiness claim.

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
| Collision source-term timing | Independent analytic regression and operator-ordering documentation | Yes | Implemented 2026-05-24; focused CPU test passed, broader gates pending |
| Instability-threshold policy | Explicit MKS24 policy, alternative-policy tests, archived input choice | Yes | Implemented 2026-05-24; CPU policy test and active paper-smoke manifest passed |
| Safety diagnostic truthfulness | Strict tests independent of backup correction | Yes | Implemented 2026-05-24; focused CPU test passed, broader gates pending |
| CPU numerical accuracy | Convergence, asymptotic, cap, directional, restart, MPI/AMR checks | Yes | Partial 2026-05-25; archived `accuracy-v2` bundle covers N-001 through N-007 locally, the AMR workflow passed, and the full CPU suite passed (`218 passed, 15 skipped`); focused CGL tests now pass (`30 passed`) with corrected RK forcing work, restartable AMR-corrected applied CGL pressure work, and the F-033 algebraic hard-wall policy; Frontier Kokkos-Serial/MPI jobs `4659053`/`4659141` pass the one-rank/four-rank AMR oracle at `rtol = atol = 1.0e-12`, although direct invocation in the local shell remains unavailable because it lacks `mpirun` |
| Frontier GPU numerical equivalence | CPU/GPU comparison suite, MPI/GPU restart checks, strict diagnostics | Yes for Frontier use | Partial 2026-05-25; immutable HIP/MPI builds archived; G-001 through G-007 pass their retained reduced gates; post-F-026 jobs `4658072`/`4658163` pass restart comparison for `lf_cpwrk`/`lf_cawrk`, `4658191`/`4658283` pass one-rank/four-rank GPU AMR comparison, and policy-confirmation job `4659663` is identical to its pre-policy reference; production-scale execution remains open |
| Paper pgen and forcing fidelity | Input/pgen review, forcing metadata, restartable reduced smoke | Yes for MKS24 claims | Partial 2026-05-25; active/passive/random reduced smoke, forcing restart, OU cadence, multi-cycle RK-companion source-work checks, passive flow-decoupling checks, and Frontier G-007 reduced-paper matrix passed; paper-grade statistical calibration open |
| Paper observables and analysis | Synthetic analysis tests and archived reduced-run products | Yes for MKS24 claims | Partial 2026-05-25; reduced histories, restartable RK-integrated applied forcing work/global active residual, RKL2-applied heat-flux contractions with fine-side AMR ownership, AMR-corrected RK-applied total/anisotropic CGL pressure-work ledgers, operator-face cap counters, windowed snapshot PDF/spectral/transfer/alignment/local-strain products, deduplicated threshold-volume history curves, manifest-qualified proxies and cadence-limited estimates, generic figures, and checksum/uncertainty-qualified reference-curve comparison plumbing implemented; populated panel comparisons and production qualification remain open |
| Standard MKS24 results | Required cases, durations, manifests, figure comparisons | Yes for reproduction claim | Inputs and guarded workflow modes defined 2026-05-24; no paper-scale runs or figure comparisons executed |
| Operational workflow | Frontier scripts, budget ledger, failure recovery, storage plan | Yes for Frontier use | Partial 2026-05-25; immutable HIP/MPI and Kokkos-Serial/MPI builds through revision `3e3c206b` are archived; sequential debug validation through G011 standard-layout startup MPI-I/O sizing and failed G012 finite-rate nonlinear hard-wall validation is retained and recorded (`0.273892` node-hours total); F-033 passes local constrained-model continuation through `t = 2.0`, while a corrected Frontier GPU rerun remains required before runtime/node-hour costing or an approved production proposal |
| User documentation | Sphinx build and accurately scoped runbook | Yes | Implemented for current functionality 2026-05-25; Sphinx warnings-as-errors and repository style suite pass in an isolated validation environment; future campaign results must still be documented when executed |
| Performance suitability | Representative timing/memory/I/O evidence; no uninvestigated prohibitive bottleneck | Yes for production use | Blocked 2026-05-25; G-008 records reduced debug-scale evidence, G010b records reduced shared-MPI-I/O timing, and G011 measures startup-only standard-layout memory, per-file size, and shared-MPI-I/O timing for a preliminary storage envelope; G012's finite-rate model aborts before its analysis window, while F-033's constrained model passes locally, so representative GPU late-time runtime/node-hour costing awaits corrected reduced GPU qualification |

The final readiness report must distinguish:

- **production ready for supported CGL-LF use**;
- **MKS24 reproduction complete**;
- **validated on Frontier GPU hardware**;
- **future extensions not yet supported**.

None of these phrases should be used as a substitute for the others.

## Current Readiness Decision (2026-05-25)

This is the current final readiness report for the implemented and retained
evidence on this branch. It is a blocked release decision, not an
authorization to execute paper-production jobs.

| Claim | Decision | Retained supporting evidence | Blocking work |
| --- | --- | --- | --- |
| Local operator and reduced-workflow implementation | Supported for the exercised local scope only | Full CPU suite (`218 passed, 15 skipped`), focused CGL suite (`30 passed`), retained `accuracy-v2`, `20260525-paper-smoke-pressure-work-v1`, `20260525-paper-convergence-pressure-work-v1`, and `20260525-amr-pressure-work-v1` bundles, F-018 through F-026 regression evidence, F-033 local exact-state hard-wall evidence, and scheduler-launched MPI CPU AMR jobs `4659053`/`4659141` | Production/statistical comparison gates remain open |
| Validated on Frontier GPU hardware | Supported for the archived reduced qualification cases through F-026 only | Immutable HIP/MPI builds; G-001 through G-007 strict/restart/decomposition/reduced-paper records; post-F-026 restart jobs `4658072`/`4658163`, GPU AMR jobs `4658191`/`4658283`, managed-policy confirmation job `4659663`, G010b reduced shared-MPI-I/O timing, and G011 standard-layout startup sizing | G012 job `4662477` validates that the pre-F-033 finite-rate hard-wall attempt fails strict near `t = 0.802`; F-033 passes locally but must be rerun on Frontier GPU before paper-scale execution or representative late-time runtime/node-hour costing |
| Production ready for supported CGL-LF use | Not established | Debug policy/accounting tooling, reduced GPU qualification, MPI CPU AMR evidence, shared-MPI-I/O measurements, preliminary defined-matrix storage envelope, and RK/AMR-consistent applied hyperbolic pressure traction histories exist | F-033 corrected hard-wall Frontier GPU qualification, representative runtime/node-hour costing, an approved production proposal, and production/statistical qualification remain required |
| MKS24 reproduction complete | Not established | Guarded standard/limiter input matrices, paper analysis products, and checksum/uncertainty-qualified reference-curve comparison interface exist | Populated MKS24 reference curves, long-time statistical calibration, authorized standard/limiter runs, quantitative panel comparisons, and archived scientific interpretation remain required |

The final audit specifically rejects two tempting overclaims:

- Original G-008 reports useful startup-scale timing, memory, output, and
  checkpoint sizes but disabled MPI-I/O timers. Follow-up G010b measures
  reduced shared-binary/restart MPI-I/O and exposes that its deliberate
  `dcycle = 1` restart cadence writes `4.654` GB in only `t = 0.01`.
  G011 then directly measures the standard layout and permits a preliminary
  storage envelope under the defined cadence, but its startup-length
  duration still cannot support runtime or node-hour costing.
- G012 attempts the missing reduced nonlinear hard-wall evidence with the
  finite-rate `limiter_nu_coll = 1.0e10` interpretation, but job `4662477`
  aborts near `t = 0.802` when a split-stage state crosses the strict
  firehose emergency bound. A matching local reduced probe with uncapped
  RKL2 STS reproduces a `hard_bound=1` abort near `t = 0.946`, and
  stage-context diagnosis from the screened state locates the failure at
  `post stage=1/5`.
- F-033 replaces the paper hard-wall decks' finite-rate pressure relaxation
  with an opt-in energy-preserving algebraic threshold projection. Its
  checksum-qualified local same-state continuation passes from `t = 1.5`
  through `t = 2.0` with zero strict safety counters and active `lf_hwproj`.
  Frontier GPU requalification remains mandatory before timing or production
  claims.
- The implemented `lf_cpwrk`/`lf_cawrk` histories supply the applied
  hyperbolic CGL pressure-traction identity for enabled decks by contracting
  stage velocity with the AMR-corrected traction divergence through the same
  RK recurrence as the state. This resolves the local identity gap, but does
  not itself qualify production statistics or paper-panel comparisons.

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
| F-001 | 2026-05-24 | Phase A | Collision updates covered two full timesteps per LF cycle | `src/mhd/mhd_sts.cpp`, `inputs/tests/cgl_lf_collision_relaxation.athinput`, focused CPU run | Blocker | Pass `dt_sweep = dt/2` to each LF collision call | `test_cgl_lf_background_collision_advances_one_physical_timestep` | Implemented; focused CPU passed |
| F-002 | 2026-05-24 | Phase B | Ordinary firehose threshold differed from MKS24 production convention | `src/eos/cgl_physics.hpp`, `inputs/tests/cgl_lf_firehose_policy.athinput`, focused CPU run | Blocker | Add `cgl_firehose_threshold = oblique/parallel`; MKS24 decks must state `parallel` | `test_cgl_lf_firehose_threshold_policies_are_distinct` | Implemented; active smoke deck and manifest select `parallel` |
| F-003 | 2026-05-24 | Phase B | Hard-bound diagnostic depended on backup correction enablement | `src/diffusion/cgl_landau_fluid.cpp`, focused CPU run | High | Count emergency bounds independently of correction policy | `test_cgl_lf_strict_hard_bound_is_reported_without_backup_correction` | Implemented; focused CPU passed |
| F-004 | 2026-05-24 | Phase D | Built-in turbulence forcing ignored documented seed selection, mapped type-1 parallel/perpendicular axes inconsistently with `B0 || z`, and added momentum without nonrelativistic energy work | `src/srcterms/turb_driver.cpp`, `src/utils/random.hpp`, `inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput`, `paper-smoke` bundle | Blocker | Add restartable `rseed`, define type 1 as `z`-guide-field Alfvenic forcing, add RK-consistent ideal/CGL source work, archive forcing metadata | `test_cgl_lf_paper_active_alfvenic_smoke_injects_energy_without_parallel_force`, seed and forcing-restart tests | Implemented for reduced active smoke mechanics; focused CPU passed |
| F-005 | 2026-05-24 | Phase D | Existing `mhd/passive=true` removed anisotropic terms from central momentum fluxes but still used CGL pressure-dependent HLLE and CFL signal speeds, allowing diagnostic anisotropy to alter MHD-like flow | `src/mhd/rsolvers/hlle_cgl.hpp`, `src/mhd/mhd_newdt.cpp`, `src/pgen/tests/cgl_lf_paper.cpp`, `inputs/cgl_lf_paper/cgl_lf_paper_smoke_passive_beta10.athinput` | Blocker | Use isothermal-MHD signal speeds for passive fluxing/timestep selection; require explicit matching passive pgen/EOS selection; add reduced flow-decoupling regression and smoke deck | `test_cgl_lf_paper_passive_delta_has_no_anisotropic_flow_feedback`, `test_cgl_lf_paper_passive_delta_must_match_eos_mode` | Implemented; focused CPU passed; long-time active/passive statistics remain open |
| F-006 | 2026-05-24 | Phase D | The plan named a local MKS24 source/figure bundle that is absent from the live worktree and should not be copied into the code branch without separately establishing redistribution permission | official arXiv `2405.02418v2` source, `scripts/stage_cgl_lf_mks24_reference.py`, staged `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/manifest.json` | High | Stage the pinned official source in ignored result storage, safely extract it, record archive/file SHA-256 checksums, and attach checksum provenance during paper analysis | Staging utility executed with archive SHA-256 `6b58e1ed01585c9820659a01b7549a08a09b7ea3fcbb28dbb96205259e82da8e`; `paper-convergence-v1/analysis/reference_provenance.json` | Implemented for reference provenance; quantitative reference curves/comparisons remain open |
| F-007 | 2026-05-24 | Phase D | `AddForcing` was chained after `InitializeModes` before the integrator and inserted before each RK update; every call advanced and applied the OU force with the full mesh timestep, so an `rk2` cycle invoked three full-duration pushes | `src/srcterms/turb_driver.cpp`, `src/mesh/meshblock_pack.cpp`, `src/driver/driver.cpp`, `docs/source/modules/srcterms.md` | Blocker | Advance OU state once per cycle in `UpdateForcing`; apply a fixed acceleration after each `RKUpdate` with `beta dt` primitive-state work and energy-preserving nonrelativistic mean-momentum projection | `test_cgl_lf_paper_forcing_ou_state_advances_once_per_cycle`, `test_cgl_lf_paper_active_alfvenic_smoke_injects_energy_without_parallel_force` | Implemented; focused CPU passed; long-time paper calibration remains open |
| F-008 | 2026-05-24 | Phase D | Reduced paper forcing selected integer-index shells and type-1 anisotropic exponents; the random override also inherited the `5/3` isotropic default instead of the documented physical `abs(k)` shell with `k^-2` power | staged `MKS24.tex:463-469`, `src/srcterms/turb_driver.cpp`, `inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput`, `scripts/cgl_lf_workflow.py` | Blocker | Add opt-in physical-shell/common-spectrum forcing controls, select them explicitly in paper inputs, archive them in manifests, and complete valid type-1 phase coefficients | `test_cgl_lf_paper_physical_forcing_shell_requires_positive_unit` plus `paper-smoke` execution | Implemented for reduced forcing contract; long-time spectral calibration remains open |
| F-009 | 2026-05-24 | Phase E | CGL primitive `mhd_w_bcc` snapshots exported `p_parallel` through the legacy `eint` field but omitted `p_perp`, making pressure-anisotropy snapshot products impossible | `src/outputs/basetype_output.cpp`, binary analysis probe, standard `bin` snapshot decks | Blocker | Retain legacy `eint = p_parallel` and append explicit `p_perp` for CGL primitive outputs | `test_cgl_lf_paper_snapshot_analysis_uses_both_pressures` validates binary output and analyzer products | Implemented; focused CPU and binary analyzer regression passed |
| F-010 | 2026-05-24 | Phase E | The LF closure applied heat-flux caps at operator faces but exposed no retained cap-activity counters for paper-window interpretation | `src/diffusion/cgl_landau_fluid.cpp`, `src/outputs/history.cpp`, `paper-smoke-core-v12` analysis bundle | High | Retain cumulative face counts and parallel/perpendicular pre-cap ratios above `q_max` and `10*q_max`; convert interval differences to fractions | `test_cgl_lf_heat_flux_cap_face_activity_is_reported` and windowed `paper-analyze` smoke result | Implemented; cap oracle and end-to-end smoke passed |
| F-011 | 2026-05-25 | Phase C | Single-case operator checks did not retain the required convergence, collisional-asymptotic, rotated-field, extreme-cap, or low-field evidence | `scripts/cgl_lf_workflow.py`, `src/pgen/tests/cgl_landau_fluid.cpp`, `inputs/unit_tests/cgl_lf_rotated_decay.athinput`, `20260524-accuracy-v2` bundle | High | Add archived `accuracy` workflow tables and pgen checks for directional decay and magnetic-floor shutdown | `test_cgl_lf_low_field_faces_disable_transport_cleanly`; `accuracy-v2` run | Implemented; focused CPU and 33-case accuracy bundle passed |
| F-012 | 2026-05-25 | Phase E | Snapshot analysis omitted MKS24 local-field velocity-gradient/strain products and produced no inspectable paper diagnostic figures | `scripts/analyze_cgl_lf_paper.py`, `scripts/plot_cgl_lf_paper.py`, binary analysis probe, `paper-smoke-core-v12` | High | Add local projected flow products and generic renderer invoked by `paper-analyze`; keep reference comparison separately gated | Expanded synthetic analysis and retained-bundle rendering | Implemented for diagnostic figures; paper-panel/reference comparison remains open |
| F-013 | 2026-05-25 | Phase F | The workflow table required a reduced `paper-convergence` path but only smoke and guarded production cases were executable | `scripts/cgl_lf_workflow.py`, `20260524-paper-convergence-v1` bundle | High | Add short resolution, heat-flux-strength, and threshold-policy variants with binary snapshots and analysis-window metadata | Six-case reduced workflow followed by `paper-analyze`/`paper-summary` | Implemented for startup/product-path qualification; steady-state convergence remains open |
| F-014 | 2026-05-25 | Phase C | Truthful hard-bound monitoring exposed that `cgl_lf_limiter_heat_flux_suppression` initialized exactly on the mirror emergency bound while requesting strict admissibility | `inputs/unit_tests/cgl_lf_limiter_heat_flux_suppression.athinput`, failed `20260524-final-full` and passing `20260524-final-full-v2` bundles | High | Keep strict checking and move the oracle state to `p_perp - p_parallel = 0.9 B^2`, above the `0.5 B^2` physical threshold but below the `1.0 B^2` emergency bound | Corrected suppression oracle and regenerated `full` plots/summary | Implemented; full workflow passed |
| F-015 | 2026-05-25 | Frontier operations | The plan prescribed debug-only root/budget/sequential-submission controls but provided templates without an executable fail-closed ledger or preparation tool | `scripts/frontier/cgl_lf_frontier.py`, `scripts/frontier/build_cgl_lf_frontier.sh` | High | Add tracked immutable-build and debug-campaign tooling that reserves budget, generates an inspectable batch script, checks an empty debug queue before manual submission, records one top-level `sacct` allocation from its run manifest, and rejects production decks | `python3 scripts/frontier/cgl_lf_frontier.py self-test`; `bash -n scripts/frontier/build_cgl_lf_frontier.sh` | Implemented; sequential live Frontier debug submission/accounting exercised through job `4661434` with `0.085003` node-hours recorded |
| F-016 | 2026-05-25 | Phase G | Required style validation could not pass because the Sphinx configuration used nonconforming indentation and the active system Python lacked lint/documentation dependencies | `docs/source/conf.py`, temporary validation environment built from `docs/requirements.txt` plus `flake8`, style/Sphinx logs | Medium | Correct the configuration indentation and run style/Sphinx through an isolated dependency-complete environment | `python run_test_suite.py --style`; `sphinx-build -W --keep-going -b html source _build/html` | Implemented; style wrapper and warning-strict documentation build passed |
| F-017 | 2026-05-25 | Phase E | Heat-flux cap occupancy alone did not quantify the transport strength, while an exact face-power reduction needs a decomposition-independent discrete accounting design | `scripts/analyze_cgl_lf_paper.py`, `scripts/plot_cgl_lf_paper.py`, `scripts/cgl_lf_workflow.py` | High | Add a manifest-qualified cell-centered reconstruction of the implemented LF closure as an explicitly non-budget heat-flux smoothing proxy; retain exact applied-face accounting as an open gate | Synthetic positive perpendicular-conduction/zero-parallel proxy check through `test_cgl_lf_paper_snapshot_analysis_uses_both_pressures` | Implemented for snapshot proxy; applied accounting with fine-side AMR ownership added by F-019 |
| F-018 | 2026-05-25 | Phase E | Retained `force_pwr` sampled instantaneous power but did not preserve the exact net RK stage forcing source, including zero-net-momentum projection, in forced paper runs or through restart | `src/srcterms/turb_driver.cpp`, `src/outputs/restart.cpp`, `src/pgen/tests/cgl_lf_paper.cpp`, `scripts/analyze_cgl_lf_paper.py` | High | Add opt-in globally reduced cumulative `force_work`, checkpoint it for paper decks, and compare active-case interval work to conserved `tot-E`; do not give passive-Delta the active budget interpretation | Expanded forcing identity/restart regressions and reduced paper analysis bundle | Implemented and corrected under F-025; long-time/production residual tolerance and scientific interpretation remain open |
| F-019 | 2026-05-25 | Phase E | The snapshot heat-flux proxy did not measure the capped face flux actually advanced by an RKL2 LF sweep | `src/diffusion/cgl_landau_fluid.cpp`, `src/mhd/mhd_sts.cpp`, `src/outputs/history.cpp`, `scripts/analyze_cgl_lf_paper.py` | High | Accumulate owned-face capped-flux/temperature-jump contractions through the same RKL2 recurrence as the state, checkpoint cumulative `lf_qprwrk`/`lf_qpewrk`, and count coarse/fine interfaces from fine-side closure fluxes that feed flux correction | Decay/sign, cap-analysis, low-field, AMR-nonzero, restart-preservation, and AMR MPI-decomposition regressions in the CGL test suites | Implemented for fixed-level and conserved-prolongation AMR paths; focused CPU AMR and fresh `amr` workflow/analyzer report nonzero applied work, and scheduler-launched Frontier Kokkos-Serial/MPI jobs `4659053`/`4659141` pass the AMR decomposition comparison at `rtol = atol = 1.0e-12` |
| F-020 | 2026-05-25 | Frontier operations | The generated environment-log filter required `MPICH=` instead of matching `MPICH_*`, so the first live manifests did not preserve the requested MPI/GPU environment variables | `scripts/frontier/cgl_lf_frontier.py`, `g001-g002-active-smoke-8gpu/manifest/run_environment.txt`, `g004-collision-relaxation-1gpu/manifest/run_environment.txt` | High | Capture configured runtime variable prefixes and archive optional restart files in the preparation utility before continuing restart qualification | Expanded `cgl_lf_frontier.py self-test` for environment capture and restart provenance; live G-004 inspection | Implemented; G-004 manifest retains requested MPICH/FI/HSA settings |
| F-021 | 2026-05-25 | Frontier G-001/G-007 | Live Slurm logs reported `MPICH_GPU_SUPPORT_ENABLED = 1` but `MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED = 0` although the generated script exported the latter as `1` | `/lustre/orion/ast207/proj-shared/dfielding/CGL/logs/slurm/cgl_g001_g002_active_smoke_8gpu.4653415.log`, `/lustre/orion/ast207/proj-shared/dfielding/CGL/logs/slurm/cgl_g006b_active_smoke_16blocks_16gpu_2node.4653601.log`, `/lustre/orion/ast207/proj-shared/dfielding/CGL/logs/slurm/cgl_g007b_paper_smoke_active_random_1gpu.4653770.log` | Medium | Treat runtime-reported behavior as authoritative, inspect the active stack and application memory space, and stop requesting managed-MPI support when it is not part of the model's memory contract | Installed Cray MPICH 9.0.1 `intro_mpi` documentation; AthenaK/Kokkos memory-space trace; job `4659663` and `analysis/f021_managed_policy_comparison.json` | Implemented and confirmed live: the AMD/HIP launch path requests managed-MPI support `0`; job `4659663` reports the requested values at runtime and is identical to job `4658072` across retained histories and state/forcing tabs |
| F-022 | 2026-05-25 | Frontier G-005 | The first GPU restart continuation restored physical fields and signed applied LF work but reset cumulative LF stage/face counters, making restart-spanning cap/exposure analysis discontinuous | `g005b-restart-resumed-1gpu/analysis/g005_restart_comparison.json`, `src/outputs/restart.cpp`, `src/pgen/pgen.cpp` | High | Checkpoint and restore all cumulative LF history diagnostics using the established shared-MPI rank-zero baseline rule, then rerun restart qualification after subsequent schema additions | Expanded CPU restart assertions for ordinary, finite-collision, and paper-forcing restarts; repeated Frontier G-005 | Implemented; corrected Frontier jobs `4653516`/`4653582` pass the original 15 LF endpoint comparisons, and post-F-026 jobs `4658072`/`4658163` pass restart comparison for added `lf_cpwrk`/`lf_cawrk` with differences below `1.2e-22` at `rtol = atol = 1.0e-12` |
| F-023 | 2026-05-25 | Phase E | The analyzer retained anisotropic-pressure transfer and heat-flux products but did not expose the configuration-space CGL pressure mechanical-work split required to interpret local anisotropic stress | `scripts/analyze_cgl_lf_paper.py`, `scripts/plot_cgl_lf_paper.py`, `20260525-paper-convergence-face-work-v1` analysis bundle | High | Add snapshot-derived `integral[p_perp div(u) - Delta p (b b : grad(u))] dV`, compare its anisotropic term with direct MKS24 transfer, label passive diagnostics without asserting applied feedback, and retain explicitly cadence-limited snapshot quadrature estimates | Synthetic correlated-pressure/strain transfer and constant-power quadrature checks through `test_cgl_lf_paper_snapshot_analysis_uses_both_pressures`; regenerated retained convergence bundle and pressure-work figure | Implemented for snapshot analysis; applied operator ledger added under F-026, while reference comparison and production qualification remain open |
| F-024 | 2026-05-25 | Phase E | The staged official source provides rendered panel PDFs rather than numeric curves, and the analyzer lacked a fail-closed contract for separately obtained or digitized reference data | `scripts/analyze_cgl_lf_paper.py`, `scripts/plot_cgl_lf_paper.py`, `scripts/cgl_lf_workflow.py` | High | Add optional reference-curve manifest ingestion requiring provenance, source/data SHA-256 values, and positive pointwise uncertainties; compare mapped products and render residual-qualified overlays | Generated exact-curve and wrong-source-checksum regressions in `test_cgl_lf_paper_snapshot_analysis_uses_both_pressures`; retained-bundle workflow rendering probe | Implemented for comparison infrastructure; populated MKS24 curve data and production comparison remain open |
| F-025 | 2026-05-25 | Phase E | `force_work` directly summed each explicit-RK source-call energy change even though the final RK state reweights earlier stage source increments; retained multi-cycle reduced bundles consequently reported false active global residuals while passing | `src/srcterms/turb_driver.cpp`, `scripts/cgl_lf_workflow.py`, `20260525-paper-convergence-face-work-v1` manifest | Blocker | Advance cumulative forcing work as a companion RK state using the cycle baseline and `gam0`/`gam1`, and reject active reduced cases whose energy/work relative residual is at least `1.0e-8` | `test_cgl_lf_paper_multicycle_forcing_work_follows_rk_state_recurrence`; regenerated `20260525-paper-smoke-rk-work-v1` and `20260525-paper-convergence-rk-work-v1` bundles | Implemented locally; focused CPU suite passes (`27 passed`), all six corrected convergence cases pass with maximum relative residual `5.141362829086986e-10` |
| F-026 | 2026-05-25 | Phase E | Snapshot-reconstructed pressure work could not identify the pressure traction actually applied by hyperbolic CGL momentum fluxes, especially after FOFC replacement and AMR flux correction | `src/mhd/rsolvers/hlle_cgl.hpp`, `src/mhd/mhd_tasks.cpp`, `src/mhd/mhd_fofc.cpp`, `src/outputs/history.cpp`, `scripts/analyze_cgl_lf_paper.py` | Blocker | Retain total and anisotropic CGL pressure traction in the numerical flux, apply the same AMR flux correction as momentum, contract corrected divergence with stage velocity, and advance `lf_cpwrk`/`lf_cawrk` by the explicit-RK companion recurrence; retain zero applied work for passive-Delta | Focused CPU AMR/restart/passive checks, dedicated CGL FOFC traction check, expanded MPI AMR decomposition oracle, regenerated local pressure-work bundles, and retained Frontier post-F-026 GPU/CPU evidence | Implemented for reduced qualification: local focused CGL CPU (`27 passed`) and style (`2 passed`) succeed; Frontier GPU restart jobs `4658072`/`4658163` and GPU AMR jobs `4658191`/`4658283` pass at `rtol = atol = 1.0e-12`; MPI CPU AMR jobs `4659053`/`4659141` also pass; production interpretation remains open |
| F-027 | 2026-05-25 | Frontier MPI CPU rerun | The debug preparation utility unconditionally generated GPU-aware MPICH exports and GPU binding; launching the intentionally non-GTL-linked Kokkos-Serial MPI executable consequently aborted before simulation with `MPIDI_CRAY_init: GPU_SUPPORT_ENABLED is requested, but GTL library is not linked` | `scripts/frontier/cgl_lf_frontier.py`, job `4658519` log and manifest | High | Add explicit `--execution-target gpu/cpu` generation; retain the existing GPU-aware path by default, while the CPU path omits GPU Slurm/srun binding and forces `MPICH_GPU_SUPPORT_ENABLED=0` | Expanded `cgl_lf_frontier.py self-test` and repeated scheduler-launched MPI CPU AMR comparison | Implemented and confirmed live: corrected jobs `4659053`/`4659141` run with `MPICH_GPU_SUPPORT_ENABLED=0`, refine cleanly, retain nonzero pressure work, and pass the one-rank/four-rank AMR oracle at `rtol = atol = 1.0e-12` |
| F-028 | 2026-05-25 | Frontier G-008 sizing | Original reduced sizing evidence retained output sizes but no MPI-I/O phase timing because runtime reported `MPICH_MPIIO_TIMERS = 0` | `scripts/frontier/cgl_lf_frontier.py`, `g008_reduced_sizing_evidence.json`, jobs `4660445`/`4661434` and their retained analysis JSON | Medium | Add manifest-recorded opt-in `--mpiio-timers`, then execute reduced and standard-layout startup shared-binary/checkpoint sizing runs; treat startup timing only as storage/I/O reconnaissance | Utility self-test; jobs `4660445`/`4661434`; retained G010b/G011 analysis JSON | Implemented for storage/I/O sizing: timers captured, strict diagnostics remain clean, and G011 supplies a preliminary standard-layout storage envelope; representative runtime/node-hour costing remains open |
| F-029 | 2026-05-25 | Production input audit after G011 | Defined standard/limiter decks relied on the default `single_file_per_rank = false`, while workflow manifests did not archive the output cadence/file-layout contract used for storage sizing | `inputs/cgl_lf_paper/`, `scripts/cgl_lf_workflow.py`, `analysis/g011_standard_layout_mpiio_sizing_evidence.json` | Medium | State shared binary/restart MPI-I/O explicitly in each guarded production deck and archive output layout/cadence choices in future workflow manifests | `test_cgl_lf_paper_production_inputs_explicitly_use_shared_mpiio`; metadata probe over `paper-standard`/`paper-nulim` cases | Implemented for defined inputs and future manifests; any later layout/cadence revision requires renewed storage and I/O review |
| F-030 | 2026-05-25 | Phase E reference-panel audit | MKS24 Figure 2(b) is an unstable-volume time history, but the checksum-qualified reference comparator admitted only snapshot-derived PDF/spectrum/transfer/alignment products | `scripts/analyze_cgl_lf_paper.py`, staged `MKS24.tex:508-512`, `docs/source/modules/cgl_landau_fluid_validation.md` | High | Export ordered threshold-volume history series, collapsing Athena's duplicate terminal row by retaining its last equal-time value, and admit them through `history.*` reference products | Expanded exact-curve regression in `test_cgl_lf_paper_snapshot_analysis_uses_both_pressures` for `history.unstable_fraction` | Implemented for comparison infrastructure; populated digitized MKS24 histories and target-run comparisons remain open |
| F-031 | 2026-05-25 | Frontier G012 nonlinear hard-wall validation | The intended active beta10 hard-wall run with uncapped RKL2 STS reaches the parallel firehose limiter and then crosses the strict emergency bound in an LF split stage before its analysis window | G012 job `4662477` log, retained histories/snapshots, and `analysis/g012_hardwall_strict_failure_evidence.json`; `analysis/f031_local_screen/f031_local_timestep_screen.json` under the G012 run; `src/mhd/mhd_sts.cpp` | Blocker | Retain strict checking, reject tested RKL2 caps, diagnose an explicit-LF or revised splitting/integration policy for limiter-active nonlinear runs, and rerun reduced GPU qualification before any runtime costing or production proposal | Local reduced screen retains strict aborts for uncapped RKL2 and capped ratios `20`, `10`, `5`, `2`, and `1`; explicit LF is clean through `t = 1.6` only when branched from a clean `t = 1` checkpoint, while switching from the clean ratio-2 `t = 1.5` checkpoint also fails near `t = 1.584` | Addressed locally by F-033 constrained-model validation; G012 consumed `0.188889` node-hours, and corrected GPU qualification remains required before runtime projection |
| F-032 | 2026-05-25 | Workflow provenance after F-031 | Workflow manifests archived closure and output layout choices but omitted the STS timestep cap and integrator controls that determine limiter-active qualification | `scripts/cgl_lf_workflow.py`, `inputs/cgl_lf_paper/`, F-031 | High | Archive time integrator, STS integrator, `sts_max_dt_ratio`, and CFL number in future workflow manifests together with the F-033 hard-wall choice | Extended `test_cgl_lf_paper_production_inputs_explicitly_use_shared_mpiio` metadata assertions | Implemented for provenance; the F-033 model/timestep contract remains pending reduced Frontier GPU qualification |
| F-033 | 2026-05-25 | F-031 model correction/local same-state validation | MKS24's hard-wall `nu_lim = 1.0e10` intent is not robustly represented by finite-rate post-LF scattering: the screened failing state violates the firehose emergency bound during internal RKL2 `post stage=1/5` | `src/eos/cgl_physics.hpp`, `src/eos/cgl_mhd.cpp`, `src/outputs/history.cpp`, and `analysis/f033_hardwall_projection_local/f033_hardwall_projection_local_evidence.json` beneath G012 (SHA-256 `29ce073e1059ac38496b0f543716b58b4cc5a3ff7a37c4457cb9217affca3eab`) | Blocker | Add opt-in `limiter_hardwall` energy-preserving projection at CGL primitive recovery, replace finite-rate limiter pressure relaxation only in hard-wall mode, append persistent `lf_hwproj`, select it in hard-wall paper decks, and rerun corrected GPU qualification | `test_cgl_lf_hardwall_projects_to_selected_firehose_threshold`, `test_cgl_lf_hardwall_requires_instability_limiter`, expanded restart/input metadata assertions; focused CPU module `30 passed`; style and Sphinx gates passed | Implemented and passed locally through exact-state `t = 1.5` to `t = 2.0` with zero strict safety counters and final `lf_hwproj = 53829433`; reduced Frontier GPU requalification remains open |

### Implemented Core Decision Log

| Date | Decision | Rationale | Evidence still required |
| --- | --- | --- | --- |
| 2026-05-24 | Keep the existing LF source placement and apply collisions after each LF half-sweep with explicit duration `dt/2`. | This repairs physical elapsed time without changing the protected LF representation/task ownership; both collision calls together span one `dt`. | Finite-collision restart and `compare` workflow passed on CPU; broader MPI/GPU gates remain |
| 2026-05-24 | Preserve `cgl_firehose_threshold = oblique` as the public default and add explicit `parallel` selection for MKS24. | Existing branch inputs retain their historical convention; paper inputs cannot accidentally depend on that default. | Paper pgen/decks and paper manifests |
| 2026-05-24 | Define firehose emergency overshoot at `p_perp - p_parallel <= -1.5 B^2` (`beta Delta <= -3`) and count it independently of `backup_limiters`. | The emergency condition must remain distinct from normal finite-rate activity when the physical MKS24 activation threshold is `beta Delta <= -2`. | Reduced nonlinear validation of overshoot interpretation |
| 2026-05-24 | Reuse and repair the built-in turbulence driver for active paper smoke rather than copy donor forcing logic into the pgen. | Shared forcing now records a deterministic seed, has explicit `B0 || z` Alfvenic orientation, preserves restart state, advances OU state once per cycle, and supplies RK-consistent energy work. | Calibrate published injection/spectral normalization and qualify larger runs |
| 2026-05-24 | Qualify reduced `passive_delta` only as isothermal-MHD flow with diagnostic CGL/LF thermodynamics, requiring matching `mhd/passive = true`. | Passive fluxes and CFL speeds are independent of diagnostic anisotropy; a paired-anisotropy regression shows identical density, velocity, and magnetic evolution under identical forcing. | Long-time active/passive statistics, execution of defined standard decks, and spatial paper observables |
| 2026-05-24 | Stage official arXiv source version `2405.02418v2` into ignored result storage with SHA-256 provenance rather than vendoring paper contents. | A pinned source archive supports exact reference interpretation while leaving third-party redistribution outside the code branch. | Analysis bundles now retain the pinned checksum provenance; implement quantitative reference curves/comparisons |
| 2026-05-24 | Add explicit physical-shell/common-spectrum turbulence settings for MKS24 inputs while preserving legacy driver defaults. | The paper specifies physical `abs(k)` bounds and `k^-2` power in its anisotropic domain; index shells and type-specific default exponents do not encode that setup. | Measure long-time injected spectra and steady-state statistics in standard configurations |
| 2026-05-24 | Extend CGL primitive outputs with `p_perp` while preserving legacy `eint = p_parallel`. | Offline anisotropy PDFs, spectra, and transfer diagnostics require both pressures; retaining the legacy label avoids breaking existing readers. | Focused CPU and binary snapshot analyzer verification passed; standard-run statistics and comparison remain open |
| 2026-05-24 | Measure heat-flux-cap occupancy from the LF operator-face closure rather than duplicating it in the offline analyzer. | The applied cap and limiter-selected diffusivity already exist in shared C++ physics; cumulative face counters preserve that exact decision and can be differenced over analysis windows. | Strong-cap oracle and windowed reduced smoke pass; signed applied face contraction with fine-side AMR ownership subsequently implemented under F-019; production statistics remain open |
| 2026-05-25 | Treat `|B| <= bfloor` faces as having no usable LF direction and verify zero LF face transport in a strict kinematic test. | The closure cannot define `bhat` below its configured floor; explicit shutdown is finite and auditable rather than directionally arbitrary. | Focused CPU low-field regression passed; nonlinear near-floor cases remain outside the paper baseline |
| 2026-05-25 | Render generic archived paper diagnostics separately from paper-panel comparison claims. | Inspectable local products are needed before production, but without reference data or standard runs they cannot establish MKS24 figure reproduction. | Binary probe and smoke rendering passed; reference comparison and production statistics remain open |
| 2026-05-25 | Define `paper-convergence` as a short reduced nonlinear qualification tier rather than a substitute for paper-scale convergence. | It exercises resolution, heat-flux, threshold, snapshot, and plotting paths locally while its `tlim = 0.02` metadata prevents steady-state interpretation. | Corrected `20260525-paper-convergence-rk-work-v1` passed and rendered six generic figures; longer reduced/standard statistics remain open |
| 2026-05-25 | Keep strict emergency monitoring enabled in limiter-suppression validation and move its physical limiter state below the emergency bound. | A limiter-action oracle must not intentionally begin in a state the independent safety contract defines as invalid. | Corrected case retains strong flux suppression and `full-v2` passes |
| 2026-05-25 | Keep Frontier submission manual while automating fail-closed preparation and accounting policy checks. | The debug-QOS one-job constraint forbids job chaining; generating an explicit script and retained reservation/manifest improves reproducibility without silently submitting work. | Utility self-test passes; immutable HIP/MPI and Kokkos-Serial/MPI builds plus sequential debug accounting through job `4661434` retain `0.085003` node-hours recorded |
| 2026-05-25 | Request GPU-aware MPI without requesting managed-memory MPI support for the AMD/HIP CGL-LF baseline. | Installed Cray MPICH 9.0.1 documents managed-memory MPI default `0` for non-NVIDIA GPUs, and AthenaK's `DvceArray*` storage resolves through Kokkos HIP to `HIPSpace`, not `HIPManagedSpace`; `HSA_XNACK=1` remains a separately recorded GPU-runtime configuration. | Job `4659663` records generated and reported `MPICH_GPU_SUPPORT_ENABLED=1` with `MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0`; its histories and state/forcing tabs are identical to job `4658072` |
| 2026-05-25 | Use an isolated dependency-complete Python environment for local style and Sphinx gate execution. | The checkout declares Sphinx dependencies but the active system module lacks them and `flake8`; installing validation-only packages under `/tmp` avoids changing the product dependency surface. | Repository style suite and warning-strict Sphinx build passed; no generated documentation output is retained as source |
| 2026-05-25 | Emit heat-flux transport strength initially as a snapshot-reconstructed proxy rather than claim an applied discrete energy budget. | Archived fields and closure choices are sufficient to reconstruct the continuous LF smoothing measure; exact finite-volume accounting must first define a decomposition-independent applied-face reduction. | Synthetic sign check passed; F-019 later adds the signed applied contraction with fine-side AMR ownership; production/reference interpretation and full budget closure remain open |
| 2026-05-25 | Track cumulative net applied forcing work at the source update only when explicitly requested by a driven input. | The before/after conserved-energy difference across the full source task, including zero-net-momentum projection, is exact at that stage and restartable; instantaneous output-time `force_pwr` cannot close an active energy residual. Opt-in persistence avoids altering legacy forcing restart records. | Focused CPU/restart coverage and corrected `20260525-paper-smoke-rk-work-v1`/`20260525-paper-convergence-rk-work-v1` bundles; production tolerance and heat-flux transfer decomposition remain open |
| 2026-05-25 | Retain applied LF heat-flux transport as a signed RKL2 companion recurrence with fine-side ownership at coarse/fine interfaces, not as a total-energy closure term. | The applied capped face flux and STS coefficients are available in the LF update; on AMR interfaces the fine-side closure flux is the quantity restricted into the coarse update, so it owns the contraction without double counting. | Focused CPU/restart tests and fresh fixed-level bundles pass; fresh AMR workflow/analyzer report total work `6.511527338026022e-4`; Frontier Kokkos-Serial/MPI jobs `4659053`/`4659141` pass the AMR decomposition regression |
| 2026-05-25 | Checkpoint every cumulative LF diagnostic history value, not only signed applied work. | Restart-spanning paper analysis differentiates face/cap/exposure counters; resetting their baseline after a valid restart silently corrupts intervals even when evolved fields agree. | Initial GPU G-005 exposed F-022; expanded CPU restart checks and corrected jobs `4653516`/`4653582` verify the correction |
| 2026-05-25 | Treat G-008/G010b/G011 debug measurements as storage/I/O reconnaissance, not a paper-scale runtime cost model. | Startup-length runs expose elapsed time, memory, output volume, and checkpoint size; G011 adds exact standard-layout shared-file timing and size, but no retained run measures late-time timestep evolution or statistical-production duration. | G010b retains reduced timing and its excessive cadence warning; G011 job `4661434` retains standard-layout binary/restart net-write means of `1772.938`/`1252.430` MiB/s and the `328.915` GB eight-case raw storage envelope; runtime/node-hour costing remains open |
| 2026-05-25 | Stop runtime costing after G012 and treat limiter-active RKL2 STS as unqualified for the screened finite-rate nonlinear state. | The finite-rate hard-wall reduced nonlinear run and local reduced branches with unlimited STS or caps `20`, `10`, `5`, `2`, and `1` all cross the strict firehose emergency bound in an LF split stage; reducing the cap through the documented conservative setting delays but does not resolve the observed failure. | F-031 is addressed locally by the F-033 constrained model; job `4662477` consumes `0.188889` node-hours for a cumulative `0.273892`, and corrected GPU requalification remains required before further costing |
| 2026-05-25 | Archive the STS timestep/integrator contract in paper workflow manifests. | F-031 makes the split-step cap a qualification-critical model choice; a result cannot be interpreted or repeated if the manifest omits it. | F-032 metadata regression records guarded-deck `sts_max_dt_ratio = -1.0`; retain it with the F-033 hard-wall model and GPU-qualify that combined contract |
| 2026-05-25 | Represent the MKS24 hard-wall limit explicitly with `limiter_hardwall = true` in hard-wall paper decks rather than finite-rate pressure relaxation at `limiter_nu_coll = 1.0e10`. | The pinned paper source describes anisotropy effectively pinned at thresholds; F-031 locates a same-state internal LF-stage overshoot under finite-rate scattering, while the energy-preserving constrained state passes locally through `t = 2.0`. Finite-rate `20`/`200` decks remain unchanged. | F-033 CPU and exact-state local validation pass; corrected reduced Frontier GPU execution and subsequent timing/proposal qualification remain required |
| 2026-05-25 | Implement CGL pressure mechanical work first as a snapshot-derived local decomposition with a transfer cross-check and sparse-time quadrature estimate. | Retained fields determine `p_perp div(u) - Delta p (b b : grad(u))` without changing the integrator; trapezoidal snapshot-time integration aids interval interpretation, but neither result is an applied budget identity. | F-023 synthetic sign/transfer/quadrature oracles and corrected `20260525-paper-convergence-rk-work-v1` products; exact production interval closure remains open |
| 2026-05-25 | Admit paper reference curves only through an external checksum- and uncertainty-qualified manifest. | The staged arXiv source supplies panel PDFs but not numeric curve tables; separating source data from code and rejecting absent uncertainties permits later author-data or documented-digitization comparisons without overstating current evidence. | F-024 exact-curve and source-checksum-rejection regressions plus retained reduced-bundle rendering probe; actual MKS24 curve acquisition/digitization and production comparisons remain open |
| 2026-05-25 | Expose threshold-volume histories as reference-comparable curves while collapsing duplicate terminal timestamps. | Figure 2(b) is a time-series result rather than a snapshot statistic, and Athena repeats an unchanged terminal history row; retaining the final equal-time value gives ordered coordinates without changing the measured endpoint. | F-030 synthetic exact-history comparison passes; actual digitized history curves and standard-run comparison remain open |
| 2026-05-25 | Advance `force_work` through the same explicit-RK recurrence as the source-modified conserved state. | Adding every source-call energy increment overcounts work retained in a multi-stage final state because later RK combinations reweight earlier source contributions; active reduced cases must fail when the ledger does not close conserved energy. | F-025 multi-cycle CPU regression; corrected smoke and convergence bundles pass with maximum active relative residual `5.141362829086986e-10`; production residual qualification remains open |
| 2026-05-25 | Define applied CGL pressure work from the isolated numerical pressure traction after momentum AMR flux correction, then integrate it with the explicit-RK companion recurrence. | A snapshot-gradient estimate is not the applied finite-volume operator and a face-jump contraction is ambiguous across coarse/fine interfaces; contracting each active cell's corrected pressure-traction momentum RHS with its stage velocity uses exactly the operator consumed by the update. | F-026 focused AMR/restart/passive/FOFC regressions and fresh pressure-work bundles pass locally; Frontier GPU restart and GPU/CPU AMR decomposition comparisons pass for `lf_cpwrk`/`lf_cawrk`; production interpretation remains open |
| 2026-05-25 | Generate distinct Frontier debug launch environments for GPU-aware qualification and MPI CPU decomposition validation. | A Kokkos-Serial MPI executable is not linked against Cray's GPU transport library and must not inherit `MPICH_GPU_SUPPORT_ENABLED=1` or GPU binding from the HIP/MPI launch path; retaining the target in each manifest avoids ambiguous evidence. | F-027 failure job `4658519` identifies the pre-launch abort; utility self-test verifies CPU and GPU generated scripts separately; corrected CPU-target jobs `4659053`/`4659141` pass |

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
| `paper-smoke` | Reduced active/passive-Delta CGL-LF turbulence with Alfvenic/random forcing, retained forcing metadata and tables |
| `full` | Eighteen operator/linear-wave/eigenmode cases plus plots |
| `plot` | Regenerate existing operator-validation plots |
| `summarize` | Regenerate a result summary from retained outputs |

This workflow layer must remain usable throughout the reproduction work.
It should not silently change meaning. In particular, `full` is presently a
closure/operator validation tier, and `paper-smoke` is only a reduced
active/passive turbulence gate; neither is an MKS24 production reproduction tier.

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

### Defect 1: Collision Relaxation Advanced Two Full Timesteps per LF Cycle

**Implementation status (2026-05-24):** Corrected in the working tree by
passing the explicit LF half-sweep duration to each post-sweep collision
update. The independent uniform-anisotropy CPU regression passes; broader
phase gates remain required.

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

### Defect 2: Baseline Firehose Activation Did Not Match MKS24

**Implementation status (2026-05-24):** Corrected in the working tree with
the `cgl_firehose_threshold = oblique|parallel` policy. Focused policy tests
pass; the reduced active and passive paper-smoke decks and the defined
paper-standard/limiter-scan inputs explicitly select `parallel`. Standard
paper execution and analysis remain open.

#### Baseline Behavior

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

### Defect 3: Hard-Bound Monitoring Depended on Corrective Backup Mode

**Implementation status (2026-05-24):** Corrected in the working tree.
`lf_hardbd` now reports emergency-bound crossings independently of
`backup_limiters`, with a focused strict negative test passing.

#### Baseline Behavior

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

**Implementation status (2026-05-24):** Partially advanced. A new
`paper-smoke` workflow runs reduced active-Delta Alfvenic/random and
passive-Delta Alfvenic forcing cases, archives forcing/model choices, retains
state and force tables, and checks clean LF safety and positive RK-integrated
forcing work. Focused CPU regressions verify one OU update per physical cycle,
the multi-cycle RK2 companion-work identity, and passive flow independence from diagnostic
initial anisotropy. It does not supply long-time paper-grade active/passive
statistics, standard-resolution results, or figure/reference comparisons. Workflow
manifests now expose dirty worktree status so development evidence cannot be
mistaken for an immutable production-revision bundle.

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

### Gap 5: Analysis Infrastructure Exists but Paper Comparisons Are Missing

The paper pgen now provides reduced volume-integrated histories for energy,
`C_B2`, pressure anisotropy, beta, physical threshold volumes, effective
collision frequency, instantaneous forcing power, and, for paper decks,
restartable cumulative RK-integrated applied forcing work. The analyzer
compares active-case cumulative applied work with conserved `tot-E` over the
selected history interval; active reduced workflows reject a relative residual
of `1.0e-8` or larger, while passive-Delta retains work without receiving this
active energy-budget interpretation. The corrected
`20260525-paper-convergence-rk-work-v1` bundle passes all six cases with a
maximum relative active residual of `5.141362829086986e-10`. Current binary output supplies both CGL
pressures, and `scripts/analyze_cgl_lf_paper.py` now produces snapshot
PDFs, shell spectra, field-projected anisotropy gradients, pressure-transfer
partitions, alignment PDFs, and per-case manifest-window ensemble averages
with synthetic checks. It also produces local-field velocity-gradient and
`b b : grad u` products, while `scripts/plot_cgl_lf_paper.py` renders generic
archived diagnostics. Snapshot products now also split the CGL pressure
mechanical power as
`integral[p_perp div(u) - Delta p (b b : grad(u))] dV`; the analyzer
compares its anisotropic term with the direct pressure-transfer integral and
labels passive-Delta values as non-feedback diagnostics. Multi-snapshot
ensembles additionally retain explicitly cadence-limited trapezoidal
pressure-work and reconstructed heat-flux integral estimates. Cumulative
`lf_q*` history counters now retain exact
operator-face heat-flux-cap activity for interval summaries. For new bundle
manifests that archive complete LF closure policy choices, the analyzer also
emits a cell-centered reconstruction of
`integral[-q_parallel b.grad(T_parallel) - q_perp b.grad(T_perp)] dV`,
separating regularized, unlimited, and cap-active contributions. The JSON
identifies it as a snapshot proxy, not an applied face-flux reduction or
energy-budget identity. Runs additionally retain restartable `lf_qprwrk` and
`lf_qpewrk` histories: signed RKL2-applied owned-face contractions of capped
parallel/perpendicular heat flux with temperature jumps. On AMR interfaces,
fine-side closure faces own the contractions because those fluxes feed coarse
flux correction. This quantity is not a total-energy or time-integrated
local-work closure. Enabled paper and AMR decks now additionally retain
`lf_cpwrk` and `lf_cawrk`, the total and anisotropic applied hyperbolic CGL
pressure-traction work. These histories isolate the numerical pressure
traction in HLLE or FOFC, apply momentum's AMR flux correction, contract its
cell divergence with stage velocity, and advance it through the explicit-RK
companion recurrence. Passive-Delta records zero applied pressure work.
The fresh `20260525-paper-smoke-pressure-work-v1` and
`20260525-paper-convergence-pressure-work-v1` bundles pass; the fresh AMR
bundle retains nonzero values (`lf_cpwrk = 8.080431744293371e-08`,
`lf_cawrk = 9.674389094505345e-08`). Other `lf_*` counters remain
safety/stage counters. The
analyzer also accepts optional external reference-curve manifests that
require checksums and pointwise uncertainty values before calculating
comparison residuals or producing overlay figures.
Figure-level reproduction still lacks:

- populated panel-to-product reference manifests for standard cases;
- quantitative comparison to machine-readable or explicitly digitized MKS24
  curves; the staged official source contains PDFs but no numeric curve
  tables;
- production interpretation of the applied pressure/heat-flux work ledgers;
- production qualification of active global budget residuals and the
  snapshot/applied local-work products.

These products must be implemented before figure-level reproduction is
possible.

## MKS24 Reproduction Contract

### What Paper Is Being Reproduced

The pinned reference is official arXiv source version `2405.02418v2`. Stage
it rather than vendor it:

```bash
python3 scripts/stage_cgl_lf_mks24_reference.py
```

A validated local staging run produced:

```text
build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/manifest.json
build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/source/MKS24.tex
archive SHA-256: 6b58e1ed01585c9820659a01b7549a08a09b7ea3fcbb28dbb96205259e82da8e
MKS24.tex SHA-256: 6d6e748fd1883c5d33167be653d67aed2f84a9b364267d9e90ea075df184af4c
```

The paper is:

```text
Self-organization in collisionless, high-beta turbulence
S. Majeski, M. W. Kunz, and J. Squire
```

Future agents should read the staged source directly before modifying
reproduction inputs or interpreting plots. At minimum, reread:

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
| Figure 2: pressure-density distributions and unstable-volume history | Active and passive, Alfvenic, representative beta values | PDFs of `delta p_parallel`, `delta p_perp`, density; volume fraction beyond thresholds versus time | Reduced threshold-volume and steady-window PDF machinery plus `history.unstable_fraction` reference comparison implemented; populated curve data and standard statistics missing |
| Figures 3-4: compressive/density behavior | Active/passive; beta and forcing variants | Density PDFs and density/compressive spectra | Steady-window density PDF/spectrum machinery implemented; target runs/comparison missing |
| Figures 5-6: turbulent spectra and rate-of-strain suppression | Active/passive, Alfvenic/random | Velocity/magnetic/pressure spectra; parallel/perpendicular velocity-gradient spectra | Local-field velocity-gradient and strain products plus generic rendering implemented; target runs/comparison missing |
| Figure 7: pressure-anisotropy organization and transfer | Active/passive, Alfvenic/random, beta 10 | `grad_parallel Delta p`, `grad_perp Delta p` spectra; `T_Delta p(k_perp)/T_total` | Steady-window gradient spectra, transfer partition/closure, and snapshot anisotropic-work cross-check implemented; reference comparison missing |
| Figures 8-9: alignment diagnostic | Active/passive; beta/forcing variants | Scale-dependent PDFs and peak locations of local-field/rate-of-strain eigenvector alignment | Steady-window selected-shell alignment PDFs, rendered peak summaries, and synthetic finite check implemented; comparison missing |
| Figure 10: kinetic comparison context | Not reproduced directly by CGL-LF | Documentation of non-reproducible kinetic comparator and any qualitative CGL-only comparison | Out of direct scope |
| Figure 11: scale-separation behavior | Selected resolutions | Resolution/convergence or scale-separation comparison products | Operator accuracy and short reduced resolution product paths implemented; steady nonlinear comparison missing |
| Figure 12: heat-flux sensitivity | Heat-flux parameter variants | Spectra/transfer/diagnostic products under heat-flux changes | Short reduced parameter path, reconstructed smoothing proxy, and signed applied face contraction with fine-side AMR ownership implemented; steady nonlinear comparison and exact interval/production budget missing |
| Figure 13: limiter-frequency sensitivity | Beta 100, Alfvenic, several `nu_lim` | Energy spectra, `beta Delta` PDFs, strain PDFs, pressure-stress transfer | Inputs, threshold policy, beta-Delta PDF and transfer foundations implemented; case runs/full products missing |

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

**Implementation status (2026-05-25):** The local `accuracy` workflow now
archives parallel/perpendicular collisionless decay at three resolutions and
three STS timestep ratios plus explicit references, both LF coefficient modes
at collision rates below, comparable to, and above `c_parallel |k_parallel|`,
and STS/explicit extreme-cap cases. The retained
`20260524-accuracy-v2` bundle passed with decreasing resolution errors, a
maximum finite-collision relative error below `3.0e-4`, and parallel and
perpendicular pre-cap ratios of `3.62e5` in both extreme-cap methods.
Field-aligned decay in the `x`, `y`, `z`, and periodic oblique directions
passed with relative errors no greater than `8.8e-4`; a strict low-field
case confirmed zero LF face transport below `bfloor`. Production forcing-energy
residual qualification remains open. The final broad local CPU gate passed with
`218 passed, 15 skipped, 59 deselected`.

Direct invocation of the required MPI CPU suite remains unavailable in the
current shell because it does not provide `mpirun`. A scheduler-capable
Frontier Kokkos-Serial/MPI build from revision `3e3c206b` instead executed
the same one-rank/four-rank AMR oracle through `srun`: jobs `4659053` and
`4659141` pass at `rtol = atol = 1.0e-12`, refine cleanly, retain nonzero
`lf_cpwrk`/`lf_cawrk`, and retain zero strict failure counters. The missing
local launcher is therefore no longer an unmet MPI numerical-evidence gate.

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

**Implementation status (2026-05-24):** Active/passive-Delta `turbulence` mode is
registered in the default binary through `src/pgen/tests/cgl_lf_paper.cpp`.
The reduced beta-10 inputs and `paper-smoke` workflow pass for active
Alfvenic/random and passive Alfvenic forcing with clean strict diagnostics,
retained forcing metadata, and restart-tested deterministic forcing.
Once-per-cycle OU evolution, stage-consistent RK coupling, and reduced
passive flow decoupling pass focused CPU regressions. Standard decks,
including the beta scan and beta-100 limiter-frequency suite, are now defined
for guarded production workflows but have not been executed. Paper-grade
forcing calibration and controlled analogue modes are not yet qualified.

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

**Status (2026-05-25):** All listed smoke, standard baseline, and beta-100
limiter-frequency input definitions are present. The standard and limiter
decks now explicitly select shared binary/restart MPI-I/O and their future
workflow manifests archive the output cadence/layout choices used for sizing
and the timestep/STS/CFL contract exposed by F-031. They remain production
definitions only; their present unlimited STS setting is not qualified for
limiter-active nonlinear evidence, and the screened `sts_max_dt_ratio = 1.0`
setting or same-state explicit-LF switch is not a qualified replacement.

### Pgen-Level Acceptance Criteria

- A small active turbulence case initializes and evolves with clean strict
  safety diagnostics. **Passed for the reduced CPU smoke case on
  2026-05-24.**
- A small passive case evolves with the documented passive feedback behavior.
  **Passed for reduced isothermal-MHD flow decoupling from diagnostic
  anisotropy on CPU on 2026-05-24.**
- Restart does not corrupt the CGL representation or forcing setup.
  **Passed for the reduced active forcing case on 2026-05-24.**
- Once-per-cycle OU correlation and RK source coupling are correct under the
  selected time integrator. **Passed for the reduced active CPU case on
  2026-05-24; long-time calibration against published forcing statistics is
  still open.**
- Threshold and LF closure choices are visible from input and manifest output.
  **Passed for threshold, passive-mode, and forcing choices in the reduced
  `paper-smoke` manifest; standard workflow metadata support is implemented
  but awaits authorized execution.**
- No pgen-local duplicated threshold or LF formula can drift from shared
  physics helpers. **Satisfied for active initialization; it calls the shared
  CGL/LF modules and does not encode limiter thresholds or heat-flux
  formulas.**

### Recommended Commit

```text
Add current-API CGL-LF paper problem generator and reduced turbulence inputs
```

## Phase E: Add Paper Diagnostics and Analysis Products

### Objective

Produce the observables needed to make quantitative comparisons with MKS24,
while keeping raw simulation output manageable and archived analyses
reproducible.

**Implementation status (2026-05-25):** The paper pgen now emits reduced
volume-integrated histories for mass, energies, `B^2`, `B^4`, pressure
anisotropy, beta, threshold volumes, effective collision frequency, and
forcing power. `paper-smoke` summaries report beta, `C_B2`, and threshold
fractions; paper histories also retain restartable cumulative RK-integrated
applied forcing work, and active reduced workflows require the
conserved-energy residual to close below the defined relative tolerance.
CGL primitive snapshots now export both pressures, and
`scripts/analyze_cgl_lf_paper.py` generates per-snapshot PDFs, shell spectra,
projected pressure/velocity-gradient spectra, `b b : grad u` strain PDFs,
pressure-transfer partitions with a real-space closure check, selected-shell
alignment PDFs, and per-case ensemble averages inside the archived analysis
window. It additionally reports snapshot CGL pressure mechanical work as
separate `p_perp div(u)` and `-Delta p (b b : grad(u))` terms, compares the
anisotropic term with direct pressure transfer, and preserves whether an
archived case applies that feedback to the flow. For retained multi-snapshot
windows it records trapezoidal time-integral estimates of the pressure and
heat-flux proxy rates, explicitly excluding them from applied budget claims.
`scripts/plot_cgl_lf_paper.py` renders available archived diagnostic figures
and is invoked by `paper-analyze`. LF history now retains exact operator
face counts for parallel/perpendicular heat-flux-cap ratios above `1` and
`10`, which `paper-analyze` differences over that window. The analyzer now
also reconstructs a cell-centered heat-flux smoothing-power proxy from
snapshots and archived closure choices, reporting regularized, unlimited, and
cap-active contributions while explicitly excluding it from discrete
face-flux or budget claims. LF `.mhd.hst` output additionally retains
restartable signed RKL2-applied capped-face heat-flux contractions, with
coarse/fine interfaces owned by the fine-side closure faces; these
contractions need not equal that snapshot proxy. Its
synthetic numerical test and a two-snapshot window-selection/rendering probe
have passed locally. For decks with `cgl_lf_record_pressure_work = true`,
LF history also retains `lf_cpwrk` and `lf_cawrk`, explicit-RK-applied total
and anisotropic CGL pressure-traction work after the same AMR flux correction
as momentum. The analyzer can now validate external reference-curve
checksums and uncertainties and render comparison overlays once curve data
are supplied. Populated
paper-panel/reference comparisons,
long-time residual qualification, and production statistics remain open.

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

Implemented paper user-history columns are `volume`, `mass`, `kinetic`, `magnetic`,
`therm_cgl`, `b2`, `b4`, `delta_p`, `abs_dp`, `beta`, `mirror_vol`,
`fire_vol`, `hard_vol`, `nu_eff`, `force_pwr`, `force_work`, and retained
forcing-orientation quantities. `force_work` is the exact cumulative
RK-integrated net forcing source, including zero-net-momentum projection, for
paper inputs that enable its restartable driver counter; active-case analysis
differences it against conserved `tot-E`. Normal LF MHD
history additionally retains
`lf_qface`, `lf_qprcap`, `lf_qpr10`, `lf_qpecap`, and `lf_qpe10` for
operator-face cap occupancy, plus restartable `lf_qprwrk` and `lf_qpewrk`
for signed RKL2-applied capped-face transport contractions, using fine-side
closure-face ownership at coarse/fine interfaces. Inputs enabling
`cgl_lf_record_pressure_work` also retain `lf_cpwrk` and `lf_cawrk` for
explicit-RK-applied total and anisotropic hyperbolic CGL pressure traction;
these fields use the corrected momentum-flux divergence on AMR meshes and
are zero for passive-Delta. Offline snapshots
additionally provide a reconstructed heat-flux smoothing proxy and a CGL
pressure-work decomposition.

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

The implemented analyzer is current-format code rather than a direct donor
copy. It reads retained binary `mhd_w_bcc` products through the repository
binary reader and writes `diagnostics.json`; CGL snapshots use legacy
`eint = p_parallel` together with the added `p_perp` field.

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

**Status:** Reduced history JSON summaries and manifest-selected steady-window
averages are implemented through the paper analyzer and workflow integration.
The analyzer also reports interval heat-flux-cap fractions from exact
operator-face counters, and the paper plotter renders retained history
summaries when available. It reports active-case conserved-energy change
against restartable cumulative RK-integrated applied forcing work and gates
reduced cases on the residual, while excluding passive-Delta from that
total-energy interpretation. A closure-reconstructed
heat-flux smoothing proxy is reported from manifest-qualified snapshots; LF
history also reports RKL2-applied capped-face transport contractions with
fine-side ownership at coarse/fine interfaces. Snapshot analysis additionally
reports the CGL pressure
mechanical-work split, its direct-transfer comparison, and a sparse-snapshot
quadrature estimate where cadence allows. It now reports interval
`lf_cpwrk`/`lf_cawrk` applied pressure-traction work where enabled; production
residual qualification and paper comparison remain open.

#### E2. PDF Products

Compute:

- pressure-density fluctuation PDFs for Figure 2-style comparisons;
- `beta Delta` PDFs for limiter-frequency comparisons;
- selected strain-distribution PDFs for Figure 13-style interpretation.

Record binning, normalization, selected snapshot times, and averaging windows.

**Status:** Per-snapshot density, parallel/perpendicular pressure-fluctuation,
`beta Delta`, and `b b : grad u` PDFs plus common-bin steady-window averages
and generic plots are implemented. Reference-panel targeting remains open.

#### E3. Spectra

Compute bin-averaged spectra in the convention described by MKS24:

- velocity and magnetic fluctuations;
- density/compressive quantities;
- `Delta p`;
- field-parallel and perpendicular gradients of `Delta p`;
- field-parallel and perpendicular components of velocity gradients.

Implement and test isotropic and `k_perp` shell binning as required. Specify
normalization conventions in documentation and output metadata.

**Status:** `k_perp` shell products for velocity, magnetic fluctuation,
density, `Delta p`, projected `Delta p` gradients, and local-field
parallel/perpendicular velocity gradients are implemented with a synthetic
single-mode binning check and generic plot output. Reference normalization
qualification remains open.

#### E4. Pressure-Anisotropy Transfer Function

Implement the paper-defined `T_Delta p(k_perp)` diagnostic. This is one of the
highest-risk analysis tasks because implementation choices can alter sign and
normalization. Requirements:

- reproduce the formula from the paper in code comments or analysis docs;
- write synthetic-field tests with known zero or sign behavior;
- compare integral transfer against a real-space anisotropic work measure where
  possible;
- archive normalization conventions in summary metadata.

**Status:** A shell-partition implementation of the stated pressure-stress
transfer and its direct real-space closure error are emitted by the analyzer;
the snapshot-derived anisotropic-stress work term is emitted with a
transfer-minus-work comparison residual. Synthetic zero-transfer/closure and
correlated negative-work/zero-residual checks pass, including passive-only
interpretation metadata and constant-power time quadrature. Sparse-snapshot
quadrature is retained as an estimate, not applied budget closure.
Normalization must still be qualified against production fields/reference
interpretation.

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

**Status:** Selected-`k_perp` stretching-eigenvector alignment PDFs and
rendered selected-shell peak curves are implemented, and synthetic output is
checked for finite bounded histograms. Paper-panel aggregation and physically
interpretable production validation remain open.

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

**Implementation status (2026-05-25):** `paper-smoke` executes locally and
archives reduced model/forcing choices. `paper-standard` and `paper-nulim`
now enumerate the defined production inputs, retain model/resolution/domain/
analysis-window metadata and in-place binary/checkpoint references, and
refuse execution unless `--authorize-paper-execution` is passed.
`paper-analyze --output-dir <bundle>` regenerates current history/snapshot
diagnostics with a synthetic test and generic diagnostic figures;
`paper-summary` regenerates retained summaries without plotting. Local
operator convergence is available through `accuracy`; `paper-convergence` is
implemented as a short startup/product-path qualification tier. Steady
nonlinear convergence and reference-comparison modes remain open.

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
| `paper-smoke` | Reduced active/passive forcing-mechanics and limiter-policy startup checks | Workstation/local |
| `paper-convergence` | Selected resolution/timestep/threshold studies | Local cluster or modest HPC |
| `paper-standard` | Published-resolution selected MKS24 matrix; defined and guarded by explicit authorization | HPC |
| `paper-nulim` | Figure 13 limiter-frequency suite; defined and guarded by explicit authorization | HPC |
| `paper-analyze` | Regenerate implemented diagnostics/figures and optionally evaluate checksum/uncertainty-qualified reference curves; populated paper-panel comparisons remain open | Analysis host |
| `paper-summary` | Regenerate human/machine summaries without plotting | Any host |

The exact CLI naming may change, but avoid overloading `full` with runs that
are orders of magnitude more expensive.

### Manifest Schema Evolution

Schema version 2 now records development-run worktree status and model/forcing
choices for the reduced `paper-smoke` workflow. Every eventual paper bundle
must record:

```text
schema_version
workflow
status
created_utc
git_revision
git_worktree_dirty
git_worktree_status
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

A tracked implementation is now provided at:

```text
scripts/frontier/cgl_lf_frontier.py
scripts/frontier/build_cgl_lf_frontier.sh
```

`cgl_lf_frontier.py` is intentionally a preparation/accounting utility, not a
submission launcher. It requires the prescribed CGL root for real use,
reserves conservative node-hours, rejects walltimes greater than two hours
and paper-standard/limiter-production decks under `debug`, emits an
inspectable batch script, queries `squeue` before the user manually submits,
and records exactly one completed top-level allocation from `sacct` together
with the source/executable/input provenance held in its manifest. It also
archives an optional `--restart-file` for explicit continuation qualification
and supports releasing an unsubmitted reservation. GPU qualification remains
the default `--execution-target gpu` path. For a separately archived
MPI-enabled Kokkos-Serial executable, use `--execution-target cpu`; that path
records the target in the manifest, disables GPU-aware MPICH, and omits GPU
binding so CPU MPI evidence cannot accidentally be launched through the HIP
transport contract. For an explicitly scoped reduced sizing run that writes
shared MPI-I/O products, add `--mpiio-timers`; this records the choice in the
manifest and exports `MPICH_MPIIO_TIMERS=1` so Cray MPICH emits phase timings
at file close. Do not enable that measurement option silently in correctness
comparisons or paper-production definitions. Reduced sizing jobs `4660445`
and `4661434` use that opt-in path and retain reported shared-binary and
checkpoint write timings; the latter uses the exact standard mesh/block
layout only for startup sizing. Neither is a production run definition.
Validate its policy logic offline with:

```bash
python3 scripts/frontier/cgl_lf_frontier.py self-test
```

The shell/Python pattern below records the original design sketch; the tracked
utility supersedes it for actual campaign accounting and fills the required
provenance fields rather than leaving them blank.

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
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
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
- retain `MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0` for this AMD/HIP
  `HIPSpace` baseline; it does not use managed MPI buffers;
- record all environment settings in each manifest;
- validate correctness and memory behavior under the selected `HSA_XNACK`
  setting without presenting it as managed-MPI support;
- do not change communication or memory-tuning variables between comparison
  runs without treating that as a new run configuration;
- enable `MPICH_MPIIO_TIMERS=1` only through a manifest-recorded
  `--mpiio-timers` sizing run whose output path actually invokes MPI-I/O;
- do not load performance instrumentation packages in the baseline correctness
  build unless profiling is the stated purpose of the run.

Live qualification jobs from `4653415` through post-F-026 job `4658283`
reported `MPICH_GPU_SUPPORT_ENABLED = 1` but
`MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED = 0`, despite the earlier generated
scripts requesting the latter as `1`. Installed Cray MPICH 9.0.1
documentation states that this managed-memory setting defaults to `0` for
non-NVIDIA GPUs. AthenaK defines its CGL-LF device arrays through
`Kokkos::DefaultExecutionSpace::memory_space`, and the HIP backend maps that
space to `HIPSpace`, not `HIPManagedSpace`. The generated GPU launch now
requests the observed AMD/HIP baseline explicitly and makes no managed-MPI
claim. Policy-confirmation job `4659663` records both the requested and
runtime-reported values as GPU support `1` and managed-memory support `0`; its
histories and state/forcing tables are identical to pre-policy reference job
`4658072` at `rtol = atol = 1.0e-12`.
The initial run-environment log filter also failed to retain `MPICH_*`
exports; the tracked utility now captures these prefixes explicitly for
subsequent runs and retains the resolved policy confirmation.

### Improved Frontier Build Script Template

The tracked `scripts/frontier/build_cgl_lf_frontier.sh` uses the repository's
existing `built_in_pgens` configuration so the CGL-LF test and paper pgen are
available in one executable. Copy the reviewed script beneath
`$CGL_ROOT/scripts/build/` for a campaign and archive its output. The
following listing documents its required behavior.

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

The `prepare` action of `scripts/frontier/cgl_lf_frontier.py` generates a
concrete version of this correctness and reduced-validation script in each
run manifest directory. It intentionally uses `debug` and therefore is not a
paper-production submission template.

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
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
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
  env | LC_ALL=C sort | grep -E '^(MPICH_|FI_|HSA_|OMP_|ROCR_|HIP_|CRAY_)' || true
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

Retained live evidence as of 2026-05-25:

- The Frontier HIP/MPI executable built from committed revision `2a1b1fef`
  has SHA-256
  `d06be3e8dbed8f0b2e00f031b0f9274ae52b77c1b6c7e6861380d368becffbb8`.
- Combined G-001/G-002 job `4653415` completed on one node with eight GPU
  ranks, clean strict LF diagnostics, finite forcing/applied-heat-flux
  histories, and `0.003056` accounted node-hours.
- G-003 job `4653417` completed on one node with one GPU rank. Its final
  primitive/tab and history products agree with the eight-rank case under
  `rtol = atol = 1.0e-12`, with maximum absolute difference
  `2.0816681711721685e-17`, and it consumed `0.002778` node-hours.
- The executable rebuilt from revision `827394fb` has SHA-256
  `a82b1f6b14a75d991559038602ac8f5b549186914aa259eeba73618ed8945ef5`.
  G-004 job `4653457` passed the independent analytic finite-collision GPU
  oracle at `1.0e-12` relative tolerance with clean strict diagnostics and
  consumed `0.002500` node-hours. Its environment manifest retains the
  earlier requested/runtime managed-memory discrepancy later addressed under
  F-021.
- G-005 reference job `4653460` passed on the same executable. Continuation
  job `4653461` reproduced final physical fields within
  `2.8865798640254070e-15` and signed applied-work values within `4.0e-20`,
  but is rejected as gate evidence because `lf_nstage`/`lf_qface` lost
  pre-checkpoint cumulative baselines.
- The corrected executable from revision `2bc8d0e0` has SHA-256
  `e6333575e42397f462410bc41b1d72dad979a0254b93c65db8e2eb2f9780c824`.
  Corrected G-005 jobs `4653516` and `4653582` pass the finite-collision
  oracle, strict diagnostics, final-field restart comparison
  (`2.8865798640254070e-15` maximum absolute difference), and all 15
  then-existing cumulative LF endpoint comparisons
  (`3.9916427639358903e-20` maximum absolute difference). This closes F-022
  for that schema.
- G-006 jobs `4653586` and `4653601` execute the identical 16-MeshBlock
  active smoke case on one node/eight ranks and two nodes/sixteen ranks.
  Final state/forcing tables agree within `8.673617379884035e-19`, all
  retained history rows agree within `3.552713678800501e-15`, and strict
  unsafe counters remain zero.
- G-007 retains an active Alfvenic, passive Alfvenic, and active random
  reduced-paper matrix through jobs `4653586`, `4653731`, and `4653770`.
  Histories and forcing products are finite, strict unsafe counters are zero,
  and the random case has nonzero parallel and perpendicular forcing
  components.
- G-008 baseline evidence is retained in `g008_reduced_sizing_evidence.json`.
  The qualified reduced Athena steps used 3 to 5 seconds, reported `sacct`
  `MaxRSS` from `452492K` to `588224K`, and wrote checkpoints from `287185`
  to `817024` bytes. It is partial because runtime logs set
  `MPICH_MPIIO_TIMERS = 0`; G010b supplies the follow-up timing measurement
  rather than retroactively qualifying this baseline.
- The post-F-026 HIP/MPI executable built from revision `3e3c206b` has
  SHA-256
  `2abb31218cb68dbae0883589b12b8297b3d647ccd40798781c5f16011009d81e`.
  Active-random restart jobs `4658072`/`4658163` pass comparison at
  `rtol = atol = 1.0e-12`; the added `lf_cpwrk` and `lf_cawrk` endpoint
  differences are `1.1745964698252993e-22` and
  `1.1994118882018901e-22`, respectively, with zero strict failures.
- GPU AMR jobs `4658191`/`4658283` pass the post-F-026 one-rank/four-rank
  comparison at `rtol = atol = 1.0e-12`; their maximum MHD endpoint
  difference is `1.7763568394002505e-15`, both refine, and
  `lf_cpwrk`/`lf_cawrk` remain nonzero.
- The Kokkos-Serial/MPI executable built from revision `3e3c206b` has
  SHA-256
  `ab74f93d6ca3c8731c0b06d0f5c3abf669b577673dc6a21cecde2bb653870c91`.
  Initial CPU-target attempt `4658519` is retained as an F-027 launcher
  failure; corrected `--execution-target cpu` jobs `4659053`/`4659141` pass
  the one-rank/four-rank AMR comparison at `rtol = atol = 1.0e-12`, with
  maximum MHD endpoint difference `1.9539925233402755e-14` and zero strict
  failures.
- F-021 policy-confirmation GPU job `4659663` reports
  `MPICH_GPU_SUPPORT_ENABLED = 1` and
  `MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED = 0`, matching its generated
  launch environment. Its retained histories and state/forcing tabs are
  identical to job `4658072`.
- G010b reduced MPI-I/O sizing job `4660445` uses an archived
  `96 x 96 x 192` active deck with SHA-256
  `240c3b55a9ca0ba600a194f60caba0e272b2cf33e8426d6d42d65ef5c9f836ea`,
  reports `MPICH_MPIIO_TIMERS = 1`, and advances cleanly through
  `t = 0.01` with zero strict failures. Its two shared binary files are
  approximately `63.7` MB each and report mean net write bandwidth
  `1132.879` MiB/s; its fifteen approximately `301.8` MB restart files
  report mean net write bandwidth `1020.535` MiB/s. The retained
  `analysis/g010b_mpiio_sizing_evidence.json` records `4.654` GB total
  output from the deliberately excessive `dcycle = 1` measurement cadence,
  so these results constrain rather than define a production cadence.
- G011 standard-layout sizing job `4661434` uses an archived startup-only
  `192 x 192 x 384` active deck with SHA-256
  `ca643cdc659df8bea66f1c5aacfc03296605582cbc55a7396e64430762ae4af1`.
  It runs on one node/eight GPU ranks through `t = 0.01`, reports
  `MPICH_MPIIO_TIMERS = 1`, and has zero strict failures. Its maximum
  same-layout shared binary/restart file sizes are `509632020` and
  `1838131787` bytes, with mean net write rates `1772.938` and `1252.430`
  MiB/s, respectively. The retained
  `analysis/g011_standard_layout_mpiio_sizing_evidence.json` calculates a
  `328.915` GB raw eight-case storage envelope and a `255.735` GB
  sequential-retention allocation with 20 percent margin; it explicitly
  withholds runtime/node-hour projection.
- G012 nonlinear validation job `4662477` uses an archived reduced
  `96 x 96 x 192` active beta10 Alfvenic deck with SHA-256
  `517eb1e7c19cf45a41bdcb72faf0fd73fae2162d0c17598a483629c0277b86fc`.
  Unlike startup sizing, it selects the pre-F-033 finite-rate hard-wall
  `limiter_nu_coll = 1.0e10` model and requests `t = 2.0`. Its retained
  history is finite and has zero recorded strict counters through
  `t = 0.8001995575541284`, with forcing-work relative residual
  `2.832e-13`; the job then fails near `t = 0.802` with
  `hard_bound=1` during an LF split stage. Retained snapshots approach the
  physical parallel firehose threshold, reaching minimum
  `(p_perp - p_parallel)/B^2 = -1.0000009517035457` by
  `t = 0.7503056885529213`, before the unretained failing intermediate state
  crosses the emergency bound. This is F-031 and invalidates runtime
  extrapolation from this run.
- Follow-up local F-031 screening is retained beneath G012 in
  `analysis/f031_local_screen/f031_local_timestep_screen.json` with SHA-256
  `08ed73b001e4f7c84c930b287a68f0ec25e032da400b104f1bf32c059d01e7be`.
  At `32 x 32 x 64`, a fresh uncapped RKL2 reduced probe fails near
  `t = 0.946`; RKL2 continuations using `sts_max_dt_ratio = 20`, `10`, `5`,
  `2`, and `1` fail near `t = 1.038`, `1.146`, `1.197`, `1.572`, and
  `1.584`, respectively, despite zero strict counters in their last retained
  history rows. Thus the documented conservative RKL2 cap does not provide a
  qualified replacement. An explicit-LF continuation from the clean ratio-10
  `t = 1.0` checkpoint completes through `t = 1.6` with zero stored strict
  counters, but an explicit-LF switch from the clean ratio-2 `t = 1.5`
  checkpoint used by the failed ratio-1 branch also aborts near `t = 1.584`
  with `hard_bound=1`. The explicit result is therefore trajectory-dependent
  diagnostic evidence, not a demonstrated same-state RKL2 repair or a
  from-`t = 0`, GPU, or paper-statistics qualification.
- F-033 evidence is retained beneath G012 in
  `analysis/f033_hardwall_projection_local/f033_hardwall_projection_local_evidence.json`
  with SHA-256
  `29ce073e1059ac38496b0f543716b58b4cc5a3ff7a37c4457cb9217affca3eab`.
  The stage-context diagnostic locates the former same-state unlimited RKL2
  abort at `t = 1.500000`, `sweep=post stage=1/5`. Commit `099e3943` adds
  opt-in energy-preserving `limiter_hardwall`, a persistent `lf_hwproj`
  application counter (including primitive-refresh support/ghost ranges, not
  normalized occupancy), and selects it in hard-wall paper decks while
  leaving finite-rate `20`/`200` definitions unchanged. A deterministic
  restart-schema upgrade of the clean F-031 `t = 1.5` checkpoint changes only
  the added input parameter and appended diagnostic baseline; the corrected
  unlimited-RKL2 continuation passes `t = 1.5` to `1.6` to `2.0` with zero
  `lf_dfloor`, `lf_pfloor`, `lf_nonfin`, `lf_nonpos`, and `lf_hardbd`, and
  final `lf_hwproj = 53829433`. This is local same-state qualification;
  corrected Frontier GPU execution remains required.
- The retained ledger records `0.273892` node-hours used with no active
  reservation, including `0.188889` node-hours for failed G012. Paper-scale
  execution, reference-panel comparison, representative runtime/node-hour
  costing, and an approved production proposal remain blocked or open as
  identified above.

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

G012 demonstrates that the finite-rate limiter-active nonlinear model cannot
provide this performance measurement: it aborts on a strict split-stage
emergency-bound crossing before its measurement window. The retained local
F-031 screen rejects RKL2 caps through ratio `1` in the tested continuation,
including the documented conservative cap, and gives trajectory-dependent
partial explicit transport evidence: the same-state explicit switch also
fails. F-033 supplies a locally tested constrained hard-wall correction, but
does not convert old G012 timing into representative timing or establish GPU
behavior. Do not infer a production rate from failed G012 or local runs;
launch only a corrected reduced Frontier GPU qualification before renewed
sizing or costing.

### Preliminary Defined-Matrix Storage Envelope (Not Authorization)

Job `4661434` is a startup-only, one-node/eight-GPU measurement at the exact
`192 x 192 x 384` mesh and `32 x 32 x 64` MeshBlock layout used by the eight
currently defined `paper-standard`/`paper-nulim` inputs. It reaches only
`t = 0.01`; it validates standard-layout launch, host-memory evidence, shared
MPI-I/O file size/timing, and strict diagnostics, not scientific production or
runtime/node-hour cost.

Those eight guarded inputs now explicitly select shared `bin` and `rst`
outputs. Their common definition writes binary snapshots every `0.25` and
restarts every `1.0` through `tlim = 10.0`. Retained reduced runs confirm
inclusive initial/final time-cadenced output, giving 41 binary snapshots and
11 restart files per case. Using the larger of the two exact G011 files of
each kind:

| Quantity | Calculation basis | Bytes | Decimal GB |
| --- | --- | ---: | ---: |
| One shared binary snapshot | G011 maximum measured file | `509632020` | `0.510` |
| One shared restart file | G011 maximum measured file | `1838131787` | `1.838` |
| Raw generated output for one defined case | `41 bin + 11 rst` | `41114362477` | `41.114` |
| Raw generated output for all eight defined cases | `8 x one case` | `328914899816` | `328.915` |
| Retained output after acceptance | all binary snapshots plus final two restarts per case | `196569411152` | `196.569` |
| Sequential peak before pruning accepted checkpoints | seven retained cases plus one raw current case | `213112597235` | `213.113` |
| Sequential retained-policy allocation with 20 percent margin | `1.2 x sequential peak` | `255735116682` | `255.735` |
| No-pruning allocation with 20 percent margin | `1.2 x all raw output` | `394697879779` | `394.698` |

The candidate retention policy is to run cases sequentially, retain all
analysis snapshots and the final two restarts for every accepted case, and
remove earlier checkpoints only after acceptance and archive verification are
recorded. The authoritative retained calculation is
`analysis/g011_standard_layout_mpiio_sizing_evidence.json` beneath the G011
run directory.

This resolves the storage-volume portion of a proposal for the currently
defined cadence. It does **not** establish a production node count, late-time
timestep behavior, walltime, node-hour cost, or whether additional cases are
needed after reference-panel mapping. Those items require representative
runtime evidence and scientific review before any production submission can
be proposed for approval; the F-033 corrected hard-wall model must first pass
reduced Frontier GPU qualification because G012 failed before the intended
nonlinear measurement window.

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

For the eight currently defined inputs and their explicit shared-I/O cadence,
the preliminary G011 storage envelope above defines the candidate sequential
retention plan and a `255.735` GB allocation including margin. A reviewed
production proposal must revisit that envelope if the case matrix, cadence,
file layout, or retention policy changes, and must still supply the missing
runtime/node-hour basis after F-031 is resolved and requalified.

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
- collision activity and forcing-power proxy;
- manifest extensions.

Gate:

- smoke workflows;
- finite history validation;
- summary regeneration from existing bundle.

**Status (2026-05-25):** Implemented for reduced global histories, exact
operator-face heat-flux-cap fractions, smoke summary generation, and
manifest-window archived-run averages and generic diagnostic plots. Snapshot
analysis now adds a manifest-qualified heat-flux smoothing proxy; paper
history retains restartable cumulative RK-integrated applied forcing work and
active-case global energy residuals with a reduced-workflow closure gate. LF
history now retains RKL2-applied capped-face
heat-flux contractions with fine-side ownership on AMR interfaces. Snapshot
analysis now retains CGL
pressure-work components, their transfer comparison, and cadence-limited
time-quadrature estimates for pressure/heat-flux proxy rates. Optional
external reference-curve manifests now produce checksum/uncertainty-qualified
residuals and overlay figures for supported snapshot products and
threshold-volume histories. Populated MKS24 panel comparisons, exact
time-integrated local-work closure, and production residual qualification
remain open.

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

Use the actual repository commands present at that time. On 2026-05-25 the
nonconforming indentation in `docs/source/conf.py` was corrected, and an
isolated validation environment populated from `docs/requirements.txt` plus
`flake8`, `pytest`, and `numpy` was used because the active system Python does
not provide these gate dependencies. In that environment,
`python run_test_suite.py --style` passed (`2 passed`), as did
`sphinx-build -W --keep-going -b html source _build/html`.

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
scripts/stage_cgl_lf_mks24_reference.py
build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/manifest.json
build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/source/MKS24.tex
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

Decision implemented on 2026-05-24:

- apply collisions after each LF sweep with the sweep duration `dt/2`, so the
  two symmetric updates cover one physical evolution cycle.

Evidence:

- `test_cgl_lf_background_collision_advances_one_physical_timestep` checks
  the analytic relaxation law independently;
- the finite-collision `compare` and restart workflows pass in the focused
  CPU validation.

### Threshold Default

Decision implemented on 2026-05-24:

- preserve `cgl_firehose_threshold = oblique` as the compatibility default;
- require MKS24 inputs to state `cgl_firehose_threshold = parallel`.

Evidence:

- the policy regression distinguishes both thresholds;
- reduced paper inputs and workflow manifests state `parallel` explicitly.

### Emergency-Bound Definition Under Parallel Threshold

Decision implemented on 2026-05-24:

- report physical firehose occupancy at the selected threshold and reserve
  the emergency hard bound for `beta Delta <= -3`;
- count hard-bound violations independently of whether backup correction is
  enabled.

Evidence:

- `test_cgl_lf_strict_hard_bound_is_reported_without_backup_correction`
  protects truthful safety reporting;
- reduced paper smoke retains separate limiter occupancy and hard-bound
  histories.

### Pgen Integration With Forcing

Decision implemented on 2026-05-24:

- reuse the common `TurbulenceDriver` and extend it with deterministic
  restartable state, repaired RK coupling, explicit Alfvenic orientation, and
  opt-in MKS24 physical-shell/common-spectrum settings.

Evidence:

- focused tests cover seed selection, restart continuation, OU cadence,
  stage work, forcing orientation, and physical-shell input validity;
- reduced paper-smoke manifests archive the selected forcing contract, while
  long-time statistical calibration remains open.

### Reference Data for Quantitative Paper Comparison

Decision implemented for comparison ingestion on 2026-05-25:

- accept either machine-readable original data or explicitly digitized curves
  through an external manifest with source/data SHA-256 checksums and positive
  pointwise uncertainties;
- interpolate matching analyzed PDF, spectrum, transfer, or alignment products
  at supplied reference coordinates and retain normalized residuals.

Data still required:

- obtain machine-readable original data from authors/archive if available; or
- digitize plotted MKS24 reference curves with documented uncertainty.

### Paper Source and Figure Provenance on Frontier

Decision implemented for staging on 2026-05-24:

- use `scripts/stage_cgl_lf_mks24_reference.py` to download or copy the
  official arXiv `2405.02418v2` source archive into ignored result storage;
- retain its generated `manifest.json`, containing SHA-256 checksums for the
  archive, `MKS24.tex`, and extracted figures, beneath the Frontier CGL run
  root or in its archived provenance record.

Context:

- the local `docs/arXiv-2405.02418v2/` directory was reported as untracked
  when this plan was authored and is absent from the current live checkout;
- a pushed code branch alone will therefore not deliver that reference bundle
  to Frontier;
- simulations need not read the paper source, but figure comparison,
  provenance, and agent interpretation do require stable access to the exact
  reference version.

Remaining criterion:

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

## Immediate Next Blocking Tasks

Collision timing, threshold selection, truthful emergency monitoring, reduced
forcing cadence/source coupling, passive-Delta flow decoupling, pinned
paper-source staging, and explicit paper forcing-shell selection now have
local implementation evidence. The next critical path is:

1. **Calibrate forcing and use paper provenance.**
   Carry the staged official-source manifest into comparison products and
   populate documented, checksum- and uncertainty-qualified MKS24 reference
   curves; calibrate long-time injection/spectral behavior against the
   documented setup.

2. **Complete remaining comparison products.**
   Run paper-panel reference comparisons through the implemented manifest
   interface on top of the cap, forcing-work residual, applied
   heat-flux/pressure-traction ledgers, snapshot heat-flux/pressure-work
   products, and steady-window PDF/spectral/transfer/alignment/strain
   products and generic figure renderer.

3. **GPU-qualify the constrained nonlinear hard-wall model.**
   G012 job `4662477` and the retained F-031 local screen reject the
   finite-rate hard-wall attempt, unlimited RKL2, RKL2 caps `20`, `10`, `5`,
   `2`, and `1`, and a same-state explicit-LF switch for the screened
   continuation. F-033 implements opt-in energy-preserving algebraic
   threshold projection for paper hard-wall decks, retains strict monitoring,
   and passes the exact-state local `t = 1.5` to `t = 2.0` continuation with
   zero safety counters and active `lf_hwproj`. Run a corrected reduced
   Frontier GPU qualification with archived model, timestep, history, and
   accounting evidence. Only after that pass may G011's file-volume evidence
   be combined with renewed runtime/node-hour costing for a proposed
   production campaign. Scheduler-launched Frontier MPI-AMR decomposition
   evidence and the AMD/HIP launch-policy confirmation are retained; F-033
   GPU qualification now blocks representative late-time runtime/node-hour
   evidence.
