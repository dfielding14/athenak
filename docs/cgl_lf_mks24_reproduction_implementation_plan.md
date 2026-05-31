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
products, Figure 2(a)-coordinate joint pressure-density PDFs,
manifest-qualified cell-centered heat-flux smoothing proxies, and
snapshot pressure-work/heat-flux time-quadrature estimates, and generic
diagnostic figure rendering. A local `accuracy`
workflow now retains collisionless resolution/timestep, finite-collision,
extreme-cap, rotated-field, and low-field evidence. Paper forcing now retains
restartable cumulative RK-integrated applied source work for active global
energy-residual analysis. LF histories now also retain restartable, RKL2-applied capped-face
heat-flux contractions with fine-side ownership on AMR interfaces. Optional
pressure-work histories now retain RK-applied total and anisotropic CGL
traction work after the same AMR correction as momentum. Optional
reference-data manifests now provide checksum- and uncertainty-qualified
curve and sampled-surface comparison plumbing. A deterministic vector-path
extractor now populates
the eight MKS24 Figure 2(b) unstable-volume histories from the checksum-pinned
official PDF, and guarded standard definitions cover each associated
active/passive and Alfvenic/random beta-10/beta-100 case. Long-time energy
qualification, steady nonlinear convergence, remaining dimensional or
unmatched-panel reference data, and authorized standard comparison remain
open. Two guarded active beta-10
heat-flux definitions now cover the nonnominal Figure 12 parameter variants;
their nominal active/passive comparisons reuse standard definitions. A
tick-calibrated vector extractor and an `alignment_peak.cos_theta` analysis
product now populate the dimensionless Figure 9 standard-case and Figure 12
heat-flux alignment panels for future comparisons. A normalization audit found
that dimensional Figure 12 spectra
and Figure 13(a),(c) curves cannot be admitted until MKS24's reported
`p0 = 100` code-unit scale is transformed to AthenaK's `v_A = 1`
normalization. MKS24's stated
`T_total ~= E_K (2 pi u_rms/L_perp)` transfer normalization is now
implemented, so dimensionless Figure 7 lower-panel and Figure 13(d) transfer
ratios join the Figure 13(b) `beta Delta` PDFs as vector-extracted curves with
checksum provenance. Two guarded Figure 11 scale-separation definitions now
cover `n_perp = 96` and `384` around the existing standard `n_perp = 192`
case, and a checksum-bound extractor populates its dimensionless lower-panel
alignment curves while withholding its dimensional upper spectra. A
checksum-bound Figure 8 raster extractor now decodes its labeled linear
colorbar, verifies the paper's per-shell probability normalization within a
declared tolerance, and emits selected-shell dimensionless alignment PDFs
with conservative uncertainty. A Figure 4(a) vector extractor now maps its
explicit `rho/<rho>` coordinate onto `pdf.density_fluctuation` through
`x -> x - 1`; its companion `E_rho` spectrum remains excluded without a
qualified spectral ordinate definition. An opt-in deterministic,
local-field-conditioned three-point structure-function estimator now emits
`eddy_anisotropy.velocity_perp` and `eddy_anisotropy.magnetic_perp`, and a
checksum-bound Figure 5(b) extractor emits its four normalized eddy-scale
curves with conservative uncertainty. Figure 2(a)'s missing joint diagnostic
now emits and renders its plotted AthenaK-coordinate distribution, and the
reference-data contract accepts checksum-qualified sampled surfaces for those
joint PDFs. A checksum-bound Figure 2(a) extractor now decodes its labeled
tiled raster and colorbar, retains donor-pressure-coordinate surface samples
for audit, and emits admitted AthenaK-coordinate samples using the
source-qualified beta-10 pressure scale and two-dimensional PDF Jacobian.
Figure 3's
Fourier-projected compressive-flow spectrum and the Figure 4(b)/6(a)
normalized-density and separate thermal/magnetic-pressure spectral fields now
emit with explicit normalization metadata and generic rendering; their
published spectral ordinates remain excluded pending qualified transforms, and
two guarded `paper-compressive` definitions now supply Figure 3's missing
active random beta-1 and beta-100 sonic-correlation cases without authorizing
execution. A follow-up source/archive audit confirms that MKS24 and its cited
donor source define shell spectra and plotted units but do not state the
discrete Fourier normalization required for absolute ordinate conversion; no
linked public donor code or numeric-data package was identified. Those
spectral and strain panels therefore remain excluded pending external
normalization evidence or machine-readable reference data. The
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

### 2026-05-30 Frontier-root audit and clean restart decision

A read-only audit of `/lustre/orion/ast207/proj-shared/dfielding/CGL` was
performed on 2026-05-29 EDT (2026-05-30 UTC). The retained filesystem agrees
with the recorded qualification and Stage I ledgers:

- archived debug qualification use is `0.851670` node-hours with no active
  reservation;
- archived Stage I use is `9.962778` node-hours with no active Stage I
  reservation;
- the pre-replacement `R16` tree retains submitted inputs, rank-local restart
  ancestry, per-segment inspections, and diagnostic-prefix bundles through
  the accepted `t = 7.5` prefix plus rejected clean `t = 8.5` segment;
- the root also hosts exploratory non-MKS24 comparison campaigns beneath
  `runs/beta25-accel05-*`. Job `4743020` was cancelled and cleaned while
  retaining provenance. Two-node debug job `4743106` was running during the
  audit from revision `89835fea`; Slurm now records it as `COMPLETED 0:0`,
  while its top-level manifest remains stale at `state = running`.

Those exploratory jobs are not MKS24 evidence and are not represented in the
MKS24 ledgers. Do not modify or prune their directories, including the stale
completed-job manifest. The shared root requires explicit campaign isolation
and a broader pre-submission queue audit.

The first clean recovery path was a new MKS24 execution epoch:

```text
archived history: E01-pre-modal-driver
archived campaign: E02-modal-driver
archived prefix:   runs/mks24-stage-i/E02-modal-driver/
new planning cap: 4000 incremental node-hours
```

The fresh planning cap is a user-authorized accounting reset for future work.
It does not erase the historical `0.851670 + 9.962778 = 10.814448` consumed
node-hours. Do not reuse the old `R16` restart state, the old Stage I ledger
as a mutable ledger for new work, or the old storage/runtime estimate as an
`E02` reservation. Epoch-aware paths, cross-epoch rejection, shared-root
preflight, and an initial `4.0` node-hour `R16`-pilot-only reservation were
implemented locally in `scripts/frontier/cgl_lf_stage_i.py` with focused
offline regression coverage. The immutable replacement-driver HIP/MPI build
and qualification ladder are now retained. Fresh `E02` R16 jobs `4743735`
through `4744158` are accepted through exact `t = 10.0`, using `6.145556`
node-hours. Their completed lineage bundle passes `t = 8`--`10` production
analysis. Standard-layout R02 jobs `4744198` and `4744205` are accepted
through native snapshot boundary `t = 0.25`, use `0.473333` node-hours, and
project `19.000000` node-hours per standard `t = 10` case from the
authenticated continuation rate. High-resolution R17 jobs `4744210` and
`4744230` are accepted through authenticated continuation `t = 0.10`, use
`4.235556` node-hours, and project `42.488889` node-hours per simulated time
unit. F-071 authorizes only the frozen sixteen-case `R02`--`R17` matrix under
sequential inspection inside a measured `900.000000` node-hour E02 envelope.
R02 job `4744249` then reaches exact native restart boundary `t = 1.0` in
`1.452778` node-hours. F-072 updates the continuation-aware projection to
`702.195370` node-hours and records the two-hour Frontier normal-QOS limit.
R02 job `4744518` continues the inspected lineage through exact `t = 1.5` in
`1.052222` node-hours. F-073 updates the projection to `725.464938`
node-hours and authorizes only a conservative continuation through `t = 2.0`.
R02 job `4744913` continues the inspected lineage through exact `t = 2.0` in
`1.068889` node-hours. F-074 updates the projection to `730.081697`
node-hours and authorizes only a conservative continuation through `t = 2.5`.
R02 job `4745305` continues the inspected lineage through exact `t = 2.5` in
`1.061944` node-hours. F-075 updated the projection to `728.164877`
node-hours and historically authorized continued frozen-matrix execution one
inspected segment at a time, with R17 last. Continuation
`R02/s07_rankio_t2p5_t3` was prepared from the inspected `s06` terminal
restart siblings and submitted as job `4745498` on 2026-05-30.

The formal source-to-code execution signoff is retained in
`docs/cgl_lf_mks24_stage_i_protocol_review.md`. Its independent review found
that `E02` does not exactly implement the published forcing semantics: random
roles inherit a generic projection and implicit parabolic-spectrum default,
while planar roles do not generally enforce `grad_perp dot u_perp = 0` for
retained `k_z != 0` modes. Job `4745498` was therefore cancelled after `498`
seconds and recorded `aborted`, consuming `0.138333` node-hours and bringing
E02 cumulative use to `15.628610`. Preserve E02 products as pipeline and cost
evidence and keep E02 preparation fail closed. F-078 closes explicit-policy
implementation and corrected-build qualification for `E03-forcing-policy`.
The reviewed token is retained, reconciliation passed, and fresh
`R02/s00_rankio_t0_t0p1` job `4745922` was submitted from `t = 0`, formally
inspected, and recorded `accepted` at exact `t = 0.1` for `0.188889`
node-hours. Authenticated `R02/s01_rankio_t0p1_t0p25` job `4746154` is
also formally inspected and recorded `accepted` through exact `t = 0.25` for
`0.289167` node-hours. Corrected E03 Stage I use is `0.478056` node-hours.
Authenticated `R02/s02_rankio_t0p25_t1` job `4746182` is formally inspected
and recorded `accepted` through exact `t = 1.0` for `1.465556` node-hours.
Corrected E03 Stage I use is now `1.943612` node-hours. Authenticated
`R02/s03_rankio_t1_t1p5` job `4746356` is also formally and independently
inspected and recorded `accepted` through exact `t = 1.5` for `1.048611`
node-hours. Corrected E03 Stage I use is now `2.992223` node-hours. Retained
recost evidence SHA-256
`b6fd6e8fcd939ca1baa06411a3cd6f04359b3bc2ff5b2a0097e2cc87b1763fc7`
historically authorized only `R02/s04_rankio_t1p5_t2`. Authenticated
`R02/s04_rankio_t1p5_t2` job `4746435` is formally and independently
inspected and recorded `accepted` through exact `t = 2.0` for `1.063611`
node-hours. Corrected E03 Stage I use is now `4.055834` node-hours. Retained
recost evidence SHA-256
`c382b78ab9e21466648adf9d7ea38d2b407f0e986579f128bdfe6bc44bef79dd`
historically authorized only `R02/s05_rankio_t2_t2p5`. Authenticated
`R02/s05_rankio_t2_t2p5` job `4746663` is formally and independently
inspected and recorded `accepted` through exact `t = 2.5` for `1.062500`
node-hours. Corrected E03 Stage I use is now `5.118334` node-hours. Retained
recost evidence SHA-256
`39115585cf7ca866e3aade1e53c69aab6ee9217361d75db7946748905e1c0e5c`
historically authorized only `R02/s06_rankio_t2p5_t3`. Authenticated
`R02/s06_rankio_t2p5_t3` job `4746773` is formally and independently
inspected and recorded `accepted` through exact `t = 3.0` for `1.077500`
node-hours. Corrected E03 Stage I use is now `6.195834` node-hours. Retained
recost evidence SHA-256
`d2c5e27553426beb03c8c1882bd6473b682029c6ba67acd32bb06c5fb9a3ff01`
historically authorized only `R02/s07_rankio_t3_t3p5`. Authenticated
`R02/s07_rankio_t3_t3p5` job `4746953` is formally and independently
inspected and recorded `accepted` through exact `t = 3.5` for `1.151111`
node-hours. Corrected E03 Stage I use is now `7.346945` node-hours. Retained
recost evidence SHA-256
`468e787db48d98661cd3ffe125829a78f135bb7ca4c9104c766a00b693f586a3`
historically authorized only `R02/s08_rankio_t3p5_t4` next under its reviewed
one-segment `4500`-second threshold. Authenticated
`R02/s08_rankio_t3p5_t4` job `4747015` is formally and independently
inspected and recorded `accepted` through exact `t = 4.0` for `1.200000`
node-hours. Corrected E03 Stage I use is now `8.546945` node-hours. Retained
recost evidence SHA-256
`dd5b0f5b82cf81005c1a481ee83d1127aa43861dcdd8c22765170e28b0ad6c99`
historically authorized only `R02/s09_rankio_t4_t4p5` next with Slurm
walltime `01:40:00`, Athena timeout `01:30:00`, and a reviewed one-segment
`4800`-second threshold retaining `600` seconds on both timeout margins.
Authenticated `R02/s09_rankio_t4_t4p5` job `4747087` is formally and
independently inspected and recorded `accepted` through exact `t = 4.5` for
`1.208056` node-hours. Corrected E03 Stage I use is now `9.755001`
node-hours. Retained recost evidence SHA-256
`8b10b1786db520740b687d989207f3cda6e0123f9bf2b62e75a5df3849511561`
historically authorized only `R02/s10_rankio_t4p5_t5` next with Slurm walltime `01:40:00`,
Athena timeout `01:30:00`, and a reviewed one-segment `4800`-second threshold
retaining `600` seconds on both timeout margins. The threshold is scoped only
to `s10` and does not ratchet automatically.
Authenticated `R02/s10_rankio_t4p5_t5` job `4747146` is formally and
independently inspected and recorded `accepted` through exact `t = 5.0` for
`1.226389` node-hours. Corrected E03 Stage I use is now `10.981390`
node-hours. Retained recost evidence SHA-256
`5256fb7fd9b75f175a16c2c16ee8a995b8d6d9281a92b70a979ea62171835bf2`
historically authorized only `R02/s11_rankio_t5_t5p5` next with Slurm walltime `01:40:00`,
Athena timeout `01:30:00`, and a reviewed one-segment `4800`-second threshold
retaining `600` seconds on both timeout margins. The threshold is scoped only
to `s11` and does not ratchet automatically.
After clean `s11` inspection, accounting, reconciliation, and recost, propose quarter units if elapsed time exceeds `4800` seconds, remaining cap headroom falls below `300` seconds, or the reviewed trend estimator projects the next half-unit at or above `4800` seconds. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation.
Authenticated `R02/s11_rankio_t5_t5p5` job `4747202` is formally and
independently inspected and recorded `accepted` through exact `t = 5.5` for
`1.268889` node-hours. Corrected E03 Stage I use is now `12.250279`
node-hours. Retained recost evidence SHA-256
`323d7cb3cbe108eee944045c189cff75ffc07d06833c49bcfa8315f07fad8d51`
authorizes only `R02/s12_rankio_t5p5_t5p75` on one node from the authenticated
`s11` terminal siblings, with Slurm walltime `01:05:00`, Athena timeout
`00:55:00`, and a reviewed one-segment
`2700`-second threshold retaining `600` seconds on both timeout margins. The
threshold is scoped only to `s12` and does not ratchet automatically. Any
scientific, provenance, scheduler, storage, budget, or reconciliation failure
blocks successor preparation.

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

At the 2026-05-29 EDT audit, the live branch baseline was:

```text
branch: feature/cgl-landau-fluid
HEAD:   89835fea Add near-isothermal MHD and pure-CGL turbulence comparison
```

The live worktree also contained an existing modified `kokkos` entry and an
untracked `vis/python/__pycache__/` directory. Treat both as user or generated
state; do not rewrite or stage them as part of MKS24 planning work.

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

| Gate | Required evidence | Release-blocking? | Current status |
| --- | --- | --- | --- |
| Physical equations and normalization | Equation-to-code review, closure/cap tests, documented code-unit mapping | Yes | Partial; collisionless mapping reviewed |
| Collision source-term timing | Independent analytic regression and operator-ordering documentation | Yes | Implemented 2026-05-24; focused CPU test passed, broader gates pending |
| Instability-threshold policy | Explicit MKS24 policy, alternative-policy tests, archived input choice | Yes | Implemented 2026-05-24; CPU policy test and active paper-smoke manifest passed |
| Safety diagnostic truthfulness | Strict tests independent of backup correction | Yes | Implemented 2026-05-24; focused CPU test passed, broader gates pending |
| CPU numerical accuracy | Convergence, asymptotic, cap, directional, restart, MPI/AMR checks | Yes | Partial 2026-05-30; archived `accuracy-v2` bundle covers N-001 through N-007 locally, the AMR workflow passed, and the pre-replacement full CPU suite passed (`218 passed, 15 skipped`). The focused post-merge serial CGL gate passed (`39 passed`) with F-056 layout-independent cap-fraction coverage, fixed-grid/AMR modal-driver restart identity, and CGL-LF interaction regressions. Scheduled one-node MPI CPU job `4743666` passed both current one-rank/four-rank decomposition regressions (`2 passed`). Production/statistical comparison gates remain open. |
| Frontier GPU numerical equivalence | CPU/GPU comparison suite, MPI/GPU restart checks, strict diagnostics | Yes for Frontier use | Corrected-build reduced qualification closed 2026-05-30 for revision `9e07542281e4e6d125582f253df3ad2e3b8b154d`. E03 `g024`--`g031` qualify explicit planar/random policies, post-refresh modal restart identity, one-rank/eight-rank decomposition identity, passive semantics, reduced nonlinear hard-wall behavior, and standard-layout ranked startup/sizing. Long-time production statistics remain open. |
| Paper pgen and forcing fidelity | Input/pgen review, forcing metadata, restartable reduced smoke | Yes for MKS24 claims | Corrected-build reduced qualification closed 2026-05-30. F-076 remains the historical E02 stop line. F-077 adds explicit restart-retained `mks24_random_unprojected` and `mks24_alfvenic_perpendicular` policies plus explicit paper-deck `spectrum = power_law`; E03 `g024`/`g025` directly qualify the force-family contracts and `g026`/`g027` qualify modal restart identity across an OU refresh. Fresh canonical production began from E03 `t = 0`; accepted R02 jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, and `4747202` now reach exact `t = 5.5`. |
| Paper observables and analysis | Synthetic analysis tests and archived reduced-run products | Yes for MKS24 claims | Partial 2026-05-25; reduced histories, restartable RK-integrated applied forcing work/global active residual, RKL2-applied heat-flux contractions with fine-side AMR ownership, AMR-corrected RK-applied total/anisotropic CGL pressure-work ledgers, operator-face cap counters, windowed snapshot PDF/spectral/transfer/alignment/local-strain products, Figure 2(a)-coordinate joint pressure-density PDFs and admitted checksum-qualified sampled surfaces, Figure 3 compressive-flow spectra, Figure 4(b) normalized-density spectra, Figure 6(a) thermal/magnetic-pressure spectra, MKS24-normalized transfer and alignment-peak curve comparisons, opt-in local-field eddy-anisotropy structure functions, deduplicated threshold-volume history curves, manifest-qualified proxies and cadence-limited estimates, generic figures, checksum/uncertainty-qualified reference-data comparison plumbing, Figure 2(b) histories, Figure 4(a) normalized-density PDFs, Figure 5(b) normalized eddy scales, Figure 7 lower-panel transfer ratios, Figure 8 selected-shell alignment PDFs, Figure 9, Figure 11 lower-panel, and Figure 12 alignment curves, and Figure 13(b),(d) dimensionless curves are implemented; guarded Figure 3 case definitions are present under `paper-compressive`; dimensional/unmatched-panel conversion, standard-run comparison, and remaining panel references remain open |
| Standard MKS24 results | Required cases, durations, manifests, figure comparisons | Yes for reproduction claim | Only the admissible fresh E03 R02 prefix has run through exact `t = 5.5`; no corrected-epoch case is complete and long-time mapped results remain not run. The frozen sixteen-case mapped Stage I manifest remains the target matrix. Preserve E02 jobs `4743735` through `4744158` only as pipeline and cost evidence. |
| Operational workflow | Frontier scripts, budget ledger, failure recovery, storage plan | Yes for Frontier use | Corrected E03 entry gate passed. The debug-helper ledger records `1.847784` node-hours after `g031`; corrected E03 qualification accounts for `0.485834`. E02 records `15.628610` historical pipeline/cost node-hours and remains fail closed. Token SHA-256 `e3fec9f35da42121b902f41ef752021f375aab38bc34c5ff75ae8780bfbca635` is retained. Fresh E03 R02 jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, and `4747202` were formally inspected and recorded `accepted` through exact `t = 5.5` for `12.250279` cumulative node-hours. F-083/F-084 hardening is promoted at canonical revision `9689c269bf329542815a1b2b137881126964b05c`; hardened reconciliation passed with `12/12/12` ledger rows/manifests/reservations, no active reservation, and no transaction. |
| User documentation | Sphinx build and accurately scoped runbook | Yes | Implemented for current functionality 2026-05-25; Sphinx warnings-as-errors and repository style suite pass in an isolated validation environment; future campaign results must still be documented when executed |
| Performance suitability | Representative timing/memory/I/O evidence; no uninvestigated prohibitive bottleneck | Yes for production use | Partial 2026-05-31. E03 `g030b` reaches reduced nonlinear `t = 2.0` in `1516` allocated seconds with fourteen MPI-I/O records. E03 `g031` launches the exact `192 x 192 x 384` ranked layout, fits one node, and reaches startup-only `t = 0.01` in `77` seconds while retaining `509685800` terminal snapshot bytes and `1381369688` terminal checkpoint bytes. E03 R02 continuation job `4747202` measures `2.537778` node-hours per simulated time unit across `t = 5.0`--`5.5`; retained corrected recost evidence projects at most `786.323334` matrix node-hours. The clean `4568`-second segment leaves only `232` seconds below its reviewed `s11` cap, so independent review authorizes only `R02/s12_rankio_t5p5_t5p75` on one node from authenticated `s11` terminal siblings with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a one-segment `2700`-second ceiling retaining `600` seconds on both timeout margins. It does not ratchet automatically. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation. Archive and catalog the current committed controller state, reconcile, prepare only `R02/s12_rankio_t5p5_t5p75` on one node from authenticated `s11` terminal siblings under the full F-094 profile, then inspect, account, reconcile, and recost. |

The final readiness report must distinguish:

- **production ready for supported CGL-LF use**;
- **MKS24 reproduction complete**;
- **validated on Frontier GPU hardware**;
- **future extensions not yet supported**.

None of these phrases should be used as a substitute for the others.

## Current Readiness Decision (2026-05-30)

This is the current readiness report for the implemented and retained
evidence on this branch. It is a blocked release decision, not a completed
reproduction claim.  Paper-production execution is governed separately by
the subsequently recorded sequential Stage I protocol and its `4000`
node-hour ceiling in `docs/cgl_lf_weak_guide_manuscript_plan.md`. The
2026-05-30 revision resets that ceiling incrementally for `E03-forcing-policy`
while preserving all historical usage records.

| Claim | Decision | Retained supporting evidence | Blocking work |
| --- | --- | --- | --- |
| Local operator and reduced-workflow implementation | Supported for the exercised replacement-driver scope; production comparison still required | Full CPU suite (`218 passed, 15 skipped`), retained historical focused suites and bundles, F-018 through F-026 regression evidence, F-033 local exact-state hard-wall evidence, current focused serial CGL gate (`39 passed`), and scheduled current MPI CPU job `4743666` (`2 passed`); the retained current gate covers modal restart, CGL-LF turbulence-driving AMR interaction, and F-056 layout-fraction accounting | Complete production/statistical comparison gates |
| Validated on Frontier GPU hardware | Supported for corrected reduced E03 qualification through standard-layout ranked startup/sizing | Immutable corrected revision `9e07542281e4e6d125582f253df3ad2e3b8b154d` and executable SHA-256 `68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c`; retained E03 `g024`--`g031` explicit-policy, restart, decomposition, passive, nonlinear-hard-wall, and sizing evidence; archived earlier GPU records remain historical context | Long-time corrected mapped production runtime and statistical qualification remain open; G012 remains the retained failure of the superseded finite-rate hard-wall interpretation |
| Production ready for supported CGL-LF use | Not established | Corrected E03 debug qualification, retained token, token-gated accounting tooling, committed epoch isolation and shared-root controls, immutable corrected build archival, retained focused serial/MPI CPU gates, historical E02 pipeline/cost evidence, accepted fresh E03 R02 jobs `4745922`, `4746154`, `4746182`, `4746356`, `4746435`, `4746663`, `4746773`, `4746953`, `4747015`, `4747087`, `4747146`, and `4747202` through exact `t = 5.5`, and archived F-079--F-094 controller provenance plus hardened clean reconciliation exist | Archive and catalog the current committed controller state, reconcile, prepare only the authenticated R02 continuation from authenticated `s11` terminal siblings through `t = 5.75` under the reviewed non-ratcheting `s12` quarter profile, inspect, account, reconcile, and recost; any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation; then execute representative mapped production and complete production/statistical qualification |
| MKS24 reproduction complete | Not established | Guarded standard/limiter/heat-flux/compressive/scale-separation input matrices and paper analysis products including Figure 2(a)-coordinate joint PDFs with source-qualified sampled-surface references, Figure 3/4(b)/6(a) spectral fields, checksum/uncertainty-qualified reference-data comparison interface, populated Figure 2(b) histories, Figure 4(a) normalized-density PDFs, Figure 5(b) normalized eddy-scale curves, dimensionless Figure 7 transfer ratios, Figure 8 selected-shell alignment PDFs, Figure 9, Figure 11 lower-panel, and Figure 12 alignment curves, and Figure 13(b),(d) curves exist. E02 R16 pipeline analysis is retained but is not direct reproduction evidence. | Execute the frozen mapped Stage I matrix from fresh E03 `t = 0`; resolve remaining dimensional or unmatched spectral transforms; complete quantitative panel comparisons and archived scientific interpretation |

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
  Corrected reduced Frontier job G013 (`4673404`) subsequently passes from
  `t = 0` through `t = 2.0` under unlimited RKL2 with zero strict safety
  counters and active `lf_hwproj`; this qualifies reduced GPU correctness,
  not paper production or long-time statistics.
- Replacement-driver debug job `4743463` (`g014`) failed before physics
  startup because its archived sizing deck lacked the
  `mhd/limiter_hardwall` key targeted by a runtime override. The failure was
  accounted as `0.002778` node-hours. F-061 adds a local prepare-time
  absent-key rejection and retains a corrected explicit qualification deck.
  Job `4743465` (`g015`) then reaches `t = 0.01` with zero strict counters,
  finite readable rank-local snapshots, and complete eight-rank restart sets.
  Evidence JSON SHA-256:
  `e2f4a047f601e6d1539d77a85a09d65e4df8b338c22eb6638441cdd6af75cd70`.
  This is startup/layout evidence only.
- Rank-local modal restart-load smoke `g016` (job `4743470`) reaches
  `t = 0.02` cleanly from the complete `g015` terminal checkpoint set.
  Uninterrupted comparator `g017` (job `4743472`) is retained as
  `inconclusive`: the forced exact `g015` stop at `t = 0.01` clipped its final
  timestep, so `g016` and uninterrupted `g017` intentionally follow different
  integration sequences. Natural-cycle comparison `g018` (job `4743473`)
  resumes the complete `g017` `.00013.rst` set and matches the uninterrupted
  terminal histories, assembled rank-local float32 snapshot, and turbulence
  force slice within declared retained-format tolerances. Evidence JSON
  SHA-256:
  `4d898cf4aa4427b81164e4eb688f011adbbb277874b1ca7453a7a2c1e60b580f`.
- One-rank comparator `g019` (job `4743474`) follows the same uninterrupted
  25-cycle trajectory as eight-rank `g017` through `t = 0.02`, with zero
  strict LF failure counters and complete rank-local products. Comparing
  assembled physical grids rather than rank-local shard boundaries gives
  maximum terminal MHD-history, user-history, float32 snapshot-field, and
  turbulence-force-slice differences `4.44e-17`, `7.11e-15`, `3.73e-9`, and
  `8.88e-16`. Evidence JSON SHA-256:
  `509ec34d4f1f27dfdbb9b507d2ca663d7233d168a0b83775b0ab532dc283c752`.
- Reduced passive-Alfvénic smoke `g020` (job `4743480`) reaches `t = 0.01`
  with zero strict counters, finite histories and force slice,
  perpendicular-only forcing (`force_prl2 = 0`, `force_prp2 = 128`), exactly
  zero passive `lf_cpwrk`/`lf_cawrk`, and three retained checkpoints. Reduced
  active-random smoke `g021` (job `4743483`) retains
  `driving_type = 0`, `rseed = 314159`, activates all three force components,
  produces finite nonzero `lf_cpwrk`/`lf_cawrk`, and retains four checkpoints.
  Evidence JSON SHA-256 values:
  `4abeddeb25db489b3334be240b837e703d035585d0814fc79d004e273cc34b3f`
  and `aaee51f186ac0f223fa4a1c35b238bb1c2e22f2a1ac9075ad1e28ffea53e8498`.
- Reduced nonlinear hard-wall run `g022` (job `4743484`) reaches `t = 2.0`
  in `1554` allocated seconds with zero strict LF failure counters, terminal
  `lf_hwproj = 10055979902`, finite forcing/LF/pressure-work ledgers,
  forcing-work relative residual `1.863e-12`, nine readable shared
  snapshots, five retained shared checkpoints, and fourteen MPI-I/O
  net-write records. Its terminal modal checkpoint is `227140784` bytes,
  smaller than the historical pre-modal G013 payload. Evidence JSON
  SHA-256:
  `339a3309276a93be7bb55c196f367ead9754f7c2ebdf94afdb9aeedd37f026cf`.
  This closes reduced nonlinear correctness only; it does not establish
  production runtime or long-time statistics.
- Standard-layout rank-local sizing run `g023` (job `4743662`) launches the
  guarded `192 x 192 x 384` layout with eight GPU ranks, advances through
  startup-only `t = 0.01` in `79` allocated seconds, and uses `0.021944`
  node-hours. Strict LF counters remain zero, forcing-work relative residual
  is `1.053e-11`, and assembled initial/final snapshots are readable at
  shape `384 x 192 x 192`. It retains two complete eight-rank binary groups
  with terminal aggregate size `509682248` bytes and two complete eight-rank
  modal checkpoint groups with aggregate size `1381366136` bytes. Evidence
  JSON SHA-256:
  `effba10b616246e1c9f6863c6a068f674551fa9ac4b4a02bdb92ead0f963c14f`.
  This closes standard-layout startup and retained-output sizing only.
- The focused post-merge local CGL gate is retained separately. The serial
  suites pass `39` tests in `23.89` seconds; evidence JSON SHA-256 is
  `3fc869c0d7455f1221d7838de900c0db7426673b2730cd954f264b8934fc0c29`.
  Scheduled one-node MPI CPU job `4743666` passes both one-rank/four-rank
  decomposition regressions in `2.81` pytest seconds and consumes
  `0.002222` node-hours outside the debug-helper ledger; evidence JSON
  SHA-256:
  `397aced4e90aa7e77ab17da05a522ad4c99c9d0e15bd48d86ce20aa195fb36d3`.
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
| F-006 | 2026-05-24 | Phase D | The plan named a local MKS24 source/figure bundle that is absent from the live worktree and should not be copied into the code branch without separately establishing redistribution permission | official arXiv `2405.02418v2` source, `scripts/stage_cgl_lf_mks24_reference.py`, staged `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/manifest.json` | High | Stage the pinned official source in ignored result storage, safely extract it, record archive/file SHA-256 checksums, and attach checksum provenance during paper analysis | Staging utility executed with archive SHA-256 `6b58e1ed01585c9820659a01b7549a08a09b7ea3fcbb28dbb96205259e82da8e`; `paper-convergence-v1/analysis/reference_provenance.json` | Implemented for source provenance; F-034/F-036/F-038/F-039/F-040/F-041/F-043/F-044 populate currently admitted reference-safe curves, while remaining panel curves/comparisons remain open |
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
| F-024 | 2026-05-25 | Phase E | The staged official source provides rendered panel PDFs rather than numeric curves, and the analyzer lacked a fail-closed contract for separately obtained or digitized reference data | `scripts/analyze_cgl_lf_paper.py`, `scripts/plot_cgl_lf_paper.py`, `scripts/cgl_lf_workflow.py` | High | Add optional reference-curve manifest ingestion requiring provenance, source/data SHA-256 values, and positive pointwise uncertainties; compare mapped products and render residual-qualified overlays | Generated exact-curve and wrong-source-checksum regressions in `test_cgl_lf_paper_snapshot_analysis_uses_both_pressures`; retained-bundle workflow rendering probe | Implemented for comparison infrastructure; F-034/F-038/F-039/F-040/F-041/F-043/F-044 provide currently admitted reference-safe panel curves, while dimensional/unmatched panels and production comparison remain open |
| F-025 | 2026-05-25 | Phase E | `force_work` directly summed each explicit-RK source-call energy change even though the final RK state reweights earlier stage source increments; retained multi-cycle reduced bundles consequently reported false active global residuals while passing | `src/srcterms/turb_driver.cpp`, `scripts/cgl_lf_workflow.py`, `20260525-paper-convergence-face-work-v1` manifest | Blocker | Advance cumulative forcing work as a companion RK state using the cycle baseline and `gam0`/`gam1`, and reject active reduced cases whose energy/work relative residual is at least `1.0e-8` | `test_cgl_lf_paper_multicycle_forcing_work_follows_rk_state_recurrence`; regenerated `20260525-paper-smoke-rk-work-v1` and `20260525-paper-convergence-rk-work-v1` bundles | Implemented locally; focused CPU suite passes (`27 passed`), all six corrected convergence cases pass with maximum relative residual `5.141362829086986e-10` |
| F-026 | 2026-05-25 | Phase E | Snapshot-reconstructed pressure work could not identify the pressure traction actually applied by hyperbolic CGL momentum fluxes, especially after FOFC replacement and AMR flux correction | `src/mhd/rsolvers/hlle_cgl.hpp`, `src/mhd/mhd_tasks.cpp`, `src/mhd/mhd_fofc.cpp`, `src/outputs/history.cpp`, `scripts/analyze_cgl_lf_paper.py` | Blocker | Retain total and anisotropic CGL pressure traction in the numerical flux, apply the same AMR flux correction as momentum, contract corrected divergence with stage velocity, and advance `lf_cpwrk`/`lf_cawrk` by the explicit-RK companion recurrence; retain zero applied work for passive-Delta | Focused CPU AMR/restart/passive checks, dedicated CGL FOFC traction check, expanded MPI AMR decomposition oracle, regenerated local pressure-work bundles, and retained Frontier post-F-026 GPU/CPU evidence | Implemented for reduced qualification: local focused CGL CPU (`27 passed`) and style (`2 passed`) succeed; Frontier GPU restart jobs `4658072`/`4658163` and GPU AMR jobs `4658191`/`4658283` pass at `rtol = atol = 1.0e-12`; MPI CPU AMR jobs `4659053`/`4659141` also pass; production interpretation remains open |
| F-027 | 2026-05-25 | Frontier MPI CPU rerun | The debug preparation utility unconditionally generated GPU-aware MPICH exports and GPU binding; launching the intentionally non-GTL-linked Kokkos-Serial MPI executable consequently aborted before simulation with `MPIDI_CRAY_init: GPU_SUPPORT_ENABLED is requested, but GTL library is not linked` | `scripts/frontier/cgl_lf_frontier.py`, job `4658519` log and manifest | High | Add explicit `--execution-target gpu/cpu` generation; retain the existing GPU-aware path by default, while the CPU path omits GPU Slurm/srun binding and forces `MPICH_GPU_SUPPORT_ENABLED=0` | Expanded `cgl_lf_frontier.py self-test` and repeated scheduler-launched MPI CPU AMR comparison | Implemented and confirmed live: corrected jobs `4659053`/`4659141` run with `MPICH_GPU_SUPPORT_ENABLED=0`, refine cleanly, retain nonzero pressure work, and pass the one-rank/four-rank AMR oracle at `rtol = atol = 1.0e-12` |
| F-028 | 2026-05-25 | Frontier G-008 sizing | Original reduced sizing evidence retained output sizes but no MPI-I/O phase timing because runtime reported `MPICH_MPIIO_TIMERS = 0` | `scripts/frontier/cgl_lf_frontier.py`, `g008_reduced_sizing_evidence.json`, jobs `4660445`/`4661434` and their retained analysis JSON | Medium | Add manifest-recorded opt-in `--mpiio-timers`, then execute reduced and standard-layout startup shared-binary/checkpoint sizing runs; treat startup timing only as storage/I/O reconnaissance | Utility self-test; jobs `4660445`/`4661434`; retained G010b/G011 analysis JSON | Implemented for storage/I/O sizing: timers captured, strict diagnostics remain clean, and G011 supplies the measured file-size basis used by the proposal; representative standard runtime remains open |
| F-029 | 2026-05-25; superseded by F-054 | Production input audit after G011 | Defined standard/limiter decks relied on the default `single_file_per_rank = false`, while workflow manifests did not archive the output cadence/file-layout contract used for storage sizing | `inputs/cgl_lf_paper/`, `scripts/cgl_lf_workflow.py`, `analysis/g011_standard_layout_mpiio_sizing_evidence.json` | Medium | State the selected binary/restart layout explicitly in each guarded production deck and archive output layout/cadence choices in future workflow manifests | Renamed `test_cgl_lf_paper_production_inputs_explicitly_use_rank_local_io`; metadata probe over `paper-standard`/`paper-nulim` cases; F-054 rank-local retained-output evidence | The initial shared-file correction was implemented, then F-054 intentionally changed production snapshots/restarts to explicit `single_file_per_rank = true`. Any later layout/cadence revision requires renewed storage and I/O review. |
| F-030 | 2026-05-25 | Phase E reference-panel audit | MKS24 Figure 2(b) is an unstable-volume time history, but the checksum-qualified reference comparator admitted only snapshot-derived PDF/spectrum/transfer/alignment products | `scripts/analyze_cgl_lf_paper.py`, staged `MKS24.tex:508-512`, `docs/source/modules/cgl_landau_fluid_validation.md` | High | Export ordered threshold-volume history series, collapsing Athena's duplicate terminal row by retaining its last equal-time value, and admit them through `history.*` reference products | Expanded exact-curve regression in `test_cgl_lf_paper_snapshot_analysis_uses_both_pressures` for `history.unstable_fraction` | Implemented for comparison infrastructure; F-034 supplies the digitized histories and target-run comparisons remain open |
| F-031 | 2026-05-25 | Frontier G012 nonlinear hard-wall validation | The intended active beta10 hard-wall run with uncapped RKL2 STS reaches the parallel firehose limiter and then crosses the strict emergency bound in an LF split stage before its analysis window | G012 job `4662477` log, retained histories/snapshots, and `analysis/g012_hardwall_strict_failure_evidence.json`; `analysis/f031_local_screen/f031_local_timestep_screen.json` under the G012 run; `src/mhd/mhd_sts.cpp` | Blocker | Retain strict checking, reject tested RKL2 caps, diagnose an explicit-LF or revised splitting/integration policy for limiter-active nonlinear runs, and rerun reduced GPU qualification before any runtime costing or production proposal | Local reduced screen retains strict aborts for uncapped RKL2 and capped ratios `20`, `10`, `5`, `2`, and `1`; explicit LF is clean through `t = 1.6` only when branched from a clean `t = 1` checkpoint, while switching from the clean ratio-2 `t = 1.5` checkpoint also fails near `t = 1.584` | Addressed by F-033; corrected G013 job `4673404` passes reduced GPU qualification, while G012 remains the retained superseded-model failure |
| F-032 | 2026-05-25 | Workflow provenance after F-031 | Workflow manifests archived closure and output layout choices but omitted the STS timestep cap and integrator controls that determine limiter-active qualification | `scripts/cgl_lf_workflow.py`, `inputs/cgl_lf_paper/`, F-031 | High | Archive time integrator, STS integrator, `sts_max_dt_ratio`, and CFL number in future workflow manifests together with the F-033 hard-wall choice | Extended `test_cgl_lf_paper_production_inputs_explicitly_use_rank_local_io` metadata assertions; G013 submitted input and manifest | Implemented for provenance and exercised in corrected G013 GPU qualification |
| F-033 | 2026-05-25 | F-031 model correction/local and Frontier GPU validation | MKS24's hard-wall `nu_lim = 1.0e10` intent is not robustly represented by finite-rate post-LF scattering: the screened failing state violates the firehose emergency bound during internal RKL2 `post stage=1/5` | `src/eos/cgl_physics.hpp`, `src/eos/cgl_mhd.cpp`, `src/outputs/history.cpp`; local `analysis/f033_hardwall_projection_local/f033_hardwall_projection_local_evidence.json` beneath G012 (SHA-256 `29ce073e1059ac38496b0f543716b58b4cc5a3ff7a37c4457cb9217affca3eab`); G013 `analysis/g013_f033_hardwall_gpu_qualification_evidence.json` (SHA-256 `2381e89a24a4ff9c2d81e4e0e91d049c3b62ac0577f044fee1816e2ccd2d6108`) | Blocker | Add opt-in `limiter_hardwall` energy-preserving projection at CGL primitive recovery, replace finite-rate limiter pressure relaxation only in hard-wall mode, append persistent `lf_hwproj`, select it in hard-wall paper decks, and rerun corrected GPU qualification | `test_cgl_lf_hardwall_projects_to_selected_firehose_threshold`, `test_cgl_lf_hardwall_requires_instability_limiter`, expanded restart/input metadata assertions; focused CPU module `30 passed`; style/Sphinx gates; Frontier G013 job `4673404` | Implemented and qualified for reduced scope: local exact-state pass through `t = 2.0`, then G013 GPU pass from `t = 0` through `t = 2.0` with zero strict safety counters, terminal `lf_hwproj = 14763255866`, and forcing-work relative residual `9.627e-14`; production/statistical claims remain open |
| F-034 | 2026-05-25 | Figure 2(b) reference/matrix audit | The published unstable-volume panel contains eight active/passive, Alfvenic/random beta-10/beta-100 histories, but the guarded standard workflow mapped only four of those series and retained no actual numeric panel curves | pinned `source/fig2b.pdf` SHA-256 `828287244aa0ec72aad6e802268fbcb8d77af65f11b03fcf9634eb9cf741ab0e`; `scripts/digitize_cgl_lf_mks24_fig2b.py`; ignored `digitized_fig2b_v1/curves.json` SHA-256 `66996af7ddd48a39c1a8fc46e1fd7ba6206ced5625c2d30a123cc3df67df3df1` | High | Add the four missing guarded standard decks and workflow entries; extract the eight visible vector histories only from the pinned PDF into an external manifest with explicit digitization uncertainty and clipped-vertex handling | Expanded production-input/workflow-case assertion; checksum-validated eight-curve manifest ingestion probe | Implemented for Figure 2(b) reference and input coverage; authorized standard executions/comparisons and remaining figure reference curves remain open |
| F-035 | 2026-05-25 | Figure 12 case-matrix audit | The pinned source specifies active beta-10 Alfvenic heat-flux sensitivity cases with `abs(k_parallel)` decreased and increased by a factor of 100, but the guarded production matrix contained only the nominal active/passive definitions | staged `MKS24.tex:642-655`, `inputs/cgl_lf_paper/cgl_lf_paper_heat_flux_beta10_{strong,weak}.athinput`, `scripts/cgl_lf_workflow.py`, `scripts/frontier/cgl_lf_frontier.py` | High | Add a guarded `paper-heat-flux` workflow with the two nonnominal definitions, exclude those decks from the debug-only utility, and revise the proposed storage/node-hour envelope | Expanded production-input/workflow-case assertion and Frontier offline self-test | Implemented for defined Figure 12 execution scope; F-036 supplies reference extraction, while authorized execution and quantitative comparison remain open |
| F-036 | 2026-05-25 | Figure 12 reference/product audit | The top Figure 12 panel is peak alignment cosine versus `k_perp`, but the comparator admitted only selected-shell alignment distributions and no Figure 12 numeric curves existed | pinned `source/fig12.pdf` SHA-256 `aeeb79fef6a24de93d6fa44d2a6c62a254bcb65c2651dadaa5e95d2204c7001b`; `scripts/digitize_cgl_lf_mks24_fig12.py`; ignored `digitized_fig12_v1/curves.json` SHA-256 `bcdb17cbf2396c7e373c83fd0bef49a4b8a92498d740deba4f3bf7229e3260e0` | High | Add `alignment_peak.cos_theta`, forward and archive explicit alignment-shell selection through `paper-analyze`, and vector-extract both Figure 12 panels with panel-specific uncertainty rules and clipped-endpoint handling | Expanded exact-curve regression for alignment peaks; checksum-validated eight-curve manifest ingestion probe | Implemented provisionally; F-038 supersedes admission of the dimensionful lower spectrum while retaining the alignment product |
| F-037 | 2026-05-25 | Figure 13 reference/product audit | Figure 13(a)-(b) appeared to map directly to existing `spectra.velocity` and `pdf.beta_delta` products, while panels (c)-(d) displayed additional strain/transfer normalizations not exposed by the analyzer; the published panels are separate PDFs | pinned `source/fig13{a,b,c,d}.pdf` SHA-256 values retained in `digitized_fig13_v1/digitization_metadata.json`; `scripts/digitize_cgl_lf_mks24_fig13.py`; ignored `digitized_fig13_v1/curves.json` SHA-256 `410e208f968f6e3b1007baff9e97ae8aeb80b03f2e3039b819fae7dc2e9494eb` | High | Admit backward-compatible multi-source digitized provenance, vector-extract provisional panel-(a)/(b) paths, and exclude normalized panels until matching products exist | Expanded multi-source checksum rejection regression; checksum-qualified seven-curve manifest ingestion probe | Implemented provisionally; F-038 supersedes admission of dimensionful panel (a) while retaining dimensionless panel (b) |
| F-038 | 2026-05-25 | Paper/reference normalization audit | MKS24 states `p0 = 100` in code units for its beta-100 Figure 13 runs and plots panel (c) against `bhat bhat : grad(u) / <B^2>`, while AthenaK's internally consistent `B0 = 1`, `v_A = 1` beta-100 decks use `p0 = 50`; an absolute ordinate transform for dimensional panels is not established | pinned `source/MKS24.tex:467,687`, `source/fig13{a,c,d}.pdf`, `src/pgen/tests/cgl_lf_paper.cpp:167-174`, `scripts/digitize_cgl_lf_mks24_fig{12,13}.py`; ignored `digitized_fig12_v2/curves.json` SHA-256 `8bc302c082ed971eda657f5aac283fd423ab182dbd676faa7da72ca3436a0b9c` and `digitized_fig13_v2/curves.json` SHA-256 `816013b3e4af9f5755fde5d64228721c86f14206a69210e3b2c38db35c1aa12e` | Blocker | Do not alter normalized decks without donor-state evidence; emit only dimensionless Figure 12 alignment and Figure 13(b) `beta Delta` comparisons, and exclude dimensional spectra/strain/transfer until their transform and `T_total` are qualified | Regenerated checksum-bound `v2` manifests; extractor compilation/style and manifest-ingestion validation | Implemented for fail-closed admission of dimensional panels; F-039 subsequently admits dimensionless normalized-transfer ratios |
| F-039 | 2026-05-25 | Dimensionless transfer-reference audit | MKS24 defines `T_Delta p/T_total` with `T_total ~= E_K (2 pi u_rms/L_perp)` in its diagnostics section, but the analyzer retained only unnormalized transfer, preventing direct admission of the dimensionless Figure 7 lower panel and Figure 13(d) | pinned `source/MKS24.tex:488-494,591`, `source/fig7.pdf` SHA-256 `5411da17391d19fceb26f617f91b4fb1153c7638ccca8afee6c8e0d013651afa`, `source/fig13d.pdf` SHA-256 `bf1c063fa5a54c2dbbe9dce417490b5bf821fcd2e1cb1926d78c923b90387a24`, `scripts/analyze_cgl_lf_paper.py`; ignored `digitized_fig7_v1/curves.json` SHA-256 `7c2f02bca4725674b9e77292a7f33c40de211bdcd48fe2e29a4b2e8fe81c1d09` and `digitized_fig13_v3/curves.json` SHA-256 `78f084d55d851561bf7a649becbfa907626f0890605b5e792fe0352cc205c924` | High | Emit the MKS24 Kolmogorov-rate-normalized transfer product and vector-extract only the dimensionless transfer panels with absolute ordinate uncertainty | Synthetic exact normalized-transfer manifest regression; checksum-bound Figure 7/13 generation and manifest-ingestion validation | Implemented for comparison infrastructure/reference curves; authorized case execution and quantitative simulation-to-paper comparisons remain open |
| F-040 | 2026-05-25 | Figure 9 dimensionless-reference audit | The comparator already exposes `alignment_peak.cos_theta` and the standard workflow already defines the eight active/passive Alfvenic/random beta-10/beta-100 cases plotted in Figure 9, but no numeric Figure 9 reference curves were admitted | pinned `source/MKS24.tex:610-614`, `source/fig9.pdf` SHA-256 `9d0dbb05514a8ce7a48eddb178b05a47ec512e58d2a2ff99d899123f8e458984`; `scripts/digitize_cgl_lf_mks24_fig9.py`; ignored `digitized_fig9_v1/curves.json` SHA-256 `07417294ca44e7aa7e8c499931e12ed7fce9eec55e4691e3172a8684470163b8` | High | Vector-extract only the dimensionless peak-alignment paths and map them to the existing standard case definitions with absolute ordinate uncertainty | Checksum-bound Figure 9 generation and eight-curve manifest-ingestion validation | Implemented for dimensionless alignment references; authorized standard execution and quantitative comparison remain open |
| F-041 | 2026-05-25 | Figure 11 scale-separation audit | The dimensionless lower-panel `alignment_peak.cos_theta` product is already comparable, but the guarded matrix lacked the plotted `n_perp = 96` and `384` active Alfvenic beta-10 cases and retained no Figure 11 numeric references | pinned `source/MKS24.tex:628-634`, `source/fig11.pdf` SHA-256 `efab7ae64d5837ee621249acfedf79cd77bb221e37504f8c41a79bf5ba30f5f8`; `inputs/cgl_lf_paper/cgl_lf_paper_scale_separation_beta10_nperp{96,384}.athinput`; `scripts/digitize_cgl_lf_mks24_fig11.py`; ignored `digitized_fig11_v1/curves.json` SHA-256 `9922eb24245cdf0e7dfe355864a610ee8ebf437ccf593ea40a2e7df6016c037a` | High | Define a guarded scale-separation production workflow reusing the standard `n_perp = 192` case, exclude it from debug preparation, vector-extract only lower-panel alignment paths, and revise storage/runtime estimates for two additional unique runs | Production-input/workflow/authorization assertions; Frontier offline self-test; checksum-bound Figure 11 generation and three-curve manifest-ingestion validation | Implemented for definition/reference infrastructure; upper-panel dimensional spectrum, approved execution, and quantitative comparison remain open |
| F-042 | 2026-05-25 | Remaining panel-reference admission audit | No additional plotted panel can be honestly admitted through the current one-dimensional curve contract without new information: Figure 8 is dimensionless but stores two `alignment` heatmaps and its colorbar as RGB raster images, while Figure 4(a) labels raw `rho` rather than the analyzer's normalized `rho/<rho>-1`; remaining vector spectra/strain panels retain unresolved dimensional or product mappings | pinned `source/MKS24.tex:496,535-552,605-608`; `source/fig8.pdf` SHA-256 `0e3e292406bc27ce24b2ef449bd95cb8addf240f02b5ac60f60286459a299d8b` containing RGB images `1616 x 1332`, `1610 x 1332`, and colorbar `96 x 1332`; `source/fig4a.pdf` SHA-256 `18a5c7119f2cb316e0283b354d25479d006aa3f27df33031035b16810937c024` | High | Withhold Figure 8 surface/slice curves until a reviewed color-scale calibration or author numeric data establishes probability values and uncertainty; withhold raw-density/dimensional panels until explicit transforms and matching products are established | Pinned-source caption/stream/object inspection and SHA-256 audit; no comparison data emitted by design | Audited and documented; the remaining reference path now requires external data/calibration or additional reviewed product/transform work, plus authorized production comparisons |
| F-043 | 2026-05-25 | Figure 8 raster-calibration follow-up | F-042's Figure 8 withholding condition can be resolved locally: the pinned PDF embeds exact RGB heatmaps and a labeled linear colorbar, the caption states per-`k_perp` unit normalization, and the existing comparator already admits selected-shell `alignment.<shell>` PDFs | pinned `source/fig8.pdf` SHA-256 `0e3e292406bc27ce24b2ef449bd95cb8addf240f02b5ac60f60286459a299d8b`; `scripts/digitize_cgl_lf_mks24_fig8.py`; ignored `digitized_fig8_v1/curves.json` SHA-256 `a533778f7777d755a682e898ae44dcb7356be8f74a9818d6b6cc7b6fd6048522`; `digitization_metadata.json` SHA-256 `25933691ade9051bd37c9830405005437279171761d5974d5ece6c43e9bc38a3` | High | Supersede F-042 only for Figure 8 by decoding the calibrated colorbar, sampling shells `8,16,32,64,128`, requiring pre-correction unit-integral error below `0.03`, renormalizing each PDF, and retaining absolute density uncertainty `0.05` | Focused extractor compile/lint/run; analyzer-ingestion probe admits ten curves with exact-fixture residual zero and records maximum raw integral error `0.0264475132735984` | Implemented for Figure 8 reference scope; raw-density/dimensional transforms, authorized execution, and quantitative simulation comparisons remain open |
| F-044 | 2026-05-25 | Figure 4(a) normalized-density correction | F-042 incorrectly treated Figure 4(a) as raw density; the pinned vector PDF explicitly labels its horizontal coordinate `rho/<rho>`, which maps exactly to the existing `pdf.density_fluctuation = rho/<rho> - 1` product | pinned `source/MKS24.tex:539-548`; `source/fig4a.pdf` SHA-256 `18a5c7119f2cb316e0283b354d25479d006aa3f27df33031035b16810937c024`; `scripts/digitize_cgl_lf_mks24_fig4a.py`; ignored `digitized_fig4a_v1/curves.json` SHA-256 `00c0fadeb15a82758be543a95b45d883ebea32388961756d0146eef3f8ca4cb9`; `digitization_metadata.json` SHA-256 `c04f32ddc9c65fecd387bc3e7a81a58f384b69f6a406dbf4c2f142f1b22d25d0` | High | Supersede F-042 for Figure 4(a), vector-extract eight standard-case density PDFs through the exact abscissa shift, retain five-percent relative PDF uncertainty, and keep Figure 4(b)'s unqualified `E_rho` spectrum excluded | Focused extractor compile/lint/run; analyzer-ingestion probe admits eight `pdf.density_fluctuation` curves and 38 boundary-clipped vertices are recorded | Implemented for Figure 4(a) reference scope; F-047 subsequently adds the Figure 2(a) joint product, while unqualified raster/spectral panels, dimensional transforms, authorized execution, and quantitative comparisons remain open |
| F-045 | 2026-05-25 | Corrected remaining-panel boundary after F-043/F-044 | Figure 5(b) is dimensionless but was not representable through existing products: MKS24 defines characteristic parallel/perpendicular eddy lengths through local-field-conditioned structure functions, while the analyzer exposed PDFs, spectra, transfer, and alignment products only. Figure 2(a) remains a joint distribution; Figures 3, 4(b), 5(a), 6, 7 upper, 11 upper, 12 lower, and 13(a),(c) remain dimensional or otherwise unqualified; Figure 10 is a hybrid-kinetic contextual comparison outside direct CGL-LF reproduction. | pinned `source/MKS24.tex:478-496,509-512,530-586,621-650,670-680`; pre-F-046 `scripts/analyze_cgl_lf_paper.py`; extraction manifests through F-044 | High | Do not emit a Figure 5(b) manifest until a validated local-field-conditioned structure-function product and case comparison contract exist; retain F-038 dimensional exclusions and unmatched/out-of-scope panel boundaries | Pinned-source diagnostic/caption audit and analyzer product inventory; no reference data emitted by design | Superseded for Figure 5(b) by F-046; remaining dimensional/unmatched boundaries and production comparison gates remain open |
| F-046 | 2026-05-25 | Figure 5(b) conditioned-structure-function implementation | MKS24 Figure 5(b) plots dimensionless characteristic parallel/perpendicular eddy lengths for `u_perp` and `B_perp`, but no analyzer product or numeric reference manifest existed after F-045 | pinned `source/MKS24.tex:557-567`; cited Squire structure-function definition/`L_perp` normalization; `source/fig5b.pdf` SHA-256 `8886c425749f97af6f50b7960db071179758e6ea68b79ca547dfe1b163c2e0fe`; `scripts/analyze_cgl_lf_paper.py`; `scripts/digitize_cgl_lf_mks24_fig5b.py`; ignored `digitized_fig5b_v1/curves.json` SHA-256 `45365392e6ff8b173e829639d3409e5c2d4af05c261c9bf817e43cceab5fcde8`; `digitization_metadata.json` SHA-256 `506b08f970468bc51d42153701690d24304458184961067f1b39efed25d68064` | High | Add opt-in deterministic three-point structure functions conditioned within 15 degrees of the local three-point field, invert equal-power perpendicular/parallel scales, archive sampling controls in `paper-analyze`, and vector-extract the four normalized panel curves | Synthetic finite-eddy oracle; comparator product-selection regression; checksum-bound extraction and four-curve manifest-ingestion probe recording eight clipped endpoints | Implemented for Figure 5(b) analysis/reference scope; authorized standard execution and quantitative simulation-to-paper comparison remain open |
| F-047 | 2026-05-25 | Figure 2(a) joint-product completion and reference-boundary audit | MKS24 Figure 2(a) is a two-dimensional PDF in coordinates `<p> delta rho/<rho>` and `delta p_parallel`/`delta p_perp`, while the analyzer retained only marginal PDFs. Its embedded raster and pressure-scaled ordinate also cannot be quantitatively admitted through the one-dimensional curve interface or the unresolved donor pressure convention. | pinned `source/MKS24.tex:503-512`; `source/fig2a.pdf` SHA-256 `277545ff1a86b5ca46e97fa034ce189f72f149ecbb515470c30996494c1de129`; `scripts/analyze_cgl_lf_paper.py`; `scripts/plot_cgl_lf_paper.py` | High | Emit common-bin per-snapshot and ensemble joint PDFs in the displayed AthenaK coordinates, render the two pressure panels with the `5/3` guide, and retain the raster as excluded comparison context pending a qualified pressure/color-density transform | Synthetic finite joint-PDF and exact `5/3` coordinate oracles plus binary-snapshot product regression; rendered-analysis probe | Implemented for Figure 2(a) analysis-product scope; F-050 supplies surface comparison, F-051 decodes raw samples, and F-052 supplies admitted transformed surfaces; authorized execution and comparison remain open |
| F-048 | 2026-05-25 | Figures 3/4(b)/6(a) spectral-product completion | MKS24 requires `E_{khat dot u}` for Figure 3, density-fluctuation spectra for Figure 4(b), and separate `p_parallel`, `p_perp`, and magnetic-pressure spectra for Figure 6(a), but the analyzer previously exposed only total velocity, raw density, magnetic-field, and `Delta p` spectra. Figure 4(b) is vector-extractable, yet its raw `E_rho` ordinate is not qualified against an explicitly normalized AthenaK density product. | pinned `source/MKS24.tex:478-486,522-580`; `fig3a.pdf` SHA-256 `7d2cf6720f390c21f2710eafe86c921c4002859623cd18b6015e79e64afe5932`; `fig3b.pdf` SHA-256 `bce1acabfb0fd9a29dfac893fcf48b75549c34f6c18c516b760b3aefbad13eb0`; `fig4b.pdf` SHA-256 `67764db729df84e6b245b88ffd476a8910962d1fd9224727090615168e8c2e3c`; `fig6a.pdf` SHA-256 `f6afdc345e647dcd71b4982c43bea1b8809ce21f4d1c22ac40bf48bf6a126dff`; analyzer/plotter | High | Emit full-wavevector-projected `spectra.compressive_velocity`, `spectra.density_fluctuation`, and separate `spectra.p_parallel`, `.p_perp`, and `.magnetic_pressure` products with normalization metadata and rendering; retain raw spectral references excluded until their transforms and required Figure 3 sonic-correlation cases are qualified | Synthetic longitudinal/non-compressive Fourier projection oracle plus binary-snapshot product/metadata regression; rendered-analysis probe | Implemented for missing spectral analysis-product scope; F-049 supplies missing guarded case definitions, while reference admission, authorized execution, and quantitative comparisons remain open |
| F-049 | 2026-05-25 | Figure 3 compressive-case matrix audit | MKS24 Figure 3 requires active random beta-1 and active random beta-100 sonic-correlation runs absent from the guarded case matrix; F-048 supplies the plotted compressive field but not those simulation definitions. | pinned `source/MKS24.tex:522-536`; `fig3a.pdf` SHA-256 `7d2cf6720f390c21f2710eafe86c921c4002859623cd18b6015e79e64afe5932`; `fig3b.pdf` SHA-256 `bce1acabfb0fd9a29dfac893fcf48b75549c34f6c18c516b760b3aefbad13eb0`; `inputs/cgl_lf_paper/cgl_lf_paper_compressive_*.athinput`; `scripts/cgl_lf_workflow.py` | High | Add a guarded `paper-compressive` workflow reusing four standard cases and defining the two missing unique cases; exclude it from debug preparation and revise cost/storage scope | Production-input/workflow/authorization regression and Frontier offline self-test | Implemented for Figure 3 execution-definition scope; spectral transforms, approved execution, and quantitative comparison remain open |
| F-050 | 2026-05-25 | Figure 2(a) surface-comparison contract | F-047 emits the panel's two-dimensional joint PDFs, but the external reference manifest and renderer accepted only one-dimensional curves, so even qualified raster or author-provided surface data could not be compared. | pinned `source/MKS24.tex:503-512`; `source/fig2a.pdf` SHA-256 `277545ff1a86b5ca46e97fa034ce189f72f149ecbb515470c30996494c1de129`; `scripts/analyze_cgl_lf_paper.py`; `scripts/plot_cgl_lf_paper.py` | High | Extend the checksum/uncertainty-qualified manifest with optional sampled surfaces mapped only to `pressure_density_joint.parallel`/`.perpendicular`, bilinear interpolation, and residual rendering; do not admit untransformed raster data | Exact sampled-surface binary-snapshot manifest comparison and rendered-surface regression | Implemented for comparison-interface scope; F-051 decodes donor-coordinate raster surfaces and F-052 admits transformed surfaces; authorized execution and quantitative comparison remain open |
| F-051 | 2026-05-25 | Figure 2(a) raster-decoding and unit-boundary audit | F-050 supplies a surface comparison contract, and the pinned panel contains a numeric colorbar and axes, but its two plotted coordinates remain in donor pressure units while AthenaK has no qualified donor-pressure conversion. | pinned `source/fig2a.pdf` SHA-256 `277545ff1a86b5ca46e97fa034ce189f72f149ecbb515470c30996494c1de129`; `scripts/digitize_cgl_lf_mks24_fig2a.py`; ignored `digitized_fig2a_v1/excluded_surfaces.json` SHA-256 `3c2c4e8984d538acda37c6c36c49dfa11eae73a7ecadd03a05074b5660cbc428`; `digitization_metadata.json` SHA-256 `5354a7206f80a0630e9f86db1dd3a9c5efd552c59e2cc2a8341bf78ebd89f541` | High | Reconstruct the tiled raster, decode its 256-color labeled bar, mask only verified black overlays, and emit sampled surfaces in donor pressure coordinates with explicit uncertainty; withhold an analyzer manifest until conversion `x_A = s x`, `y_A = s y`, `z_A = z/s^2`, `sigma_z_A = sigma_z/s^2` is qualified | Checksum-bound extraction; palette/overlay invariant; `6953` parallel and `6952` perpendicular retained surface samples | Superseded for Figure 2(a) admission by F-052; retained as raw-raster provenance |
| F-052 | 2026-05-25 | Figure 2(a) pressure-transform admission | F-051 left donor-coordinate surfaces excluded, but the pinned source determines the beta-10 pressure mapping: MKS24 holds `B0` fixed while varying `p0` with `beta0`, and its beta-100 limiter footnote states `p0 = 100` in code units; thus its beta-10 Figure 2(a) run has `p0 = 10`, while the matched AthenaK beta-10 deck uses `p0 = 5`. | pinned `source/MKS24.tex` SHA-256 `6d6e748fd1883c5d33167be653d67aed2f84a9b364267d9e90ea075df184af4c` at lines `467,687`; `source/fig2a.pdf` SHA-256 `277545ff1a86b5ca46e97fa034ce189f72f149ecbb515470c30996494c1de129`; `inputs/cgl_lf_paper/cgl_lf_paper_standard_active_alfvenic_beta10.athinput`; ignored `digitized_fig2a_v2/surfaces.json` SHA-256 `3e2f961b1bdc8bb46180be1e4b38448a997d306738a834898caf3b5b0bd89845`; `digitization_metadata.json` SHA-256 `89575ae0dd3953a72c903fca377c3709923e93199b522fd2bab4178b0a862016` | High | Require both pinned source checksums; retain raw donor samples and emit analyzer-compatible surfaces with `s = 0.5`, `x_A = s x`, `y_A = s y`, `z_A = z/s^2`, and `sigma_z_A = sigma_z/s^2`; do not treat this pressure-only mapping as a qualified spectral transform | Checksum-bound extraction and transform metadata; `6953` parallel and `6952` perpendicular admitted surface samples; manifest-ingestion probe | Implemented for Figure 2(a) reference admission; authorized execution and quantitative simulation comparison remain open |
| F-053 | 2026-05-25; followed up 2026-05-26 | Remaining spectral/strain normalization evidence audit | MKS24 defines its shell-binned Fourier spectra but does not specify the discrete Fourier normalization behind absolute plotted ordinates. Its cited Squire et al. donor source additionally states spectra are plotted in `phi^2 L_perp` units and defines rate-of-strain fields, but likewise omits the discrete transform convention needed to map digitized absolute curves to AthenaK products. A targeted publication/archive and public-repository audit identified article/preprint records but no linked diagnostic source or machine-readable panel-data deposit. | pinned `source/MKS24.tex:478-496`; cited donor source `AthenaMagnetoimmutability_JPP.tex:976-1002`; Cambridge/official preprint-record and public-repository audit performed 2026-05-25; official journal-page text check plus exact arXiv-ID Zenodo and exact-title GitHub repository searches repeated 2026-05-26 | Blocker | Keep Figures 3, 4(b), 5(a), 6, Figure 7 upper, Figure 11 upper, Figure 12 lower, and Figure 13(a),(c) spectral/strain ordinates excluded unless author numeric data, donor diagnostic code, or an explicit discrete-normalization statement establishes a checksum/provenance-qualified mapping | Manual source-line audit and targeted public-source/data/repository search; the retrieved journal-page text contained no `Data availability` heading and the exact-query repository checks returned no matching data/code record on 2026-05-26, so no speculative curve manifest was emitted | Documented external-evidence boundary; the existing dimensionless and Figure 2(a) pressure-only admissions are unchanged, and authorized simulations/comparisons remain open |
| F-054 | 2026-05-26; superseded for continuation by F-057 | Stage I `R16` retained-output and bundle qualification | The first production segment used shared MPI-I/O retained snapshots/restarts and incurred an approximately twenty-minute output delay at its `t = 1` boundary. After switching to rank-local products, the first accepted-prefix bundling attempt also selected the obsolete clean diagnostic branch and incorrectly required sampled restart histories to repeat a boundary row. | Frontier jobs `4674731`, `4675322`, `4676557`, `4676696`, `4679803`, `4681476`, `4683291`, and `4683918`; retained segment inspections under `runs/mks24-stage-i/R16`; commits `1f52784d`, `476f9dbd`, `9c5c1b40`, and `bad8ba05`; diagnostic prefix bundles through `R16_rankio_t7p5_prefix_cadence_20260526` | High | Use `single_file_per_rank = true` for production snapshots/restarts, retain complete rank sets, select only the explicit inspected parent lineage for bundling, and permit only history-boundary gaps consistent with retained cadence and timestep records | Focused CPU regression `test_cgl_lf_stage_i_bundle_selects_terminal_restart_lineage`; formal segment inspections with zero strict failure counters through `t = 7.5`; prefix assembly containing 31 deduplicated rank-local snapshots; generic analysis retaining source provenance and correctly selecting zero pre-window samples while passing its synthetic check | Implemented for the archived pre-replacement retained-output path only; F-057 prevents continuing that lineage with the replacement driver. |
| F-055 | 2026-05-26 | Tier 3 mapped-case authorization consistency audit | The Tier 3 minimum-case table still required an active-Alfvenic beta-1 production run after the later frozen Stage I manifest excluded that separate definition as lacking a displayed-result mapping; the required published beta-1 Figure 3 role is the mapped active-random case `R10`. | `inputs/cgl_lf_paper/mks24_stage_i_manifest.json`; `docs/cgl_lf_weak_guide_manuscript_plan.md`; `MKS24.tex:522-536`; implementation-plan Tier 3 table | High | Replace the stale generic minimum table with the sixteen mapped canonical-case groups and state explicitly that the separate active-Alfvenic beta-1 definition requires a new published-panel mapping and revised authorization before execution | Manifest/source/document consistency audit; existing production-input regression excludes `standard_active_alfvenic_beta1` from Stage I while retaining `paper_compressive_active_random_beta1` | Implemented in the Tier 3 execution scope; no new production case was authorized, and F-057 records the completed-and-rejected continuation after the driver contract changed. |
| F-056 | 2026-05-29 | Post-merge diagnostic audit | `lf_qface` and heat-flux-cap counters counted every block-local LF face while `lf_qprwrk` and `lf_qpewrk` used single-owner face accounting, so equivalent MeshBlock tilings reported different cap fractions. | `src/diffusion/cgl_landau_fluid.cpp`; one-block versus two-block reduced paper smoke with an activated heat-flux cap | High | Apply the existing owned-face rule to `lf_qface`, `lf_qprcap`, `lf_qpr10`, `lf_qpecap`, and `lf_qpe10`, matching the work-ledger topology convention. | `test_cgl_lf_heat_flux_cap_fractions_are_meshblock_layout_independent`; focused CGL CPU and MPI/AMR reruns | Implemented locally; production bundles generated before this correction must not be used for cap-fraction comparison without regeneration. |
| F-057 | 2026-05-29 | Turbulent-driver replacement versus in-flight Stage I provenance | Job `4686032` completed after submission under immutable executable `e129a7ff8be6`, but the merged turbulent-driver replacement changes the forced-state representation: restart output now records modal OU coefficients with `TurbulenceRestartMetadata`, and restart input requires that new metadata and modal state. The existing ranked `R16` restart lineage therefore cannot be continued as evidence for the replacement executable. | merged changes in `src/srcterms/turb_driver.{hpp,cpp}`, `src/outputs/restart.cpp`, and `src/pgen/pgen.cpp`; job `4686032`; `s08_rankio_t7p5_t8p5/manifest/segment_inspection.json`; Stage I ledger updated 2026-05-29 | Blocker | Preserve pre-replacement segments as archival pipeline/timing evidence only, record `s08` consumed time without admitting it to a current reproduction lineage, and hold further production until the replacement driver has a qualified immutable build and fresh lineage. | Formal `s08` inspection reached `t = 8.5` with zero strict LF failure counters; accounting recorded `1.282778` node-hours as `rejected`, increasing Stage I actual use to `9.962778` with zero active reservation. | Production paused; the replacement driver requires new qualification and a fresh Stage I lineage from `t = 0` before any MKS24 comparison claim. |
| F-058 | 2026-05-30 | Shared Frontier-root operational audit | `/lustre/orion/ast207/proj-shared/dfielding/CGL` is not an MKS24-only root: exploratory `runs/beta25-accel05-*` campaigns coexist with Stage I archives. During the audit, job `4743106` was actively writing beneath that root. Slurm now records it as `COMPLETED 0:0`, but its top-level manifest remains stale at `state = running`. The prior Stage I helper's submission preflight queried `batch` but rejected only lines containing `cgl_mks24_`, so it did not prevent unintended overlap with other CGL-root campaigns. | `squeue -u "$USER"` on 2026-05-29 EDT; exploratory `prepared_run.json` manifests for cancelled job `4743020` and stale-record job `4743106`; `scripts/frontier/cgl_lf_stage_i.py:production_queue_output`, `shared_root_campaign_conflicts`, and `check_submit` | High | Keep exploratory campaigns out of MKS24 evidence ledgers, fail closed on any queued user job, inspect active top-level CGL-root campaign records, and require explicit acknowledgement only after an isolation review. Never modify or prune another campaign while preparing MKS24 work. | Read-only root inventory, queue audit, and live Slurm-state follow-up; focused offline regression `test_cgl_lf_stage_i_isolates_epoch_and_checks_all_shared_root_jobs`. | Implemented and committed as `462b9dbd`; explicitly acknowledge the reviewed stale `4743106` record and perform the live preflight again immediately before each submission. |
| F-059 | 2026-05-30 | Fresh Stage I execution-epoch boundary | The prior Stage I helper stored runs beneath `runs/mks24-stage-i/<case>` and used one mutable Stage I ledger. Starting a fresh replacement-driver `R16` lineage in that namespace would leave archival and current segments discoverable together and make recosting/accounting harder to audit, even though digest checks reject some invalid combinations. | archived `runs/mks24-stage-i/R16`; new `runs/mks24-stage-i/E02-modal-driver`; `accounting/mks24_stage_i_E02_modal_driver_*`; `scripts/frontier/cgl_lf_stage_i.py` layout and bundle discovery; F-056/F-057 | Blocker | Preserve historical work as `E01-pre-modal-driver`; use explicit `E02-modal-driver` run, reservation, summary, and ledger paths; reject cross-epoch restarts and bundle discovery; reset the future planning cap incrementally to `4000` node-hours while continuing to report historical use; permit only an initial `4.0` node-hour `R16` pilot until recosting review. | Focused offline epoch-isolation, queue-preflight, continuation-rejection, and bundle-lineage regressions. | Implemented and committed as `462b9dbd`. The full Stage I matrix reservation remains blocked pending replacement-driver qualification and pilot recosting. |
| F-060 | 2026-05-30 | Retained source-provenance and debug-preflight audit | The campaign-root working checkout was intentionally replaced by Git bundles on 2026-05-29, but both Frontier preparation helpers still required a source checkout beneath the root. The build script now correctly refuses such a checkout, so real preparation would fail or encourage restoring mutable campaign-root source state. The debug helper also copied only one selected restart file, which is insufficient for a rank-local MPI continuation because AthenaK rewrites the `rank_00000000` path for each rank, and its preflight still queried only `batch/debug` rather than every user job and active top-level root record. | `source-archives/README.md`; `scripts/frontier/build_cgl_lf_frontier.sh`; `scripts/frontier/cgl_lf_frontier.py`; `scripts/frontier/cgl_lf_stage_i.py`; `src/main.cpp:235-262` | High | Permit a clean external source checkout, require and checksum an in-root Git bundle containing each prepared revision, record bundle provenance in manifests, archive complete rank-local restart sets for debug continuation qualification, and apply the all-user-job plus reviewed-shared-root-record preflight to debug jobs as well as Stage I. | Debug helper self-test including bundle, rank-local restart-set, all-user-job, and shared-root-record checks; focused Stage I bundle-provenance regression; retained bundle `source-archives/athenak-feature-cgl-through-eab6e12b.bundle` with SHA-256 `dbbf139dfd41d50de1816858d0372ac15cb75e6b79a795c240406ba4ca530483`. | Implemented and committed as `d210cdd5` plus `eab6e12b`; live `g015` manifest retains the bundle checksum and immutable executable revision. |
| F-061 | 2026-05-30 | Replacement-driver G014/G015 reduced startup qualification | Debug job `4743463` (`g014`) reached the allocation but failed input parsing before physics startup: the archived MPI-I/O sizing deck predated `limiter_hardwall`, and AthenaK command-line overrides may change only parameters already declared in the deck. The helper accepted that invalid combination and spent an avoidable `0.002778` node-hours. | `runs/qualification-462b9dbd-e02/g014-modal-rankio-startup-8gpu`; Slurm job `4743463`; `scripts/frontier/cgl_lf_frontier.py`; corrected shared-root input `inputs/archived/cgl_lf_modal_rankio_startup_active_beta10_96x96x192.athinput` SHA-256 `f41e8a3fc98bc2d9881ad9281fa9af76c2ecb2043bec95a19f2dfda015d39ec5`; corrected G015 evidence JSON SHA-256 `e2f4a047f601e6d1539d77a85a09d65e4df8b338c22eb6638441cdd6af75cd70` | High | Parse block-qualified input keys during debug `prepare`, reject malformed or absent override targets before reserving node-hours, retain an explicit reduced modal-driver hard-wall/rank-local deck, and rerun the startup qualification. | Expanded debug helper self-test plus live negative prepare probe leave no run tree or reservation; corrected job `4743465` (`g015`) completes `t = 0.01` in 13 cycles with zero strict safety counters, finite readable eight-rank snapshots, fifteen complete eight-rank restart groups, and forcing-work relative residual `6.066e-12`. | Implemented and committed as `7fae0bcf`; retained current bundle `athenak-feature-cgl-through-7fae0bcf.bundle` has SHA-256 `15235ba07685388523d618c9c6d59df305b05d1ffbaded71795640ba954330c9`. Qualified for reduced startup scope only. |
| F-062 | 2026-05-30 | Replacement-driver G016-G018 rank-local modal restart identity qualification | `g016` proves that all eight archived `g015` restart siblings load and continue cleanly, but its intended identity comparator `g017` follows a different timestep sequence: stopping `g015` at exact `t = 0.01` clips the final step. Treating that divergence as a restart failure or pass would be incorrect. | `runs/qualification-462b9dbd-e02/g016-modal-rankio-resumed-t0p01-t0p02-8gpu`; `g017-modal-rankio-reference-t0-t0p02-8gpu`; `g018-modal-rankio-natural-restart-identity-8gpu`; jobs `4743470`, `4743472`, and `4743473`; `g018` evidence JSON SHA-256 `4d898cf4aa4427b81164e4eb688f011adbbb277874b1ca7453a7a2c1e60b580f` | High | Retain `g016` as restart-load smoke and `g017` as an inconclusive comparison-design attempt. For identity, branch from one natural-cycle `g017` checkpoint and compare a resumed continuation against the uninterrupted parent terminal state. Use `rtol = atol = 1e-12` for histories and turbulence-force tabs; use `rtol = 1e-6`, `atol = 1e-7` for retained float32 snapshot fields. | `g018` resumes complete `.00013.rst` rank set, reaches `t = 0.02`, retains zero strict safety counters and complete eight-rank products, and passes every evidence check. Maximum MHD-history, user-history, snapshot-field, and force-tab absolute differences are `6.10e-17`, `1.82e-13`, `3.73e-9`, and `2.04e-14`. | Reduced rank-local modal restart identity qualified. GPU/MPI decomposition, passive/random smoke, nonlinear hard-wall, and standard-layout rank-local I/O remain open. |
| F-063 | 2026-05-30 | Replacement-driver G019 one-rank versus eight-rank GPU/MPI decomposition qualification | Rank-local output filenames and shard boundaries differ by MPI decomposition, so comparing individual files would test storage layout rather than the physical solution. | `runs/qualification-462b9dbd-e02/g017-modal-rankio-reference-t0-t0p02-8gpu`; `g019-modal-rankio-reference-t0-t0p02-1gpu`; jobs `4743472` and `4743474`; `g019` evidence JSON SHA-256 `509ec34d4f1f27dfdbb9b507d2ca663d7233d168a0b83775b0ab532dc283c752` | High | Run the corrected reduced deck uninterrupted through `t = 0.02` with one GPU rank and compare it against eight-rank `g017`. Require the same timestep sequence, finite histories, zero strict LF failure counters, complete products for each layout, assembled-grid float32 snapshot agreement at `rtol = 1e-6`, `atol = 1e-7`, and history/turbulence-force agreement at `rtol = atol = 1e-12`. | `g019` completes the same 25-cycle sequence as `g017`; terminal MHD-history, user-history, assembled snapshot-field, and force-tab maximum absolute differences are `4.44e-17`, `7.11e-15`, `3.73e-9`, and `8.88e-16`. | Reduced one-rank versus eight-rank GPU/MPI decomposition qualified. Passive/random smoke, nonlinear hard-wall, and standard-layout rank-local I/O remain open. |
| F-064 | 2026-05-30 | Replacement-driver G020/G021 reduced forcing-family smoke | The modal-driver build had active-Alfvénic startup, restart, and decomposition evidence, but still needed replacement-driver coverage for the mapped passive branch and isotropic-random forcing orientation. | `runs/qualification-462b9dbd-e02/g020-modal-smoke-passive-alfvenic-1gpu`; `g021-modal-smoke-active-random-1gpu`; jobs `4743480` and `4743483`; archived debug-only input SHA-256 values `28f4a8bee361abb56976c5da3224033fb76b276f3303e8763d39de498e7f89ff` and `c0c1bec0c0bacc74b5e297201b4316942924afc40d14a4e34c3e3502636f469d`; evidence JSON SHA-256 values `4abeddeb25db489b3334be240b837e703d035585d0814fc79d004e273cc34b3f` and `aaee51f186ac0f223fa4a1c35b238bb1c2e22f2a1ac9075ad1e28ffea53e8498` | High | Rerun the established one-GPU reduced matrix with current tracked smoke physics copied into explicitly debug-only archived decks. Require zero strict counters and finite histories/tabs for both; require perpendicular-only force and zero passive pressure work for passive-Alfvénic `g020`; retain `driving_type = 0`, `rseed = 314159`, all three force components, and finite nonzero active pressure work for active-random `g021`. | Both runs reach `t = 0.01` cleanly. `g020` retains three checkpoints with `force_prl2 = 0`, `force_prp2 = 128`, and `lf_cpwrk = lf_cawrk = 0`. `g021` retains four checkpoints, all three force components, and finite nonzero `lf_cpwrk`/`lf_cawrk`. | Reduced replacement-driver passive-Alfvénic and active-random startup smoke qualified. Long-time statistics, nonlinear hard-wall, and standard-layout rank-local I/O remain open. |
| F-065 | 2026-05-30 | Replacement-driver G022 corrected nonlinear hard-wall qualification | Replacement-driver startup, modal restart, decomposition, and forcing-family smoke passed, but the merged driver still needed the corrected F-033 reduced nonlinear trajectory rerun before production recosting. | `runs/qualification-462b9dbd-e02/g022-modal-f033-hardwall-active-beta10-t2-96x96x192-8gpu`; job `4743484`; retained input SHA-256 `443b52f67a44cf49b74a248222b73db6e0de3dc62d2e79c450c26342763f2fc1`; evidence JSON SHA-256 `339a3309276a93be7bb55c196f367ead9754f7c2ebdf94afdb9aeedd37f026cf` | High | Rerun the corrected hard-wall active-beta10 `96x96x192` trajectory through `t = 2.0` with the immutable modal-driver build. Require zero strict counters, active algebraic projection, finite forcing/LF/pressure-work ledgers, forcing-work residual below `1e-10`, readable snapshots, retained checkpoints, and MPI-I/O records without production extrapolation. | Job `4743484` completes in `1554` allocated seconds and uses `0.431667` node-hours. Terminal `lf_hwproj = 10055979902`; strict counters remain zero; forcing-work relative residual is `1.863e-12`; nine snapshots, five checkpoints, and fourteen net-write records are retained. The terminal modal restart is `227140784` bytes. | Reduced replacement-driver nonlinear hard-wall correctness qualified. Standard-layout rank-local I/O sizing and long-time production statistics remain open. |
| F-066 | 2026-05-30 | Replacement-driver G023 standard-layout rank-local startup and retained-output sizing qualification | Reduced replacement-driver correctness passed through G022, but the production handoff still lacked a current exact-standard-layout launch using the intended rank-local snapshot and modal-restart protocol. Historical G011 measured shared MPI-I/O before the modal-driver replacement and cannot establish current retained-product sizing. | `runs/qualification-462b9dbd-e02/g023-modal-rankio-sizing-active-beta10-192x192x384-8gpu`; job `4743662`; retained input SHA-256 `cee4c1c97ceb1be2f5d72488c15ca4a02f8f95e5402c58343b461e85c77e5c96`; evidence JSON SHA-256 `effba10b616246e1c9f6863c6a068f674551fa9ac4b4a02bdb92ead0f963c14f` | High | Launch the guarded `192 x 192 x 384` layout with one node/eight GPU ranks and intended rank-local binary/restart settings. Require startup through `t = 0.01`, zero strict counters, finite ledgers, forcing-work residual below `1e-10`, complete initial/final eight-rank snapshot and modal-restart groups, and readable assembled endpoints without extrapolating production runtime. | Job `4743662` completes in `79` allocated seconds and uses `0.021944` node-hours. It retains two complete binary groups with terminal aggregate size `509682248` bytes and two complete modal checkpoint groups with aggregate size `1381366136` bytes. The assembled endpoint reader passes at shape `384 x 192 x 192`; forcing-work relative residual is `1.053e-11`; strict counters remain zero. | Replacement-driver Frontier debug startup and retained-output sizing qualification closed. F-068 supplies initial rank-local R16 runtime and storage measurements; long-time production statistics remain open. |
| F-067 | 2026-05-30 | Post-merge focused local serial and scheduled MPI CPU retention | The `E02` entry gate explicitly requires retained current local evidence after the modal-driver merge and F-056 correction. Frontend direct invocation lacks `mpirun`, so the MPI pytest files require a scheduled wrapper rather than an unrecorded login-node workaround. | `runs/local-validation-462b9dbd-e02`; serial evidence JSON SHA-256 `3fc869c0d7455f1221d7838de900c0db7426673b2730cd954f264b8934fc0c29`; scheduled MPI CPU job `4743666`; MPI evidence JSON SHA-256 `397aced4e90aa7e77ab17da05a522ad4c99c9d0e15bd48d86ce20aa195fb36d3` | High | Run the focused serial CGL suites against the current source. Build a current Kokkos-Serial/MPI executable, schedule one isolated one-node debug allocation, translate test `mpirun -np` invocations into nested `srun --overlap` steps, and retain executable, wrapper, batch, log, source, test-count, and Slurm accounting provenance. | Serial suites pass `39` tests in `23.89` seconds. Scheduled job `4743666` passes both one-rank/four-rank MPI decomposition regressions in `2.81` pytest seconds, completes `0:0` in `8` allocated seconds, and uses `0.002222` node-hours separately from the debug-helper ledger. | Focused post-merge local serial/MPI entry gate closed. |
| F-068 | 2026-05-30 | Fresh E02 R16 segmented pilot, restart continuation, and preliminary completion recost | Replacement-driver debug and local gates closed, but the initial `4.0` node-hour production-path cap was intentionally only a pilot authorization. A current rank-local production measurement and one authenticated modal-restart continuation were required before extending even the first case. | `runs/mks24-stage-i/E02-modal-driver/R16/s00_rankio_t0_t2`; `s01_rankio_t2_t3p5`; jobs `4743735` and `4743933`; accepted-prefix bundle `R16_rankio_t3p5_prefix_20260530`; recost evidence JSON SHA-256 `5095d8833da4ca65eac0552d8dcbbfea87f696f6fd3932e25277dcb89c5534a5` | High | Run fresh `R16` from `t = 0` with rank-local products; inspect and record each segment before continuation; archive all terminal restart siblings; assemble transient-prefix bundles; require zero strict counters, finite ledgers, clean continuation, and synthetic analyzer pass. Use the developed continuation cost only to size an R16-completion envelope, not to authorize the mapped matrix. | `s00` reaches exact `t = 2.0` in `0.995278` node-hours. Authenticated restart continuation `s01` reaches exact `t = 3.5` in `0.870833` node-hours. The prefix retains fifteen snapshot groups, uses `1.866111` node-hours total, and stores `1820024680` raw snapshot/restart bytes. The transient `t = 0`--`2` forcing-work residual is `1.803e-12`; default `t = 8`--`10` analysis correctly selects zero prefix samples while its synthetic check passes. Developed-state extrapolation projects `5.639721` node-hours through `t = 10`. | Expand only the `R16` completion envelope to `7.0` node-hours, leaving `1.360279` projected margin. Continue sequential inspection through `t = 10`; do not reserve the remaining matrix until completed R16 production-window analysis and final recost review. |
| F-069 | 2026-05-30 | Fresh E02 R16 completion, production-window analysis, and standard-layout timing-pilot authorization | F-068 safely bounded R16 completion, but low-resolution R16 cannot determine standard-layout cost from cell count alone: the `192 x 192 x 384` cases also reduce CFL timestep, and `R17` adds another scaling risk. | accepted R16 jobs `4743735`, `4743933`, `4743988`, `4744019`, `4744056`, `4744120`, and `4744158`; completed bundle `R16_rankio_t10_20260530` manifest SHA-256 `9955085f3557a15c9775722bd1574dbd3636adbbb8df4c4ba5e4e3b9b37d6aa8`; production diagnostics SHA-256 `b3547ff697391dc7d56588028da2d7dd595e9b1ad33fbf7887e0e98188f9e6ab`; completion-recost evidence JSON SHA-256 `d42a0432858127a94ab4aee5aab0c0b79572b70730085352aec53147670adad3` | High | Inspect and record every continuation through exact `t = 10`; assemble the complete ranked lineage; analyze `t = 8`--`10`; calculate actual logical and retained storage; bracket matrix runtime with both cell-count-only and CFL-aware scaling; authorize only the smallest standard-layout timing measurement needed to collapse that uncertainty. | R16 completes in `6.145556` node-hours inside its `7.0` envelope and retains `4686854064` raw segment-tree bytes. The analyzer selects nine snapshots and 101 history rows, passes its synthetic test, and closes forcing work to relative residual `8.806e-13`. Matrix runtime brackets `1087.763412` to `2956.012436` node-hours; cell-scaled logical raw storage is `798538119152` bytes before margin. | Preserve completed R16. Authorize only a `2.0` node-hour `R02` standard-layout timing pilot next; do not authorize the remaining matrix until its measured review. |
| F-070 | 2026-05-30 | Fresh E02 R02 standard-layout timing pilot and bounded R17 high-resolution-probe authorization | F-069 isolated standard-layout runtime uncertainty, but a measured standard result still cannot choose between the `8x` cell-only and `16x` CFL-aware R02-to-R17 node-hour factors. | `runs/mks24-stage-i/E02-modal-driver/R02/s00_rankio_t0_t0p1`; `s01_rankio_t0p1_t0p25`; jobs `4744198` and `4744205`; R02 recost evidence JSON SHA-256 `81c9f6e56cca9eb5a69449a1166bf5d3e3cc8a0a4102b6b4fd0c3f73f3c2123b` | High | Run R02 from `t = 0` through `t = 0.1`, inspect and record the complete eight-rank modal-restart set, then continue from that authenticated restart through native snapshot boundary `t = 0.25`. Require zero strict counters and complete ranked products. Use only the developed continuation rate for the standard-case estimate; keep R17 bracketed until an eight-node/64-rank high-resolution probe preserves the standard 27-meshblock-per-rank load. | `s00` reaches exact `t = 0.1` in `0.188333` node-hours. Authenticated continuation `s01` reaches exact `t = 0.25` in `0.285000` node-hours with zero strict counters and complete products. Its `0.15` simulated-time interval projects `19.000000` node-hours per standard `t = 10` case. Measured terminal groups are `509682424` snapshot bytes and `1381366304` restart bytes. The mapped-matrix bracket narrows to `424.145556`--`576.145556` node-hours; cell-scaled logical raw storage is `798538318560` bytes before margin. | Preserve completed R02 timing evidence. Authorize only an `8.0` node-hour `R17` high-resolution timing probe next; do not authorize the remaining matrix until its measured review. |
| F-071 | 2026-05-30 | Fresh E02 R17 high-resolution timing pilot and measured mapped-matrix authorization | F-070 left only the R17 scaling factor unresolved. A fresh R17 segment alone still includes startup effects, so an authenticated continuation was required before replacing the bracket with a measured matrix envelope. | `runs/mks24-stage-i/E02-modal-driver/R17/s00_rankio_t0_t0p05`; `s01_rankio_t0p05_t0p1`; jobs `4744210` and `4744230`; R17 recost evidence JSON SHA-256 `79147715729a2c9017400e170de2191a3f7d6b0bf411acee09a948ada464b822` | High | Run R17 from `t = 0` through exact `t = 0.05` on eight nodes/64 GPU ranks, inspect and record the complete restart set, then continue through exact `t = 0.10`. Require 27 meshblocks per rank, zero strict counters, complete ranked products, and a developed high-resolution rate. Recalculate logical raw, acceptance-retained, and sequential-peak storage from measured ranked groups. | Both segments pass formal inspection with zero strict counters and complete products. The developed `t = 0.05`--`0.10` interval uses `2.124444` node-hours, projecting `424.888889` node-hours for full R17. Reusing completed prefixes yields a `697.019444` node-hour mapped E02 projection. Measured logical raw matrix storage is `798559758560` bytes; sequential peak with 20-percent margin is `747543662189` bytes. | Authorize only the frozen `R02`--`R17` mapped matrix under sequential inspection inside a `900.000000` node-hour E02 envelope. Resume R02 first and complete R17 last. |
| F-072 | 2026-05-30 | First mapped R02 production continuation and normal-QOS walltime guard | F-071 authorizes sequential mapped execution, but the first `02:15:00` R02 continuation preparation reached Slurm `--test-only` before revealing Frontier normal QOS rejects requests above 120 minutes. The helper should reject that limit before archiving restart siblings. | released `R02/s02_rankio_t0p25_t1`; accepted `R02/s03_rankio_t0p25_t1`; job `4744249`; R02 continuation evidence JSON SHA-256 `7ba577644cd65eefb4c9e3a782f284c36be9e86894af1ec17e6242eba1c17b06` | High | Preserve the released unsubmitted manifest, reprepare at `02:00:00` with Athena `01:50:00`, run from authenticated `t = 0.25` restart siblings through native restart boundary `t = 1.0`, require zero strict counters and complete ranked products, and encode the observed normal-QOS walltime cap during `prepare`. | Slurm rejects only the unsubmitted `s02` request. Job `4744249` completes `s03` at exact `t = 1.0` in `1.452778` node-hours with zero strict counters, expected active hard-wall projection, and complete products. The longer developed rate is `1.937037` node-hours per simulated unit, updating the mapped projection to `702.195370` node-hours and leaving `197.804630` margin. | Continue R02 conservatively from `t = 1.0` to `1.5`, inspect, and recost before lengthening later segments. Reject normal-QOS preparations above two hours. |
| F-073 | 2026-05-30 | Second mapped R02 production continuation and late-time recost | F-072 deliberately limited the next standard-layout interval to one half-unit so a late-time developed rate could be measured before longer production segments. | accepted `R02/s04_rankio_t1_t1p5`; job `4744518`; R02 continuation evidence JSON SHA-256 `0319284b72628f2456087f378588515726f0fa8586b20660d37da2a20a1ba6ba` | High | Run from authenticated `t = 1.0` restart siblings through native restart boundary `t = 1.5`, require zero strict counters and complete ranked products, and apply its measured rate conservatively to unfinished standard-layout work before authorizing another segment. | Job `4744518` completes exact `t = 1.5` in `1.052222` node-hours with zero strict counters, expected active hard-wall projection, and complete products. The late-time rate is `2.104444` node-hours per simulated unit, updating the mapped projection to `725.464938` node-hours and leaving `174.535062` margin. | Continue R02 conservatively from `t = 1.5` to `2.0`, inspect, and recost before lengthening later segments. |
| F-074 | 2026-05-30 | Third mapped R02 production continuation and late-time recost | F-073 retained a one-half-unit standard-layout interval so measured rate drift could be reviewed before lengthening later production segments. | accepted `R02/s05_rankio_t1p5_t2`; job `4744913`; R02 continuation evidence JSON SHA-256 `d0c7e0f9cdd8031a9d88c195adee9822943a358eb55193b0e1146bf3b95a08e4` | High | Run from authenticated `t = 1.5` restart siblings through native restart boundary `t = 2.0`, require zero strict counters and complete ranked products, and apply its measured rate conservatively to unfinished standard-layout work before authorizing another segment. | Job `4744913` completes exact `t = 2.0` in `1.068889` node-hours with zero strict counters, expected active hard-wall projection, and complete products. The latest rate is `2.137778` node-hours per simulated unit, updating the mapped projection to `730.081697` node-hours and leaving `169.918303` margin. | Continue R02 conservatively from `t = 2.0` to `2.5`, inspect, and recost before lengthening later segments. |
| F-075 | 2026-05-30 | Fourth mapped R02 production continuation and stable late-time recost | F-074 retained one more half-unit standard-layout interval so the latest late-time drift could be reviewed before authorizing sustained matrix execution. | accepted `R02/s06_rankio_t2_t2p5`; job `4745305`; R02 continuation evidence JSON SHA-256 `f823ef6993045d5d90f016f7052518a353f348f7f83b9efdc19d32ee097f0c59` | High | Run from authenticated `t = 2.0` restart siblings through native restart boundary `t = 2.5`, require zero strict counters and complete ranked products, and apply its measured rate conservatively to unfinished standard-layout work before authorizing sustained sequential execution. | Job `4745305` completes exact `t = 2.5` in `1.061944` node-hours with zero strict counters, expected active hard-wall projection, and complete products. The latest rate is `2.123888` node-hours per simulated unit, updating the mapped projection to `728.164877` node-hours and leaving `171.835123` margin. | Continue only the frozen matrix one inspected segment at a time and complete R17 last. |
| F-076 | 2026-05-30 | Phase A forcing-contract stop-the-line audit | The tracked paper decks describe MKS24 `k^-2` forcing and the production campaign treated the modal-driver qualification as sufficient, but an independent source-to-code review found that the E02 driver does not exactly implement the published forcing semantics. Random paper roles omit `spectrum = power_law` and inherit the generic projected-force policy; planar roles apply a full-`abs(k)^2` projection and therefore do not generally enforce published `grad_perp dot u_perp = 0` when retained `k_z != 0` modes are present. | pinned `MKS24.tex:469`; `src/srcterms/turb_driver.cpp:126-132,692-767`; E02 executable revision `462b9dbd53e085dea46c2478b567781576d7d03e`; cancelled `R02/s07_rankio_t2p5_t3` job `4745498`; stop-line evidence JSON SHA-256 `de85c443c3bc9bab795f7f1903e5de224d7d86b553ac40e9dd7249986c3df211`; `docs/cgl_lf_mks24_stage_i_protocol_review.md` | Blocker | Cancel active E02 work, account the allocation, preserve existing E02 output only as pipeline and cost evidence, fail closed on further E02 preparation, add explicit restart-retained MKS24 random and Alfvenic forcing policies while preserving generic driver defaults, make paper deck spectrum selection explicit, requalify, and start a fresh epoch from `t = 0`. | Historical stop-line snapshot: job `4745498` is cancelled after `498` seconds and recorded `aborted`, consuming `0.138333` node-hours; E02 cumulative use is `15.628610` with no active reservation. F-078 later closes source correction and replacement qualification. | E02 remains permanently fail closed: no additional E02 preparation or submission is permitted. |
| F-077 | 2026-05-30 | Corrected `E03-forcing-policy` worktree and fail-closed qualification boundary | F-076 requires a clean epoch because E02 physics is not exact. A fresh executable also needs stronger operational provenance than a prose-only requalification instruction. | `src/srcterms/turb_driver.cpp`; `src/srcterms/turb_driver.hpp`; `src/outputs/restart.cpp`; all tracked paper decks; `scripts/frontier/cgl_lf_stage_i.py`; `scripts/analyze_cgl_lf_paper.py`; `scripts/cgl_lf_workflow.py`; `inputs/cgl_lf_paper/mks24_stage_i_manifest.json`; focused CPU and local-root regressions; independent replay/recovery exploit probes; staged full-field policy-smoke input SHA-256 `33a712e85dc2bf522daa56bd7d9e08bb0a32f4e2a2a6d2b49d1e2918922c2a0a` and `0c8bd9f66981ecf6852135db999a0cb2192d866ce128ec3c68944b1045905fda` | Blocker | Preserve generic turbulence defaults for non-paper users; add explicit restart-retained `mks24_random_unprojected` and `mks24_alfvenic_perpendicular` policies; enforce planar `sol_fraction = 1`; add explicit `spectrum = power_law`; retain authoritative `time/restart_time`; isolate E03 accounting; require a reviewed immutable-build qualification token; atomically submit under lock; authenticate prepared artifacts across submission, runtime, record, reconciliation, and bundling; retain the schema-v2 panel inventory with checksum-bound references and pending scientific criteria. | Historical pre-F-078 snapshot: corrected source, decks, E03 helper, panel inventory, reference bindings, and local regressions were implemented in the worktree. The independent helper audit closed targeted rejection probes for the prior scheduler-timestamp, replay-baseline, reservation-preservation, and cleared-submit-evidence exploits. Deterministic replay passed before and after reservation replacement. At that boundary the canonical E03 summary reported zero actual use, zero active reservations, and a pending qualification token. Staged one-rank policy-smoke decks retained full force fields for direct HIP FFT checks. | Closed by F-078 after clean committed bundle archival, immutable HIP/MPI build retention, reviewed Frontier qualification, and approval-token creation. |
| F-078 | 2026-05-30 | Corrected E03 Frontier qualification closure and approval-token transition | F-077 intentionally fails canonical preparation closed until one immutable forcing-correct build passes the smallest discriminating Frontier ladder and its evidence is reviewed. | corrected revision `9e07542281e4e6d125582f253df3ad2e3b8b154d`; executable SHA-256 `68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c`; source-bundle SHA-256 `c39d55809989d20aa5438711803f4fd43237fa9284c7e15183e84fdd693d3687`; evidence JSON SHA-256 values `g024=7e0edfc2e45c141a58df57c15dfda77c937dec97a52cec23c4987b0efadd5916`, `g025=6479fee28f36872e9113e8df487405c4dadc1969df5eace20ffac2ed2b791912`, `g026=31bd80bc7f07a30e5ac40cd67443830f3ec3a0a93bd66ce7a5da49145aef394d`, `g027=ca73511829c2ce6389053295a4febde057e507da09245c0b7791f56d975219b4`, `g028=76ce1796332d3642dbe9610733fd6303c7acfb21d2e8eabe1cffd39a799ea2fc`, `g029c=1702a5f298f2874771616c6128811e0a5d11681c9f9570c1d941201bf95d8fec`, `g030b=593b0c85fb30a5dbbd44c39be71b883b6831129590c6444c988629dbbd579b05`, canonical `g031=a16415c5f9c557a0dad35936c0bb0e89658da9e3b78a4d86852ead4e006f0d98`, supplementary independent `g031=e000a789d10e749b476d09a1be5e6526f0f0dc2dbe0e46c47b0ef09bc0f56ff6`; approval token SHA-256 `e3fec9f35da42121b902f41ef752021f375aab38bc34c5ff75ae8780bfbca635` | High | Run explicit-policy smokes, uninterrupted/restarted/decomposition comparisons, passive semantics, reduced nonlinear hard wall, and exact-standard-layout ranked sizing against one immutable executable. Retain checksum-addressed evidence and require an independent review before creating the E03 approval token. | Jobs `4745621`, `4745643`, `4745651`, `4745756`, `4745825`, `4745843`, `4745848`, and `4745890` pass as `g024` through `g031`. `g027` crosses an OU refresh and matches the uninterrupted reference; `g030b` reaches exact `t = 2.0` with zero strict counters, terminal `lf_hwproj = 10084905222`, and forcing-work relative residual `1.8681072370071585e-12`; `g031` reaches exact `t = 0.01` with 27 meshblocks per GPU rank and complete ranked products. Corrected E03 debug use is `0.485834` node-hours; the full debug ledger is `1.847784`. Token retention and reconciliation passed; fresh `R02/s00_rankio_t0_t0p1` job `4745922` is accepted through exact `t = 0.1` for `0.188889` node-hours. | F-079/F-080 is archived; prepare the authenticated R02 continuation and keep R17 last. |
| F-079 | 2026-05-30 | Preserve reviewed shared-root acknowledgements in the printed atomic-submit follow-up | The first R02 `check-submit` preview enforced the reviewed stale exploratory-record acknowledgement but printed a follow-up `submit` command without that flag. Submitting the printed command verbatim would fail closed, but the omission is an operator usability defect. | `scripts/frontier/cgl_lf_stage_i.py:3037-3055`; atomic R02 submission job `4745922`; parser-round-trip regression in `tst/test_suite/cgl/test_cgl_landau_fluid_cpu.py`; isolated formatter probe | Low | Render every reviewed `--allow-shared-root-campaign` value in deterministic deduplicated equals style, shell-quote the complete option token, and round-trip unsorted duplicates, whitespace, and option-like values through the real parser. | Job `4745922` retained the required acknowledgement and is valid. The preview formatter and regression are corrected before any later submission. | Closed by F-080; prepare the authenticated R02 continuation. |
| F-080 | 2026-05-30 | Authenticate recorded historical helper bytes from each segment's retained immutable source bundle | Applying F-079 after recording R02 exposed that a single mutable live-helper checksum cannot authenticate both an already-recorded segment and a future helper revision. Repeatedly swapping live helper bytes is not an admissible lifecycle. | `scripts/frontier/cgl_lf_stage_i.py:1742-1854,2127-2248`; historical-bundle and route-level state-gating regression plus Stage I helper subset in `tst/test_suite/cgl/test_cgl_landau_fluid_cpu.py`; canonical `reconcile`; direct canonical recorded-parent, submitted-drift-rejection, accepted-lineage, and uncommitted-preparation probes; retained F-080 evidence JSON SHA-256 `46fc1c4054e75f4224302be1895c1b9eaae88097544512ac067c53376e542e34` | High | Keep `prepared` and `submitted` manifests strict against the live helper path and checksum. For `recorded` manifests only, authenticate the historical helper blob by revision and canonical repository-relative path from the checksum-bound retained source bundle; retain the normalized batch-script checksum while avoiding regeneration from future helper code. Require the retained bundle to contain the historical utility revision and fail closed on bundle tampering. | Closed by final controller transition `05cb4c324bfd8feec72ebdeb33b1961c9fde70bf`, archived as `athenak-feature-cgl-through-05cb4c324.bundle` with SHA-256 `1381918e471730d8c9639014475566f93fc87bac66b62738dd8316fc69c03570`. The offline Stage I helper subset passes (`10 passed`); post-archive canonical reconciliation is clean with one accepted R02 ledger row and no active reservation; direct recorded-parent, submitted-drift-rejection, accepted-lineage, and uncommitted-preparation probes pass. | Prepare the authenticated R02 continuation and keep R17 last. |
| F-081 | 2026-05-30 | Fresh corrected E03 R02 native-boundary continuation and bounded recost | F-080 closes the controller transition, but sustained corrected production still requires one authenticated developed-state timing interval before extending R02. | accepted `R02/s01_rankio_t0p1_t0p25`; job `4746154`; record-time inspection SHA-256 `f0c4229d9a2cd86bed76e87efd32204cc504fbd2d8585d0f5d4468bd25328d6f`; recorded manifest SHA-256 `d836755911182b554fc3e49c729d7f0e5e1e08a90d180816ae0cd752d81c122b`; retained R02 recost evidence JSON SHA-256 `eb6071cf53d453b2f75707003616ac9cbdfe025d62969d00e34e7948bee3310c` | High | Continue from the authenticated E03 `t = 0.1` restart siblings through native boundary `t = 0.25`; require exact terminal time, complete eight-rank snapshot and restart groups, zero strict LF failure counters, finite histories, forcing-work review, accounting, and reconciliation before extension. | Job `4746154` reaches exact `t = 0.25` in `0.289167` node-hours. Strict LF failure counters remain zero; terminal `lf_hwproj = 28701760`; complete eight-rank terminal products are retained; sampled-history forcing-work relative residual is `1.7741014899016423e-11`; reconciliation closes with two accepted ledger rows and no active reservation. The developed rate is `1.927778` node-hours per simulated unit, provisionally projecting `700.868889` matrix node-hours and `199.131111` margin inside the `900.000000` envelope. | Historical authorization was only `R02/s02_rankio_t0p25_t1`; that segment is now accepted under F-085. Keep R17 last. |
| F-082 | 2026-05-30 | Require committed live helper state during historical authentication and reserve the case-level analysis namespace | Writing the F-081 retained recost beneath `R02/analysis/` exposed that reconciliation treated a legitimate case-level analysis directory as an orphaned segment. The same review tightened F-080: historical-bundle authentication must never mask an uncommitted live helper. | commits `9480e62764528a3f40066d22a192f0e99b369891` and `ef1e42fa088203ac9ef6ec8e47e668db4fb95a3c`; `scripts/frontier/cgl_lf_stage_i.py`; Stage I helper subset in `tst/test_suite/cgl/test_cgl_landau_fluid_cpu.py`; launch bundle `athenak-feature-cgl-through-ef1e42fa.bundle` SHA-256 `80c211c41a8de32688c3be577d8c727ec3268c347c2a7c3b2d5143e1c34593ec` | High | Require the live helper to resolve to committed repository bytes before either strict-current or recorded-historical utility authentication. Exclude only reserved case-level `analysis/` from orphaned-segment discovery and reject `--segment analysis` so the namespace cannot be captured by a run. Re-run focused Stage I tests and canonical reconciliation before preparing another segment. | Focused Stage I tests pass (`10 passed`); canonical reconciliation is clean; job `4746182` prepares and submits from the verified launch bundle with production utility revision `ef1e42fa088203ac9ef6ec8e47e668db4fb95a3c` and SHA-256 `957f536c107ff9042ec38ed8fc407dda5a09eecb015b3e60d6e2301ac553f4b4`. | Preserve committed live-helper enforcement and the reserved analysis namespace at every later lifecycle transition. |
| F-083 | 2026-05-30 | Render follow-up lifecycle commands with the canonical absolute helper path | A relative helper invocation can select the wrong checkout when operators follow a retained command from another working directory. | commits `e7008b410`, `0ec89bcd9`, `e3bfa966a`, and `57cb7ca47`; promoted canonical revision `9689c269bf329542815a1b2b137881126964b05c`; `scripts/frontier/cgl_lf_stage_i.py`; public runbook | High | Render and validate follow-up lifecycle commands through the committed canonical absolute helper path while preserving checksum-bound historical helper authentication. | Promoted after accepted job `4746356` was recorded. The helper subset passes (`10 passed, 35 deselected`), Sphinx warnings-as-errors passes, and hardened reconciliation is clean. | Keep absolute-path command rendering in every later lifecycle transition. |
| F-084 | 2026-05-30 | Bind canonical case layout, R17 ordering, scheduler evidence, and replay arithmetic in the production controller | The retained protocol requires one-node standard cases, an eight-node R17, R17-last execution, and cumulative budget controls. Those invariants must be machine-enforced rather than prose-only. | promoted canonical revision `9689c269bf329542815a1b2b137881126964b05c`; helper SHA-256 `1c633ebb58294938a0a0609742ae8f8d1f88cd15242ffe649578796edeb39375`; `scripts/frontier/cgl_lf_stage_i.py`; Stage I helper regressions | Blocker | Require one node for R02--R16 and eight nodes for R17; require accepted `t = 10` R02--R16 lineages before R17; permanently block lower-case continuation after R17 starts; replay reservation, manifest, scheduler, ledger, cumulative-ceiling, and positive finite threshold bindings during lifecycle checks. | Promoted after accepted job `4746356` was recorded. Focused helper tests pass (`10 passed, 35 deselected`); canonical syntax and diff checks pass; hardened reconciliation reports four recorded manifests, no active reservation, and no issue. | Preserve F-084 enforcement for the remainder of E03. Execute R17 last. |
| F-085 | 2026-05-30 | Fresh corrected E03 R02 continuation through exact `t = 1.0` and longer-interval recost | F-081 bounds one longer R02 continuation so corrected developed-state timing can replace the short `t = 0.1`--`0.25` estimate before later segments. | accepted `R02/s02_rankio_t0p25_t1`; job `4746182`; record-time inspection SHA-256 `482c7f0f9552e319cc50e92f6d84f17565a90f9562b06f1d67ef1c57d240a2dd`; recorded manifest SHA-256 `73aa132eaff0d5bf92610c2ab153e1012d0b0695acbd4a9c572e97f9f13c8fed`; retained R02 recost evidence JSON SHA-256 `7d11e8a0004e24417167aeb1f186f8b1f68c8eb9d6206467a708a7cad16cf252` | High | Continue from authenticated E03 `t = 0.25` restart siblings through exact `t = 1.0`; require formal inspection, complete ranked products, zero strict LF failure counters, finite histories, forcing-work closure, accounting, and reconciliation before extension. | Job `4746182` reaches exact `t = 1.0` in `1.465556` node-hours. Strict LF failure counters remain zero; terminal `lf_hwproj = 65252911379`; complete two-group eight-rank snapshots and terminal eight-rank restart siblings are retained; sampled-history forcing-work relative residual is `8.616579960442532e-12`; reconciliation closes with three accepted ledger rows and no active reservation. The first retained `t = 1` recost artifact is preserved as superseded arithmetic history. Corrected conservative projection is `704.604815` matrix node-hours with `195.395185` margin. | Historical authorization was only `R02/s03_rankio_t1_t1p5`; that segment is now accepted under F-086. Keep R17 last. |
| F-086 | 2026-05-30 | Fresh corrected E03 R02 continuation through exact `t = 1.5` and corrected recost | F-085 authorized only one half-unit continuation so developed-state rate drift and corrected matrix arithmetic could be reviewed before another extension. | accepted `R02/s03_rankio_t1_t1p5`; job `4746356`; record-time inspection SHA-256 `5a66f4e037fe08f176b51c7f949db0377bc9270e666ae1e1c69dda85d7f96633`; recorded manifest SHA-256 `2172897379f87f7b2538ea27b6c12d0052a7bc37b1066068f9f70842d507e93e`; retained corrected R02 recost evidence JSON SHA-256 `b6fd6e8fcd939ca1baa06411a3cd6f04359b3bc2ff5b2a0097e2cc87b1763fc7` | High | Continue from authenticated E03 `t = 1.0` restart siblings through exact `t = 1.5`; require formal and independent inspection, complete ranked products, zero strict LF failure counters, finite histories, forcing-work closure, accounting, corrected recosting, and reconciliation before extension. | Job `4746356` reaches exact `t = 1.5` in `1.048611` node-hours. Strict LF failure counters remain zero; terminal `lf_hwproj = 106136472818`; complete two-group eight-rank snapshots and terminal eight-rank restart siblings are retained; sampled-history forcing-work relative residual is `4.978727845741857e-13`; reconciliation closes with four accepted ledger rows and no active reservation. The half-unit rate is `2.097222` node-hours per simulated unit. Corrected authorization-facing projection is `724.645556` matrix node-hours with `175.354444` margin. | Historical authorization was only `R02/s04_rankio_t1p5_t2`; that segment is now accepted under F-087. Keep R17 last. |
| F-087 | 2026-05-30 | Fresh corrected E03 R02 continuation through exact `t = 2.0` and late-time recost | F-086 authorized only one half-unit continuation so the corrected developed-state rate could be reviewed before another extension. | accepted `R02/s04_rankio_t1p5_t2`; job `4746435`; record-time inspection SHA-256 `fee92e6d8b484d1050266f6d5600577d19099df02335299114da59dd42050606`; recorded manifest SHA-256 `ff5ee6b8d408c0768ebf02bf4e5b0951609bd5cdd7f6591ebd2c924df942f9e6`; retained corrected R02 recost evidence JSON SHA-256 `c382b78ab9e21466648adf9d7ea38d2b407f0e986579f128bdfe6bc44bef79dd` | High | Continue from authenticated E03 `t = 1.5` restart siblings through exact `t = 2.0`; require formal and independent inspection, complete ranked products, zero strict LF failure counters, finite histories, forcing-work closure, accounting, recosting, and reconciliation before extension. | Job `4746435` reaches exact `t = 2.0` in `1.063611` node-hours. Strict LF failure counters remain zero; terminal `lf_hwproj = 123445535090`; complete two-group eight-rank snapshots and terminal eight-rank restart siblings are retained; sampled-history forcing-work relative residual is `3.034159691118515e-12`; reconciliation closes with five accepted ledger rows and no active reservation. The half-unit rate is `2.127222` node-hours per simulated unit. Authorization-facing projection is `728.845556` matrix node-hours with `171.154444` margin. | Historical authorization was only `R02/s05_rankio_t2_t2p5`; that segment is now accepted under F-088. Keep R17 last. |
| F-088 | 2026-05-30 | Fresh corrected E03 R02 continuation through exact `t = 2.5` and late-time recost | F-087 authorized only one half-unit continuation so the corrected developed-state rate could be reviewed before another extension. | accepted `R02/s05_rankio_t2_t2p5`; job `4746663`; record-time inspection SHA-256 `57fb0b9f1a19c0125aec0619802f2c41bd91ac97c51f32e70ceae9a434355da4`; recorded manifest SHA-256 `8051e68e17775aa4e7aac27cac3e65967bec668186e8a609b41dce7110899482`; retained corrected R02 recost evidence JSON SHA-256 `39115585cf7ca866e3aade1e53c69aab6ee9217361d75db7946748905e1c0e5c` | High | Continue from authenticated E03 `t = 2.0` restart siblings through exact `t = 2.5`; require formal and independent inspection, complete ranked products, zero strict LF failure counters, finite histories, forcing-work closure, accounting, recosting, and reconciliation before extension. | Job `4746663` reaches exact `t = 2.5` in `3825` seconds (`1.062500` node-hours). Strict LF failure counters remain zero; terminal `lf_hwproj = 171943591926`; complete two-group eight-rank snapshots and terminal eight-rank restart siblings are retained; sampled-history forcing-work relative residual is `3.0678457367645077e-12`; reconciliation closes with six accepted ledger rows and no active reservation. The half-unit rate is `2.125000` node-hours per simulated unit. Authorization-facing projection is `728.534445` matrix node-hours with `171.465555` margin. | Historical authorization was only `R02/s06_rankio_t2p5_t3`; that segment is now accepted under F-089. Keep R17 last. |
| F-089 | 2026-05-30 | Fresh corrected E03 R02 continuation through exact `t = 3.0` and hardened staged recost | F-088 authorized only one half-unit continuation so late-time rate drift could be reviewed before another extension. | accepted `R02/s06_rankio_t2p5_t3`; job `4746773`; record-time inspection SHA-256 `4b32798174222c9b9ca165285397ed9a02a411bb3e3dbf09e1f10611eebc70db`; recorded manifest SHA-256 `06b582e73d7b513e6380e5728bc10232949ab73b32bf57407791bdcc28dbfdae`; retained independent validation SHA-256 `200f60f69d6e42440211c557fb0e580ba9e3706761c18596568c511d4997c3ef`; retained corrected R02 recost evidence JSON SHA-256 `d2c5e27553426beb03c8c1882bd6473b682029c6ba67acd32bb06c5fb9a3ff01` | High | Continue from authenticated E03 `t = 2.5` restart siblings through exact `t = 3.0`; require formal and independent inspection, complete ranked products, zero strict LF failure counters, finite histories, forcing-work closure, accounting, reconciliation, staged recost review, and byte-preserving promotion before extension. | Job `4746773` reaches exact `t = 3.0` in `3879` seconds (`1.077500` node-hours). Strict LF failure counters remain zero; terminal `lf_hwproj = 184131904300`; complete two-group eight-rank snapshots and terminal eight-rank restart siblings are retained; sampled-history forcing-work relative residual is `2.8133353495278692e-12`; reconciliation closes with seven accepted ledger rows and no active reservation. The half-unit rate is `2.155000` node-hours per simulated unit. Authorization-facing projection is `732.734445` matrix node-hours with `167.265555` margin. Linnaeus independently approved byte-preserving promotion of the staged recost artifact. | Historical authorization was only `R02/s07_rankio_t3_t3p5`; that segment is now accepted under F-090. Keep R17 last. |
| F-090 | 2026-05-31 | Fresh corrected E03 R02 continuation through exact `t = 3.5`, recost stop-line review, and bounded threshold override | F-089 authorized only one half-unit continuation. The accepted clean interval exceeded the prior `4000`-second operational guard, so another segment required an explicit reviewed threshold decision rather than an implicit extension. | accepted `R02/s07_rankio_t3_t3p5`; job `4746953`; record-time inspection SHA-256 `4714af847c1ce10517258ea34f6151eb2a60d7fbd1bf1717d2de38d2e6822933`; recorded manifest SHA-256 `33ba23f4986b62461142b046c7d8e9568a8ae03bd50e67e8f4e8d1f55be1b642`; retained independent validation SHA-256 `aa37d73e4efcd2ea67cb8538f3e52fc49437207c558291c7e9c83ae20e4e4486`; retained corrected R02 recost evidence JSON SHA-256 `468e787db48d98661cd3ffe125829a78f135bb7ca4c9104c766a00b693f586a3` | High | Continue from authenticated E03 `t = 3.0` restart siblings through exact `t = 3.5`; require formal and independent inspection, complete ranked products, zero strict LF failure counters, finite histories, forcing-work closure, accounting, reconciliation, staged recost review, and byte-preserving promotion. If elapsed time exceeds the prior `4000`-second guard, stop for explicit review. | Job `4746953` reaches exact `t = 3.5` in `4144` seconds (`1.151111` node-hours). Strict LF failure counters remain zero; terminal `lf_hwproj = 192268745176`; complete two-group eight-rank snapshots and terminal eight-rank restart siblings are retained; sampled-history forcing-work relative residual is `2.6602098521945525e-13`; reconciliation closes with eight accepted ledger rows and no active reservation. The half-unit rate is `2.302222` node-hours per simulated unit. Authorization-facing projection is `753.345556` matrix node-hours with `146.654444` margin. The fixed `4000`-second guard rejected automatic extension. Linnaeus independently approved a staged-only recost artifact with an explicit `4500`-second cap for `s08` only, retaining `300` seconds before the Athena timeout and prohibiting automatic ratcheting. | Historical authorization was only `R02/s08_rankio_t3p5_t4`; that segment is now accepted under F-091. Keep R17 last. |
| F-091 | 2026-05-31 | Fresh corrected E03 R02 continuation through exact `t = 4.0`, independent recost provenance correction, and bounded timeout-profile revision | F-090 authorized only one half-unit continuation under a `4500`-second cap scoped to `s08`. The accepted clean interval remained above the baseline `4000`-second guard, so no next segment was authorized until a deliberate timeout-profile review closed. | accepted `R02/s08_rankio_t3p5_t4`; job `4747015`; record-time inspection SHA-256 `be705544f6c8f4a71f5501130b1e317555cdf1cac0663fba267acec82d967649`; recorded manifest SHA-256 `475ded71708ebe51bf1d696b9227ece14026e1b3884d1a024953928cf6ba7043`; retained independent validation SHA-256 `3f4d392437b3fe4a346656d3a4536575b7621c7d6049d58d524c196a84b3426c`; retained corrected R02 recost evidence JSON SHA-256 `dd5b0f5b82cf81005c1a481ee83d1127aa43861dcdd8c22765170e28b0ad6c99` | High | Continue from authenticated E03 `t = 3.5` restart siblings through exact `t = 4.0`; require formal and independent inspection, complete ranked products, zero strict LF failure counters, finite histories, forcing-work closure, accounting, reconciliation, staged recost review, and byte-preserving promotion. Review whether the next half-unit should receive a deliberate timeout revision or switch to quarter-unit increments. | Job `4747015` reaches exact `t = 4.0` in `4320` seconds (`1.200000` node-hours). Strict LF failure counters remain zero; terminal `lf_hwproj = 203117871743`; complete two-group eight-rank snapshots and terminal eight-rank restart siblings are retained; sampled-history forcing-work relative residual is `1.3875872397103518e-13`; reconciliation closes with `9/9/9` ledger rows/manifests/reservations, no active reservation, and no transaction. The half-unit rate is `2.400000` node-hours per simulated unit. Authorization-facing projection is `767.034445` matrix node-hours with `132.965555` margin. Linnaeus blocked the first staged artifact because it mislabeled the immediate predecessor threshold; the regenerated staged artifact distinguishes the `4000`-second baseline, `4500`-second predecessor cap, and `4800`-second `s09` cap. Linnaeus then independently approved byte-preserving promotion. | Archive and catalog the current committed controller state, reconcile, then prepare only `R02/s09_rankio_t4_t4p5` with Slurm walltime `01:40:00`, Athena timeout `01:30:00`, and a reviewed one-segment `4800`-second threshold retaining `600` seconds on both timeout margins. Do not ratchet automatically. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation. Inspect, account, reconcile, and recost before any further extension. Keep R17 last. |
| F-092 | 2026-05-31 | Fresh corrected E03 R02 continuation through exact `t = 4.5` and explicit one-segment timeout-profile renewal | F-091 authorized only one half-unit continuation under the reviewed `4800`-second cap scoped to `s09`. The accepted clean interval satisfies that cap, but no later segment inherits it automatically. | accepted `R02/s09_rankio_t4_t4p5`; job `4747087`; record-time inspection SHA-256 `0373fa3505beaadbd818db5271dae13bf69a41b1683bb561d558d8d71eff7581`; recorded manifest SHA-256 `c240bbec83817f50bbb3463f00fa8ce150c033af857e0859cc4d2313c4cda439`; retained independent validation SHA-256 `d4af66b36e115b6dde40fe41b157089b0b178ee87efb24ffa1a4e97e19fcc811`; retained corrected R02 recost evidence JSON SHA-256 `8b10b1786db520740b687d989207f3cda6e0123f9bf2b62e75a5df3849511561` | High | Continue from authenticated E03 `t = 4.0` restart siblings through exact `t = 4.5`; require formal and independent inspection, complete ranked products, zero strict LF failure counters, finite histories, forcing-work closure, accounting, reconciliation, staged recost review, and byte-preserving promotion. Review the next interval explicitly rather than inheriting the `s09` cap. | Job `4747087` reaches exact `t = 4.5` in `4349` seconds (`1.208056` node-hours). Strict LF failure counters remain zero; terminal `lf_hwproj = 207989256307`; complete two-group eight-rank snapshots and terminal eight-rank restart siblings are retained; sampled-history forcing-work relative residual is `1.9837412357887464e-12`; reconciliation closes with `10/10/10` ledger rows/manifests/reservations, no active reservation, and no transaction. The half-unit rate is `2.416111` node-hours per simulated unit. Authorization-facing projection is `769.290001` matrix node-hours with `130.710000` margin. Linnaeus independently approved byte-preserving promotion after replaying the full ten-segment evidence chain. Independent policy reviews renew the same timeout profile for `s10` only. | Archive and catalog the current committed controller state, reconcile, then prepare only `R02/s10_rankio_t4p5_t5` with Slurm walltime `01:40:00`, Athena timeout `01:30:00`, and a reviewed one-segment `4800`-second threshold retaining `600` seconds on both timeout margins. Do not ratchet automatically. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation. Inspect, account, reconcile, and recost before any further extension. Keep R17 last. |
| F-093 | 2026-05-31 | Fresh corrected E03 R02 continuation through exact `t = 5.0` and final explicitly reviewed half-unit probe | F-092 authorized only one half-unit continuation under the reviewed `4800`-second cap scoped to `s10`. The accepted clean interval satisfies that cap with narrowing headroom, so independent policy reviews authorize only one final half-unit probe before a mandatory fresh recost. | accepted `R02/s10_rankio_t4p5_t5`; job `4747146`; record-time inspection SHA-256 `8fb8069dc1f75448a8f79a38de98141c5d6dfed83bc3530eb633f8c009655014`; recorded manifest SHA-256 `8ce396a5f9e061b87d713636d257af0b8e57109eaffe48bac3945f028ae05744`; retained independent validation SHA-256 `082fa6f3fb04b1e46f50799c128bc7d4119595ad7c080a3ab0c30b8617c13daf`; retained corrected R02 recost evidence JSON SHA-256 `5256fb7fd9b75f175a16c2c16ee8a995b8d6d9281a92b70a979ea62171835bf2` | High | Continue from authenticated E03 `t = 4.5` restart siblings through exact `t = 5.0`; require formal and independent inspection, complete ranked products, zero strict LF failure counters, finite histories, forcing-work closure, accounting, reconciliation, staged recost review, and byte-preserving promotion. Review whether the next interval remains a defensible half unit or must split to quarters. | Job `4747146` reaches exact `t = 5.0` in `4415` seconds (`1.226389` node-hours). Strict LF failure counters remain zero; terminal `lf_hwproj = 220146873111`; complete two-group eight-rank snapshots and terminal eight-rank restart siblings are retained; sampled-history forcing-work relative residual is `5.300795515586916e-12`; reconciliation closes with `11/11/11` ledger rows/manifests/reservations, no active reservation, and no transaction. The half-unit rate is `2.452778` node-hours per simulated unit. Authorization-facing projection is `774.423334` matrix node-hours with `125.576666` margin. Linnaeus independently approved byte-preserving promotion after replaying the evidence chain. Linnaeus and Godel independently authorize one final half-unit probe under the same timeout profile for `s11` only. | Archive and catalog the current committed controller state, reconcile, then prepare only `R02/s11_rankio_t5_t5p5` with Slurm walltime `01:40:00`, Athena timeout `01:30:00`, and a reviewed one-segment `4800`-second threshold retaining `600` seconds on both timeout margins. Do not ratchet automatically. After clean `s11` inspection, accounting, reconciliation, and recost, propose quarter units if elapsed time exceeds `4800` seconds, remaining cap headroom falls below `300` seconds, or the reviewed trend estimator projects the next half-unit at or above `4800` seconds. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation. Keep R17 last. |
| F-094 | 2026-05-31 | Fresh corrected E03 R02 continuation through exact `t = 5.5`, hard-stop wording correction, and reviewed quarter-unit fallback | F-093 authorized one final half-unit continuation under the reviewed `4800`-second cap scoped to `s11`. The accepted clean interval leaves only `232` seconds of cap headroom, below the reviewed `300`-second floor. Independent policy review therefore requires quarter-unit proposals and clarifies that scientific, provenance, scheduler, storage, budget, or reconciliation failures block preparation rather than merely shorten the interval. | accepted `R02/s11_rankio_t5_t5p5`; job `4747202`; record-time inspection SHA-256 `928c281788b8ba359eb499aa607529b6c8eba2f8236e374f808d57fe239d9c1e`; recorded manifest SHA-256 `77192bba5dcfc378627d32b1ef316b5c08f77777fb003a505f1afe6d31322ac5`; retained independent validation SHA-256 `cfa3b6372fbfcaa65db60c33e826f67ff3ba0a301192095f17520e35c94a52a8`; retained corrected R02 recost evidence JSON SHA-256 `323d7cb3cbe108eee944045c189cff75ffc07d06833c49bcfa8315f07fad8d51`; audited one-use recost generator SHA-256 `0cfb8aaef913c3e677b147637704dd23d3b49b1710d32b072e4aad12a85e467e` | High | Continue from authenticated E03 `t = 5.0` restart siblings through exact `t = 5.5`; require formal and independent inspection, complete ranked products, zero strict LF failure counters, finite histories, forcing-work closure, accounting, reconciliation, staged recost review, and byte-preserving promotion. Define an explicit reviewed timing estimator before any quarter-unit successor. | Job `4747202` reaches exact `t = 5.5` in `4568` seconds (`1.268889` node-hours). Strict LF failure counters remain zero; terminal `lf_hwproj = 224321212896`; complete two-group eight-rank snapshots and terminal eight-rank restart siblings are retained; sampled-history forcing-work relative residual is `5.093197346120908e-12`; reconciliation closes with `12/12/12` ledger rows/manifests/reservations, no active reservation, and no transaction. The half-unit rate is `2.537778` node-hours per simulated unit. Authorization-facing projection is `786.323334` matrix node-hours with `113.676666` margin. The reviewed estimator selects `4721` seconds for the next half-unit and `2361` seconds for the proposed quarter-unit segment. Linnaeus and Newton independently approved byte-preserving promotion after replaying the evidence chain. | Archive and catalog the current committed controller state, reconcile, then prepare only `R02/s12_rankio_t5p5_t5p75` on one node from authenticated `s11` terminal siblings with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and a reviewed one-segment `2700`-second threshold retaining `600` seconds on both timeout margins. Do not ratchet automatically. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation. Inspect, account, reconcile, and recost before any further extension. Keep R17 last. |

### Implemented Core Decision Log

| Date | Decision | Rationale | Evidence still required |
| --- | --- | --- | --- |
| 2026-05-31 | Accept corrected E03 R02 continuation job `4747202` through exact `t = 5.5` and authorize only the reviewed quarter-unit fallback through `t = 5.75`. | Job `4747202` passes formal and independent inspection, records `1.268889` node-hours, and brings corrected E03 use to `12.250279`. Retained corrected F-094 evidence SHA-256 `323d7cb3cbe108eee944045c189cff75ffc07d06833c49bcfa8315f07fad8d51` projects at most `786.323334` matrix node-hours inside the `900.000000` envelope. The accepted `4568`-second half-unit leaves `232` seconds of cap headroom, so the reviewed s12 quarter profile uses Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and threshold `2700` seconds, retaining `600` seconds on both timeout margins. It is scoped only to `s12` and cannot ratchet automatically. | Archive and catalog the current committed controller state, reconcile, then prepare only `R02/s12_rankio_t5p5_t5p75` on one node from authenticated `s11` terminal siblings under the full F-094 profile. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation. Inspect, account, reconcile, and recost before any extension. Finish R02, execute R03--R16 sequentially, and execute R17 last. |
| 2026-05-31 | Accept corrected E03 R02 continuation job `4747146` through exact `t = 5.0` and authorize only one final reviewed half-unit probe through `t = 5.5`. | Job `4747146` passes formal and independent inspection, records `1.226389` node-hours, and brings corrected E03 use to `10.981390`. Retained corrected F-093 evidence SHA-256 `5256fb7fd9b75f175a16c2c16ee8a995b8d6d9281a92b70a979ea62171835bf2` projects at most `774.423334` matrix node-hours inside the `900.000000` envelope. The reviewed `s11` profile uses Slurm walltime `01:40:00`, Athena timeout `01:30:00`, and threshold `4800` seconds, retaining `600` seconds on both timeout margins. It is scoped only to `s11` and cannot ratchet automatically. | Archive and catalog the current committed controller state, reconcile, then prepare only `R02/s11_rankio_t5_t5p5`. Inspect, account, reconcile, and recost before any extension. After clean `s11` inspection, accounting, reconciliation, and recost, propose quarter units if elapsed time exceeds `4800` seconds, remaining cap headroom falls below `300` seconds, or the reviewed trend estimator projects the next half-unit at or above `4800` seconds. Any scientific, provenance, scheduler, storage, budget, or reconciliation failure blocks successor preparation. Finish R02, execute R03--R16 sequentially, and execute R17 last. |
| 2026-05-31 | Accept corrected E03 R02 continuation job `4747087` through exact `t = 4.5` and renew the reviewed timeout profile only for the `t = 4.5`--`5.0` continuation. | Job `4747087` passes formal and independent inspection, records `1.208056` node-hours, and brings corrected E03 use to `9.755001`. Retained corrected F-092 evidence SHA-256 `8b10b1786db520740b687d989207f3cda6e0123f9bf2b62e75a5df3849511561` projects at most `769.290001` matrix node-hours inside the `900.000000` envelope. The reviewed `s10` profile uses Slurm walltime `01:40:00`, Athena timeout `01:30:00`, and threshold `4800` seconds, retaining `600` seconds on both timeout margins. It is scoped only to `s10` and cannot ratchet automatically. | Archive and catalog the current committed controller state, reconcile, then prepare only `R02/s10_rankio_t4p5_t5`. Preserve formal inspection, accounting, projection review, shared-root queue audit, hardened reconciliation, and the one-segment threshold scope before any further extension. Finish R02, execute R03--R16 sequentially, and execute R17 last. |
| 2026-05-31 | Accept corrected E03 R02 continuation job `4747015` through exact `t = 4.0` and authorize only the reviewed `t = 4.0`--`4.5` continuation next under a non-ratcheting timeout-profile revision. | Job `4747015` passes formal and independent inspection, records `1.200000` node-hours, and brings corrected E03 use to `8.546945`. Retained corrected F-091 evidence SHA-256 `dd5b0f5b82cf81005c1a481ee83d1127aa43861dcdd8c22765170e28b0ad6c99` projects at most `767.034445` matrix node-hours inside the `900.000000` envelope. The reviewed `s09` profile uses Slurm walltime `01:40:00`, Athena timeout `01:30:00`, and threshold `4800` seconds, retaining `600` seconds on both timeout margins. It is scoped only to `s09` and cannot ratchet automatically. | Archive and catalog the current committed controller state, reconcile, then prepare only `R02/s09_rankio_t4_t4p5`. Preserve formal inspection, accounting, projection review, shared-root queue audit, hardened reconciliation, and the one-segment threshold scope before any further extension. Finish R02, execute R03--R16 sequentially, and execute R17 last. |
| 2026-05-31 | Accept corrected E03 R02 continuation job `4746953` through exact `t = 3.5` and authorize only the reviewed `t = 3.5`--`4.0` continuation next under a one-segment `4500`-second ceiling. | Job `4746953` passes formal and independent inspection, records `1.151111` node-hours, and brings corrected E03 use to `7.346945`. Retained corrected F-090 evidence SHA-256 `468e787db48d98661cd3ffe125829a78f135bb7ca4c9104c766a00b693f586a3` projects at most `753.345556` matrix node-hours inside the `900.000000` envelope. Its `4144` seconds exceed the prior `4000`-second guard; reviewed `4500` seconds retains `300` seconds before the Athena timeout and is scoped only to `s08`. | Historical authorization was only `R02/s08_rankio_t3p5_t4`; that segment is now accepted under F-091. |
| 2026-05-30 | Accept corrected E03 R02 continuation job `4746773` through exact `t = 3.0` and authorize only the conservative `t = 3.0`--`3.5` continuation next. | Job `4746773` passes formal and independent inspection, records `1.077500` node-hours, and brings corrected E03 use to `6.195834`. Retained corrected F-089 evidence SHA-256 `d2c5e27553426beb03c8c1882bd6473b682029c6ba67acd32bb06c5fb9a3ff01` projects at most `732.734445` matrix node-hours inside the `900.000000` envelope. | Historical authorization was only `R02/s07_rankio_t3_t3p5`; that segment is now accepted under F-090. |
| 2026-05-30 | Accept corrected E03 R02 continuation job `4746663` through exact `t = 2.5` and authorize only the conservative `t = 2.5`--`3.0` continuation next. | Job `4746663` passes formal and independent inspection, records `1.062500` node-hours, and brings corrected E03 use to `5.118334`. Retained corrected F-088 evidence SHA-256 `39115585cf7ca866e3aade1e53c69aab6ee9217361d75db7946748905e1c0e5c` projects at most `728.534445` matrix node-hours inside the `900.000000` envelope. | Historical authorization was only `R02/s06_rankio_t2p5_t3`; that segment is now accepted under F-089. |
| 2026-05-30 | Accept corrected E03 R02 continuation job `4746435` through exact `t = 2.0` and authorize only the conservative `t = 2.0`--`2.5` continuation next. | Job `4746435` passes formal and independent inspection, records `1.063611` node-hours, and brings corrected E03 use to `4.055834`. Retained corrected F-087 evidence SHA-256 `c382b78ab9e21466648adf9d7ea38d2b407f0e986579f128bdfe6bc44bef79dd` projects at most `728.845556` matrix node-hours inside the `900.000000` envelope. | Historical authorization was only `R02/s05_rankio_t2_t2p5`; that segment is now accepted under F-088. |
| 2026-05-30 | Accept corrected E03 R02 continuation job `4746356` through exact `t = 1.5`, promote F-083/F-084 controller hardening, and authorize only the conservative `t = 1.5`--`2.0` continuation next. | Job `4746356` passes formal and independent inspection, records `1.048611` node-hours, and brings corrected E03 use to `2.992223`. Retained corrected F-086 evidence SHA-256 `b6fd6e8fcd939ca1baa06411a3cd6f04359b3bc2ff5b2a0097e2cc87b1763fc7` projects at most `724.645556` matrix node-hours inside the `900.000000` envelope. Promoted canonical revision `9689c269bf329542815a1b2b137881126964b05c` passes focused tests, documentation build, syntax check, diff check, and hardened reconciliation. | Historical authorization was only `R02/s04_rankio_t1p5_t2`; that segment is now accepted under F-087. |
| 2026-05-30 | Accept the corrected E03 R02 continuation through exact `t = 1.0` and authorize only the conservative `t = 1.0`--`1.5` continuation next. | Job `4746182` passes formal and independent inspection, records `1.465556` node-hours, and brings corrected E03 use to `1.943612`. The retained F-085 artifact SHA-256 `7d11e8a0004e24417167aeb1f186f8b1f68c8eb9d6206467a708a7cad16cf252` is preserved as superseded arithmetic history; corrected conservative projection at that boundary is `704.604815` matrix node-hours inside the `900.000000` envelope. | Historical authorization was only `R02/s03_rankio_t1_t1p5`; that segment is now accepted. |
| 2026-05-30 | Preserve committed live-helper enforcement and reserve case-level `analysis/` for retained evidence. | The F-081 recost placement exposed an orphan-discovery ambiguity. Commits `9480e62764528a3f40066d22a192f0e99b369891` and `ef1e42fa088203ac9ef6ec8e47e668db4fb95a3c` close that ambiguity and prevent historical-bundle authentication from masking dirty live helper bytes. | Preserve the focused Stage I regression and verified source bundle at later lifecycle transitions. |
| 2026-05-30 | Accept the corrected E03 R02 continuation through exact `t = 0.25` and authorize only the bounded `t = 0.25`--`1.0` continuation next. | Job `4746154` passes formal and independent inspection, records `0.289167` node-hours, and brings corrected E03 use to `0.478056`. Retained F-081 evidence SHA-256 `eb6071cf53d453b2f75707003616ac9cbdfe025d62969d00e34e7948bee3310c` projects `700.868889` provisional matrix node-hours inside the `900.000000` envelope. | Historical authorization was only `R02/s02_rankio_t0p25_t1`; that segment is now accepted. |
| 2026-05-30 | Accept the fresh corrected E03 R02 pilot and retain the controller-provenance transition before continuation. | Job `4745922` reaches exact `t = 0.1`, passes formal inspection with zero strict LF safety counters and complete initial/final eight-sibling ranked snapshot and restart groups, and records `0.188889` node-hours. The record-time `segment_inspection.json` and embedded `scientific_inspection` payload SHA-256 are `f4b0a31304e31610ec9ec83133f79b747fb0d8cf0f9b3b013febc1cddc2ae784`; forcing-work relative residual is `1.920820308941694e-11`. F-080 is archived by transition `05cb4c324bfd8feec72ebdeb33b1961c9fde70bf` and retained evidence SHA-256 `46fc1c4054e75f4224302be1895c1b9eaae88097544512ac067c53376e542e34`. | Prepare an authenticated R02 continuation from the accepted checkpoint. |
| 2026-05-30 | Isolate corrected work beneath `E03-forcing-policy`, require a canonical immutable-build qualification token before preparation, and retain schema-v2 panel status as fail-closed until numeric criteria review. | The replacement forcing policy, restart timestamp, provenance controller, and panel-reference bindings are part of one auditable execution contract. Before F-078, the canonical E03 summary reported an absent token and therefore blocked production. | Closed by F-078 after immutable build archival, reviewed `g024`--`g031` qualification, token retention, and clean reconciliation. |
| 2026-05-30 | Stop every new E02 preparation or submission, preserve its products only as pipeline and cost evidence, correct the forcing contract, and start a fresh epoch from `t = 0`. | The independent F-076 source-to-code audit found that E02 random roles inherit generic projected forcing and paper decks rely on an implicit parabolic-spectrum default, while planar roles do not generally enforce `grad_perp dot u_perp = 0` for retained `k_z != 0` modes. Job `4745498` was cancelled and recorded `aborted`; the stop-line evidence JSON has SHA-256 `de85c443c3bc9bab795f7f1903e5de224d7d86b553ac40e9dd7249986c3df211`. | Explicit restart-retained MKS24 random and Alfvenic forcing policies, explicit paper-deck `spectrum = power_law`, fresh qualification, immutable executable archival, and a new execution epoch. |
| 2026-05-30 | Continue only the frozen matrix one inspected segment at a time, with R17 last. | R02 job `4745305` reaches exact `t = 2.5` cleanly in `1.061944` node-hours. Its latest late-time rate updates the mapped projection to `728.164877` node-hours, leaving `171.835123` node-hours of margin. | Preserve formal inspection, accounting, projection review, and shared-root queue audit for every segment. |
| 2026-05-30 | Continue R02 conservatively from `t = 2.0` to `2.5` under the existing F-071 envelope. | R02 job `4744913` reaches exact `t = 2.0` cleanly in `1.068889` node-hours. Its latest late-time rate updates the mapped projection to `730.081697` node-hours, leaving `169.918303` node-hours of margin. | Continue R02 from `t = 2.0` to `2.5`, inspect, and recost before lengthening later segments. |
| 2026-05-30 | Continue R02 conservatively from `t = 1.5` to `2.0` under the existing F-071 envelope. | R02 job `4744518` reaches exact `t = 1.5` cleanly in `1.052222` node-hours. Its late-time rate updates the mapped projection to `725.464938` node-hours, leaving `174.535062` node-hours of margin. | Continue R02 from `t = 1.5` to `2.0`, inspect, and recost before lengthening later segments. |
| 2026-05-30 | Continue R02 conservatively and reject Stage I normal-QOS preparations above two hours. | R02 job `4744249` reaches exact `t = 1.0` cleanly in `1.452778` node-hours. The longer developed rate updates the mapped projection only modestly to `702.195370` node-hours. Slurm preflight rejects the released unsubmitted `02:15:00` attempt because normal QOS is capped at 120 minutes. | Continue R02 from `t = 1.0` to `1.5`, inspect, and recost before lengthening later segments. |
| 2026-05-30 | Authorize only the frozen `R02`--`R17` mapped matrix under a `900.000000` node-hour E02 envelope and sequential inspection. | R17 jobs `4744210` and `4744230` preserve 27 meshblocks per rank, pass strict inspection, and close the developed high-resolution rate. The continuation-aware mapped projection is `697.019444` node-hours with `202.980556` node-hours of envelope margin. | Resume R02 first, inspect each segment, retain accepted snapshots plus final two restart groups, and complete R17 last to control sequential peak storage. |
| 2026-05-30 | Preserve completed `R02` timing evidence and authorize only an `8.0` node-hour eight-node `R17` high-resolution timing probe next. | R02 jobs `4744198` and `4744205` reach native cadence `t = 0.25` cleanly. Their authenticated continuation projects `19.000000` node-hours per standard `t = 10` case, narrowing the matrix bracket to `424.145556`--`576.145556` node-hours. The residual spread is isolated to R17. | Run and inspect one bounded R17 timing probe with 64 GPU ranks so each rank retains the standard 27-meshblock load, then revise mapped-matrix runtime/storage authorization; no broader matrix execution is authorized yet. |
| 2026-05-30 | Preserve completed `R16` and authorize only a `2.0` node-hour `R02` standard-layout timing pilot next. | R16 completes cleanly in `6.145556` node-hours, but cell-count-only and CFL-aware mapped-matrix projections span `1087.763412` to `2956.012436` node-hours. `R02` holds active-Alfvenic beta-10 physics fixed while changing from `96 x 96 x 192` to `192 x 192 x 384`, making it the smallest discriminating production-path measurement. | Run and inspect one bounded R02 timing segment, then revise mapped-matrix runtime/storage authorization; no broader matrix execution is authorized yet. |
| 2026-05-30 | Expand only the fresh `E02` `R16` completion envelope from the initial `4.0` node-hour pilot cap to `7.0` node-hours. | Accepted jobs `4743735` and `4743933` establish clean rank-local execution and modal-restart continuation through `t = 3.5`. Their developed-state rate projects `5.639721` node-hours through `t = 10`, leaving `1.360279` node-hours of margin under the revised cap. | Complete and analyze `R16` through `t = 10`, then review the measured matrix runtime and storage reservation; no other mapped case is authorized yet. |
| 2026-05-30 | Start MKS24 production over as isolated epoch `E02-modal-driver`, with a fresh incremental `4000` node-hour planning ceiling and immutable `E01-pre-modal-driver` history. | F-056 invalidates old cap-fraction products for comparison, and F-057 makes old forcing restarts incompatible with the modal-driver executable. A clean namespace prevents archival segments from contaminating future discovery, bundles, and cost decisions. | Commit the locally implemented epoch-aware Stage I helper and separate ledger paths; complete replacement-driver qualification, fresh `R16` pilot from `t = 0`, and recosted Stage I reservation. |
| 2026-05-30 | Treat `/lustre/orion/ast207/proj-shared/dfielding/CGL` as a shared campaign root, not as an MKS24-exclusive directory. | Exploratory beta-25 comparison jobs coexist with MKS24 archives, and the prior Stage I helper filtered queue conflicts by `cgl_mks24_` job names only. MKS24 work must not modify, prune, or unintentionally overlap another live campaign. | Commit the locally implemented all-user-job and active-campaign-record preflight; verify the root is quiescent or explicitly reviewed before the next Stage I submission. |
| 2026-05-29 | Stop admitting the pre-replacement `R16` lineage after clean inspection of job `4686032`; account its cost as `rejected` for current reproduction evidence. | The merged turbulent-driver implementation changes the forcing/restart-state contract, so an old-executable terminal restart cannot qualify continuation under the new executable. | Replacement-driver forcing/restart/work-accounting qualification, immutable production build provenance, and a fresh Stage I lineage from `t = 0`. |
| 2026-05-24 | Keep the existing LF source placement and apply collisions after each LF half-sweep with explicit duration `dt/2`. | This repairs physical elapsed time without changing the protected LF representation/task ownership; both collision calls together span one `dt`. | Finite-collision restart and `compare` workflow passed on CPU; broader MPI/GPU gates remain |
| 2026-05-24 | Preserve `cgl_firehose_threshold = oblique` as the public default and add explicit `parallel` selection for MKS24. | Existing branch inputs retain their historical convention; paper inputs cannot accidentally depend on that default. | Paper pgen/decks and paper manifests |
| 2026-05-24 | Define firehose emergency overshoot at `p_perp - p_parallel <= -1.5 B^2` (`beta Delta <= -3`) and count it independently of `backup_limiters`. | The emergency condition must remain distinct from normal finite-rate activity when the physical MKS24 activation threshold is `beta Delta <= -2`. | Reduced nonlinear validation of overshoot interpretation |
| 2026-05-24 | Reuse and repair the built-in turbulence driver for active paper smoke rather than copy donor forcing logic into the pgen. | Shared forcing now records a deterministic seed, has explicit `B0 || z` Alfvenic orientation, preserves restart state, advances OU state once per cycle, and supplies RK-consistent energy work. | Calibrate published injection/spectral normalization and qualify larger runs |
| 2026-05-24 | Qualify reduced `passive_delta` only as isothermal-MHD flow with diagnostic CGL/LF thermodynamics, requiring matching `mhd/passive = true`. | Passive fluxes and CFL speeds are independent of diagnostic anisotropy; a paired-anisotropy regression shows identical density, velocity, and magnetic evolution under identical forcing. | Long-time active/passive statistics, execution of defined standard decks, and spatial paper observables |
| 2026-05-24 | Stage official arXiv source version `2405.02418v2` into ignored result storage with SHA-256 provenance rather than vendoring paper contents. | A pinned source archive supports exact reference interpretation while leaving third-party redistribution outside the code branch. | Analysis bundles retain pinned checksum provenance; F-034/F-036/F-037 populate Figures 2(b), 12, and 13(a)-(b), while remaining panel curves/comparisons remain open |
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
| 2026-05-25 | Treat G-008/G010b/G011 debug measurements as storage/I/O reconnaissance, not a paper-scale runtime cost model. | Startup-length runs expose elapsed time, memory, output volume, and checkpoint size; G011 adds exact standard-layout shared-file timing and size, but no retained run measures late-time timestep evolution or statistical-production duration. | G010b retains reduced timing and its excessive cadence warning; G011 file sizes plus F-041/F-049 case additions feed the estimated `880.368` GB sequential allocation for seventeen unique runs; only a G013-based scaling estimate, not measured production runtime, is available |
| 2026-05-25 | Stop runtime costing after G012 and treat limiter-active RKL2 STS as unqualified for the screened finite-rate nonlinear state. | The finite-rate hard-wall reduced nonlinear run and local reduced branches with unlimited STS or caps `20`, `10`, `5`, `2`, and `1` all cross the strict firehose emergency bound in an LF split stage; reducing the cap through the documented conservative setting delays but does not resolve the observed failure. | F-031 is addressed by F-033; G012 consumes `0.188889` node-hours and remains invalid for extrapolation, while corrected G013 supplies new reduced GPU evidence |
| 2026-05-25 | Archive the STS timestep/integrator contract in paper workflow manifests. | F-031 makes the split-step cap a qualification-critical model choice; a result cannot be interpreted or repeated if the manifest omits it. | F-032 metadata regression records guarded-deck `sts_max_dt_ratio = -1.0`; G013 archives and passes that setting with the F-033 hard-wall model |
| 2026-05-25 | Represent the MKS24 hard-wall limit explicitly with `limiter_hardwall = true` in hard-wall paper decks rather than finite-rate pressure relaxation at `limiter_nu_coll = 1.0e10`. | The pinned paper source describes anisotropy effectively pinned at thresholds; F-031 locates a same-state internal LF-stage overshoot under finite-rate scattering, while the energy-preserving constrained state passes locally through `t = 2.0`. Finite-rate `20`/`200` decks remain unchanged. | F-033 CPU/exact-state validation and corrected G013 reduced Frontier GPU execution pass; standard costing, proposal approval, and production/statistical qualification remain required |
| 2026-05-25 | Implement CGL pressure mechanical work first as a snapshot-derived local decomposition with a transfer cross-check and sparse-time quadrature estimate. | Retained fields determine `p_perp div(u) - Delta p (b b : grad(u))` without changing the integrator; trapezoidal snapshot-time integration aids interval interpretation, but neither result is an applied budget identity. | F-023 synthetic sign/transfer/quadrature oracles and corrected `20260525-paper-convergence-rk-work-v1` products; exact production interval closure remains open |
| 2026-05-25 | Admit paper reference curves only through an external checksum- and uncertainty-qualified manifest. | The staged arXiv source supplies panel PDFs but not numeric curve tables; separating source data from code and rejecting absent uncertainties permits author-data or documented-digitization comparisons without overstating simulation evidence. | F-024 exact-curve and source-checksum-rejection regressions plus F-034/F-036/F-037/F-039/F-040/F-041/F-043/F-044/F-046 checksum-bound extraction for Figures 2(b), 4(a), 5(b), 7, 8, 9, 11, 12, and 13; remaining panels and production comparisons remain open |
| 2026-05-25 | Expose threshold-volume histories as reference-comparable curves while collapsing duplicate terminal timestamps. | Figure 2(b) is a time-series result rather than a snapshot statistic, and Athena repeats an unchanged terminal history row; retaining the final equal-time value gives ordered coordinates without changing the measured endpoint. | F-030 synthetic exact-history comparison passes; F-034 now supplies eight digitized Figure 2(b) histories, while standard-run comparison remains open |
| 2026-05-25 | Complete Figure 2(b) coverage with four additional guarded standard definitions and vector-derived external reference histories. | The official panel shows active and passive Alfvenic/random cases at beta 10 and 100; omitting any series prevents the implemented comparison interface from evaluating the displayed result. | F-034 adds nine standard definitions in total and external `digitized_fig2b_v1/curves.json`; no production execution is authorized by this definition/reference work |
| 2026-05-25 | Define the nonnominal Figure 12 heat-flux sensitivity cases in a separate guarded workflow. | The pinned paper source varies the active beta-10 Alfvenic heat-flux parameter by factors of 100 around the nominal case; two additional definitions avoid duplicating nominal active/passive output while making the intended campaign auditable. | F-035 adds two `paper-heat-flux` decks and debug exclusion; F-036 supplies reference curves; authorized production runs and comparison remain open |
| 2026-05-25 | Expose peak-alignment versus physical `k_perp` as a reference product and extract both Figure 12 vector panels. | Figure 12 compares scale-dependent peak alignment, not a single shell's PDF, and its lower spectrum panel already maps to an existing analyzed product. | F-036 adds explicit analysis-shell forwarding and external `digitized_fig12_v1/curves.json`; authorized runs and scientific comparison remain open |
| 2026-05-25 | Extract only directly comparable Figure 13 panels and validate multiple source PDFs in one digitized manifest. | Panels (a)-(b) use existing velocity-spectrum and `beta Delta` PDF products, whereas panels (c)-(d) include plotted normalizations that cannot be represented honestly by the current raw strain/transfer products. | F-037 adds external `digitized_fig13_v1/curves.json` for seven curves and records panels (c)-(d) as excluded product gaps; authorized limiter execution and comparison remain open |
| 2026-05-25 | Supersede dimensional reference-curve admission until the paper-to-AthenaK code-unit transform is qualified. | MKS24's beta-100 footnote reports `p0 = 100`, while the normalized AthenaK decks use `p0 = 50`; this does not affect dimensionless `beta Delta` or alignment products, but it prevents honest raw spectral/strain/transfer comparisons without a declared conversion. | F-038 emits corrected `digitized_fig12_v2`/`digitized_fig13_v2` manifests containing only dimensionless admitted curves; earlier `v1` outputs remain archived as provisional evidence |
| 2026-05-25 | Implement MKS24's stated normalized transfer comparator and admit dimensionless transfer panels. | The paper defines `T_total ~= E_K (2 pi u_rms/L_perp)`, so Figure 7's lower panel and Figure 13(d) do not require the unresolved dimensional ordinate transform. | F-039 adds `pressure_transfer.transfer_normalized_by_total` and external `digitized_fig7_v1`/`digitized_fig13_v3` manifests; production comparisons remain open |
| 2026-05-25 | Admit Figure 9 peak-alignment references through existing standard-case mappings. | Its `|cos(theta)|` ordinate is dimensionless and all eight legend cases already exist in the guarded standard matrix, so it is not blocked by the dimensional normalization audit. | F-040 adds external `digitized_fig9_v1/curves.json`; production comparisons remain open |
| 2026-05-25 | Define Figure 11 scale-separation execution scope and admit its lower alignment panel. | The lower panel is dimensionless and directly uses `alignment_peak.cos_theta`, but comparing its resolution dependence requires explicit `n_perp = 96` and `384` inputs in addition to the standard `192` case. | F-041 adds guarded `paper-scale-separation`, external `digitized_fig11_v1/curves.json`, and revised estimated resource scope; execution/comparison remain open |
| 2026-05-25 | Stop admitting additional PDF panels at the reviewed curve/reference boundary. | Figure 8 has physically comparable axes but only raster heatmap/colorbar content, so deriving normalized probability slices requires a reviewed color mapping and uncertainty contract; raw-density and dimensional vector panels still lack qualified paper-to-AthenaK transforms or matching products. | F-042 retains the source/object/checksum audit without generating a speculative manifest; machine-readable author data, reviewed raster calibration, or explicit normalization work is required before expanding reference coverage |
| 2026-05-25 | Supersede the Figure 8 portion of F-042 with a calibrated selected-shell raster extraction contract. | The pinned PDF's labeled linear colorbar decodes its embedded RGB surfaces to absolute probability density, and sampled slices satisfy the published per-`k_perp` unit-normalization constraint within a conservative tolerance. | F-043 adds external `digitized_fig8_v1/curves.json` for ten `alignment.<shell>` curves at shells `8,16,32,64,128`; raw-density and dimensional-panel exclusions and production comparison gates remain open |
| 2026-05-25 | Supersede the Figure 4(a) portion of F-042 through its explicitly normalized density coordinate. | The panel labels its abscissa `rho/<rho>`, which differs from the implemented `rho/<rho>-1` PDF only by an exact unit translation and needs no donor code-unit scale. | F-044 adds external `digitized_fig4a_v1/curves.json` for eight `pdf.density_fluctuation` curves; Figure 4(b), dimensional panels, and production comparison gates remain open |
| 2026-05-25 | Retain Figure 5(b) outside admitted references until its conditioned-structure-function observable is implemented. | Unlike the corrected Figure 4(a) PDF and Figure 8 alignment slices, Figure 5(b)'s dimensionless axes did not map to an existing analyzer curve: its characteristic eddy lengths are derived from local-field-conditioned structure functions. | F-045 documented the product boundary; F-046 supersedes that defer-only state with an opt-in product and checksum-bound extraction, without authorizing production execution |
| 2026-05-25 | Admit Figure 5(b) only through an explicitly sampled local-field structure-function contract. | The paper and cited method define three-point second-order structure functions, local three-point magnetic-field conditioning, equal-power inversion, and `L_perp` length normalization; the computationally heavier estimate should be opt-in and reproducible rather than silently imposed on every snapshot analysis. | F-046 adds archived `--eddy-samples`, `--eddy-bins`, and `--eddy-seed`, two `eddy_anisotropy.*` products, and external `digitized_fig5b_v1/curves.json`; production comparison remains open |
| 2026-05-25 | Implement Figure 2(a)'s joint product without treating its raster as admitted reference data. | The panel's plotted coordinates define a calculable joint pressure-density PDF in each AthenaK snapshot, but its absolute pressure-scale/color-density mapping is not established by the pinned source and cannot be represented by the current curve-only manifest. | F-047 adds `pressure_density_joint` products and rendering; raster extraction and quantitative comparison remain open until a reviewed transform or author numeric data exists |
| 2026-05-25 | Complete missing spectral analysis fields without admitting unqualified spectral references. | The plotted field definitions for Figures 3, 4(b), and 6(a) can be emitted from retained snapshots independently of the unresolved donor ordinate normalization; a Fourier-projected compressive flow product is necessary before Figure 3 can be compared at all. | F-048 adds compressive, normalized-density, and separate thermal/magnetic-pressure spectra plus rendering and projection oracles; spectral reference transforms remain open and F-049 supplies the missing guarded Figure 3 definitions |
| 2026-05-25 | Define Figure 3's missing compressive-case execution scope without launching it. | The panel legends add an active random beta-1 run and a beta-100 random run with `t_corr = L_parallel/v_th`; the latter is `0.2` in the normalized AthenaK beta-100 deck, while four other panel cases already exist in `paper-standard`. | F-049 adds guarded `paper-compressive`, two definition-only decks, debug exclusion, and a revised seventeen-unique-run proposal; execution and comparison remain open |
| 2026-05-25 | Extend admitted reference infrastructure to sampled two-dimensional surfaces without admitting Figure 2(a)'s unresolved raster. | F-047 provides joint-PDF analysis grids, but a curve-only manifest cannot compare them once a valid transform or author dataset becomes available. | F-050 adds `surfaces` entries with positive `z_uncertainty`, bilinear `pressure_density_joint.*` comparison, and residual plotting; F-051 decodes excluded raw samples while donor-pressure conversion remains open |
| 2026-05-25 | Decode Figure 2(a)'s calibrated raster only into excluded donor-coordinate samples. | The pinned PDF provides a numeric linear colorbar and numeric pressure axes, resolving the RGB/raster mapping, but both axes are dimensionful and any later pressure scale conversion must also transform the joint-PDF density and its uncertainty by the two-dimensional Jacobian. | F-051 adds `scripts/digitize_cgl_lf_mks24_fig2a.py` and external `digitized_fig2a_v1/excluded_surfaces.json`; no analyzer-compatible manifest is emitted until `x_A = s x`, `y_A = s y`, `z_A = z/s^2`, and `sigma_z_A = sigma_z/s^2` are qualified |
| 2026-05-25 | Admit Figure 2(a) surfaces through a pressure-only beta-10 transform while retaining broader dimensional exclusions. | The pinned MKS24 source fixes the donor beta-10 initial pressure as `p0 = 10` by combining fixed `B0`, linear `p0(beta0)`, and the beta-100 `p0 = 100` code-unit statement; the matched AthenaK deck explicitly uses `p0 = 5`. This establishes `s = 0.5` for both Figure 2(a) pressure axes and the joint-PDF Jacobian, but not the ordinate conventions of unrelated spectra or strain panels. | F-052 requires pinned `MKS24.tex` and `fig2a.pdf`, emits external `digitized_fig2a_v2/surfaces.json`, and retains `digitized_fig2a_v2/donor_surfaces.json` for audit |
| 2026-05-25 | Keep remaining absolute spectral and strain ordinates excluded at an explicit external-evidence boundary. | MKS24 and its cited donor manuscript define shell spectra and the donor's plotted `phi^2 L_perp` units, but neither states the discrete Fourier convention needed to convert absolute ordinates; targeted public-source/data discovery found no linked diagnostic implementation or machine-readable deposit. | F-053 records the source and discovery audit; a later admission requires author/archive numeric data, donor diagnostic code, or an explicit discrete-normalization statement |
| 2026-05-26 | Retain the F-053 exclusion after a targeted external-data follow-up while production is queued. | The pinned source still supplies no discrete Fourier convention; retrieved official journal-page text contained no `Data availability` heading, and exact arXiv-ID Zenodo and exact-title GitHub repository searches returned no matching numeric-data or donor-code record. | Leave the absolute spectral/strain panels blocked pending author/archive numeric data, donor diagnostic code, or an explicit discrete-normalization statement |
| 2026-05-25 | Replace shared retained output with rank-local products and bundle only authenticated restart ancestry. | The shared-I/O `R16` diagnostic pause made the original retained-output path impractical, while unrelated clean partial branches and sampled history boundaries must not contaminate or block an accepted lineage. | F-054 records the `R16` production jobs, lineage/cadence regressions, and diagnostic prefix bundles; full `t = 10` acceptance and comparisons remain open |
| 2026-05-26 | Align the Tier 3 minimum execution matrix with the frozen source-mapped Stage I manifest. | The published beta-1 compressive role is the active-random `R10` case; retaining the earlier active-Alfvenic beta-1 requirement would direct execution of an unmapped definition outside the recorded authorization. | F-055 replaces the stale generic case row with canonical mapped groups; any later active-Alfvenic beta-1 execution requires a demonstrated panel role and revised authorization |
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
| Figure 2(a) surface-reference extractor | `scripts/digitize_cgl_lf_mks24_fig2a.py` |
| Figure 2(b) reference extractor | `scripts/digitize_cgl_lf_mks24_fig2b.py` |
| Figure 4(a) normalized-density PDF extractor | `scripts/digitize_cgl_lf_mks24_fig4a.py` |
| Figure 5(b) normalized eddy-anisotropy extractor | `scripts/digitize_cgl_lf_mks24_fig5b.py` |
| Figure 7 lower-panel reference extractor | `scripts/digitize_cgl_lf_mks24_fig7.py` |
| Figure 8 selected-shell alignment-PDF extractor | `scripts/digitize_cgl_lf_mks24_fig8.py` |
| Figure 9 alignment reference extractor | `scripts/digitize_cgl_lf_mks24_fig9.py` |
| Figure 11 lower-panel alignment reference extractor | `scripts/digitize_cgl_lf_mks24_fig11.py` |
| Figure 12 reference extractor | `scripts/digitize_cgl_lf_mks24_fig12.py` |
| Figure 13(b),(d) dimensionless reference extractor | `scripts/digitize_cgl_lf_mks24_fig13.py` |
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
paper-standard/limiter-scan/heat-flux/compressive/scale-separation inputs
explicitly select `parallel`. Standard paper execution and analysis remain
open.

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
analyzer also accepts optional external reference-data manifests that require
checksums and pointwise uncertainty values before calculating curve or sampled
surface comparison residuals or producing overlay figures.
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
The working AthenaK decks choose `rho0 = B0 = v_A = 1`, so the beta-100
initial pressure is `p0 = 50`. MKS24 separately states in the Figure 13(c)
discussion that `p0 = 100` in its code units for the beta-100 limiter runs.
That statement does not by itself establish the donor density, magnetic
normalization, or a transform for dimensional plotted ordinates. Therefore
the normalized AthenaK decks are retained, but raw dimensional reference
curves are not admitted for comparison until the observable conversion is
explicitly qualified.

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
| Figure 2: pressure-density distributions and unstable-volume history | Active and passive, Alfvenic/random, beta 10 and 100 | PDFs of `delta p_parallel`, `delta p_perp`, density; volume fraction beyond thresholds versus time | Reduced threshold-volume, steady-window PDFs, panel-(a)-coordinate joint pressure-density products, source-qualified color-calibrated sampled surfaces, and all eight vector-extracted panel-(b) histories are implemented with checksummed provenance; quantitative comparison still requires authorized standard execution |
| Figures 3-4: compressive/density behavior | Active/passive; beta and forcing variants | Density PDFs and density/compressive spectra | Steady-window density PDFs plus `spectra.compressive_velocity` and normalized `spectra.density_fluctuation` products/rendering are implemented under F-048; F-049's guarded `paper-compressive` workflow supplies Figure 3's missing active random beta-1 and beta-100 sonic-correlation definitions while reusing four standard cases; Figure 4(a)'s explicit `rho/<rho>` coordinate is mapped to eight admitted curves under F-044; Figure 4(b)'s raw `E_rho` normalization and Figure 3 reference transform remain unqualified, and required case execution/comparison remains incomplete |
| Figures 5-6: turbulent spectra and rate-of-strain suppression | Active/passive, Alfvenic/random | Velocity/magnetic/pressure spectra; parallel/perpendicular velocity-gradient spectra; local-field eddy scales | Local-field velocity-gradient/strain products, separate `p_parallel`, `p_perp`, and AthenaK `B^2/2` spectra/rendering, plus opt-in conditioned three-point `eddy_anisotropy.velocity_perp`/`.magnetic_perp` products are implemented; Figure 5(b)'s four normalized eddy-scale curves are admitted under F-046, while Figure 5(a)/6 spectral panels retain reference-transform gaps; target runs/comparison missing |
| Figure 7: pressure-anisotropy organization and transfer | Active/passive, Alfvenic/random, beta 10 | `grad_parallel Delta p`, `grad_perp Delta p` spectra; `T_Delta p(k_perp)/T_total` | Steady-window gradient spectra, transfer partition/closure, snapshot anisotropic-work cross-check, MKS24-normalized transfer product, and dimensionless lower-panel reference curves are implemented; dimensional upper-panel reference transform, target runs, and comparison remain open |
| Figures 8-9: alignment diagnostic | Active/passive; beta/forcing variants | Scale-dependent PDFs and peak locations of local-field/rate-of-strain eigenvector alignment | Steady-window selected-shell alignment PDFs, reference-comparable peak curves, rendered peak summaries, and synthetic finite check implemented; Figure 8 active/passive beta-10 Alfvenic raster panels emit ten calibrated `alignment.<shell>` curves with explicit uncertainty under F-043; Figure 9 eight peak-alignment reference curves are vector-extracted for the standard case matrix |
| Figure 10: kinetic comparison context | Not reproduced directly by CGL-LF | Documentation of non-reproducible kinetic comparator and any qualitative CGL-only comparison | Out of direct scope |
| Figure 11: scale-separation behavior | Active Alfvenic beta 10 at `n_perp = 96, 192, 384` | Kinetic-energy spectra and peak-alignment scale separation | Guarded `paper-scale-separation` definitions cover `96` and `384` while reusing the standard `192` case; three dimensionless lower-panel `alignment_peak.cos_theta` curves are admitted; dimensional upper-spectrum transform, approved runs, and comparison remain open |
| Figure 12: heat-flux sensitivity | Heat-flux parameter variants | Spectra/transfer/diagnostic products under heat-flux changes | Two nonnominal active beta-10 guarded definitions exist and nominal active/passive cases reuse standard definitions; dimensionless top-panel curves are admitted through `alignment_peak.cos_theta`; the dimensional lower spectrum is excluded pending the F-038 unit transform; execution, steady nonlinear comparison, and exact interval/production budget remain missing |
| Figure 13: limiter-frequency sensitivity | Beta 100, Alfvenic, several `nu_lim` | Energy spectra, `beta Delta` PDFs, strain PDFs, pressure-stress transfer | Inputs and foundations exist; dimensionless panel-(b) `beta Delta` PDFs and panel-(d) normalized transfer curves are vector-extracted for comparison; panel (a) spectra and panel (c) strain are excluded pending the F-038 unit transform; case runs/comparisons are missing |

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
frontend shell because it does not provide `mpirun`. Historical
scheduler-capable Frontier Kokkos-Serial/MPI jobs `4659053` and `4659141`
executed the one-rank/four-rank AMR oracle through `srun`, passing at
`rtol = atol = 1.0e-12`. After the modal-driver merge, F-067 closes the
current retained gate explicitly: scheduled one-node job `4743666` runs both
MPI pytest files against a current Kokkos-Serial/MPI executable by translating
their `mpirun -np` invocations into nested `srun --overlap` steps. Both tests
pass, so the missing frontend launcher is not an unmet numerical-evidence
gate.

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
  cgl_lf_paper_standard_passive_alfvenic_beta100.athinput
  cgl_lf_paper_standard_active_random_beta10.athinput
  cgl_lf_paper_standard_active_random_beta100.athinput
  cgl_lf_paper_standard_passive_random_beta10.athinput
  cgl_lf_paper_standard_passive_random_beta100.athinput
  cgl_lf_paper_heat_flux_beta10_strong.athinput
  cgl_lf_paper_heat_flux_beta10_weak.athinput
  cgl_lf_paper_compressive_active_random_beta1.athinput
  cgl_lf_paper_compressive_active_random_beta100_sonic.athinput
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

**Status (2026-05-25):** All listed smoke, Figure 2(b)-complete standard
baseline, Figure 11 scale-separation, Figure 12 heat-flux, and beta-100
limiter-frequency input definitions are present. The
standard matrix now maps all eight active/passive, Alfvenic/random beta-10
and beta-100 histories in that panel, while retaining the active Alfvenic
beta-1 case needed elsewhere. Figure 12 reuses nominal active/passive
standard cases and adds the two active nonnominal heat-flux definitions.
Figure 3 reuses four standard cases and adds active random beta-1 plus
active random beta-100 sonic-correlation definitions through
`paper-compressive`.
Figure 11 adds active Alfvenic beta-10 definitions at `n_perp = 96` and
`384` while reusing the standard `n_perp = 192` case. The standard,
heat-flux, compressive, limiter, and scale-separation
decks now explicitly select rank-local binary/restart products with
`single_file_per_rank = true`, and their future
workflow manifests archive the output cadence/layout choices used for sizing
and the timestep/STS/CFL contract exposed by F-031. They remain production
definitions only; the F-033 hard-wall model has reduced-GPU qualification
under unlimited STS, but no paper-standard execution or statistical
qualification has occurred.

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
pressure-transfer partitions with a real-space closure check and the
MKS24-normalized `T_Delta p/T_total` ratio using
`T_total ~= E_K (2 pi u_rms/L_perp)`, selected-shell alignment PDFs, and
per-case ensemble averages inside the archived analysis window. It
additionally reports snapshot CGL pressure mechanical work as
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
as momentum. The analyzer can now validate external reference-data
checksums and uncertainties and render curve or sampled-surface comparison
overlays once qualified data are supplied. Populated
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
owned operator-face cap occupancy, plus restartable `lf_qprwrk` and
`lf_qpewrk` for signed RKL2-applied capped-face transport contractions.
Both diagnostics count a shared face once and use fine-side closure-face
ownership at coarse/fine interfaces. Inputs enabling
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
and generic plots are implemented. Figure 2(a)-coordinate joint
pressure-density PDFs, two-panel rendering, and optional sampled-surface
reference comparisons are now implemented. F-051 decodes the raster/colorbar
into donor-coordinate samples, and F-052 applies the source-qualified
`s = 0.5` pressure conversion and PDF Jacobian to emit admitted surfaces.
Figure 13(b)
dimensionless `beta Delta`
reference curves are populated; Figure 13(c) dimensional conversion remains
open under F-038. Figure 13(d)'s normalized transfer is addressed under
E4/F-039.

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
normalized density fluctuation, separate parallel/perpendicular/magnetic
pressures, `Delta p`, projected `Delta p` gradients, and local-field
parallel/perpendicular velocity gradients are implemented. The Figure 3
compressive product projects `u_k` onto the full Fourier wavevector before
`k_perp` binning. Synthetic checks now cover both the existing single mode and
exact zero compressive power for a transverse mode; the JSON records field
and normalization definitions and the plotter renders the new products.
Reference normalization qualification and authorized execution of the
guarded Figure 3 case matrix remain open. F-053 establishes that the
available manuscript definitions and donor plotted units are insufficient for
absolute ordinate conversion without donor diagnostic code, machine-readable
reference data, or an explicit discrete Fourier-normalization statement.

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
MKS24's stated `T_total ~= E_K (2 pi u_rms/L_perp)` normalization is emitted
as `pressure_transfer.transfer_normalized_by_total`, with exact synthetic
manifest comparison coverage. Figure 7 lower-panel and Figure 13(d)
dimensionless transfer references are populated under F-039; production-field
evaluation and scientific comparison remain open.

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
reference-comparable peak curves are implemented, with explicit
`paper-analyze --alignment-shells` selection for wide plotted ranges.
Synthetic output is checked for finite bounded histograms and exact
peak-curve comparison. Eight dimensionless Figure 9 peak-alignment references
are populated under F-040 for the standard case matrix; three dimensionless
Figure 11 lower-panel peak curves are populated under F-041 for the
scale-separation matrix. Physically interpretable production validation
remains open.

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
archives reduced model/forcing choices. `paper-standard`, `paper-nulim`,
`paper-heat-flux`, `paper-compressive`, and `paper-scale-separation` now
enumerate the defined
production inputs, retain
model/resolution/domain/analysis-window metadata and in-place
binary/checkpoint references, and refuse execution unless
`--authorize-paper-execution` is passed.
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
| `paper-heat-flux` | Figure 12 nonnominal heat-flux sensitivity suite; defined and guarded by explicit authorization | HPC |
| `paper-compressive` | Figure 3 compressive-case suite, including beta-1 and sonic-correlation definitions; defined and guarded by explicit authorization | HPC |
| `paper-scale-separation` | Figure 11 `n_perp = 96, 192, 384` suite; defined and guarded by explicit authorization | HPC |
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

The following user-specified constraints governed the archived debug
qualification ladder and remain the default for future debug-only
requalification unless revised explicitly:

```text
Project/account: AST207
Project run root: /lustre/orion/ast207/proj-shared/dfielding/CGL
Allowed QOS: debug only
Maximum testing budget: 1000 Frontier node-hours total across this effort
Concurrency: one debug-QOS job in any state at a time
```

The archived `E01-pre-modal-driver` Stage I protocol separately authorized
sequential paper production on `batch` with Frontier's default `normal` QOS.
The 2026-05-30 clean restart decision gave `E02-modal-driver` a fresh
incremental `4000` node-hour planning ceiling while retaining earlier use as
history. F-071 through F-075 then measured E02 R16/R02/R17 costs and
historically authorized sequential mapped work inside a `900.000000`
node-hour envelope. F-076 supersedes that authorization: preserve E02 as
pipeline and cost evidence only. F-078 closes corrected E03 qualification.
The reviewed E03 token is retained, `reconcile` passed, and fresh
`R02/s00_rankio_t0_t0p1` job `4745922` was submitted from `t = 0`, formally
inspected, and recorded `accepted` at exact `t = 0.1` for `0.188889`
node-hours. F-079/F-080 controller provenance is archived. Authenticated
`R02/s01_rankio_t0p1_t0p25` job `4746154` is also accepted through exact
`t = 0.25` for `0.289167` node-hours. Authenticated
`R02/s02_rankio_t0p25_t1` job `4746182` is accepted through exact `t = 1.0`
for `1.465556` node-hours. Authenticated `R02/s03_rankio_t1_t1p5` job
`4746356` is accepted through exact `t = 1.5` for `1.048611` node-hours,
bringing corrected E03 use to `2.992223`. F-083/F-084 controller hardening is
promoted at canonical revision `9689c269bf329542815a1b2b137881126964b05c`,
helper SHA-256
`1c633ebb58294938a0a0609742ae8f8d1f88cd15242ffe649578796edeb39375`.
Authenticated `R02/s04_rankio_t1p5_t2` job `4746435` is accepted through
exact `t = 2.0` for `1.063611` node-hours, bringing corrected E03 use to
`4.055834`. Retained corrected F-087 recost evidence JSON SHA-256 is
`c382b78ab9e21466648adf9d7ea38d2b407f0e986579f128bdfe6bc44bef79dd`.
Authenticated `R02/s05_rankio_t2_t2p5` job `4746663` is accepted through
exact `t = 2.5` for `1.062500` node-hours, bringing corrected E03 use to
`5.118334`. Retained corrected F-088 recost evidence JSON SHA-256 is
`39115585cf7ca866e3aade1e53c69aab6ee9217361d75db7946748905e1c0e5c`.
Authenticated `R02/s06_rankio_t2p5_t3` job `4746773` is accepted through
exact `t = 3.0` for `1.077500` node-hours, bringing corrected E03 use to
`6.195834`. Retained corrected F-089 recost evidence JSON SHA-256 is
`d2c5e27553426beb03c8c1882bd6473b682029c6ba67acd32bb06c5fb9a3ff01`.
Authenticated `R02/s07_rankio_t3_t3p5` job `4746953` is accepted through
exact `t = 3.5` for `1.151111` node-hours, bringing corrected E03 use to
`7.346945`. Retained corrected F-090 recost evidence JSON SHA-256 is
`468e787db48d98661cd3ffe125829a78f135bb7ca4c9104c766a00b693f586a3`.
Authenticated `R02/s08_rankio_t3p5_t4` job `4747015` is accepted through
exact `t = 4.0` for `1.200000` node-hours, bringing corrected E03 use to
`8.546945`. Retained corrected F-091 recost evidence JSON SHA-256 is
`dd5b0f5b82cf81005c1a481ee83d1127aa43861dcdd8c22765170e28b0ad6c99`.
Authenticated `R02/s09_rankio_t4_t4p5` job `4747087` is accepted through
exact `t = 4.5` for `1.208056` node-hours, bringing corrected E03 use to
`9.755001`. Retained corrected F-092 recost evidence JSON SHA-256 is
`8b10b1786db520740b687d989207f3cda6e0123f9bf2b62e75a5df3849511561`.
Authenticated `R02/s10_rankio_t4p5_t5` job `4747146` is accepted through
exact `t = 5.0` for `1.226389` node-hours, bringing corrected E03 use to
`10.981390`. Retained corrected F-093 recost evidence JSON SHA-256 is
`5256fb7fd9b75f175a16c2c16ee8a995b8d6d9281a92b70a979ea62171835bf2`.
Authenticated `R02/s11_rankio_t5_t5p5` job `4747202` is accepted through
exact `t = 5.5` for `1.268889` node-hours, bringing corrected E03 use to
`12.250279`. Retained corrected F-094 recost evidence JSON SHA-256 is
`323d7cb3cbe108eee944045c189cff75ffc07d06833c49bcfa8315f07fad8d51`.
Hardened reconciliation closes with `12/12/12` ledger
rows/manifests/reservations, no active reservation, and no transaction.
Archive and catalog the current committed controller state, reconcile, then
prepare only `R02/s12_rankio_t5p5_t5p75` on one node from the authenticated `s11` terminal
siblings with Slurm walltime `01:05:00`, Athena timeout `00:55:00`, and its
reviewed one-segment `2700`-second threshold retaining `600` seconds on both
timeout margins. The threshold is scoped only to `s12` and cannot ratchet
automatically. Any scientific, provenance, scheduler, storage, budget, or
reconciliation failure blocks successor preparation. Inspect, account,
reconcile, and recost before any further extension.

### Corrected E03 Successor-Segment Review Policy

After accepted `s11`, quarter-unit proposals use a reviewed timing envelope,
not an automatic ratchet. Let `E11` be the top-level Slurm allocation
`ElapsedRaw` in integer seconds. Retain the latest six developed half-unit
intervals and compute:

```text
rolling_six = [3879, 4144, 4320, 4349, 4415, E11]
ols_next_half = (24938 + 10 * E11) / 15
local_next_half = max(E11, 2 * E11 - 4415)
next_half = ceil(max(ols_next_half, local_next_half))
projected_quarter = ceil(next_half / 2)
```

The selected `s12` profile is admissible only when `projected_quarter <= 2700`.
For accepted `s11`, `E11 = 4568`, so `ols_next_half = 4707.866667`,
`local_next_half = 4721`, `next_half = 4721`, and
`projected_quarter = 2361`. The `2700`-second threshold is a reviewed
authorization and post-run recost boundary, not a controller-enforced kill
deadline. If an otherwise clean segment completes after its threshold but
before its Athena timeout, retain and inspect the endpoint and block automatic
extension pending recost. Any scientific, provenance, scheduler, storage,
budget, or reconciliation failure is a hard stop and does not authorize a
shorter segment.

All Frontier simulations for this CGL-LF project must place their output,
logs, manifests, restart products, and usage accounting beneath:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL
```

Do not write CGL-LF simulation products into the older `divb` directory or a
personal scratch path merely because a prior script used those locations.

This root is shared by multiple CGL campaigns. During the 2026-05-29 EDT
audit it contained exploratory beta-25 comparison runs in addition to MKS24
archives, including then-active debug job `4743106`. Slurm now records that
job as `COMPLETED 0:0`, but its top-level manifest remains stale at
`state = running`. Before any new MKS24 prepare or submit action:

1. Inspect `squeue -u "$USER"` without filtering only for MKS24 job names.
2. Inspect active campaign manifests beneath the CGL root.
3. Refuse to modify or prune any live campaign tree.
4. Avoid overlapping another root-writing CGL campaign unless an explicit
   review documents independent directories, accounting, and operational
   safety.
5. Record exploratory allocations in a separate ledger; they are neither
   debug-qualification evidence nor MKS24 Stage I evidence merely because
   they share the root.

### Debug-QOS Restriction and Its Consequence

The Frontier guide states that the `debug` QOS is for short non-production
debug tasks, that production work and job chaining using that QOS are
prohibited, that each user may have only one `debug` job in any state at a
time, and that requested walltime cannot exceed two hours.

Therefore:

- every debug-only qualification job submitted by agents under this plan must
  include:

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

This distinction matters: debug-queue qualification can establish
**production quality**, but paper production must use the separately reviewed
OLCF-compliant production path. The archived `E01-pre-modal-driver` campaign
used `batch` with default `normal` QOS. Reusing that scheduler path for
`E03-forcing-policy` requires the reviewed token, measured reservation,
clean reconciliation, and shared-root review recorded in this plan.

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

The live Frontier debug qualification ledger is governed by:

```text
1000 node-hours
```

Its current recorded use is `1.847784` node-hours, including historical
qualification plus corrected E03 `g024`--`g031`. The pre-replacement
historical component is `0.851670`; the clean-restart decisions do not erase
it.

The earlier `E02-modal-driver` reset used a fresh incremental project planning
ceiling:

```text
4000 node-hours
```

This block is historical and superseded by F-076. It does not authorize any
new E02 preparation or submission. The epoch-aware helper initially carried a separate `4.0`
node-hour production-path envelope for an `R16` pilot only, then a `7.0`
node-hour completion envelope after F-068. Accepted replacement-driver jobs
`4743735` through `4744158` now complete exact `t = 10.0` using `6.145556`
node-hours and pass production-window analysis. F-069 then authorized a
separate `2.0` node-hour R02 standard-layout timing pilot. Jobs `4744198` and
`4744205` complete that pilot through native cadence `t = 0.25` using
`0.473333` node-hours. R17 jobs `4744210` and `4744230` complete the
high-resolution timing pilot through authenticated continuation `t = 0.10`
using `4.235556` node-hours. F-071 replaces the bracket with a measured
`900.000000` node-hour mapped-matrix envelope. F-072 records R02 job `4744249`
through exact `t = 1.0`, updates
the projection to `702.195370` node-hours, and requires conservative
two-hour-or-shorter normal-QOS segments. F-073 records R02 job `4744518`
through exact `t = 1.5`, updates the projection to `725.464938` node-hours,
and permits only the next conservative segment through `t = 2.0`. F-074
records R02 job `4744913` through exact `t = 2.0`, updates the projection to
`730.081697` node-hours, and permits only the next conservative segment
through `t = 2.5`. F-075 records R02 job `4745305` through exact `t = 2.5`,
updates the projection to `728.164877` node-hours, and permits continued
frozen-matrix execution one inspected segment at a time, with R17 last. The
old `1068.888889` Stage I reservation is planning history, not an `E02`
entitlement. F-076 cancels job `4745498`, records E02 cumulative use
`15.628610`, and requires a forcing-correct `E03-forcing-policy` epoch from
`t = 0`.

According to OLCF accounting, a job is charged according to requested nodes
multiplied by actual time from entering the running state until exit; unused
cores within requested nodes do not reduce the charge, and failed jobs still
consume allocation. Historical, exploratory, debug-qualification, and MKS24
production ledgers must remain separate while still being reconcilable to
actual project usage.

#### Before Every Submission

An agent must:

1. Inspect the relevant epoch ledger, the historical ledgers, and the
   shared-root campaign records.
2. Calculate the conservative reservation for the proposed job:

   ```text
   reserved_node_hours = requested_nodes * requested_walltime_hours
   ```

3. Refuse every new E02 preparation or submission. For E03, require the
   reviewed qualification-approval token and refuse submission if the new
   epoch reservation would exceed its currently approved envelope or the
   fresh project ceiling:

   ```text
   e03_consumed_actual + e03_active_reserved + proposed_reserved > approved_e03_envelope
   e03_consumed_actual + e03_active_reserved + proposed_reserved > 4000
   ```

4. Confirm there is no other queued user job in `PD`, `R`, or any other state
   before issuing `sbatch`.
5. Confirm there is no unrelated live root-writing CGL campaign unless an
   explicit concurrency review permits overlap.
6. Refuse successor preparation after any scientific, provenance, scheduler,
   storage, budget, or reconciliation failure. A shorter segment is a reviewed
   performance response, not a recovery path for invalid evidence.
7. Record the proposed job, scientific purpose, acceptance criterion, node
   count, time limit, and reservation in the campaign notes before it runs.

Recommended pre-submission commands:

```bash
squeue -u "$USER" -o "%.18i %.12P %.30j %.8u %.2t %.10M %.6D %R"
tail -n 20 /lustre/orion/ast207/proj-shared/dfielding/CGL/accounting/frontier_node_hours.csv
tail -n 20 /lustre/orion/ast207/proj-shared/dfielding/CGL/accounting/mks24_stage_i_E03_forcing_policy_node_hours.csv
test -f /lustre/orion/ast207/proj-shared/dfielding/CGL/accounting/mks24_stage_i_E03_forcing_policy_qualification_approval.json
python3 /autofs/nccs-svm1_home2/dfielding/athenak-df/scripts/frontier/cgl_lf_stage_i.py \
  --root /lustre/orion/ast207/proj-shared/dfielding/CGL reconcile
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

Use `showusage` as an external allocation cross-check when appropriate. The
live debug ledger remains authoritative for its `1000` node-hour testing
constraint; the E03 epoch ledger must be authoritative for mapped-production
authorization.

#### Archived Initial Debug Budget Envelope

The following table governed the archived debug ladder. It remains useful as
historical rationale but does not allocate `E02` spending:

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

#### Archived Fresh `E02` Initial Envelope

The superseded E02 epoch started narrowly:

| Activity | Archived initial `E02` upper envelope | Historical expansion rule |
| --- | ---: | --- |
| Replacement-driver local and Frontier debug requalification | `50` node-hours | Stop after the smallest tests establish modal restart, decomposition, hard-wall, and rank-local I/O behavior |
| Completed Stage I `R16` | `6.145556` node-hours actual inside `7.0` reserved | F-069 completed low-resolution case and production-window analysis |
| Completed Stage I `R02` standard-layout timing pilot | `0.473333` node-hours actual inside `2.0` reserved | F-070 completes native-cadence standard-layout timing and storage measurement |
| Completed Stage I `R17` high-resolution timing pilot | `4.235556` node-hours actual | F-071 completes fresh and authenticated developed eight-node timing measurements |
| Stage I frozen mapped matrix | `900.000000` node-hours | F-071 continuation-aware projection is `697.019444`; historical E02 authorization only, superseded by F-076 |
| Accepted Stage I `R02` mapped continuation through `t = 1.0` | `1.452778` node-hours actual | F-072 longer developed interval updates mapped projection to `702.195370`; continue conservatively |
| Accepted Stage I `R02` mapped continuation from `t = 1.0` through `t = 1.5` | `1.052222` node-hours actual | F-073 late-time developed interval updates mapped projection to `725.464938`; continue conservatively through `t = 2.0` |
| Accepted Stage I `R02` mapped continuation from `t = 1.5` through `t = 2.0` | `1.068889` node-hours actual | F-074 late-time developed interval updates mapped projection to `730.081697`; continue conservatively through `t = 2.5` |
| Accepted Stage I `R02` mapped continuation from `t = 2.0` through `t = 2.5` | `1.061944` node-hours actual | F-075 stable late-time developed interval updates mapped projection to `728.164877`; continue only the frozen matrix one inspected segment at a time |
| Stage II extension | `0` | Remains unauthorized until Stage I acceptance and recosting |

#### Recommended Accounting Helper

A tracked implementation is now provided at:

```text
scripts/frontier/cgl_lf_frontier.py
scripts/frontier/build_cgl_lf_frontier.sh
```

`cgl_lf_frontier.py` is intentionally a preparation/accounting utility, not a
submission launcher. It requires the prescribed CGL root for real use,
reserves conservative node-hours, rejects walltimes greater than two hours
and paper-standard/limiter/heat-flux/compressive/scale-separation production
decks under
`debug`, emits an
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
Real preparations may use a clean source checkout outside the campaign root,
but they must pass `--source-bundle` naming a retained Git bundle beneath
`$CGL_ROOT`. The helper clones that bundle transiently, verifies that it
contains the prepared revision, records the bundle checksum, and archives a
complete sibling set when a rank-local restart is supplied. Its submission
preflight rejects any queued user job and requires explicit acknowledgement
of an active top-level CGL-root campaign record only after an isolation review.
Preparation also parses the selected input deck and rejects any malformed or
absent `block/name=value` override target before reserving node-hours.
This utility remains the debug-qualification helper. Production work uses the
separate `scripts/frontier/cgl_lf_stage_i.py` utility. That production helper
now carries explicit `E03-forcing-policy` run/reservation/summary/ledger paths,
requires the reviewed qualification token, rejects cross-epoch restart and
bundle ancestry, authenticates explicit restart physical-time markers, uses
an all-user-job queue preflight and active top-level CGL-root campaign-record
checks, submits atomically under lock, and retains the F-071 measured
`900.000000` node-hour frozen-matrix envelope only as an E03 planning prior.
Preserve the committed F-060/F-061 controls, keep E02 fail closed, and
re-estimate later segment sizes from corrected E03 measurements. F-072's
normal-QOS preparation rejection above two hours remains active before
restart archival.
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
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}
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
  sequential-retention allocation with 20 percent margin for the pre-F-034
  eight-definition inventory; F-035 applies those measured per-file sizes to
  the expanded thirteen-unique-run baseline proposal, yielding `403.162` GB
  sequential storage with 20 percent margin. F-041 subsequently extends the
  proposed matrix to fifteen unique runs and estimates `781.694` GB
  sequential storage with 20 percent margin after adding the mixed-resolution
  Figure 11 cases. F-049 adds two standard-layout Figure 3 definitions and
  revises the candidate to seventeen unique runs with an estimated
  `880.368` GB sequential allocation with margin. G011 itself explicitly withholds
  runtime/node-hour projection.
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
  corrected Frontier GPU execution is recorded by G013 below.
- G013 F-033 corrected nonlinear validation job `4673404` uses immutable
  revision `9e1986dfebfadfae0fcb4dedd7fa15326e6ae4c0`, HIP/MPI executable
  SHA-256
  `f9857fd48b1e07da7eecc1a03d8fbe335089237ee435a6469c5c2c9ab941f4d7`,
  and a reduced G012-derived input whose only model change is
  `limiter_hardwall = true` (input SHA-256
  `443b52f67a44cf49b74a248222b73db6e0de3dc62d2e79c450c26342763f2fc1`).
  It completes the one-node/eight-GPU unlimited-RKL2 run through `t = 2.0`
  in `2080` allocated seconds. Its 101 retained MHD history rows are finite
  with zero maxima for `lf_dfloor`, `lf_pfloor`, `lf_nonfin`, `lf_nonpos`,
  and `lf_hardbd`, terminal `lf_hwproj = 14763255866`, and forcing-work
  relative residual `9.627e-14`. Nine shared binary snapshots and five
  shared restarts are retained; initial/terminal snapshots are finite and
  readable through the implemented binary reader. Full spectral/strain
  analysis was deliberately not claimed for this correctness gate after an
  exploratory FFT/eigendecomposition pass was interrupted. Authoritative
  evidence is retained at
  `analysis/g013_f033_hardwall_gpu_qualification_evidence.json` with SHA-256
  `2381e89a24a4ff9c2d81e4e0e91d049c3b62ac0577f044fee1816e2ccd2d6108`.
- F-034 Figure 2(b) reference extraction is retained outside git under
  `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/digitized_fig2b_v1`.
  The tracked extractor refuses a source figure not matching SHA-256
  `828287244aa0ec72aad6e802268fbcb8d77af65f11b03fcf9634eb9cf741ab0e`
  and emits all eight visible active/passive threshold-volume histories with
  absolute digitization uncertainty `0.0025`; six vertices clipped above the
  panel range are explicitly omitted. `curves.json` has SHA-256
  `66996af7ddd48a39c1a8fc46e1fd7ba6206ced5625c2d30a123cc3df67df3df1`
  and validates through the implemented reference-manifest reader. These are
  published-panel references only, not AthenaK comparison results.
- F-036 Figure 12 reference extraction is retained outside git under
  `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/digitized_fig12_v1`.
  The tracked extractor refuses a source figure not matching SHA-256
  `aeeb79fef6a24de93d6fa44d2a6c62a254bcb65c2651dadaa5e95d2204c7001b`
  and emits four `alignment_peak.cos_theta` and four
  `spectra.grad_parallel_delta_p` curves with tick-calibrated vector
  coordinates, omitting 16 x-clipped endpoints and recording 30 recommended
  alignment shells. `curves.json` has SHA-256
  `bcdb17cbf2396c7e373c83fd0bef49a4b8a92498d740deba4f3bf7229e3260e0`;
  `digitization_metadata.json` has SHA-256
  `bd902377c45183f507a91cebf8eb86afd115a74bc69762f38475811f7c521ef6`.
  This provisional manifest predates F-038 and must not be used for
  dimensional lower-panel comparison.
- F-037 Figure 13 partial reference extraction is retained outside git under
  `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/digitized_fig13_v1`.
  The tracked extractor verifies all four panel PDFs and emits seven directly
  comparable curves: three panel-(a) `spectra.velocity` series and four
  panel-(b) `pdf.beta_delta` series. It records five-percent relative
  uncertainty and omits seven panel-(a) and twelve panel-(b) boundary-clipped
  vertices. `curves.json` has SHA-256
  `410e208f968f6e3b1007baff9e97ae8aeb80b03f2e3039b819fae7dc2e9494eb`;
  `digitization_metadata.json` has SHA-256
  `7070a0b5cfbb0e303e4d4655b7733fab00618b4d091f063c067132c75beaa582`.
  This provisional manifest predates F-038 and must not be used for
  dimensional panel-(a) comparison.
- F-038 corrected dimensionless-only reference manifests are retained outside
  git under `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/`.
  `digitized_fig12_v2/curves.json` emits four
  `alignment_peak.cos_theta` curves and has SHA-256
  `8bc302c082ed971eda657f5aac283fd423ab182dbd676faa7da72ca3436a0b9c`;
  its `digitization_metadata.json` has SHA-256
  `8dc1d61bc62af1ea11e386a756a27e932f6168805bf2bdd4e71f3c05e5e9d5e8`.
  `digitized_fig13_v2/curves.json` emits four `pdf.beta_delta` curves and
  has SHA-256
  `816013b3e4af9f5755fde5d64228721c86f14206a69210e3b2c38db35c1aa12e`;
  its `digitization_metadata.json` has SHA-256
  `67be32ec6fcac377541cf6345f6b913766cfd0a8f7179bdd74043f2d50685a46`.
  The `v2` metadata records Figure 12's dimensional lower spectrum and
  Figure 13(a),(c),(d) as excluded at the F-038 admission boundary; F-039
  supersedes that exclusion for the dimensionless panel-(d) transfer ratio
  only. These manifests are reference data, not simulation outcomes.
- F-039 dimensionless transfer reference manifests are retained outside git
  under `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/`.
  MKS24 defines the plotted transfer normalization in `source/MKS24.tex:591`
  as `T_total ~= E_K (2 pi u_rms/L_perp)`. The tracked Figure 7 extractor
  verifies `source/fig7.pdf` SHA-256
  `5411da17391d19fceb26f617f91b4fb1153c7638ccca8afee6c8e0d013651afa`
  and `digitized_fig7_v1/curves.json` emits three lower-panel
  `pressure_transfer.transfer_normalized_by_total` curves with SHA-256
  `7c2f02bca4725674b9e77292a7f33c40de211bdcd48fe2e29a4b2e8fe81c1d09`;
  its metadata SHA-256 is
  `7ed7869b9e1bf43556dff2d97a05647a698e092592cae8d2087bd55296e269cc`.
  The updated Figure 13 extractor verifies `source/fig13d.pdf` SHA-256
  `bf1c063fa5a54c2dbbe9dce417490b5bf821fcd2e1cb1926d78c923b90387a24`
  and `digitized_fig13_v3/curves.json` emits four `pdf.beta_delta` plus four
  normalized-transfer curves with SHA-256
  `78f084d55d851561bf7a649becbfa907626f0890605b5e792fe0352cc205c924`;
  its metadata SHA-256 is
  `24f5f6f771fbb894ef2c62d4d75663345686499315116f5e573ede2c4ef690d5`.
  Both extractors record absolute transfer ordinate uncertainty `0.025`;
  dimensional Figure 7 upper-panel and Figure 13(a),(c) comparison remains
  excluded pending the paper-to-AthenaK transform.
- F-040 Figure 9 dimensionless alignment reference extraction is retained
  outside git under
  `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/digitized_fig9_v1/`.
  The tracked extractor verifies `source/fig9.pdf` SHA-256
  `9d0dbb05514a8ce7a48eddb178b05a47ec512e58d2a2ff99d899123f8e458984`,
  emits eight `alignment_peak.cos_theta` curves for the guarded standard
  active/passive Alfvenic/random beta-10/beta-100 cases, omits 16
  horizontal-boundary vertices, and records 27 recommended alignment shells
  with absolute ordinate uncertainty `0.005`. `curves.json` has SHA-256
  `07417294ca44e7aa7e8c499931e12ed7fce9eec55e4691e3172a8684470163b8`;
  `digitization_metadata.json` has SHA-256
  `55db154b26d3b040258f79c6ee56042fe0e215c8f6ea9a39fc9155ecf177d242`.
  Manifest ingestion admits all eight curves; these are references only, not
  simulation comparisons.
- F-041 Figure 11 dimensionless lower-panel alignment extraction is retained
  outside git under
  `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/digitized_fig11_v1/`.
  The tracked extractor verifies `source/fig11.pdf` SHA-256
  `efab7ae64d5837ee621249acfedf79cd77bb221e37504f8c41a79bf5ba30f5f8`,
  emits three `alignment_peak.cos_theta` curves for active Alfvenic beta-10
  cases at `n_perp = 96`, `192`, and `384`, omits four horizontal-boundary
  vertices, and records 56 recommended alignment shells with absolute
  ordinate uncertainty `0.005`. `curves.json` has SHA-256
  `9922eb24245cdf0e7dfe355864a610ee8ebf437ccf593ea40a2e7df6016c037a`;
  `digitization_metadata.json` has SHA-256
  `b59b0db90ec786724095ab645ec0e84f915c5abbb3db97b2bd379c9745d0665c`.
  The dimensional upper-panel kinetic spectra remain excluded under F-038;
  these lower-panel curves are references only, not simulation comparisons.
- F-043 Figure 8 selected-shell alignment-PDF extraction is retained outside
  git under
  `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/digitized_fig8_v1/`.
  The tracked extractor verifies `source/fig8.pdf` SHA-256
  `0e3e292406bc27ce24b2ef449bd95cb8addf240f02b5ac60f60286459a299d8b`,
  decodes its embedded RGB heatmaps using the plotted linear colorbar,
  and emits ten `alignment.<shell>` curves for active/passive Alfvenic
  beta-10 cases at shells `8,16,32,64,128`. Its maximum decoded
  pre-renormalization unit-integral error is `0.0264475132735984`, beneath
  the declared `0.03` tolerance, with absolute PDF-density uncertainty
  `0.05`. `curves.json` has SHA-256
  `a533778f7777d755a682e898ae44dcb7356be8f74a9818d6b6cc7b6fd6048522`;
  `digitization_metadata.json` has SHA-256
  `25933691ade9051bd37c9830405005437279171761d5974d5ece6c43e9bc38a3`.
  Manifest ingestion admits all ten curves; these are references only, not
  simulation comparisons.
- F-044 Figure 4(a) normalized-density PDF extraction is retained outside git
  under
  `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/digitized_fig4a_v1/`.
  The tracked extractor verifies `source/fig4a.pdf` SHA-256
  `18a5c7119f2cb316e0283b354d25479d006aa3f27df33031035b16810937c024`,
  maps the plotted `rho/<rho>` coordinate exactly onto
  `pdf.density_fluctuation` through `x -> x - 1`, and emits eight standard-
  case curves with five-percent relative ordinate uncertainty. It records
  38 omitted plot-boundary vertices. `curves.json` has SHA-256
  `00c0fadeb15a82758be543a95b45d883ebea32388961756d0146eef3f8ca4cb9`;
  `digitization_metadata.json` has SHA-256
  `c04f32ddc9c65fecd387bc3e7a81a58f384b69f6a406dbf4c2f142f1b22d25d0`.
  Manifest ingestion admits all eight curves; these are references only, not
  simulation comparisons.
- F-047 adds the missing Figure 2(a)-coordinate product in tracked analysis:
  `pressure_density_joint.parallel` and `.perpendicular` retain common-bin
  two-dimensional PDFs in `<p> delta rho/<rho>` versus `delta p` coordinates,
  and the plotter emits case-specific two-panel diagnostics with the `5/3`
  guide. The pinned `source/fig2a.pdf` SHA-256 is
  `277545ff1a86b5ca46e97fa034ce189f72f149ecbb515470c30996494c1de129`.
  Its pressure-scaled raster is not admitted reference data until a reviewed
  donor-pressure transform exists.
- F-050 extends the tracked reference-data interface to optional sampled
  surfaces for `pressure_density_joint.parallel` and `.perpendicular`, with
  data/source checksums, positive `z_uncertainty`, bilinear interpolation,
  and residual rendering. An exact synthetic surface fixture validates the
  route; it does not by itself convert or admit the Figure 2(a) raster.
- F-051 Figure 2(a) raw surface extraction is retained outside git under
  `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/digitized_fig2a_v1/`.
  The tracked extractor verifies the pinned PDF checksum, reconstructs its
  `1406 x 1250` tiled raster, decodes the 256-color labeled bar, and rejects
  panel-overlay pixels other than the two observed black RGB values. It emits
  `6953` parallel and `6952` perpendicular donor-pressure-coordinate samples
  with absolute density uncertainty `0.5`. `excluded_surfaces.json` has
  SHA-256 `3c2c4e8984d538acda37c6c36c49dfa11eae73a7ecadd03a05074b5660cbc428`;
  `digitization_metadata.json` has SHA-256
  `5354a7206f80a0630e9f86db1dd3a9c5efd552c59e2cc2a8341bf78ebd89f541`.
  This output is intentionally not an analyzer manifest: a qualified pressure
  scale `s` would require `x_A = s x`, `y_A = s y`, and
  `z_A = z/s^2`, with the same `1/s^2` scaling for `z_uncertainty`,
  before comparison.
- F-052 supersedes that exclusion for Figure 2(a) only. The pinned
  `source/MKS24.tex` SHA-256
  `6d6e748fd1883c5d33167be653d67aed2f84a9b364267d9e90ea075df184af4c`
  states fixed `B0` and `p0` proportional to `beta0`, and reports
  beta-100 `p0 = 100` in code units; therefore the beta-10 panel has donor
  `p0 = 10`, while the matched AthenaK input uses `p0 = 5`.
  `build-cgl-implementation/cgl_lf_reference/arXiv-2405.02418v2/digitized_fig2a_v2/`
  retains the raw donor samples and admits transformed samples with
  `s = 0.5`; `surfaces.json` has SHA-256
  `3e2f961b1bdc8bb46180be1e4b38448a997d306738a834898caf3b5b0bd89845`,
  `donor_surfaces.json` has SHA-256
  `2430af5c3d3eaf025c200f2e7ef4ab719b1053e9e01437ca1ee91918c235dd5e`,
  and `digitization_metadata.json` has SHA-256
  `89575ae0dd3953a72c903fca377c3709923e93199b522fd2bab4178b0a862016`.
  The admitted parallel and perpendicular surfaces retain `6953` and `6952`
  samples respectively, with transformed absolute density uncertainty `2.0`.
  This does not qualify unrelated spectral or strain ordinate transforms.
- F-053 audits the remaining spectral/strain boundary. MKS24
  `source/MKS24.tex:478-496` defines shell-binned spectra and transfer
  filtering; the cited donor source
  `AthenaMagnetoimmutability_JPP.tex:976-1002` defines corresponding spectra,
  rate-of-strain fields, and plotted `phi^2 L_perp` units. Neither passage
  states the discrete Fourier normalization needed for absolute ordinate
  conversion. A targeted public publication/archive/repository audit on
  2026-05-25 found article/preprint records but no linked donor diagnostic
  code or numeric panel-data package. A 2026-05-26 follow-up checked the
  official journal-page text and repeated exact arXiv-ID Zenodo and
  exact-title GitHub repository searches; retrieved page text contained no
  `Data availability` heading and the exact queries returned no matching
  numeric-data/donor-code record.
  No further manifest was emitted.
- F-048 adds the missing spectral analysis fields in tracked analysis:
  `spectra.compressive_velocity` implements the Figure 3 full-wavevector
  projection, `spectra.density_fluctuation` provides a normalized Figure
  4(b)-oriented density field, and `spectra.p_parallel`, `.p_perp`, and
  `.magnetic_pressure` provide Figure 6(a)-oriented fields. Their emitted
  definitions and synthetic transverse-mode oracle do not qualify raw
  published spectral ordinates. F-049 now defines the previously missing
  sonic-correlation and beta-1 random Figure 3 production cases behind an
  explicit authorization guard.
- The retained ledger records `0.851670` node-hours used with no active
  reservation, including `0.188889` node-hours for failed G012 and
  `0.577778` node-hours for passed G013. Paper-scale
  execution, remaining unmatched/dimensional reference transforms, measured
  standard runtime, approved production execution, and scientific comparison
  remain blocked or open as identified above.

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
emergency-bound crossing before its measurement window. G013 establishes
corrected reduced-GPU behavior for the constrained F-033 model through
`t = 2.0`, with `2080` allocated seconds and retained MPI-I/O measurements,
but its half-linear-resolution correctness tier is not by itself a
standard-matrix production cost model. Do not infer a production rate from
failed G012; use G013 only as input to a separately reviewed standard-matrix
runtime/node-hour proposal.

### Archived Shared-I/O Storage Envelope and Archived `E02` Startup Sizing

Job `4661434` is a startup-only, one-node/eight-GPU measurement at the exact
`192 x 192 x 384` mesh and `32 x 32 x 64` MeshBlock layout used by the
guarded standard-layout
`paper-standard`/`paper-nulim`/`paper-heat-flux`/`paper-compressive`
definitions. It reaches only
`t = 0.01`; it validates historical standard-layout launch, host-memory
evidence, shared MPI-I/O file size/timing, and strict diagnostics, not
scientific production or runtime/node-hour cost. F-034 expanded the standard definitions to cover the
published Figure 2(b) series; F-035 adds the two nonnominal Figure 12
heat-flux cases. F-049 adds two guarded standard-layout Figure 3
compressive-case definitions. F-041 adds two guarded Figure 11
scale-separation inputs at
`96 x 96 x 192` and `384 x 384 x 768`; their file sizes below are estimated
from cell count until an approved measurement is retained. The standard active Alfvenic beta-100 hard-wall definition
and the limiter-scan hard-wall definition are model-identical apart from
labels/comments, so the proposed campaign reuses one accepted result in both
analysis roles. Figure 11 similarly reuses the accepted standard active
Alfvenic beta-10 `n_perp = 192` result.  The later frozen Stage I manifest
also excludes the unmapped active-Alfvenic beta-1 definition, yielding sixteen
mapped simulations: fourteen standard-layout cases plus `n_perp = 96` and
`384` scale-separation cases.

Historical E02 replacement-driver `g023` (job `4743662`) supplies a retained
startup-only rank-local standard-layout counterpart. Corrected E03 `g031`
(job `4745890`) now launches the exact `192 x 192 x 384` layout with 27
meshblocks per GPU rank, reaches `t = 0.01` in `77` allocated seconds with
zero strict counters, retains terminal aggregate snapshot size `509685800`
bytes, and retains modal checkpoint-group size `1381369688` bytes. This
closes corrected E03 startup, one-node memory-fit, and ranked retained-output
sizing only; it is not a production runtime or cadence measurement. Fresh
historical E02 R16 jobs `4743735` through `4744158`
complete rank-local production through `t = 10.0`; F-069 records
`4686854064` raw retained segment-tree bytes (`4.365` GiB), a deduplicated
logical cadence of `4514126544` bytes, and `6.145556` actual node-hours.
F-070 records R02 timing through native cadence `t = 0.25`, measured terminal
snapshot and restart groups of `509682424` and `1381366304` bytes, and
`0.473333` actual node-hours. F-071 records R17 timing through authenticated
continuation `t = 0.10`, measured maximum ranked snapshot and restart groups
of `4077462784` and `11052866880` bytes, and `4.235556` actual node-hours.

The following G011 calculations describe the archived shared-file assumption.
The tracked guarded paper decks now explicitly set
`single_file_per_rank = true` for `bin` and `rst` outputs, and the replacement
driver adds authoritative modal forcing state to the restart contract.
The common cadence still writes binary snapshots every `0.25` and restarts
every `1.0` through `tlim = 10.0`. Retained reduced runs confirm inclusive
initial/final time-cadenced output, giving 41 binary snapshots and 11 restart
points per executed case. Using the larger of the two exact historical G011
shared files of each kind for the fourteen mapped standard-layout runs, and
cell-count-scaled estimates for the additional one-eighth-size
`n_perp = 96` and eight-times-size `n_perp = 384` cases:

| Quantity | Calculation basis | Bytes | Decimal GB |
| --- | --- | ---: | ---: |
| One shared binary snapshot | G011 maximum measured file | `509632020` | `0.510` |
| One shared restart file | G011 maximum measured file | `1838131787` | `1.838` |
| Raw generated output for one standard-layout case | `41 bin + 11 rst` | `41114362477` | `41.114` |
| Estimated raw output for added `n_perp = 96` case | `one standard case / 8` | `~5139295310` | `~5.139` |
| Estimated raw output for added `n_perp = 384` case | `8 x one standard case` | `~328914899816` | `~328.915` |
| Estimated raw output for sixteen mapped cases | `14 standard + n96 + n384` | `~909655269804` | `~909.655` |
| Estimated retained output after acceptance | all snapshots plus final two restarts per case | `~543637277717` | `~543.637` |
| Estimated sequential peak before pruning accepted checkpoints | retain 14 standard and n96; run n384 last | `~675982766381` | `~675.983` |
| Estimated sequential retained-policy allocation with 20 percent margin | `1.2 x sequential peak` | `~811179319658` | `~811.179` |
| Estimated no-pruning allocation with 20 percent margin | `1.2 x all raw output` | `~1091586323765` | `~1091.586` |

These are archived shared-file lower-envelope estimates for the sixteen
mapped executions. The `E01-pre-modal-driver` Stage I protocol conservatively
retained its earlier `880.368` GB sequential-retention and `1140.924` GB
no-pruning allocations while a segment was submitted; that retained margin
did not authorize execution of the excluded beta-1 inventory definition.
Neither value is an `E02-modal-driver` storage reservation.

The historical candidate retention policy is to run cases sequentially,
retain all analysis snapshots and the final two restarts for every accepted
case, and remove earlier checkpoints only after acceptance and archive
verification are recorded. The standard-layout file-size measurement remains
authoritative in `analysis/g011_standard_layout_mpiio_sizing_evidence.json`
beneath the G011 run directory as shared-file reconnaissance. F-069 measures
completed rank-local R16 products: terminal snapshot and restart groups are
`63758976` and `172728048` bytes, respectively, and the seven accepted
segments retain `4686854064` raw bytes. At the same cadence, G023 terminal
groups project `36091999664` bytes (`33.613` GiB) for one standard-layout
case. F-070 replaces that startup-only projection with measured R02 terminal
groups: `509682424` snapshot bytes and `1381366304` restart bytes project
`36092008728` bytes (`33.613` GiB) per standard case. F-071 replaces
cell-scaled R17 storage with measured ranked groups. The mapped matrix projects
`798559758560` logical raw bytes (`743.717` GiB), `523477249904` accepted
retention bytes, and a sequential peak of `622953051824` bytes when R17 runs
last. Reserve `747543662189` bytes after a 20-percent sequential-retention
margin; the no-pruning alternative with the same margin is `958271710272`
bytes.
Executing all eighteen guarded definitions independently would duplicate the
equivalent hard-wall beta-100 role and execute the unmapped active-Alfvenic
beta-1 inventory entry; the frozen Stage I manifest does neither.

### Archived Candidate Standard-Run Cost Proposal

The only archived corrected nonlinear late-time GPU timing is G013: one node/eight
GPUs advances the `96 x 96 x 192` hard-wall beta-10 case through `t = 2.0` in
`2080` allocated seconds. A planning estimate for the standard mesh applies
an eightfold cell-count factor and a fivefold simulated-duration factor:

| Quantity | Planning calculation | Estimate |
| --- | --- | ---: |
| Standard `t = 10` time per unique case on one node | `2080 s x 8 x 5` | `83200 s = 23.111111` node-hours |
| Fourteen mapped standard/limiter/heat-flux/compressive cases | `14 x 83200 s` | `323.555556` node-hours |
| Added `n_perp = 96` scale-separation case | `83200 s / 8` | `2.888889` node-hours |
| Added `n_perp = 384` scale-separation case | `8 x 83200 s` | `184.888889` node-hours |
| Sixteen mapped cases including Figure 3 and Figure 11 | `323.555556 + 2.888889 + 184.888889` | `511.333333` node-hours |
| Recalculated twofold mapped-case planning estimate | `2 x 511.333333` | `1022.666667` node-hours |
| Conservative Stage I reservation retained by the submitted manifest | frozen protocol ceiling | `1068.888889` node-hours |
| Existing recorded debug qualification plus retained reservation | `0.851670 + 1068.888889` | `1069.740559` node-hours |

F-070 supplies the current replacement-driver standard-layout measurement.
Authenticated R02 continuation from `t = 0.1` to `0.25` consumes `0.285000`
node-hours, or `1.900000` node-hours per simulated time unit. It projects
`19.000000` node-hours per standard `t = 10` case and `266.000000` node-hours
for the fourteen mapped standard-layout cases. Keeping R17 bracketed at
`8x`--`16x` the R02 node-hours yields a provisional mapped-matrix range of
`424.145556`--`576.145556` node-hours including completed R16. This remains a
historical F-070 planning bracket.

F-071 supplies the developed high-resolution measurement. Authenticated R17
continuation from `t = 0.05` to `0.10` consumes `2.124444` node-hours, or
`42.488889` node-hours per simulated time unit, projecting `424.888889`
node-hours through high-resolution `t = 10`. Reusing completed R02 and R17
prefixes yields a continuation-aware mapped E02 projection of `697.019444`
node-hours. The measured `900.000000` node-hour envelope leaves `202.980556`
node-hours of projected margin for case-to-case and late-time variation.

F-072 extends R02 from `t = 0.25` to native restart boundary `t = 1.0`.
Job `4744249` consumes `1.452778` node-hours, or `1.937037` node-hours per
simulated time unit over the longer developed interval. Applying that rate to
the unfinished standard-layout work updates the continuation-aware mapped
projection to `702.195370` node-hours and leaves `197.804630` node-hours of
F-071 envelope margin. Continue with conservative normal-QOS segments at or
below two requested hours and recost after inspection.

F-073 extends the inspected R02 lineage from `t = 1.0` to native restart
boundary `t = 1.5`. Job `4744518` consumes `1.052222` node-hours, or
`2.104444` node-hours per simulated time unit over the late-time interval.
Applying that rate conservatively to the unfinished standard-layout work
updates the continuation-aware mapped projection to `725.464938` node-hours
and leaves `174.535062` node-hours of F-071 envelope margin. Continue only
from `t = 1.5` to `2.0`, then inspect and recost before lengthening later
segments.

F-074 extends the inspected R02 lineage from `t = 1.5` to native restart
boundary `t = 2.0`. Job `4744913` consumes `1.068889` node-hours, or
`2.137778` node-hours per simulated time unit over the latest late-time
interval. Applying that rate conservatively to the unfinished standard-layout
work updates the continuation-aware mapped projection to `730.081697`
node-hours and leaves `169.918303` node-hours of F-071 envelope margin.
Continue only from `t = 2.0` to `2.5`, then inspect and recost before
lengthening later segments.

F-075 extends the inspected R02 lineage from `t = 2.0` to native restart
boundary `t = 2.5`. Job `4745305` consumes `1.061944` node-hours, or
`2.123888` node-hours per simulated time unit over the latest late-time
interval. Applying that rate conservatively to the unfinished standard-layout
work updates the continuation-aware mapped projection to `728.164877`
node-hours and leaves `171.835123` node-hours of F-071 envelope margin.
F-076 supersedes the continuation authorization after the source-to-code
forcing audit. Prohibit every new E02 preparation or submission.

This archived E02 estimate is not a corrected standard-mesh benchmark. It
assumes linear scaling with cell count and simulated duration, equivalent
timestep behavior among cases, and no late-time performance change. The 2x
factor is a historical planning allowance, not a measured error bar. The
active compliant E03 proposal is narrower: the reviewed token is retained,
`reconcile` passed, and fresh `R02/s00_rankio_t0_t0p1` job `4745922` is
formally inspected and recorded `accepted` through exact `t = 0.1` for
`0.188889` node-hours. F-079/F-080 controller provenance is archived.
Authenticated `R02/s01_rankio_t0p1_t0p25` job `4746154` is also accepted
through exact `t = 0.25`. Authenticated `R02/s02_rankio_t0p25_t1` job
`4746182` is also accepted through exact `t = 1.0`. Authenticated
`R02/s03_rankio_t1_t1p5` job `4746356` is accepted through exact `t = 1.5`;
authenticated `R02/s04_rankio_t1p5_t2` job `4746435` is accepted through
exact `t = 2.0`; authenticated `R02/s05_rankio_t2_t2p5` job `4746663` is
accepted through exact `t = 2.5`;
authenticated `R02/s06_rankio_t2p5_t3` job `4746773` is accepted through
exact `t = 3.0`;
authenticated `R02/s07_rankio_t3_t3p5` job `4746953` is accepted through
exact `t = 3.5`;
authenticated `R02/s08_rankio_t3p5_t4` job `4747015` is accepted through
exact `t = 4.0`;
authenticated `R02/s09_rankio_t4_t4p5` job `4747087` is accepted through
exact `t = 4.5`;
authenticated `R02/s10_rankio_t4p5_t5` job `4747146` is accepted through
exact `t = 5.0`;
authenticated `R02/s11_rankio_t5_t5p5` job `4747202` is accepted through
exact `t = 5.5`;
archive and catalog the current committed controller state, reconcile, then
continue only with `R02/s12_rankio_t5p5_t5p75` on one node from its authenticated `s11`
restart siblings under the full F-094 profile: Slurm walltime `01:05:00`, Athena
timeout `00:55:00`, and a one-segment `2700`-second threshold retaining `600`
seconds on both timeout margins. The threshold is scoped only to `s12` and
cannot ratchet automatically. Any scientific, provenance, scheduler, storage,
budget, or reconciliation failure blocks successor preparation. Inspect, account,
reconcile, and recost before any further extension.

The debug-QOS validation launcher remains prohibited for production work.
The archived `E01-pre-modal-driver` Stage I protocol supplied the non-debug
authorization path: it froze the sixteen mapped cases, retained the
conservative `1068.888889` node-hour planning basis, and used sequential
`batch`/default-`normal` submissions with recorded inspection/accounting. Its pre-replacement
execution is accepted only through the `R16` prefix at `t = 7.5`; job
`4686032` reached `t = 8.5` cleanly but was recorded `rejected` after the
replacement turbulent driver changed forcing/restart state.  Stage I has
consumed `9.962778` node-hours with no active segment reservation.  The
historical mapped lower-envelope storage estimates are `811.179` GB under the
sequential retained-policy margin or `1091.586` GB without pruning. Fresh
`E02-modal-driver` R16 reaches exact `t = 10.0` and passes its
production-window analyzer. F-069 historically authorizes only a `2.0` node-hour R02
standard-layout timing pilot next. F-070 records its accepted completion and
authorizes only an `8.0` node-hour R17 high-resolution timing probe. F-071
records its accepted completion and historically authorizes only the frozen
mapped matrix inside the measured `900.000000` node-hour envelope. F-076
retains these records only as pipeline and cost evidence.
Any additional panel-reference-driven case beyond the frozen manifest
requires a revised cost before execution.

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

This tier is scientifically required for a completed reproduction claim.
OLCF prohibits production work in the `debug` QOS. The archived
`E01-pre-modal-driver` Stage I protocol in
`docs/cgl_lf_weak_guide_manuscript_plan.md` permitted only the frozen sixteen
mapped cases, submitted sequentially through
`scripts/frontier/cgl_lf_stage_i.py` on `batch` with default `normal` QOS.
Its accepted `R16` prefix ends at `t = 7.5`; no complete accepted case or
paper-panel result exists in that archived epoch. The replacement driver
invalidates continuation of that lineage. Fresh `E02-modal-driver` R16 jobs
`4743735` through `4744158` are accepted through exact `t = 10.0` after helper
isolation, requalification, and immutable-build archival. R02 jobs `4744198`
and `4744205` close the standard-layout timing gate through native cadence
`t = 0.25`. R17 jobs `4744210` and `4744230` close high-resolution timing
through authenticated continuation `t = 0.10`. F-071 authorizes only the
frozen mapped matrix under sequential inspection. R02 job `4744249` then
reaches exact `t = 1.0`; F-072 requires conservative two-hour-or-shorter
normal-QOS continuation segments. R02 job `4744518` then reaches exact
`t = 1.5`; F-073 authorizes only a conservative continuation through
`t = 2.0`. R02 job `4744913` then reaches exact `t = 2.0`; F-074 authorizes
only a conservative continuation through `t = 2.5`. R02 job `4745305` then
reaches exact `t = 2.5`; F-075 historically authorizes continued
frozen-matrix execution one inspected segment at a time, with R17 last. F-076
supersedes that authorization after the source-to-code forcing audit.
Preserve E02 as pipeline and cost evidence. The retained E03 approval token
binds fresh `R02/s00_rankio_t0_t0p1` job `4745922`, accepted from `t = 0`
through exact `t = 0.1`, and authenticated `R02/s01_rankio_t0p1_t0p25` job
`4746154`, accepted through exact `t = 0.25`. Authenticated
`R02/s02_rankio_t0p25_t1` job `4746182` is accepted through exact `t = 1.0`.
Authenticated `R02/s03_rankio_t1_t1p5` job `4746356` is accepted through
exact `t = 1.5`. Authenticated `R02/s04_rankio_t1p5_t2` job `4746435` is
accepted through exact `t = 2.0`. Authenticated `R02/s05_rankio_t2_t2p5`
job `4746663` is accepted through exact `t = 2.5`. Authenticated
`R02/s06_rankio_t2p5_t3` job `4746773` is accepted through exact `t = 3.0`.
Authenticated `R02/s07_rankio_t3_t3p5` job `4746953` is accepted through
exact `t = 3.5`.
Authenticated `R02/s08_rankio_t3p5_t4` job `4747015` is accepted through
exact `t = 4.0`.
Authenticated `R02/s09_rankio_t4_t4p5` job `4747087` is accepted through
exact `t = 4.5`.
Authenticated `R02/s10_rankio_t4p5_t5` job `4747146` is accepted through
exact `t = 5.0`.
Authenticated `R02/s11_rankio_t5_t5p5` job `4747202` is accepted through
exact `t = 5.5`.
Archive and catalog the
current committed controller state, reconcile, then prepare only
`R02/s12_rankio_t5p5_t5p75` on one node from the authenticated `s11` terminal siblings under
the full F-094 profile: Slurm walltime `01:05:00`, Athena timeout `00:55:00`,
and a one-segment `2700`-second threshold retaining `600` seconds on both
timeout margins. The threshold is scoped only to `s12` and cannot ratchet
automatically. Any scientific, provenance, scheduler, storage, budget, or
reconciliation failure blocks successor preparation. Inspect, account, reconcile, and recost
before any further extension. Complete R02, then
R03--R16 sequentially, and execute R17 last.

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
| Standard mapped comparisons | `R02`-`R09`: active/passive Alfvenic/random roles at mapped beta-10 and beta-100 values |
| Figure 3 compressive additions | `R10`: active random beta `1`; `R11`: active random beta `100` with sonic-correlation forcing |
| Heat-flux sensitivity | `R12`-`R13`, reusing nominal active/passive beta-10 roles from `R02` and `R06` |
| Limiter-frequency scan | `R14`-`R15`, reusing the hard-wall and passive beta-100 roles from `R03` and `R07` |
| Figure 11 scale separation | `R16`-`R17`, reusing the `n_perp = 192` beta-10 role from `R02` |

The separately defined `standard_active_alfvenic_beta1` input is not part of
this minimum matrix: the frozen Stage I manifest records it as unmapped and
unauthorized because the identified beta-1 displayed role is the random
compressive `R10` case. Executing that separate definition would require a
new source-to-panel mapping and revised cost/authorization record.

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

Like Tier 3, this is not an authorized `debug`-QOS workload.  The frozen
Stage I production manifest includes the mapped Figure 11 scale-separation
cases, which may be run only through the sequential production protocol.
Other convergence additions still require a revised mapped scope and cost.

Selected cases should be repeated at:

```text
96x96x192
192x192x384
384x384x768, resources permitting
```

F-041 defines the active Alfvenic beta-10 Figure 11 realization of this
three-resolution suite behind the same explicit paper-execution guard and
admits its dimensionless lower-panel alignment references. The archived
`E01-pre-modal-driver` Stage I protocol maps those three runs and retained an
accepted `R16` prefix only; it does not qualify the dimensional upper-panel
kinetic spectrum or establish a completed comparison.

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

For the frozen sixteen-case Stage I manifest, the historical G011
shared-file sizes plus the F-041 cell-count-scaled resolution estimates define
an archived lower-envelope sequential-retention allocation of `811.179` GB
including margin. The `E01-pre-modal-driver` protocol retained the more
conservative earlier `880.368` GB allocation. Neither is an
`E02-modal-driver` reservation: the tracked decks now select rank-local
snapshot/restart products, and the replacement driver adds authoritative
modal state to restarts. F-069 records `4686854064` raw snapshot/restart
bytes (`4.365` GiB) across the completed segmented R16 tree and a
`4514126544`-byte deduplicated logical cadence. G023 terminal groups project
`36091999664` bytes (`33.613` GiB) for one standard-layout case at the same
cadence. F-070 R02 terminal groups refine that standard-layout projection to
`36092008728` bytes (`33.613` GiB). F-071 R17 ranked groups replace cell
scaling with a measured high-resolution estimate: the mapped matrix projects
`798559758560` logical raw bytes (`743.717` GiB), and running R17 last yields
a `622953051824`-byte sequential peak before margin. Reserve
`747543662189` bytes after the 20-percent sequential-retention margin. Any
change in case matrix, cadence, file layout, retention policy, or required
panel coverage must revise the storage/runtime proposal.

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
external reference-data manifests now produce checksum/uncertainty-qualified
residuals and overlay figures for supported curves, sampled
`pressure_density_joint.*` surfaces, and threshold-volume histories.
Populated MKS24 panel comparisons, exact
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

- `paper-smoke`, `paper-standard`, `paper-nulim`, `paper-heat-flux`,
  `paper-compressive`, `paper-scale-separation`, and analysis modes;
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

This paragraph originally withheld paper-production submission pending the
debug qualification ladder and resource estimate.  The subsequent tracked
Stage I protocol in `docs/cgl_lf_weak_guide_manuscript_plan.md` supplies that
execution update: it freezes a deduplicated sixteen-case mapped matrix,
records a `4000` node-hour project ceiling, and permits only sequential
Frontier `batch`/default-`normal` submissions through
`scripts/frontier/cgl_lf_stage_i.py`.  It does not permit paper production
through Frontier `debug`.  The archived pre-replacement Stage I execution
accepted `R16` only as a prefix through `t = 7.5`; job `4686032` reached
`t = 8.5` cleanly but is rejected for continuation after the turbulent-driver
forcing/restart-state contract changed. Fresh replacement-driver `E02`
execution has since completed R16 from `t = 0` through exact `t = 10.0` in
jobs `4743735` through `4744158`. R02 jobs `4744198` and `4744205` close the
standard-layout timing gate through native cadence `t = 0.25`. R17 jobs
`4744210` and `4744230` close high-resolution timing through authenticated
continuation `t = 0.10`. F-071 permits only the frozen `R02`--`R17` mapped
matrix inside the `900.000000` node-hour envelope. R02 job `4744249`
subsequently reaches exact `t = 1.0`; F-072 encodes the two-hour normal-QOS
preparation limit. R02 job `4744518` subsequently reaches exact `t = 1.5`;
F-073 authorizes only a conservative continuation through `t = 2.0`. R02
job `4744913` subsequently reaches exact `t = 2.0`; F-074 authorizes only a
conservative continuation through `t = 2.5`. R02 job `4745305` subsequently
reaches exact `t = 2.5`; F-075 historically authorizes continued
frozen-matrix execution one inspected segment at a time, with R17 last. F-076
supersedes that authorization: prohibit every new E02 preparation or
submission.

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
scripts/digitize_cgl_lf_mks24_fig2a.py
scripts/digitize_cgl_lf_mks24_fig2b.py
scripts/digitize_cgl_lf_mks24_fig4a.py
scripts/digitize_cgl_lf_mks24_fig7.py
scripts/digitize_cgl_lf_mks24_fig8.py
scripts/digitize_cgl_lf_mks24_fig9.py
scripts/digitize_cgl_lf_mks24_fig11.py
scripts/digitize_cgl_lf_mks24_fig12.py
scripts/digitize_cgl_lf_mks24_fig13.py
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

Figure 2(b), Figure 4(a) normalized-density PDFs, Figure 5(b) normalized
eddy-anisotropy curves, dimensionless Figure 7
lower-panel transfer, dimensionless Figure 9, Figure 11 lower-panel, and
Figure 12 alignment, and dimensionless Figure 13(b),(d) data are now
extracted from pinned vector PDFs by
`scripts/digitize_cgl_lf_mks24_fig2b.py` and
`scripts/digitize_cgl_lf_mks24_fig4a.py` and
`scripts/digitize_cgl_lf_mks24_fig5b.py` and
`scripts/digitize_cgl_lf_mks24_fig7.py` and
`scripts/digitize_cgl_lf_mks24_fig9.py` and
`scripts/digitize_cgl_lf_mks24_fig11.py` and
`scripts/digitize_cgl_lf_mks24_fig12.py` and
`scripts/digitize_cgl_lf_mks24_fig13.py`, with checksum-qualified external
manifests and explicit line-extraction uncertainty. Under F-043,
`scripts/digitize_cgl_lf_mks24_fig8.py` additionally extracts ten selected-
shell `alignment.<shell>` PDF curves from the pinned RGB heatmaps using
their labeled linear colorbar, a `0.03` unit-integral guard, and absolute
density uncertainty `0.05`. Figure 9, Figure 11 lower-panel, and Figure 12
comparisons use `alignment_peak.cos_theta` and require explicitly selected
alignment shells spanning the reference coordinates; Figure 8 comparisons
require shells `8,16,32,64,128`. Figure 7 lower-panel and Figure 13(d) use
`pressure_transfer.transfer_normalized_by_total` under F-039. Under F-038,
Figure 7's upper spectra, Figure 11's upper kinetic spectrum, Figure 12's
lower spectrum, and Figure 13(a),(c) remain excluded because their
dimensionful ordinates do not have a qualified paper-to-AthenaK transform.
Under F-051, `scripts/digitize_cgl_lf_mks24_fig2a.py` first decoded Figure
2(a)'s labeled tiled raster into donor-pressure-coordinate samples with
absolute density uncertainty `0.5`. Under F-052, the extractor also requires
the pinned `MKS24.tex` source and emits `surfaces.json`: the source-derived
beta-10 pressure scale `s = 0.5` maps both coordinates and transforms density
and its uncertainty by `1/s^2`.
Under F-047/F-050, Figure 2(a) now maps to the emitted
`pressure_density_joint` diagnostic and a sampled-surface comparison route,
and under F-052 its labeled raster/colorbar is admitted as converted sampled
surface comparison data. Under F-048, Figures 3, 4(b), and 6(a)
  have matching emitted spectral-field products, but their raw plotted
  ordinates remain excluded because they do not yet map to qualified supported
  normalizations. F-053 confirms that the manuscripts' shell definitions and
  plotted units do not supply the missing discrete Fourier normalization and
  that the 2026-05-26 follow-up still found no `Data availability` heading in
  retrieved official page text or matching exact-query public
  numeric-data/donor-code record.
  Under F-049, Figure 3's sonic-correlation and beta-1 random
  cases are guarded definitions only; they are not simulation evidence.
Data still required:

- obtain machine-readable original data from an author/archive for remaining
  quantitative spectral/strain panels; or
- obtain donor diagnostic code or an explicit discrete Fourier-normalization
  statement sufficient to qualify a digitized ordinate conversion.

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
11. Before any Frontier submission, confirm the launcher's allowed QOS,
    sequential single-job rule, run-root placement, and remaining node-hour
    budget. Debug-only qualification uses `debug`; canonical E03 production
    uses `batch` with Frontier's default `normal` QOS.
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

1. **Preserve the reviewed E03 approval token and clean reconciliation.**
   F-078 closes the corrected Frontier ladder for revision
   `9e07542281e4e6d125582f253df3ad2e3b8b154d`, executable SHA-256
   `68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c`,
   and source-bundle SHA-256
   `c39d55809989d20aa5438711803f4fd43237fa9284c7e15183e84fdd693d3687`.
   Jobs `g024`--`g031` are retained and independently reviewed. Token SHA-256
   `e3fec9f35da42121b902f41ef752021f375aab38bc34c5ff75ae8780bfbca635`
   was retained atomically and the initially empty E03 ledger reconciled
   cleanly before the first reservation.

2. **Preserve the accepted fresh E03 R02 pilot from `t = 0`.**
   Use
   `/autofs/nccs-svm1_home2/dfielding/athenak-df/scripts/frontier/cgl_lf_stage_i.py`
   beneath `E03-forcing-policy`.
   Job `4745922` submitted `R02/s00_rankio_t0_t0p1` with absolute
   `time/tlim = 0.1`, one node, eight GPU ranks, Slurm walltime `01:30:00`,
   and Athena timeout `01:15:00`. Preparation and atomic submission required
   an empty user queue and explicit acknowledgement of the reviewed stale
   `beta25-accel05-gamma10001-purecgl-256` campaign record without modifying
   it. Formal inspection accepted exact `t = 0.1`, zero strict LF safety
   counters, forcing-work relative residual `1.920820308941694e-11`, and
   complete initial/final eight-sibling ranked snapshot and restart groups.
   Accounting records `0.188889` accepted node-hours. The record-time
   `segment_inspection.json` and embedded `scientific_inspection` payload
   SHA-256 are
   `f4b0a31304e31610ec9ec83133f79b747fb0d8cf0f9b3b013febc1cddc2ae784`.
   F-079 corrects the nonblocking preview-rendering omission found during this
   submission. Job `4745922` retained the reviewed acknowledgement and is
   valid. F-079/F-080 controller provenance is archived. Authenticated
   continuation `R02/s01_rankio_t0p1_t0p25` job `4746154` passes formal and
   independent inspection through exact `t = 0.25`, records `0.289167`
   node-hours, and leaves corrected E03 cumulative use at `0.478056`.
   Retained F-081 recost evidence JSON SHA-256 is
   `eb6071cf53d453b2f75707003616ac9cbdfe025d62969d00e34e7948bee3310c`.
   Post-record lifecycle hardening is retained by F-082. Authenticated
   `R02/s02_rankio_t0p25_t1` job `4746182` passes formal and independent
   inspection through exact `t = 1.0`, records `1.465556` node-hours, and
   leaves corrected E03 cumulative use at `1.943612`. Retained F-085 recost
   evidence JSON SHA-256 is preserved as superseded arithmetic history:
   `7d11e8a0004e24417167aeb1f186f8b1f68c8eb9d6206467a708a7cad16cf252`.
   Authenticated `R02/s03_rankio_t1_t1p5` job `4746356` passes formal and
   independent inspection through exact `t = 1.5`, records `1.048611`
   node-hours, and leaves corrected E03 cumulative use at `2.992223`.
   Retained corrected F-086 recost evidence JSON SHA-256 is
   `b6fd6e8fcd939ca1baa06411a3cd6f04359b3bc2ff5b2a0097e2cc87b1763fc7`.
   Authenticated `R02/s04_rankio_t1p5_t2` job `4746435` passes formal and
   independent inspection through exact `t = 2.0`, records `1.063611`
   node-hours, and leaves corrected E03 cumulative use at `4.055834`.
   Retained corrected F-087 recost evidence JSON SHA-256 is
   `c382b78ab9e21466648adf9d7ea38d2b407f0e986579f128bdfe6bc44bef79dd`.
   Authenticated `R02/s05_rankio_t2_t2p5` job `4746663` passes formal and
   independent inspection through exact `t = 2.5`, records `1.062500`
   node-hours, and leaves corrected E03 cumulative use at `5.118334`.
   Retained corrected F-088 recost evidence JSON SHA-256 is
   `39115585cf7ca866e3aade1e53c69aab6ee9217361d75db7946748905e1c0e5c`.
   Authenticated `R02/s06_rankio_t2p5_t3` job `4746773` passes formal and
   independent inspection through exact `t = 3.0`, records `1.077500`
   node-hours, and leaves corrected E03 cumulative use at `6.195834`.
   Retained corrected F-089 recost evidence JSON SHA-256 is
   `d2c5e27553426beb03c8c1882bd6473b682029c6ba67acd32bb06c5fb9a3ff01`.
   Authenticated `R02/s07_rankio_t3_t3p5` job `4746953` passes formal and
   independent inspection through exact `t = 3.5`, records `1.151111`
   node-hours, and leaves corrected E03 cumulative use at `7.346945`.
   Retained corrected F-090 recost evidence JSON SHA-256 is
   `468e787db48d98661cd3ffe125829a78f135bb7ca4c9104c766a00b693f586a3`.
   Authenticated `R02/s08_rankio_t3p5_t4` job `4747015` passes formal and
   independent inspection through exact `t = 4.0`, records `1.200000`
   node-hours, and leaves corrected E03 cumulative use at `8.546945`.
   Retained corrected F-091 recost evidence JSON SHA-256 is
   `dd5b0f5b82cf81005c1a481ee83d1127aa43861dcdd8c22765170e28b0ad6c99`.
   Authenticated `R02/s09_rankio_t4_t4p5` job `4747087` passes formal and
   independent inspection through exact `t = 4.5`, records `1.208056`
   node-hours, and leaves corrected E03 cumulative use at `9.755001`.
   Retained corrected F-092 recost evidence JSON SHA-256 is
   `8b10b1786db520740b687d989207f3cda6e0123f9bf2b62e75a5df3849511561`.
   Authenticated `R02/s10_rankio_t4p5_t5` job `4747146` passes formal and
   independent inspection through exact `t = 5.0`, records `1.226389`
   node-hours, and leaves corrected E03 cumulative use at `10.981390`.
   Retained corrected F-093 recost evidence JSON SHA-256 is
   `5256fb7fd9b75f175a16c2c16ee8a995b8d6d9281a92b70a979ea62171835bf2`.
   Authenticated `R02/s11_rankio_t5_t5p5` job `4747202` passes formal and
   independent inspection through exact `t = 5.5`, records `1.268889`
   node-hours, and leaves corrected E03 cumulative use at `12.250279`.
   Retained corrected F-094 recost evidence JSON SHA-256 is
   `323d7cb3cbe108eee944045c189cff75ffc07d06833c49bcfa8315f07fad8d51`.
   F-083/F-084 controller hardening is promoted at canonical revision
   `9689c269bf329542815a1b2b137881126964b05c`, helper SHA-256
   `1c633ebb58294938a0a0609742ae8f8d1f88cd15242ffe649578796edeb39375`;
   focused tests, Sphinx
   warnings-as-errors, syntax, diff, and hardened reconciliation pass.
   Archive and catalog the current committed controller state, reconcile,
   then prepare only `R02/s12_rankio_t5p5_t5p75` on one node from the authenticated `s11`
   terminal siblings with Slurm walltime `01:05:00`, Athena timeout `00:55:00`,
   and its reviewed one-segment `2700`-second threshold retaining `600`
   seconds on both timeout margins. The threshold is scoped only to `s12` and
   cannot ratchet automatically. Any scientific, provenance, scheduler,
   storage, budget, or reconciliation failure blocks successor preparation.
   Inspect, account, reconcile,
   and recost. Finish R02, execute R03--R16 sequentially, and execute R17
   last.

3. **Recost corrected E03 and continue the frozen matrix sequentially.**
   Preserve E01 and E02 as historical pipeline and cost evidence only. No E02
   restart is admissible. Use the measured E03 R02 pilot and authenticated
   continuations to revise segment increments conservatively. Run only one
   production segment at a time, retain accepted snapshots plus final two
   restart groups, inspect and account every segment, and complete `R17`
   last. The historical `900.000000` node-hour mapped envelope and E02 timing
   rates are planning priors, not permission to skip corrected-E03 recosting.

4. **Complete reference coverage and long-time calibration in parallel.**
   Figure 2(b)'s eight unstable-volume histories, Figure 7's dimensionless
   lower-panel transfer curves, Figure 4(a)'s normalized-density PDFs,
   Figure 5(b)'s normalized eddy-anisotropy curves,
   Figure 8's selected-shell alignment PDFs, Figure 9's, Figure 11's
   lower-panel, and Figure 12's dimensionless alignment curves, and Figure
   13(b),(d)'s dimensionless curves are populated through pinned reference
   extractors.
   F-047/F-050 now supply Figure 2(a)'s joint diagnostic and sampled-surface
   comparison interface in AthenaK coordinates; F-051 decodes its raw
   color-calibrated raster into donor-coordinate samples, and F-052 applies
   the pinned-source beta-10 pressure conversion and PDF Jacobian to admit
   those surfaces for comparison. F-048
   supplies the missing Figure 3/4(b)/6(a) spectral fields,
   and F-049 supplies guarded Figure 3 case definitions, but neither admits
   raw-ordinate spectral reference transforms nor execution evidence. Resolve
   the F-038/F-053 paper-to-AthenaK normalization boundary before admitting
   dimensional spectra or strain: obtain author/archive numeric data, donor
   diagnostic code, or an explicit discrete Fourier-normalization statement
   before digitizing absolute spectral/strain ordinates. Revise the proposal
   when scope changes, and
   calibrate long-time injection/spectral behavior against the documented
   setup.

5. **Complete remaining comparison products from accepted corrected-epoch
   outputs.**
   Run paper-panel reference comparisons through the implemented manifest
   interface on top of the cap, forcing-work residual, applied
   heat-flux/pressure-traction ledgers, snapshot heat-flux/pressure-work
   products, and steady-window PDF/joint-PDF/spectral/transfer/alignment/strain
   products and generic figure renderer.
