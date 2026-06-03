# AthenaK MHD-PIC Sun & Bai Reproduction and Release Work Plan

## Document Status

| Item | Value |
| --- | --- |
| Date | 2026-05-31 |
| Scope | Dependency-ordered execution companion for the AthenaK MHD-PIC release effort |
| Governing document | [`PIC_PRODUCTION_READINESS_PLAN.md`](PIC_PRODUCTION_READINESS_PLAN.md) |
| Bulk-artifact root | `/lustre/orion/ast207/proj-shared/dfielding/PIC` |
| Frontier testing cap | `10000` cumulative reserved-plus-consumed node-hours |
| Current release status | Not production-ready |
| External review | Pending until a named reviewer is assigned |

This document is a concise operational companion to the canonical
[`PIC_PRODUCTION_READINESS_PLAN.md`](PIC_PRODUCTION_READINESS_PLAN.md). It does
not replace that document, relax any gate, or independently authorize Frontier
execution. If the two documents differ, the canonical readiness plan controls.

## Objective

Deliver a rock-solid AthenaK MHD-PIC implementation that:

1. Reliably reproduces the non-relativistic parallel-shock acceleration
   benchmark in Section 5.4 of Sun & Bai (2023), including the qualitative
   morphology and downstream particle spectra.
2. Quantitatively reproduces the complete Sun & Bai test suite:
   Sections 5.1 through 5.7 and the relevant expanding-box appendices.
3. Is reliable under the execution modes needed for later science:
   conservation, AMR, SMR, boundaries, particle migration, load balancing,
   checkpoint/restart, MPI decomposition, Frontier HIP/GPU execution,
   interruption recovery, and output reconstruction.
4. Produces archived raw artifacts, metrics, and figures that can be
   independently regenerated from a frozen clean candidate.
5. Separately implements and qualifies the selected optional extensions:
   Hall Bell, ion-neutral-damped CRSI, and adaptive-delta-f physical-damping
   CRPAI.

## Explicit Non-Goal

This work plan does **not** include development or execution of a driven
turbulent MHD-PIC box. That is a later science application. The purpose of this
plan is to make the underlying MHD-PIC implementation trustworthy before that
application is attempted.

## Release Milestones

The project should advance through four explicit milestones. Completion of an
earlier milestone must not be reported as completion of a later one.

| Milestone | Meaning | Minimum evidence |
| --- | --- | --- |
| M1: Shock-qualified candidate | Section 5.4 is reproduced by the frozen candidate | Calibrated coarse, fine, and AMR shock runs; restart parity; MPI/HIP evidence; regenerated figures and metrics |
| M2: Sun-Bai reproduction bundle | The complete paper test suite is reproduced | Sections 5.1-5.7 and Appendices A/B pass with archived raw artifacts and deterministic analysis |
| M3: Production-ready `paper_mhd_pic_vl2_tsc` | The additive VL2/TSC paper mode is suitable as a supported release profile | M2 plus nonlinear, reliability, portability, scaling, documentation, and archive-integrity closure |
| M4: Production-ready `extended_mhd_pic` | The previously selected optional extensions are supported | M3 plus qualified Hall Bell, ion-neutral-damped CRSI, and adaptive-delta-f physical-damping CRPAI |

## Current Position

The codebase contains substantial local-readiness work, but none of the four
release milestones is complete.

### Locally Demonstrated Building Blocks

| Area | Current bounded evidence | Important limitation |
| --- | --- | --- |
| Paper-mode identity | Active additive `paper_mhd_pic_vl2_tsc` selects the VL2/TSC paper candidate, rejects direct-current CT induction, and requires coherent coupled feedback; `paper_mhd_pic` is retained as archival chronology only | Broader campaign matrix remains open |
| Relativistic pusher | Host convergence and a registered one-rank Frontier GPU gyro oracle pass | Full CPU/MPI/GPU portability matrix remains open |
| Conservative coupling | Host and registered one-rank Frontier GPU coupling oracles pass | Multi-rank and paper-campaign matrix remains open |
| True delta-f | State, deposition, restart, and bounded source-local mechanics exist | Long-horizon CRSI/CRPAI paper reproduction remains open |
| Expanding box | Bounded host CPAW, gyro-history, comoving-flux, div(B), and restart checks exist | Full paper Section 5.7 and Frontier qualification remain open |
| AMR lifetime | Serial refine/derefine, restart, boundary, ASan, and UBSan successors pass | Multi-rank MPI migration, HIP lifetime, and scientific AMR equivalence remain open |
| Shock generator | Paper-surface injection and startup-cohort removal are implemented; a bounded actual-particle audit passes; the first registered HIP pressure slice exposed a strided host-to-device copy fault; the packed-transfer repair passes bounded regressions; the rebuilt HIP/MPI v2 retry completed four serialized pressure-engineering calibration slices | Aggregate calibration publication, human pressure review, prerequisite slices, and qualifying shock campaigns remain open |
| Particle provenance | Schema-7 provenance and independently reconstructed weighted spectra pass bounded local tests, including a two-rank Orion-local replay | Frontier HIP parity and shock-campaign binding remain open |
| Frontier control plane | Serialized ledger, root restriction, registered F0/F1 slices, and bounded F2 engineering evidence exist | Remaining registered science and stress slices remain open |

### Section 5.4 Preparation Already Present

The active prepared paper deck is
[`inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput`](../../inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput).
The sibling deck without the `_vl2_tsc` suffix is retained as archival
pre-successor chronology only.
The source-local preparation record is
[`q011_parallel_shock_section54_paper_preparation_successor_v10_2026-06-02.json`](readiness/q011_parallel_shock_section54_paper_preparation_successor_v10_2026-06-02.json).
The bounded injection audit is
[`q011_injection_distribution_runtime_local_successor_v8_2026-06-02.json`](readiness/q011_injection_distribution_runtime_local_successor_v8_2026-06-02.json).

The first registered Frontier pressure slice failed closed before producing
physics evidence: the HIP executable attempted a strided host-to-device
particle subview copy during injection. The packed particle-transfer repair
passes bounded serial-host append, recenter and startup-cohort compaction
regressions. The rebuilt clean HIP/MPI candidate then completed the registered
v2 pressure-pilot retry for `problem/ps_p0 = 1.0`, `0.05`, `0.10`, and `0.20`.
All four slices are reconciled with immutable raw descriptors. Aggregate
publication is paused after a failed-closed parser-compatibility discovery and
two intentionally cancelled worker reruns exposed additional publication
hardening requirements before any public bundle was accepted. Human pressure
review, qualifying-plan hardening, and prerequisite slices remain open.

The prepared contract includes:

- a reflecting left wall and periodic transverse direction;
- upstream flow and magnetic field parallel to `x`;
- `M_A = 30`, `gamma = 5/3`, `eta = 10^-3`, and `C/U_A0 = 10000`;
- continuous injection at the ideal shock surface;
- isotropic monoenergetic injection with `p/m = sqrt(10) u0`;
- gas mass, momentum, and energy subtraction;
- removal of particles born before `45 Omega_0^-1`;
- coarse, fine, and three-level AMR variants with cell sizes `12`, `6`, and
  `3 c/omega_pi`;
- mesh `bin` outputs for density, magnetic-field magnitude, and current;
- particle `pvtk` outputs for full particle payloads;
- restart checkpoints.

This is preparation, not reproduction evidence.

## Operating Rules

These rules apply to every phase.

1. Use `/lustre/orion/ast207/proj-shared/dfielding/PIC` as the exact root for
   every Frontier simulation artifact. Do not use Kronos.
2. Run no large or qualifying campaign from a dirty worktree.
3. Bind every qualifying run to:
   - clean git commit;
   - clean source-bundle checksum;
   - executable checksum;
   - input-deck checksum;
   - analyzer checksum;
   - scheduler script;
   - redacted environment allowlist;
   - artifact root;
   - node-hour ledger mutation.
4. Submit at most one AthenaK PIC job at a time across both QOS classes.
5. Use the registered submission wrapper and locked ledger. Do not submit
   qualifying runs with direct `sbatch`.
6. Use `debug` only for eligible short non-production work. Use `normal` on
   `batch` for full shock, nonlinear, saturation, and controlled scaling runs.
7. Refuse submission when consumed node-hours, active reservations, and the
   proposed maximum reservation would exceed `10000` node-hours.
8. Freeze estimators, fit windows, exclusions, seeds, and tolerances before
   inspecting qualifying outputs.
9. Archive failures and rejected attempts. Do not silently discard
   discrepancies or broaden thresholds after inspecting results.
10. Keep paper reproduction and optional-extension results separate in
    manifests, figures, and claims.
11. Run the opt-in Orion live-state preflight with
    `PIC_RUN_LIVE_PREFLIGHT=1` immediately before paired control-plane install,
    policy promotion, reservation, and launch. Require exact mirrored policy
    and ledger state, an accepted empty-queue snapshot, no pending marker, and
    no active reservation before proceeding.

## Phase Overview

| Phase | Objective | Depends on | Exit result |
| --- | --- | --- | --- |
| 0 | Curate and freeze the candidate | Existing local-readiness work | Clean reviewable release candidate |
| 1 | Close local physics and shock-calibration gaps | Phase 0 | Frozen Section 5.4 physical contract and analyzers |
| 2 | Close registered Frontier prerequisite slices | Phases 0-1 | MPI/HIP, AMR, restart, telemetry, and resilience evidence |
| 3 | Reproduce analytical and bounded Sun-Bai tests | Phase 2 | Sections 5.1-5.3, 5.5-5.7, and Appendices A/B pass |
| 4 | Execute the full Section 5.4 shock campaign | Phases 1-3 | M1: shock-qualified candidate |
| 5 | Freeze the complete paper bundle | Phase 4 | M2: Sun-Bai reproduction bundle |
| 6 | Complete paper-mode production hardening | Phase 5 | M3: production-ready `paper_mhd_pic_vl2_tsc` candidate |
| 7 | Qualify selected optional extensions | Phase 6 | M4: production-ready `extended_mhd_pic` candidate |
| 8 | Perform terminal review and sign-off | Phases 6-7 as selected | Signed release disposition |

## Phase 0: Curate and Freeze the Candidate

### Objective

Turn the current dirty worktree into a reviewable, immutable candidate before
using additional Frontier allocation.

### Required Work

1. Inventory all modified and untracked files.
2. Separate source changes, tests, documentation, readiness sidecars, generated
   report files, and non-versioned artifacts.
3. Split retained source work into narrow commits with matching tests.
4. Review the final diff against nearby AthenaK patterns.
5. Build the final host candidate and run the broad local suite.
6. Create the canonical clean source bundle and executable pin.
7. Bind the canonical paper release profile explicitly to double precision
   with `Athena_SINGLE_PRECISION=OFF`. Single-precision shock runs require a
   separate qualification campaign.
8. Regenerate and review
   `tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json`,
   then freeze its exact checksums for all prepared PIC decks and publication
   analyzers from archived source bytes.
9. Preserve predecessor evidence as chronology without promoting it to release
   evidence.

### Required Verification

- `git diff --check`
- touched-file C++ lint
- Python AST or syntax checks
- publication test suite
- readiness-registry suite
- qualification-manifest suite
- warning-as-error documentation build
- host Debug and Release builds
- MPI and OpenMP compile coverage
- Orion live-state preflight with `PIC_RUN_LIVE_PREFLIGHT=1` immediately before
  paired install and policy promotion

### Exit Gate

- The candidate worktree is clean.
- The reviewed source commit, recursive source bundle, executable, and
  analysis-tool checksums are archived.
- No generated simulation output is accidentally versioned.
- Every later run can point to a single immutable candidate identity.

### Execution Status: In Progress

The first curation pass exposed two freeze-blocking provenance gaps during
independent review: clean-candidate creation accepted an operator-selected
prepared-artifact subset, and Q-006, Q-007, and Q-011 analyzers reopened
retained evidence by pathname after verification. A second independent review
then found that the copied private snapshot remained owner-mutable under the
same Unix account. A third stable-boundary review then found that the handoff
captured its topology map from an unverified scan. A fourth adversarial review
found that unknown post-handoff members could fall back into the writable
topology-only staging tree and that open ordinary inodes could still be
rewritten in place between hashing and semantic inspection. The successor now
requires
the canonical prepared-artifact inventory path, independently derives the
complete archived PIC-deck and publication-analyzer closure, keeps verified
tree topology in memory, and exposes consumed regular payloads lazily through
sealed descriptor-backed snapshot members. It compares the handed-off
topology map back to the verified tree, rejects unknown active-snapshot
members, copies mutation-checked ordinary validator inputs into sealed
descriptors before multi-pass semantic inspection, and rejects sealed
descriptors that were not routed by the active retained-tree snapshot. A fifth
schema-alias audit then found Python numeric-equality gaps at release
boundaries: boolean, floating-point, or string schema aliases could survive
selected profile, receipt, provenance, invocation, recovery, reservation, and
promotion checks. The committed repair requires exact integer schema versions
at those boundaries and adds adversarial regressions. The transition policy
also clears
historical registered-science slices and marks the schema-v4 clean-candidate
freeze pending. The historical nonqualifying transition record is
[`phase0_curated_candidate_successor_2026-05-31.json`](readiness/phase0_curated_candidate_successor_2026-05-31.json).
The current local successor chronology is
[`phase0_curated_candidate_successor_v6_2026-06-01.json`](readiness/phase0_curated_candidate_successor_v6_2026-06-01.json);
it remains nonqualifying until the final clean freeze.

The first exact clean-tree Phase 0 code-boundary validation matrix passed at
`d0ad350562799e609b4d22f45ce9ca314b38da82`: 727 publication tests with two
intentional skips, 364 Frontier control-plane tests, 102 focused hardening
tests, warning-as-error documentation rendering, syntax sweeps and
`git diff --check`. The nonqualifying receipt is
[`phase0_exact_boundary_validation_receipt_2026-05-31.json`](readiness/phase0_exact_boundary_validation_receipt_2026-05-31.json).
Independent rereview rejected that boundary after finding additional
schema-alias paths in policy promotion, qualification registries, pre-submit
coercions, ledger records and Bell retained geometry. The expanded repair is
now pending a fresh exact clean-tree matrix and independent rereview. Phase 0
remains open until the stable final diff is independently rereviewed,
installed as a paired Orion
and Project Home control-plane generation, and frozen with its final source,
executable, deck, and analyzer checksums. No later phase is authorized yet.
The historically installed `6f3458ca5c412b866f579d67dde51d73497b65629f6e7a35e10a3af84d7ecd8b`
control-plane generation is chronology-only: do not launch any new Frontier
job through it. Launches remain prohibited until the reviewed successor is
installed, paired, promoted, and bound to a registered slice.

## Phase 1: Close Local Physics and Shock-Calibration Gaps

### Objective

Resolve the remaining locally decidable Section 5.4 questions before spending
Frontier node-hours on shock pilots.

### Required Shock Work

1. Audit the gas pressure and unit normalization used by
   [`pic_parallel_shock.cpp`](../../src/pgen/tests/pic_parallel_shock.cpp).
2. Verify the artificial-light-speed conversion and all uses of `C/U_A0`.
3. Verify the shock-surface speed:

   $$
   u_{\rm sh}' = \frac{\Gamma - 1}{2} u_0.
   $$

4. Verify the upstream-relative swept-mass budget and `eta = 10^-3`
   injection rate.
5. Verify gas mass, momentum, and energy subtraction for injected particles.
6. Calibrate the macro-particle mass needed to reproduce approximately:
   - `40` downstream particles per finest cell;
   - `640` downstream particles per coarsest cell.
7. Verify the one-time removal of particles injected before
   `45 Omega_0^-1` across uninterrupted and restarted runs.
8. Review output cadence and storage volume for the full runs.
9. Freeze the spectrum estimator and downstream selection:
   - spatial downstream region;
   - source and birth-time filters;
   - energy or momentum binning;
   - fit interval for the late-time slope;
   - snapshot selection tolerance around `t = 500` and `1200`.
10. Freeze AMR-versus-fine residual tolerances before looking at campaign
    results.

### Required Component Work

1. Close the remaining AMR deposition-policy review:
   - use the paper-faithful smooth-interface policy for reproduction;
   - keep any conservative alternative explicitly labeled as an AthenaK
     extension.
2. Review stage ordering for particle push, moment deposition, MHD feedback,
   boundary synchronization, migration, and CT updates.
3. Confirm that the independent offline weighted-spectrum reconstruction uses
   schema-7 particle provenance consistently.
4. Create the qualifying Section 5.4 analyzer. The analyzer must regenerate:
   - shock position;
   - upstream field amplification;
   - density, magnetic, and current morphology;
   - downstream weighted spectra;
   - late-time power-law slope;
   - AMR-versus-uniform residuals;
   - load-balance and performance diagnostics.

### Exit Gate

- Every locally decidable Section 5.4 mapping question has a reviewed answer.
- The deck and analyzer are checksum-frozen.
- Fit ranges, exclusions, tolerances, and seed handling are preregistered.
- No known physical-calibration question is deferred into an expensive run.

## Phase 2: Close Registered Frontier Prerequisite Slices

### Objective

Demonstrate that the candidate behaves correctly under the execution modes
required by the paper suite and the full shock campaign.

### Required Registered Slices

| Slice | Purpose | Required evidence |
| --- | --- | --- |
| AMR MPI/HIP lifetime | Exercise repeated refine, derefine, migration, and ownership refresh | Stable particle identity, no memory faults, restart parity |
| Coupled boundary MPI/HIP | Exercise periodic, reflecting, outflow, and required inflow paths | Conservation, sign, symmetry, and decomposition parity |
| Q016 HIP parity | Verify provenance and spectra on Frontier | Schema-7 payload parity and independent spectrum reconstruction |
| Restart publication | Exercise per-rank manifests and completion markers | Valid last-known-good recovery and deterministic rejection of corrupt artifacts |
| Scheduler-pretimeout continuation | Exercise controlled interruption | Restarted endpoint agrees with uninterrupted control |
| Frontier filesystem failure drills | Exercise partial or failed writes | Fail-closed behavior with preserved recovery point |
| Q017 telemetry | Measure timers, load, and memory | Retained stage timers, rank-load metrics, tracked GPU memory, and output overhead |
| GPU-aware MPI A/B | Freeze communication choice | Correctness parity and measured performance comparison |

### Required Shock-Specific Pilots

1. One short uniform coarse setup run.
2. One short uniform fine setup run.
3. One short AMR run with refinement and derefinement.
4. One short AMR uninterrupted-versus-restart comparison.
5. One short load-balance tuning matrix for
   `pic_load_balance_cost_per_particle`.

These are pilots, not paper reproduction evidence.

### Exit Gate

- Required MPI, HIP, AMR, restart, and telemetry slices pass on the clean
  candidate.
- The selected load-balance cost is frozen with measured evidence.
- Total GPU memory and output overhead are measured well enough to size the
  full paper runs.
- No unresolved runtime fault remains in a code path needed by Section 5.4.

## Phase 3: Reproduce Analytical and Bounded Sun-Bai Tests

### Objective

Run the cheaper and more diagnostic paper tests before committing substantial
resources to the long Section 5.4 campaign.

### Required Paper Matrix

| Paper target | Required result | Main release gates |
| --- | --- | --- |
| Section 5.1 gyro-motion | Relativistic Boris convergence and artificial-`C` dependence on CPU and Frontier GPU | Q-003, Q-025, Q-026 |
| Section 5.2 Bell instability | Real and imaginary dispersion relation in 1D, 2D, and 3D, including signed phase | Q-003, Q-004, Q-005, Q-023, Q-025, Q-026 |
| Section 5.3 oscillation | Analytic frequency agreement on uniform, SMR, and audited AMR grids | Q-004, Q-006, Q-009, Q-025, Q-026 |
| Section 5.5 CRSI | True-delta-f polarization-resolved spectra and fitted linear growth | Q-003, Q-004, Q-007, Q-023, Q-025, Q-026 |
| Section 5.6 CRPAI | Prolate and oblate growth rates and polarization handedness | Q-003, Q-004, Q-007, Q-023, Q-025, Q-026 |
| Section 5.7 driven CRPAI | Expanding/compressing-box distribution and spectral evolution | Q-007, Q-008, Q-023, Q-025, Q-026 |
| Appendix A | CPAW amplitude and phase under expansion and compression | Q-008, Q-025, Q-026 |
| Appendix B | Expanding-box gyro gamma and phase history | Q-003, Q-008, Q-025, Q-026 |

### Required Evidence for Each Test

1. Frozen input deck and analyzer.
2. Clean-candidate executable binding.
3. Raw output and stdout/stderr.
4. Machine-readable measured metrics.
5. Analytical or manuscript reference values.
6. Frozen tolerance and pass/fail result.
7. Restart parity where the mode carries restart state.
8. MPI and HIP parity appropriate to the mode.
9. Independently regenerated plot and metric record.

### Exit Gate

- Every non-shock Sun-Bai benchmark passes or has a reviewed discrepancy that
  is demonstrated to be a documented AthenaK correction.
- The paper-mode pusher, feedback, AMR, delta-f, and expanding-box paths are
  qualified before the full shock allocation is used.

## Phase 4: Execute the Full Section 5.4 Shock Campaign

### Objective

Reproduce the non-relativistic parallel-shock acceleration benchmark with the
frozen, qualified candidate.

### Paper Contract

Use the prepared Section 5.4 contract:

| Quantity | Required value |
| --- | --- |
| Geometry | Reflecting left wall, periodic transverse direction, parallel upstream `B_x` |
| Domain | `(48 x 3.12) * 10^3 c/omega_pi` |
| Alfven Mach number | `M_A = 30` |
| Gas adiabatic index | `gamma = 5/3` |
| Injection efficiency | `eta = 10^-3` |
| Injection momentum | `p/m = sqrt(10) u0` |
| Artificial light speed | `C = 10^4 U_A0` |
| Startup-cohort exclusion | Remove particles injected before `45 Omega_0^-1` |
| AMR cell sizes | `12`, `6`, and `3 c/omega_pi` |
| Refinement thresholds | Refine at `1.0`, derefine below `0.1` |
| Comparison snapshots | Near `t = 500` and `1200 Omega_0^-1` |

### Required Run Matrix

| Run | Purpose | Required comparison |
| --- | --- | --- |
| Short uniform coarse | Early setup sanity | Shock structure, injection budget, conservation |
| Short uniform fine | Resolution response | Compare with coarse trend |
| Short AMR | Refinement and load balance | Compare resolved region with fine run |
| Short AMR restart | Continuation correctness | Compare with uninterrupted AMR control |
| Full uniform coarse | Paper comparison | Morphology, amplification, spectra |
| Full uniform fine | High-resolution paper comparison | Morphology, amplification, spectra |
| Full three-level AMR | Fiducial paper comparison | Fine-run overlap plus AMR efficiency |
| Fixed-seed repeats | Statistical qualification | Apply preregistered eight-seed policy |

### Required Raw Outputs

At each required snapshot, archive:

| Output | Format | Purpose |
| --- | --- | --- |
| Gas density | `bin` | Shock morphology, cavities, filaments |
| Magnetic-field magnitude | `bin` | Upstream amplification |
| Current diagnostic | `bin` | Upstream CR-current and Bell morphology comparison |
| Full particles | `pvtk` | Provenance-filtered phase space and weighted spectra |
| Restart checkpoint | `rst` | Recovery and continuation |
| Stdout telemetry | Text plus manifest | Timers, particle counts, load balance, memory |

### Required Scientific Results

The archived analysis must show:

1. A corrugated shock surface.
2. Upstream cavities and filaments associated with Bell growth.
3. Upstream magnetic-field amplification of approximately `2-4` in the
   relevant regions near `t = 500 Omega_0^-1`.
4. Downstream particle spectra at `t = 500` and `1200 Omega_0^-1`.
5. Increasing maximum particle energy with time.
6. A late-time spectrum trending toward:

   $$
   f(\varepsilon) \propto \varepsilon^{-3/2}.
   $$

7. AMR morphology and spectra that substantially overlap the fine uniform
   result within preregistered tolerances.
8. Measured AMR runtime and memory behavior with an honest statement of the
   benefit or lack of benefit.
9. Load-balancing behavior that addresses the particle concentration near the
   shock and downstream region.

### Exit Gate

- The Section 5.4 figures and metrics regenerate from frozen raw artifacts.
- Independent recomputation agrees with the primary analyzer.
- Coarse, fine, and AMR differences are documented quantitatively.
- Checkpoint/restart, MPI, and HIP behavior is demonstrated for the qualified
  shock path.
- M1 is complete.

## Phase 5: Freeze the Complete Sun-Bai Reproduction Bundle

### Objective

Combine the complete paper suite into one reviewable reproduction artifact
without implying broader production readiness.

### Required Bundle

1. Clean candidate manifest and recursive source-bundle checksum.
2. Executable checksum and environment allowlist.
3. Sections 5.1 through 5.7 metric records.
4. Appendix A and B metric records.
5. Section 5.4 coarse, fine, and AMR raw outputs.
6. Reconstructed paper-comparison figures.
7. Independent metric recomputation results.
8. Statistical qualification records and discrepancy ledger.
9. Node-hour accounting and artifact inventories.
10. Known limitations:
    - no thermal-pool injection model;
    - no implied oblique-shock qualification;
    - no optional-extension substitution for paper mode;
    - no driven turbulent-box application claim.

### Exit Gate

- Every mandatory `sun_bai_2023_reproduction` child gate is closed.
- The immutable aggregate bundle is labeled as paper reproduction, not as
  production-ready software.
- M2 is complete.

## Phase 6: Complete Paper-Mode Production Hardening

### Objective

Convert the paper-faithful candidate into a supported `paper_mhd_pic_vl2_tsc`
release profile suitable for later science applications. Retain
`paper_mhd_pic` as archival pre-VL2 chronology only.

### Required Work

1. Complete nonlinear non-Hall Bell qualification:
   - amplification;
   - wavelength evolution;
   - cavities and filaments;
   - energy transfer;
   - saturation behavior;
   - sensitivity matrix.
2. Complete nonlinear CRSI qualification:
   - growth and saturation;
   - spectra;
   - pitch-angle evolution;
   - scattering and diffusion;
   - weight-validity envelope;
   - reduced full-f controls.
3. Complete driven and undriven CRPAI nonlinear qualification:
   - branch evolution;
   - anisotropy;
   - scattering;
   - saturation;
   - sensitivity matrix;
   - reduced full-f controls.
4. Close independent matched comparisons required by the canonical plan.
5. Complete performance qualification:
   - push, deposition, migration, sorting, and coupling timers;
   - communication time;
   - output overhead;
   - memory high-water marks;
   - AMR overhead;
   - load-balance effectiveness;
   - representative multi-node scaling.
6. Complete resilience:
   - corrupt checkpoint;
   - missing checkpoint member;
   - interrupted writer;
   - scheduler pretimeout;
   - node-loss or equivalent fail-closed drill;
   - Frontier filesystem failure handling;
   - last-known-good recovery.
7. Complete documentation and usability:
   - supported equations and modes;
   - parameter reference;
   - small paper-mode quick start;
   - Frontier runbook;
   - output definitions;
   - troubleshooting;
   - unsupported-composition errors.
8. Complete archive-integrity, restore, licensing, and Orion-retention review.

### Exit Gate

- Accuracy, convergence, conservation, restart, AMR, MPI, GPU, resilience,
  scaling, documentation, and archive-integrity evidence is complete.
- Every enabled paper-mode composition is documented and qualified.
- Unsupported combinations fail clearly before execution.
- M3 is complete.

## Phase 7: Qualify Selected Optional Extensions

### Objective

Complete the selected `extended_mhd_pic` release scope without conflating
extension behavior with Sun-Bai paper mode.

### Required Extensions

| Extension | Required work | Release gate |
| --- | --- | --- |
| Hall Bell | Review physical Bai mapping; freeze coefficient grid and tolerances; run linear, nonlinear, shock-front, MPI, and GPU qualification | Q-029 |
| Ion-neutral-damped CRSI | Freeze matched reduced-map applicability; extract comparison data and thresholds; run damped CRSI, nonlinear, MPI, and GPU qualification | Q-032 |
| Adaptive-delta-f physical-damping CRPAI | Freeze reference mapping and calibration metrics; measure effective scattering, spectra, anisotropy, saturation, restart, MPI, and GPU behavior | Q-033 |

### Exit Gate

- Each extension has an `implemented_and_qualified` disposition.
- Extension figures, manifests, and limitations are separate from paper-mode
  figures and claims.
- M4 is complete.

## Phase 8: Terminal Review and Sign-Off

### Objective

Produce a defensible release decision after all selected child gates close.

### Required Work

1. Freeze the aggregate claim manifest with disposition
   `pending_terminal_review`.
2. Assign a named external reviewer.
3. Review:
   - source candidate;
   - paper reproduction bundle;
   - production-hardening evidence;
   - extension bundles if included;
   - discrepancy ledger;
   - unsupported-mode register;
   - node-hour ledger;
   - Orion-only retention risk;
   - fresh-directory restore drill.
4. Record the terminal Q-014 disposition.
5. Sign and archive the final manifest.

### Exit Gate

- The signed terminal manifest states exactly which milestone is complete.
- No `pending external review` placeholder is promoted without an assigned
  reviewer disposition.
- Remaining limitations are explicit and scoped.

## Section 5.4 Analyzer Checklist

Before the first qualifying shock submission, the analyzer must fail closed
unless all required bindings and products are present.

| Check | Required behavior |
| --- | --- |
| Candidate identity | Require clean git commit and clean executable checksum |
| Root containment | Require every artifact below the authorized Orion root |
| Snapshot timing | Select only snapshots within frozen tolerances of `t=500` and `1200` |
| Particle provenance | Use only `shock_injected` particles with required birth-time and downstream filters |
| Mesh products | Require density, magnetic magnitude, and current `bin` files |
| Particle products | Require complete `prtcl_all` `pvtk` payloads |
| Restart products | Require validated restart manifests and completion markers |
| Spectra | Independently regenerate weighted spectra from raw particle payloads |
| Morphology | Generate density, magnetic, current, and CR-distribution panels |
| AMR comparison | Compute preregistered AMR-versus-fine residual metrics |
| Performance | Parse Q017 telemetry, load distribution, runtime, and memory |
| Failures | Preserve rejected attempts and refuse partial promotion |

## Completion Checklist

| Deliverable | Status at creation of this document |
| --- | --- |
| Clean curated candidate | Open |
| Frozen Section 5.4 physical calibration | Open |
| Frozen Section 5.4 qualifying analyzer | Open |
| Registered Frontier AMR MPI/HIP stress | Open |
| Registered Frontier restart and filesystem resilience matrix | Open |
| Registered Q016 HIP provenance and spectra parity | Open |
| Complete Section 5.1-5.3 analytical paper qualification | Open |
| Complete Section 5.5-5.7 and Appendix A/B qualification | Open |
| Full Section 5.4 coarse/fine/AMR paper campaign | Open |
| M1 shock-qualified candidate | Open |
| M2 Sun-Bai reproduction bundle | Open |
| M3 production-ready `paper_mhd_pic_vl2_tsc` | Open |
| M4 production-ready `extended_mhd_pic` | Open |
| Named external review and terminal Q-014 sign-off | Pending external review |

## Immediate Next Actions

Execute these actions in order:

1. Commit and push the integrated ninth authenticated-storage,
   exact-historical-anchor and candidate-lifecycle repair tree, and run the
   full clean-worker validation matrix from that commit.
2. Capture authenticated mirrored storage-preflight evidence, install the
   paired controller, promote the one-use consumed-slice retirement successor,
   freeze and independently revalidate a fresh worker-built candidate, and
   promote the candidate-only successor.
3. Provision the fixed sibling publication-acceptance authority through the
   reviewed one-time transition. Finish and verify the immutable four-slice
   pressure aggregate publication and pressure-review packet on workers.
4. Record the required human
   `problem/ps_p0` selection receipt.
5. Close the Section 5.4 normalization, macro-particle-mass, gas-subtraction,
   spectrum-window, snapshot-tolerance, and AMR-residual reviews.
6. Materialize and freeze the qualifying Section 5.4 plan and analyzer.
7. Run registered Frontier AMR MPI/HIP, Q016 HIP, restart-resilience, Q017
   telemetry, and load-balance prerequisite slices.
8. Reproduce the cheaper analytical Sun-Bai cases.
9. Run the full registered Section 5.4 coarse, fine, and AMR campaign.
10. Freeze the paper reproduction bundle.
11. Complete paper-mode production hardening.
12. Qualify the selected optional extensions.
13. Perform named external review and terminal sign-off.

Do not launch the full Section 5.4 campaign before Actions 1 through 8 close.
