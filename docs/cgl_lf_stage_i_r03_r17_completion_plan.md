# CGL-LF MKS24 Stage I Rapid Completion Plan: R03 Through R17

## 1. Purpose

This document is the execution-focused plan for completing the corrected
AthenaK CGL-LF MKS24 Stage I matrix from the present `R03` stop line through
accepted `R17` completion and the final campaign analysis bundle.

The goal is to move quickly. The implementation is already scientifically
credible, the corrected forcing policies are qualified, `R02` is complete,
and independent use by Stephen Majeski provides encouraging external evidence.
The remaining work is a production campaign with one known restart-metadata
hardening issue, not an open-ended physics-development program.

This plan preserves fail-closed controller behavior while using parallel
subagents aggressively for software hardening, packet review, scientific
monitoring, recosting, analysis, provenance, and documentation. After the R03
recovery gate, it deliberately promotes a bounded-concurrency controller
transition: shared-root metadata mutations remain serial, but independent
Frontier case lanes may run concurrently.

### 1.1 Acceleration Transition Status

The bounded-concurrency controller transition is implemented in the live
working tree and passes the focused Stage I helper subset (`84 passed, 35
deselected`). The implementation:

- permits up to four active distinct `R03` through `R16` case lanes;
- preserves at most one prepared packet globally while allowing already
  submitted lanes to remain active;
- authenticates allowed scheduler overlap against retained submitted
  reservation manifests and rejects unknown user jobs;
- authorizes `1/2/4`-node profiles for `R04` through `R15`, `1/2` for `R16`,
  one node for `R03`, and eight exclusive nodes for `R17`;
- retains summed reservation accounting, durable replay validation, root-lock
  serialization, and R17-last enforcement.

This is staged source code, not a production promotion. Before a shared-root
writer uses the transition, commit, independently review, archive, checksum,
catalog, and verify the helper revision through the ordinary production
provenance gates. Keep the active retained campaign root unchanged until that
promotion packet is complete.

## 2. Governing Documents And Boundaries

Read these documents before mutating the production root:

1. `docs/cgl_lf_phase_i_handoff.md`
2. `docs/cgl_lf_mks24_stage_i_protocol_review.md`
3. `docs/cgl_lf_weak_guide_manuscript_plan.md`
4. `docs/cgl_lf_mks24_reproduction_implementation_plan.md`
5. `inputs/cgl_lf_paper/mks24_stage_i_manifest.json`

This completion plan is subordinate to the retained controller and frozen
matrix. It is deliberately shorter and more operational than the historical
evidence ledgers. Update it when a decision changes the execution strategy,
but do not duplicate every accepted-segment hash already retained by the
controller.

The production roots are:

```text
source repository:
  /autofs/nccs-svm1_home2/dfielding/athenak-df

shared campaign root:
  /lustre/orion/ast207/proj-shared/dfielding/CGL
```

The active source branch is:

```text
feature/cgl-landau-fluid
```

Production packets must pin the validated frozen E03 source tree rather than
the moving live checkout:

```text
/autofs/nccs-svm1_home2/dfielding/athenak-cgl-e03-9e075422
```

The current complete-history source bundle retained by the submitted R03
packet is:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL/source-archives/athenak-feature-cgl-through-dbe7e5004.bundle

SHA-256:
  b8437f066f8391a696efaaaf0de531430a9dac27c95dfc38aa7328dd49fb19fe
```

The Stage I controller is:

```text
scripts/frontier/cgl_lf_stage_i.py
```

The retained recost lifecycle companion is:

```text
scripts/frontier/cgl_lf_stage_i_checkpoint.py
```

Stage II guide-field and background-collisionality extensions remain outside
this plan. Do not start them before the accepted R17 campaign bundle and
Stage I scientific report exist.

## 3. Current Snapshot

### 3.1 Closed Work

The corrected `E03-forcing-policy` executable is qualified for Frontier
production:

```text
revision:
  9e07542281e4e6d125582f253df3ad2e3b8b154d

executable SHA-256:
  68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c
```

Retained `g024` through `g031` evidence qualifies:

- explicit `mks24_alfvenic_perpendicular` forcing;
- explicit `mks24_random_unprojected` forcing;
- modal restart identity across an OU refresh;
- one-rank versus eight-rank decomposition identity;
- passive-Delta semantics;
- nonlinear hard-wall behavior;
- standard-layout ranked startup, memory fit, and output sizing.

Corrected-E03 `R02` is accepted from fresh `t = 0` through exact `t = 10`.
Its 29 accepted segments are recorded, reconciled, assembled once as a
whole-case bundle, and analyzed successfully:

```text
runs/mks24-stage-i/E03-forcing-policy/bundles/R02
```

The F-112 recost publication projected:

```text
mapped-matrix envelope:
  900.000000 node-hours

projected use:
  874.3677783333333 node-hours

projected margin:
  25.632221666666737 node-hours

incremental project ceiling:
  4000.000000 node-hours
```

### 3.2 R03 Diagnostic Stop Line

Fresh `R03/s00_rankio_t0_t0p5` was submitted as Frontier job `4762472`.
Slurm reports `COMPLETED 0:0`, but AthenaK reached its wall-clock guard at:

```text
t = 0.31282347945569927
```

instead of:

```text
t = 0.5
```

The ranked outputs are present and the strict LF failure counters remain zero.
The current controller does not authorize the terminal restart siblings as
continuation sources because their explicit parameter-dump marker is rounded:

```text
restart_time = 0.312823
```

The hardened inspector correctly rejects the `4.79e-7` mismatch against its
`1e-10` restart-marker tolerance. Do not relax that check and do not enable a
production marker bypass.

The checkpoint payload itself is not rounded. A read-only audit found the exact
binary double `0.31282347945569927` in every one of the eight terminal
`.00001.rst` siblings. AthenaK writes `pm->time` separately into the restart
header in `src/outputs/restart.cpp` and reloads that binary value into
`Mesh::time` in `src/mesh/build_tree.cpp`. The R03 stop is therefore an
audit-metadata defect with a recoverable checkpoint candidate, not evidence of
a scientifically invalid checkpoint. Recovery still requires an explicit
fail-closed controller transition before any continuation.

The current shared-root reconciliation is internally consistent:

```text
ledger rows:         29
manifests:           30
reservations:        30
active reservations: 1
transactions:         0
issues:              []
```

The one active reservation belongs to unrecorded job `4762472`.

### 3.3 External Scientific Evidence

The retained correspondence:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL/Majeski_AthenaK_turbulence_driver.pdf

SHA-256:
  8f30cd53f630bdc2e4af181fcdb64c092e02279bf7d952e46d1dc0ce764f674b
```

records Stephen Majeski's independent AthenaK turbulence-driver run. His
beta-10 mixed-forcing calculation shows the expected:

- magnetic-field-strength distribution;
- `beta Delta` distribution;
- instability-threshold occupancy;
- Kolmogorov kinetic and magnetic spectra;
- active-run perpendicular pressure balance;
- absence of that enforced balance in the passive run;
- suppression of `bb:grad(u)` without suppression of `u_parallel`.

This is strong directional evidence that the implemented physics behaves as
expected. It justifies an ambitious production cadence. It does not replace
the retained Stage I acceptance checks or paper-panel comparisons.

## 4. Definition Of Completion

Stage I is complete only when all of the following are true:

1. Job `4762472` is retained and accounted either as a binary-authenticated
   `clean_partial` under the reviewed recovery transition or as an `aborted`
   diagnostic segment followed by a fresh launch.
2. Every mapped case `R03` through `R16` reaches exact `t = 10` through
   authenticated accepted segments.
3. `R17` reaches exact `t = 10` last, on eight nodes, through authenticated
   accepted segments.
4. Every mapped case `R02` through `R17` has one accepted whole-case bundle.
5. `bundle-campaign` succeeds for the frozen 16-case matrix.
6. The campaign-level paper analyzer succeeds and regenerates the panel-status
   table.
7. Comparison-ready panels are reviewed against their retained references.
8. Reference-blocked panels remain explicitly disclosed unless author data,
   archive data, donor diagnostics, or a qualified conversion unblock them.
9. The final ledger, reservation store, transactions, source archives,
   controller provenance, scheduler evidence, and analysis products reconcile.
10. The durable handoff is refreshed with the exact R17 completion boundary.

## 5. Execution Philosophy

### 5.1 Move Quickly Where The Evidence Is Strong

Use the already qualified corrected-E03 executable for production unless a
measured trigger requires a controlled executable transition. Avoid inserting
a new production binary into the middle of Stage I merely to improve a
recoverability path that the current controller already rejects safely.

The known restart marker issue affects the text audit metadata emitted at
non-round, wall-clock-limited partial endpoints. It does not round the binary
checkpoint payload and it does not invalidate exact-target segments. `R02`
demonstrated that exact bounded targets can be inspected, accepted, and
continued reliably.

### 5.2 Keep Shared-Root Mutations Serial

Only one lead execution agent may mutate the shared production root or issue a
Stage I submission. Every controller mutation remains serialized under the
canonical root lock. This does not require Frontier compute to remain serial.

Keep the R03 recovery lane exclusive until its binary-aware disposition,
continuation or fallback prefix, recost, and controller transition are closed.
After the bounded-concurrency transition in Section 9.1 is promoted, run
independent `R03` through `R16` case lanes concurrently. Preserve one in-flight
segment per case lineage. Drain all lower-resolution lanes before starting
`R17`, which remains an exclusive final lane.

Parallel subagents should do everything else:

- code review and patch development;
- regression testing;
- immutable-build preparation;
- concurrent-lane packet preparation;
- readiness-packet review;
- read-only scheduler and output monitoring;
- scientific health review;
- budget and storage recosting;
- case-bundle analysis;
- panel-status review;
- documentation and source-archive preparation.

### 5.3 Prefer Objective Stop Lines

Pause the production lane only for:

- nonzero strict LF failure counters;
- nonfinite synchronized histories;
- missing ranked products;
- forcing-work closure failure;
- conserved-mass or normalized-`divB` regression;
- checkpoint authentication failure;
- controller reconciliation issue;
- unresolved scheduler state;
- budget or storage exhaustion;
- evidence that a segment target repeatedly cannot fit the wall-clock guard.

Do not invent new manual review ceremonies for already automated checks.

## 6. Parallel Subagent Operating Model

### 6.1 Roles

Use a lead agent plus parallel experts with explicit ownership boundaries.

Every subagent works under the same quality charter:

- make surgical changes that match the existing codebase;
- reuse existing helpers, controller contracts, and analysis conventions;
- distinguish physics evidence from operational evidence;
- reason from the CGL, MHD, Landau-fluid, forcing, and limiter contracts;
- keep claims no broader than retained evidence;
- return concise artifacts that another expert can audit quickly.

| Role | Primary responsibility | May mutate shared production root? |
| --- | --- | --- |
| Lead execution agent | Controller lifecycle, submit decisions, integration, final authorization | Yes, exclusively |
| Restart-hardening engineer | Binary-aware controller transition, precision patch, tests, sidecar build, qualification design | No |
| Concurrency engineer | Multi-reservation policy, authenticated queue allowlist, allocation-profile replay, focused tests | No |
| Lane scheduler | Wave composition, per-case packet queue, lane-cap recommendation, completion drain | Read-only |
| Controller reviewer | Adversarial review of manifests, reservations, transactions, helper behavior | Read-only |
| Frontier operations reviewer | Queue, Slurm, allocation, storage, and timing review | Read-only |
| Plasma-physics reviewer | LF counters, limiter behavior, forcing semantics, invariants, scientific interpretation | Read-only |
| Segment evidence reviewer | Independent inspection regeneration and endpoint review | Read-only |
| Recost and budget reviewer | Throughput estimates, matrix projection, reservation recommendation | Read-only |
| Case-analysis agent | Whole-case bundle and analyzer review after each case reaches `t = 10` | Only case-analysis paths when authorized |
| Panel-comparison agent | Cross-case figure readiness and retained-reference comparison | Read-only until campaign analysis is authorized |
| Scientific-criteria agent | Uncertainty-aware pass/fail criteria for admitted reference comparisons | Repository documentation only when authorized |
| Reference-unblock agent | Author tables, Fourier normalization, donor diagnostics, and archive-data requests | Read-only |
| Provenance agent | Commit packet, source bundle, checksums, catalog update, durable handoff draft | Repository documentation only when authorized |

Use more than one independent reviewer for transitions that alter executable
provenance, controller behavior, the Stage I reservation ceiling, or the R17
launch profile. Use one reviewer for ordinary accepted exact-boundary
continuations once the profile is stable.

### 6.2 Agent Waves

Run these work streams concurrently whenever possible:

| While the lead lane is doing this | Parallel subagent work |
| --- | --- |
| Dispositioning job `4762472` | Binary-aware controller patch, precision patch, fixture regression, salvaged R03 `.5` packet, fresh `.25` fallback packet, updated recost draft |
| Running the accepted R03 prefix | Concurrency patch, multi-reservation regression, authenticated-queue regression, `1/2/4`-node scaling packets, wave composition |
| Running one R03-R16 wave | Per-lane output monitoring, completed-lane inspection drafts, next-packet preparation, provisional recost, storage review, completed-case analysis |
| Draining a wave barrier | Authoritative recost publication, allocation-profile update, storage audit, lane-cap ratchet review, next-wave packet review |
| Running the final R04-R16 lanes | Final R17 readiness review, campaign-bundle dry-run, 64-rank storage audit, final-analysis staging |
| Running R17 | Campaign-bundle assembly preparation, panel comparison review, final handoff draft |

### 6.3 Mutation Ownership

The lead agent owns these commands against the shared root:

```text
prepare
check-submit
submit
recover-submit
clear-submit-pending
inspect-segment
record
cancel
recover-transactions
bundle-case
bundle-campaign
summary
```

The command-safety classification is:

| Class | Commands |
| --- | --- |
| Read-only controller actions | `validate-matrix`, `check-submit`, `reconcile` |
| Shared-root controller writes | `init`, `approve-qualification`, `prepare`, `submit`, `mark-submitted`, `recover-submit`, `clear-submit-pending`, `inspect-segment`, `record`, `bundle-case`, `bundle-campaign`, `cancel`, `summary`, `recover-transactions` |
| Read-only system inspection | `squeue`, `sacct`, `find`, `du`, `sha256sum`, `git status`, `git diff` |

Subagents may run read-only commands such as:

```text
squeue
sacct
find
du
sha256sum
git status
git diff
reconcile
```

Subagents may develop code or documentation only in disjoint repository write
sets assigned by the lead agent.

Analysis tools also write outputs. Give an analysis subagent an explicit
directory lease or a scratch output directory before it runs `paper-analyze`.
Only the lead agent invokes controller writers against retained shared-root
paths.

## 7. Binary-Aware Recovery And Writer Precision Strategy

### 7.1 Production Fast Lane

Keep Stage I on the already qualified corrected-E03 executable. Run one focused
controller-hardening sprint, salvage job `4762472` if every binary-aware gate
passes, and continue R03 from its exact binary checkpoint. If any salvage gate
fails, account the job as `aborted` and relaunch fresh immediately.

The diagnostic run measured:

```text
6600 seconds / 0.31282347945569927 = 21098.161 seconds per simulated time unit
```

The inferred targets are:

```text
0.31282347945569927 -> 0.50: approximately 3949 seconds
0.25 simulated time units:          approximately 5275 seconds
0.50 simulated time units:         approximately 10549 seconds
```

The preferred next packet is the salvaged exact continuation to `t = 0.5`.
After that, use absolute output-aligned quarter-unit endpoints unless measured
evidence supports larger exact increments. The fresh fallback starts at
`t = 0` and targets `t = 0.25`.

For the qualified E03 executable:

- continue only from inspector-authenticated complete restart sibling sets;
- accept ordinary exact-target segments through the existing strict path;
- accept a non-round `clean_partial` only through the reviewed binary-aware
  path;
- reduce the next exact target immediately after any wall-clock partial.

### 7.2 Controller Salvage Hardening

Extend:

```text
scripts/frontier/cgl_lf_stage_i.py
```

with a narrow restart-header authenticator bound to the qualified E03
executable ABI and retained build manifest. It must:

1. parse the binary `Mesh::time` field from each selected restart sibling;
2. require one complete ranked sibling set and exact sibling agreement;
3. require the terminal binary time to match final synchronized history within
   the existing strict tolerance;
4. retain and parse the text `time/restart_time` marker;
5. require either full-precision marker agreement or exact agreement with the
   legacy default-precision serialization of the authenticated binary time;
6. reject unknown executable revisions, unknown ABI layouts, absent markers,
   sibling disagreements, nonfinite values, and arbitrary tolerance bypasses.

Changing the live helper also invalidates the strict live-helper authentication
recorded by the already submitted R03 manifest. Add a narrow, checksum-bound,
tested transition for this one historical submitted segment: retain the
original prepared helper bytes, authenticate them from their retained bundle,
authenticate the promoted helper transition, and preserve the full audit
chain. Do not introduce a general submitted-manifest bypass.

### 7.3 Parallel Writer Repair Lane

Develop the restart-marker precision repair concurrently. The minimum patch is
local to restart metadata:

```text
src/outputs/restart.cpp
```

Format `pm->time` with:

```cpp
std::setprecision(std::numeric_limits<Real>::max_digits10)
```

and store the formatted marker through `ParameterInput::SetString`.

Prefer this scoped repair over changing `ParameterInput::SetReal` globally.
The global helper has many callers and would broaden provenance churn without
improving the immediate Stage I physics path. If a later audit identifies
multiple restart-critical real-valued fields, introduce a reusable exact-real
setter as a separate reviewed change.

Do not switch the E03 production executable merely to gain the writer repair.
The binary-aware controller path is the immediate recovery tool; the writer
repair is qualification-ready code for a deliberate future executable
transition.

### 7.4 Required Hardening Tests

Add and run:

1. A controller fixture with final time `0.31282347945569927` proving:
   - a binary-exact sibling set with legacy text marker `0.312823` is accepted
     only by the binary-aware path;
   - full-precision text and binary agreement is accepted;
   - a text marker inconsistent with legacy serialization is rejected;
   - binary sibling disagreement is rejected;
   - terminal-history disagreement is rejected;
   - an unknown ABI or executable revision is rejected;
   - `record --result clean_partial` succeeds only after valid inspection.
2. A historical submitted-manifest fixture proving the checksum-bound helper
   transition permits exactly the reviewed R03 recovery shape and rejects a
   broad bypass.
3. A CGL-LF restart regression that parses emitted `time/restart_time` and
   requires full-precision round-trip agreement with terminal history.
4. The modal forcing restart regression across an OU refresh.
5. Focused CPU, MPI CPU, turbulence-driver CPU, and turbulence-driver MPI CPU
   suites.
6. A compact Frontier one-node, eight-rank restart qualification:
   - uninterrupted comparator;
   - arbitrary non-round wall-clock partial;
   - resumed continuation through an OU refresh.

### 7.5 Executable Promotion Trigger

Do not promote the repaired executable into Stage I production by default.
The current E03 qualification token is immutable after segment preparation,
and a mid-campaign binary replacement would require an explicit controller
and provenance transition.

Promote the repair only if one of these occurs:

1. the ABI-bound binary-aware recovery path cannot support the needed partial
   restarts;
2. R17 cannot reliably complete useful exact-target increments;
3. a new audit finds another restart-critical precision defect that affects
   accepted exact endpoints;
4. the team deliberately chooses a versioned qualification transition because
   its operational benefit exceeds the transition cost.

If promotion becomes necessary, create a reviewed new execution epoch rather
than rewriting E03 history. The current controller does not implement
cross-epoch campaign bundling. Before any E04 promotion, explicitly implement,
test, qualify, and document an accepted-case import path that preserves
completed E03 `R02` and deliberately assembles E03 and E04 case bundles.

## 8. Immediate Recovery: Job 4762472 And R03

### 8.1 Run The Bounded Salvage Decision

The lead agent should retain the job `4762472` log, rank-local outputs,
restart siblings, and scheduler evidence while the controller and test agents
complete Section 7.2 through 7.4. The read-only evidence already establishes:

```text
Slurm state:                 COMPLETED 0:0
charged elapsed:             6613 seconds
charged node-hours:          1.836944
Athena terminal time:        0.31282347945569927
strict LF counters:          zero
terminal restart siblings:   eight complete .00001.rst files
binary Mesh::time siblings:  exact 0.31282347945569927
text restart_time siblings:  rounded 0.312823
```

Promote the narrow controller transition only after focused tests, independent
review, committed helper bytes, retained source archive, checksum cataloging,
and publication audit. Then inspect the existing segment through the
binary-aware path.

If inspection passes, record job `4762472` as `clean_partial` and state that:

- AthenaK stopped cleanly on its wall-clock guard;
- the complete rank-local restart sibling set is binary-authenticated;
- text markers match the exact legacy serialization of the binary time;
- strict LF counters are zero;
- same-executable E03 continuation is authorized.

If any requirement fails, record job `4762472` as `aborted`, state the precise
failed gate, and prohibit continuation from its outputs. Do not turn the
focused salvage sprint into an open-ended delay.

After either disposition, reconcile until:

```text
active reservations: 0
transactions:        0
issues:              []
```

Do not delete the diagnostic outputs.

### 8.2 Preferred R03 Continuation

After an authenticated `clean_partial` record, prepare a uniquely named R03
continuation from the retained terminal sibling set:

```text
start:           t = 0.31282347945569927
target:          t = 0.5
nodes:           1
Slurm walltime:  02:00:00
Athena timeout:  01:50:00
expected time:   approximately 3949 seconds
```

Require the ordinary controller preflight:

- no queued user jobs;
- free strict root lock;
- no pending transactions;
- authenticated source bundle;
- reviewed stale shared-root campaign acknowledgement;
- clean reconciliation;
- committed helper bytes.

After the accepted exact `t = 0.5` prefix:

1. record and reconcile immediately;
2. generate an updated timing and budget projection;
3. choose the largest measured successor target that retains useful wall-clock
   margin;
4. remain on absolute quarter-unit endpoints unless measurements justify a
   larger exact increment.

### 8.3 Fail-Fast Fresh Fallback

If the binary-aware transition or R03 inspection fails, prepare a new uniquely
named fresh `R03` segment from `t = 0` with:

```text
target:          t = 0.25
nodes:           1
Slurm walltime:  02:00:00
Athena timeout:  01:50:00
expected time:   approximately 5275 seconds
```

Use the same preflight and recost lifecycle. Continue from accepted
output-aligned quarter-unit endpoints.

## 9. Concurrent Production Cadence For R03 Through R16

### 9.1 Bounded-Concurrency Transition

The historical promoted helper intentionally enforces an earlier conservative
policy:

- exactly one active prepared or submitted reservation globally;
- rejection of submission whenever any user job is queued;
- exactly one node for every canonical `R02` through `R16` segment;
- exactly eight nodes for `R17`.

The surgical transition is implemented and locally verified ahead of the R03
recovery gate so it does not remain on the critical path. After the first
accepted R03 continuation or fallback prefix and reviewed recost, promote the
reviewed controller revision that:

1. permits multiple active reservations only when they belong to distinct
   `R03` through `R16` case lanes;
2. preserves at most one prepared or submitted segment per case lineage;
3. serializes every metadata mutation under the existing canonical root lock;
4. accounts the sum of all active reservation node-hours against the reviewed
   Stage I envelope and the `4000` node-hour project ceiling;
5. permits queued or running Stage I jobs only when each one matches an
   authenticated submitted reservation and exact retained manifest;
6. continues to reject unknown user jobs, unreviewed shared-root campaign
   records, orphaned run directories, pending transactions, and ambiguous
   scheduler state;
7. allows reviewed multi-node allocation profiles for `R04` through `R16`;
8. keeps `R17` fixed at eight nodes, requires accepted `t = 10` predecessors,
   and permanently locks out lower-resolution preparation after R17 starts;
9. replays historical one-node records without rewriting their provenance.

Complete focused controller fixtures proving:

- two distinct cases may be prepared, submitted, inspected, and recorded in
  either completion order;
- a second active segment for the same case is rejected;
- summed reservations and the lane cap fail closed;
- an exact authenticated Stage I queue set is accepted while one unknown user
  job is rejected;
- a pending transaction or orphaned run directory still blocks mutation;
- approved `1`, `2`, and `4` node profiles replay correctly where authorized,
  while unapproved or decomposition-infeasible profiles are rejected;
- R17 remains exclusive and last.

The local coverage tranche proves node-profile authorization, duplicate case
rejection, the four-lane cap, one-prepared-packet enforcement, R17 exclusivity,
authenticated submitted-job queue overlap, unknown-job rejection, and durable
preservation of an unrelated submitted lane while a completed lane is
recorded. Commit, push, archive, checksum, catalog, and independently review
this transition before the first concurrent wave.

### 9.2 Case-Lane Lifecycle

For every segment in every active case lane:

1. The lane scheduler identifies the latest accepted exact endpoint.
2. The timing reviewer proposes the largest evidence-backed next target and
   reviewed node profile.
3. The controller reviewer audits the bounded readiness packet.
4. The lead agent runs `prepare`, `check-submit`, and `submit` serially under
   the root lock.
5. While Frontier lanes run concurrently, parallel agents monitor output
   growth, scheduler state, storage, and the next packet for each lane.
6. After any lane completes, the lead agent runs `inspect-segment`.
7. The segment evidence reviewer independently checks endpoint time, ranked
   inventory, strict LF counters, forcing-work closure, and restart siblings.
8. The lead agent runs `record` and `reconcile`.
9. The recost reviewer updates the provisional case-family throughput model.
10. The lead agent prepares the next packet for that case immediately when its
    retained evidence is coherent and the wave policy allows refill.

Controller writes remain serial. Frontier allocations, read-only monitoring,
packet preparation, and case analysis should overlap aggressively.

The retained recost publication companion currently requires zero active
reservations. Preserve that fail-closed rule initially: update provisional
models continuously, then execute the authoritative staged verify, independent
audit, promote, verify-promoted, publication-audit, archive, and reconcile
lifecycle at planned drained-wave barriers. Harden the companion for
checksum-bound in-flight reservation snapshots only if barrier cost becomes a
measured bottleneck.

### 9.3 Target And Node Selection

Choose the largest exact increment whose measured or family-inferred runtime
remains comfortably below the `6600`-second Athena guard. Use a review
threshold of approximately `5700` to `6000` seconds. Ratchet upward promptly
when an accepted prefix supports it and downward immediately after an early
wall-clock termination.

Treat node count as a measured optimization variable:

- keep the salvaged or fallback R03 lane on one node while it closes the
  immediate recovery and cost-calibration gate;
- qualify `1`, `2`, and `4` node profiles promptly for the standard-layout
  `R04` through `R15` cases;
- qualify `1` and `2` node profiles for lower-resolution `R16`, and expand to
  `4` only if meshblock decomposition and measured efficiency justify it;
- permit a larger R04-R16 profile only after explicit decomposition,
  throughput, I/O, budget, and replay review;
- keep R17 at its separately qualified eight-node profile.

Optimize time-to-R17, not node-hour minimization in isolation. Prefer a larger
profile when it materially reduces elapsed time while retaining useful
node-hour efficiency, balanced meshblock ownership, complete rank-local
products, acceptable I/O pressure, and budget headroom. Prefer another
concurrent case lane over weak strong-scaling gains.

Do not mechanically reuse the R02 half-unit profile for beta-100, random,
compressive, heat-flux-extreme, or finite-limiter cases.

### 9.4 Initial Aggressive Profiles

These are first-packet proposals, not permanent limits. Recost provisionally
after each first accepted prefix and authoritatively at wave barriers.

| Case | Runtime family | Initial exact increment | Initial node study | Reason |
| --- | --- | ---: | --- | --- |
| `R03` | active Alfvenic beta-100 hard wall | Salvage to `0.5`, then `0.25` | `1` | Directly inferred from job `4762472` |
| `R04` | active random beta-10 | `0.25` | `1`, `2`, `4` scaling leader | Establish standard-layout multi-node profile and random-forcing cost |
| `R05` | active random beta-100 | `0.25` | Reuse reviewed standard winner | Beta-100 and random forcing both merit measured calibration |
| `R06` | passive Alfvenic beta-10 | `0.50` | Reuse reviewed standard winner | Closest passive counterpart to completed R02 |
| `R07` | passive Alfvenic beta-100 | `0.25` | Reuse reviewed standard winner | Reuse R03 beta-100 bound until measured faster |
| `R08` | passive random beta-10 | `0.25` | Reuse reviewed standard winner | Reuse R04 random-family calibration, then ratchet |
| `R09` | passive random beta-100 | `0.25` | Reuse reviewed standard winner | Reuse beta-100 random-family calibration |
| `R10` | compressive active random beta-1 | `0.25` | Recheck reviewed standard winner | New compressive family |
| `R11` | compressive active random beta-100 sonic-correlation | `0.25` | Recheck reviewed standard winner | New beta-100 sonic-correlation family |
| `R12` | stronger LF heat flux | `0.25` | Recheck reviewed standard winner | Extreme closure coefficient needs measured calibration |
| `R13` | weaker LF heat flux | `0.25` | Recheck reviewed standard winner | Distinct closure-cost profile |
| `R14` | beta-100 finite limiter `nu_lim = 20` | `0.25` | Reuse reviewed standard winner | New finite-limiter profile |
| `R15` | beta-100 finite limiter `nu_lim = 200` | `0.25` | Reuse reviewed standard winner | New finite-limiter profile |
| `R16` | beta-10 `96 x 96 x 192` scale separation | `1.50` | `1`, `2`; test `4` only if useful | Retained low-resolution timing supports larger increments |

Complete each case through exact `t = 10`, bundle it once, and start analysis
as soon as that individual case closes.

### 9.5 Family-Based Learning

Use early completed cases to accelerate later cases:

- `R03` calibrates beta-100 Alfvenic hard-wall cost for `R07`, `R14`, and
  `R15`.
- `R04` calibrates corrected random-forcing cost and the standard-layout node
  profile for `R05`, `R08`, `R09`, `R10`, and `R11`.
- `R06` calibrates passive beta-10 savings.
- `R12` and `R13` establish whether LF-strength extremes alter wall-clock cost
  or the preferred node profile materially.
- `R16` verifies the low-resolution Figure 11 lane and its separate node
  profile before R17 starts.

Use measured evidence, not inherited assumptions, but avoid repeating pilot
work once a family profile is stable.

### 9.6 Wave Rollout

After the R03 recovery prefix and concurrency transition close:

1. Run compact standard-layout `1/2/4`-node and R16 `1/2`-node scaling packets
   in qualification namespaces or reviewed fresh-prefix pilot records. Do not
   fork accepted production lineages merely to benchmark scaling.
2. Start an initial four-lane packet wave using `R03`, `R04`, `R12`, and
   `R16`.
3. Let that first packet wave drain, then promote the authoritative recost,
   allocation-profile table, and storage audit.
4. Raise the lane cap from four toward six when measured scheduler behavior,
   filesystem pressure, storage growth, and budget headroom remain healthy.
5. Refill a lane immediately after its accepted segment is recorded if the
   provisional model, storage monitor, and reservation sum remain healthy.
6. Fill subsequent lanes from `R05` through `R11` and `R13` through `R15`,
   adding `R06` promptly as the passive beta-10 calibration lane.
7. Drain periodically for authoritative recost barriers and always before a
   lane-cap or allocation-profile increase.
8. Drain every lower-resolution lane, publish the final predecessor recost,
   and only then authorize R17.

## 10. R17 High-Resolution Completion

### 10.1 Entry Gate

The controller correctly requires R17 to run last. Before preparing `R17`,
require:

- accepted exact `t = 10` lineages for `R02` through `R16`;
- no active Stage I reservation;
- no pending transaction;
- clean reconciliation;
- eight-node allocation;
- fresh storage and node-hour recost;
- independently reviewed 64-rank readiness packet.

Once R17 starts, do not return to a lower-cost mapped case.

### 10.2 Retained Timing Evidence

Historical E02 timing pilots remain valid as cost evidence:

```text
R17/s00: t = 0.00 -> 0.05, 8 nodes, 00:15:50, 2.111111 node-hours
R17/s01: t = 0.05 -> 0.10, 8 nodes, 00:15:56, 2.124444 node-hours
```

The developed rate is:

```text
42.488889 node-hours per simulated time unit
```

Use this only as the initial E03 estimate. Recost from accepted corrected-E03
R17 prefixes.

### 10.3 Initial R17 Profile

Start corrected-E03 R17 fresh from `t = 0` with:

```text
nodes:           8
ranks:           64
initial target:  t = 0.25
expected time:   approximately 4780 elapsed seconds
Slurm walltime:  02:00:00
Athena timeout:  01:50:00
```

The quarter-unit initial prefix is output-aligned and should preserve more than
`1800` seconds against the Athena timeout. After the first accepted prefix,
keep quarter-unit exact increments if measured elapsed time preserves the
guard. Increase only if the output cadence and measured corrected-E03 rate
justify it.

If output cadence, filesystem behavior, or developed corrected-E03 runtime
cost is materially different, reduce the next exact increment immediately.
Do not continue from an unauthenticated wall-clock partial.

### 10.4 R17 Storage Discipline

Before R17 submission:

1. inventory free project storage;
2. refresh the measured no-pruning projection from retained corrected-E03
   output and verify available retention capacity against the larger of that
   projection and `958271710272` bytes;
3. estimate per-snapshot and per-restart 64-rank size from retained pilots;
4. confirm the required retained cadence for the `t = 8` through `10`
   analysis window;
5. retain all controller-required restart siblings;
6. avoid unnecessary debug outputs;
7. monitor storage after every segment;
8. bundle and analyze promptly after exact `t = 10`.

The historical `747543662189`-byte sequential peak assumes retention of only
the final two restart groups. The current controller authenticates retained
inventory and does not yet implement that compaction lifecycle. Do not use the
smaller number as an operational gate unless a reviewed, tested,
controller-supported compaction transition exists.

Do not delete accepted evidence merely to recover space. Change retention only
through an explicit reviewed storage plan.

Concurrent `R03` through `R16` waves also require storage discipline. Before
the first wave and at every drained-wave barrier:

1. measure retained bytes by completed case and active lane;
2. project the maximum simultaneous in-flight snapshot and restart growth;
3. include multi-node rank-local inventory growth in the projection;
4. cap or drain lanes if filesystem pressure or free-quota headroom degrades;
5. preserve every authenticated artifact unless a controller-supported
   compaction transition is implemented and tested.

## 11. Budget Strategy

### 11.1 Record Reality Early

Job `4762472` consumes `6613 / 3600 = 1.836944` node-hours. Account for it,
then update the Stage I projection from measured R03 throughput.

The F-112 projection has limited margin inside the current `900` node-hour
Stage I reservation. R03's observed rate implies that the old matrix estimate
is already likely to exceed that reservation before contingency. The
incremental project ceiling is `4000` node-hours. Recompute the matrix now and
promote a bounded reviewed reservation with useful headroom, likely in the
`1100` to `1200` node-hour range if the refreshed arithmetic confirms it. Do
not submit the next production segment under a stale `900` node-hour envelope.

### 11.2 Reservation Transition

Before the next R03 production submission:

1. have the recost agent produce a measured case-family projection;
2. have an independent reviewer audit arithmetic and contingency;
3. select a revised Stage I reservation with useful headroom;
4. patch the controller constant and documentation surgically;
5. run focused controller and recost regressions;
6. commit, push, archive, checksum, and catalog the transition;
7. resume R03 under the retained `4000` node-hour project ceiling.

Before concurrent fan-out:

1. refresh the matrix model from the accepted R03 recovery prefix;
2. include projected multi-node node-hours plus the summed maximum active-wave
   reservation;
3. promote the bounded-concurrency controller transition from Section 9.1;
4. run the expanded controller and recost regressions;
5. commit, push, archive, checksum, catalog, and independently review the
   transition;
6. start the compact scaling packets and first concurrent wave.

This is an accounting-control update, not a reason to reopen qualified
physics.

### 11.3 Recost Frequency

Update the provisional throughput and node-efficiency model:

- after accounting for job `4762472`;
- after every scaling packet;
- after the first accepted prefix of each new runtime family;
- after any early wall-clock termination;
- after every completed case.

Publish an authoritative recost under the retained zero-active-reservation
companion contract:

- after the R03 recovery prefix;
- after the first drained concurrent-wave barrier;
- before raising the lane cap or adding a larger node profile;
- at later drained-wave barriers when the projected envelope moves materially;
- before R17;
- after the first corrected-E03 R17 prefix;
- at final R17 completion.

## 12. Scientific Health Review

### 12.1 Per-Segment Mechanical Checks

Require:

- finite synchronized MHD and user histories;
- conserved mass to roundoff;
- controlled normalized `divB` wherever available;
- zero `lf_dfloor`;
- zero `lf_pfloor`;
- zero `lf_nonfin`;
- zero `lf_nonpos`;
- zero `lf_hardbd`;
- zero hard-bound volume;
- complete expected ranked snapshots;
- complete expected ranked restart siblings;
- exactly one terminal restart group whose binary time matches final history,
  with text metadata accepted only by the strict current rule or the reviewed
  legacy-serialization rule;
- active-Delta `Delta tot-E - Delta force_work` closure within the existing
  relative tolerance;
- passive-Delta forcing work retained without mislabeling it as an active-CGL
  energy closure;
- expected hard-wall projection activity where applicable.

Archive interval LF cap fractions, signed applied LF work, applied CGL
pressure work, magnetic-energy history, beta, kinetic and magnetic energies,
threshold volumes, and effective collision rate. Interpret `lf_hwproj`
correctly: it counts hard-wall projection activity and should generally be
nonzero in hard-wall production. It is not a failure counter.

### 12.2 Per-Case Physics Review

After each case reaches exact `t = 10`, review:

- density, velocity, magnetic-field, `p_parallel`, and `p_perp` histories;
- beta and `Delta p`;
- mirror and firehose occupancy;
- limiter and hard-wall activity;
- LF heat-flux cap fractions;
- forcing-work and pressure-work accounting;
- steady-window selection over `t = 8` through `10`;
- PDFs, spectra, transfer, alignment, and local-strain products relevant to
  that case.

The physics reviewer should explicitly compare active and passive behavior:

- active runs should show pressure-anisotropy feedback;
- passive runs should not enforce perpendicular pressure balance;
- active runs should suppress `bb:grad(u)` without simply suppressing
  `u_parallel`;
- random and Alfvenic forcing families should retain their qualified forcing
  semantics;
- finite-limiter runs should remain distinct from hard-wall and passive roles.

Audit the archived forcing controls once per completed case and monitor their
histories continuously:

- Alfvenic decks use `mks24_alfvenic_perpendicular`, retain `f_z = 0`, and
  remain perpendicular-divergence-free even when retained modes have
  `k_z != 0`.
- Random decks use `mks24_random_unprojected` and retain unconstrained
  three-component forcing with nonzero parallel content.
- Both forcing families retain physical `|k| / pi` shell `[1, 3]`,
  isotropic-total-`|k|` power spectrum `k^-2`, `dedt = 0.32`, and restartable
  modal state.
- The ordinary correlation time is `tcorr = 2`; `R11` intentionally uses
  sonic-correlation `tcorr = 0.2`.

### 12.3 Cross-Case Panel Review

Run cross-case comparison work concurrently as dependencies become available:

| Panels | Required cases |
| --- | --- |
| Figure 2(a) | `R02` |
| Figures 2(b), 4(a), 9 | `R02` through `R09` |
| Figure 5(b) | `R02`, `R04` |
| Figure 7 lower | `R02`, `R04`, `R06` |
| Figure 8 | `R02`, `R06` |
| Figure 11 lower | `R16`, `R02`, `R17` |
| Figure 12 upper | `R12`, `R02`, `R06`, `R13` |
| Figures 13(b), 13(d) | `R14`, `R15`, `R03`, `R07` |

The remaining dimensional or normalization-blocked panels should proceed in a
parallel reference-resolution lane. Ask for author or archive data early.
Do not block simulation production while that correspondence proceeds.

All admitted checksum-backed comparison panels still need reviewed numeric
acceptance criteria. Assign a scientific-criteria subagent immediately and
finish those uncertainty-aware criteria before each dependency set closes.
Contact Stephen Majeski immediately for raw curves, FFT-normalization details,
or donor diagnostic code so reference work cannot become the final campaign
stall.

## 13. Case Bundling And Analysis

### 13.1 Bundle Each Completed Case Immediately

After one case reaches accepted exact `t = 10`:

1. run `bundle-case`;
2. authenticate the lineage and merged histories;
3. run the paper analyzer once for that whole-case bundle;
4. retain diagnostics and generated figures;
5. review the exact `t = 8` through `10` steady window;
6. update panel dependencies;
7. start the next production case without waiting for every scientific
   interpretation to finish.

This overlaps compute and analysis cleanly.

### 13.2 Final Campaign Bundle

After accepted R17 completion:

1. run `bundle-case` for R17;
2. run `bundle-campaign --required-final-time 10`;
3. run the campaign paper analyzer;
4. regenerate panel status;
5. retain the final campaign manifest, diagnostics, figures, and reference
   provenance;
6. audit all 16 case lineages;
7. refresh the durable handoff.

## 14. Command Skeleton

Use exact reviewed paths and acceptance prose from the current retained
packet. The following is a skeleton, not a substitute for packet review.

```bash
ROOT=/lustre/orion/ast207/proj-shared/dfielding/CGL
REPO=/autofs/nccs-svm1_home2/dfielding/athenak-df
FROZEN_SOURCE=/autofs/nccs-svm1_home2/dfielding/athenak-cgl-e03-9e075422
STAGE="$REPO/scripts/frontier/cgl_lf_stage_i.py"
MATRIX="$REPO/inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
SOURCE_BUNDLE="$ROOT/source-archives/athenak-feature-cgl-through-dbe7e5004.bundle"

# Rotate SOURCE_BUNDLE only after the replacement complete-history bundle is
# committed, cataloged, checksum-validated, and independently reviewed.

python3 "$STAGE" --root "$ROOT" reconcile

# Preferred route after binary-aware inspection and clean_partial accounting.
python3 "$STAGE" --root "$ROOT" prepare \
  --case-id R03 \
  --segment s01_rankio_t0p312823_t0p5 \
  --acceptance-criterion "<reviewed exact-boundary acceptance prose>" \
  --executable "<qualified executable>" \
  --build-manifest "<qualified build manifest>" \
  --source-dir "$FROZEN_SOURCE" \
  --source-bundle "$SOURCE_BUNDLE" \
  --matrix "$MATRIX" \
  --restart-file "<authenticated terminal R03 .00001 restart sibling>" \
  --nodes 1 \
  --walltime 02:00:00 \
  --athena-walltime 01:50:00 \
  --override time/tlim=0.5

python3 "$STAGE" --root "$ROOT" check-submit \
  --manifest "<prepared manifest>" \
  --allow-shared-root-campaign beta25-accel05-gamma10001-purecgl-256

python3 "$STAGE" --root "$ROOT" submit \
  --manifest "<prepared manifest>" \
  --allow-shared-root-campaign beta25-accel05-gamma10001-purecgl-256

python3 "$STAGE" --root "$ROOT" inspect-segment \
  --manifest "<submitted manifest>" \
  --required-time "<exact target>"

python3 "$STAGE" --root "$ROOT" record \
  --manifest "<submitted manifest>" \
  --job-id "<slurm job id>" \
  --result accepted \
  --notes "<reviewed endpoint summary>"

python3 "$STAGE" --root "$ROOT" reconcile
```

For job `4762472`, prefer `record --result clean_partial` after the reviewed
binary-aware transition and successful inspection. If any salvage gate fails,
use `record --result aborted` with explicit non-continuation notes and prepare a
fresh uniquely named R03 `t = 0` to `0.25` fallback without `--restart-file`.

After the bounded-concurrency transition is promoted, issue `prepare`,
`check-submit`, and `submit` serially for each selected `R03` through `R16`
lane, but allow authenticated submitted allocations to overlap on Frontier.
Use:

```bash
  --nodes "<reviewed lane allocation profile>"
```

Keep R03 on one node. Do not hard-code one node for `R04` through `R16`.
Retain one active segment per case, the reviewed total lane cap, summed
reservation accounting, and exact manifest binding for every queued or running
Stage I job.

## 15. Durable Checkpoints

Commit, push, bundle, checksum, and catalog at these boundaries:

1. disposition of job `4762472` and authorization of the salvaged R03
   continuation or fresh quarter-unit fallback;
2. bounded-concurrency, allocation-profile, and Stage I reservation transition;
3. every drained-wave authoritative recost and lane-cap change;
4. any executable qualification or execution-epoch transition;
5. completion of each mapped case;
6. R17 readiness authorization;
7. R17 completion;
8. final campaign bundle and scientific report.

For ordinary exact-boundary continuations within one stable profile, retain
controller evidence and recost publications without turning every segment
into a large documentation rewrite.

## 16. Final Checklist

### Immediate

- [ ] Implement, test, review, archive, and catalog the binary-aware controller
      recovery transition.
- [ ] Inspect job `4762472` through the binary-aware path.
- [ ] Record job `4762472` as `clean_partial` if every salvage gate passes;
      otherwise record it `aborted` and activate the fresh fallback.
- [ ] Reconcile to zero active reservations and `issues = []`.
- [ ] Prepare and review salvaged R03 `t = 0.31282347945569927` to `0.5`, or
      fresh fallback `t = 0` to `0.25`.
- [ ] Start the scoped restart-marker precision patch in parallel.
- [ ] Recast and promote the reviewed Stage I budget reservation before the
      next submission.
- [x] Implement and run the first focused regression tranche for the
      bounded-concurrency controller transition.
- [ ] Independently review, archive, catalog, and promote the
      bounded-concurrency controller transition.
- [ ] Qualify standard-layout `1/2/4`-node and R16 `1/2`-node scaling packets.
- [ ] Request raw curves, FFT-normalization details, or donor diagnostics from
      Stephen Majeski in parallel.

### Standard Matrix

- [ ] Launch the initial four-lane `R03`, `R04`, `R12`, `R16` wave.
- [ ] Promote the first drained-wave recost, storage audit, and node-profile
      table.
- [ ] Raise the lane cap toward six if measured evidence supports it.
- [ ] Complete and bundle `R03`.
- [ ] Complete and bundle `R04`.
- [ ] Complete and bundle `R05`.
- [ ] Complete and bundle `R06`.
- [ ] Complete and bundle `R07`.
- [ ] Complete and bundle `R08`.
- [ ] Complete and bundle `R09`.
- [ ] Complete and bundle `R10`.
- [ ] Complete and bundle `R11`.
- [ ] Complete and bundle `R12`.
- [ ] Complete and bundle `R13`.
- [ ] Complete and bundle `R14`.
- [ ] Complete and bundle `R15`.
- [ ] Complete and bundle `R16`.

### High Resolution

- [ ] Recost and independently review R17.
- [ ] Verify eight-node storage readiness.
- [ ] Complete and bundle `R17`.

### Final Analysis

- [ ] Assemble the 16-case campaign bundle.
- [ ] Run final paper analysis.
- [ ] Review comparison-ready panels.
- [ ] Disclose or resolve reference-blocked panels.
- [ ] Reconcile the final ledger, archives, and provenance.
- [ ] Refresh the durable handoff at accepted R17 completion.

## 17. Completion Decision

The campaign should proceed aggressively. The code has broad retained tests,
corrected forcing-policy qualification, one complete corrected-E03 production
case, a controller that fails closed, and encouraging independent physics
evidence from Stephen Majeski. The shortest robust path is:

1. promote the narrow binary-aware controller transition and account for R03
   as an authenticated `clean_partial`, or fail fast to the fresh fallback;
2. continue salvaged R03 to exact `t = 0.5`, then use output-aligned
   quarter-unit endpoints;
3. promote bounded multi-case concurrency and reviewed multi-node profiles;
4. run R03 through R16 as four-to-six concurrent case lanes while serializing
   shared-root mutations and preserving one in-flight segment per case;
5. use wave barriers for authoritative recost, storage review, and lane-cap
   ratchets while analysis continues in parallel;
6. execute R17 last with eight-node, recosted exact increments;
7. assemble and review the final campaign bundle.
