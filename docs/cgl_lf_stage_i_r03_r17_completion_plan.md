# CGL-LF MKS24 Stage I Rapid Completion Plan: R03 Through R17

## 1. Purpose

This document is the execution-focused plan for completing the corrected
AthenaK CGL-LF MKS24 Stage I matrix from the drained first accelerated wave
through accepted `R17` completion and the final campaign analysis bundle.

The goal is to move quickly. The implementation is already scientifically
credible, the corrected forcing policies are qualified, `R02` is complete,
and independent use by Stephen Majeski provides encouraging external evidence.
The remaining work is a production campaign plus the current F116/F117 release
gate, not an open-ended physics-development program.

This plan preserves fail-closed controller behavior while using parallel
subagents aggressively for software hardening, packet review, scientific
monitoring, recosting, analysis, provenance, and documentation. The first
accelerated wave demonstrated the bounded-concurrency model: shared-root
metadata mutations remain serial, but independent Frontier case lanes may run
concurrently.

### 1.1 Acceleration Transition Status

The bounded-concurrency controller transition is active and was exercised by
the first accelerated wave across four independent lanes and ten concurrent
nodes. The combined implementation:

- permits up to four active distinct `R03` through `R16` case lanes;
- preserves at most one prepared packet globally while allowing already
  submitted lanes to remain active;
- authenticates every queued `cgl_` job against retained submitted reservation
  manifests, permits well-formed unrelated non-CGL jobs, and rejects unbound
  or malformed Stage I jobs;
- authorizes `1/2/4`-node profiles for `R04` through `R15`, `1/2` for `R16`,
  one node for `R03`, and eight exclusive nodes for `R17`;
- retains summed reservation accounting, prospective durable replay
  validation, root-lock serialization, and retained R17-last enforcement.

The next release boundary is not yet complete. Before the next shared-root
production mutation, finish the reviewed release commit and publish:

- `F-116`, the independently reviewed current-source-authority supersession
  that selects the final committed seven-tool release and complete-history
  source bundle; then
- `F-117`, the independently reviewed recost and next-wave recommendation,
  including the fresh R12 rerun.

Until both publications are promoted and verified, retain the drained canonical
campaign state and do not claim the next wave is authorized.

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

The active source branch and pre-F116 source boundary are:

```text
branch:
  feature/cgl-landau-fluid

pre-F116 local baseline:
  c4ddb25d574816f469c4fc61f756de5b9cf82d25

currently published bridge head:
  36140ea825cb853b298714c27720440fdab60b9e
```

F116 must replace the pre-release baseline with the exact final committed and
pushed release HEAD. Do not copy either hash above into a post-F116 production
packet without recomputing and verifying the release identity.

Production packets must pin the validated frozen E03 source tree rather than
the moving live checkout:

```text
/autofs/nccs-svm1_home2/dfielding/athenak-cgl-e03-9e075422
```

The retained bridge bundle for the currently published production history is:

```text
/lustre/orion/ast207/proj-shared/dfielding/CGL/source-archives/athenak-feature-cgl-through-36140ea82.bundle

SHA-256:
  2c2f57a166877387244dd5bb6bdf87beb12492ea075a7431939b78e5df7307a0
```

This bridge terminates at revision
`36140ea825cb853b298714c27720440fdab60b9e`. It remains historical bridge
evidence, not the final F116-selected current-source bundle. F116 must bind the
final release HEAD and its new complete-history bundle before F117 or the next
production wave is promoted.

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

The historical F-112 recost publication projected:

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

### 3.2 Current Campaign Snapshot

The first accelerated wave is drained, recorded, and independently reviewed.
It used four overlapping case lanes and ten concurrent Frontier nodes:

| Case | Current retained production state |
| --- | --- |
| `R02` | Fully complete through exact `t = 10`, bundled, and analyzed |
| `R03` | Exact `t = 0.5` accepted from job `4766828` |
| `R04` | Exact `t = 0.25` accepted from job `4766847` |
| `R12` | `s00_rankio_t0_t0p25`, job `4766856`, retained as `clean_partial` at exact `t = 0.1371931229426507`; inventory-only and explicitly not continuation-authorizing |
| `R16` | Exact `t = 1.5` accepted from job `4766866` |
| `R05-R11`, `R13-R15`, `R17` | Not started |

The next R12 production lineage is a fresh rerun:

```text
case:          R12
segment:       s01_rankio_t0_t0p12
start:         t = 0
target:        t = 0.12
parent:        null
restart:       null
nodes:         4
ranks:         32
Slurm:         02:00:00
Athena:        01:50:00
```

Do not use job `4766856`, its restart siblings, or its `clean_partial` result as
continuation authority. Preserve it as exact historical inventory and evidence.
That job reached `t = 0.1371931229426507` in `6657` seconds. The measured rate
projects the fresh `t = 0.12` target at `5822.7` seconds, making `0.12` the
largest natural `0.02`-aligned target below the `5940`-second 90% runtime
limit.

The canonical campaign root is currently drained:

```text
ledger rows:                  34
reservations:                 36 total; 34 recorded and 2 cancelled
active reservations:         0
manifests:                    36
transactions:                0
E03 cumulative node-hours:   36.577224
```

The next production action is gated on completed F116 and F117 publication, not
on additional interpretation of the first-wave jobs.

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

1. The inventory-only R12 job `4766856` remains non-authorizing, and the fresh
   `R12/s01_rankio_t0_t0p12` lineage starts from `t = 0` with null parent and
   restart.
2. Every mapped case `R03` through `R16` reaches exact `t = 10` through
   authenticated accepted segments.
3. `R17` reaches exact `t = 10` last, on eight nodes, through authenticated
   accepted segments.
4. Every mapped case `R02` through `R17` has one accepted whole-case bundle.
5. `bundle-campaign` succeeds for the frozen 16-case matrix.
6. The campaign-level paper analyzer succeeds and regenerates the panel-status
   table.
7. Every admitted comparison panel has preregistered quantitative criteria and
   is either `passed` or explicitly `blocked_out_of_scope`; `not_run`,
   `pending_review`, `failed`, and `inconclusive` do not satisfy final release.
8. Every case passes the preregistered steady-state, statistical-adequacy, and
   family-specific physics gates, or an inconclusive case is extended and
   reevaluated under the same criteria.
9. Retained terminal and late-time restart states pass the independent sampled
   normalized CT-`divB` audit, with the historical limitation disclosed.
10. The `R16`/`R02`/`R17` lane passes the preregistered resolution-convergence
    criteria.
11. Reference-blocked panels remain explicitly disclosed unless author data,
   archive data, donor diagnostics, or a qualified conversion unblock them.
12. The final ledger, reservation store, transactions, source archives,
   controller provenance, scheduler evidence, and analysis products reconcile.
13. The durable handoff is refreshed with the exact R17 completion boundary.

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

The first accelerated wave established that independent `R03` through `R16`
case lanes may overlap while controller mutations remain serialized. Continue
with rolling bounded-concurrency waves, preserve one in-flight segment per case
lineage, and use up to the reviewed ten-node wave envelope. Drain all
lower-resolution lanes before starting `R17`, which remains exclusive and last.

The production filesystem threat model is explicit:

- the canonical Frontier hierarchy intentionally traverses the trusted
  project-owned mode-`2770` `/lustre/orion/ast207/proj-shared` directory;
  production tools must not reject this deployment solely for group
  writability;
- tools must reject world-writable or otherwise untrusted authority profiles
  and bind the exact public-root, parent, lock, and target identities through
  every operation;
- tools defend against concurrent namespace, profile, link, and content
  mutation observable before or after a raw filesystem syscall;
- if authority is lost during an ambiguous syscall, tools durably classify the
  resulting state, fail closed, and perform no rollback or second namespace
  mutation;
- user-space code does not claim it can prevent a hostile actor from
  interposing inside the kernel syscall itself.

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
| Closing the current F116/F117 release gate | Final seven-tool validation, adversarial security review, source-bundle construction, F116 reviews, F117 recost/review, next-wave packet preparation |
| Preparing the next wave | R03/R04/R16 continuation packets, fresh R12 packet, queue and storage review, provisional later-wave composition |
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

## 7. Closed Recovery Decisions And Retained Hardening

### 7.1 Production Lineage Decisions

Keep Stage I on the already qualified corrected-E03 executable. The R03 recovery
decision is closed: job `4766828` reached exact `t = 0.5`, was accepted, and is
the authoritative R03 continuation boundary.

The R12 decision is also closed. Job `4766856` stopped cleanly at
`t = 0.1371931229426507`, but its retained evidence is inventory-only and does
not authorize continuation. The next R12 packet is
`s01_rankio_t0_t0p12`, fresh from `t = 0`, with null parent and restart. Its
four-node, 32-rank, `02:00:00` Slurm / `01:50:00` Athena profile targets
`t = 0.12`: the measured `6657` seconds to `t = 0.1371931229426507` projects
`5822.7` seconds to `t = 0.12`, the largest natural `0.02` target below the
`5940`-second 90% limit.

For the qualified E03 executable:

- continue exact accepted lineages only from inspector-authenticated complete
  restart sibling sets;
- accept ordinary exact-target segments through the existing strict path;
- retain non-authorizing partials as inventory without forcing continuation;
- reduce the next exact target immediately after any wall-clock partial;
- prefer a fresh rerun whenever retained evidence does not independently
  authorize the intended continuation.

### 7.2 Retained Historical Recovery Boundary

The binary-aware restart-header authenticator and checksum-bound historical
helper transition remain important retained protections. They must continue to:

1. parse the binary `Mesh::time` field from each selected restart sibling;
2. require one complete ranked sibling set and exact sibling agreement;
3. require terminal binary time to match final synchronized history;
4. reject unknown executable revisions, unknown ABI layouts, absent markers,
   sibling disagreements, nonfinite values, and arbitrary tolerance bypasses;
5. preserve historical submitted-manifest evidence without introducing a
   general manifest bypass.

These protections do not turn the R12 `clean_partial` into continuation
authority. Fresh R12 execution is the selected production strategy.

### 7.3 Writer Precision Repair

Retain the scoped restart-marker precision repair local to restart metadata:

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
The binary-aware controller path remains retained historical protection; the
writer repair is qualification-ready code for a deliberate future executable
transition.

### 7.4 Required Retained Hardening Tests

Retain and rerun as part of the release suite:

1. A controller fixture with a non-round final time proving:
   - a binary-exact sibling set with legacy text marker `0.312823` is accepted
     only by the binary-aware path;
   - full-precision text and binary agreement is accepted;
   - a text marker inconsistent with legacy serialization is rejected;
   - binary sibling disagreement is rejected;
   - terminal-history disagreement is rejected;
   - an unknown ABI or executable revision is rejected;
   - `record --result clean_partial` succeeds only after valid inspection.
2. A historical submitted-manifest fixture proving the checksum-bound helper
   transition permits only its reviewed historical shape and rejects a broad
   bypass.
3. A CGL-LF restart regression that parses emitted `time/restart_time` and
   requires full-precision round-trip agreement with terminal history.
4. The modal forcing restart regression across an OU refresh.
5. Focused CPU, MPI CPU, turbulence-driver CPU, and turbulence-driver MPI CPU
   suites.
6. A compact Frontier one-node, eight-rank restart qualification when an
   executable transition requires it:
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

## 8. Current Release Gate And Next Wave

### 8.1 Complete F116 And F117

The campaign is drained at the first-wave barrier. Before the next production
submission:

1. finish the final seven-tool implementation, focused adversarial tests, full
   release suite, and independent security and plasma review;
2. commit and push the final release revision;
3. construct and verify its complete-history source bundle while retaining the
   `36140ea82` bridge bundle as historical evidence;
4. publish and verify F116 current-source authority;
5. publish and verify F117 recost and next-wave recommendation;
6. reconcile again with zero active reservations and zero transactions.

F116 and F117 are release gates, not completed work. F116 changes current source
selection only. F117 must preserve job `4766856` as non-authorizing inventory
and recommend the fresh R12 rerun.

F116 and F117 authenticate exactly the retained seven-tool production-control
vector. The standalone scientific-acceptance utility and criteria may be
present in the complete-history source bundle, but they are non-authorizing
until their separate independent plasma/statistical and restart-format reviews
are approved and bind their exact final digests.

### 8.2 Next Ten-Node Wave

Subject to promoted F116/F117 authority and ordinary preflight, the intended
next rolling wave is:

| Case | Segment | Start | Target | Nodes | Parent/restart |
| --- | --- | ---: | ---: | ---: | --- |
| `R03` | `s03_rankio_t0p5_t0p75` | `0.5` | `0.75` | `1` | accepted job `4766828` lineage |
| `R04` | `s02_rankio_t0p25_t1p25` | `0.25` | `1.25` | `4` | accepted job `4766847` lineage |
| `R12` | `s01_rankio_t0_t0p12` | `0` | `0.12` | `4` | null parent; null restart |
| `R16` | `s01_rankio_t1p5_t4p5` | `1.5` | `4.5` | `1` | accepted job `4766866` lineage |

This preserves the demonstrated ten-node envelope while allowing all four
allocations to overlap. Issue `prepare`, `check-submit`, and `submit` serially;
the Frontier allocations may run concurrently.

### 8.3 Fresh R12 Rule

The fresh R12 rerun is not a continuation workaround or waiver:

- use segment `s01_rankio_t0_t0p12`;
- start at `t = 0`;
- target exact `t = 0.12`;
- set parent and restart to null;
- use four nodes, 32 ranks, Slurm `02:00:00`, and Athena `01:50:00`;
- keep later R12 continuation increments on the natural `0.02` cadence,
  subject to authoritative recosting; all other cases retain the quarter-unit
  increment rule;
- preserve `R12/s00` job `4766856` unchanged as inventory-only evidence;
- prohibit any later planner, recost, or controller path from treating job
  `4766856` as continuation authority.

## 9. Concurrent Production Cadence For R03 Through R16

### 9.1 Bounded-Concurrency Transition

The first accelerated wave demonstrated the bounded-concurrency policy with
R03, R04, R12, and R16 allocations overlapping across ten nodes. The next
release must retain the reviewed controller behavior that:

1. permits multiple active reservations only when they belong to distinct
   `R03` through `R16` case lanes;
2. preserves at most one prepared or submitted segment per case lineage;
3. serializes every metadata mutation under the existing canonical root lock;
4. accounts the sum of all active reservation node-hours against the reviewed
   Stage I envelope and the `4000` node-hour project ceiling;
5. permits queued or running Stage I jobs only when each one matches an
   authenticated submitted reservation and exact retained manifest;
6. continues to reject unbound, malformed, forged, or duplicate Stage I CGL
   queue rows, unreviewed shared-root campaign records, orphaned run
   directories, pending transactions, and ambiguous scheduler state;
7. allows reviewed multi-node allocation profiles for `R04` through `R16`;
8. keeps `R17` fixed at eight nodes, requires accepted `t = 10` predecessors,
   and permanently locks out lower-resolution preparation after R17 starts;
9. replays historical one-node records without rewriting their provenance.
10. validates current, prior, and payload reservation snapshots against the
    promoted Stage I envelope and the `4000` node-hour project ceiling before
    any ledger append, reservation-store rewrite, or manifest rewrite, while
    permitting the durable recovery journal and counting a recorded journal
    row exactly once during replay;
11. retains authenticated initial and final queue snapshots, with the final
    query immediately before the durable ambiguity barrier and real `sbatch`;
12. rechecks the R17-last invariant during every lifecycle replay and
    reconciliation.

Complete focused controller fixtures proving:

- two distinct cases may be prepared, submitted, inspected, and recorded in
  either completion order;
- a second active segment for the same case is rejected;
- summed reservations and the lane cap fail closed;
- an exact authenticated Stage I queue set and well-formed unrelated non-CGL
  jobs are accepted while an unbound Stage I CGL job is rejected;
- a queue change between preflight and real `sbatch` is rejected before the
  ambiguity barrier;
- prepared and recorded replay budget overruns fail before any ledger append,
  reservation-store rewrite, or manifest rewrite;
- R17 submitted, recorded, and reconciliation paths fail closed without
  complete exact-`t = 10` predecessors;
- a pending transaction or orphaned run directory still blocks mutation;
- approved `1`, `2`, and `4` node profiles replay correctly where authorized,
  while unapproved or decomposition-infeasible profiles are rejected;
- R17 remains exclusive and last.

The first wave supplies production evidence for authenticated submitted-job
overlap and durable preservation of unrelated submitted lanes while completed
lanes are recorded. The F116/F117 release gate must retain focused coverage for
node-profile authorization, duplicate-case rejection, the four-lane cap,
one-prepared-packet enforcement, unbound-CGL-job rejection, fresh-R12 lineage
selection, and R17 exclusivity before the next concurrent wave.

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

- keep R03 on one node under its measured accepted profile;
- use the reviewed four-node standard-layout profile for R04 and the fresh R12
  rerun, then apply measured family evidence to later standard cases;
- use the reviewed one-node profile for lower-resolution R16;
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

### 9.4 Aggressive Profiles

The first-wave rows report retained results; the unstarted rows remain initial
proposals, not permanent limits. Recost provisionally after each accepted
prefix and authoritatively at wave barriers.

| Case | Runtime family | Current or initial exact increment | Reviewed node profile | Reason |
| --- | --- | ---: | --- | --- |
| `R03` | active Alfvenic beta-100 hard wall | Accepted exact `0.5`; next target `0.75` | `1` | Accepted job `4766828` establishes the current lineage |
| `R04` | active random beta-10 | Accepted exact `0.25`; next target `1.25` | `4` | Accepted job `4766847` establishes the four-node standard profile |
| `R05` | active random beta-100 | `0.25` | Reuse reviewed standard winner | Beta-100 and random forcing both merit measured calibration |
| `R06` | passive Alfvenic beta-10 | `0.50` | Reuse reviewed standard winner | Closest passive counterpart to completed R02 |
| `R07` | passive Alfvenic beta-100 | `0.25` | Reuse reviewed standard winner | Reuse R03 beta-100 bound until measured faster |
| `R08` | passive random beta-10 | `0.25` | Reuse reviewed standard winner | Reuse R04 random-family calibration, then ratchet |
| `R09` | passive random beta-100 | `0.25` | Reuse reviewed standard winner | Reuse beta-100 random-family calibration |
| `R10` | compressive active random beta-1 | `0.25` | Recheck reviewed standard winner | New compressive family |
| `R11` | compressive active random beta-100 sonic-correlation | `0.25` | Recheck reviewed standard winner | New beta-100 sonic-correlation family |
| `R12` | stronger LF heat flux | Fresh `0 -> 0.12`; then `0.02`-aligned continuation increments; prior partial is inventory-only | `4` | Start `s01_rankio_t0_t0p12` with null parent/restart |
| `R13` | weaker LF heat flux | `0.25` | Recheck reviewed standard winner | Distinct closure-cost profile |
| `R14` | beta-100 finite limiter `nu_lim = 20` | `0.25` | Reuse reviewed standard winner | New finite-limiter profile |
| `R15` | beta-100 finite limiter `nu_lim = 200` | `0.25` | Reuse reviewed standard winner | New finite-limiter profile |
| `R16` | beta-10 `96 x 96 x 192` scale separation | Accepted exact `1.5`; next target `4.5` | `1` | Accepted job `4766866` supports a larger next increment |

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

The initial four-lane R03/R04/R12/R16 wave is complete and drained. It used ten
concurrent nodes and established the accepted R03, R04, and R16 prefixes plus
the inventory-only R12 partial. Continue as follows:

1. Complete and verify the F116/F117 release gate.
2. Launch the Section 8.2 ten-node wave, including fresh R12 from `t = 0`.
3. Keep the controller cap at four lanes and the reviewed wave envelope at ten
   nodes. Any increase requires a separate reviewed transition.
4. Refill a lane immediately after its accepted segment is recorded if the
   provisional model, storage monitor, and reservation sum remain healthy.
5. Fill subsequent unstarted lanes in this order:
   `R06 -> R10 -> R11 -> R13 -> R05 -> R07 -> R14 -> R15 -> R08 -> R09`.
6. Drain periodically for authoritative recost barriers and always before a
   lane-cap or allocation-profile increase.
7. Drain every lower-resolution lane, publish the final predecessor recost,
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
4. confirm the required retained cadence for the `t = 4` through `10`
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

The drained canonical campaign currently records `36.577224` cumulative E03
node-hours. The first accelerated wave contributes measured R03, R04, R12, and
R16 throughput for the F117 recost.

The F-112 projection had limited margin inside the historical `900` node-hour
Stage I reservation. The measured R03 rate gives the conservative refreshed
transition arithmetic:

```text
measured conservative projection: 1217.254885 node-hours
promoted Stage I envelope:         1400.000000 node-hours
remaining planning headroom:        182.745115 node-hours
incremental project ceiling:       4000.000000 node-hours
```

Retain the reviewed `1400` node-hour Stage I envelope unless the F117 measured
recost requires a reviewed change. Increase it again only when measured scaling
or family timing requires another reviewed transition.

### 11.2 Reservation Transition

Before the next concurrent wave:

1. have F117 consume the accepted first-wave results and the inventory-only R12
   partial without granting it continuation authority;
2. include projected multi-node node-hours plus the summed maximum active-wave
   reservation;
3. retain `required_storage_safety_bytes = 1099511627776` (exactly `1 TiB`),
   which exceeds the reviewed `65,998,006,704`-byte F117 wave projection and
   preserves the separately required R17 storage margin;
4. have an independent reviewer audit arithmetic, contingency, and the fresh
   R12 profile;
5. promote F116 current source authority before F117;
6. verify both publications and reconcile the canonical root;
7. resume under the retained `4000` node-hour project ceiling.

This is an accounting-control update, not a reason to reopen qualified
physics.

### 11.3 Recost Frequency

Update the provisional throughput and node-efficiency model:

- after every accepted or inventory-only measured packet;
- after the first accepted prefix of each new runtime family;
- after any early wall-clock termination;
- after every completed case.

Publish an authoritative recost under the retained zero-active-reservation
companion contract:

- at the current F117 first-wave barrier;
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
- steady-window selection over retained `t = 4` through `10`, with
  stationarity comparison between `t = 4` through `7` and `t = 7` through
  `10`;
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

### 12.4 Preregistered Scientific Acceptance

Use the standalone Stage I scientific-acceptance utility and its independently
reviewed criteria artifact. Do not make final scientific acceptance depend on
uncommitted analyzer work. Evaluate each completed case concurrently from
immutable whole-case products over:

```text
full window:   t = 4 through 10
early window:  t = 4 through 7
late window:   t = 7 through 10
```

The preregistered default gates are:

- endpoint-clipped trapezoidal time-weighted means;
- full-window independent forcing-correlation blocks `>= 3` and half-window
  blocks `>= 1.5`; never claim more independent blocks than the physical
  window duration divided by the authenticated forcing correlation time;
- deterministic moving-block-bootstrap uncertainties;
- reject or mark inconclusive histories with physical-time gaps larger than
  the preregistered cadence/forcing-correlation threshold;
- early/late scalar drift `z <= 3` and relative change `<= 25%`;
- occupancy drift `z <= 3` and absolute change `<= 0.002`;
- forcing-power relative drift `<= 10%`;
- reference-panel normalized-residual RMS `<= 2` and maximum absolute
  normalized residual `<= 5`;
- generic meaningful accumulated activity `> 1e-6` and normalized activity
  `> 1e-8`.

An undersampled but otherwise healthy case is `inconclusive`, not failed.
Extend only that case under the retained extension policy; do not relax the
criteria after results are visible.

Require the following family gates:

- finite-limiter `R14`/`R15`: exact zero hard-wall projections, positive
  late-time effective collisionality, threshold occupancy in both halves, and
  statistically resolved `nu_eff(R15) > nu_eff(R14)`;
- hard-wall cases: hard-wall activity in both halves, nonzero threshold
  occupancy, and exact zero hard-bound volume;
- LF-strength `R12`/`R02`/`R06`/`R13`: nonzero LF face activity, meaningful
  LF work, and statistically resolved retained-response differences;
- active/passive pairs `R02/R06`, `R03/R07`, `R04/R08`, and `R05/R09`:
  exact-zero passive CGL pressure-work terms, meaningful active work, and at
  least one Holm-corrected standardized contrast `>= 0.5`;
- forcing families: Alfvenic parallel-forcing fraction `<= 1e-10` and random
  parallel-forcing fraction lower 95% bound `> 0.05`.

Audit sampled normalized CT-`divB` from authenticated native restart face
fields at every retained terminal and late-time state. Require:

```text
abs(divB) * min(dx) / max(abs(B), bfloor) < 1e-12
```

Qualify the independent restart parser against a diagnostic-enabled test. The
historical E03 claim is limited to sampled retained states because those runs
did not retain a full-time-history `max_ndiv` diagnostic.

Before inspecting R17 results, retain the resolution-convergence interval:

```text
4 <= k_perp / pi <= 24
alignment shells = [4, 6, 8, 12, 16, 24]
```

Require `R02`/`R17` alignment difference `<= 0.05`, spectral log-RMS
difference `<= 0.15`, and global scalar differences `<= 10%` or within two
combined standard errors. Where `R16`/`R02` disagreement is resolved, require:

```text
distance(R02, R17) <= 0.75 * distance(R16, R02)
```

## 13. Case Bundling And Analysis

### 13.1 Bundle Each Completed Case Immediately

After one case reaches accepted exact `t = 10`:

1. run `bundle-case`;
2. authenticate the lineage and merged histories;
3. run the paper analyzer once for that whole-case bundle;
4. retain diagnostics and generated figures;
5. review the exact `t = 4` through `10` steady window and its `t = 4..7`
   versus `t = 7..10` stationarity split;
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
SOURCE_BUNDLE="$ROOT/source-archives/athenak-feature-cgl-through-36140ea82.bundle"

# Replace the historical bridge SOURCE_BUNDLE only after F116 publishes and
# verifies the final complete-history bundle.

python3 "$STAGE" --root "$ROOT" reconcile

# Fresh R12 rerun; no parent or restart-file is permitted.
python3 "$STAGE" --root "$ROOT" prepare \
  --case-id R12 \
  --segment s01_rankio_t0_t0p12 \
  --acceptance-criterion "<reviewed exact-boundary acceptance prose>" \
  --executable "<qualified executable>" \
  --build-manifest "<qualified build manifest>" \
  --source-dir "$FROZEN_SOURCE" \
  --source-bundle "<F116-selected final complete-history bundle>" \
  --matrix "$MATRIX" \
  --nodes 4 \
  --walltime 02:00:00 \
  --athena-walltime 01:50:00 \
  --override time/tlim=0.12

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

Do not add a parent or `--restart-file` to the fresh R12 packet. Job `4766856`
remains retained inventory only and is not a continuation source.

After F116 and F117 are promoted and verified, issue `prepare`, `check-submit`,
and `submit` serially for each selected `R03` through `R16` lane, but allow
authenticated submitted allocations to overlap on Frontier.
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

1. F116 current-source-authority publication and final complete-history bundle;
2. F117 first-wave recost and next-wave recommendation, including fresh R12;
3. every later drained-wave authoritative recost and lane-cap change;
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

- [x] Complete, bundle, and analyze R02 through exact `t = 10`.
- [x] Run and drain the first accelerated four-lane, ten-node wave.
- [x] Accept R03 job `4766828` at exact `t = 0.5`.
- [x] Accept R04 job `4766847` at exact `t = 0.25`.
- [x] Retain R12 job `4766856` as inventory-only `clean_partial` at
      `t = 0.1371931229426507` and prohibit continuation from it.
- [x] Accept R16 job `4766866` at exact `t = 1.5`.
- [x] Reconcile the first-wave barrier to zero active reservations and zero
      transactions.
- [x] Implement the scoped restart-marker precision patch in parallel.
- [x] Implement and run the first focused regression tranche for the
      bounded-concurrency controller transition.
- [ ] Finish the final seven-tool release implementation, validation, and
      independent audits.
- [ ] Commit, push, bundle, independently review, publish, and verify F116.
- [ ] Generate, independently review, publish, and verify F117.
- [ ] Prepare and launch the next ten-node wave, including fresh
      `R12/s01_rankio_t0_t0p12` with null parent/restart.
- [ ] Request raw curves, FFT-normalization details, or donor diagnostics from
      Stephen Majeski in parallel.

### Standard Matrix

- [x] Launch and drain the initial four-lane `R03`, `R04`, `R12`, `R16` wave.
- [ ] Promote the first drained-wave F117 recost, storage audit, and
      next-wave profile table.
- [ ] Keep the rolling campaign at the reviewed four-lane cap.
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
- [ ] Validate and independently review the preregistered scientific-acceptance
      criteria.
- [ ] Pass steady-state, statistical-adequacy, and family-specific case gates.
- [ ] Pass sampled retained-state normalized CT-`divB` audit.
- [ ] Pass the `R16`/`R02`/`R17` resolution-convergence gate.
- [ ] Require every admitted comparison panel to pass or be explicitly
      blocked out of scope.
- [ ] Disclose or resolve reference-blocked panels.
- [ ] Reconcile the final ledger, archives, and provenance.
- [ ] Refresh the durable handoff at accepted R17 completion.

## 17. Completion Decision

The campaign should proceed aggressively. The code has broad retained tests,
corrected forcing-policy qualification, one complete and analyzed
corrected-E03 production case, a successful ten-node accelerated wave, a
controller that fails closed, and encouraging independent physics evidence
from Stephen Majeski. The shortest robust path is:

1. finish and publish F116 current source authority, then F117 first-wave
   recost and next-wave recommendation;
2. launch the next four-lane, ten-node wave with R03, R04, fresh R12, and R16;
3. run R03 through R16 as rolling bounded-concurrency waves while serializing
   shared-root mutations and preserving one in-flight segment per case;
4. use drained barriers for authoritative recost, storage review, and
   allocation review while analysis continues in parallel;
5. drain every lower-resolution lane and execute R17 exclusively and last with
   eight-node, recosted exact increments;
6. assemble and review the final campaign bundle.
