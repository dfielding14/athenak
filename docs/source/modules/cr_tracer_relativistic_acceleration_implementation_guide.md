---
orphan: true
---

# CR Tracer Relativistic Acceleration Implementation Guide

## Document Status

| Item | Value |
| --- | --- |
| Purpose | Implementation, review, validation, and handoff guide for passive relativistic cosmic-ray tracer acceleration and deceleration |
| Working branch | `feature/CR_tracers_relativistic_acceleration` |
| Required prerequisite branch | `feature/CR_tracers_followup_architecture` |
| Prerequisite tip when this guide was drafted | `64a4d1be` |
| Initial physical scope | Passive test particles in a Newtonian ideal-MHD background |
| Initial coupling scope | One-way MHD-to-particle field sampling only; no particle feedback onto MHD |
| Relationship to the existing CR tracer plans | Extends the existing fixed-energy CR tracer architecture and accuracy plans without weakening their regression contract |
| Relationship to the local `PIC` branch | Bounded source-level prior art only; do not merge or transplant the coupled PIC implementation wholesale |
| Document role | Living implementation gate. Update it when evidence changes the plan. Do not treat it as a static checklist to rush through. |

This guide is deliberately more procedural than the earlier CR tracer plans.
The added physics changes the particle state, the pusher, the restart contract,
the timestep logic, and the meaning of several diagnostics. A patch can compile
and still be physically wrong. Agents must therefore work in small, reviewable
increments, record design choices, run independent reviews repeatedly, and
pause at the reflection gates below before widening the implementation.

## Control Header

Update this table at every accepted checkpoint. Do not leave it stale while
implementation proceeds.

| Control | Current value at guide creation | Update rule |
| --- | --- | --- |
| Frozen feature base | `64a4d1be` | Change only through an accepted integration-target decision |
| Current implementation phase | `Phase 0: freeze baseline and open the ledger` | Advance only after the current reflection gate records `PROCEED` |
| Intended eventual integration target | Unresolved | Set in `DR-000` before production code edits |
| Allowed write manifest | Guide and ledger only until `CP-0 Baseline` records `PROCEED` | Open only the next phase's edit set; each dependent decision through `DR-018` must be accepted first |
| Last accepted checkpoint | Guide creation only | Record checkpoint ID, commit, reviewers, verdict, and evidence paths |
| Open blocking findings | None recorded yet | Keep a live list; do not proceed while any blocking finding remains |
| GPU qualification status | Not started | Do not claim GPU readiness until run on the accepted SHA |

### Qualification Verdicts

Track these separately:

- `CPU/MPI QUALIFIED`: CPU, MPI, restart, diagnostics, AMR, and documentation
  gates pass on the accepted SHA.
- `GPU QUALIFIED`: required accelerator gates pass on the same accepted SHA.
- `MERGE READY`: `CPU/MPI QUALIFIED` and `CP-4 Solver Coupling` are mandatory;
  accepted `DR-000` states whether `GPU QUALIFIED` is additionally required;
  all required verdicts pass; and residual risk is recorded.

Unavailable backends remain explicit residual risk. They are not silently
converted into passing evidence. Record required backend verdicts in accepted
`DR-000` before `CP-7`. An unavailable required backend forces `HOLD`; an
intentionally optional backend remains explicit residual risk.

### Integration-Target Risk Must Be Resolved Early

A read-only `git merge-tree` comparison against the current
`origin/development` already reports core integration overlaps. At guide
creation, content conflicts were reported in:

- `src/mesh/mesh.cpp`;
- `src/mesh/mesh.hpp`;
- `src/outputs/outputs.cpp`;
- `src/outputs/outputs.hpp`;
- `src/pgen/pgen.cpp`.

Additional overlapping edits that appeared auto-mergeable at guide creation
included:

- `src/CMakeLists.txt`;
- `src/mesh/mesh_refinement.cpp`;
- `src/outputs/basetype_output.cpp`;
- `src/pgen/pgen.hpp`.

This does not mean the relativistic feature should be rebased immediately. It
means implementation agents must record the intended target branch, expected
stacking order, conflict strategy, and retest obligations before the branch
grows. Rerun the merge-tree audit after target-branch updates and before final
handoff.

## Executive Recommendation

The first implementation should add an **opt-in passive relativistic
ideal-MHD CR tracer mode** while preserving the current fixed-energy magnetic
tracer mode as a compatibility control.

The new mode should:

1. First validate the relativistic kernel against prescribed analytical
   electric and magnetic fields in a narrowly scoped test harness.
2. After prescribed-field acceptance, gather the fluid velocity and magnetic
   field at the particle location using a documented common interpolation
   policy.
3. Construct the particle-sampled pointwise ideal-MHD field from those gathered
   quantities, using an explicitly chosen convention:

   ```text
   cE_particle = -u_fluid,particle x B_particle
   ```

   Do not use `E` and `cE` interchangeably. The exact normalization, including
   whether the code variable represents physical `E`, `cE`, or another
   normalized quantity, must be settled and documented before implementation.
4. Advance a relativistic momentum-like particle state with an explicit
   electromagnetic pusher whose kernel order and solver-coupled temporal order
   are qualified separately.
5. Derive the Lorentz factor and physical drift velocity consistently from the
   authoritative momentum state.
6. Keep the implementation passive: the new mode samples MHD fields but does
   not deposit force, current, energy, or momentum back onto the gas.
7. Preserve old fixed-energy tests and add a separate relativistic validation
   ladder rather than rewriting existing tests until they pass.

### Physics Correction: Do Not Curl The Particle Force

The particle acceleration does **not** require taking a curl of the
interpolated velocity or magnetic field. Write the adopted normalized
relativistic Lorentz equations explicitly before coding. For example, if the
stored electromagnetic quantity is `cE`, the exact physical-units form is:

```text
cE_particle = -u_fluid,particle x B_particle
d(p/m)/dt = (q / (m c)) [cE_particle + v(p) x B_particle]
dK/dt = (q / c) cE_particle dot v(p)
```

The adopted code-unit convention may absorb factors of physical `c`, but the
decision record must show where each factor went. Replacing physical `c` with a
configurable reduced `C` is a reduced-light-speed model choice, not a units
conversion. `DR-001` must state whether and where that replacement enters the
state relation, Lorentz-force prefactor, sampled-field mapping, and work
diagnostic. Do not copy the example blindly.

The curl appears in the MHD induction equation,

```text
dB/dt = -curl(E)
```

and AthenaK's constrained-transport MHD path already owns that operation. Do
not take a curl inside the particle pusher unless a later, separately reviewed
physical model actually requires a spatial derivative.

### Recommended State Model

Use a momentum-like state as the authoritative evolving CR quantity and derive
the Lorentz factor from it. A useful candidate is the mass-normalized momentum

```text
u_cr = p / m = gamma * v
gamma = sqrt(1 + |u_cr|^2 / c_model^2)
v = u_cr / gamma
```

where `c_model` is physical `c` unless an accepted reduced-light-speed model
replaces it with configurable `C`.
The exact variable names, units, restart representation, and whether `gamma`
is cached or always derived are mandatory pre-code decisions. Storing only a
Lorentz factor is insufficient because the force update needs directional
momentum.

### Required Relativistic Work-Energy Check

For accelerating cases, fixed-speed assertions are no longer appropriate.
The new validation ladder must measure the relativistic work-energy identity
in the selected normalization:

```text
Delta(gamma * m * c^2) = integral(q * E dot v dt)
```

or, if the adopted stored quantity is `cE`,

```text
Delta(gamma * m * c^2) = integral((q / c) * cE dot v dt)
```

Record the discrete form used by the pusher, output an explicit quadrature of
the sampled field work needed to audit it, and preregister a tolerance before
inspecting results. Never compute accumulated work from `Delta K` itself. If a
reduced-light-speed model is accepted, derive and document its corresponding
work identity rather than replacing symbols mechanically.

## Scope Boundary

### Required In The First Relativistic Feature

- An explicit opt-in runtime mode for passive relativistic CR tracers.
- A documented particle state and units contract.
- A documented charge-to-mass or mass-to-charge convention.
- Interpolation of fluid velocity alongside magnetic-field interpolation.
- Construction and optional output of the particle-sampled ideal-MHD electric
  field.
- A reviewed relativistic electromagnetic pusher.
- Dynamic Lorentz factor, derived velocity, acceleration, and deceleration.
- Subcycle constraints that remain valid while particle momentum changes.
- New-version restart write/read support and an intentional legacy-restart
  policy.
- Inspector support for every restart version that remains supported.
- Diagnostics with unambiguous names and units.
- Dedicated analytical, restart, MPI, AMR, and compatibility tests.
- Documentation of every user-facing input and output semantic change.

### Explicitly Out Of Scope For The First Relativistic Feature

- Particle feedback onto MHD momentum or energy.
- CR current deposition.
- CR Hall corrections to the CT electric field.
- Resistive, non-ideal, or externally imposed electric-field models unless
  introduced as a narrowly scoped test harness.
- Delta-f particle methods.
- Expanding-box transformations.
- Full MHD-PIC reproduction claims.
- SRMHD or GRMHD field transformations.
- Performance optimization that obscures the reference implementation before
  the analytical tests pass.

If an agent believes one of these exclusions must be crossed, stop and create
a decision record before changing code. Do not widen the patch opportunistically.

## Mandatory Working Protocol

This section is binding for implementation work on this branch.

### Rule 1: Establish A Baseline Before Editing

Before the first implementation patch:

1. Record the exact base commit and dirty-tree state.
2. Run the existing CR tracer particle suite on the unmodified branch.
3. Archive command lines, configuration, and results.
4. Confirm that the current fixed-energy CR tracer mode still has a clean,
   named runtime path that can remain available as a control.
5. Read the existing feature, follow-up architecture, and accuracy plans.
6. Create the initial decision ledger entries listed below, even if some remain
   `proposed`.
7. Run a read-only merge-tree comparison against the intended integration
   target and record every overlap.
8. Reconcile the historical CR tracer plans into the baseline ledger below.

Do not begin implementation by changing the particle data layout. First freeze
the existing behavior and document the migration target.

### Rule 2: Use Small Commits With One Claim Each

Each commit should make one reviewable claim. Good examples:

- add a unit-tested relativistic state helper;
- add runtime parsing for an opt-in mode without changing push behavior;
- add fluid-velocity gather and gather-only tests;
- add the pusher behind the new mode;
- add restart versioning and compatibility tests;
- add diagnostics and their analysis checks.

Avoid commits that combine state migration, pusher changes, restart rewrites,
output renames, test loosening, and performance refactors. A reviewer must be
able to identify the single intended behavior change in each diff.

### Rule 3: Record Every Nontrivial Choice

When two or more credible implementation choices exist, write a decision
record before or during the experiment. Do not silently choose the quickest
path. The record must identify:

- the question;
- constraints;
- considered options;
- selected option;
- evidence supporting the choice;
- rejected alternatives and why they were not chosen;
- migration and compatibility implications;
- tests that would falsify the decision;
- a revisit trigger.

The ledger is not bureaucratic overhead. It is how later agents can distinguish
deliberate design from accidental behavior.

At each checkpoint, list every choice made since the previous checkpoint,
including retained defaults and inherited behavior. Reviewers must identify
plausible unlogged alternatives. Record `none` explicitly if no new choice was
made.

### Rule 4: Spawn Independent Review Agents At Every Done Claim

Every time an implementation agent thinks a phase, milestone, or final patch is
done, it must spawn independent sub-agents to challenge that conclusion before
moving on. At minimum:

1. One **physics and numerics reviewer** checks equations, units, pusher
   algebra, limiting cases, and acceptance criteria.
2. One **integration and validation reviewer** checks code surfaces, restart
   migration, MPI/AMR behavior, diagnostics, tests, and documentation.

For restart, output, GPU, or performance milestones, spawn an additional
specialist reviewer for that area. Provide reviewers the commit hash, the
decision ledger, the exact commands run, and the evidence artifacts. Ask for
skeptical findings, not approval.

A done claim includes any commit presented as correct, any request to open the
next edit set, any checkpoint verdict, and any handoff or merge recommendation.
Spawn reviewers before further feature edits. Phase-boundary reviews broaden
the review surface; they do not replace per-claim review. Phase-specific
specialists are additional to the two-reviewer minimum. A substitution may
change reviewer expertise but must not reduce the two-independent-reviewer
minimum; record its rationale and reviewer concurrence.

An accepted done-claim review is a named checkpoint log containing reviewers,
prompts, evidence paths, findings, dispositions, rechecks, and an explicit
`PROCEED` verdict. A message saying only that review is complete does not open
the next edit set.

After review:

1. Record each finding.
2. Fix or explicitly defer each nonblocking finding with reviewer concurrence
   and a named residual-risk or follow-up entry.
3. Rerun affected tests.
4. Ask the independent reviewer to recheck the revised state.
5. Do not declare the milestone complete until blocking findings are closed.

The reviewer assigns initial severity. A blocking finding cannot be deferred
within the current scope. Reclassification requires reviewer concurrence and
an accepted decision record. A blocking finding that is moved out of the branch
requires an accepted split-scope decision and a `HOLD` or `STOP` verdict.

Use fresh reviewers for the final cold review. Do not substitute self-review,
an implementation agent's own summary, or a review timeout for independent
evidence. If sub-agent tooling is unavailable, record `HOLD` and do not open
the next edit set. Equivalent review means a documented review by a separate
agent or person with no implementation ownership, including prompt, findings,
dispositions, and recheck.

### Rule 5: Do Not Convert Failures Into Passing Tests By Relaxation

Never respond to a failed analytical or regression test by loosening its
tolerance, renaming it as a smoke test, or removing an assertion without a
written decision record. Determine whether the failure is caused by:

- an implementation defect;
- an invalid oracle;
- a unit mismatch;
- a temporal-centering mismatch;
- a gather-policy mismatch;
- an unsupported regime;
- floating-point limitations;
- a test setup error.

Change a threshold only after documenting the measured convergence behavior and
the reason the original threshold was wrong.

### Rule 6: Pause At Reflection Gates

The reflection gates below are mandatory. At each gate, stop adding features,
review accumulated evidence, ask whether the model and plan still make sense,
and decide whether to proceed, revise the plan, split the branch, or stop.
Every reflection record must state what changed in the plan, ledger, control
header, checkpoint sequence, or write manifest. If nothing changed, record why
the reviewed evidence and adversarial findings did not justify a change.

### Rule 7: Preserve Compatibility Controls

Do not silently change the meaning of the existing `boris` mode, existing
particle outputs, or legacy restart records. New semantics need explicit mode
names, output names, restart metadata, or an intentional migration decision.

### Rule 8: Prefer A Clear Reference Path Before Optimization

Implement and validate a readable device-capable reference path first. Optimize
only after the analytical suite passes and the profiler identifies a real
bottleneck. Every optimization needs before/after tests and a differential
check against the reference path.

### Rule 9: Use Named Checkpoints With A Verdict

Every phase checkpoint must end with:

- assumptions falsified;
- scope drift found;
- unresolved risks;
- evidence added;
- independent reviewer verdicts;
- an explicit `PROCEED`, `HOLD`, or `STOP` verdict;
- the next allowed edit set.

An agent may not treat a long checklist with unknown items as a passing
checkpoint.

Any missing artifact, unknown result, failed required test, post-hoc threshold,
or unresolved blocking finding forces `HOLD`. Link every minimum-evidence item
to a concrete artifact and preregistered acceptance criterion.

## Global Stop Conditions

Stop the current phase immediately, record the reason, and issue a `HOLD` or
`STOP` checkpoint verdict when any of these occurs:

- unexpected dirty-tree changes;
- a changed base SHA or target branch without an accepted `DR-000` update;
- an unresolved decision record that affects the next edit;
- edits outside the approved phase write manifest;
- a legacy behavior change without an accepted migration decision;
- restart incompatibility or ambiguous conversion;
- diagnostics with undefined relativistic meaning;
- a missing analytical oracle or unregistered tolerance;
- a pusher result that passes smoke tests but fails convergence or work closure;
- an unexplained CPU/MPI/AMR/GPU disagreement;
- an unresolved independent-review objection;
- a request to claim accelerator readiness without accelerator evidence;
- a request to widen into feedback, SRMHD, GRMHD, resistivity, orbital
  advection, or shearing-box behavior without reopening scope.

Do not bury a stop condition in a later TODO list. Resolve it, revise the plan,
split the branch, or document why the branch is blocked.

## Required Record Templates

### Decision Record Template

Create a ledger near this guide or in a clearly linked implementation handoff
document. Use one entry for each decision.

```md
## DR-XXX: Short decision title

- Status: proposed | accepted | superseded | rejected
- Date:
- Author:
- Commit:
- Question:
- Constraints:
- Options considered:
  1. Option A:
  2. Option B:
  3. Option C:
- Selected option:
- Evidence and reasoning:
- Why the other options were not selected:
- Compatibility and migration impact:
- Validation required:
- Failure signals:
- Revisit trigger:
- Independent reviewers:
- Follow-up actions:
```

### Reflection Record Template

Write a reflection record at every reflection gate.

```md
## RG-XXX: Short reflection title

- Date:
- Commit range reviewed:
- Evidence reviewed:
- Assumptions that still hold:
- Assumptions weakened or falsified:
- New risks:
- Blocking findings:
- Decision: proceed | revise plan | split scope | stop
- Plan changes:
- Tests added or changed:
- Independent reviewer findings and responses:
- Next smallest reviewable increment:
```

### Test Evidence Template

```md
## Evidence: short test name

- Commit:
- Build configuration:
- Runtime command:
- Input checksum:
- Analysis command:
- Expected result:
- Measured result:
- Acceptance threshold:
- Pass or fail:
- Artifact paths:
- Notes:
```

### Baseline Reconciliation Ledger Template

Classify relevant historical-plan items before implementation. The older plans
remain useful evidence, but they contain a mixture of completed work, retained
requirements, and intentionally excluded future work.

```md
| Historical item | Source plan | Classification | Code evidence | Test evidence | Action |
| --- | --- | --- | --- | --- | --- |
| Example item | plan path and section | implemented | file and line | command and artifact | retain |
```

Allowed classifications:

- `implemented`;
- `retained`;
- `superseded`;
- `deferred`;
- `not applicable`.

Do not rewrite history by deleting old plans or describing unimplemented work
as inherited behavior.

### Checkpoint Log Template

```md
## CP-X: checkpoint title

- Commit:
- Allowed edit set reviewed:
- Evidence reviewed:
- Reviewer verdicts:
- Assumptions falsified:
- Scope drift found:
- Unresolved risks:
- Blocking findings:
- Evidence added:
- Choices made since previous checkpoint, including retained defaults:
- Decision records still proposed and why they do not authorize current edits:
- Verdict: PROCEED | HOLD | STOP
- Next allowed edit set:
```

## Mandatory Checkpoint Sequence

These checkpoints supplement the reflection gates. The checkpoint is the
evidence bundle and verdict; the reflection gate is the deliberate reassessment
of direction.

| Checkpoint | Must be complete before | Required independent review | Minimum evidence |
| --- | --- | --- | --- |
| `CP-0 Baseline` | Any production-code edit | Branch, SHA, dirty state, historical-plan reconciliation, merge-tree, and affected-file inventory review | Frozen base, target strategy, baseline tests, write manifest |
| `CP-1 Physics` | State-layout or pusher edit | Physics and numerics review | Accepted normalized equations, units, charge semantics, pusher spike, oracle plan |
| `CP-2 Data Contract` | Restart or output-format edit | Architecture and migration review | State schema, authoritative-versus-derived values, restart matrix, output dictionary |
| `CP-3 Minimal Serial` | MPI or AMR qualification | Physics and test-adversary review | Zero-field, uniform-`B`, prescribed-`E`, crossed-field, high-`gamma`, work-closure, and convergence metrics |
| `CP-4 Solver Coupling` | Any production ideal-MHD solver-coupling qualification claim | Physics, MHD integration, and test-adversary review | Gather tests, pointwise ideal-field checks, temporal-model statement, prescribed-versus-coupled comparison, and time-dependent manufactured convergence before any dynamic-MHD order claim |
| `CP-5 Restart And Diagnostics` | Broad infrastructure qualification | Restart and output specialist review | Legacy policy, new-format round trips, inspector decode, matched-checkpoint policy, output validation |
| `CP-6 AMR And MPI` | GPU handoff | Migration and decomposition review | `1`, `2`, `4`, and `8` rank evidence where supported, empty ranks, SMR, forced AMR, boundary cases |
| `CP-7 Final Cold Review` | Merge recommendation; opens only after `CP-4 Solver Coupling` records `PROCEED` | Fresh sub-agents with no implementation ownership | Full diff, all records, all evidence, merge-tree refresh, docs overlay, GPU status, residual-risk ledger |

Every `CP-*` entry must name the next permitted edit set. No downstream edit
set opens until the prerequisite checkpoint and reflection gate both record
`PROCEED`. Rerun affected checkpoints after rebases, conflict resolution,
schema changes, or pusher changes.

`CP-7 Final Cold Review` applies to one candidate SHA and its evidence bundle.
Any subsequent source, test, documentation, ledger, acceptance-criterion, or
evidence change invalidates `CP-7` and dependent qualification verdicts until
the changed candidate is rechecked.

### Phase Entry Map

Open only the named edit set after all listed prerequisites record `PROCEED`.
The solver-coupled portion of Phase 4 may be opened experimentally after its
listed prerequisites, but it cannot receive a production qualification claim
until `CP-4 Solver Coupling` records `PROCEED`.

| Phase or edit set | Opened by |
| --- | --- |
| `0` baseline and ledger | Guide creation; documentation edits only |
| `1` state helpers | `CP-0 Baseline`, `CP-1 Physics`, and `RG-001` |
| `2` runtime contract | accepted Phase 1 done-claim review and `RG-002` |
| `3` gathers and prescribed fields | accepted Phase 2 done-claim review and `RG-003` |
| `4a` prescribed-field relativistic push | accepted Phase 3 done-claim review and `RG-004` |
| `4b` experimental solver coupling | `CP-3 Minimal Serial` and `RG-005` |
| `5` subcycling and outer timestep | `CP-3 Minimal Serial` and `RG-005` |
| `6` restart | `CP-2 Data Contract`, accepted Phase 5 done-claim review, and `RG-006` |
| `7` diagnostics | `CP-2 Data Contract`, accepted Phase 6 done-claim review, and `RG-007` |
| `8` MPI and AMR | `CP-5 Restart And Diagnostics` and `RG-008` |
| `9` GPU and handoff | `CP-6 AMR And MPI` and `RG-009` |
| merge recommendation | `CP-7 Final Cold Review` and `RG-010` |

`CP-4 Solver Coupling` must record `PROCEED` before `CP-7 Final Cold Review`
opens and before the branch can receive `MERGE READY`.

## Phase Write Manifest

This is the default file boundary. Reopen it through a decision record when a
phase genuinely needs another surface. An unexpected edit outside the current
row is a stop signal, not a cleanup task to hide in the same commit.

| Phase | Default allowed write surfaces | Surfaces to inspect but avoid editing casually |
| --- | --- | --- |
| `0` baseline | Guide, ledger, test evidence artifacts | All affected sources; integration target |
| `1` state helpers | `src/athena.hpp`, `src/particles/particles.hpp`, narrowly scoped helper implementation, helper tests | Restart, outputs |
| `2` runtime contract | `src/particles/particles.hpp`, `src/particles/particles.cpp`, initialization inputs, parser tests, user docs | Existing input decks |
| `3` gathers and prescribed fields | `src/particles/particles_pushers.cpp`, narrowly scoped pgen or test harness, gather tests; conditionally `src/pgen/pgen.cpp`, `src/pgen/pgen.hpp`, and `src/CMakeLists.txt` if a registered harness is accepted | MHD edge-field implementation |
| `4` relativistic push | `src/particles/particles_pushers.cpp`, kernel tests, prescribed-field decks | Legacy `boris` path |
| `5` subcycling and outer timestep | `src/particles/particles_pushers.cpp`, `src/particles/particles_tasks.cpp`, `src/pgen/part_random.cpp`, timestep tests; conditionally `src/mesh/mesh.cpp` after ownership review | MHD timestep ownership |
| `6` restart | `src/outputs/rst_prtcl.cpp`, `src/particles/particles.cpp`, `src/pgen/part_random.cpp`, inspector, restart tests | Mesh restart format |
| `7` diagnostics | Particle-output files, analysis scripts, output tests, metadata docs; conditionally `src/outputs/outputs.cpp`, `src/outputs/outputs.hpp`, `src/outputs/basetype_output.cpp`, and `src/outputs/derived_variables.cpp` after audit | Legacy column semantics |
| `8` MPI and AMR | Migration tests and inputs; production sources only for isolated defects | Generic particle exchange and mesh-refinement infrastructure |
| `9` GPU and handoff | Device-specific tests, docs, evidence artifacts; optimized source only after profiling | Public `gh-pages` overlay |

## Current Code Baseline

The implementation agents must verify these observations against their working
commit before editing. They describe the prerequisite branch when this guide
was drafted.

### Existing Particle State

The current particle record contains three integer fields and fourteen real
fields. The real fields include position, velocity, the `IPM` control value,
sampled magnetic field, and displacement diagnostics. The current CR tracer
path evolves velocity while preserving speed under magnetic rotation.

Relevant surfaces:

| File | Current role | Relativistic review obligation |
| --- | --- | --- |
| `src/athena.hpp` | Particle data indices | Decide whether to replace velocity storage, add momentum storage, cache `gamma`, add sampled electric field, and version semantics |
| `src/particles/particles.hpp` | Runtime controls and particle-module interface | Add explicit mode, units, light-speed, validation controls, and any helper contracts |
| `src/particles/particles.cpp` | Input parsing, initialization support, restart sizing | Remove hard-coded assumptions only through a versioned design; validate mode combinations |
| `src/particles/particles_pushers.cpp` | Gather, subcycling, magnetic rotation, displacement update | Add reviewed relativistic state conversion, fluid-velocity gather, ideal electric field, pusher, and changing-momentum timestep logic |
| `src/particles/particles_tasks.cpp` | Task ordering | Prove that field time-centering and particle push scheduling match the intended method |
| `src/pgen/part_random.cpp` | Random-particle initialization and restart loading | Add explicit relativistic initialization and legacy conversion rules |
| `src/bvals/bvals_part.cpp` | Particle exchange | Verify that added real fields migrate correctly across blocks, ranks, and AMR rebuilds |
| `src/outputs/rst_prtcl.cpp` | Particle restart writing | Replace fixed record assumptions with a documented versioned contract |
| `src/outputs/df_prtcl.cpp` | Particle spectrum diagnostics | Redefine or add momentum, `gamma`, and relativistic kinetic-energy outputs without semantic ambiguity |
| `src/outputs/pos_prtcl.cpp` | Particle-position output | Preserve stable semantics and add regression coverage |
| `src/outputs/dxhist_prtcl.cpp` | Displacement-history output | Preserve stable semantics if velocity remains physical; otherwise migrate explicitly |
| `src/outputs/track_prtcl.cpp` | Tracked-particle output | Repair or isolate confirmed pre-existing repeated-buffer and byte-offset defects before accepting it as evidence |
| `src/outputs/vtk_prtcl.cpp` | Particle VTK output | Repair or isolate confirmed pre-existing integer-field naming gaps and audit new-state coverage before accepting it as evidence |
| `scripts/particles/cr_tracer_inspect.py` | Restart inspection | Decode each supported record version and report explicit field names |
| `tst/test_suite/particles/` | Particle tests | Preserve fixed-energy controls and add independent relativistic validation |
| `tst/run_test_suite.py` | Test discovery and naming | Follow existing `_cpu`, `_mpicpu`, and `_gpu` conventions |
| `inputs/particles/` | CR tracer decks | Add explicit opt-in examples with documented normalization |

### Existing `IPM` Hazard

In the prerequisite CR tracer path, the magnetic rotation divides by
`pr(IPM,p)`, so `IPM` behaves like an inverse gyro coefficient or a
mass-to-charge-like control in the normalized equations. Negative values also
have existing special behavior in the magnetic-only path. In the local `PIC`
comparison branch, `IPM` is instead read as `q_over_m`.

Do not silently reinterpret `IPM`. Before supporting physical charge signs,
decide whether to:

1. keep a mass-to-charge convention and rename it;
2. migrate to charge-to-mass and convert legacy values;
3. introduce separate explicit fields;
4. initially support only positive-charge passive tracers while retaining a
   compatibility mode for legacy sentinel behavior.

This is a mandatory decision record because an unnoticed reciprocal or sign
change can produce plausible-looking but incorrect trajectories.

### Existing Restart Hazard

The prerequisite branch contains fixed-size particle restart assumptions:
three integer values plus fourteen real values, for seventeen values per
particle record. Those assumptions appear in particle initialization, restart
writing, and the Python inspector.

The relativistic work must not scatter a new magic record size through the
code. Define a versioned schema, document each field, serialize integers as
integers, validate magic/version/checksum, and record checkpoint ID, time,
cycle, topology policy, and shard manifest. Preserve a compatibility reader
only where conversion is safe. Reject precision-unsafe legacy conversion and
add automated round-trip, stale-sidecar, missing-sidecar, nonintegral-count,
large-tag, and legacy-load tests.

The current restart path reads a shard derived from the current MPI rank and
restores saved `PGID` values. It does not automatically support changing rank
count. `DR-009` must choose either manifest-based old-shard discovery plus
redistribution onto the current topology, or deterministic rejection of rank-
count changes. Test the selected policy; do not promise transparent support by
default.

### Existing MHD Electric-Field Context

AthenaK MHD already implements constrained transport for the induction
equation, with an edge-centered electric-field-like variable and ideal-MHD
`-(v x B)` semantics. The passive particle feature should not casually reuse
that edge-centered field:

- its staggering differs from the current cell-centered particle gather;
- its task-stage availability and time centering need analysis;
- interpolation from edges is a separate numerical design;
- the current particle push occurs before the time integrator stages.

The recommended first path is to gather fluid velocity and magnetic field
consistently at the particle position and construct a pointwise reconstructed
ideal field there. It is not the edge-centered CT electric field: interpolation
and crossing do not generally commute. Record this choice, its limitations, and
grid-convergence evidence. If exact reuse of the CT electromotive force is
later desired, implement it as a separate reviewed design.

The current particle push runs before the MHD integrator stages. A reference
path that gathers at a spatial orbit midpoint from grid fields frozen at
`t^n` can establish second-order kernel behavior for static or genuinely
time-centered fields. It does not establish second-order solver-coupled
temporal accuracy for evolving MHD fields. Qualify that separately with a
time-dependent manufactured solution. Until then, qualify frozen-background
coupling only for `time/evolution = static` or explicitly frozen prescribed-
field harnesses. Treat kinematic advection and dynamic MHD as evolving
backgrounds requiring separate temporal qualification.

### Local `PIC` Branch As Bounded Prior Art

The local `PIC` branch already demonstrates several useful ideas:

- interpolate fluid velocity alongside magnetic field;
- compute sampled ideal `cE = -u x B`;
- split an electromagnetic update into half electric kick, magnetic rotation,
  and half electric kick;
- store sampled electric fields;
- collect energy and momentum-change diagnostics.

It also contains coupled-feedback, deposition, expanding-box, species, and
other broad changes that are outside this feature. Its pusher is
nonrelativistic velocity-based at the time of this guide. Use it for comparison
and review questions, not as a wholesale merge source.

### Candidate Data-Layout Strategy To Evaluate

One compatibility-oriented option is to keep the existing fourteen real fields
as a stable prefix, preserve `IPVX`, `IPVY`, and `IPVZ` as physical velocity,
and append momentum, optional cached `gamma`, and sampled electric-field
fields. That can reduce churn in migration, displacement, and legacy output
paths because MPI and AMR exchanges already iterate over runtime `nrdata`.

This is a candidate, not an accepted design. A duplicated velocity and momentum
representation creates stale-cache risk. Evaluate it against a momentum-only
authoritative layout, decide which quantities are authoritative versus derived,
and record invariant checks and restart implications in `DR-002`.

### Additional Pre-Existing Risks To Audit Separately

- The particle timestep constraint must be refreshed after momentum changes;
  subcycling alone is not automatically a sufficient outer-MHD-step policy.
- Mesh restart and particle restart are separate artifacts. Decide whether the
  new mode requires matched time and cycle metadata and synchronized
  checkpoints.
- The inherited migration and AMR paths wrap particle coordinates
  periodically. Treat the first relativistic mode as `strictly_periodic`
  unless a separately reviewed boundary implementation is accepted.
- The legacy `lin_legacy` gather uses staggered magnetic components and does not
  automatically provide a clean matching fluid-velocity gather. Reject it
  initially unless a separate decision record defines and validates a matching
  fluid gather.
- Repair or isolate `track_prtcl.cpp` before using tracked output as acceptance
  evidence. Direct source inspection found that each loop writes `data[0]`
  instead of the current particle slice and that record offsets omit the float
  byte size.
- Repair or isolate `vtk_prtcl.cpp` before using VTK output as acceptance
  evidence. Direct source inspection found that it iterates all integer fields
  while naming only `PGID` and `PTAG`, leaving `PSP` output structurally
  suspect.
- Review resistivity, orbital-advection, shearing-box, SRMHD, and GRMHD modes
  explicitly. The initial feature should reject unsupported combinations rather
  than assuming Newtonian ideal-MHD semantics apply.

## Mandatory Pre-Code Decisions

Open all of the following records in Phase 0. Accept each record before its
first dependent production edit. Experiments and small prototypes are allowed,
but production code changes must not outrun accepted decisions. A checkpoint
must list records that remain proposed and explain why they do not authorize
the current edit set.

| ID | Decision | Questions that must be answered | Initial recommendation |
| --- | --- | --- | --- |
| `DR-000` | Base, target, and stacking strategy | What base SHA is frozen? What eventual integration target is intended? What merge-tree overlaps exist? When should rebasing happen? | Record the current stacked base and a deliberate conflict strategy before implementation |
| `DR-001` | Physical normalization and reduced-light model | Does the particle equation use `E`, `cE`, `q/m`, `m/q`, physical `c`, or reduced `C`? Is `C` a deliberate model approximation? Where does it enter the state relation, force prefactor, sampled-field mapping, and work diagnostic? How do MHD code units map to the pusher? | Write the exact physical equations first, then any reduced-light model explicitly; test dimensional consistency before coding |
| `DR-002` | Authoritative particle state | Store velocity, momentum, mass-normalized momentum, or another variable? Is `gamma` stored or derived? | Store momentum-like state; derive `gamma` and velocity; cache only if justified |
| `DR-003` | Relativistic pusher | Standard relativistic Boris, Vay, Higuera-Cary, or another reviewed method? | Run a focused comparison spike and select only after preregistered criteria are reviewed |
| `DR-004` | Runtime compatibility mode | New pusher name or silent replacement of existing `boris`? | Add a new explicit opt-in mode and retain legacy `boris` as a fixed-energy control |
| `DR-005` | Charge and species semantics | What exactly does `IPM` mean? Are negative charges supported now? What replaces legacy negative-value behavior? | Introduce explicit semantics and avoid overloaded sentinel values |
| `DR-006` | Gather policy | Linear, trilinear, TSC, legacy policy, or selectable policy? Gather `u` and `B` separately then cross, or gather a precomputed field? | Gather `u` and `B` with the same documented policy, then construct the explicitly normalized pointwise ideal field at the particle |
| `DR-007` | Temporal centering | Frozen fields over each particle substep, midpoint fields, old/new MHD fields, or stage-coupled fields? | Begin with a clearly documented frozen-background midpoint gather; measure convergence before widening scope |
| `DR-008` | Dimensional semantics | 3-D only initially, or intentional 1-D/2-D position plus 3-vector velocity and field behavior? | Either restrict the first mode to 3-D or explicitly validate a 2D3V contract |
| `DR-009` | Restart schema | How are new records versioned? How do legacy velocity records load? What metadata is required? | Add a versioned extensible schema and explicit legacy conversion policy |
| `DR-010` | Diagnostics migration | Which old outputs retain meaning? Which get new names? Which new values are required? | Add unambiguous relativistic outputs; do not silently redefine old columns |
| `DR-011` | Subcycling constraints | Bound cell crossing, block crossing, relativistic gyro-angle, electric kick, and fractional momentum or energy change? | Add conservative bounds and validate each bound independently |
| `DR-012` | Failure policy and validity domain | What invalid inputs and runtime states abort immediately? Are electric-dominated sampled states supported, rejected, or warned? What maximum `|u_fluid| / c_model` is intended? | Fail fast for nonfinite state, invalid `c_model`, unsupported dimensions, inconsistent mode combinations, superluminal reconstruction, and states outside the accepted model-speed domain |
| `DR-013` | Electric-field scope | Only ideal Newtonian MHD initially, or additional field models? If the initial mode assumes a magnetically dominated background, how is `|cE| < c_model |B|` enforced? | Keep the production feature ideal and passive; require or deliberately reject `|u_fluid,perp| < c_model` as appropriate; use test-only prescribed fields only if needed for validation |
| `DR-014` | Position integration | Drift-kick-drift, kick-drift-kick, existing position update, or another second-order split? | Choose a method consistent with the field gather time and prove second-order behavior |
| `DR-015` | Work accumulation | What discrete `q E dot v`, `(q / c) cE dot v`, or explicitly derived reduced-model quantity is accumulated and compared to kinetic-energy change? Is it serialized? | Define a pusher-consistent explicit sampled-field quadrature and preregister closure tolerances |
| `DR-016` | Outer timestep refresh | How is the mesh-level particle timestep constraint updated after acceleration? | Add a reviewed refresh point and prove that the next MHD step consumes current particle limits |
| `DR-017` | Particle boundary semantics | Does the first relativistic mode support only periodic coordinates, or are reflecting and other physical boundaries implemented deliberately? | Start with `strictly_periodic` and fail fast otherwise; widen only through separate review |
| `DR-018` | Compatibility matrix | Which EOS, static/kinematic/dynamic evolution, ion-neutral, turbulence-forcing, resistive, shearing-box, SRMHD, GRMHD, and dynamical-GR combinations are supported? | List every combination explicitly; fail fast for every unqualified opt-in mode |

## Relativistic Pusher Selection Spike

Do not choose a pusher by familiarity alone. Implement or script a small,
isolated comparison before wiring the production path.

Evaluate a justified candidate set that includes at least two credible methods.
Record why any named method is excluded. Use implementation, scripting, or
literature-backed elimination to the depth needed for a defensible selection.
Named candidates include:

1. Standard relativistic Boris.
2. Vay's relativistic method.
3. Higuera-Cary's structure-preserving second-order method.

The comparison must include:

- uniform `B`, zero `E`: energy conservation and gyro-phase error;
- crossed uniform `E` and `B`: correct drift velocity;
- exact crossed-field force cancellation;
- high-`gamma` behavior;
- timestep convergence;
- reversibility where expected;
- phase-space-volume behavior or literature-backed elimination;
- a long-time resonance stress case;
- boosted-frame force-cancellation behavior;
- finite-state behavior near configured limits;
- GPU-suitable algebra and branch structure;
- operation count and clarity;
- reference values produced independently from the AthenaK kernel.

Do not copy formulas from memory. Read the primary references, transcribe the
chosen algorithm into the decision record, and ask a physics reviewer to check
the algebra before integration.

## Planned Implementation Phases

Each phase ends with tests, independent reviews, a reflection record, and a
small commit. Do not start the next phase merely because the code compiles.

### Phase 0: Freeze Baseline And Open The Ledger

#### Work

1. Record base commit, compiler, build configuration, MPI mode, and dirty state.
2. Run existing particle tests and the fixed-energy CR tracer accuracy ladder.
3. Identify all hard-coded particle record sizes and all consumers of particle
   velocity, `IPM`, sampled `B`, displacement, spectra, and restart fields.
4. Open `DR-000` through `DR-018`.
5. Complete the pusher selection spike.
6. Write the first reflection record.
7. Complete the historical-plan reconciliation ledger.
8. Run and archive the merge-tree audit against the selected target.

#### Required Evidence

- Existing fixed-energy suite results.
- Search inventory for affected symbols and record sizes.
- Pusher comparison artifact.
- Reviewed normalized equations.
- Initial compatibility matrix.
- Historical-plan reconciliation ledger.
- Merge-tree result and accepted conflict strategy.

#### Independent Review

Spawn:

- a physics reviewer for normalized equations and pusher candidates;
- an integration reviewer for code-surface inventory and compatibility risks;
- a test reviewer for baseline coverage and missing oracles.

#### Stop Conditions

Stop before production coding if:

- units are still implicit;
- `IPM` semantics remain ambiguous;
- legacy compatibility behavior is not stated;
- integration target and conflict strategy are unresolved;
- the selected pusher has not passed a crossed-field drift test;
- no independent reviewer has checked the algebra.

### Reflection Gate RG-001: Is The Feature Still Passive And Bounded?

Reassess whether the first branch can remain a one-way MHD-to-particle feature.
If a desired test requires feedback, separate it from this branch unless there
is a documented reason to widen scope.

### Phase 1: Add Relativistic State Helpers Without Changing Runtime Behavior

#### Work

Add small device-capable helper functions for the chosen state model:

- `gamma(momentum, c_model)`;
- `velocity(momentum, gamma)` or equivalent;
- kinetic energy;
- finite-state checks;
- optional conversion from legacy velocity to new momentum;
- optional conversion from momentum to legacy velocity for compatibility
  outputs.

Keep the existing runtime pusher behavior unchanged. Unit-test the helpers
before integrating them into particle movement.

#### Required Tests

- `gamma >= 1`;
- zero-momentum limit;
- nonrelativistic limit;
- high-momentum limit;
- round-trip conversion within stated tolerance;
- no reconstructed `|v| >= c_model`;
- rejection of invalid or nonfinite `c_model`, state, and conversion input;
- host/device parity if helpers are exercised on both paths.

#### Reflection Gate RG-002: Is The Chosen State Still The Right One?

Review helper complexity, restart implications, and diagnostic semantics. If a
cached `gamma` adds synchronization risk without measured value, remove it.

### Phase 2: Add An Explicit Opt-In Runtime Contract

#### Work

1. Add the new runtime mode without activating it by default.
2. Add light-speed and normalization inputs with fail-fast validation.
3. Add initialization parsing for momentum, `gamma`, energy, or speed only
   after recording precedence rules and the directional-momentum contract for
   scalar initialization inputs.
4. Preserve legacy mode semantics and old input decks.
5. Add user-facing documentation for supported and rejected combinations.
6. Implement the accepted `DR-018` compatibility matrix and fail fast for
   every unqualified opt-in combination.

#### Required Tests

- legacy inputs run unchanged;
- new mode requires valid `c_model`;
- ambiguous initialization combinations fail clearly;
- `|v| >= c_model` initialization fails clearly;
- unsupported dimensional modes fail clearly until validated;
- legacy sentinel or negative-charge cases follow the recorded policy.
- every unsupported EOS, evolution, source-term, frame, and relativity
  combination fails clearly according to `DR-018`.

#### Independent Review

Spawn an integration reviewer specifically to look for silent behavior changes
in old decks and restart paths.

### Reflection Gate RG-003: Is Compatibility Explicit?

Read the input documentation as a new user. Verify that a user cannot
accidentally activate relativistic semantics or misunderstand `IPM`,
`c_model`, physical `c`, reduced `C`,
momentum, velocity, or output units.

### Phase 3: Gather Fluid Velocity And Construct Sampled Ideal Electric Field

#### Work

1. Add a narrowly scoped prescribed-field harness for analytical kernel
   validation. Keep it mechanically distinct from the production ideal-MHD
   source path.
2. Extend the existing particle gather framework to interpolate fluid velocity
   using the same stated spatial policy as magnetic field.
3. Construct particle-sampled ideal `E` or `cE` from gathered values.
4. Keep the solver-coupled electromagnetic update disabled or diagnostic-only
   at first.
5. Add optional sampled electric-field output fields only after deciding the
   particle record and output contracts.
6. Verify boundary, ghost-zone, AMR, and lower-dimensional behavior.

#### Required Gather-Only Tests

- uniform `u_fluid`, uniform `B`: exact sampled values;
- spatially varying manufactured `u_fluid` with uniform `B`: expected
  interpolation order;
- spatially varying manufactured `B` with uniform `u_fluid`: expected
  interpolation order;
- `u_fluid || B`: sampled ideal field is zero;
- known perpendicular vectors: sign and component checks;
- pointwise `E dot B` near roundoff for ideal construction;
- serial/MPI parity;
- block-boundary continuity;
- static-refinement and forced-AMR migration behavior;
- every enabled gather policy separately.

Treat pointwise `E dot B` as a debug invariant, not independent validation
evidence. It is nearly automatic when the field is reconstructed as a cross
product, and it does not detect every sign error.

#### Stop Conditions

Stop before enabling acceleration if:

- sampled-field sign is not checked analytically;
- the gather policy differs silently between `u_fluid` and `B`;
- `E dot B` is unexpectedly nonzero;
- MPI or block-boundary behavior differs without explanation.

### Reflection Gate RG-004: Is Pointwise Ideal Construction Adequate?

Review whether gathering `u_fluid` and `B` then crossing them remains the right
first model. Record the tradeoff versus gathering a separately discretized
electric field or CT electromotive force. Do not change staggering or time
centering without a new decision record and tests.

### Phase 4: Enable The Relativistic Electromagnetic Push

#### Work

1. Integrate the selected pusher behind the explicit new runtime mode, first
   against the prescribed-field harness.
2. Use the reviewed momentum state consistently.
3. Derive drift velocity from momentum and `gamma` for position updates.
4. Preserve the legacy fixed-energy path unchanged.
5. Update sampled fields and displacement diagnostics deliberately.
6. Add debug-mode consistency checks that can expose nonfinite state,
   impossible `gamma`, and superluminal reconstructed velocity early.
7. After the prescribed-field analytical suite passes, enable the
   solver-coupled ideal-MHD gather path and rerun the applicable ladder.
8. Keep prescribed-field and solver-coupled evidence separate in reports.

#### Required Kernel Tests

- zero-field straight-line drift;
- pure uniform electric field with independent closed-form `p(t)`, `gamma(t)`,
  and `x(t)` oracles;
- magnetic-only orbit with constant `gamma`;
- uniform magnetic field with independent radius and phase oracles;
- crossed-field drift;
- crossed-field zero-force equilibrium;
- production-path cancellation when particle velocity matches fluid velocity;
- acceleration and deceleration in controlled fields;
- work-energy closure for every accelerating or decelerating case;
- time-dependent manufactured field for solver-coupled temporal-order
  qualification;
- timestep convergence for momentum and trajectory;
- time reversibility where the selected method promises it;
- charge-sign behavior if supported;
- nonrelativistic agreement with a low-speed reference;
- high-`gamma` stability;
- comparison against an independent high-accuracy reference integrator.

Accumulate work through an explicit quadrature of sampled fields. Never
construct the work diagnostic from `Delta K`, because that would make closure
tautological.

Some controlled tests may need a narrowly scoped prescribed-field harness
because ideal MHD cannot express every useful field combination. Keep that
harness clearly separated from the production ideal-MHD mode.

#### Independent Review

Spawn:

- a physics reviewer to independently derive and check the pusher;
- a test reviewer to challenge every analytical oracle;
- an integration reviewer to look for accidental changes to legacy mode,
  displacement, periodic-boundary behavior, nonperiodic rejection, and task
  scheduling.

### Reflection Gate RG-005: Does The Pusher Merit Integration?

Do not proceed on the basis of one orbit plot. Review quantitative errors,
convergence rates, drift accuracy, limit behavior, and reviewer findings. If
the selected algorithm performs poorly, revisit `DR-003` rather than tuning
tests around it. Classify prescribed-field evidence and solver-coupled
ideal-MHD evidence separately.

### Phase 5: Redesign Subcycling For Changing Momentum

#### Work

The old subcycling policy was designed for fixed-speed magnetic trajectories.
The relativistic mode needs bounds that remain valid as momentum and `gamma`
change.

Review and implement limits for:

- maximum cell crossing;
- maximum block crossing;
- relativistic gyro-angle per substep;
- maximum electric kick;
- maximum fractional momentum change;
- maximum fractional kinetic-energy change, if useful;
- extreme or nonfinite substep requests;
- global caps that prevent runaway subcycle counts.
- mesh-level timestep refresh after acceleration.

State whether fields are frozen over a full MHD step, regathered at each
particle substep, or handled another way. Specify these independently:

- which temporal grid state is sampled;
- whether fields are spatially regathered at every particle substep;
- whether constraints are recomputed adaptively after electric kicks or bounded
  conservatively over the full outer step.

The current fixed-energy path computes the subcycle count once before its loop.
Do not assume that remains conservative after momentum-changing electric kicks.
Test the accepted policy directly.

#### Required Tests

- increasing electric field triggers more substeps where intended;
- high `gamma` changes the gyro bound correctly;
- stricter limits converge toward the reference solution;
- loose limits fail or warn according to policy;
- subcycle caps fail clearly rather than hanging;
- MPI ranks choose behavior consistently;
- AMR migration remains correct when particles cross blocks during many
  substeps.
- the next outer MHD step consumes a current particle timestep constraint.

### Reflection Gate RG-006: Is The Temporal Model Still Defensible?

Review whether frozen MHD fields over the particle step are adequate for the
target science. If not, stop and design a stage-coupled follow-up rather than
quietly mixing stage fields into the passive pusher.

### Phase 6: Version Restart And Inspection Contracts

#### Work

1. Define a named, versioned particle restart schema.
2. Document every integer and real field.
3. Write new-format records with explicit metadata.
4. Add a new-format reader.
5. Add an intentional legacy-format reader or a clear rejection path.
6. Update `scripts/particles/cr_tracer_inspect.py`.
7. Decide whether restart preserves momentum only, caches `gamma`, or stores
   derived values for inspection without treating them as authoritative.
8. Require or reject paired mesh and particle checkpoints explicitly. If paired
   checkpoints are required, validate matching time and cycle metadata.
9. Choose and implement either manifest-based rank-count-change redistribution
   or deterministic rank-count-change rejection.

#### Required Tests

- legacy fixed-energy restart round trip;
- new relativistic restart round trip;
- legacy-to-new conversion according to policy;
- restart continuation versus uninterrupted run;
- MPI-rank-count change across restart;
- rank-count-change support or rejection exactly follows `DR-009`;
- restart after AMR migration;
- inspector decode of every supported version;
- mismatched mesh and particle checkpoint rejection if paired checkpoints are
  required;
- stale or missing sidecar rejection;
- nonintegral particle counts and precision-unsafe large legacy tags;
- corrupted metadata and truncated-record rejection;
- clear diagnostic for unsupported future versions.

#### Independent Review

Spawn a restart specialist reviewer. Ask them to search again for every hard-
coded record size and decode assumption across C++ and Python.

### Reflection Gate RG-007: Is The Migration Contract Honest?

Check that legacy restarts are either faithfully converted or explicitly
rejected. Do not label a lossy or ambiguous conversion as backward compatible.

### Phase 7: Add Unambiguous Relativistic Diagnostics

#### Work

Add or revise diagnostics only through `DR-010`. Candidate outputs include:

- momentum components;
- momentum magnitude;
- Lorentz factor;
- derived velocity components and speed;
- relativistic kinetic energy;
- sampled `B`;
- sampled ideal `E` or `cE`;
- `E dot B`;
- instantaneous `q E dot v`, `(q / c) cE dot v`, or the explicitly derived
  reduced-model quantity, according to `DR-001`;
- accumulated work and work-energy residual;
- displacement diagnostics;
- subcycle count and active limiting constraint;
- full-orbit applicability diagnostics such as `r_L / dx` and gyro-angle
  resolution;
- finite-state failure counters.

The existing spectrum code currently uses speed and nonrelativistic
`0.5 * speed^2` energy proxies. Preserve old proxies for legacy mode or version
their names. Do not silently call a relativistic quantity by an old ambiguous
name.

Name the current traps explicitly in `DR-010`: legacy `pspec p` and `psamp p`
mean speed, legacy energy means `0.5 * v^2`, and legacy `psamp mass` emits the
overloaded `IPM` value.

Audit every particle-facing output path, including `ppd`, `pmom`, `df`,
`pspec`, `pspec2`, `psamp`, position, displacement history, tracked-particle,
and VTK outputs. Classify each as preserved, extended, versioned, repaired in a
separate prerequisite, or intentionally unsupported for the new mode.

#### Required Tests

- diagnostic values match independently reconstructed values;
- magnetic-only run preserves `gamma` and kinetic energy;
- acceleration and deceleration change kinetic energy with the expected sign;
- accumulated work closes against relativistic kinetic-energy change within
  preregistered tolerance;
- restart does not alter diagnostic meaning;
- spectra bins and labels match documented quantities;
- legacy output retains its documented meaning.

### Reflection Gate RG-008: Can A Scientist Interpret The Outputs Correctly?

Have an independent reviewer read output metadata and analysis scripts without
the implementation patch open. If they cannot determine units and semantics,
the output contract is not done.

### Phase 8: MPI, Block, SMR, And AMR Qualification

#### Work

Exercise the validated kernel through the full particle-migration path:

- multiple mesh blocks;
- block boundaries;
- static mesh refinement;
- forced refinement and derefinement;
- MPI rank migration;
- rank-count-change support or rejection according to `DR-009`;
- particles near periodic boundaries;
- displacement tracking across migrations;
- changing momentum during migration.

Treat `src/bvals/bvals_part.cpp`, `src/mesh/mesh_refinement.cpp`, and
`src/particles/particles.cpp::RemapAfterAMR()` as required review surfaces even
if generic loops already move all real fields.

#### Required Tests

- serial versus MPI decomposition invariance;
- decomposition invariance for `1`, `2`, `4`, and `8` ranks, including an
  empty-rank case where supported;
- one versus multiple blocks;
- uniform versus static-refined mesh where the physical field is uniform;
- forced refine/derefine cycles;
- restart before and after migration;
- periodic wrapping;
- fail-fast rejection of nonperiodic boundaries unless `DR-017` explicitly
  widens scope;
- no stale state after AMR reconstruction;
- particle count and unique identity preservation;
- finite-state diagnostics remain clean.

### Reflection Gate RG-009: Is Infrastructure Noise Hiding Physics Error?

Compare kernel-only tests with end-to-end tests. If errors appear only after
MPI, AMR, or restart, isolate infrastructure defects before adding new
physics. If errors exist in both, return to the pusher or gather phase.

### Phase 9: GPU, Performance, Documentation, And Handoff

#### Work

1. Run device builds and tests on the intended GPU backend.
2. Confirm that helper functions and pusher kernels remain device-capable.
3. Profile before optimizing.
4. Add optimization only in isolated commits with differential tests.
5. Update user documentation, input examples, test documentation, and handoff
   notes.
6. Archive final evidence and independent-review responses.
7. Follow repository test naming conventions for `_cpu`, `_mpicpu`, and
   `_gpu` coverage.
8. Archive equivalent manual qualification evidence when the intended
   integration target is not covered by repository Actions. At guide creation,
   the checked-in workflow triggers only for `main`.

#### Required Tests

- CPU debug suite;
- CPU release suite;
- MPI debug and release suite;
- GPU smoke and analytical suite;
- GPU MPI smoke where available;
- restart and inspector suite;
- AMR stress suite;
- deep AMR `divB` regression;
- legacy fixed-energy regression suite;
- documentation build;
- strict Sphinx validation with warnings as errors in a temporary
  `origin/gh-pages` documentation overlay before public integration;
- deterministic analysis rerun from archived artifacts.

### Reflection Gate RG-010: Final Scope And Scientific Readiness Review

Before merging, ask:

- Did the branch remain passive?
- Are old semantics preserved or intentionally migrated?
- Are units explicit?
- Are the pusher and gather independently validated?
- Can restarts be trusted?
- Are diagnostics scientifically interpretable?
- Did independent reviewers close all blocking findings?
- Are deferred features clearly listed?
- Is the next follow-up branch obvious and bounded?

Do not describe the branch as full MHD-PIC, SRMHD, or production scientific
readiness unless a separate qualification plan supports that claim.

## Validation Matrix

The matrix below is the minimum required ladder. Add cases when failures expose
gaps.

| Layer | Case | Main assertion | Compatibility role |
| --- | --- | --- | --- |
| Existing regression | Legacy magnetic-only uniform orbit | Existing fixed-energy behavior remains unchanged | Mandatory control |
| Existing regression | Legacy CR tracer accuracy suite | Existing interpolation and subcycling claims still hold | Mandatory control |
| Helper unit | Momentum to `gamma` and velocity | Correct limits and finite state | New |
| Helper unit | Legacy velocity conversion | Explicit migration behavior | New |
| Gather unit | Uniform `u_fluid` and `B` | Exact sampled ideal field | New |
| Gather unit | Manufactured spatial fields | Expected interpolation order | New |
| Gather unit | `u_fluid || B` | Sampled ideal field vanishes | New |
| Kernel analytic | Zero field | Straight drift | New |
| Kernel analytic | Uniform `B`, zero `E` | Constant `gamma`, orbit phase convergence | New |
| Kernel analytic | Crossed uniform `E` and `B` | Correct drift | New |
| Kernel analytic | Controlled acceleration | Momentum and energy increase correctly | New |
| Kernel analytic | Controlled deceleration | Momentum and energy decrease correctly | New |
| Kernel analytic | Accelerating and decelerating work closure | `Delta(gamma m c_model^2)` matches explicit sampled-field work quadrature in the accepted model | New |
| Kernel analytic | High `gamma` | Stable finite trajectory | New |
| Kernel differential | Selected pusher versus high-accuracy oracle | Error and convergence stay within preregistered bounds | New |
| Subcycle | Bound activation | Correct constraint controls each regime | New |
| Restart | New-format continuation | Matches uninterrupted run | New |
| Restart | Legacy policy | Converts or rejects exactly as documented | Compatibility |
| Restart | Rank-count policy | Redistributes through a manifest or rejects clearly according to `DR-009` | Compatibility |
| Restart | Integer safety and synchronization | Integer fields, checkpoint metadata, sidecar, checksum, and shard manifest validate | New |
| Inspector | Every supported record version | Correct field decode and metadata | Compatibility |
| MPI | Decomposition variants | Particle results agree within stated tolerance | New |
| SMR | Uniform physical fields | Refinement topology does not change solution materially | New |
| AMR | Forced refine/derefine | Particle identity and solution remain valid | New |
| Boundary | Periodic wrapping and nonperiodic rejection | Correct position and momentum behavior within `DR-017` scope | New |
| Output | Momentum, `gamma`, energy, sampled fields | Values and labels are unambiguous | New |
| GPU | Device execution | Results match CPU within stated tolerance | New |
| Docs | Sphinx build | New guide and user docs render cleanly | New |

## Failure Injection Checklist

Before declaring the feature complete, deliberately test:

- configured `c_model <= 0`;
- missing required model-speed input in relativistic mode;
- nonfinite initialization;
- initialization at or above `|v| = c_model`;
- inconsistent momentum and `gamma` inputs;
- nonfinite sampled fields;
- unsupported negative charge or legacy sentinel combinations;
- unsupported dimensional mode;
- corrupted restart version;
- truncated restart record;
- mismatched mesh and particle checkpoint metadata;
- unreasonable subcycle demand;
- zero magnetic field;
- zero fluid velocity;
- `u_fluid || B`;
- electric-dominated or out-of-domain sampled state under the accepted
  model-speed policy;
- `IPM == 0`;
- reciprocal overflow or underflow in the accepted `IPM` convention;
- particle exactly on or near a block boundary;
- particle crossing multiple blocks under aggressive acceleration;
- restart after a migration-heavy step;
- forced AMR transitions;
- MPI decomposition changes;
- single-precision builds if supported.
- precision-unsafe legacy integer tags and nonintegral legacy counts.

Every injected failure should either produce a documented valid result or fail
fast with an actionable message.

## Full-Orbit Applicability Trigger

This feature remains a full-orbit tracer method. Report `r_L / dx`, gyro-angle
resolution, and any other accepted applicability metric in diagnostics and
qualification artifacts. Define a revisit trigger for regimes where the
full-orbit method is unresolved or prohibitively expensive. Treat a
guiding-center approximation as a separately reviewed follow-up, not an
optimization hidden inside this branch.

## Suggested Commit Sequence

Use this as a planning aid, not an excuse to skip reviews.

1. `docs: add relativistic CR design ledger and normalized equations`
2. `tests: add standalone relativistic state and pusher reference harness`
3. `particles: add device-capable relativistic state helpers`
4. `particles: parse opt-in passive relativistic CR mode`
5. `particles: gather fluid velocity and diagnose ideal electric field`
6. `tests: add manufactured gather and ideal-field checks`
7. `particles: add reviewed relativistic electromagnetic pusher`
8. `tests: add analytic pusher convergence and drift suite`
9. `particles: extend relativistic subcycling constraints`
10. `outputs: version relativistic particle restart schema`
11. `scripts: decode versioned relativistic CR restart records`
12. `outputs: add unambiguous relativistic diagnostics`
13. `tests: add restart MPI AMR and failure-injection qualification`
14. `docs: finalize relativistic CR user guide and handoff evidence`

At each commit boundary, rerun the smallest relevant suite immediately and
apply Rule 4 before further feature edits. At each phase boundary, rerun the
broader suite and broaden the independent-review surface.

## Required Independent Review Prompts

Use prompts like these. Include the commit hash and evidence paths.

### Physics And Numerics Review

```text
Review this relativistic passive CR tracer milestone skeptically. Independently
derive the normalized equations, inspect the state conversion and pusher
algebra, check the ideal-MHD electric-field sign and units, identify unsupported
regimes, and challenge the analytical tests. Do not assume passing tests prove
the implementation. Return blocking findings, nonblocking risks, and missing
oracles with file and line references.
```

### Integration And Migration Review

```text
Review this relativistic CR tracer milestone for hidden integration defects.
Search every particle-record consumer, restart reader/writer, inspector,
diagnostic, migration path, AMR path, MPI path, boundary path, and legacy input.
Look for silent semantic changes, hard-coded record sizes, sign inversions,
stale derived values, and missing compatibility tests. Return blocking findings,
nonblocking risks, and exact follow-up checks.
```

### Test-Adversary Review

```text
Try to falsify the claimed milestone. Inspect the tests and acceptance
thresholds, propose failure injections, identify cases where a test can pass
while the physics is wrong, and distinguish smoke tests from analytical
validation. Return a ranked list of missing or weak tests.
```

### Final Handoff Review

```text
Assume another agent will rush the merge. Audit the full branch against the
implementation guide and decision ledger. Verify that every completed gate has
evidence, every reviewer finding has a disposition, every deferred feature is
explicit, and every user-facing semantic change is documented. Return a merge
recommendation only after listing all residual risks.
```

## Handoff Artifact Contract

The final branch handoff must include:

- exact branch and commit range;
- base commit and prerequisite branch;
- decision ledger;
- reflection records;
- normalized equations and units table;
- selected-pusher derivation reference;
- affected-file inventory;
- build configurations;
- complete command log for qualification tests;
- machine-readable analytical metrics;
- restart schema documentation;
- restart compatibility matrix;
- input examples;
- output-field dictionary;
- CPU/MPI/GPU evidence;
- AMR and migration evidence;
- independent sub-agent prompts, findings, responses, and rechecks;
- current merge-tree status and accepted integration strategy;
- public Pages overlay documentation status;
- `CPU/MPI QUALIFIED`, `GPU QUALIFIED`, and `MERGE READY` verdicts with accepted
  SHAs and any unavailable-backend residual risk;
- known limitations;
- deferred follow-up branches.

## Definition Of Done

The feature is not done when it compiles, when a single orbit looks reasonable,
or when one agent says the patch is finished. It is done only when:

1. The passive physical scope is explicit and preserved.
2. Units, state variables, charge semantics, and light speed are documented.
3. The selected pusher has independent analytical validation.
4. Existing fixed-energy behavior remains available and passes its regression
   suite.
5. Acceleration and deceleration are quantitatively validated.
6. Relativistic work closes against the explicit sampled-field quadrature from
   `DR-015` within the preregistered tolerance.
7. Gather, subcycle, outer-timestep refresh, restart, inspector, output, MPI,
   AMR, and boundary paths
   pass their dedicated tests.
8. Invalid states fail clearly.
9. Every phase has a reflection record.
10. Every nontrivial choice has a decision record.
11. Independent sub-agents have challenged each done claim and rechecked fixes.
12. Documentation renders successfully.
13. The final handoff names residual risk honestly.
14. `CP-4 Solver Coupling` records `PROCEED`.
15. `CP-7 Final Cold Review` records `PROCEED`.
16. The accepted SHA has an explicit `MERGE READY` verdict.

## Primary Numerical References

- A. V. Higuera and J. R. Cary, "Structure-preserving second-order
  integration of relativistic charged particle trajectories in
  electromagnetic fields," *Physics of Plasmas* 24, 052104 (2017),
  <https://arxiv.org/abs/1701.05605>.
- J.-L. Vay, "Simulation of beams or plasmas crossing at relativistic
  velocity," *Physics of Plasmas* 15, 056701 (2008),
  <https://doi.org/10.1063/1.2837054>.

Read the primary papers before implementing formulas. The references guide
pusher selection; they do not replace AthenaK-specific unit, scheduling,
restart, interpolation, and compatibility decisions.
