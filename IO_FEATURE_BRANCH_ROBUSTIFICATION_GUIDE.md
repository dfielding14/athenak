# IO Output Formats And Sharding Robustification Guide

## Document Status

| Item | Value |
| --- | --- |
| Target branch | `feature/io-output-formats-and-sharding` |
| Reviewed branch HEAD | `a273ede97fa93f926df754a8a40047f717f67ac0` |
| Clean base | `origin/main` at `886dd2a1437e45a3a30b3eeebf2adfa838328f73` |
| Deferred Pages baseline inspected during the prior pass | `origin/gh-pages` at `4833aa9341e19861297e330ff02aabfd8001935c` |
| Existing reconstruction guide | `IO_FEATURE_BRANCH_GUIDE.md` |
| Existing selective-integration plan | `IO_FEATURE_BRANCH_INTEGRATION_PLAN.md` |
| Durable compatibility contract | `IO_FORMAT_COMPATIBILITY.md` |
| Decision record to extend | `IO_FEATURE_BRANCH_DECISION_LOG.md` |
| Audit ledger to extend | `IO_FEATURE_AUDIT_LEDGER.md` |
| Status of this document | Execution guide for a new post-review hardening pass |

This guide governs the next implementation pass on
`feature/io-output-formats-and-sharding`. It does not reopen the original Gotham
port indiscriminately. It starts from the reconstructed branch that already exists,
preserves the useful compatibility work, and addresses the locally actionable
robustness, maintainability, testing, and documentation issues found during a fresh
whole-branch review.

The previous audit ledger records a strong local baseline:

- `204` serial IO tests passed, including a CPU smoke run of the GPU-selectable test;
- `59` MPI IO tests passed on one physical node;
- `148` Python reader tests passed;
- repository style checks passed;
- all `27` immutable fixture artifacts passed checksum verification; and
- the deferred Pages package built successfully with warnings treated as errors
  against the recorded `origin/gh-pages` baseline.

Those results remain valuable evidence. They are not permission to rush this pass.
This guide exists because a fresh source audit found issues that the previous audit
cycle missed. Every implementing agent must treat that fact as a process lesson:
passing tests and a prior independent audit are necessary, but neither substitutes
for careful reasoning, explicit decisions, and repeated adversarial review.

## Executive Mandate

Make the IO feature branch suitable for deliberate merge review by:

1. completing end-to-end large-checkpoint safety;
2. checking all branch-added MPI operations consistently;
3. hardening serial positioned IO;
4. defining and enforcing the supported coarsened-binary contract;
5. making file publication guarantees consistent and accurately documented;
6. fully qualifying the advertised diagnostic-variable surface;
7. resolving derived-variable ghost-zone semantics;
8. eliminating avoidable spherical-slice filename collisions;
9. measuring and, if justified, improving node-restart manifest scaling;
10. refactoring touched code so it matches the surrounding codebase more closely;
11. keeping one usable, maintainable Python reader/converter surface;
12. automating the deferred GitHub Pages staging workflow without publishing early;
13. deciding deliberately which large process records belong in the eventual merge;
14. completing real CUDA and scheduler-backed multi-node qualification; and
15. recording every consequential choice, rejected alternative, qualification
    boundary, and reflection-driven change of direction.

## Completion Standard

Do not describe this branch as complete, merge-ready, production-ready, or robust
until all of the following are true:

1. Every locally actionable `ROB-*` item in this guide is either:
   - implemented and independently re-audited;
   - disproven with reproducible evidence and independently accepted as invalid; or
   - explicitly deferred by the user after the risk is documented.
2. Every checkpoint has:
   - a pre-edit reconnaissance record;
   - any required decision-log entries;
   - focused tests;
   - a read-only independent audit by an agent that did not implement the change;
   - correction of blocking findings;
   - a second independent re-audit after correction;
   - a reflection entry; and
   - main-agent verification in the integration worktree.
3. The full serial, MPI, Python-reader, fixture, and style matrices pass again.
4. The GPU-selectable IO regression has executed on a CUDA-capable build.
5. Node-sharded output and native direct restart routing have executed under a
   scheduler on multiple physical nodes, including a genuinely empty or non-owning
   node.
6. The deferred Pages package has been staged into a detached, refreshed
   `origin/gh-pages` worktree by a deterministic helper, reviewed, and built with
   warnings treated as errors.
7. The live `gh-pages` branch remains untouched until after the code branch merges.
8. The decision log, audit ledger, compatibility contract, examples, Python tools,
   source comments, and deferred Pages text agree with the final implementation.

A P1 or P2 correctness issue cannot close merely because an agent declines to
implement it. It must be fixed, disproven with reproducible evidence and independent
acceptance, or explicitly deferred by the user with a narrowed merge-readiness
claim.

## Non-Negotiable Slow-Work Protocol

### Core Rule

Optimize for correctness, recoverability, and reviewability. Do not optimize this
pass for speed.

The branch touches MPI collectives, restart recovery, on-disk layouts, Python
analysis APIs, and public documentation. An agent that compresses discovery,
implementation, test execution, and approval into one pass is not following this
guide.

### Required Phase Order For Every Checkpoint

Perform each checkpoint in this order:

1. **Orient.** Record branch, `HEAD`, dirty state, relevant prior decisions, and the
   exact files under consideration.
2. **Inspect.** Read the current implementation, tests, examples, compatibility
   contract, and deferred documentation before editing.
3. **Delegate reconnaissance.** Spawn at least one bounded read-only subagent to
   inspect the risk surface independently.
4. **Decide.** Write decision-log entries for every material choice with more than
   one plausible implementation.
5. **Implement narrowly.** Change one coherent behavior family only.
6. **Run focused verification.** Build and test the affected surface before
   expanding scope.
7. **Delegate adversarial audit.** Spawn a different read-only subagent that did not
   implement the change and did not author the pre-edit recommendation.
8. **Correct findings.** Address blocking findings one by one. Do not batch
   unrelated cleanup into the correction.
9. **Re-run verification.** Repeat focused checks after every correction.
10. **Delegate re-audit.** Ask an independent agent to confirm that the correction
    closes the finding without introducing a new regression.
11. **Reflect.** Record whether the checkpoint exposed a flaw in the plan, an
    incorrect assumption, a missing test family, or a need to reorder later work.
12. **Commit coherently.** Only after the stop gate closes, produce a reviewable
    commit or commit series.

If any phase is skipped, the checkpoint remains open.

### Prohibited Rush Patterns

Do not:

- edit code before reading the existing tests and compatibility contract;
- treat a green focused test as proof that the design is complete;
- let the implementing agent approve its own checkpoint;
- accept a subagent summary without inspecting its evidence;
- use one broad subagent prompt such as "review everything" in place of bounded
  audits;
- combine restart arithmetic, MPI error handling, Python refactors, and Pages work
  into one commit;
- change a public file layout without a compatibility decision and fixture plan;
- hide an unsupported configuration behind comments or assumptions;
- claim multi-node behavior from multiple ranks on one physical node;
- claim GPU behavior from a CPU execution of a `_gpu` test;
- update the live `gh-pages` branch before the code feature branch merges;
- remove historical records silently; or
- resolve a difficult design question by copying the old Gotham or rejected
  extraction branch implementation verbatim.

### Mandatory Pause Conditions

Stop editing and reassess if any of these occur:

1. A focused correction changes an on-disk schema.
2. A correction affects legacy shared or per-rank behavior.
3. A correction requires touching a subsystem outside the checkpoint file list.
4. A new failure suggests a pre-existing baseline defect in a touched path.
5. A subagent reports disagreement with an accepted decision.
6. Serial and MPI behavior diverge unexpectedly.
7. A test passes only after weakening validation.
8. A new abstraction starts collecting unrelated responsibilities.
9. A docs update requires claims that the test suite does not prove.
10. External qualification reveals behavior that local one-node tests did not model.

At a pause condition:

1. append a reflection entry to `IO_FEATURE_AUDIT_LEDGER.md`;
2. append a new accepted, rejected, pending-evidence, or superseding decision to
   `IO_FEATURE_BRANCH_DECISION_LOG.md`;
3. update this guide if the checkpoint order or scope needs revision; and
4. request an independent design audit before continuing.

### Plan-Amendment Protocol

Reflection points are allowed to change the plan. They are not allowed to erase why
the plan changed.

When evidence requires a different direction:

1. stop the active checkpoint;
2. record the triggering evidence;
3. identify every later checkpoint affected by the change;
4. append a superseding decision rather than editing an old decision silently;
5. update this guide's checkpoint text and live status board;
6. ask a read-only subagent to review the revised sequence;
7. rerun any earlier verification invalidated by the new direction; and
8. resume only after the revised stop gate is explicit.

## Worktree And Branch Discipline

### Freeze The Approved Guide Before Execution

This guide may remain untracked while it is being drafted and reviewed. Do not begin
RCP-00 until the user-approved revision is tracked in a process-only commit or frozen
as an immutable artifact by another explicit user-approved method.

Record:

1. path;
2. commit SHA when tracked;
3. SHA-256 checksum;
4. line count;
5. byte count; and
6. approval date.

If the guide changes during an audit, stop that audit, record the superseding
revision, and restart the audit against one immutable checksum.

At the start of every work session and checkpoint, recompute the guide SHA-256,
line count, and byte count. Compare them with the approved frozen revision.
Stop immediately if they differ. Resume only after the changed guide is
explicitly approved, frozen, and re-audited.

Start every work session and checkpoint with:

```bash
git fetch --prune origin
wc -l -c IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md
shasum -a 256 IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md
git status --short --branch
git rev-parse HEAD origin/main origin/gh-pages
git worktree list --porcelain
git diff --check origin/main...HEAD
```

Rules:

1. Work on `feature/io-output-formats-and-sharding`.
2. Do not reset, discard, or overwrite user changes.
3. Keep implementation workers in separate worktrees with disjoint write scopes.
4. Use temporary build directories outside the repository.
5. Use isolated runtime directories for executable tests.
6. Keep generated Python caches, pytest caches, build products, and output files out
   of commits.
7. Inspect `git status --short` after every build or test batch.
8. Run `git diff --check origin/main...HEAD` after every checkpoint correction.
9. Record whether remote refs were refreshed successfully. If remote-ref refresh
   fails, record the failure and stop RCP-00 unless the user explicitly authorizes
   bounded offline work. Offline work must not close RCP-00, RCP-07, RCP-10, or
   any merge-readiness claim. Refresh refs and revalidate the affected diffs
   before closure.
10. Classify every dirty or untracked path by owner and intended disposition before
    editing. Stop if ownership is ambiguous.
11. Before committing, inspect:

```bash
git diff --stat
git diff --name-status
git diff --check
git status --short
```

12. Before final review, inspect the full branch diff again:

```bash
git diff --stat origin/main...HEAD
git diff --check origin/main...HEAD
```

## Source-Of-Truth Hierarchy

When records disagree, resolve the contradiction before editing. Use this hierarchy:

1. current executable source and current tests;
2. frozen legacy fixtures produced from `origin/main`;
3. `IO_FORMAT_COMPATIBILITY.md`;
4. accepted non-superseded entries in `IO_FEATURE_BRANCH_DECISION_LOG.md`;
5. current checkpoint records in `IO_FEATURE_AUDIT_LEDGER.md`;
6. this robustification guide;
7. the prior integration plan and historical guide;
8. deferred Pages content; and
9. Gotham history or the rejected extraction branch as historical evidence only.

This hierarchy is not permission to ignore documentation drift. If a lower-ranked
record disagrees with a higher-ranked source, correct or supersede the stale record
as part of the same checkpoint.

## Durable Records

### Extend The Existing Decision Log

Continue using `IO_FEATURE_BRANCH_DECISION_LOG.md`. The current log ends at
`D-068`. New hardening decisions should begin at `D-069`.

Add a decision entry whenever:

- two or more plausible designs exist;
- a format guarantee is added, narrowed, or rejected;
- a public Python API is added, removed, split, or retained;
- a practical resource cap becomes configurable or remains fixed;
- a configuration is supported, rejected, or deferred;
- a durability promise is chosen;
- a subagent recommendation is rejected;
- a checkpoint is reordered after reflection;
- a large refactor is intentionally bounded; or
- process artifacts are retained, archived, or removed.

Use this template:

```markdown
### D-XXX: Short Decision Title

| Field | Value |
| --- | --- |
| Status | Accepted, Rejected, Pending evidence, Superseded, or Accepted with deferred optimization |
| Context | What problem or ambiguity forced a choice |
| Options considered | Enumerate the plausible alternatives |
| Decision | State the selected behavior precisely |
| Reason | Explain why this option is preferred |
| Tradeoffs | State what becomes harder, less general, or deferred |
| Compatibility consequence | State whether file layouts, APIs, tests, or docs change |
| Evidence | Files, tests, benchmark, external qualification, or auditor report |
| Reversal path | Explain how a future implementation could choose differently |
| Supersedes | Prior decision ID if applicable |
| Follow-up | Checkpoint and required action |
```

### Extend The Existing Audit Ledger

Continue using `IO_FEATURE_AUDIT_LEDGER.md`. For every checkpoint append:

```markdown
### YYYY-MM-DD: RCP-XX Short Checkpoint Title

| Field | Record |
| --- | --- |
| Snapshot | Branch, integrated HEAD, guide checksum, and dirty-state summary |
| Scope | Files and behavior reviewed |
| Pre-edit inspection | Relevant code, tests, docs, and prior decisions read |
| Decisions | New or superseded D-XXX entries |
| Implementation | Files changed and behavioral summary |
| Focused verification | Commands, exit codes, retained log paths, and exact results |
| Independent audit | Agent or report ID, scope, reviewed SHA, findings, and disposition |
| Corrections | Finding IDs, fixes, and correction commit SHA |
| Re-audit | Agent or report ID, finding-to-disposition mapping, and result; state explicitly when not applicable |
| Reflection | What changed in the plan or what assumption was validated |
| Remaining risk | Explicit unresolved risk |
| Status | Open, blocked, locally closed, or externally qualified |
```

### Record Reflection Points Explicitly

Do not treat reflection as an informal thought. At every `R-*` reflection point,
append:

```markdown
#### Reflection R-X.Y: Short Title

| Field | Record |
| --- | --- |
| Trigger | Checkpoint completion, failed test, external result, or audit disagreement |
| Assumptions rechecked | List the assumptions reviewed |
| Evidence | Tests, diff inspection, benchmark, or subagent report |
| Direction decision | Continue, revise, split, reorder, defer, or stop |
| Plan changes | Checkpoint edits and new decisions required |
| New risks | Newly visible risks |
```

## Subagent Operating Protocol

### Tool Discovery

Before spawning subagents, discover the available multi-agent tools with
`tool_search`. Use the available spawn, message, wait, and close operations. Do not
assume a particular tool name without discovery.

If no independent-agent mechanism is available, stop before implementation. Record
the missing capability as a blocker. Do not simulate independence by asking the
implementing agent to perform an additional self-review. Resume only when independent
threads or subagents can be launched and their report IDs retained.

### Required Roles

Use bounded subagents deliberately:

| Role | Responsibility |
| --- | --- |
| Reconnaissance auditor | Read-only pre-edit inspection of one risk surface |
| Implementation worker | Narrow code or test changes in a disjoint worktree and file set |
| Adversarial auditor | Read-only post-edit review looking for defects, missing tests, and contract drift |
| Re-auditor | Read-only verification that specific corrections close prior findings |
| Scope auditor | Checks that a checkpoint did not import unrelated refactors |
| Documentation auditor | Compares code, examples, compatibility contract, staged Pages text, and navigation |
| Qualification auditor | Reviews CUDA or scheduler-backed evidence and confirms that it proves the claimed topology |
| Final red-team auditor | Performs a fresh whole-branch review without relying on prior closure claims |

### Required Independence

1. The implementing agent cannot close its own checkpoint.
2. The same subagent should not both implement and approve a change.
3. After correcting an auditor finding, use a different auditor where feasible.
4. The main agent must inspect diffs and rerun critical checks even when a subagent
   reports success.
5. If two auditors disagree, record the disagreement and resolve it through evidence,
   not seniority or convenience.

### Required Subagent Report Format

Every subagent report must include:

1. assigned scope;
2. files inspected or edited;
3. commands run;
4. findings ordered by severity;
5. tests or builds run and exact results;
6. unresolved risks;
7. recommended decision-log entries;
8. recommended audit-ledger entry; and
9. explicit confirmation that unrelated files were not modified.

Also retain:

1. auditor identity, agent ID, thread ID, or report ID;
2. reviewed guide checksum;
3. reviewed integrated commit SHA or working-tree checksum;
4. finding IDs and severity;
5. correction commit SHA when applicable;
6. re-auditor identity and finding-to-disposition mapping;
7. exact command exit codes;
8. retained log paths for long runs; and
9. working-tree state before and after validation.

When a clean audit requires no correction, record:

```text
Re-audit: Not applicable; no correction was required.
```

### Coordinator-Only Integration Protocol

The main agent owns the integration worktree. Implementation workers do not merge
their own changes and do not edit coordinator-owned records.

Coordinator-owned files:

- `IO_FEATURE_BRANCH_DECISION_LOG.md`
- `IO_FEATURE_AUDIT_LEDGER.md`
- `IO_FORMAT_COMPATIBILITY.md` unless documentation ownership is delegated
- this robustification guide

For every implementation worker:

1. create or record a dedicated worker branch and worktree;
2. record the worker base SHA;
3. assign a disjoint write set;
4. require a changed-file list and focused test results;
5. require a worker commit or an explicit patch artifact;
6. inspect the worker diff before integration;
7. integrate only through the main worktree;
8. record conflict resolution and helper-file ownership;
9. rerun focused tests on the integrated snapshot;
10. record the integrated commit SHA; and
11. run post-edit audits against the integrated SHA, never a worker worktree.

### Read-Only Audit Prompt Template

```text
Read-only audit. Do not edit files, create commits, switch refs, or clean the tree.

Repository:
  /Users/dbf75/.codex/worktrees/e948/athenak-DF

Target branch:
  feature/io-output-formats-and-sharding

Checkpoint:
  <RCP-XX>

Bounded scope:
  <files, functions, tests, docs, and exact questions>

Required output:
  1. Findings first, ordered by severity.
  2. File and line references for each finding.
  3. Commands run.
  4. Missing tests.
  5. Decision-log entries that should be added or revised.
  6. Unresolved risks and explicit non-findings.

Do not assume prior closure claims are correct. Reconstruct the reasoning from the
current source and tests.
```

### Implementation Worker Prompt Template

```text
Implement only the assigned checkpoint in your separate worktree.

You are not alone in this codebase. Do not revert or overwrite edits made by other
agents. Do not touch files outside your assigned ownership list. If a required
change crosses that boundary, stop and report it.

Target branch baseline:
  feature/io-output-formats-and-sharding

Checkpoint:
  <RCP-XX>

Owned files:
  <explicit file list>

Required behavior:
  <bounded implementation contract>

Before editing:
  - read the relevant source, tests, IO_FORMAT_COMPATIBILITY.md, and prior decisions;
  - record alternative designs that require a decision-log entry.

Before returning:
  - run focused tests;
  - run git diff --check;
  - report changed files, commands, exact results, unresolved risks, and proposed
    decision-log/audit-ledger text.

Do not merge, cherry-pick, publish gh-pages content, or claim checkpoint closure.
```

### Final Red-Team Prompt Template

```text
Perform a fresh read-only red-team review of the complete branch.
Do not trust prior closure summaries. Do not edit files.

Compare:
  feature/io-output-formats-and-sharding
against:
  origin/main

Audit:
  - restart layout arithmetic and routing;
  - MPI error handling and collective participation;
  - publication guarantees and directory handling;
  - cbin support boundaries;
  - diagnostic semantics and ghost zones;
  - sphslice naming and reconstruction;
  - Python public APIs, limits, and malformed-input handling;
  - tests, fixtures, examples, and external qualification evidence;
  - deferred gh-pages staging helper and documentation consistency;
  - scope discipline and process-artifact disposition.

Return findings first, ordered by severity, with file and line references,
commands run, and explicit residual risks.
```

## Finding Register

The fresh review identified the following work items. Do not silently drop any row.
If a row is rejected or deferred, add a decision-log entry explaining why.

| ID | Severity | Finding | Primary checkpoint |
| --- | --- | --- | --- |
| ROB-001 | P1 | Restart writer and reader arithmetic can overflow before 64-bit chunked IO receives the request | RCP-01 |
| ROB-002 | P1 | Branch-added MPI collectives and communicator operations do not consistently check return codes | RCP-02 |
| ROB-003 | P1 qualification | CUDA execution and true multi-node scheduler qualification remain open | RCP-09 |
| ROB-004 | P2 | Serial positioned IO ignores seek failure and does not reject negative tell results | RCP-01 |
| ROB-005 | P2 | Supported `cbin` semantics for lower dimensions and AMR remain unclear while source comments admit uncertainty | RCP-03 |
| ROB-006 | P2 | Publication guarantees differ across formats; directory creation errors are ignored; crash-durability language is not explicit | RCP-02 |
| ROB-007 | P2 | Advertised generic diagnostics lack full analytic coverage and explicit zero-density/projection-bound behavior | RCP-04 |
| ROB-008 | P2 | Derived variables populate active zones while ordinary outputs can request ghost zones | RCP-04 |
| ROB-009 | P2 | `sphslice` radius filename token uses low-precision `%g` and can collide | RCP-04 |
| ROB-010 | P2 scaling | Node-restart manifest and replicated-header validation run independently on every rank | RCP-08 |
| ROB-011 | P3 | Output registration retains duplicated parsing debt and an overgrown inline PDF parser | RCP-05 |
| ROB-012 | P3 | Python reader/converter modules are monolithic and repeat strict-validation helpers; practical caps are fixed internally | RCP-06 |
| ROB-013 | P3 | Deferred Pages content is prepared, but post-merge staging remains manual and drift-prone | RCP-07 |
| ROB-014 | P3 | Large process records need a deliberate merge-packaging decision | RCP-10 |
| ROB-015 | P3 | Touched APIs and source comments need cleanup so guarantees and boundaries are explicit | All, close in RCP-10 |
| ROB-016 | P3 | Restart startup logs a missing file and then continues; touched-path cleanup should fail clearly and remove stale local comments | RCP-01 |
| ROB-017 | P1 | Five-digit output-number buffers truncate dump numbers at `100000`, allowing deterministic filename reuse and overwrite | RCP-02 |
| ROB-018 | P2 | `cbin` Kokkos launch and normalization ranges still multiply iteration counts in `int` | RCP-03 |
| ROB-019 | P2 | Strict node-restart manifest parsing caps record counts but not individual line or total file size | RCP-01 |
| ROB-020 | P2 | The MPI inventory must include branch-added and inherited-but-touched MPI calls, including `.bin`, `.cbin`, and base output paths | RCP-02 |
| ROB-021 | P2 | Multiple configured outputs can resolve to one deterministic public namespace, including modern PDF blocks with the same explicit `id` | RCP-02 |
| ROB-022 | P3 | Writer-side PDF and `sphslice` geometry allocations need deliberate preflight limits and diagnostics | RCP-04 |

## Starting Evidence Map

These locations are starting points, not exhaustive proof. Re-run searches and read
the surrounding implementation before editing.

| Finding | Starting source evidence |
| --- | --- |
| ROB-001 | `src/outputs/restart.cpp:354`, `src/outputs/restart.cpp:423`, `src/pgen/pgen.cpp:235`, `src/pgen/pgen.cpp:332` |
| ROB-002 | `src/globals.cpp:39`, `src/globals.cpp:42`, `src/driver/driver.cpp:37`, `src/outputs/pdf.cpp:278`, `src/outputs/restart.cpp:245`, `src/outputs/restart.cpp:754` |
| ROB-003 | `IO_FEATURE_AUDIT_LEDGER.md:509`, `IO_FORMAT_COMPATIBILITY.md:43`, `tst/test_suite/io/test_node_sharding_mpicpu.py:475` |
| ROB-004 | `src/outputs/io_wrapper.cpp:617`, `src/outputs/io_wrapper.cpp:725`, `src/outputs/io_wrapper.cpp:818` |
| ROB-005 | `src/outputs/outputs.cpp:79`, `src/outputs/coarsened_binary.cpp:575`, `src/outputs/coarsened_binary.cpp:587` |
| ROB-006 | `src/outputs/binary.cpp:95`, `src/outputs/binary.cpp:140`, `src/outputs/coarsened_binary.cpp:120`, `src/outputs/coarsened_binary.cpp:409`, `src/outputs/pdf.cpp:120`, `src/outputs/pdf.cpp:377`, `src/outputs/spherical_slice.cpp:331`, `src/outputs/spherical_slice.cpp:484`, `src/outputs/restart.cpp:88`, `src/outputs/restart.cpp:722` |
| ROB-007 | `src/outputs/derived_variables.cpp:1296`, `src/outputs/derived_variables.cpp:1340`, `src/outputs/derived_variables.cpp:1404`, `deferred_docs/gh-pages/io-output-formats-and-sharding/overlay/docs/source/modules/outputs.md:129` |
| ROB-008 | `src/outputs/basetype_output.cpp:862`, `src/outputs/basetype_output.cpp:944`, `src/outputs/derived_variables.cpp:1282` |
| ROB-009 | `src/outputs/spherical_slice.cpp:472` |
| ROB-010 | `src/restart_manifest.cpp:248`, `src/restart_manifest.cpp:360`, `src/main.cpp:269`, `IO_FEATURE_BRANCH_DECISION_LOG.md:364` |
| ROB-011 | `src/outputs/outputs.cpp:129`, `src/outputs/outputs.cpp:208`, `src/outputs/outputs.cpp:226`, `src/outputs/outputs.cpp:301` |
| ROB-012 | `vis/python/bin_convert.py:95`, `vis/python/read_pdf.py:29`, `vis/python/read_sphslice.py:21`, `tst/test_suite/io/test_python_io_readers_cpu.py` |
| ROB-013 | `deferred_docs/gh-pages/io-output-formats-and-sharding/MANIFEST.md:41`, `deferred_docs/gh-pages/io-output-formats-and-sharding/VALIDATION.md:29` |
| ROB-014 | `IO_FEATURE_BRANCH_GUIDE.md`, `IO_FEATURE_BRANCH_INTEGRATION_PLAN.md`, `IO_FEATURE_BRANCH_DECISION_LOG.md`, `IO_FEATURE_AUDIT_LEDGER.md`, `IO_FORMAT_COMPATIBILITY.md` |
| ROB-015 | `src/outputs/io_wrapper.hpp:42`, `src/globals.hpp:25`, `src/restart_manifest.hpp:37`, `src/outputs/outputs.cpp:18` |
| ROB-016 | `src/main.cpp:291` |
| ROB-017 | `src/outputs/outputs.cpp:125`, `src/outputs/binary.cpp:120`, `src/outputs/coarsened_binary.cpp:382`, `src/outputs/pdf.cpp:335`, `src/outputs/pdf.cpp:443`, `src/outputs/restart.cpp:224`, `src/outputs/spherical_slice.cpp:472`, `vis/python/read_pdf.py:22` |
| ROB-018 | `src/outputs/coarsened_binary.cpp:285`, `src/outputs/coarsened_binary.cpp:294`, `src/outputs/coarsened_binary.cpp:302`, `src/outputs/coarsened_binary.cpp:334` |
| ROB-019 | `src/restart_manifest.cpp:36`, `src/restart_manifest.cpp:161`, `src/restart_manifest.cpp:345`, `src/restart_manifest.cpp:385`, `src/restart_manifest.cpp:406` |
| ROB-020 | `src/outputs/binary.cpp:222`, `src/outputs/coarsened_binary.cpp:216`, `src/outputs/coarsened_binary.cpp:505`, `src/outputs/basetype_output.cpp:924` |
| ROB-021 | `src/outputs/pdf.cpp:93`, `src/outputs/pdf.cpp:378`, `src/outputs/pdf.cpp:445`, `src/outputs/outputs.cpp:97` |
| ROB-022 | `src/outputs/spherical_slice.cpp:75`, `src/outputs/spherical_slice.cpp:178`, `src/outputs/outputs.cpp:457` |

## Checkpoint Overview

| Checkpoint | Purpose | Required stop gate |
| --- | --- | --- |
| RCP-00 | Reopen the baseline carefully and freeze current facts | Baseline record, issue register, and audit assignments are accepted |
| RCP-01 | Complete restart arithmetic and serial positioned-IO safety | Overflow, serial seek/tell, resume, and startup failure tests pass after independent re-audit |
| RCP-02 | Centralize MPI failure handling and filesystem publication primitives | All branch-added MPI calls and publication paths pass independent review |
| RCP-03 | Settle and enforce the `cbin` support contract | AMR and lower-dimensional behavior are either qualified or rejected explicitly |
| RCP-04 | Harden diagnostic semantics, derived ghost zones, and `sphslice` naming | Every advertised diagnostic and naming boundary has automated evidence |
| RCP-05 | Refactor output registration and touched C++ interfaces | Behavior-preserving cleanup passes regression and scope audits |
| RCP-06 | Refactor Python tooling while preserving one public surface | API, malformed-input, budget, fixture, and CLI tests pass |
| RCP-07 | Automate deferred Pages staging and align documentation | Detached refreshed Pages staging and warnings-as-errors build pass |
| RCP-08 | Prepare manifest-scaling instrumentation, collect evidence during RCP-09, and optimize only if justified | Benchmark evidence and a documented keep-or-refactor decision exist |
| RCP-09 | Execute external CUDA and multi-node qualification, including RCP-08 measurements | Real environment evidence closes both external gates |
| RCP-10 | Package process records and run final whole-branch red-team review | Final independent audits, full regression matrix, and documentation agreement pass |

## Live Checkpoint Board

Copy this table into the top of the new robustification section in
`IO_FEATURE_AUDIT_LEDGER.md` and update it incrementally. Do not mark several rows
closed retroactively at the end of the pass.

| Checkpoint | Status | Decisions | Pre-edit auditors | Implementation commit(s) | Focused tests | Post-edit auditors | Reflection | Remaining risk |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RCP-00 | Not started |  |  |  |  |  |  |  |
| RCP-01 | Not started |  |  |  |  |  |  |  |
| RCP-02 | Not started |  |  |  |  |  |  |  |
| RCP-03 | Not started |  |  |  |  |  |  |  |
| RCP-04 | Not started |  |  |  |  |  |  |  |
| RCP-05 | Not started |  |  |  |  |  |  |  |
| RCP-06 | Not started |  |  |  |  |  |  |  |
| RCP-07 | Not started |  |  |  |  |  |  |  |
| RCP-08 | Not started |  |  |  |  |  |  |  |
| RCP-09 | Not started |  |  |  |  |  |  |  |
| RCP-10 | Not started |  |  |  |  |  |  |  |

Allowed statuses:

- `Not started`
- `Reconnaissance`
- `Decision pending`
- `Implementing`
- `Focused verification`
- `Independent audit`
- `Correction required`
- `Re-audit`
- `Reflection required`
- `Locally closed`
- `Externally qualified`
- `Blocked`

Rules:

1. Update the board when the phase changes.
2. Link decision IDs and audit-ledger headings.
3. Never use `Locally closed` when an independent re-audit is still pending.
4. Never use `Externally qualified` for a CPU-only GPU smoke test or a one-node MPI
   run.
5. Never mark `RCP-10` closed while another row remains open.

## Seeded Decision Queue

These IDs are reservations for the first expected hardening decisions. Append them
in evidence order. If a new decision appears earlier, use the next available ID and
keep the queue current rather than forcing an inaccurate ordering.

| Suggested ID | Checkpoint | Decision question |
| --- | --- | --- |
| D-069 | RCP-00 | Accept, revise, or expand the fresh `ROB-*` register and checkpoint sequence |
| D-070 | RCP-01 | Use a shared restart-layout descriptor, mirrored checked helpers, or narrower arithmetic utilities |
| D-071 | RCP-01 | Select the portable serial large-file seek/tell implementation |
| D-072 | RCP-02 | Select the shared MPI error-reporting helper boundary and teardown policy |
| D-073 | RCP-02 | Define namespace atomicity versus crash durability for each output family |
| D-074 | RCP-02 | Define checked directory creation, stale-temporary behavior, and generation collision handling |
| D-075 | RCP-03 | Define the complete supported and rejected `cbin` matrix for dimensions, ghost zones, slices, and AMR |
| D-076 | RCP-04 | Define diagnostic projection clamps, zero-density behavior, and axis-origin behavior |
| D-077 | RCP-04 | Populate derived ghost zones or reject unsupported derived `ghost_zones=true` output |
| D-078 | RCP-04 | Select collision-resistant `sphslice` radius naming |
| D-079 | RCP-05 | Bound the C++ parser and interface refactor |
| D-080 | RCP-06 | Extract private Python reader utilities or retain local helpers |
| D-081 | RCP-06 | Select the public configuration mechanism for practical reader limits |
| D-082 | RCP-07 | Select Pages staging language, anchor policy, idempotence markers, and drift handling |
| D-083 | RCP-08 | Keep replicated manifest validation, centralize it, use node leaders, or split a scaling branch |
| D-084 | RCP-10 | Retain, archive, summarize, or move process artifacts to pull-request history |
| D-085 | RCP-01 | Select bounded node-manifest line and total-file limits |
| D-086 | RCP-02 | Preserve five-digit-minimum output numbering with unbounded rendering or reject an explicit maximum |
| D-087 | RCP-02 | Define construction-time output namespace collision detection and compatibility boundaries |
| D-088 | RCP-03 | Select checked `cbin` kernel-range and normalization-range representation |
| D-089 | RCP-04 | Select writer-side PDF and `sphslice` allocation preflight policy |
| D-090 | RCP-08 | Preregister production-scale benchmark topology and the threshold that triggers optimization |

For every row:

1. record options before implementation;
2. record evidence and tradeoffs;
3. state the reversal path; and
4. update the queue if reflection changes the required choice.

## RCP-00: Reopen The Baseline And Freeze Facts

### Purpose

Prevent the new hardening pass from inheriting unverified assumptions from the prior
closure record. Re-establish exactly what exists before editing.

### Required Local Inspection

Run:

```bash
git status --short --branch
git rev-parse HEAD origin/main origin/gh-pages
git merge-base HEAD origin/main
git diff --stat origin/main...HEAD
git diff --check origin/main...HEAD
git worktree list --porcelain
```

Read:

- `IO_FEATURE_BRANCH_DECISION_LOG.md`;
- `IO_FEATURE_AUDIT_LEDGER.md`;
- `IO_FORMAT_COMPATIBILITY.md`;
- this guide;
- `src/outputs/restart.cpp`;
- `src/pgen/pgen.cpp`;
- `src/outputs/io_wrapper.*`;
- `src/globals.*`;
- `src/outputs/outputs.*`;
- `src/outputs/coarsened_binary.cpp`;
- `src/outputs/derived_variables.cpp`;
- `src/outputs/spherical_slice.cpp`;
- `vis/python/bin_convert.py`;
- `vis/python/read_pdf.py`;
- `vis/python/read_sphslice.py`;
- the complete `tst/test_suite/io/` inventory; and
- `deferred_docs/gh-pages/io-output-formats-and-sharding/`.

### Required Subagents

Spawn three read-only auditors in parallel:

1. **Runtime baseline auditor:** restart arithmetic, `IOWrapper`, MPI calls, and
   publication.
2. **Format/tooling baseline auditor:** `cbin`, PDF, `sphslice`, Python APIs, reader
   caps, examples, and fixtures.
3. **Documentation/process baseline auditor:** compatibility contract, deferred
   Pages package, audit ledger, decision log, and process-artifact inventory.

Do not edit while these auditors run. Use the waiting period to inspect the same
files locally and compare your reasoning after reports arrive.

### Required Decisions

Append new decisions beginning at `D-069` for:

1. the locally actionable issue register;
2. the checkpoint sequence;
3. whether publication durability means atomic namespace publication or
   crash-durable persistence;
4. whether `cbin` should support AMR and lower-dimensional meshes in this branch or
   reject them explicitly;
5. whether Python budget overrides belong in function arguments, CLI arguments,
   environment variables, a configuration object, or a combination; and
6. whether manifest-validation scaling is a merge blocker or evidence-driven
   optimization checkpoint.

### Reflection Point R-00

Ask:

1. Did any auditor find a correctness defect not present in this guide?
2. Are any planned refactors too broad for the feature branch?
3. Does the checkpoint order minimize risk, or should a format-contract choice move
   earlier?
4. Are any existing prior "resolved" findings invalidated by the fresh audit?

### Stop Gate

Do not edit implementation code until:

- the branch snapshot is recorded;
- auditor reports are recorded;
- new findings are added to the register;
- initial D-069+ decisions are appended; and
- the main agent has explicitly accepted or rejected each auditor recommendation.

## RCP-01: Restart Arithmetic, Serial Positioned IO, And Startup Failure

### Purpose

Make the large-checkpoint safety promise true from layout calculation through
physical IO, not merely inside `IOWrapper`.

### Primary Files

- `src/outputs/restart.cpp`
- `src/pgen/pgen.cpp`
- `src/outputs/io_wrapper.cpp`
- `src/outputs/io_wrapper.hpp`
- `src/main.cpp`
- `src/restart_manifest.cpp`
- `src/restart_manifest.hpp`
- focused IO wrapper and restart tests under `tst/test_suite/io/`

### Problem Details

The wrapper uses `IOWrapperSizeT = std::uint64_t` and chunks MPI byte transfers, but
restart callers still use expressions such as:

```cpp
nout1*nout2*nout3*nhydro*sizeof(Real)
```

where `nout1`, `nout2`, `nout3`, and component counts are `int`. The product can
overflow before assignment to `IOWrapperSizeT`.

Restart paths also narrow Kokkos view sizes:

```cpp
int mbcnt = mbptr.size();
```

A single large view can therefore overflow before entering the chunking layer.

Serial positioned IO calls `std::fseek` without checking its result, and
`GetPosition()` returns `ftell()` results without rejecting negative values.

The touched restart-startup path logs a missing restart file but continues into the
open path, producing delayed and less useful failure behavior.

The node-manifest parser bounds payload and segment counts, but still reads signature,
scalar, payload, segment, and trailing lines into unconstrained strings. A malformed
manifest can force large allocations before semantic validation.

### Pre-Edit Reconnaissance

Map every restart layout component:

1. replicated parameter/header bytes;
2. MeshBlock metadata bytes;
3. problem-generator restart bytes;
4. Hydro cell-centered bytes;
5. MHD cell-centered bytes;
6. MHD face-centered bytes;
7. radiation bytes;
8. forcing bytes;
9. Z4c bytes;
10. ADM bytes;
11. per-MeshBlock bytes;
12. per-rank offsets;
13. per-node offsets;
14. manifest-declared payload bytes; and
15. reader virtual shared-file offsets.

Do not implement until the writer and reader calculations are compared field by
field in a written audit note.

### Required Design Decision

Choose deliberately between:

1. a shared restart-layout descriptor used by writer and reader;
2. mirrored checked helper functions with an explicit parity test; or
3. a smaller checked-arithmetic utility applied separately to existing code.

Prefer the shared layout descriptor if it reduces duplicated arithmetic without
forcing an unrelated restart rewrite. Record why the chosen level of abstraction
fits the codebase.

Also decide whether serial large-file support uses:

- `fseeko`/`ftello`;
- a platform abstraction already present in the codebase; or
- another portable checked wrapper.

Record portability assumptions.

Select and record production-appropriate limits for:

1. total manifest bytes;
2. signature line bytes;
3. scalar-record line bytes;
4. payload-record line bytes;
5. segment-record line bytes;
6. trailing-record line bytes; and
7. generated payload path bytes.

Limits must accommodate expected production topology while rejecting malformed
inputs before unbounded allocation.

### Implementation Requirements

1. Use checked multiply and add operations for every restart size and offset.
2. Promote operands before multiplication.
3. Never narrow Kokkos view sizes to `int` before IO.
4. Pass `IOWrapperSizeT` counts through writer and reader helpers.
5. Check serial seek failures before reading or writing.
6. Reject negative tell results explicitly.
7. Keep MPI offset-range checks in the wrapper.
8. Make missing restart files fail immediately with a clear message.
9. Remove stale local comments and unused touched-path variables.
10. Preserve legacy shared and per-rank restart bytes unless a documented
    compatibility decision says otherwise.
11. Preserve node manifest-only entry semantics and direct span routing.
12. Bound individual node-manifest lines and total manifest bytes before
    unconstrained string growth.
13. Report which manifest record exceeded its limit.

### Required Tests

Add or extend:

1. checked-layout unit coverage with artificial extents above `INT_MAX`;
2. multiplication overflow negatives;
3. addition overflow negatives;
4. single-view count above `INT_MAX` behavior without allocating enormous memory;
5. serial seek failure behavior;
6. serial negative-position handling through a focused harness or injected wrapper;
7. frozen shared restart resume;
8. frozen per-rank restart resume;
9. generated shared restart resume;
10. generated per-rank restart resume;
11. generated node-manifest resume;
12. forced-small-chunk node restart write and read;
13. changed-rank-count node resume on the local node;
14. malformed manifest byte-count rejection; and
15. oversized manifest signature, scalar, payload, segment, trailing line, and total
    file rejection; and
16. missing restart path failure.

### Required Subagents

1. Spawn a pre-edit restart-layout auditor.
2. Assign implementation only after the layout map and design decision exist.
3. Spawn a different post-edit arithmetic auditor.
4. Spawn a separate compatibility auditor to compare frozen shared/per-rank behavior.
5. After corrections, spawn a re-auditor focused only on the previously reported
   defects.

### Reflection Point R-01

Ask:

1. Did the shared layout work reveal duplicated logic elsewhere?
2. Is the abstraction narrow enough?
3. Did any legacy fixture change unexpectedly?
4. Are wrapper-level chunk tests still representative after caller changes?
5. Should any new helper move into a general checked-IO utility for RCP-02?

### Stop Gate

Do not proceed until caller-level large-count safety, serial positioned IO, legacy
resume, node resume, bounded manifest parsing, and startup failure behavior are
independently accepted.

## RCP-02: MPI Failure Handling And Filesystem Publication

### Purpose

Make branch-added MPI and filesystem behavior fail predictably rather than silently.

### Primary Files

- `src/globals.cpp`
- `src/globals.hpp`
- `src/driver/driver.cpp`
- `src/outputs/io_wrapper.cpp`
- `src/outputs/io_wrapper.hpp`
- `src/outputs/binary.cpp`
- `src/outputs/coarsened_binary.cpp`
- `src/outputs/pdf.cpp`
- `src/outputs/spherical_slice.cpp`
- `src/outputs/restart.cpp`
- `src/file_sharding.hpp`
- relevant tests under `tst/test_suite/io/`

### Pre-Edit MPI Inventory

Produce an explicit inventory of every branch-added and inherited-but-touched MPI
call. Start from the branch diff rather than a hand-maintained subset:

```bash
git diff --name-only origin/main...HEAD -- src |
  while IFS= read -r file; do
    rg -n "MPI_[A-Za-z0-9_]+\\(" "$file" || true
  done
```

At minimum, verify `src/globals.cpp`, `src/driver/driver.cpp`,
`src/outputs/io_wrapper.cpp`, `src/outputs/binary.cpp`,
`src/outputs/coarsened_binary.cpp`, `src/outputs/basetype_output.cpp`,
`src/outputs/pdf.cpp`, `src/outputs/spherical_slice.cpp`, and
`src/outputs/restart.cpp`.

Classify every result as:

1. branch-added and in scope;
2. inherited but touched by this branch and requiring disposition;
3. inherited and outside the checkpoint scope; or
4. fatal-shutdown handling that must remain simple.

Classify each call:

| Category | Examples | Required behavior |
| --- | --- | --- |
| Communicator lifecycle | `MPI_Comm_split_type`, `MPI_Comm_free` | Check and report failure |
| Node metadata | `MPI_Exscan`, `MPI_Bcast`, `MPI_Allreduce` | Check and report failure |
| Output reductions | PDF and timing `MPI_Reduce` | Check and report failure |
| Publication synchronization | restart `MPI_Barrier`, `MPI_Gather` | Check and report failure |
| File IO | `MPI_File_*` | Preserve wrapper checks and improve consistency |
| Fatal shutdown | `MPI_Abort` | Keep simple and avoid recursive reporting failure |

### Required Design Decisions

Record:

1. whether to expose one common `CheckMpi` helper in a shared header or keep a
   narrowly scoped utility;
2. how MPI errors are rendered using `MPI_Error_string`;
3. whether communicator cleanup failure is fatal during normal teardown;
4. whether publication should be:
   - namespace-atomic only through temporary file plus rename; or
   - crash-durable through file sync and containing-directory sync;
5. whether `.bin` and `.cbin` should adopt atomic rename publication now;
6. how stale temporary files are handled;
7. how output directories are created and validated; and
8. whether restart generation names need stronger collision resistance or
   exclusive-create behavior;
9. whether output sequence tokens preserve five-digit minimum width while rendering
   larger values without truncation, or whether an explicit maximum is rejected
   before publication;
10. how Python readers discover widened sequence tokens while retaining old names;
    and
11. whether construction builds a deterministic namespace map that rejects two
    output blocks resolving to the same public targets.

### Implementation Requirements

1. Check every branch-added MPI return code.
2. Give errors enough context to identify operation, file family, communicator
   scope, and checkpoint stage.
3. Add a shared checked directory helper:
   - accept an existing directory;
   - reject a conflicting non-directory path;
   - fail on permission errors;
   - retain existing permission choices unless deliberately changed.
4. Add a reusable publication helper where it reduces duplication.
5. Keep temporary files private until complete.
6. Clean up failed temporary files where practical.
7. Document whether a rename is merely namespace-atomic or also crash-durable.
8. If crash durability is promised, sync payloads, manifests, and containing
   directories in the correct order.
9. Reconcile `.bin`, `.cbin`, modern PDF, `sphslice`, and node restart guarantees.
10. Do not silently change legacy PDF behavior without a compatibility decision.
11. Eliminate fixed-buffer truncation for output numbers at and above `100000`, or
    reject an explicitly documented maximum before any write.
12. Preserve existing filenames below `100000`.
13. Audit Python reader filename matching for widened output numbers.
14. Detect deterministic target namespace collisions during construction, including
    modern PDF blocks with the same explicit `id` but different first variables,
    weights, bins, or scales.
15. Audit collision behavior for `.bin`, `.cbin`, modern PDF, `sphslice`, and
    restart output while preserving legacy behavior unless a decision authorizes a
    change.

### Required Tests

Add or extend:

1. conflicting output-directory path rejection;
2. unwritable output-directory rejection where portable;
3. stale temporary file handling;
4. `.bin` temporary publication and cleanup if adopted;
5. `.cbin` temporary publication and cleanup if adopted;
6. PDF temporary publication and cleanup regression;
7. `sphslice` temporary publication and cleanup regression;
8. restart payload-before-manifest publication regression;
9. restart manifest rejection for unpublished payloads;
10. wrapper harness coverage for MPI collective count mismatch and forced chunks;
11. output timing MPI label regression; and
12. a diff-driven source audit gate proving branch-added and inherited-but-touched
    MPI calls are dispositioned;
13. output numbers `99999`, `100000`, and `100001` across `.bin`, `.cbin`, PDF,
    `sphslice`, and restart output;
14. Python reader compatibility with the widened sequence-token contract if
    adopted; and
15. duplicate target namespace rejection for PDF blocks that differ by variables,
    weights, bins, or scales and for other deterministic collisions found by audit.

### Required Subagents

1. Spawn a pre-edit MPI auditor and a separate filesystem-publication auditor.
2. If code changes split cleanly, assign disjoint implementation workers:
   - worker A: common MPI helpers and node communicator paths;
   - worker B: publication and directory helpers;
   - worker C: focused tests.
3. Spawn a post-edit collective-participation auditor.
4. Spawn a post-edit filesystem-semantics auditor.
5. Re-audit corrections with agents that did not implement them.

### Reflection Point R-02

Ask:

1. Did a shared helper reduce duplication without obscuring control flow?
2. Did publication hardening change any filenames or reader assumptions?
3. Is crash durability actually required, or should documentation narrow the claim?
4. Are temporary-file semantics consistent enough for users and tooling?
5. Does any format still fail differently for rank and node shards without reason?

### Stop Gate

Do not proceed until the MPI inventory is fully dispositioned, publication semantics
are explicit, directory errors fail clearly, focused regressions pass, and two
independent audits accept the result.

## RCP-03: Settle The Coarsened-Binary Contract

### Purpose

Turn `cbin` from an implementation with historical uncertainty into an explicit,
tested format contract.

### Primary Files

- `src/outputs/coarsened_binary.cpp`
- `src/outputs/outputs.cpp`
- `src/outputs/outputs.hpp`
- `vis/python/bin_convert.py`
- `inputs/io/`
- `tst/inputs/`
- `tst/test_suite/io/`
- `IO_FORMAT_COMPATIBILITY.md`
- deferred Pages overlays and insertions

### Questions That Must Be Answered

1. Is full-volume `cbin` supported on AMR meshes?
2. If yes, how are reduced-grid logical locations and refinement levels represented?
3. If no, how is AMR rejected before output publication?
4. Is `cbin` supported for 1D and 2D meshes?
5. Should singleton axes remain singleton instead of participating in the
   `coarsen_factor` lower bound?
6. Are ghost-expanded extents supported when divisible by the factor?
7. Is sliced `cbin` still deliberately rejected in every shard mode?
8. Do shared, rank, and node files use identical per-record semantics?
9. Are reader reconstruction and ATHDF/XDMF conversion correct for every supported
   case?

### Required Design Decisions

Record one explicit support matrix. Replace every `Pending D-075 decision` status
before implementation:

| Configuration | Initial status | Evidence required | Decision ID | Test owner |
| --- | --- | --- | --- | --- |
| Uniform 3D full-volume shared | Required supported path | Writer, direct read, ATHDF/XDMF conversion, malformed metadata rejection | D-075 | RCP-03 |
| Uniform 3D full-volume rank | Required supported path | Writer, shard assembly, numerical equality, conversion | D-075 | RCP-03 |
| Uniform 3D full-volume node | Required supported path | Writer, explicit-empty inventory, shard assembly, numerical equality, conversion | D-075 | RCP-03 and RCP-09 |
| Uniform 1D full-volume | Pending D-075 decision | Positive writer/read/conversion evidence or construction-time rejection | D-075 | RCP-03 |
| Uniform 2D full-volume | Pending D-075 decision | Positive writer/read/conversion evidence or construction-time rejection | D-075 | RCP-03 |
| Ghost-zone full-volume | Pending D-075 decision | Divisibility, placement, readback, and conversion evidence or construction-time rejection | D-075 | RCP-03 |
| Sliced output | Required rejected path unless redesigned | Construction-time rejection and no publication | D-075 | RCP-03 |
| AMR full-volume | Pending D-075 decision | Correct reduced-grid location/refinement evidence or construction-time rejection | D-075 | RCP-03 |
| AMR sliced output | Required rejected path unless redesigned | Construction-time rejection and no publication | D-075 | RCP-03 |

### Implementation Requirements

1. Remove uncertainty comments only after replacing them with enforced behavior and
   accurate explanatory comments.
2. Validate unsupported modes during construction, before data loading.
3. Keep reader validation aligned with writer guarantees.
4. Avoid changing legacy supported bytes unless documented and tested.
5. Add user-facing error messages that name the unsupported configuration and the
   supported alternative.
6. Update examples and deferred docs to state the exact matrix.
7. Use checked arithmetic for coarsen cubes, total elements, Kokkos launch ranges,
   indexing ranges, and moment-normalization ranges.
8. Reject a range that cannot be represented before launching a Kokkos kernel.

### Required Tests

For every supported row:

1. writer execution;
2. canonical `bin_convert.py` direct read;
3. shard assembly where applicable;
4. ATHDF/XDMF conversion where supported;
5. shared/rank/node numerical equality; and
6. malformed metadata rejection.

For every rejected row:

1. construction-time failure;
2. no public file publication; and
3. stable diagnostic message.

Also add representability-boundary tests for checked `cbin` products and Kokkos
launch ranges without allocating production-scale arrays.

### Required Subagents

1. Spawn an AMR-aware format auditor before design.
2. Spawn a reader/writer parity auditor after implementation.
3. Spawn a docs auditor after the support matrix changes.

### Reflection Point R-03

Ask:

1. Is supporting AMR now worth the additional format complexity?
2. Would an explicit uniform-grid-only contract be more honest and maintainable?
3. Did lower-dimensional support expose assumptions shared with `.bin`?
4. Did any docs claim exceed the implemented matrix?

### Stop Gate

Do not proceed until `cbin` has a precise support matrix, construction-time guards,
checked Kokkos range arithmetic, positive and negative tests, reader parity, and
matching documentation.

## RCP-04: Diagnostic Semantics, Ghost Zones, And Spherical-Slice Naming

### Purpose

Qualify the complete user-visible diagnostic surface and eliminate ambiguous output
behavior.

### Primary Files

- `src/outputs/basetype_output.cpp`
- `src/outputs/derived_variables.cpp`
- `src/outputs/spherical_slice.cpp`
- `src/outputs/outputs.hpp`
- `inputs/io/`
- `tst/inputs/`
- `tst/test_suite/io/test_output_formats_cpu.py`
- `tst/test_suite/io/test_output_formats_mpicpu.py`
- `tst/test_suite/io/test_output_formats_gpu.py`
- documentation and compatibility records

### Advertised Diagnostic Families

Audit every name:

| Family | Names |
| --- | --- |
| Cartesian coordinates | `coord_x`, `coord_y`, `coord_z` |
| Spherical coordinates | `coord_r`, `coord_theta`, `coord_phi`, `coord_costheta`, `coord_abscostheta` |
| Cylindrical coordinates | `coord_cyl_R`, `coord_cyl_phi`, `coord_cyl_z` |
| Velocity projections | `vel_sph_r`, `vel_sph_theta`, `vel_sph_phi`, `vel_cyl_R`, `vel_cyl_phi` |
| Spherical mass flux | `mdot_sph`, `mdot_sph_out`, `mdot_sph_in` |
| Vertical mass flux | `mdot_vert`, `mdot_vert_out`, `mdot_vert_in` |
| Spherical energy flux | `edot_sph`, `edot_sph_out`, `edot_sph_in`, `edot_sph_kin`, `edot_sph_th`, `edot_sph_mag` |
| Vertical energy flux | `edot_vert`, `edot_vert_out`, `edot_vert_in` |
| Passive scalars | `hydro_u_s_N`, `hydro_w_s_N`, `mhd_u_s_N`, `mhd_w_s_N` |

### Required Design Decisions

Record:

1. whether `z/r` is clamped to `[-1, 1]` before `acos` and cosine outputs;
2. the zero-density policy for velocity and energy diagnostics:
   - use a physical floor already established by the fluid module;
   - emit a defined sentinel;
   - fail at runtime when a non-finite diagnostic is detected;
   - reject only statically incompatible module combinations during construction;
     or
   - another explicitly justified runtime policy;
3. behavior at `r = 0` and cylindrical radius `R = 0`;
4. whether derived variables support `ghost_zones=true`;
5. if not, whether unsupported combinations fail during construction;
6. whether a future ghost-zone-safe derived-kernel design belongs in this branch;
7. a stable, collision-resistant `sphslice` radius token format; and
8. whether two configured spherical slices may share `id` safely;
9. writer-side preflight limits for PDF total bins and retained arrays; and
10. writer-side preflight limits for `sphslice` angular geometry, interpolation
    arrays, sparse buffers, and dense output buffers.

### Implementation Requirements

1. Clamp floating-point projection inputs where mathematically bounded.
2. Define zero-density behavior explicitly.
3. Preserve the existing two-fluid rejection for ambiguous generic fluid
   diagnostics.
4. Preserve the `sphslice` rejection of derived arrays until ghost-zone-safe
   sampling exists.
5. Either populate derived ghost zones correctly or reject unsupported derived
   `ghost_zones=true` combinations before data loading.
6. Replace `%g` radius naming with deterministic high-precision naming.
7. Test collision resistance for nearby radii.
8. Keep existing `file_type=sph` distinct from `file_type=sphslice`.
9. Estimate writer-side PDF and `sphslice` allocation bytes before allocation.
10. Reject unreasonable writer geometry with clear diagnostics.
11. Make writer limits configurable or document why a fixed limit is appropriate.

### Required Analytic Tests

Construct deterministic uniform or simple-gradient states and test:

1. every coordinate diagnostic;
2. all spherical and cylindrical velocity projections;
3. inflow, outflow, and signed mass flux;
4. Hydro total, kinetic, and thermal energy flux;
5. MHD magnetic and total energy flux;
6. vertical sign handling above and below the midplane;
7. behavior at `r = 0`;
8. behavior at cylindrical `R = 0`;
9. zero-density policy;
10. floating-point clamp behavior;
11. passive scalar indexing;
12. nonideal EOS rejection where required;
13. ion-neutral two-fluid rejection;
14. derived ghost-zone acceptance or rejection;
15. nearby `sphslice` radii producing distinct filenames;
16. shared/rank/node `sphslice` numerical equality; and
17. writer-side PDF and `sphslice` capacity rejection without uncontrolled
    allocation; and
18. CUDA execution of representative derived-variable PDF axes.

### Required Subagents

1. Spawn a numerical-semantics auditor before editing.
2. Spawn a test-design auditor to propose analytic oracles independently.
3. Use a separate post-edit GPU-readiness auditor.
4. Spawn a docs-to-code auditor after final semantics are chosen.

### Reflection Point R-04

Ask:

1. Are these diagnostics general output fields or PDF-focused fields?
2. Should any advertised diagnostic be removed until its semantics are strong?
3. Are ghost-zone restrictions consistent with existing derived-variable behavior?
4. Does high-precision radius naming preserve usability?
5. Does CUDA qualification require a different test state?

### Stop Gate

Do not proceed until every advertised diagnostic has either positive analytic
coverage or a documented rejection boundary, `sphslice` naming cannot collide under
the tested precision cases, writer-side PDF and `sphslice` capacity is preflighted
before allocation, and independent auditors accept the numerical contract.

## RCP-05: Refactor Output Registration And Touched C++ Interfaces

### Purpose

Improve maintainability without changing established behavior accidentally.

### Primary Files

- `src/outputs/outputs.cpp`
- `src/outputs/outputs.hpp`
- `src/outputs/basetype_output.cpp`
- `src/outputs/io_wrapper.hpp`
- `src/globals.hpp`
- `src/restart_manifest.hpp`
- related focused tests

### Known Debt

`Outputs::Outputs` contains:

- duplicated ghost-zone parsing;
- duplicated `gid` parsing and validation;
- duplicated variable/id parsing;
- an inherited tracked-particle inconsistency;
- an overgrown inline PDF parser; and
- stale top-of-file file-type documentation.

Touched headers also under-document:

- count and offset units;
- collective-participation requirements;
- communicator lifetime;
- serial-versus-MPI mode behavior;
- restart manifest schema and invariants; and
- supported versus rejected output combinations.

### Required Design Decisions

Record:

1. helper boundaries for common output parsing;
2. whether format construction uses helper functions, a small factory, or the
   existing chain with extracted parsers;
3. whether `OutputParameters` fields receive in-class defaults;
4. how tracked-particle parsing behavior is preserved or repaired;
5. whether shared checked arithmetic/publication helpers belong in a new header;
6. whether raw owning pointers introduced by this branch should become RAII owners;
   and
7. how far cleanup extends before it becomes an unrelated refactor.

### Implementation Requirements

1. Remove duplicated parsing.
2. Preserve parser behavior unless a decision and test authorize a change.
3. Extract PDF parsing into a focused helper.
4. Initialize touched parameter fields safely.
5. Replace stale or personal comments with codebase-style explanations.
6. Add comments for non-obvious invariants, not line-by-line narration.
7. Document collective requirements in the wrapper and node-communicator headers.
8. Document manifest schema and direct-loading invariants.
9. Keep the refactor separate from functional changes where practical.

### Required Tests

1. parser regression matrix for every output family;
2. tracked-particle construction regression;
3. sharding exclusivity regression;
4. legacy PDF parser compatibility;
5. modern PDF parser matrix;
6. invalid `cbin` factors and modes;
7. invalid `sphslice` combinations;
8. output timing and final-output policy regression; and
9. full build after refactor.

### Required Subagents

1. Spawn a pre-edit code-structure auditor.
2. Spawn a behavior-preservation auditor after refactor.
3. Spawn a scope auditor that rejects unrelated cleanup.

### Reflection Point R-05

Ask:

1. Did the refactor reduce cognitive load materially?
2. Did any extracted helper conceal a format-specific rule?
3. Are defaults safe for all constructors?
4. Did the diff become harder to review than the original?
5. Should any cleanup be reverted or deferred?

### Stop Gate

Do not proceed until functional behavior is unchanged except for explicitly logged
decisions, parser regressions pass, comments are current, and a scope auditor accepts
the diff.

## RCP-06: Refactor Python Tooling Without Breaking Users

### Purpose

Keep one authoritative, strict, usable Python surface while reducing duplication and
making practical limits configurable.

### Primary Files

- `vis/python/bin_convert.py`
- `vis/python/read_pdf.py`
- `vis/python/read_sphslice.py`
- new private helper module only if justified
- `vis/python/examples/read_io_outputs.py`
- Python tests under `tst/test_suite/io/`
- visualization documentation

### Compatibility Rules

1. `vis/python/bin_convert.py` remains the sole supported binary converter module.
2. Do not restore `bin_convert_new.py`.
3. Preserve established public functions used by repository consumers.
4. Keep legacy file read compatibility.
5. Keep strict malformed-input behavior.
6. Preserve memory preflight before materialization.
7. Treat private-helper extraction as an internal refactor unless an additive public
   API is deliberately approved.

### Required Design Decisions

Record:

1. whether to add `io_reader_common.py` or keep helpers local;
2. which validation utilities are genuinely common:
   - checked products and sums;
   - bounded reads;
   - retained-memory budgets;
   - canonical shard discovery;
   - strict identifier validation;
   - inventory normalization;
3. whether practical limits are controlled through:
   - optional function keyword arguments;
   - a `ReaderLimits` object;
   - CLI flags;
   - environment variables; or
   - documented constants only;
4. backward-compatible defaults;
5. how test modules should be split without reducing coverage; and
6. whether decomposition is valuable enough to justify churn now.

### Implementation Requirements

1. Preserve public signatures unless an additive migration is logged.
2. Keep safe defaults equivalent to current caps unless evidence supports a change.
3. Avoid global mutable configuration where possible.
4. Ensure callers can read legitimate larger outputs without editing source.
5. Preserve strict rejection of malformed files before large allocations.
6. Keep aggregate reconstruction incremental.
7. Split test files by format or behavior family if it improves reviewability.
8. Update examples and visualization docs for any additive controls.

### Required Tests

1. frozen shared and per-rank fixtures;
2. `.bin` and `.cbin` read and assembly;
3. ATHDF/XDMF conversion;
4. legacy PDF;
5. modern dense PDF;
6. rank/node sparse PDF;
7. shared/rank/node `sphslice`;
8. malformed headers, payloads, inventories, aliases, and extents;
9. reduced-cap pre-materialization guards;
10. larger-limit override success;
11. invalid override rejection;
12. CLI compatibility;
13. example execution;
14. `py_compile`;
15. targeted `flake8`; and
16. import smoke for every public helper.

### Required Subagents

1. Spawn a public-API auditor before editing.
2. Spawn a memory-budget auditor before and after editing.
3. Spawn a reader-contract auditor after extraction.
4. Spawn a test-organization auditor if test files are split.

### Reflection Point R-06

Ask:

1. Did helper extraction reduce duplication without making lifecycle accounting
   harder to follow?
2. Are limits configurable without weakening safe defaults?
3. Did any public caller change unnecessarily?
4. Are tests easier to locate and maintain?
5. Is additional decomposition useful, or would it be churn?

### Stop Gate

Do not proceed until public APIs remain compatible, limits are documented, strict
preflight behavior remains proven, Python checks pass, and independent auditors
accept the refactor.

## RCP-07: Automate Deferred Pages Staging

### Purpose

Keep public documentation deferred until code merge while making later Pages
integration deterministic, reviewable, and drift-aware.

### Primary Files

- `deferred_docs/gh-pages/io-output-formats-and-sharding/MANIFEST.md`
- `deferred_docs/gh-pages/io-output-formats-and-sharding/VALIDATION.md`
- staged overlays and insertions
- a new staging helper under an appropriate script directory
- `IO_FORMAT_COMPATIBILITY.md`

### Required Behavior

Create a helper that:

1. accepts a detached refreshed `gh-pages` worktree path;
2. refuses a dirty target worktree;
3. verifies the target is detached at the intended refreshed `origin/gh-pages`
   baseline;
4. verifies expected target paths exist;
5. records and verifies the baseline Git blob ID for every whole-page replacement;
6. verifies insertion anchors;
7. fails clearly if Pages drift makes an insertion ambiguous;
8. defaults to strict mode, refusing to overwrite a replacement target whose blob
   differs from the recorded baseline;
9. supports reviewed-drift mode only by generating a three-way reconciliation
   packet and performing no replacement write until an auditor approves the merged
   page;
10. applies approved overlays idempotently;
11. inserts reviewed fragments without duplicating content;
12. preserves unrelated live reference material;
13. enforces an exact modified-file allowlist;
14. reports every copied, inserted, skipped, or rejected file;
15. emits a reviewable `git diff --stat` and `git diff`;
16. runs or prints the strict Sphinx HTML and link-check commands;
17. runs contradiction searches; and
18. never commits, pushes, or edits the live `gh-pages` checkout automatically.

### Required Design Decisions

Record:

1. script language;
2. insertion-anchor strategy;
3. idempotence markers;
4. baseline-SHA handling;
5. drift policy:
   - strict mode refuses replacement-target blob drift by default;
   - reviewed-drift mode emits a three-way reconciliation packet without writing;
   - insertion-only drift may proceed only when anchors remain unique and an
     auditor accepts the resulting diff;
6. whether build execution is automatic or an explicit flag; and
7. how the helper proves it did not touch unrelated pages.

### Two-Phase Pages Lifecycle

1. **Pre-merge preview:** apply and validate in a detached refreshed Pages worktree
   while the code branch is still under review. Do not publish.
2. **Post-code-merge restaging:** after the code branch merges, fetch the then-current
   `origin/gh-pages`, create a fresh detached worktree, rerun the helper, resolve
   drift deliberately, rerun all checks, and open a separate Pages review.

Preserve the deferred bundle until the Pages change itself merges.

### Required Validation Procedure

Use a detached worktree:

```bash
git fetch --prune origin
git worktree add --detach /tmp/athenak-gh-pages-io-docs origin/gh-pages
```

Run the helper only against the detached worktree. Then run:

```bash
cd /tmp/athenak-gh-pages-io-docs/docs
make clean html SPHINXOPTS="-W --keep-going"
make linkcheck SPHINXOPTS="-W --keep-going"
```

Also inspect:

```bash
git -C /tmp/athenak-gh-pages-io-docs status --short
git -C /tmp/athenak-gh-pages-io-docs diff --stat
git -C /tmp/athenak-gh-pages-io-docs diff
```

### Required Tests

1. clean application;
2. second idempotent application;
3. missing target rejection;
4. missing anchor rejection;
5. duplicate marker rejection;
6. unrelated content preservation;
7. drift detection;
8. contradiction search;
9. exact modified-file allowlist;
10. warnings-as-errors Sphinx HTML build;
11. warnings-as-errors link check;
12. rendered-page or browser inspection of modified public pages; and
13. confirmation that live `gh-pages` remains untouched.

### Required Subagents

1. Spawn a live-Pages structure auditor before writing the helper.
2. Spawn a helper-safety auditor after implementation.
3. Spawn a docs-to-code auditor after detached staging.
4. Spawn a navigation and Sphinx-build auditor before checkpoint closure.

### Reflection Point R-07

Ask:

1. Does the helper fail safely when Pages moves?
2. Are insertion anchors stable enough?
3. Does staged documentation describe only qualified behavior?
4. Should any page remain manual because automation would be brittle?
5. Is the later publication procedure clear to a reviewer who did not implement it?

### Stop Gate

Do not proceed until detached staging is deterministic, idempotent, drift-aware,
reviewable, warnings-as-errors clean, and independently audited. Treat this as
provisional pre-merge validation: rerun detached staging after RCP-09 and any RCP-08B
refactor, then repeat it again after the code branch merges against the then-current
Pages baseline before opening the separate Pages review.

## RCP-08: Prepare, Measure, And Decide Manifest Validation Scaling

### Purpose

Prepare a benchmark plan, collect evidence during RCP-09 multi-node qualification,
then decide whether replicated manifest and payload-header validation should remain
per rank or move to centralized validation plus metadata broadcast.

### Primary Files

- `src/restart_manifest.cpp`
- `src/restart_manifest.hpp`
- `src/main.cpp`
- scheduler qualification scripts or instructions
- benchmark records

### Background

The current implementation validates strict manifest inventory and replicated
headers independently on every rank. The prior decision log accepted this as correct
and deferred centralized distribution because a late protocol change would enlarge
the blast radius without local multi-node hardware.

Do not refactor this path merely because a centralized design sounds cleaner.
Measure first.

### Required Measurements

Collect:

1. payload count;
2. rank count;
3. node count;
4. replicated header bytes;
5. manifest bytes;
6. per-rank validation wall time;
7. aggregate filesystem read amplification;
8. startup time before and after manifest parsing;
9. direct local-span load time; and
10. filesystem and scheduler environment.

Run at multiple rank and node counts representative of production during RCP-09.

Before consuming scheduler time, preregister:

1. node counts;
2. ranks per node;
3. physical hostnames or the command that will capture them;
4. payload-size targets;
5. replicated-header-size targets;
6. repeat count;
7. filesystem;
8. measured percentiles;
9. acceptable startup-overhead threshold; and
10. the threshold that triggers centralized validation work or a separate scaling
    branch.

### Two-Stage Execution Rule

RCP-08 begins before RCP-09 and closes after RCP-09 baseline measurements:

1. **RCP-08A: prepare.** Define measurements, add narrowly scoped instrumentation if
   needed, and obtain benchmark-design audit approval.
2. **RCP-09 baseline execution.** Run scheduler-backed qualification and collect the
   RCP-08 measurements.
3. **RCP-08B: decide.** Review evidence and record a keep, refactor, or
   separate-branch decision.
4. **If refactored:** repeat the affected restart portion of RCP-09 and rerun focused
   restart regressions before closing either checkpoint.

### Required Design Decision

Choose:

1. keep replicated validation and document measured acceptability;
2. parse and validate on rank 0, then broadcast structured metadata;
3. validate headers on node leaders, then distribute node-local metadata; or
4. stage a separate scaling branch if the optimization is material but too risky for
   this feature branch.

Record correctness, complexity, performance, and qualification tradeoffs.

### Required Tests If Refactored

1. malformed manifest rejection;
2. mismatching replicated header rejection;
3. missing payload rejection;
4. payload marker rejection;
5. traversal and alias rejection;
6. changed rank count;
7. multiple nodes;
8. empty/non-owning node;
9. forced chunking;
10. collective error handling; and
11. no `.assembled` staging.

### Required Subagents

1. Spawn a benchmark-design auditor before running measurements.
2. Spawn a restart-protocol auditor before any refactor.
3. Spawn a performance-evidence auditor after measurements.
4. If refactored, spawn a different correctness auditor and re-auditor.

### Reflection Point R-08

Ask:

1. Is validation overhead material at production scale?
2. Would a protocol change create more risk than value?
3. Should optimization remain deferred with measured evidence?
4. Is a separate scaling branch more reviewable?

### Stop Gate

Before starting RCP-09, do not advance past RCP-08A without an accepted benchmark
plan and any required instrumentation. After RCP-09 baseline execution, do not close
RCP-08 without benchmark evidence and a recorded keep, refactor, or separate-branch
decision. If the path is refactored, rerun the affected RCP-09 restart qualification
before closing RCP-08 or RCP-09.

## RCP-09: External CUDA And Multi-Node Qualification

### Purpose

Close the two environment-dependent gates that local workstation testing cannot
prove.

### CUDA Qualification

Use a CUDA-capable AthenaK build. Record:

1. host;
2. compiler;
3. CUDA toolkit;
4. Kokkos configuration;
5. CMake command;
6. build command;
7. selected GPU;
8. test command;
9. exact result;
10. generated output inventory; and
11. any warnings.

Run at minimum:

- `tst/test_suite/io/test_output_formats_gpu.py`;
- representative N-D PDF axes using derived coordinates and velocity projection;
- scalar weighting;
- one Hydro case;
- one MHD case if supported by the harness; and
- any new analytic diagnostic GPU coverage from RCP-04.

### Multi-Node Qualification

Run under the target scheduler on at least two physical nodes with more than one rank
per node. Include:

1. full-volume node-sharded `.bin`;
2. sliced node-sharded `.bin`;
3. full-volume node-sharded `.cbin`;
4. node-sharded modern PDF;
5. node-sharded `sphslice`;
6. node-sharded restart publication;
7. native direct restart resume;
8. changed rank count;
9. changed rank distribution across nodes where feasible;
10. a genuinely empty or non-owning node;
11. output timing labels;
12. expected shard inventory;
13. reader assembly and equality checks;
14. no `.assembled` staging;
15. parallel-filesystem behavior; and
16. manifest-validation scaling measurements for RCP-08B.

### Qualification Evidence Template

```markdown
### External Qualification: Short Title

| Field | Record |
| --- | --- |
| Environment | Host, scheduler, filesystem, compiler, CUDA, MPI, Kokkos |
| Branch snapshot | Exact SHA |
| Guide snapshot | Exact SHA-256 checksum and line count |
| Scheduler evidence | Job ID plus physical rank-to-host mapping |
| CUDA evidence | Device inventory and selected device when applicable |
| Build | Exact commands |
| Runtime topology | Nodes, ranks, ranks per node, empty/non-owning node arrangement |
| Commands | Exact execution and reader commands with exit codes |
| Retained logs | Paths to scheduler, test, and timing logs |
| Output inventory | Expected and observed files |
| Results | Test results and numerical comparisons |
| Timing | Relevant measurements |
| Working tree | State before and after validation |
| Failures | Any failure and disposition |
| Auditor | Independent qualification reviewer |
| Status | Passed, failed, or incomplete |
```

### Required Subagents

1. Spawn a qualification-plan auditor before consuming scheduler time.
2. Spawn an evidence auditor after CUDA execution.
3. Spawn a separate topology-evidence auditor after multi-node execution.
4. If any failure occurs, spawn a bounded root-cause auditor before editing.

### Reflection Point R-09

Ask:

1. Did external execution invalidate any local assumption?
2. Is empty-node behavior correct and understandable?
3. Do timings support the feature's intended scaling benefit?
4. Is manifest validation acceptable, or must the flow return to RCP-08B for a
   refactor decision?
5. Does any format require a deployment-specific note?

### Stop Gate

Do not claim merge readiness until both external qualifications pass and independent
auditors agree that the evidence proves the intended behavior.

## RCP-10: Process Packaging And Final Whole-Branch Review

### Purpose

Prepare a reviewable final branch without discarding useful history or carrying
unnecessary process bulk accidentally.

### Process Artifact Decision

Review:

- `IO_FEATURE_BRANCH_GUIDE.md`;
- `IO_FEATURE_BRANCH_INTEGRATION_PLAN.md`;
- `IO_FEATURE_BRANCH_DECISION_LOG.md`;
- `IO_FEATURE_AUDIT_LEDGER.md`;
- `IO_FORMAT_COMPATIBILITY.md`;
- this robustification guide; and
- deferred Pages records.

For each artifact decide:

1. retain in the final code branch;
2. retain but move under a development-documentation archive;
3. summarize into a smaller durable record;
4. preserve in the pull request rather than merge; or
5. remove only after its durable content is migrated.

Do not remove or relocate records silently. Add a decision-log entry.

### Final Full Validation Matrix

Re-run and record:

1. serial build;
2. MPI build;
3. full serial IO pytest matrix;
4. full MPI IO pytest matrix;
5. Python reader pytest matrix;
6. actual CUDA IO regression;
7. scheduler-backed multi-node matrix;
8. fixture checksum verification;
9. frozen shared restart resume;
10. frozen per-rank restart resume;
11. generated node-manifest restart resume;
12. promoted examples;
13. Python `py_compile`;
14. targeted `flake8`;
15. repository style checks;
16. `git diff --check origin/main...HEAD`;
17. contradiction searches;
18. detached Pages helper application;
19. detached Pages warnings-as-errors Sphinx build; and
20. clean working-tree inspection.

### Contradiction Searches

At minimum:

```bash
rg -n "\\.assembled|StageNodeRestart|CopyFileRange|bin_convert_new" src vis inputs
rg -n "single_file_per_node|sphslice|final_output_policy|output_timing|AKPDFV2" \
  src vis inputs deferred_docs IO_FORMAT_COMPATIBILITY.md
rg -n "TODO|FIXME|DBF|not sure|probably won't work" \
  src/outputs src/main.cpp src/pgen/pgen.cpp src/restart_manifest.* src/globals.*
```

Classify every match. Do not remove unrelated baseline comments casually, but do not
leave uncertainty in touched feature code.

### Required Final Subagents

Spawn separate read-only auditors:

1. restart and MPI red-team auditor;
2. file-format and compatibility auditor;
3. numerical diagnostic auditor;
4. Python public-API and malformed-input auditor;
5. tests, examples, fixture, and qualification-evidence auditor;
6. deferred Pages and documentation auditor;
7. scope and process-packaging auditor; and
8. final whole-branch red-team auditor with permission to revisit earlier decisions.

Do not ask one agent to collapse these lanes into a superficial summary.

### Reflection Point R-10

Ask:

1. What did the fresh red-team audit find that earlier audits missed?
2. Are any "locally complete" claims still relying on external assumptions?
3. Does the branch remain one coherent IO feature rather than a cleanup branch?
4. Are all user-facing guarantees documented?
5. Are all documented guarantees tested?
6. Can a new user read every new output format with shipped helper tools?
7. Can a reviewer understand why rejected alternatives were rejected?
8. Is the later Pages application procedure deterministic and isolated?

### Stop Gate

Do not call the branch merge-ready until:

- all blocking red-team findings are resolved;
- every correction is re-audited;
- the full matrix passes;
- external gates pass;
- docs match code;
- process records have an intentional disposition; and
- the final audit-ledger closure entry names all residual risks explicitly.

## Suggested Verification Command Inventory

The implementing agent must inspect the current repository harness and adjust paths
if needed. Do not copy commands blindly. A starting inventory is:

### Serial Build

```bash
cmake -S . -B /tmp/athenak-io-robust-build \
  -DCMAKE_BUILD_TYPE=Debug \
  -DAthena_ENABLE_MPI=OFF
cmake --build /tmp/athenak-io-robust-build -j 4
```

### MPI Build

```bash
cmake -S . -B /tmp/athenak-io-robust-build-mpi \
  -DCMAKE_BUILD_TYPE=Debug \
  -DAthena_ENABLE_MPI=ON
cmake --build /tmp/athenak-io-robust-build-mpi -j 4
```

### Fixture Integrity

```bash
cd tst/fixtures/io/origin_main_886dd2a1
shasum -a 256 -c SHA256SUMS
```

### Python Syntax And Style

```bash
/Users/dbf75/.uv/envs/interactive/.venv/bin/python -m py_compile \
  vis/python/bin_convert.py \
  vis/python/read_pdf.py \
  vis/python/read_sphslice.py
```

Run targeted `flake8` over changed Python files and pytest modules using the
repository's established invocation.

### IO Pytest Inventory

Inspect and run all applicable modules under:

```text
tst/test_suite/io/
```

The suite currently includes:

- `test_chunked_io_mpicpu.py`
- `test_io_examples_cpu.py`
- `test_io_finalization_timing_cpu.py`
- `test_io_finalization_timing_mpicpu.py`
- `test_io_wrapper_harness_mpicpu.py`
- `test_node_sharding_mpicpu.py`
- `test_output_formats_cpu.py`
- `test_output_formats_gpu.py`
- `test_output_formats_mpicpu.py`
- `test_python_io_readers_cpu.py`
- `test_writer_hardening_cpu.py`
- `test_writer_hardening_mpicpu.py`

### Canonical Matrix Registration

During RCP-00, record the exact canonical pytest commands in
`IO_FEATURE_AUDIT_LEDGER.md`. Do not leave the final matrix as an informal selection
of "applicable" files.

The prior local baseline recorded these minimum floors:

| Matrix | Prior floor |
| --- | --- |
| Python reader module | `148 passed` |
| Full serial IO plus CPU smoke of GPU-selectable regression | `204 passed` |
| Full MPI IO | `59 passed` |
| Repository style | `2 passed` |
| Immutable fixtures | `27` artifacts verified |

Rules:

1. Preserve or increase test coverage as behavior expands.
2. If a count decreases, stop and explain every removed or deselected case.
3. Record exact commands, environment variables, working directories, exit codes,
   and retained log paths.
4. Add explicit timeouts for MPI tests that exercise mismatch, failure, or deadlock
   boundaries.
5. Keep CUDA execution separate from CPU smoke evidence.
6. Keep scheduler-backed multi-node qualification separate from one-node MPI
   evidence.

### Diff Hygiene

```bash
git diff --check origin/main...HEAD
git status --short --branch
```

## Test-Failure Triage Protocol

When a test fails:

1. stop expanding scope;
2. record the exact command and failure;
3. determine whether the failure is:
   - intended new rejection;
   - implementation defect;
   - incorrect test oracle;
   - stale docs/example expectation;
   - environment issue;
   - pre-existing baseline issue; or
   - external-qualification limitation;
4. spawn a bounded root-cause auditor for ambiguous failures;
5. add or update a decision entry if behavior is ambiguous;
6. fix the smallest coherent cause;
7. rerun the focused test;
8. rerun neighboring regressions;
9. inspect the diff;
10. request re-audit; and
11. only then resume the checkpoint.

Do not weaken validation or delete a failing test merely to restore green status.

## Commit Discipline

Prefer reviewable commits organized by coherent behavior:

1. restart layout arithmetic and serial positioned IO;
2. MPI checking and filesystem publication helpers;
3. `cbin` contract enforcement;
4. diagnostic semantics, ghost-zone boundaries, and `sphslice` naming;
5. output registration and header-comment refactor;
6. Python helper extraction and configurable limits;
7. tests and examples;
8. deferred Pages staging helper and documentation refresh;
9. benchmark or scaling-path changes if justified;
10. process-record packaging and final audit updates.

Rules:

1. Do not hide functional changes inside refactor commits.
2. Do not mix deferred Pages publication with code changes.
3. Keep audit-ledger and decision-log updates close to the commits they explain.
4. Run focused tests before each commit.
5. Run broader tests before pushing.
6. Inspect every staged diff before committing.

## Final Deliverables

The completed branch should contain:

1. end-to-end checked restart arithmetic;
2. checked serial positioned IO;
3. consistent MPI error handling;
4. explicit directory and publication helpers;
5. sequence-token-safe output naming beyond dump `99999`;
6. bounded node-manifest parsing for individual records and total file size;
7. deterministic output namespace mapping with duplicate-target rejection;
8. an enforced `cbin` support matrix;
9. checked `cbin` Kokkos range arithmetic;
10. fully qualified diagnostic semantics;
11. an explicit derived ghost-zone contract;
12. collision-resistant `sphslice` radius naming;
13. writer-side PDF and `sphslice` capacity preflight before allocation;
14. measured manifest-validation scaling with a recorded keep-or-refactor decision;
15. cleaner output registration and current source comments;
16. maintainable Python readers with preserved public compatibility;
17. configurable practical reader limits with safe defaults if approved;
18. expanded CPU, MPI, GPU, malformed-input, fixture, and scheduler-backed tests;
19. runnable examples for supported workflows;
20. a deterministic deferred Pages staging helper;
21. refreshed Pages-compatible content that remains unpublished until code merge;
22. an updated compatibility contract;
23. an updated decision log with alternatives and reversal paths;
24. an updated audit ledger with independent audits and reflection points; and
25. a deliberate disposition for large process records.

## Final Instruction To Implementing Agents

Proceed slowly. A quick green test run is not completion. A prior audit is not
completion. A confident implementation summary is not completion.

For every checkpoint:

1. inspect before editing;
2. write decisions before coding through ambiguity;
3. implement narrowly;
4. test the affected behavior;
5. spawn independent subagents to challenge the result;
6. correct every blocking finding;
7. spawn re-auditors after correction;
8. reflect on whether the plan still makes sense;
9. update the plan when evidence changes the direction; and
10. only then close the checkpoint.

The branch should become robust because each claim has evidence and each design
choice is explicit, not because the implementation was completed quickly.
