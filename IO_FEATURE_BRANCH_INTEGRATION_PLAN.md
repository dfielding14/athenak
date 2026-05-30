# IO Output Formats And Sharding Selective Integration Plan

## Document Status

| Item | Value |
| --- | --- |
| Target branch | `feature/io-output-formats-and-sharding` |
| Clean base | `origin/main` at `886dd2a1437e45a3a30b3eeebf2adfa838328f73` |
| Current local state | Reconstructed implementation preserved in `ca00581b`, `3049aa92`, and `fcc534fe`; CP-01 planning and deferred documentation corrections remain uncommitted |
| Comparison reference | Local remote-tracking ref `origin/feature/single-file-per-node-outputs` at `47462c5da45d3c763b37fb21505fac5fd3498805` |
| Comparison rule | Use the remote branch as read-only evidence. Do not merge or cherry-pick it wholesale. |
| Decision record | `IO_FEATURE_BRANCH_DECISION_LOG.md` |
| Existing implementation audit record | `IO_FEATURE_AUDIT_LEDGER.md` |
| Compatibility contract | `IO_FORMAT_COMPATIBILITY.md` |

This document governs the next implementation pass. The local branch already
contains a substantial clean reconstruction of Gotham IO behavior. The goal is
not to replace that work with `origin/feature/single-file-per-node-outputs`.
The goal is to preserve the stronger local design, selectively port useful
remote hardening, close the remaining scalability gaps, and leave a
well-tested branch that can be reviewed and merged intentionally.

## Executive Goal

Produce one reviewable IO feature branch that contains:

1. `single_file_per_node` sharding for supported binary, coarsened-binary, PDF,
   spherical-slice, and restart products.
2. Versioned N-dimensional PDF output with legacy one- and two-dimensional PDF
   compatibility.
3. Shared, rank-sharded, and node-sharded spherical-slice output with strict
   reconstruction validation.
4. A transactional per-node restart format with direct distributed loading.
5. Chunked MPI byte IO that does not rely on one `int`-sized MPI count.
6. Opt-in timing diagnostics and explicit final-output policy.
7. One canonical `vis/python/bin_convert.py` module with the useful newer
   convenience APIs folded in without breaking established callers.
8. Runnable examples, negative tests, frozen compatibility fixtures, deferred
   GitHub Pages content, and an auditable evidence trail.

## Non-Negotiable Rules

### Integration Rules

1. Do not merge `origin/feature/single-file-per-node-outputs`.
2. Do not cherry-pick its four commits wholesale.
3. Port behavior by intent after inspecting the exact remote implementation.
4. Keep the local implementation as the authority for public file formats,
   compatibility guarantees, runtime policy, and Python reader strictness.
5. Keep shared and existing per-rank behavior working while adding node support.
6. Do not modify or publish the live `gh-pages` branch from this feature branch.
7. Do not claim production readiness until GPU and real multi-node
   qualifications are complete.

### Worktree Rules

1. Start every checkpoint with `git status --short --branch`.
2. Do not discard or overwrite existing local changes.
3. Record the exact branch, `HEAD`, comparison ref, and merge base before edits.
4. Use isolated build directories and isolated test run directories.
5. Keep temporary comparison snapshots under `/tmp` and remove them after use.
6. Run `git diff --check` after each implementation checkpoint.
7. Append evidence to `IO_FEATURE_AUDIT_LEDGER.md`.
8. Append decisions to `IO_FEATURE_BRANCH_DECISION_LOG.md`.

## Current Baseline Assessment

### Local Implementation To Preserve

The local tree is the preferred baseline for:

| Area | Local strength |
| --- | --- |
| PDF format | Versioned `AKPDFV2` payloads plus byte-compatible legacy text PDFs |
| PDF reader | Reads legacy text, current `AKPDFV2`, and the remote branch's transitional unversioned dense/sparse payloads |
| Spherical-slice reader | Strict metadata, payload-size, duplicate-ownership, and missing-angle validation |
| Spherical-slice semantics | Rejects derived fields until ghost-zone-safe interpolation exists |
| Generic diagnostics | Rejects ambiguous generic flow/flux diagnostics for two-fluid ion-neutral runs |
| Node setup | Lazily initializes node communicators only when required |
| Restart publication | Atomically publishes generation-qualified payloads before atomically publishing a text manifest |
| Runtime policy | Provides opt-in `<time>/output_timing` and `<time>/final_output_policy = all|restart_only|none` |
| Converter consolidation | Exposes one canonical `vis/python/bin_convert.py` and removes `bin_convert_new.py` redundancy |
| Tests and docs | Uses the repository pytest harness, frozen fixtures, promoted examples, and a deferred Pages overlay |

### Remaining Local Gaps

The next pass must address these gaps before merge readiness:

| ID | Gap | Why it matters |
| --- | --- | --- |
| GAP-001 | Per-node restart input reconstructs `manifest.rst.assembled` on rank 0 before using the legacy reader | Correct but not scalable for large restart files |
| GAP-002 | `IOWrapper` lacks generalized chunked MPI byte operations and forced-small-chunk tests | Large IO operations still depend on MPI `int` count limits |
| GAP-003 | Local canonical `bin_convert.py` lacks useful newer helpers retained by the remote consolidation | One canonical module should include the useful modern API surface |
| GAP-004 | Local `tst/scripts/utils/athena.py` still constructs `['mpiexec -n', ...]` | The older test helper invokes MPI incorrectly |
| GAP-005 | Empty-node sliced-output behavior is not qualified on multiple physical nodes | One-node local MPI tests cannot prove true node-shard behavior |
| GAP-006 | The compatibility record overstates the restart loader as native | Documentation must match the implementation at every checkpoint |
| GAP-007 | Resolved during CP-01: earlier planning artifacts used a superseded prefixed branch name | Branch identity should remain consistent |

## Subagent Operating Protocol

### Why Subagents Are Required

Restart IO, MPI collectives, file-format compatibility, Python API
compatibility, and deferred documentation are separate failure surfaces. Use
subagents to audit them independently. The main implementing agent remains
responsible for integration, conflict resolution, verification, and final
judgment.

### When To Spawn Subagents

Before using subagents, discover the available multi-agent tools with
`tool_search`. Spawn bounded workers or auditors only after assigning an
explicit scope.

Use subagents in these cases:

| Situation | Required action |
| --- | --- |
| A checkpoint changes MPI collectives or restart loading | Spawn one implementation reviewer and one read-only correctness auditor |
| A checkpoint changes a public Python API | Spawn a read-only API compatibility auditor |
| A checkpoint changes file layouts or reader behavior | Spawn a format-contract auditor |
| A checkpoint adds or modifies tests | Spawn a test-coverage auditor after implementation |
| Deferred Pages content is updated | Spawn a docs-to-code consistency auditor |
| Final merge readiness is being considered | Spawn a final independent audit that did not implement the changes |

### Subagent Scope Rules

1. Give each subagent a bounded scope and a concrete deliverable.
2. Prefer read-only audits for cross-cutting review.
3. If an implementation worker edits files, assign disjoint file ownership and
   use a separate worktree.
4. Require every subagent to report:
   - files inspected or edited;
   - commands run;
   - findings ordered by severity;
   - tests run and their results;
   - unresolved risks;
   - recommended decision-log entries.
5. Do not accept a subagent report as proof by itself. Re-run critical checks in
   the main worktree.
6. Record the audit and the main-agent disposition in
   `IO_FEATURE_AUDIT_LEDGER.md`.

### Prompt Template For Read-Only Auditors

```text
Audit only. Do not edit files, create commits, or change refs.

Repository:
  /Users/dbf75/.codex/worktrees/e948/athenak-DF

Target branch:
  feature/io-output-formats-and-sharding

Reference branch:
  origin/feature/single-file-per-node-outputs

Scope:
  <bounded files and behavior>

Questions:
  1. What correctness defects or regressions remain?
  2. What useful behavior exists in the reference branch but not locally?
  3. What behavior must not be ported?
  4. What tests are missing?

Return findings first, ordered by severity, with file and line references.
Include commands run and a short recommended decision-log update.
```

### Prompt Template For Implementation Workers

```text
Implement only the assigned checkpoint in a separate worktree.
Do not merge, cherry-pick, or copy the reference branch wholesale.
Do not touch files outside the assigned ownership list.

Preserve:
  - existing shared and per-rank behavior;
  - local file-format contracts;
  - local strict reader validation;
  - transactional restart publication;
  - lazy node communicator setup.

Reference behavior may be inspected from:
  origin/feature/single-file-per-node-outputs

Before returning:
  - run focused tests;
  - run git diff --check;
  - report changed files, commands, results, and unresolved risks.
```

## Checkpoint Overview

| Checkpoint | Purpose | Must pass before continuing |
| --- | --- | --- |
| CP-00 | Freeze facts and preserve the reconstructed baseline | Branch identity, worktree inventory, and baseline verification recorded |
| CP-01 | Correct planning records and finalize contracts | Decision log, compatibility contract, and branch naming agree |
| CP-02 | Add chunked MPI byte IO | Serial/MPI wrapper behavior and forced-small-chunk tests pass |
| CP-03 | Replace staged restart assembly with native distributed reads | No `.assembled` file is needed; restart regressions pass |
| CP-04 | Harden node-sharded writers and empty-node behavior | Binary, coarsened-binary, and spherical-slice tests pass |
| CP-05 | Complete canonical Python converter API | One converter module exposes tested legacy and useful modern APIs |
| CP-06 | Expand regression and qualification coverage | CPU, MPI, negative, fixture, GPU, and multi-node requirements are explicit |
| CP-07 | Reconcile examples and deferred Pages content | Docs match code and build with warnings as errors in a temporary Pages worktree |
| CP-08 | Perform final independent audits and prepare commits | All blocking findings are resolved or explicitly deferred outside scope |

## CP-00: Freeze Facts And Preserve The Reconstructed Baseline

### Purpose

Preserve the existing local reconstruction before layering selective ports on
top of it. The current implementation is largely uncommitted. Do not start MPI
or restart rewrites until its contents and verification evidence are captured.

### Required Actions

1. Capture:

   ```bash
   git status --short --branch
   git rev-parse HEAD origin/main origin/feature/single-file-per-node-outputs
   git merge-base HEAD origin/feature/single-file-per-node-outputs
   git diff --stat origin/main
   git ls-files --others --exclude-standard
   git diff --check origin/main
   ```

2. Confirm that the active branch is
   `feature/io-output-formats-and-sharding`.
3. Review every modified and untracked file before staging.
4. Re-run the existing local serial and MPI IO suites against isolated builds.
5. Verify frozen fixture checksums.
6. Preserve the local reconstruction in logical commits before adding remote
   hardening. Keep tests with the behavior they verify.

### Recommended Baseline Commit Shape

1. `outputs: add runtime timing and explicit final-output policy`
2. `outputs: add N-dimensional PDFs and generic diagnostics`
3. `outputs: add spherical-slice format and strict readers`
4. `outputs: add node-sharded products and transactional restarts`
5. `tools: consolidate binary conversion and add IO readers`
6. `tests: add IO fixtures, examples, and deferred documentation package`

Adjust the split if the actual diff shows tighter dependencies. Do not split a
commit in a way that leaves the tree uncompilable.

### Subagent Assignment

Spawn one read-only baseline auditor. Ask it to verify that the staged baseline
matches the documented local contracts and that unrelated files are not
included.

### Stop Gate

- [ ] Branch identity recorded.
- [ ] Existing local files inventoried.
- [ ] Baseline tests re-run.
- [ ] Fixture checksums pass.
- [ ] Baseline commits created without unrelated changes.
- [ ] Baseline audit recorded in `IO_FEATURE_AUDIT_LEDGER.md`.

## CP-01: Correct Planning Records And Finalize Contracts

### Purpose

Make the planning artifacts truthful before code changes continue.

### Required Actions

1. Replace stale branch references with
   `feature/io-output-formats-and-sharding`.
2. Correct the current restart-loader description:
   - publication is transactional;
   - manifest validation is strict;
   - loading currently stages a temporary `.assembled` shared file;
   - CP-03 will replace that staging path with native distributed reads.
3. Update `IO_FORMAT_COMPATIBILITY.md` only when implementation evidence
   changes.
4. Append every accepted, rejected, pending, or superseded decision to
   `IO_FEATURE_BRANCH_DECISION_LOG.md`.

### Subagent Assignment

Spawn one read-only consistency auditor for:

- `IO_FEATURE_BRANCH_GUIDE.md`
- `IO_FEATURE_BRANCH_INTEGRATION_PLAN.md`
- `IO_FEATURE_BRANCH_DECISION_LOG.md`
- `IO_FEATURE_AUDIT_LEDGER.md`
- `IO_FORMAT_COMPATIBILITY.md`
- `deferred_docs/gh-pages/io-output-formats-and-sharding/`

### Stop Gate

- [x] Branch naming is consistent.
- [x] Restart-loader claims match the current code.
- [x] Deferred docs do not claim behavior that has not landed.
- [x] Decision log entries D-001 through D-024 are reviewed.

## CP-02: Add Chunked MPI Byte IO

### Purpose

Generalize MPI IO so large reads, writes, and broadcasts are not limited by one
`int`-sized MPI count. This is useful independently and is a prerequisite for
native distributed restart loading.

### Reference Behavior To Inspect

Inspect the remote branch versions of:

- `src/outputs/io_wrapper.cpp`
- `src/outputs/io_wrapper.hpp`
- `src/parameter_input.cpp`
- `src/parameter_input.hpp`
- `tst/scripts/restart/per_node_restart.py`

Do not copy those files wholesale.

### Required Design

1. Add checked multiplication for `size * count`.
2. Add chunked byte helpers for:
   - sequential read;
   - offset read;
   - collective offset read;
   - sequential write;
   - offset write;
   - collective offset write;
   - broadcast.
3. Keep zero-byte collective participants inside collective calls.
4. Preserve serial behavior.
5. Add a test-only environment override such as
   `ATHENAK_TEST_MAX_MPI_BYTES` to force tiny chunks.
6. Keep the production default bounded by the MPI count limit.
7. Avoid unrelated `ParameterInput` behavior changes. If a larger parameter
   header limit is needed, justify and test it separately.

### Files Expected To Change

- `src/outputs/io_wrapper.cpp`
- `src/outputs/io_wrapper.hpp`
- `src/parameter_input.cpp` only if broadcast migration is needed
- focused tests under `tst/test_suite/io/`

### Subagent Assignments

1. Spawn one bounded implementation worker for `io_wrapper.*`.
2. Spawn one read-only MPI auditor after integration.
3. Require the auditor to inspect zero-byte participants, overflow checks,
   serial fallbacks, communicator selection, and collective call symmetry.

### Verification

- [x] Serial build passes.
- [x] MPI build passes.
- [x] Existing shared and per-rank IO tests pass.
- [x] Forced-small-chunk shared restart write/read passes.
- [ ] Forced-small-chunk node restart write/read passes after CP-03.
- [x] `git diff --check` passes.

### Stop Gate

Do not proceed to CP-03 until chunked wrapper behavior is independently audited.

## CP-03: Replace Staged Restart Assembly With Native Distributed Reads

### Purpose

Remove the rank-0 `.assembled` reconstruction bottleneck while preserving the
local transactional manifest and strict validation contract.

### Design Principle

Keep both good intents:

1. Keep the local text manifest, generation-qualified payloads, atomic payload
   publication, atomic manifest publication, and strict inventory validation.
2. Adapt the remote branch's native distributed loading strategy, coalesced
   read spans, node-collective reads, and chunked transfer handling.

Do not retain both implementations permanently. The final path must not require
`manifest.rst.assembled`.

### Required Architecture

Prefer extracting restart-manifest logic from `src/main.cpp` into a dedicated
module rather than expanding `main.cpp` further. The final design should expose
a structured manifest object used by both path normalization and restart
loading.

The native flow should:

1. Accept a public manifest path.
2. Optionally accept a node payload path only when it can be normalized
   unambiguously to its public manifest.
3. Reject malformed shard directories, path traversal, absolute payload paths,
   mixed generations, missing nodes, duplicate nodes, incomplete manifests,
   inconsistent byte counts, and incomplete segment coverage.
4. Parse and validate the public text manifest once in a controlled location.
5. Broadcast validated structured metadata where required.
6. Read parameter and mesh metadata without assembling a full shared payload.
7. Route MeshBlock field reads directly to node payload shards.
8. Coalesce contiguous payload requests into spans.
9. Open payload shards collectively over the relevant node communicator.
10. Include zero-byte collective participants.
11. Use CP-02 chunked IO for large reads.
12. Preserve shared and existing per-rank restart behavior.
13. Preserve ordinary restart counter advancement so resumed runs cannot
    overwrite terminal checkpoints.

### Explicitly Rejected Remote Behavior

Do not port:

- the remote binary restart manifest schema;
- eager node communicator initialization at process startup;
- stale-payload auto-detection that can silently reinterpret an ordinary shared
  restart;
- unrelated mesh-diagnostics refactors;
- unrelated particle restart changes.

### Files Expected To Change

Likely files:

- `src/main.cpp`
- a new dedicated restart-manifest module if justified
- `src/mesh/build_tree.cpp`
- `src/mesh/mesh.cpp`
- `src/mesh/mesh.hpp`
- `src/outputs/io_wrapper.cpp`
- `src/outputs/io_wrapper.hpp`
- `src/parameter_input.cpp`
- `src/parameter_input.hpp`
- `src/pgen/pgen.cpp`
- `src/pgen/pgen.hpp`
- `src/outputs/restart.cpp` only if manifest metadata must be extended
- focused tests under `tst/test_suite/io/`

### Subagent Assignments

1. Spawn one read-only restart-design auditor before implementation.
2. Spawn one bounded implementation worker in a separate worktree after the
   design is reviewed.
3. Spawn a different read-only auditor after integration.
4. Ask the final auditor to prove that no `.assembled` path remains in the
   production node-restart flow.

### Required Tests

- [ ] Shared restart round trip.
- [ ] Existing per-rank restart round trip.
- [ ] Per-node restart by manifest path.
- [ ] Per-node restart by normalized node-payload path if that API is accepted.
- [ ] No `.assembled` file is created during node restart.
- [ ] Terminal node checkpoint resume advances numbering without overwrite.
- [ ] Corrupted manifest rejection: traversal, absolute path, incomplete marker,
      wrong byte count, missing node, duplicate node, mixed generation,
      overlapping segment, missing segment coverage.
- [ ] Missing payload rejection.
- [ ] Truncated payload rejection.
- [ ] Forced-small-chunk restart read and write.
- [ ] Zero-payload rank participation.
- [ ] Real multi-node qualification with an empty or non-owning node.

### Stop Gate

Do not describe the restart reader as native until:

- [ ] production restart loading no longer creates `.assembled`;
- [ ] direct distributed reads pass focused MPI tests;
- [ ] the new implementation has an independent correctness audit.

## CP-04: Harden Node-Sharded Writers And Empty-Node Behavior

### Purpose

Carry forward useful remote hardening for binary, coarsened-binary, and
spherical-slice output without regressing local format contracts.

### Binary And Coarsened-Binary Requirements

1. Define the contract for a node with no selected sliced output:
   - either skip the shard deliberately and teach readers the contract; or
   - write an explicit valid empty shard.
2. Ensure every rank advances output counters consistently.
3. Ensure collective calls include zero-byte ranks.
4. Preserve full-volume node-sharded `.cbin`.
5. Keep sliced node-sharded `.cbin` outside the promoted workflow until the
   existing zero-width extent defect is repaired separately.
6. Add optional phase-level stats only behind an explicit opt-in mechanism.

### Spherical-Slice Requirements

1. Preserve the local reader format and strict validation.
2. Add checked writes and checked close handling to the writer.
3. Validate node-merged angular ownership before publication:
   - index range;
   - duplicate ownership;
   - complete ownership where required;
   - deterministic ordering.
4. Keep derived spherical-slice fields rejected until ghost-zone-safe sampling
   exists.
5. Keep legacy `file_type=sph` VTK output unchanged.

### Subagent Assignments

1. Spawn one implementation worker for binary and coarsened-binary hardening.
2. Spawn one implementation worker for spherical-slice writer hardening.
3. Use separate worktrees and disjoint file ownership.
4. Spawn one read-only MPI auditor across both diffs before integration closes.

### Verification

- [ ] Shared, rank, and node `.bin` reconstruction equality.
- [ ] Full-volume shared, rank, and node `.cbin` reconstruction equality.
- [ ] Sliced `.bin` empty-owner tests.
- [ ] Shared, rank, and node spherical-slice reconstruction equality.
- [ ] Spherical-slice malformed payload and ownership rejection.
- [ ] Explicit exclusion remains documented for sliced node-sharded `.cbin`.
- [ ] Real multi-node empty-node qualification is recorded.

## CP-05: Complete The Canonical Python Converter API

### Purpose

Keep one supported converter module while restoring useful modern convenience
APIs from the remote consolidation.

### Required API Decisions

Keep:

- `read_binary`
- `read_coarsened_binary`
- `read_all_ranks_binary`
- `read_all_ranks_coarsened_binary`
- `read_binary_as_athdf`
- `read_all_ranks_binary_as_athdf`
- `read_all_ranks_coarsened_binary_as_athdf`
- `read_single_rank_binary_as_athdf`
- `read_coarsened_binary_as_athdf`
- `write_athdf`
- `write_xdmf_for`
- `convert_file`
- the CLI with `--assemble-shards`

Add:

- `read_rank_binary_as_athdf`
- an additive indexed-meshblock option for
  `read_single_rank_binary_as_athdf`, if it can preserve existing callers

Do not add by default:

- a new `bin_convert_new.py`;
- a compatibility alias named `bin_convert_new.py`;
- the remote ad hoc `athinput()` parser unless an actual in-repository consumer
  is identified and its behavior is tested.

### Compatibility Rules

1. Existing positional and keyword call forms must continue to work.
2. New indexed behavior must be additive and unambiguous.
3. Shared, rank-sharded, and node-sharded binary assembly must retain strict
   metadata validation.
4. Empty shards must remain valid where the file-format contract permits them.
5. CLI output naming must remove only the final `.bin` or `.cbin` suffix.

### Subagent Assignments

1. Spawn one read-only Python API auditor before implementation.
2. Spawn one bounded implementation worker for
   `vis/python/bin_convert.py` and its focused tests.
3. Spawn a different read-only auditor after integration to compare public
   symbols and call signatures against both the base and remote variants.

### Verification

- [ ] Public API inventory test passes.
- [ ] Legacy converter call forms pass.
- [ ] Indexed meshblock read passes if added.
- [ ] `read_rank_binary_as_athdf()` passes on a rank shard and a node shard.
- [ ] Shared `.bin` and `.cbin` CLI conversion passes.
- [ ] `--assemble-shards` rank and node conversion passes.
- [ ] `rg -n "bin_convert_new" .` finds no supported code or docs dependency.

## CP-06: Expand Regression And Qualification Coverage

### Purpose

Translate useful remote scenarios into the repository's pytest harness and add
the qualification runs that local single-node development cannot prove.

### Test Placement Rules

1. Put promoted tests under `tst/test_suite/io/`.
2. Put focused input decks under `tst/inputs/`.
3. Put user-facing examples under `inputs/io/`.
4. Use isolated run directories.
5. Keep frozen baseline fixtures under
   `tst/fixtures/io/origin_main_886dd2a1/`.
6. Do not import the remote `tst/scripts/` tests verbatim.

### Scenarios To Port From The Remote Branch

- [ ] Forced-small MPI byte chunks using `ATHENAK_TEST_MAX_MPI_BYTES`.
- [ ] Restart from public manifest.
- [ ] Restart from normalized node-payload path if supported.
- [ ] Relative and absolute restart path handling.
- [ ] Stale node-payload collision beside a valid shared restart.
- [ ] Detailed restart phase statistics if retained.
- [ ] Detailed output stats for skipped empty node shards if retained.
- [ ] Empty-node sliced binary behavior.
- [ ] Converter convenience API behavior.

### Existing Local Scenarios To Preserve

- [ ] Legacy one-dimensional PDF byte compatibility.
- [ ] Legacy two-dimensional PDF byte compatibility.
- [ ] Modern `AKPDFV2` dense and sparse readback.
- [ ] Transitional unversioned dense and sparse PDF readback.
- [ ] Three- and four-dimensional PDFs.
- [ ] Linear, log, and symlog axes.
- [ ] Volume, mass, and variable weighting.
- [ ] Spherical-slice strict malformed-input rejection.
- [ ] Derived spherical-slice rejection.
- [ ] Two-fluid generic diagnostic rejection.
- [ ] Final-output policy and terminal-checkpoint numbering.
- [ ] Canonical converter CLI and example readback.

### Qualification Matrix

| Environment | Required evidence |
| --- | --- |
| Serial CPU | Builds, focused IO pytest modules, fixture checksums |
| MPI CPU on one node | Shared/rank/node equivalence, restart round trips, forced chunks |
| MPI CPU on multiple physical nodes | Node file counts, empty-node behavior, direct restart reads, restart numbering |
| GPU-capable build | Existing IO GPU regression plus any new derived-variable paths |
| Docs environment | Temporary `origin/gh-pages` worktree Sphinx build with warnings as errors |
| Style gate | Repository style suite |

### Stop Gate

- [ ] Every accepted behavior has a positive test.
- [ ] Every parser or manifest validation rule has a negative test.
- [ ] GPU qualification is run, not merely defined.
- [ ] Real multi-node qualification is run, not inferred from two local ranks.

## CP-07: Reconcile Examples And Deferred Pages Content

### Purpose

Keep documentation merge-ready without publishing behavior prematurely.

### Required Actions

1. Update examples after the implementation stabilizes.
2. Update deferred Pages content under
   `deferred_docs/gh-pages/io-output-formats-and-sharding/`.
3. Ensure restart docs accurately describe:
   - transactional payload and manifest publication;
   - supported restart entry points;
   - native direct loading after CP-03;
   - strict validation;
   - no `.assembled` staging in the final implementation.
4. Ensure PDF docs describe:
   - legacy text compatibility;
   - `AKPDFV2`;
   - transitional unversioned read compatibility;
   - dense shared and sparse rank/node layouts.
5. Ensure converter docs expose only canonical `bin_convert.py`.
6. Keep sliced node-sharded `.cbin` explicitly outside the promoted workflow.
7. Apply the deferred overlay only in a temporary detached Pages worktree.
8. Run:

   ```bash
   make clean html SPHINXOPTS="-W --keep-going"
   ```

### Subagent Assignment

Spawn one read-only docs-to-code auditor. Require it to compare every documented
parameter, filename, reader entry point, and qualification claim against the
current implementation and tests.

### Stop Gate

- [ ] Deferred docs match code.
- [ ] Sphinx warnings-as-errors build passes in a temporary Pages worktree.
- [ ] Live `gh-pages` remains untouched.

## CP-08: Final Independent Audits And Commit Preparation

### Purpose

Confirm that the final branch contains the intended IO feature and no accidental
remote-branch baggage.

### Required Independent Audits

Spawn separate read-only auditors for:

1. C++ file formats and output registration.
2. MPI collectives, chunking, node communicators, and restart loading.
3. Python API compatibility and malformed-file handling.
4. Test completeness, examples, and fixture integrity.
5. Deferred Pages accuracy and publication instructions.
6. Final scope review against `origin/main`,
   `origin/gotham-1.0`, and
   `origin/feature/single-file-per-node-outputs`.

### Final Scope Review Questions

1. Is any remote binary restart-manifest assumption still present?
2. Does any production node-restart path create `.assembled`?
3. Are all node communicators initialized lazily?
4. Does any supported code or documentation mention `bin_convert_new.py` as a
   public API?
5. Are legacy PDF fixtures still byte-identical?
6. Are transitional remote PDF payloads still readable?
7. Are shared and per-rank restart paths preserved?
8. Are terminal checkpoint counters advanced normally?
9. Are detailed stats opt-in?
10. Are sliced node-sharded `.cbin` claims still excluded?
11. Are unrelated remote changes absent?
12. Do decision-log and audit-ledger entries match the final implementation?

### Final Verification Commands

Record exact commands and results in `IO_FEATURE_AUDIT_LEDGER.md`.

```bash
git status --short --branch
git diff --check origin/main
shasum -a 256 -c tst/fixtures/io/origin_main_886dd2a1/SHA256SUMS
rg -n "bin_convert_new|\\.assembled" \
  IO_FEATURE_BRANCH_GUIDE.md \
  IO_FEATURE_BRANCH_INTEGRATION_PLAN.md \
  IO_FEATURE_BRANCH_DECISION_LOG.md \
  IO_FEATURE_AUDIT_LEDGER.md \
  IO_FORMAT_COMPATIBILITY.md \
  deferred_docs src tst vis inputs
```

Also run the repository's serial CPU, MPI CPU, GPU, docs, and style commands
recorded during earlier checkpoints.

### Stop Gate

- [ ] All blocking findings resolved.
- [ ] Remaining scope exclusions are explicit.
- [ ] GPU qualification passed.
- [ ] Real multi-node qualification passed.
- [ ] Deferred Pages validation passed.
- [ ] Final independent audit reports are recorded.
- [ ] Commit sequence is coherent and reviewable.

## Recommended Final Commit Shape

After the reconstructed baseline is preserved, use narrow commits for the
selective integration:

1. `io: add chunked MPI byte operations`
2. `restart: load transactional node payloads directly`
3. `outputs: harden empty node shards and spherical-slice publication`
4. `tools: restore canonical binary converter convenience APIs`
5. `tests: cover chunked restart IO and node-shard edge cases`
6. `docs: reconcile IO contracts and deferred Pages overlay`

Keep behavior and its focused tests together where practical. Avoid one large
integration commit.

## Definition Of Done

The branch is complete only when:

1. The local `AKPDFV2`, legacy PDF, spherical-slice, timing, final-policy, and
   transactional publication contracts are preserved.
2. Per-node restart loading reads node payloads directly without rank-0 shared
   staging.
3. Chunked MPI byte IO is used for large transfers and exercised with forced
   tiny chunks.
4. Shared, per-rank, and per-node workflows remain test-backed.
5. One canonical `bin_convert.py` contains required legacy and useful modern
   behavior.
6. Tests cover valid, malformed, empty, truncated, and resumed cases.
7. GPU and real multi-node qualification have actually run.
8. Deferred Pages content matches the final implementation and builds cleanly.
9. `IO_FEATURE_AUDIT_LEDGER.md` contains the evidence.
10. `IO_FEATURE_BRANCH_DECISION_LOG.md` contains the final accepted, rejected,
    pending, and superseded decisions.
