# IO Output Formats And Sharding Decision Log

## Purpose

This file records consequential decisions for
`feature/io-output-formats-and-sharding`. It is the durable record for why the
combined IO branch keeps some behavior, ports some behavior selectively, and
rejects other behavior.

Update this file whenever:

- a public API changes;
- a file layout changes;
- a compatibility promise changes;
- a remote-branch behavior is accepted or rejected;
- a test or qualification requirement is added or removed;
- a pending question is resolved;
- an earlier decision is superseded.

Do not silently rewrite history. Add a new entry that names the superseded
decision.

## Status Values

| Status | Meaning |
| --- | --- |
| Accepted | The branch should implement and preserve this decision. |
| Rejected | The branch should not implement this behavior. |
| Pending evidence | More design work, consumer evidence, or qualification is required. |
| Superseded | A later decision replaces this one. |

## Decision Summary

| ID | Status | Decision | Implementation checkpoint |
| --- | --- | --- | --- |
| D-001 | Accepted | Use `feature/io-output-formats-and-sharding` as the branch name | CP-01 |
| D-002 | Accepted | Treat the local reconstruction as the integration baseline | CP-00 |
| D-003 | Rejected | Do not merge or wholesale cherry-pick `origin/feature/single-file-per-node-outputs` | All |
| D-004 | Accepted | Preserve local versioned `AKPDFV2` output | Baseline |
| D-005 | Accepted | Preserve byte-compatible legacy one- and two-dimensional text PDFs | Baseline |
| D-006 | Accepted | Keep read compatibility for transitional unversioned remote PDF payloads | Baseline |
| D-007 | Accepted | Preserve the local transactional text restart manifest and payload publication contract | CP-03 |
| D-008 | Accepted | Replace `.assembled` node-restart staging with native distributed reads | CP-03 |
| D-009 | Rejected | Do not adopt the remote binary restart manifest schema | CP-03 |
| D-010 | Accepted | Port generalized chunked MPI byte IO behavior | CP-02 |
| D-011 | Accepted | Keep node communicator setup lazy and opt-in | CP-02, CP-03, CP-04 |
| D-012 | Accepted | Preserve opt-in high-level output timing | Baseline |
| D-013 | Accepted | Preserve explicit final-output policy with default `all` | Baseline |
| D-014 | Accepted | Advance counters normally for every written final output | Baseline |
| D-015 | Accepted | Maintain one public `vis/python/bin_convert.py` | CP-05 |
| D-016 | Accepted | Add `read_rank_binary_as_athdf()` to canonical `bin_convert.py` | CP-05 |
| D-017 | Pending evidence | Add indexed single-meshblock reads only if existing call forms remain unambiguous | CP-05 |
| D-018 | Rejected | Do not restore `bin_convert_new.py` as a public module or shim | CP-05 |
| D-019 | Pending evidence | Do not add the remote `athinput()` helper unless a real consumer is identified | CP-05 |
| D-020 | Accepted | Port remote test scenarios into pytest rather than copying remote `tst/scripts/` files | CP-06 |
| D-021 | Accepted | Keep sliced node-sharded `.cbin` outside the promoted workflow | CP-04, CP-06 |
| D-022 | Accepted | Require actual GPU qualification before merge readiness | CP-06, CP-08 |
| D-023 | Accepted | Require actual multi-node MPI qualification before production readiness | CP-06, CP-08 |
| D-024 | Accepted | Stage documentation for later `gh-pages` integration without publishing from this branch | CP-07 |
| D-025 | Pending evidence | Decide whether detailed per-phase IO stats should be retained as a separate opt-in mode | CP-04 |
| D-026 | Pending evidence | Decide whether a node payload path is a supported restart entry point or only a tested normalization convenience | CP-03 |
| D-027 | Rejected | Do not auto-detect stale node payloads beside an ordinary shared restart and silently reinterpret the shared file | CP-03 |
| D-028 | Rejected | Do not port unrelated remote `imex2+` messaging or mesh-diagnostics refactors as part of this IO integration | CP-08 |
| D-029 | Accepted | Preserve the reconstructed baseline in three reviewable commits before selective hardening | CP-00 |
| D-030 | Accepted | Mark immutable legacy PDF fixtures as whitespace-insensitive in `.gitattributes` without changing their frozen bytes | CP-00 |
| D-031 | Accepted | Until CP-03 lands, describe node restart loading as strict manifest validation followed by transient rank-0 `.assembled` staging; support manifest-path restart only | CP-01, CP-03 |
| D-032 | Accepted | Port MPI chunking by intent with checked multiplication, checked offsets, and strict collective symmetry | CP-02 |
| D-033 | Accepted | Keep zero-byte collective participants in the communicator-wide chunk schedule with dummy buffers | CP-02 |
| D-034 | Accepted | Truncate MPI files once per communicator, check deletion errors, and synchronize before collective open | CP-02 |
| D-035 | Rejected | Do not couple `ParameterInput` to `FileShardMode` or change unrelated parameter-header limits during CP-02 | CP-02 |
| D-036 | Accepted | Require forced-small-chunk wrapper and shared-restart tests before CP-03; repeat for native node restart after CP-03 | CP-02, CP-03 |

## Detailed Decisions

### D-001: Branch Name

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Use `feature/io-output-formats-and-sharding`. |
| Reason | The user explicitly requested the renamed feature branch. Existing planning files contained stale references to the earlier prefixed name and required correction. |
| Evidence | Current `git status --short --branch` reports `feature/io-output-formats-and-sharding`. |
| Follow-up | Replace stale references during CP-01. |

### D-002: Integration Baseline

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Use the current local reconstruction as the authoritative integration baseline. |
| Reason | It has stricter file-format validation, compatibility fixtures, transactional restart publication, explicit runtime policy, promoted examples, and deferred Pages content. |
| Evidence | Fresh effective-tree comparison against `origin/feature/single-file-per-node-outputs`. |
| Follow-up | Preserve the baseline in logical commits during CP-00. |

### D-003: No Direct Remote Merge

| Field | Value |
| --- | --- |
| Status | Rejected |
| Decision | Do not merge or cherry-pick the remote extraction branch wholesale. |
| Reason | The branches contain competing restart formats, different communicator-lifecycle assumptions, different PDF schemas, different readers, and unrelated changes. |
| Evidence | Fresh comparison found 50 local-only files, 12 remote-only files, and 29 shared files with differing contents. |
| Follow-up | Port behavior selectively under the relevant checkpoints. |

### D-004 Through D-006: PDF Compatibility

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Keep local `AKPDFV2` writes, legacy text byte compatibility, and read compatibility for transitional unversioned dense/sparse payloads. |
| Reason | The local writer has an explicit version and stronger schema. The local reader already consumes remote experimental PDF artifacts without requiring the older writer. |
| Evidence | `src/outputs/pdf.cpp`; `vis/python/read_pdf.py`; frozen fixture tests. |
| Follow-up | Preserve coverage during CP-06 and document all three read paths during CP-07. |

### D-007 Through D-009: Restart Format And Loader

| Field | Value |
| --- | --- |
| Status | Accepted for D-007 and D-008; rejected for D-009 |
| Decision | Keep the local transactional text manifest and replace rank-0 `.assembled` staging with native distributed reads. Do not adopt the remote binary manifest. |
| Reason | Local publication and validation are stronger. Remote direct reads solve a real scalability problem. The right integration keeps both good properties without combining incompatible schemas. |
| Evidence | `src/outputs/restart.cpp`; `src/main.cpp`; remote `src/pgen/pgen.cpp`; remote `src/mesh/build_tree.cpp`. |
| Follow-up | Execute CP-03 with pre- and post-implementation restart audits. |

### D-010 And D-011: MPI IO Foundation

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Port generalized chunked MPI byte IO while keeping node setup lazy. |
| Reason | Large transfers should not depend on an MPI `int` count. Eager communicator setup changes ordinary runs unnecessarily. |
| Evidence | Remote `src/outputs/io_wrapper.cpp`; local `src/globals.cpp`; local `src/outputs/outputs.cpp`. |
| Follow-up | Implement CP-02 before the native restart loader. |

### D-012 Through D-014: Runtime Output Policy

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Keep opt-in output timing, `final_output_policy = all|restart_only|none`, default `all`, and normal counter advancement for every final write. |
| Reason | This exposes production controls without silently changing baseline behavior or allowing a resumed run to overwrite a terminal checkpoint. |
| Evidence | `src/driver/driver.cpp`; existing CPU and MPI policy regressions. |
| Follow-up | Preserve tests throughout all checkpoints. |

### D-015 Through D-019: Canonical Converter

| Field | Value |
| --- | --- |
| Status | Accepted except D-017 and D-019 pending evidence |
| Decision | Keep one `bin_convert.py`, add `read_rank_binary_as_athdf()`, consider additive indexed single-block reads, and omit `athinput()` unless a consumer is identified. |
| Reason | The branch should remove redundancy without losing useful modern workflows or breaking established conversion helpers. |
| Evidence | Local and remote `vis/python/bin_convert.py`; existing frozen conversion tests. |
| Follow-up | Execute CP-05 with API inventory and compatibility tests. |

### D-020 Through D-024: Tests, Qualification, And Documentation

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Port remote scenarios into pytest, keep sliced node-sharded `.cbin` excluded, require real GPU and multi-node runs, and stage Pages docs without publishing. |
| Reason | Repository-native tests are maintainable; local one-node MPI is insufficient for production node-shard claims; Pages content must not advertise code before merge. |
| Evidence | Existing `tst/test_suite/io/`; remote `tst/scripts/`; deferred Pages bundle. |
| Follow-up | Execute CP-06 and CP-07. |

### D-025: Detailed Per-Phase IO Statistics

| Field | Value |
| --- | --- |
| Status | Pending evidence |
| Decision | Decide whether to retain remote-style detailed phase statistics as a separate opt-in diagnostic mode. |
| Reason | They are useful for scaling analysis but overlap partially with the existing high-level `<time>/output_timing` records. The final interface should avoid redundant or noisy default logging. |
| Evidence needed | Proposed interface, example output, tests proving opt-in behavior, and documentation review. |
| Follow-up | Resolve during CP-04. |

### D-026 And D-027: Restart Entry Points

| Field | Value |
| --- | --- |
| Status | D-026 pending evidence; D-027 rejected |
| Decision | Consider strict normalization from a node payload path to its public manifest. Do not silently reinterpret an ordinary shared restart merely because stale node payloads exist nearby. |
| Reason | Explicit convenience can be safe. Ambient filesystem inference is risky. |
| Evidence needed | Component-based normalization design, negative tests, and stale-collision tests. |
| Follow-up | Resolve during CP-03. |

### D-028: Excluded Remote Drift

| Field | Value |
| --- | --- |
| Status | Rejected |
| Decision | Exclude unrelated remote messaging and mesh-diagnostics refactors from the IO integration. |
| Reason | They are not required by the IO feature and increase review surface. |
| Evidence | Remote diff includes an `imex2+` driver-message change and mesh-diagnostics cleanup unrelated to IO correctness. |
| Follow-up | Final scope auditor must confirm their absence during CP-08. |

### D-029: Preserve The Reconstructed Baseline

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Preserve the existing reconstruction before selective hardening in three reviewable commits: C++ IO implementation, Python readers/examples, and regression fixtures/tests. |
| Reason | The local reconstruction already passed its serial, MPI, fixture-integrity, and whitespace gates. A preserved baseline makes the subsequent chunking and restart-loader changes auditable. |
| Evidence | Commits `ca00581b`, `3049aa92`, and `fcc534fe`; fresh CP-00 verification recorded in `IO_FEATURE_AUDIT_LEDGER.md`. |
| Follow-up | Keep CP-02 and later hardening in narrow follow-up commits. |

### D-030: Preserve Intentional Fixture Whitespace

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Add a narrow `.gitattributes` rule for `tst/fixtures/io/origin_main_886dd2a1/pdf/**` with `-whitespace`. |
| Reason | Legacy PDF artifacts are frozen byte-for-byte compatibility fixtures. Editing their trailing spaces or final blank lines would invalidate the compatibility oracle, while treating their intentional bytes as source-format whitespace defects would make `git diff --check` noisy. |
| Evidence | `shasum -a 256 -c SHA256SUMS` passes from the fixture directory; `git diff --check origin/main` passes with the scoped attribute. |
| Follow-up | Keep the rule scoped to immutable legacy PDF fixtures. |

### D-031: Describe The Current Restart Loader Truthfully

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Until CP-03 lands, describe node-restart input as strict public-manifest validation followed by transient rank-0 `<manifest>.assembled` staging into the legacy shared-file loader. Support the public manifest path only. |
| Reason | Transactional publication and strict manifest validation have landed, but production restart loading is not yet native or distributed. Direct payload entry is not currently implemented. |
| Evidence | `src/main.cpp`; CP-01 documentation consistency audit; corrected `IO_FORMAT_COMPATIBILITY.md`. |
| Follow-up | Supersede this decision when CP-03 removes staging and the accepted entry-point behavior is tested. |

### D-032 Through D-036: Chunked MPI IO Contract

| Field | Value |
| --- | --- |
| Status | Accepted for D-032, D-033, D-034, and D-036; rejected for D-035 |
| Decision | Port chunking by intent with overflow-safe byte math, checked MPI offsets, communicator-wide collective chunk schedules, dummy buffers for zero-byte participants, synchronized one-rank file truncation, and forced-small-chunk tests. Do not copy unrelated `ParameterInput` coupling or header-limit changes. |
| Reason | Existing MPI operations narrow 64-bit sizes to `int`, current shared-file truncation is racy, and the reference collective helpers can deadlock if a rank returns early while peers continue. |
| Evidence | CP-02 preimplementation MPI audit of local and reference `src/outputs/io_wrapper.*`, `src/mesh/build_tree.cpp`, and reference restart tests. |
| Follow-up | Implement and independently audit CP-02 before starting native restart reads. |

## Pending Decision Queue

Resolve these before merge readiness:

| ID | Question | Required evidence | Owner checkpoint |
| --- | --- | --- | --- |
| D-017 | Can indexed single-meshblock reads be added without ambiguous legacy positional calls? | Signature design and compatibility tests | CP-05 |
| D-019 | Does any real consumer need `athinput()` in canonical `bin_convert.py`? | Repository import sweep or identified external requirement | CP-05 |
| D-025 | Should detailed per-phase stats be retained, and under which opt-in interface? | Interface proposal, tests, docs audit | CP-04 |
| D-026 | Should users be allowed to restart from an individual node payload path? | Strict normalization design and negative tests | CP-03 |

## Decision Update Template

Append new decisions in this form:

```markdown
### D-XXX: Short Title

| Field | Value |
| --- | --- |
| Status | Accepted, Rejected, Pending evidence, or Superseded |
| Decision | Concrete behavior or scope choice |
| Reason | Why this choice is correct |
| Evidence | Files, tests, audit report, or qualification run |
| Supersedes | Decision ID if applicable |
| Follow-up | Checkpoint and required action |
```
