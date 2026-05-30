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

Review entry point: the compact table below summarizes the original
selective-integration decisions through `D-067`. It is not a complete summary
of the active branch contract. `D-068` onward are controlling post-review
robustification decisions: read their full entries under **Detailed Decisions**
below and use the current **Robustification Pass** section of
`IO_FEATURE_AUDIT_LEDGER.md` for checkpoint status. Active execution is governed
by `IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md`.

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
| D-017 | Superseded by D-040 | Add indexed single-meshblock reads only if existing call forms remain unambiguous | CP-05 |
| D-018 | Rejected | Do not restore `bin_convert_new.py` as a public module or shim | CP-05 |
| D-019 | Superseded by D-040 | Do not add the remote `athinput()` helper unless a real consumer is identified | CP-05 |
| D-020 | Accepted | Port remote test scenarios into pytest rather than copying remote `tst/scripts/` files | CP-06 |
| D-021 | Superseded by D-055 | Keep sliced node-sharded `.cbin` outside the promoted workflow | CP-04, CP-06 |
| D-022 | Accepted | Require actual GPU qualification before merge readiness | CP-06, CP-08 |
| D-023 | Accepted | Require actual multi-node MPI qualification before production readiness | CP-06, CP-08 |
| D-024 | Accepted | Stage documentation for later `gh-pages` integration without publishing from this branch | CP-07 |
| D-025 | Superseded by D-039 | Decide whether detailed per-phase IO stats should be retained as a separate opt-in mode | CP-04 |
| D-026 | Superseded by D-037 | Decide whether a node payload path is a supported restart entry point or only a tested normalization convenience | CP-03 |
| D-027 | Rejected | Do not auto-detect stale node payloads beside an ordinary shared restart and silently reinterpret the shared file | CP-03 |
| D-028 | Rejected | Do not port unrelated remote `imex2+` messaging or mesh-diagnostics refactors as part of this IO integration | CP-08 |
| D-029 | Accepted | Preserve the reconstructed baseline in three reviewable commits before selective hardening | CP-00 |
| D-030 | Accepted | Mark immutable legacy PDF fixtures as whitespace-insensitive in `.gitattributes` without changing their frozen bytes | CP-00 |
| D-031 | Superseded by D-037 | Until CP-03 lands, describe node restart loading as strict manifest validation followed by transient rank-0 `.assembled` staging; support manifest-path restart only | CP-01, CP-03 |
| D-032 | Accepted | Port MPI chunking by intent with checked multiplication, checked offsets, and strict collective symmetry | CP-02 |
| D-033 | Accepted | Keep zero-byte collective participants in the communicator-wide chunk schedule with dummy buffers | CP-02 |
| D-034 | Accepted | Truncate MPI files once per communicator, check deletion errors, and synchronize before collective open | CP-02 |
| D-035 | Rejected | Do not couple `ParameterInput` to `FileShardMode` or change unrelated parameter-header limits during CP-02 | CP-02 |
| D-036 | Accepted | Require forced-small-chunk wrapper and shared-restart tests before CP-03; repeat for native node restart after CP-03 | CP-02, CP-03 |
| D-037 | Accepted and verified | Use native node-restart reads through the public manifest only | CP-03 |
| D-038 | Accepted and verified | Publish explicit valid empty node shards | CP-04 |
| D-039 | Rejected | Do not port detailed per-phase IO statistics into this branch | CP-04 |
| D-040 | Accepted and verified | Add canonical converter helpers additively and omit `athinput()` | CP-05 |
| D-041 | Accepted and verified | Validate full positioned MPI byte ranges before collective IO | CP-02 |
| D-042 | Accepted and verified | Canonicalize node-restart payload paths, reject symlink escapes, and compare every replicated payload header byte-for-byte against payload 0 | CP-03 |
| D-043 | Accepted and verified | Add optional node-inventory metadata to new shards, validate it when present, and atomically publish spherical-slice files | CP-04 |
| D-044 | Accepted and verified | Reject generated node-restart payload artifacts through lexical, repeated-separator, symlink, and temporary-file aliases | CP-03 |
| D-045 | Accepted and verified | Reject oversized declared shard inventories before allocation, aggregate assembled node MeshBlock counts, and reject sliced node `.cbin` during construction | CP-04 |
| D-046 | Accepted with deferred optimization | Retain strict manifest validation on every MPI rank for this branch; defer central validation plus structured broadcast as a scaling optimization | CP-03, CP-08 |
| D-047 | Accepted and verified | Mark node-restart payload bytes explicitly, bound declared payload inventory, and check payload publication writes and closes | CP-03 |
| D-048 | Accepted and verified | Publish modern PDF files atomically, record sparse shard inventory, and bound Python reader allocations | CP-04, CP-05 |
| D-049 | Accepted and verified | Require declared PDF V2 payloads, normalize reconstructed metadata, bound reader bulk reads, and execute frozen restart fixtures | CP-05, CP-06 |
| D-050 | Accepted and verified | Bound restart segment inventory and reject non-positive segment records | CP-03 |
| D-051 | Accepted and verified | Complete malformed-reader preflight before aggregate materialization | CP-05 |
| D-052 | Accepted and verified | Ignore generated Python bytecode and pytest caches | CP-06 |
| D-053 | Accepted and verified | Bound retained aggregates and athdf-like conversion allocations | CP-05 |
| D-054 | Accepted and verified | Restore bounded athdf-like reconstruction and enforce complete producer/reader metadata contracts | CP-05, CP-06 |
| D-055 | Accepted and verified | Reject sliced `cbin` consistently and validate emitted extents during writer construction | CP-04, CP-06 |
| D-056 | Accepted and verified | Reject non-finite PDF axis bounds before edge construction | CP-04, CP-06 |
| D-057 | Accepted and verified | Complete bounded Python reconstruction for ghost zones, cropped cells, malformed grids, partial shards, and reader live peaks | CP-05, CP-06 |
| D-058 | Accepted and verified | Detect omitted ghost counts and bound remaining Python reconstruction temporaries | CP-05, CP-06 |
| D-059 | Accepted and verified | Release incorporated sparse sibling state and bound remaining metadata-validation temporaries | CP-05, CP-06 |
| D-060 | Accepted and verified | Exercise retained ATHDF orchestration APIs and promoted node `cbin` helper workflow positively | CP-06 |
| D-061 | Accepted and verified | Require exact logical geometry and preflight PDF edge materialization | CP-05, CP-06 |
| D-062 | Accepted and verified | Use bounded absolute geometry tolerance and strictly bound legacy PDF ASCII expansion | CP-05, CP-06 |
| D-063 | Accepted and verified | Freeze generated-edge guard ordering with a fail-fast regression | CP-06 |
| D-064 | Accepted and verified | Forward retained legacy-header budgets and bound spherical metadata tokenization | CP-05, CP-06 |
| D-065 | Accepted and verified | Enforce ASCII PDF tokens and retain spherical sibling metadata budgets | CP-05, CP-06 |
| D-066 | Accepted and verified | Prove public spherical sibling-budget forwarding in automation | CP-06 |
| D-067 | Accepted and verified | Retain spherical metadata through coordinate generation | CP-05, CP-06 |
| D-068 onward | See detailed decisions and current robustification ledger | Controlling post-review robustification decisions are recorded in full below and summarized by the current `IO_FEATURE_AUDIT_LEDGER.md` robustification checkpoint board | RCP-00 onward |

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
| Status | Accepted; D-017 and D-019 resolved by D-040 |
| Decision | Keep one `bin_convert.py`, add `read_rank_binary_as_athdf()`, add indexed single-block reads only through a keyword-only selector, and omit `athinput()`. |
| Reason | The branch should remove redundancy without losing useful modern workflows or breaking established conversion helpers. |
| Evidence | Local and remote `vis/python/bin_convert.py`; existing frozen conversion tests. |
| Follow-up | Completed by D-040 with API inventory and compatibility tests. |

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
| Status | Superseded by D-039 |
| Decision | Decide whether to retain remote-style detailed phase statistics as a separate opt-in diagnostic mode. |
| Reason | They are useful for scaling analysis but overlap partially with the existing high-level `<time>/output_timing` records. The final interface should avoid redundant or noisy default logging. |
| Evidence needed | Proposed interface, example output, tests proving opt-in behavior, and documentation review. |
| Follow-up | Resolved by D-039. |

### D-026 And D-027: Restart Entry Points

| Field | Value |
| --- | --- |
| Status | D-026 superseded by D-037; D-027 rejected |
| Decision | Consider strict normalization from a node payload path to its public manifest. Do not silently reinterpret an ordinary shared restart merely because stale node payloads exist nearby. |
| Reason | Explicit convenience can be safe. Ambient filesystem inference is risky. |
| Evidence needed | Component-based normalization design, negative tests, and stale-collision tests. |
| Follow-up | Resolved by D-037. |

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
| Status | Superseded by D-037 |
| Decision | Until CP-03 lands, describe node-restart input as strict public-manifest validation followed by transient rank-0 `<manifest>.assembled` staging into the legacy shared-file loader. Support the public manifest path only. |
| Reason | Transactional publication and strict manifest validation have landed, but production restart loading is not yet native or distributed. Direct payload entry is not currently implemented. |
| Evidence | `src/main.cpp`; CP-01 documentation consistency audit; corrected `IO_FORMAT_COMPATIBILITY.md`. |
| Follow-up | Superseded by D-037 after CP-03 removed staging and tested the accepted entry-point behavior. |

### D-032 Through D-036: Chunked MPI IO Contract

| Field | Value |
| --- | --- |
| Status | Accepted for D-032, D-033, D-034, and D-036; rejected for D-035 |
| Decision | Port chunking by intent with overflow-safe byte math, checked MPI offsets, communicator-wide collective chunk schedules, dummy buffers for zero-byte participants, synchronized one-rank file truncation, and forced-small-chunk tests. Do not copy unrelated `ParameterInput` coupling or header-limit changes. |
| Reason | Existing MPI operations narrow 64-bit sizes to `int`, current shared-file truncation is racy, and the reference collective helpers can deadlock if a rank returns early while peers continue. |
| Evidence | CP-02 preimplementation MPI audit of local and reference `src/outputs/io_wrapper.*`, `src/mesh/build_tree.cpp`, and reference restart tests. |
| Follow-up | Implement and independently audit CP-02 before starting native restart reads. |

### D-037: Native Node Restart Uses The Public Manifest Only

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Replace transient `.assembled` staging with validated direct reads from generation-qualified node payloads. Keep the public manifest path as the only supported node-restart entry point; reject payload-path restart. |
| Reason | The public manifest is the transactional commit point. Payload paths are implementation artifacts, and normalizing them adds ambiguity without improving the supported workflow. |
| Evidence | CP-03 restart-design audit of local `src/main.cpp`, `src/outputs/restart.cpp`, and the reference direct-reader pattern. |
| Supersedes | D-026 and D-031 once CP-03 implementation passes its stop gate. |
| Follow-up | Extract a structured manifest module, remove production `.assembled` paths, and test direct node-payload reads during CP-03. |

### D-038: Explicit Empty Node Shards

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Publish explicit valid empty node shards for binary, coarsened-binary, and spherical-slice output. Do not skip empty nodes. |
| Reason | Explicit empty shards keep inventories auditable, preserve counter semantics, and avoid silent omissions. The local readers already have a natural header-only or zero-point representation. |
| Evidence | CP-04 preimplementation audit of `src/outputs/binary.cpp`, `src/outputs/coarsened_binary.cpp`, `src/outputs/spherical_slice.cpp`, `vis/python/bin_convert.py`, and `vis/python/read_sphslice.py`. |
| Follow-up | Add writer hardening, sibling-inventory validation, and scheduler-backed empty-node qualification during CP-04 and CP-06. |

### D-039: Reject Detailed Per-Phase IO Statistics In This Branch

| Field | Value |
| --- | --- |
| Status | Rejected |
| Decision | Do not port the reference `ATHENAK_OUTPUT_IO_STATS` phase diagnostics. Preserve the existing opt-in `<time>/output_timing=true` event timing interface. |
| Reason | The reference interface overlaps the existing slowest-rank timing records, uses a separate environment-variable control surface, and adds format-specific output noise. |
| Evidence | CP-04 preimplementation audit of local `src/driver/driver.cpp` and the reference statistics paths. |
| Supersedes | D-025. |
| Follow-up | Treat detailed distributed-read phase statistics as a separate scaling-analysis feature if future evidence justifies them. |

### D-040: Canonical Python Rank Reader Must Delegate

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Add `read_rank_binary_as_athdf()` as a thin delegate to the canonical logical-location mapper. Add indexed single-meshblock reads only as keyword-only `meshblock_index_in_file=...`. Omit `athinput()`. |
| Reason | The reference rank-reader scatter implementation misplaces later logical MeshBlocks at the root-grid origin. Its positional index overload also changes the meaning of existing calls. No target-tree consumer needs `athinput()`. |
| Evidence | CP-05 Python API audit, including a rank-1 fixture probe and repository import sweep. |
| Supersedes | Resolves D-017 and D-019. |
| Follow-up | Implement the additive APIs, harden malformed-file handling, and add signature and placement regressions during CP-05. |

### D-041: Validate Full Positioned MPI Byte Ranges

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Before any non-empty positioned MPI read or write, validate the inclusive end offset `offset + total_bytes - 1` as well as the starting offset. Add a direct MPI wrapper harness for asymmetric zero-byte collective reads and writes, range overflow, multiplication overflow, and communicator-rank chunk-limit disagreement. |
| Reason | The first CP-02 implementation checked each chunk start but could still pass a representable start plus an unrepresentable non-empty byte range into MPI. Existing application-level pytest coverage did not directly exercise that boundary. |
| Evidence | CP-02 independent post-integration MPI audit of `src/outputs/io_wrapper.cpp` and `tst/test_suite/io/test_chunked_io_mpicpu.py`. |
| Supersedes | Tightens D-032 and D-036. |
| Follow-up | Preserve the corrected wrapper checkpoint and use it for CP-03 native node restart reads. |

### D-042: Validate Node-Restart Payload Containment And Replicated Headers

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Canonicalize the manifest directory and every generation-qualified payload path, reject payload symlinks that escape the manifest directory, use payload 0 as the canonical header stream, and compare every replicated payload header byte-for-byte against payload 0 before accepting the manifest. |
| Reason | The public manifest is only a valid transactional commit point when its payload inventory cannot escape the checkpoint tree and all replicated headers describe the same restart state. |
| Evidence | CP-03 implementation in `src/restart_manifest.cpp`; expanded MPI regressions for symlink escape and replicated-header mismatch; local native node-restart suite returned `32 passed`. |
| Supersedes | Tightens D-007, D-008, and D-037. |
| Follow-up | Preserve the validated parser and repeat cross-node qualification on multiple physical nodes. |

### D-043: Validate Optional Node Inventories And Publish Spherical Slices Atomically

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Add `number of nodes` metadata to newly written node-sharded binary and coarsened-binary files; validate complete, unique sibling inventories when additive count metadata is present while continuing to accept older files without it; publish spherical slices through checked temporary files and atomic rename. |
| Reason | Explicit inventories make empty shards auditable without breaking older files, while atomic spherical-slice publication prevents readers from observing partial output. |
| Evidence | CP-04 implementation in `src/outputs/binary.cpp`, `src/outputs/coarsened_binary.cpp`, `src/outputs/spherical_slice.cpp`, `vis/python/bin_convert.py`, and `vis/python/read_sphslice.py`; focused local CPU reader matrix returned `51 passed`; focused MPI writer matrix returned `2 passed`. |
| Supersedes | Tightens D-038. |
| Follow-up | Repeat the explicit-empty-node qualification on multiple physical nodes. |

### D-044: Reject Generated Restart Payload Aliases

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Classify generated node-restart payload artifacts after lexical normalization and, when the path exists, canonical symlink resolution. Reject exact payload paths, `./` aliases, repeated separators, symlink aliases, and unpublished `.tmp` payload artifacts as direct restart entry points. Preserve ordinary shared restart filenames that merely end in `.payload.rst`. |
| Reason | The public manifest is the transactional commit point. An alias that opens a generated payload directly bypasses inventory and completion validation even when the exact spelling is rejected. |
| Evidence | Post-integration restart audit; `src/restart_manifest.cpp`; expanded alias matrix in `test_node_sharding_mpicpu.py`; focused native restart suite returned `36 passed`. |
| Supersedes | Tightens D-037 and D-042. |
| Follow-up | Preserve the manifest-only contract in deferred Pages content and repeat it during multi-node qualification. |

### D-045: Bound Reader Inventory Validation And Reject Sliced Node Cbin Early

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Compare declared shard counts with discovered sibling counts before materializing expected-ID ranges; aggregate `number_of_meshblocks` in assembled node binary/coarsened-binary objects; and reject sliced node-sharded `.cbin` during construction before output-data loading begins. |
| Reason | File-controlled counts must fail cheaply, assembled metadata must describe the assembled object, and the deliberately excluded sliced-node-`.cbin` path must fail deterministically even with one MPI rank. |
| Evidence | Post-integration writer/reader audit; `src/outputs/coarsened_binary.cpp`; `vis/python/bin_convert.py`; `vis/python/read_sphslice.py`; CPU reader-hardening suite returned `17 passed`; MPI writer-hardening suite returned `3 passed`. |
| Supersedes | Tightens D-021 and D-043. |
| Follow-up | Preserve the early rejection and bounded-validation regressions. |

### D-046: Defer Centralized Manifest Validation Broadcast

| Field | Value |
| --- | --- |
| Status | Accepted with deferred optimization |
| Decision | Keep strict node-restart manifest inventory and replicated-header validation independently on every MPI rank in this branch. Defer validate-once central parsing plus structured metadata broadcast to a future scaling-focused change. |
| Reason | The current path is correct and removes the dominant full-file `.assembled` bottleneck. Retrofitting a new distributed metadata protocol late in this integration would enlarge the restart blast radius without local multi-node hardware to qualify it. |
| Evidence | Post-integration restart audit of `src/main.cpp` and `src/restart_manifest.cpp`; local direct-read MPI regressions; explicit external multi-node qualification gate. |
| Follow-up | Measure validation overhead at production scale and implement centralized metadata distribution if it is material. |

### D-047: Mark And Bound Node-Restart Payloads

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Write an explicit node-payload marker after the replicated parameter dump, require it when loading through a node manifest, reject it when an ordinary restart path opens payload bytes directly, bound declared payload inventory at `1048576`, and check node-payload header writes and close returns before rename. |
| Reason | Pathname classification alone cannot distinguish hard links or copied payload bytes from ordinary shared restarts. The on-disk marker enforces the manifest-only contract by content while preserving legacy shared and per-rank restart bytes. Bounded inventory and checked publication prevent malformed manifests and partial payloads from consuming unbounded memory or becoming visible. |
| Evidence | Second independent restart audit; `src/main.cpp`; `src/outputs/restart.cpp`; `src/restart_manifest.*`; expanded hard-link, copied-payload, corrupt-marker, and oversized-count regressions; focused native restart suite returned `40 passed`. |
| Supersedes | Tightens D-007, D-037, D-042, and D-044. |
| Follow-up | Repeat the manifest-only marker contract during real multi-node qualification. |

### D-048: Publish And Read Modern PDF Shards Defensively

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Publish each modern PDF header and payload through a checked temporary file plus atomic rename; add `rank`/`number_of_ranks` or `node`/`number_of_nodes` header metadata for sparse shards; validate canonical shard paths and identifiers; and bound file-controlled dense allocations in the shipped Python readers. |
| Reason | Readers must not observe partially written modern files, sparse reconstruction must reject incomplete or aliased shard inventories, and malformed files must fail before allocating unreasonable memory. Atomicity applies independently to each PDF file rather than to the header/payload family as one filesystem transaction. |
| Evidence | Second independent writer/reader audits; `src/outputs/pdf.cpp`; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; `vis/python/bin_convert.py`; focused CPU reader matrix returned `70 passed`; focused MPI PDF suite returned `1 passed`. |
| Supersedes | Tightens D-004, D-006, D-043, and D-045. |
| Follow-up | Preserve the per-file publication wording and exercise rank/node sparse inventory during external multi-node qualification. |

### D-049: Close Final Reader And Fixture Audit Gaps

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Require `AKPDFV2` headers to pair with V2 payload preambles and matching cycles; clear shard-local rank/node metadata from reconstructed PDF and spherical-slice aggregates; set reconstructed spherical-slice `npoints` to full surface coverage; bound PDF metadata/payload reads, spherical-slice whole-file reads, and binary metadata plus aggregate payload construction; and resume frozen `origin/main` shared/per-rank restart fixtures in pytest. |
| Reason | Strict schemas must not silently downgrade corrupted V2 bytes to transitional parsing, reconstructed APIs must describe the aggregate rather than shard zero, file-controlled sizes must fail before unreasonable bulk reads or copies, and checksum-only restart fixtures do not prove loader compatibility. |
| Evidence | Final Python and test-coverage audits; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; `vis/python/bin_convert.py`; expanded CPU/MPI fixture, malformed-file, aggregate, and helper-CLI regressions; focused Python reader/writer matrix returned `83 passed`. |
| Supersedes | Tightens D-006, D-042, D-045, and D-048. |
| Follow-up | Preserve frozen-fixture resume coverage and repeat external GPU plus scheduler-backed multi-node qualification. |

### D-050: Bound Restart Segment Inventory

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Bound node-restart segment records before storing them, reject segment inventories larger than the declared MeshBlock count or the practical `1048576` ceiling, and reject non-positive segment counts. |
| Reason | A payload-count bound alone does not constrain malformed manifests: arbitrarily many zero-length segment records could otherwise consume memory and validation time independently on every MPI rank. Valid producer output never needs zero-length segments, and a positive segment inventory cannot contain more records than MeshBlocks. |
| Evidence | Final restart re-audit; `src/restart_manifest.cpp`; expanded malformed-manifest regressions in `test_node_sharding_mpicpu.py`; focused MPI restart/output/chunk/writer matrix returned `53 passed`. |
| Supersedes | Tightens D-047 and D-049. |
| Follow-up | Retain the bounded-segment negatives and include malformed-manifest handling in scheduler-backed multi-node qualification. |

### D-051: Complete Reader Preflight Before Aggregate Materialization

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Require sparse PDF siblings to agree on `binary_magic` and parse each payload against its local header; bound cumulative binary metadata records, positive binary variable counts, fixed MeshBlock metadata, and incremental shard payload/metadata totals; and preflight spherical-slice coordinate arrays plus embedded input-header offsets before reading or materializing aggregate objects. |
| Reason | Per-file limits are insufficient when malformed sibling metadata can downgrade validation, many individually bounded records can accumulate beyond a practical limit, metadata-only binary records can grow without variable payloads, or file-controlled coordinate and offset values reach allocation or read APIs unchecked. |
| Evidence | Final Python API audit; `vis/python/read_pdf.py`; `vis/python/bin_convert.py`; `vis/python/read_sphslice.py`; expanded malformed-reader regressions; focused Python reader/writer matrix returned `92 passed`; full serial IO matrix returned `114 passed`; full MPI IO matrix returned `59 passed`. |
| Supersedes | Tightens D-045, D-048, and D-049. |
| Follow-up | Preserve incremental preflight tests and repeat GPU qualification on a CUDA-capable build. |

### D-052: Ignore Generated Python Test Caches

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Add standard `.gitignore` rules for `__pycache__/`, `*.py[cod]`, and `.pytest_cache/`, then remove generated caches before staging. |
| Reason | Reader and style validation should not leave bytecode or pytest caches as untracked workspace noise during audit and commit preparation. |
| Evidence | Final scope audit; `.gitignore`; post-cleanup `git status --short`; repeated reader and style validation. |
| Follow-up | Keep generated test artifacts out of reviewable commits. |

### D-053: Bound Retained Aggregates And Conversion Allocations

| Field | Value |
| --- | --- |
| Status | Accepted and verified |
| Decision | Combine spherical-slice shards incrementally while retaining only small inventory summaries; retain only inventory summaries for PDF siblings; account for retained binary metadata and transient aggregate-copy peaks; reject non-positive coarsening factors; validate spherical embedded-header offsets in the public header-only API; and preflight coordinate, field, level-map, and restriction-map allocations in every athdf-like conversion helper. |
| Reason | Per-shard file bounds do not constrain retained sibling arrays, serialized metadata bytes undercount widened in-memory arrays, aggregate materialization temporarily coexists with source arrays, and conversion helpers allocate from parsed grid metadata after raw readback. |
| Evidence | Deep final Python API audit; `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; expanded reduced-limit regressions; focused Python reader/writer matrix returned `95 passed`; full serial IO matrix returned `117 passed`; full MPI IO matrix returned `59 passed`. |
| Supersedes | Tightens D-049 and D-051. |
| Follow-up | Preserve reduced-limit regressions and repeat CUDA-capable qualification before merge readiness. |

### D-054: Restore Bounded Reconstruction And Complete Metadata Contracts

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Centralize binary athdf-like MeshBlock placement; crop prolongation before materialization; restore subsampled, Cartesian, and volume-weighted restriction; validate bounded logical levels and non-negative locations before exponent arithmetic; preserve caller-provided single-rank arrays and dtype conversion; emit non-wrapping passive-scalar labels; reject duplicate binary labels; require finite PDF metadata and mandatory sparse V2 inventories; account for cumulative PDF/spherical retained arrays; and validate the documented `cbin` factor contract before writer construction or reader arithmetic. |
| Reason | The previous helper paths inherited placeholder restriction branches, could allocate complete prolonged MeshBlocks for narrow selections, and accepted malformed metadata that either hid distinct variables, published unreadable products, or reached arithmetic before contract validation. |
| Evidence | Final Python and format re-audits; `src/outputs/basetype_output.cpp`; `src/outputs/outputs.cpp`; `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; expanded reader and producer regressions; focused reader matrix returned `115 passed`; full serial IO matrix returned `146 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`. |
| Supersedes | Tightens D-051 and D-053. |
| Follow-up | Preserve malformed-metadata and reconstruction regressions; repeat CUDA-capable and scheduler-backed multi-node qualification before merge readiness. |

### D-055: Reject Incompatible Coarsened-Binary Extents During Construction

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Reject sliced `cbin` consistently in every shard mode and require every emitted full-volume extent, including optional ghost zones, to be divisible by the validated `coarsen_factor` during writer construction. |
| Reason | The prior node-only guard left shared and rank sliced output to fail through the pre-existing zero-width path, and a range-valid factor could still produce incompatible ghost-expanded extents. Construction-time validation gives one explicit contract before output-data loading. |
| Evidence | Final documentation audit; `src/outputs/coarsened_binary.cpp`; serial construction negatives for sliced and ghost-expanded extents; updated MPI sliced-output regression; full serial IO matrix returned `146 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`; detached Pages build passed with warnings treated as errors. |
| Supersedes | Tightens D-021, D-045, and D-054. |
| Follow-up | Keep sliced `cbin` outside the promoted workflow until a deliberately designed sliced coarsening contract exists. |

### D-056: Reject Non-Finite PDF Axis Bounds Before Edge Construction

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Require finite `pdf_bin_min` and `pdf_bin_max` values for every dimension before ordering checks, logarithmic-domain checks, edge construction, or PDF publication. |
| Reason | Ordered infinities passed the previous range check and could publish malformed PDF metadata containing non-finite edges. |
| Evidence | Final producer re-audit; `src/outputs/outputs.cpp`; runtime negative for `output1/bin4_max=inf`; focused producer module returned `19 passed`; full serial IO matrix returned `146 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`; auditor probe confirmed zero PDF files are published after rejection. |
| Supersedes | Tightens D-054. |
| Follow-up | Preserve the producer-side finite-bound regression alongside reader malformed-metadata coverage. |

### D-057: Complete Bounded Python Reconstruction

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Validate binary root-grid, MeshBlock, ghost-zone, logical-location, and geometry metadata before reconstruction; place emitted ghost-zone MeshBlocks by interior width and extend coordinates beyond root bounds; retain cells intersecting requested lower crop bounds; initialize uncovered partial-shard level cells to `-1`; preflight cumulative ATHDF coordinates, output arrays, and exact-restriction temporaries; preflight PDF bin-center, payload-copy, sparse-validation, and spherical-slice dump/diagnostic peaks; and reject non-shared distribution declarations in modern dense PDF headers. |
| Reason | Individually bounded arrays were still able to exceed practical live-memory caps when retained together, malformed grid metadata could silently produce empty or zero-filled products, ghost-zone placement overlapped at the wrong offsets, and lower-bound cropping discarded intersecting cells. |
| Evidence | Final Python tooling re-audit; `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; expanded reduced-cap, malformed-grid, crop, partial-level, and ghost-placement regressions; Python reader module returned `105 passed`; full serial IO matrix returned `161 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Supersedes | Tightens D-051, D-053, and D-054. |
| Follow-up | Preserve reduced-cap regressions and run CUDA-capable plus scheduler-backed multi-node qualification before merge readiness. |

### D-058: Detect Omitted Ghost Counts And Bound Remaining Python Temporaries

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Reject ghost-bearing ATHDF-like products when callers omit the matching `num_ghost`; reject duplicate logical MeshBlocks and MeshBlock geometry outside the ghost-extended root domain during direct reads; preflight prolongation index arrays before allocation; generate root centers without an intermediate Python list; conservatively preflight symlog metadata vectorization; and preflight spherical-slice final coordinate-generation temporaries while the assembled surface and last payload remain live. |
| Reason | The previous pass corrected the primary data placement and retained-array paths, but default ghost handling could still truncate valid data silently and several adjacent vectorized temporaries remained outside practical live-memory caps. |
| Evidence | Second final Python tooling re-audit; `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; focused omitted-ghost, direct-inventory, reduced-cap prolongation, symlog, and spherical-coordinate regressions; Python reader module returned `111 passed`; full serial IO matrix returned `167 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Supersedes | Tightens D-057. |
| Follow-up | Preserve reduced-cap regressions and run CUDA-capable plus scheduler-backed multi-node qualification before merge readiness. |

### D-059: Release Sparse Sibling State And Bound Remaining Validation Temporaries

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Bind binary MeshBlock geometry to each logical-location envelope; reject impossible singleton-axis ghost widths; account for NumPy coordinate, prolongation-index, and exact-restriction source temporaries; include sparse duplicate-validation work arrays in PDF and spherical-slice live peaks; account for retained PDF reference state while parsing replacement shard headers; and release incorporated PDF local headers plus PDF/spherical sparse arrays before reading each sibling. Preserve the public one-argument `read_pdf_header()` API by keeping replacement-header accounting in a private helper. Preserve the legacy single-MeshBlock ATHDF helper signature while rejecting unsupported ghost-bearing or sliced emitted extents. |
| Reason | Globally plausible geometry can still belong to a different logical block, inactive-axis malformed widths can bypass ghost validation, and individually bounded arrays do not prove bounded live peaks while NumPy source arrays, duplicate-validation arrays, or prior sibling metadata remain live. Internal accounting state is not part of the public reader API. |
| Evidence | Final Python tooling re-audit; `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; reduced-cap, weak-reference, and malformed single-MeshBlock extent regressions; Python reader module returned `125 passed`; full serial IO matrix returned `181 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Supersedes | Tightens D-057 and D-058. |
| Follow-up | Preserve the weak-reference and reduced-cap regressions; run CUDA-capable plus scheduler-backed multi-node qualification before merge readiness. |

### D-060: Exercise Retained ATHDF APIs And Node Coarsened-Binary Example

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Add positive frozen-fixture comparisons for `read_all_ranks_binary_as_athdf()`, `read_coarsened_binary_as_athdf()`, and `read_all_ranks_coarsened_binary_as_athdf()`; execute the documented node-sharded `cbin --assemble-shards` helper command in the promoted MPI example regression. |
| Reason | Negative malformed-shard coverage and manual probes are not durable substitutes for positive automation of retained public APIs and published helper commands. |
| Evidence | Final test/example re-audit; `tst/test_suite/io/test_python_io_readers_cpu.py`; `tst/test_suite/io/test_node_sharding_mpicpu.py`; Python reader module returned `125 passed`; promoted node example returned `1 passed`; full serial IO matrix returned `181 passed`; full MPI IO matrix returned `59 passed`. |
| Follow-up | Keep the frozen-fixture comparisons and promoted helper subprocess assertions in the required IO matrix. |

### D-061: Require Exact Logical Geometry And Preflight PDF Edge Materialization

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Require every binary MeshBlock geometry record to match its exact logical physical interval within floating-point tolerance. Preflight modern explicit-edge parsing before Python float-list and NumPy-array construction; preflight generated-edge and edge-validation temporaries before NumPy work; apply equivalent cumulative parsing and validation checks to legacy PDF headers; and include explicit arrays for later dimensions plus retained reference-header state in every live-peak budget. |
| Reason | A containment envelope still accepted shifted or shrunken MeshBlock geometry, and bounded final PDF arrays did not prove that parsing, generation, or monotonicity-validation temporaries stayed under the practical live-memory ceiling before materialization. |
| Evidence | Final Python tooling re-audit correction; `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; shifted-geometry and reduced-cap parser regressions; Python reader module returned `130 passed`; full serial IO matrix returned `186 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; detached Pages warnings-as-errors build passed; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Supersedes | Tightens D-059. |
| Follow-up | Preserve reduced-cap pre-materialization regressions and run CUDA-capable plus scheduler-backed multi-node qualification before merge readiness. |

### D-062: Use Bounded Absolute Geometry Tolerance And Strictly Bound Legacy PDF ASCII Expansion

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Compare binary MeshBlock physical intervals with zero relative tolerance and a storage-dtype-aware absolute tolerance capped at one eighth of the logical block width. Count ASCII tokens without materializing a split list, conservatively preflight token strings plus parsed floats before NumPy conversion, reject malformed trailing tokens, accumulate legacy payload-row bytes, and preflight the final two-dimensional legacy stack copy. |
| Reason | Default relative tolerance admitted shifted neighboring blocks on large-offset domains. Header file-size bounds and final histogram bounds did not constrain transient token-list expansion, malformed numeric suffixes, cumulative legacy row arrays, or the final `vstack` coexistence peak. |
| Evidence | Final Python tooling re-audit correction; `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; large-offset geometry, tokenization, strict legacy suffix, cumulative row, and `vstack` regressions; Python reader module returned `139 passed`; full serial IO matrix returned `195 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; detached Pages warnings-as-errors build passed; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Supersedes | Tightens D-061. |
| Follow-up | Preserve strict ASCII and reduced-cap regressions; run CUDA-capable plus scheduler-backed multi-node qualification before merge readiness. |

### D-063: Freeze Generated-Edge Guard Ordering With A Fail-Fast Regression

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Add a reduced-cap regression that replaces PDF edge generation with a failing stub and proves the retained-peak guard rejects the header before NumPy edge generation begins. |
| Reason | The implementation order was correct, but a reduced-cap symlog assertion alone would still pass if a later refactor moved the guard after materialization. |
| Evidence | Final test/example re-audit; `tst/test_suite/io/test_python_io_readers_cpu.py`; Python reader module returned `139 passed`; full serial IO matrix returned `195 passed`; repository style returned `2 passed`; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Follow-up | Keep this assertion beside explicit-edge, legacy-row, and `vstack` pre-materialization regressions. |

### D-064: Forward Retained Legacy-Header Budgets And Bound Spherical Metadata Tokenization

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Forward externally retained bytes through private legacy `.bins.pdf` header delegation. Bound cumulative spherical-slice metadata bytes and every header line before decoding, count variable tokens without splitting, and preflight token-string expansion before materializing the variable list. |
| Reason | Sparse replacement-header parsing can encounter a legacy candidate while reference arrays and the dense accumulator remain live. Spherical whole-file bounds alone do not prevent a malformed metadata line from expanding into decoded strings and a large split-token list before rejection. |
| Evidence | Final Python tooling re-audit correction; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; retained legacy-candidate, oversized spherical metadata-line, and variable-token preflight regressions; Python reader module returned `142 passed`; full serial IO matrix returned `198 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; detached Pages warnings-as-errors build passed; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Supersedes | Tightens D-062. |
| Follow-up | Preserve reduced-cap parser-boundary regressions; run CUDA-capable plus scheduler-backed multi-node qualification before merge readiness. |

### D-065: Enforce ASCII PDF Tokens And Retain Spherical Sibling Metadata Budgets

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Reject non-ASCII PDF numeric rows before applying the ASCII-sized token peak model. Reject duplicate spherical `variables:` declarations; retain a bounded private variable-metadata summary; include reference metadata and dense reconstruction state while parsing and loading sibling shards; release non-reference candidate headers promptly; and carry metadata through duplicate-ownership checks. |
| Reason | Python accepts Unicode numeric tokens that occupy more memory than the ASCII parser model. Spherical variable lists remain live across sibling reads, and previous candidate headers can remain bound unless explicitly released. |
| Evidence | Final Python tooling re-audit correction; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; Unicode numeric, duplicate variable declaration, retained reference metadata, and sibling coexistence regressions; Python reader module returned `145 passed`; full serial IO matrix returned `201 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; detached Pages warnings-as-errors build passed; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Supersedes | Tightens D-064. |
| Follow-up | Preserve parser-model and sibling-coexistence regressions; run CUDA-capable plus scheduler-backed multi-node qualification before merge readiness. |

### D-066: Prove Public Spherical Sibling-Budget Forwarding In Automation

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Extend the public two-shard spherical lifecycle regression to assert that each sibling read after the first receives a nonzero externally retained-byte budget. |
| Reason | A private `_read_header()` reduced-cap test proved parser behavior but would not catch a future regression that stopped forwarding retained state from `read_sphslice()`. |
| Evidence | Final test/example re-audit correction; `tst/test_suite/io/test_python_io_readers_cpu.py`; Python reader module returned `145 passed`; full serial IO matrix returned `201 passed`; repository style returned `2 passed`; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Follow-up | Keep the forwarding assertion in the public multi-shard lifecycle test. |

### D-067: Retain Spherical Metadata Through Coordinate Generation

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Keep the private spherical variable-metadata byte summary until after theta/phi construction and include it in the final coordinate-generation peak check. |
| Reason | The returned header still owns its variable strings while coordinate temporaries are materialized; removing the summary before that check undercounted the live peak. |
| Evidence | Final Python tooling re-audit correction; `vis/python/read_sphslice.py`; focused retained-metadata coordinate regression; Python reader module returned `146 passed`; full serial IO matrix returned `202 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; detached Pages warnings-as-errors build passed; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Supersedes | Tightens D-065. |
| Follow-up | Keep retained metadata in every spherical live-peak calculation until the corresponding strings are no longer reachable. |

### D-068: Reject Mixed Emitted MeshBlock Extents Within One Binary File

| Field | Value |
| --- | --- |
| Status | Accepted and verified locally |
| Decision | Require every nonempty `.bin` or `.cbin` file to use one uniform emitted MeshBlock extent. Reject a later record with a different extent in the shared MeshBlock decoder before exposing parsed file data or entering athdf-like reconstruction. |
| Reason | Dense reconstruction derives one block size from the first MeshBlock. Accepting a later differently shaped block can silently leave uncovered cells initialized to zero instead of reporting malformed output. |
| Evidence | Final Python tooling re-audit correction; `vis/python/bin_convert.py`; direct malformed two-MeshBlock `.bin` and `.cbin` regressions in `tst/test_suite/io/test_python_io_readers_cpu.py`; Python reader module returned `148 passed`; full serial IO matrix returned `204 passed`; full MPI IO matrix returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; detached Pages warnings-as-errors build passed; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Follow-up | Preserve direct-file mixed-extent rejection beside the cross-shard output-shape checks. |

### D-069: Accept The Expanded Robustification Register And Sequence

| Field | Value |
| --- | --- |
| Status | Accepted |
| Decision | Execute `RCP-01` through `RCP-10` in the order defined by `IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md`. Retain `ROB-001` through `ROB-022` and add `ROB-023` through `ROB-027`: ignored restart-header reposition failures, signed output-sequence domain failures, incorrect `cbin` global-`gid` filtering, converter-side allocation preflight, and preservation-aware deferred Pages staging. |
| Reason | Three independent `RCP-00` baseline auditors confirmed the seeded register and identified additional locally actionable defects. The checkpoint order still minimizes risk: restart arithmetic first, shared MPI/publication behavior second, format semantics third, diagnostics fourth, then maintainability, tooling, Pages staging, scaling evidence, external qualification, and packaging. |
| Alternatives | Reorder format or documentation work ahead of restart safety; fold all new findings into broad existing rows without durable identifiers. |
| Why not | Restart correctness and shared publication helpers are prerequisites for later work. Separate identifiers keep new evidence auditable. |
| Reversal path | Append a superseding decision if a checkpoint reflection shows that a prerequisite must move earlier or a finding belongs to a different checkpoint. Re-run affected audits before resuming. |
| Evidence | Frozen guide commit `e40621d81e1948567415a1435f531b275fe14c2e`; runtime, format/tooling, and documentation/process baseline audits; fresh local regression floors recorded in `IO_FEATURE_AUDIT_LEDGER.md`. |
| Follow-up | Track every row in the robustification checkpoint board and do not close a row without independent re-audit. |

### D-070: Use A Small Shared Restart-Layout Descriptor

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-01` implementation |
| Decision | Add a narrowly scoped shared restart-layout descriptor and checked unsigned add/multiply helpers used by both writer and reader. Keep module-specific serialization loops in their existing files. |
| Reason | Writer and reader currently duplicate field-size arithmetic with pre-assignment `int` overflow risk. A shared descriptor makes parity explicit while avoiding an unrelated restart rewrite. |
| Alternatives | Mirror checked arithmetic independently in `restart.cpp` and `pgen.cpp`; add only scattered narrow casts; redesign the restart format. |
| Why not | Mirrored logic can drift, scattered casts are difficult to audit, and a format redesign is outside scope. |
| Reversal path | If the descriptor begins to absorb serialization control flow or module ownership, reduce it to shared arithmetic utilities and retain an explicit parity test. |
| Evidence | `ROB-001`; `src/outputs/restart.cpp`; `src/pgen/pgen.cpp`; runtime baseline audit. |
| Follow-up | Add artificial large-extent and overflow harness coverage before closing `RCP-01`. |

### D-071: Use Checked POSIX Large-File Serial Positioning

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-01` implementation |
| Decision | Implement serial positioned IO with `fseeko()` and `ftello()`. Reject offsets above the representable `off_t` range before conversion, check seek failures, and reject negative positions explicitly. |
| Reason | AthenaK production platforms are POSIX systems. `fseeko()` and `ftello()` provide the intended large-file interface without narrowing through `long` on platforms where that differs from `off_t`. |
| Alternatives | Continue using unchecked `fseek()` and `ftell()`; introduce a broad platform abstraction; rely on MPI IO for all large files. |
| Why not | Existing serial paths are public behavior and need direct protection. A broader abstraction is unnecessary unless non-POSIX support is requested. |
| Reversal path | Add a platform wrapper with equivalent checked semantics if AthenaK gains a supported non-POSIX build target. |
| Evidence | `ROB-004`; `ROB-023`; `src/outputs/io_wrapper.cpp`; `src/parameter_input.cpp`. |
| Follow-up | Add a focused serial wrapper harness and check the parameter-header reposition call site. |

### D-073: Promise Namespace-Atomic Publication, Not Crash Durability

| Field | Value |
| --- | --- |
| Status | Accepted direction for `RCP-02` |
| Decision | Describe temporary-file publication followed by rename as namespace-atomic publication. Do not promise crash-durable persistence. Add checked directory creation and reconcile `.bin`, `.cbin`, modern PDF, `sphslice`, and node-restart publication under this narrower contract. Preserve legacy PDF bytes and behavior. |
| Reason | Rename prevents readers from observing partially written public files. Full crash durability would require payload sync plus directory sync ordering across filesystems and MPI paths, which is a materially larger operational contract. |
| Alternatives | Promise full crash durability now; retain direct public-path writes for `.bin` and `.cbin`; change legacy PDF publication. |
| Why not | The branch needs a truthful, consistent contract without claiming guarantees it does not implement. Legacy PDF compatibility is explicitly preserved. |
| Reversal path | Add a separate durability feature with filesystem-specific sync semantics and failure-injection qualification if production requirements demand it. |
| Evidence | `ROB-006`; baseline publication audit; existing temporary-file behavior in modern PDF, `sphslice`, and node restart. |
| Follow-up | Resolve helper boundaries and stale-temporary behavior during `RCP-02`. |

### D-075: Keep Coarsened Binary Deliberately Uniform-3D In This Branch

| Field | Value |
| --- | --- |
| Status | Accepted direction for `RCP-03` |
| Decision | Support uniform three-dimensional full-volume `.cbin` in shared, rank, and node layouts. Permit ghost-expanded full-volume extents only when divisible by the validated factor. Reject lower-dimensional meshes, AMR meshes, and sliced output during construction until each receives a deliberately designed format contract. |
| Reason | The current lower-dimensional reader arithmetic can collapse singleton axes to zero, AMR logical-location semantics remain unresolved, and sliced output is already deliberately excluded. A narrow explicit contract is safer than partially promoting ambiguous behavior. |
| Alternatives | Promote lower-dimensional and AMR output immediately; silently rely on current incidental behavior. |
| Why not | Promotion requires new reader/writer semantics, positive reconstruction evidence, and documentation beyond the risk-reduction goal of this branch. |
| Reversal path | Add support row-by-row in a later feature branch with writer, reader, conversion, and shard-equality tests. |
| Evidence | `ROB-005`; format/tooling baseline audit; `src/outputs/outputs.cpp`; `src/outputs/coarsened_binary.cpp`; `vis/python/bin_convert.py`. |
| Follow-up | Enforce the matrix, repair global-`gid` filtering, and check Kokkos ranges during `RCP-03`. |

### D-081: Add Keyword-Only Python Reader Limit Overrides

| Field | Value |
| --- | --- |
| Status | Accepted direction for `RCP-06` |
| Decision | Preserve safe defaults and expose deliberate reader/converter budget overrides through keyword-only configuration objects or keyword-only arguments. Do not use environment variables as the primary public API. Add CLI flags only where an existing CLI workflow needs equivalent control. |
| Reason | Python analysis users sometimes need larger trusted products, but hidden constants force source edits. Keyword-only overrides preserve compatibility and make elevated budgets explicit at the call site. |
| Alternatives | Keep fixed private constants; use environment variables only; add positional parameters. |
| Why not | Fixed constants are inflexible, environment-only behavior is hard to audit, and positional additions risk API ambiguity. |
| Reversal path | Consolidate keyword-only arguments into a shared immutable configuration object if local repetition becomes harder to maintain than the object boundary. |
| Evidence | `ROB-012`; `ROB-026`; format/tooling baseline audit. |
| Follow-up | Inventory each public reader and converter allocation before implementing `RCP-06`. |

### D-083: Treat Replicated Manifest Validation As Evidence-Driven Optimization

| Field | Value |
| --- | --- |
| Status | Accepted with measurement gate |
| Decision | Retain strict replicated manifest and payload-header validation on every rank until production-topology measurements justify a structured-broadcast or node-leader optimization. Treat measurement and an explicit keep-or-refactor decision as required evidence, but do not refactor preemptively. |
| Reason | Replicated validation is simple and defensible for correctness. Optimization without topology evidence risks weakening validation or adding collective asymmetry. |
| Alternatives | Centralize immediately on rank 0; validate only on node leaders; ignore scaling evidence. |
| Why not | The correct boundary depends on filesystem and topology measurements unavailable on the local one-node host. |
| Reversal path | Implement the preregistered optimization if `RCP-08` and `RCP-09` measurements exceed the threshold recorded in `D-090`. |
| Evidence | `ROB-010`; existing D-046; local topology probe reports one physical host only. |
| Follow-up | Add instrumentation and preregister thresholds before scheduler-backed qualification. |

### D-085: Bound Node-Restart Manifest Materialization

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-01` implementation |
| Decision | Reject node-restart manifests above 64 MiB before parsing. Bound signature lines at 256 bytes, scalar records at 256 bytes, payload records at 4096 bytes, segment records at 256 bytes, trailing records at 256 bytes, and generated payload paths at 1024 bytes. Apply the signature-line bound in both format probing and full parsing. |
| Reason | The existing record-count caps do not prevent a malformed file from forcing an unbounded `std::string` allocation before semantic validation. These limits exceed plausible production records while stopping pathological inputs early. |
| Alternatives | Depend only on payload and segment count caps; use substantially smaller limits; permit unbounded line reads. |
| Why not | Count caps do not bound bytes, very small limits create avoidable topology-path constraints, and unconstrained reads violate the strict parser contract. |
| Reversal path | Raise an individual limit with production evidence and regression coverage; keep the total cap and early-rejection structure. |
| Evidence | `ROB-019`; `src/restart_manifest.cpp`; runtime baseline audit. |
| Follow-up | Add oversized signature, scalar, payload, segment, trailing, path, and total-file regressions. |

### D-091: Include Restart Metadata Reconstruction In The Shared Layout Boundary

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-01` implementation |
| Decision | Keep the shared restart-layout facility narrow, but use its checked arithmetic for replicated metadata reconstruction in `Mesh::BuildTreeFromRestart()` as well as per-MeshBlock payload layout in the writer and reader. Preserve serialization order, wire bytes, and module-specific read/write loops. |
| Reason | The focused pre-edit layout audit found that `src/mesh/build_tree.cpp` reconstructs `listsize * nmb_total` into an allocation, broadcast byte count, and signed cursor. Hardening only `restart.cpp` and `pgen.cpp` would leave the header path vulnerable before payload routing begins. |
| Alternatives | Restrict the descriptor to payload arithmetic; redesign restart serialization around a broad object model; patch metadata reconstruction with scattered casts only. |
| Why not | The narrow descriptor can cover the arithmetic contract without absorbing control flow or changing persisted bytes. Scattered casts remain difficult to audit, while a broader restart rewrite increases compatibility risk. |
| Reversal path | Reduce the facility to checked arithmetic helpers if it starts to own serialization behavior. Retain explicit writer-reader parity tests and frozen resume coverage. |
| Evidence | Focused `RCP-01` pre-edit restart-layout audit; `src/outputs/restart.cpp`; `src/pgen/pgen.cpp`; `src/mesh/build_tree.cpp`. |
| Follow-up | Add artificial large-count harness coverage and preserve frozen shared and per-rank restart resumes. |

### D-092: Enforce Symmetric Restart Manifest Budgets And Early Header Rejection

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-01` closure |
| Decision | Make the strict node-manifest writer obey the same 64 MiB total-byte and 1024-byte generated payload-path limits enforced by the reader. Validate serialized restart mesh structure before root-grid division, metadata allocation, or payload routing. Expand serial positioned-IO checks through the terminal byte and through `size_t` representability. |
| Reason | A writer that can publish a manifest its own reader rejects violates the transaction contract. Likewise, checked payload arithmetic does not protect startup if malformed serialized mesh dimensions are consumed first. Serial positioned IO must prove that the entire requested span, not only its first byte, is representable. |
| Alternatives | Treat reader limits as malformed-input defenses only; validate restart structure after allocation; rely on standard-library narrowing or seek failure for terminal-range rejection. |
| Why not | Those alternatives leave locally generated unreadable products, permit unsafe arithmetic before validation, and make overflow behavior platform-dependent. |
| Reversal path | Raise reader and writer manifest budgets together with production-topology evidence and focused regressions. Broaden structural validation only if a valid historical restart fixture demonstrates a missing compatibility case. |
| Evidence | `ROB-029` through `ROB-032`; independent arithmetic audit; `src/restart_layout.hpp`; `src/mesh/build_tree.cpp`; `src/outputs/io_wrapper.cpp`; `src/outputs/restart.cpp`; focused restart/layout and serial-wrapper tests returned `21 passed`; MPI restart tests returned `56 passed`; detached pre-edit shared and per-rank restart bytes matched exactly. |
| Follow-up | Close `RCP-01` only after a correction-focused independent re-audit and the full local matrix pass. |

### D-093: Reject Malformed Restart Topologies Before Tree Mutation

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-01` closure |
| Decision | Treat restart leaf metadata as an untrusted persisted topology. Before constructing `MeshBlockTree`, require signed-safe logical levels `<=30`, representable axis capacities, valid active and inactive axes, positive finite per-block costs, finite aggregate cost, unique locations, no ancestor overlap, complete dimensionally active children, and canonical persisted Z order. After safe construction, traverse neighbors to enforce the existing 2:1 balance contract and validate that load balancing assigns one contiguous nonempty range to every rank. Apply the same signed-safe maximum-level and wide adaptive-level arithmetic to scratch construction. |
| Reason | Tree insertion, geometry creation, neighbor discovery, and load balancing assume a complete canonical leaf inventory. A malformed restart could otherwise reach signed shifts, null-child dereferences, silent payload relabeling, or partially initialized rank maps before a useful error. |
| Alternatives | Depend on post-construction MeshBlock counts; harden only individual tree methods; accept permuted valid leaf sets and reorder payload metadata implicitly. |
| Why not | Counts do not prove topology completeness, tree-local checks occur after unsafe input has entered recursion, and implicit reordering disconnects serialized payload records from their logical locations. |
| Reversal path | If a historical valid fixture demonstrates a missing topology case, refine the preflight with a fixture-backed compatibility rule. Do not weaken signed-safe levels or canonical payload ordering. |
| Evidence | `ROB-033` through `ROB-038`; correction-focused topology audits; `src/mesh/build_tree.cpp`; `tst/inputs/io_restart_metadata.athinput`; corruption and positive regressions in `tst/test_suite/io/test_io_finalization_timing_cpu.py`; focused serial restart/layout set returned `41 passed`; focused MPI restart set returned `58 passed`; full serial matrix returned `238 passed`; full MPI matrix returned `67 passed`; style returned `2 passed`; fixture checksum verification passed for all `27` artifacts. |
| Follow-up | Retain the restart-specific preflight locally. Revisit a more general tree validation API only if later non-restart paths need the same boundary. |

### D-072: Use A Narrow Shared MPI Failure Utility

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-02` implementation |
| Decision | Add a header-only shared MPI utility that renders MPI failures with `MPI_Error_string()` on a best-effort basis, falls back to the numeric code if rendering itself fails, checks actionable branch-added and inherited-but-touched MPI returns, and aborts `MPI_COMM_WORLD` for MPI-sensitive fatal paths. Treat `MPI_Comm_free(node_comm)` failure as fatal during normal teardown and clear the handle only after successful release. Keep intentional `MPI_Abort()` shutdown calls simple. |
| Reason | The diff-driven inventory found fragmented checks and unchecked returns across communicator lifecycle, metadata collectives, output reductions, timing, and restart publication. World abort prevents a node-local or rank-local failure from stranding peers in a later collective. |
| Alternatives | Keep file-local helpers; check only newly added calls; return errors from deep output paths; broaden the checkpoint into all inherited MPI startup and shutdown code. |
| Why not | File-local helpers duplicate rendering and miss collective-participation hazards. Deep error propagation is a larger interface redesign. Inherited startup and finalization cleanup is outside the narrow output checkpoint unless a touched path depends on it. |
| Reversal path | Promote the header-only helper into a compiled utility if later users require non-inline state or policy. Keep the rendering fallback and simple world-abort semantics. |
| Evidence | `ROB-039`; `ROB-040`; pre-edit MPI auditor classified `83` direct call sites across `13` touched source files. |
| Follow-up | Re-audit every diff-touched MPI call after implementation and preserve wrapper partial-IO propagation where callers deliberately inspect counts. |

### D-074: Use Checked Directories, Owned Temporaries, And Exclusive Node Generations

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-02` implementation |
| Decision | Add a narrow output-filesystem helper for checked directory creation, deterministic temporary names, owned-temporary cleanup, and rename publication. Accept an existing directory, reject a conflicting non-directory path, and report permission errors. Adopt private temporary plus rename publication for `.bin`, `.cbin`, modern PDF, `sphslice`, and every restart layout while preserving legacy PDF append bytes and behavior. For node restarts, reserve a generation token exclusively before writing and clean only attempt-owned reservation or temporary artifacts. Do not sweep ambient directories. |
| Reason | Current publication guarantees differ by format, raw `mkdir()` failures are ignored, and clock-derived node generations are not exclusive. A small shared helper removes repeated error-prone filesystem code without hiding communicator ordering. |
| Alternatives | Promise crash-durable publication; retain direct-public `.bin`, `.cbin`, and shared restart writes; sweep stale files broadly; change legacy PDF publication. |
| Why not | Crash durability requires payload and directory sync ordering that is not implemented or qualified. Direct writes expose incomplete public files. Broad sweeps can remove unrelated work. Legacy PDF compatibility is frozen byte-for-byte. |
| Reversal path | Add a separate durability feature with sync semantics and filesystem qualification if required. Raise or alter generation policy only with fault-injection coverage. |
| Evidence | `ROB-006`; `ROB-041`; `ROB-042`; pre-edit filesystem-publication audit. |
| Follow-up | Add stale-temporary, directory-conflict, rename-failure, generation-collision, and payload-before-manifest regressions. |

### D-086: Render Minimum-Five-Digit Sequences Without Truncation

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-02` implementation |
| Decision | Replace fixed `%05d` sequence buffers with one dynamic renderer. Preserve a minimum width of five ASCII digits, emit wider values without truncation, reject negative counters before publication, and use a checked advance helper. Keep persisted counters as `int` for compatibility; the maximum emittable value is `INT_MAX - 1` because every successful output advances and persists the next counter. Change Python PDF inference to accept `[0-9]{5,}` while preserving historical five-digit names. |
| Reason | Existing fixed buffers deterministically reuse the `10000` namespace at `100000` and `100001`. Widening alone is incomplete because negative values and signed overflow remain unsafe. |
| Alternatives | Reject every value above `99999`; migrate persisted counters to an unsigned or wider type; widen buffers independently at each writer. |
| Why not | A five-digit maximum is an unnecessary operational limit, a persisted type migration is broader than needed, and piecemeal buffers can drift. |
| Reversal path | Migrate persisted counters to a wider integer in a separate compatibility-reviewed change if runs need more than `INT_MAX - 1` numbered publications. |
| Evidence | `ROB-017`; `ROB-024`; independent sequence inventory across table, VTK, particle VTK, Cartesian-grid, spherical-surface, `.bin`, `.cbin`, PDF, `sphslice`, and restart writers. |
| Follow-up | Add exact `99999`, `100000`, and `100001` naming tests plus negative and exhausted-counter rejection before publication. |

### D-087: Reserve Deterministic Output Families Before Writer Construction

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-02` implementation |
| Decision | Preflight active output blocks before constructing writers. Register typed public-target templates using `{SEQ}` and `{PARTITION}` placeholders so blocks collide even when initial counters or cadences differ. Cover `.bin`, `.cbin`, modern and legacy PDF data or metadata overlap, `sphslice`, restart manifests, and inherited deterministic numbered or append families where their keys are straightforward. Permit genuinely disjoint shared, rank, and node roots and distinct `.cbin` factors. Reject unsafe generated `id` components and PDF directory components. Preserve historical basename handling except for the existing stricter node-restart safe-leaf rule. |
| Reason | Constructors currently create directories while parsing blocks and no namespace map exists. Modern PDF directories omit the first axis variable, allowing same-`id` blocks with different histogram definitions to overwrite one another. Generated `id` and PDF directory fragments can also escape intended paths. |
| Alternatives | Detect only exact current filenames; reject all non-leaf basenames globally; refactor every writer around a new path object immediately. |
| Why not | Current-filename checks miss future overlap, global basename tightening risks historical configurations without evidence, and a broad path-object rewrite exceeds this checkpoint. |
| Reversal path | Expand safe-component rules or promote shared writer path renderers after compatibility evidence. Keep family reservations before any writer constructor side effects. |
| Evidence | `ROB-021`; `ROB-043`; independent filesystem and namespace audits. |
| Follow-up | Add duplicate-family rejection and allowed-control regressions, including same explicit PDF `id` with different variables, bins, scales, or weights. |

### D-094: Keep The Coarsened-Binary Producer Active-Zone-Only

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-03`; supersedes `D-075` and resolves `D-088` |
| Decision | Produce `.cbin` only for uniform three-dimensional full-volume active-zone output in shared, rank, and node layouts. Reject `ghost_zones=true`, lower-dimensional meshes, AMR meshes, and slices during construction. Represent checked Kokkos plane strides, factor powers, coarsening ranges, and normalization ranges as `std::int64_t`; preflight retained allocations and serialized payload sizes separately as `std::size_t`. |
| Reason | Ghost-expanded output has no promoted producer contract and introduces avoidable ambiguity. The narrow active-zone producer is explicit, documented, and sufficient for the feature. Kernel ranges and host allocations have distinct representability domains and should not be conflated. |
| Alternatives | Permit divisible ghost-expanded extents; promote lower-dimensional, AMR, or sliced output now; use `int` launch ranges; use `std::size_t` directly inside Kokkos launch arithmetic. |
| Why not | The broader rows need deliberate reconstruction semantics and evidence. Signed checked 64-bit launch arithmetic matches the Kokkos kernel boundary while `std::size_t` remains the correct allocation boundary. |
| Reversal path | Promote one producer row at a time in a later branch with writer, reader, ATHDF/XDMF, sharding-equality, and documentation evidence. |
| Evidence | `src/outputs/coarsened_binary.cpp`; `src/outputs/coarsened_binary_layout.hpp`; `tst/test_suite/io/cbin_layout_harness.cpp`; focused arithmetic harness returned `8 passed`; MPI moment conversion regression passed for shared, rank, and node layouts. |
| Follow-up | Retain explicit constructor rejection tests and the reader-compatibility choice in `D-095`. |

### D-095: Preserve Broader Read-Only Legacy Coarsened-Binary Compatibility

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-03` |
| Decision | Keep `vis/python/bin_convert.py` able to read historical `.cbin` products that are broader than the active producer contract when their emitted metadata and payload are internally valid. Treat this as read-only legacy compatibility, not as permission for the C++ writer to emit those forms. Continue rejecting malformed moment groups before payload decoding. |
| Reason | Tightening the producer prevents new ambiguous files. Tightening the reader would unnecessarily strand historical analysis products and is not required for writer safety. |
| Alternatives | Reject every historical `.cbin` row excluded from the current producer; silently broaden the producer to match the reader. |
| Why not | Reader compatibility is useful and low risk when validation is strict. Producer promotion requires a stronger design and test matrix. |
| Reversal path | Narrow an individual legacy reader row only if a concrete malformed-input risk cannot be handled by validation, and document the compatibility break. |
| Evidence | `vis/python/bin_convert.py`; early malformed-moment regressions in `tst/test_suite/io/test_python_io_readers_cpu.py`; independent RCP-03 re-audit. |
| Follow-up | Document the producer-versus-reader distinction in the compatibility matrix and deferred Pages overlay. |

### D-076: Define Finite Diagnostic Semantics At Geometric Singularities

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-04` |
| Decision | Clamp mathematically bounded `z/r` inputs to `[-1, 1]` before inverse-cosine and cosine-derived outputs. At `r=0`, define radial velocity and radial mass or energy projections as zero. At cylindrical `R=0`, define cylindrical radial, azimuthal, and polar velocity projections as zero. Reject non-positive density and non-finite fluid diagnostics at runtime rather than emitting sentinels. Retain the existing construction-time rejection for ambiguous generic ion-neutral two-fluid diagnostics. |
| Reason | Roundoff may move a theoretically bounded projection slightly outside its domain. Coordinate-axis singularities have a useful finite limiting convention for output, while division by non-positive density does not have a defensible generic fluid interpretation. |
| Alternatives | Emit NaNs; apply a fluid floor silently; reject every axis-origin sample; emit a sentinel for non-positive density. |
| Why not | NaNs make analysis failures late and opaque. Silent floors alter physical meaning. Axis-origin rejection is unnecessarily strict for otherwise valid meshes. Sentinels are not self-describing in ordinary output fields. |
| Reversal path | Add module-qualified diagnostics if a fluid package establishes a different physical policy. Keep generic names conservative. |
| Evidence | `src/outputs/diagnostic_semantics.hpp`; `src/outputs/derived_variables.cpp`; `tst/test_suite/io/diagnostic_semantics_harness.cpp`; analytic Hydro and MHD writer regressions. |
| Follow-up | Preserve representative CUDA execution as an external `RCP-09` gate. |

### D-077: Reject Derived Ghost-Zone Output Until Kernels Populate It

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-04` |
| Decision | Reject `ghost_zones=true` during construction whenever an output requests derived diagnostics. Modern PDFs additionally reject ghost-zone sampling explicitly and operate on active zones only. Do not attempt a partial derived-kernel ghost-zone implementation in this branch. |
| Reason | Existing derived arrays populate active zones. Allowing ghost-zone output would publish uninitialized or stale values. A future implementation must define every derived kernel's ghost-domain behavior consistently. |
| Alternatives | Populate selected derived kernels only; allow the request and rely on callers to avoid ghosts; broadly redesign derived kernels now. |
| Why not | Partial support is difficult to explain and easy to misuse. Silent allowance emits incorrect data. A broad redesign exceeds this IO branch's bounded scope. |
| Reversal path | Promote derived ghost-zone support in a focused branch with a complete kernel inventory and boundary-value regressions. |
| Evidence | `src/outputs/basetype_output.cpp`; `src/outputs/pdf.cpp`; derived ghost-zone construction regressions. |
| Follow-up | Keep `sphslice` derived-field rejection until ghost-zone-safe interpolation exists. |

### D-078: Use Round-Trip Scientific Spherical-Slice Radius Tokens

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-04` |
| Decision | Render `sphslice` radius filename components with deterministic round-trip scientific precision, for example `r_2.5000000000000000e-01`. Permit separately configured spherical slices to share an `id` when their precise radius tokens differ. |
| Reason | The inherited `%g` representation can collapse nearby radii into one public path. Round-trip scientific tokens remain human-readable while preserving the configured floating-point distinction. |
| Alternatives | Keep `%g`; append a sequence-independent slice ordinal; reject duplicate IDs across all spherical slices. |
| Why not | `%g` can collide. Ordinals hide the physical radius and make configuration ordering observable. Blanket duplicate-ID rejection blocks safe, useful multi-radius families. |
| Reversal path | Add a versioned naming layer only if downstream tooling requires shorter tokens; retain collision resistance. |
| Evidence | `src/outputs/output_file_utils.hpp`; `src/outputs/spherical_slice.cpp`; namespace preflight; nearby-radius regression. |
| Follow-up | Preserve the explicit distinction between `file_type=sph` and `file_type=sphslice`. |

### D-089: Bound PDF And Spherical-Slice Writer Allocations

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-04` |
| Decision | Add a positive per-output `max_writer_allocation_bytes` parameter with a default of `536870912` bytes. Preflight modern PDF histogram, edge, host-mirror, copied-field, derived-field, metadata, and serialization staging bytes before allocation. Preflight spherical-slice angular geometry, interpolation arrays, sparse buffers, dense staging, and serialization buffers before allocation. |
| Reason | File-controlled dimensions can otherwise trigger unreasonable host or device allocations before a clear diagnostic. A configurable practical cap is safer than an unbounded writer and still permits deliberate larger production products. |
| Alternatives | Hard-code a cap; cap only final payload bytes; depend on allocator failures. |
| Why not | Fixed caps require source edits for legitimate large products. Final payload size omits retained and transient arrays. Allocator failures are late and backend-dependent. |
| Reversal path | Adjust the default with production evidence while keeping positive explicit overrides and early preflight. |
| Evidence | `src/outputs/pdf.cpp`; `src/outputs/spherical_slice.cpp`; reduced-cap constructor regressions and PDF load-cap-before-derived-allocation regression. |
| Follow-up | Keep reader-side configurable budgets separate in `RCP-06`. |

### D-096: Use Backend-Portable PDF Atomics And Host-Staged Reductions

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-04` correction cycle |
| Decision | Accumulate PDF bins with `Kokkos::atomic_add()` into the result view. For MPI reductions, copy result views to host mirrors, reduce host pointers, and copy reduced root results back to device. Include host mirrors in writer budget preflight. |
| Reason | `ScatterView` replication is backend-sensitive and complicates allocation accounting. Passing device pointers directly to MPI silently assumes GPU-aware MPI. The explicit atomic and host-staging path is portable across CPU and accelerator builds. |
| Alternatives | Retain `ScatterView`; require GPU-aware MPI; select behavior through backend-specific compile-time branches. |
| Why not | The branch should not depend on hidden backend allocation multipliers or deployment-specific MPI capabilities. Compile-time divergence increases the qualification matrix without a demonstrated need. |
| Reversal path | Add an opt-in GPU-aware reduction path only after measured benefit and qualification on supported production stacks. |
| Evidence | `src/outputs/pdf.cpp`; local serial and MPI PDF regressions; deferred documentation overlay. |
| Follow-up | Execute the representative CUDA PDF regression during `RCP-09`. |

### D-097: Revalidate Derived Shapes After AMR And Preflight Spherical-Slice Peaks

| Field | Value |
| --- | --- |
| Status | Accepted for reopened `RCP-04` correction cycle |
| Decision | Reallocate the shared derived-variable view whenever any meshblock, variable, or cell extent differs from the current request, rather than only on first use. For `sphslice`, retain a geometry baseline and preflight each transient peak before allocation: ownership validation, shared dense storage, local sparse storage and sorting, node metadata and gather buffers, publication sorting, structurally sized parameter-dump and counted metadata serialization, and dense payload serialization. Compute parameter-dump bytes from the parsed hierarchy without invoking its allocating serializer, stream the actual serialized header directly to the temporary file, and admit sort peaks before mutating retained shard vectors. |
| Reason | Adaptive refinement can increase the local meshblock pack after an earlier diagnostic output, making a first-allocation-only view too small. A constructor-only spherical-slice estimate also misses late node and serialization buffers, especially input-controlled parameter dumps. |
| Alternatives | Allocate derived views at a fixed mesh maximum; retain constructor-only spherical-slice checks; check string sizes only after materialization; remove spherical-slice configurability. |
| Why not | Fixed maxima consume avoidable accelerator memory. Constructor-only checks do not cover later peaks. Post-materialization checks are too late for input-controlled text. Removing configurability would block legitimate production tuning without improving the contract. |
| Reversal path | Replace the conservative peak model only with an audited equivalent that remains allocation-before-use and retains the adaptive-PDF regression. |
| Evidence | `src/outputs/derived_variables.cpp`; `src/outputs/spherical_slice.cpp`; adaptive pack-growth derived-PDF regression in `tst/test_suite/io/test_output_formats_gpu.py`; serialization-cap regression in `tst/test_suite/io/test_output_formats_cpu.py`; fresh correction audit after direct-streaming repair. |
| Follow-up | Require a fresh independent `RCP-04` acceptance audit and retain CUDA execution in `RCP-09`. |

### D-079: Keep C++ Output Registration Cleanup Bounded

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-05` |
| Decision | Keep the existing format-construction chain. Extract common variable and ID parsing, MeshBlock selection, slices, and PDF parameter parsing into focused private helpers. Add in-class defaults for touched `OutputParameters` fields. Convert the branch-added spherical-slice owner to `std::unique_ptr`. Document wrapper byte units, MPI collectivity, lazy node-communicator lifetime, and restart-manifest post-load invariants. Repair the inherited tracked-particle path only where the existing fixed-width format requires correctness: skip irrelevant `variable` parsing, use a device counter, require a dense unique tag range, compute checked byte offsets, preserve native-endian payload bytes, and retain collective ordering. |
| Reason | These changes remove duplicated parsing and make branch-touched contracts explicit without replacing the established output factory style or broadening the wire-format surface. The tracked-particle inconsistency was directly exposed by common-parser cleanup and could corrupt or misplace records if left untouched. |
| Alternatives | Introduce a registration factory or table-driven constructor; redesign all output ownership to RAII; change tracked-particle wire bytes; defer the inherited tracked-particle defect. |
| Why not | A broad factory or ownership rewrite increases review surface without advancing the IO formats. Changing tracked bytes breaks compatibility. Deferral would leave a touched constructor path knowingly incorrect and the no-variable contract inconsistent. |
| Reversal path | Consider a table-driven registration layer in a separate refactor branch after behavior snapshots exist for every inherited format. Version tracked-particle payloads explicitly if a future format redesign is needed. |
| Evidence | `src/outputs/outputs.cpp`; `src/outputs/outputs.hpp`; `src/outputs/track_prtcl.cpp`; touched interface headers; dedicated tracked-particle serial and MPI regressions. |
| Follow-up | Require behavior-preservation and scope audits before closing `RCP-05`. |

### D-090: Preregister Restart-Manifest Scaling Before Optimizing

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-08A`; `RCP-08B` remains pending external evidence |
| Decision | Keep strict per-rank manifest and replicated-header validation in this feature branch until scheduler-backed measurements exist. Add default-off `ATHENAK_RESTART_MANIFEST_TIMING=1` records for validation, startup parse, and direct local-block loading. Preregister 2-node and 4-node topology tiers, one warm-up plus five measured resumes, median/p95/maximum reporting, and explicit keep/follow-up/separate-branch thresholds. |
| Reason | The current protocol is correct and reviewable. A rank-0 or node-leader metadata distribution protocol may reduce filesystem amplification, but changing restart semantics without production topology measurements would increase risk based on speculation. |
| Alternatives | Parse and validate only on rank 0 now; validate headers on node leaders now; omit instrumentation and rely on external wall-clock timing. |
| Why not | Both centralized designs alter the restart protocol and require multi-node requalification. External wall-clock timing alone cannot isolate manifest parse/validation from routed local-span loading. |
| Reversal path | If `RCP-09` evidence exceeds the preregistered threshold, open a separate scaling branch to compare rank-0 structured broadcast and node-leader header validation, then rerun the full restart qualification lane. |
| Evidence | `IO_RESTART_MANIFEST_SCALING_QUALIFICATION.md`; `IO_EXTERNAL_IO_QUALIFICATION_PLAN.md`; `src/restart_manifest.cpp`; `src/main.cpp`; opt-in timing regression in `tst/test_suite/io/test_node_sharding_mpicpu.py`. |
| Follow-up | Benchmark-design audit approved the local instrumentation contract. Run scheduler-backed `RCP-09`, then record the `RCP-08B` keep, follow-up, or separate-branch disposition. |

### D-084: Preserve Active Qualification Records And Package Deliberately

| Field | Value |
| --- | --- |
| Status | Accepted for local `RCP-10` packaging preparation |
| Decision | Retain `IO_FORMAT_COMPATIBILITY.md`, the external qualification plans, the deferred Pages bundle, and its deterministic staging helper with the branch. Keep the large guide, integration, decision, and ledger records at their current review paths until external gates close. Before final merge, summarize durable content into the compatibility contract, user docs, and pull-request record; archive process records in-repository only if maintainers explicitly want that history. |
| Reason | The branch still has active external qualification gates, so removing or moving the evidence now would lose audit context and invalidate the frozen guide path. Long implementation diaries should not silently become permanent user-facing documentation either. |
| Alternatives | Delete process records now; retain every root-level guide indefinitely; move the frozen guide into an archive before qualification finishes. |
| Why not | Immediate deletion loses provenance. Permanent root-level diaries add maintenance bulk. Moving the frozen guide before closure breaks the checksum-path contract used throughout robustification. |
| Reversal path | A maintainer may request an in-repository development archive after qualification. Update the disposition record and preserve a migration map in the cleanup commit. |
| Evidence | `IO_FEATURE_BRANCH_PROCESS_ARTIFACT_DISPOSITION.md`; current open `RCP-08B` and `RCP-09` gates; frozen-guide checksum rule. |
| Follow-up | Revisit after external qualification and before the final merge commit or pull-request handoff. |

### D-082: Stage Deferred Pages Documentation With Bounded Detached-Worktree Edits

| Field | Value |
| --- | --- |
| Status | Accepted for local `RCP-07` |
| Decision | Use `scripts/stage_gh_pages_io_docs.py`, implemented with the Python 3 standard library, and `deferred_docs/gh-pages/io-output-formats-and-sharding/manifest.json` as the machine-readable contract. Require a canonical detached Pages worktree, frozen target and protected blobs, source SHA256 values, expected add-target absence, unique anchors, baseline section SHA256 values, marker-wrapped bounded edits, an exact eight-file allowlist, in-memory planning, pure second-transform idempotence, and atomic writes. Permit drift handling only through `--reviewed-drift <packet-dir>`, which fingerprints the target and emits an external three-way packet without target writes. |
| Reason | Literal complete-page overlays can delete unrelated live Pages guidance while still producing a successful Sphinx build. Bounded verified edits preserve unrelated content and make drift explicit before publication. |
| Alternatives | Continue manual copy-and-merge staging; retain destructive whole-page overlays; edit a checked-out `gh-pages` branch directly. |
| Why not | Manual staging is difficult to reproduce and audit. Broad replacement is unsafe under live Pages drift. Direct branch edits violate the deferred separate-review lifecycle. |
| Reversal path | If the Pages structure changes, run reviewed-drift mode against a fresh detached worktree, audit the packet, update anchors and hashes deliberately, and rerun strict staging. |
| Evidence | `scripts/stage_gh_pages_io_docs.py`; `manifest.json`; bounded fragments; `tst/test_suite/io/test_stage_gh_pages_io_docs_cpu.py`; local `22 passed`; targeted `flake8`; strict detached stage and `--verify-staged`; warnings-as-errors HTML and link-check builds; rendered inspection of all eight staged pages; accepted docs-to-code audit; unchanged local `origin/gh-pages=4833aa9341e19861297e330ff02aabfd8001935c`. |
| Follow-up | After the code branch merges, refresh `origin/gh-pages`, restage in a fresh detached worktree, and open a separate Pages review. |

### D-080: Extract Shared Private Python Reader Utilities

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-06` |
| Decision | Extract common private reader-limit, checked-arithmetic, ASCII-token, and CLI-limit helpers into `vis/python/io_reader_common.py`. Keep `vis/python/bin_convert.py` as the one supported binary converter, add package-relative imports with direct-script fallbacks, and thread retained metadata accounting through binary, coarsened-binary, PDF, and spherical-slice readback. |
| Reason | The readers need one consistent practical-limit contract while preserving both existing script-style imports and package-form imports. Local duplicate helpers made it too easy for one format to bypass a live-memory preflight. |
| Alternatives | Retain duplicated format-local helpers; expose a new public configuration module; restore `bin_convert_new.py` as a compatibility alias. |
| Why not | Duplicated helpers already diverged. A new public module would unnecessarily expand the supported API. Restoring the redundant converter would recreate the ambiguity this branch is intended to remove. |
| Reversal path | Add public configuration only if user workflows demonstrate a need beyond keyword-only `limits=` overrides and shared CLI flags. |
| Evidence | `vis/python/io_reader_common.py`; `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; `vis/python/examples/read_io_outputs.py`; `tst/test_suite/io/test_python_io_readers_cpu.py`; independent API/memory-budget correction audit; `181 passed`. |
| Follow-up | Retain malformed-input, package-import, direct-script, and reduced-budget rows in the final RCP-10 matrix. |

### D-098: Document Native Multi-Field Spherical Slices

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-07` documentation alignment |
| Decision | Preserve and document `sphslice` support for native state-backed scalar fields and native multi-field groups. Continue rejecting derived arrays until interpolation is ghost-zone-safe. Document that the origin-centered spherical surface must fit inside a 3D domain, with a positive radius strictly interior to every domain face, and add producer-reader coverage for a native hydro group. |
| Reason | The writer serializes `outvars.size()` variables and the public reader returns `(theta, phi, variable)` arrays. Native groups are therefore an implemented format capability, not an accidental parser side effect. |
| Alternatives | Reject every group that expands to more than one native variable; advertise groups without adding executable evidence. |
| Why not | Artificially narrowing the writer would remove a coherent existing behavior. Advertising an untested behavior would leave the docs-to-code audit unresolved. |
| Reversal path | Narrow the contract later only through an explicit compatibility decision and migration note if production use shows that group output is unsafe. |
| Evidence | `src/outputs/spherical_slice.cpp`; `vis/python/read_sphslice.py`; native-group and invalid-domain regressions in `tst/test_suite/io/test_output_formats_cpu.py`; deferred Pages fragments; `IO_FORMAT_COMPATIBILITY.md`. |
| Follow-up | Rerun detached Pages staging and obtain a fresh docs-to-code audit after helper-safety corrections land. |

### D-099: Prove Deferred Pages Staging Against A Canonical Index

| Field | Value |
| --- | --- |
| Status | Accepted for `RCP-07` |
| Decision | Before trusting status output, require the detached Pages worktree index to match `HEAD` exactly and reject noncanonical index flags such as `assume-unchanged` and `skip-worktree`. Require every staged target to remain a regular file with its expected mode. Include target-root metadata in reviewed-drift fingerprints. |
| Reason | Porcelain status alone can hide tracked mutations behind index flags. Byte comparison alone can accept an allowlisted symlink to identical bytes or an executable-bit change. Descendant-only fingerprints cannot prove that reviewed-drift packet generation leaves the worktree root untouched. |
| Alternatives | Keep status-only verification; document that operators must avoid index flags; ignore Markdown file-mode changes; fingerprint descendants only. |
| Why not | The helper is a publication-safety boundary. Its checks must fail closed under Git features and filesystem substitutions that a normal review workflow can accidentally preserve. |
| Reversal path | None for the fail-closed requirements. If a future Pages payload intentionally needs executable files or symlinks, extend the manifest schema with explicit per-target type and mode declarations and add matching regressions. |
| Evidence | `scripts/stage_gh_pages_io_docs.py`; exploit regressions in `tst/test_suite/io/test_stage_gh_pages_io_docs_cpu.py`; local `22 passed`; strict detached stage, verify, HTML, linkcheck, and reviewed-drift packet generation against `origin/gh-pages`. |
| Follow-up | Keep the exploit regressions in the final matrix and repeat staging against refreshed `origin/gh-pages` after code merge. |

### D-100: Make New Reader Contracts Fail Closed Without Breaking Legacy Rank Fixtures

| Field | Value |
| --- | --- |
| Status | Accepted during `RCP-10` red-team correction |
| Decision | Require explicit distribution, shard ID, sibling count, and MeshBlock count metadata for new node-sharded `.bin` and `.cbin` files. Require explicit distribution, shard ID, and sibling count metadata for every new rank- or node-sharded `sphslice` file. Require every intrinsic `AKPDFV2` dense or sparse payload to have a matching header declaration. Reject duplicate scalar metadata in binary preheaders, PDF headers, PDF dimension aliases, and spherical-slice headers. Preserve legacy inventory-free rank-sharded binary fixtures and transitional unversioned sparse PDF compatibility. |
| Reason | The red-team audit demonstrated that relocating a shared binary fixture into a `node_########` directory, stripping a sparse PDF V2 header declaration, or omitting spherical-slice shard counts could downgrade new formats into less-audited compatibility paths. New layouts have no historical inventory-free files to preserve. Legacy per-rank binary fixtures and transitional unversioned sparse PDF readers do have an intentional compatibility purpose. |
| Alternatives | Require inventory metadata for every historical rank-sharded binary file; remove transitional unversioned sparse PDF support; treat duplicate keys as last-write-wins; retain optional inventory for new formats. |
| Why not | Breaking historical rank fixtures is unnecessary. Transitional sparse support remains bounded when intrinsic V2 payloads cannot masquerade as unversioned files. Last-write-wins metadata and optional new-format inventory create ambiguous parsing contracts. |
| Reversal path | Narrow a retained compatibility row only through an explicit migration decision with affected-reader evidence. Do not relax metadata requirements for new producer layouts. |
| Evidence | `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; downgrade, duplicate-key, and positive shard regressions in `tst/test_suite/io/test_python_io_readers_cpu.py`; corrected reader suite returned `197 passed`. |
| Follow-up | Obtain a fresh independent file-format and Python-reader audit against the settled tree. |

### D-101: Account For Text Metadata And Preserve Suffix-Only Batch Conversion

| Field | Value |
| --- | --- |
| Status | Accepted during `RCP-10` red-team correction |
| Decision | Charge conservatively decoded textual metadata against `max_live_bytes` before parser-controlled materialization using a four-times byte charge. Carry retained charges through later reconstruction checks. Keep public return dictionaries free of internal accounting fields. In `make_athdf.py`, replace only the final `.bin` suffix when naming `.athdf` output. |
| Reason | The API red-team audit reproduced megabyte-scale textual metadata accepted under a `64 KiB` live budget and demonstrated that the batch wrapper corrupted identifiers containing an interior `.bin`. The live budget must cover scalable text expansion as well as array allocations. |
| Alternatives | Rely only on header-read limits; exempt a fixed text allowance from `max_live_bytes`; use interpreter-specific `sys.getsizeof()` accounting for every object; keep global string replacement in the batch helper. |
| Why not | Header-read limits do not bound live expansion. A fixed allowance weakens the advertised aggregate budget. Exact interpreter object accounting is brittle and still needs conservative transient assumptions. Global replacement changes user-visible basenames. |
| Reversal path | Replace the conservative multipliers only with an audited parser-specific accounting model that remains admission-before-materialization and preserves keyword-only budget overrides. |
| Evidence | `vis/python/io_reader_common.py`; `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; `vis/python/make_athdf.py`; large-text and interior-suffix regressions in `tst/test_suite/io/test_python_io_readers_cpu.py`; corrected reader suite returned `197 passed`. |
| Follow-up | Keep the malformed-input regressions in the full local matrix and fresh Python API audit. |

### D-102: Partition Energy Channels By Signed Flux And Reject Non-Positive PDF Mass

| Field | Value |
| --- | --- |
| Status | Accepted during `RCP-10` red-team correction |
| Decision | Define `edot_sph_out`, `edot_sph_in`, `edot_vert_out`, and `edot_vert_in` as positive and negative partitions of signed total energy flux. For MHD, include the Poynting contribution before partitioning. For PDF `weight = mass`, require positive finite conserved density before multiplying by cell volume. Preserve finite signed `weight = variable` values. |
| Reason | In MHD, the Poynting contribution can oppose the gas velocity, so gas-motion-selected energy channels contradict the `_out` and `_in` names. A mass histogram has no defensible meaning for zero, negative, or non-finite conserved density, even if ordinary conversion floors make the path uncommon. |
| Alternatives | Keep gas-motion-selected transport channels and rename or document them; classify by gas velocity only for Hydro; silently floor PDF mass density; reject signed variable weights. |
| Why not | Flux-sign partitioning is consistent across Hydro and MHD and matches the public names. Silent floors alter diagnostics. Finite signed variable weighting is a useful deliberate capability distinct from mass weighting. |
| Reversal path | Add separately named gas-motion-selected diagnostics later if a concrete analysis workflow needs them. Do not overload the signed flux partitions. |
| Evidence | `src/outputs/diagnostic_semantics.hpp`; `src/outputs/derived_variables.cpp`; `src/outputs/pdf.cpp`; direct harness, adversarial radial and vertical MHD regressions, and non-positive mass-density regressions in `tst/test_suite/io`; focused diagnostic set returned `5 passed`; broader serial producer set returned `75 passed` before the final mass-density rows were added. |
| Follow-up | Rerun the full producer matrix, refresh deferred Pages staging, and obtain a fresh numerical-semantics audit. |

### D-103: Bind Node-PDF Payload Identity And Reject Zero Restart Segments

| Field | Value |
| --- | --- |
| Status | Accepted during `RCP-10` red-team correction |
| Decision | Add `payload_rank` to every newly written node-sharded `AKPDFV2` header and require it to match the leader rank embedded in the binary preamble. Reject duplicate node payload-rank inventories. In node-restart publication, reject non-positive rank segment counts before emitting a manifest so writer and loader enforce the same positive-segment contract. |
| Reason | Separate node headers and PDF payloads could previously be swapped while each file remained individually parseable. The restart loader already rejects zero-count segments, but the publisher accepted them even though AthenaK load balancing intends each participating rank to own at least one MeshBlock. Both boundaries should fail closed before reconstruction or publication. |
| Alternatives | Encode node ID into the existing payload-rank field; checksum each PDF payload; accept swappable node payloads because ordinary writes do not move files; keep the restart writer looser than the loader. |
| Why not | The preamble field already has a stable writer-rank meaning, so adding matching header metadata is the narrow compatibility-preserving fix. Checksums are a larger format revision. Reader and writer disagreement creates avoidable malformed artifacts. |
| Reversal path | Introduce payload checksums only through a future explicit format-version decision. Keep `payload_rank` validation for V2 node shards. |
| Evidence | `src/outputs/pdf.cpp`; `src/outputs/restart.cpp`; `vis/python/read_pdf.py`; swapped-node-payload regression in `tst/test_suite/io/test_python_io_readers_cpu.py`; corrected reader suite `230 passed`; serial and MPI Debug rebuilds passed. |
| Follow-up | Include node-PDF and restart publication behavior in fresh file-format and restart/MPI correction audits. |

### D-104: Expand Deferred Pages Staging To Nine Bounded Files

| Field | Value |
| --- | --- |
| Status | Accepted during `RCP-10` documentation correction |
| Decision | Supersede only the eight-file allowlist portion of `D-082`. Add one bounded replacement for `docs/source/modules/index.md` so the Outputs landing-page row reports 13 registered formats. Require the staging helper to prove an exact nine-file payload while preserving the detached-worktree, hash, marker, idempotence, build, linkcheck, rendered-inspection, and write-free drift-packet requirements. |
| Reason | The deferred Outputs module page advertises the added IO format while the live modules landing page still says 12 registered formats. Leaving that row untouched produces an internally inconsistent public site. |
| Alternatives | Keep the stale count; remove the count from the landing page; replace the whole modules index; edit live Pages directly after code merge. |
| Why not | A stale count is user-visible drift. Removing or broadly replacing unrelated live documentation is outside this branch. A bounded hash-checked row replacement preserves the deferred publication lifecycle. |
| Reversal path | If the live Pages row changes before publication, run reviewed-drift mode and deliberately update the bounded section fingerprint before restaging. |
| Supersedes | The eight-file payload count in `D-082`; all other `D-082` requirements remain controlling. |
| Evidence | `deferred_docs/gh-pages/io-output-formats-and-sharding/manifest.json`; bounded modules-index fragment; helper suite `22 passed`; exact nine-file detached strict stage; `--verify-staged --run-builds`; empty linkcheck output; reviewed-drift packet generation outside the target; rendered inspection of the 13-format row, IO example route, `payload_rank` prose, and preserved home-page iframe. |
| Follow-up | Record fresh nine-file stage, build, linkcheck, drift-packet, and navigation evidence in the ledger. |

### D-105: Require Immutable External Qualification Packets

| Field | Value |
| --- | --- |
| Status | Accepted during `RCP-10` process-packaging correction |
| Decision | Require one durable immutable evidence-packet root per external qualification run family, including copied inputs, overrides, commands, environment, raw logs, generated inventories, packet index, and a SHA-256 manifest covering retained artifacts. Record the archive root and checksum-manifest digest in every external evidence row. |
| Reason | Scheduler scratch paths and terminal summaries can disappear or become unauditable before independent review. The external CUDA, topology, filesystem, and scaling gates need evidence that survives handoff. |
| Alternatives | Retain cluster-local paths only; paste summarized output into the ledger; archive logs without checksums or copied input decks. |
| Why not | Those choices cannot prove what was run or whether retained artifacts changed after collection. |
| Reversal path | A site-specific archival system may replace the directory layout only if it preserves immutable indexed artifacts and checksums. |
| Evidence | `IO_EXTERNAL_IO_QUALIFICATION_PLAN.md`; process-packaging auditor finding; frozen scaling deck checksum in `IO_RESTART_MANIFEST_SCALING_QUALIFICATION.md`. |
| Follow-up | Use the packet contract during external `RCP-08B` and `RCP-09`; do not claim merge readiness until independent auditors accept both lanes. |

### D-106: Close Late Reader And Restart Compatibility Apertures

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-10` correction |
| Decision | Reject non-positive restart payload block counts as well as non-positive segment counts. Preserve transitional unversioned sparse PDF compatibility for shared and historical rank layouts, but reject unversioned node shards because node sharding is a V2 producer feature. Enforce conservative binary metadata admission before callers decode or split records and legacy PDF row admission before ASCII decoding. Keep parser accounting fields internal to the canonical binary reader. Route `make_athdf.py` through `bin_convert.convert_file`. |
| Reason | A synthetic manifest could add unused zero-block header-only payloads that the producer cannot emit. Relocating an unversioned sparse PDF below `node_########/` bypassed V2 inventory and payload identity checks without preserving a real historical format. Some metadata paths still decoded before low-memory rejection, public readers leaked private bookkeeping, and the batch wrapper duplicated canonical conversion steps. |
| Alternatives | Document zero-block restart payloads; retain weak unversioned node-PDF reconstruction; treat decode-before-reject as acceptable because header byte limits exist; expose accounting fields as public API; keep duplicate batch-wrapper conversion logic. |
| Why not | Those choices preserve noncanonical layouts or implementation leakage without a user requirement. The narrower contract preserves actual historical compatibility while making new node formats fail closed. |
| Reversal path | Broaden transitional node-PDF support only with a frozen historical producer artifact and explicit migration tests. Introduce empty restart payloads only through a versioned manifest decision. |
| Evidence | `src/restart_manifest.cpp`; `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; `vis/python/make_athdf.py`; zero-payload, unversioned-node, duplicate-payload-rank, public-accounting, package-import, and low-budget regressions in `tst/test_suite/io`; focused readers `236 passed`; corrupted-manifest subset `29 passed`. |
| Follow-up | Require fresh restart/file-format and Python API audits against the settled tree. |

### D-107: Keep Deferred Pages Markers Outside Markdown Table Rows

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-07` correction |
| Decision | Supersede the row-level landing-page edit portion of `D-104`. Replace the full bounded Support Systems section in `docs/source/modules/index.md`, preserving every neighboring row while changing the Outputs count to 13. Replace the full bounded Implementation Entry Points section in `docs/source/modules/outputs.md` so the refreshed PDF row remains inside its Markdown table. Require semantic helper regression and rendered browser checks for both tables. |
| Reason | The row-level replacement inserted marker comments between Markdown rows. Sphinx and linkcheck passed, but rendered HTML split the table and dropped neighboring navigation. The retained Outputs implementation table also continued to understate PDF dimensionality. |
| Alternatives | Remove bounded markers for the row; hand-edit Pages after merge; tolerate split rendered HTML; replace the whole module index file. |
| Why not | Marker-free edits weaken idempotence, manual repair is not reproducible, broken navigation is user-visible, and whole-file replacement is broader than required. A bounded section preserves unrelated content and places comments outside row flow. |
| Reversal path | If the live Support Systems section changes, generate reviewed drift, reconcile its bounded fingerprint deliberately, and repeat rendered checks before publication. |
| Supersedes | The row-level `docs/source/modules/index.md` replacement shape in `D-104`; the nine-file payload and all other `D-104` requirements remain controlling. |
| Evidence | `manifest.json`; Support Systems and Implementation Entry Points fragments; staging-helper semantic regression; helper suite `23 passed` before the final table-width correction; strict detached HTML and linkcheck builds; rendered browser record retained under `/tmp/athenak-io-robust-logs-pages/browser.log`. |
| Follow-up | Obtain a fresh independent deferred-Pages audit and repeat staging after code merge before opening a Pages review. |

### D-108: Enforce Reader Budgets Before Materialization And Validate PDF Weights

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-10` re-audit correction |
| Decision | Derive binary, spherical-slice, and legacy-PDF textual `readline()` caps from the remaining live-memory allowance before yielding bytes for decoding. Stream declared embedded athinput regions one line at a time, decode retained records once, and count retained list and string objects. Reject non-finite dense, sparse, aggregated, and legacy PDF histogram values while preserving finite signed values. Extend `make_athdf.py` to expose canonical `.bin`/`.cbin`, shard-assembly, and reader-budget behavior. Explicitly reject non-positive root-side restart payload totals before manifest publication. |
| Reason | Header and payload byte caps bound absolute reads but did not strictly honor a smaller configured live-memory ceiling before materialization. Newline-heavy athinput content amplified Python containers beyond a byte-only estimate. Reader acceptance of NaN or infinity contradicted writer-side diagnostic validation. The batch wrapper and restart publisher should enforce the same contracts as their canonical implementations and consumers. |
| Alternatives | Treat bounded transient overshoot as acceptable; retain byte-only embedded-header accounting; reject negative PDF histogram values as well as non-finite values; document `make_athdf.py` as a narrower legacy wrapper; rely on restart topology invariants rather than an explicit payload-total guard. |
| Why not | Those choices weaken the advertised practical-limit contract, reject valid signed variable-weight diagnostics, preserve unnecessary wrapper drift, or leave writer-loader symmetry implicit. |
| Reversal path | Replace conservative text charges only with an audited admission-before-materialization model. Broaden payload validation only through an explicit format decision and compatibility evidence. |
| Evidence | `vis/python/bin_convert.py`; `vis/python/read_sphslice.py`; `vis/python/read_pdf.py`; `vis/python/make_athdf.py`; `src/outputs/restart.cpp`; focused low-budget, newline-heavy, non-finite, signed-weight, batch-wrapper, and guard-order regressions in `tst/test_suite/io`. |
| Follow-up | Rerun the full local matrix and obtain fresh Python API, restart/MPI, compatibility, and whole-branch audits against the settled committed tree. |

### D-109: Keep CUDA As The Frozen Device Lane And Require Explicit HIP Evidence For HIP Deployment

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-09` qualification-plan clarification |
| Decision | Preserve CUDA as the device-execution lane required by the frozen robustification guide. Do not infer HIP readiness from CPU smoke or CUDA evidence. If a production deployment uses HIP, require a separate immutable HIP evidence packet with the same representative regressions, environment inventory, and independent audit before making a HIP production-readiness claim. |
| Reason | This workstation exposes neither CUDA nor HIP, while the frozen guide explicitly specifies CUDA qualification. Stating the HIP boundary prevents an ambiguous device-readiness claim without expanding the frozen guide retroactively. |
| Alternatives | Add HIP as a mandatory merge gate for every deployment; claim HIP readiness from Kokkos portability; omit HIP entirely from the qualification plan. |
| Why not | A universal HIP gate expands scope without a named deployment requirement. Backend portability is not execution evidence. Silence leaves the deployment boundary unclear. |
| Reversal path | Promote HIP into a required merge gate if maintainers identify a HIP production target for this branch. |
| Evidence | `IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md`; `IO_EXTERNAL_IO_QUALIFICATION_PLAN.md`; local environment probe recorded in `IO_FEATURE_AUDIT_LEDGER.md`. |
| Follow-up | External operators must retain CUDA evidence and, when applicable, a separate HIP packet before the corresponding production-readiness claim. |

### D-088: Superseded By D-094

| Field | Value |
| --- | --- |
| Status | Superseded by `D-094` |
| Decision | Preserve this navigable record for the earlier `RCP-03` question about checked coarsened-binary kernel-range and normalization-range representation. |
| Reason | `D-094` resolved the question by selecting signed checked `std::int64_t` launch arithmetic and separate `std::size_t` allocation and serialized-payload preflights while narrowing the producer to uniform 3D active-zone-only output. |
| Alternatives | Use `int` launch ranges; use `std::size_t` directly for kernel arithmetic; broaden the producer contract before qualification. |
| Why not | Those alternatives either narrow representability incorrectly, mix signed launch arithmetic with allocation domains, or expand unsupported output semantics. |
| Reversal path | Promote one additional producer row at a time only with explicit writer, reader, conversion, equality, and documentation evidence. |
| Evidence | `D-094`; `src/outputs/coarsened_binary.cpp`; `src/outputs/coarsened_binary_layout.hpp`; `tst/test_suite/io/cbin_layout_harness.cpp`. |
| Follow-up | None beyond the controlling `D-094` qualification rows. |

### D-110: Parse Fixed Binary Grammar Without Token Expansion And Preallocate Shard Assembly

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-10` Python re-audit correction |
| Decision | Parse fixed binary signatures, scalar records, and preheader `key=value` records with bounded exact grammar rather than unconstrained token lists. Before reconstructing partitioned `.bin` or `.cbin`, account conservatively for the transient duplicate-logical-ownership dictionary, release it after validation, and copy shard rows directly into preallocated NumPy arrays. Preserve direct `make_athdf.main(file_stem=..., verbose=...)` compatibility by defaulting additive assembly and reader-limit options. |
| Reason | A malformed whitespace-heavy metadata record could expand into Python token containers beyond the configured live-memory cap before later rejection. Shard assembly budgeted duplicated NumPy payloads but materialized uncharged Python row lists and logical-owner bookkeeping. The batch wrapper's additive options accidentally became mandatory for historical direct callers. |
| Alternatives | Pre-count every fixed-grammar token split; add broad conservative charges while retaining Python row-list assembly; remove callable-wrapper compatibility; rely on byte caps alone. |
| Why not | Exact grammar avoids unnecessary materialization. Direct array assembly is simpler to bound than interpreter-specific row-list amplification. The callable wrapper is an existing low-cost compatibility surface. Byte caps alone do not prove the live-memory contract. |
| Reversal path | Replace conservative ownership accounting only with an audited lower-overhead duplicate detector whose peak is admitted before construction. Remove direct-wrapper defaults only through an explicit API migration decision. |
| Evidence | `vis/python/bin_convert.py`; `vis/python/make_athdf.py`; whitespace-amplification, ownership-peak, empty-shard, and direct-wrapper regressions in `tst/test_suite/io/test_python_io_readers_cpu.py`; focused reader module returned `228 passed`. |
| Follow-up | Rerun the complete local matrix, capture an immutable local packet, and require fresh independent Python plus whole-branch acceptance audits against the committed tree. |

### D-111: Preserve Historical Spherical-Slice Fallback And Bound Python Metadata Objects

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-10` Python correction re-audit |
| Decision | Restore narrow historical spherical-slice version-1 inference for shared and rank layouts that use `single_file_per_rank`, while requiring complete metadata for every new explicit rank/node layout. Stream post-discovery shard validation without rebuilding unbounded identifier containers. Conservatively preflight direct binary logical-owner sets, athdf-like slice-classification sets, sorted shard inventories, `make_athdf.py` batch inventories, modern and legacy PDF metadata records, and spherical-slice metadata records. Preallocate legacy dense PDF payload rows. Preserve the canonical direct-reader `mb_data[variable]` NumPy array contract. |
| Reason | Version-1 fallback is an established reader compatibility surface; new node formats still need fail-closed inventories. Text-byte charges alone do not cover high-cardinality Python metadata objects, and bounded file payloads do not justify unaccounted shard lists, tuple sets, or legacy row arrays. The ndarray result matches the canonical module documentation and avoids row-list amplification. |
| Alternatives | Break all inventory-free spherical-slice files; retain broad fallback for new node layouts; accept interpreter-object overhead as negligible; expose lists of per-MeshBlock arrays for direct readers; keep unbounded `glob` inventories in the batch wrapper. |
| Why not | Those choices either discard useful compatibility, weaken new-format validation, or leave scalable live-memory gaps. Restoring list results would conflict with the documented array contract and reintroduce avoidable container growth. |
| Reversal path | Remove historical fallback only through a versioned migration decision with fixture evidence. Replace conservative record charges only with an audited admission model. Add a list-like compatibility adapter only if a real external caller requires mutation semantics that ndarray indexing does not preserve. |
| Evidence | `vis/python/io_reader_common.py`; `vis/python/bin_convert.py`; `vis/python/read_pdf.py`; `vis/python/read_sphslice.py`; `vis/python/make_athdf.py`; high-cardinality, legacy-fallback, bounded-inventory, direct-validation, slice-classification, preallocated-row, and ndarray-contract regressions in `tst/test_suite/io/test_python_io_readers_cpu.py`; focused reader and hardening suite returned `270 passed`. |
| Follow-up | Require fresh Python API audit acceptance, rerun the canonical local matrix, refresh detached Pages hashes, and include these boundaries in the committed-tree evidence packet. |

### D-112: Use Analytical And Binary-Backed Spherical-Slice Producer Oracles

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-10` numerical re-audit |
| Decision | Keep a cycle-zero analytical shock-tube oracle that forces spherical-slice interpolation across an internal MeshBlock face, strengthen the MPI shared/rank/node comparison to force the same cross-rank face interpolation, and validate adaptive `Rebuild()` behavior against an independent same-cycle binary-snapshot trilinear oracle after redistribution. |
| Reason | Shape, finiteness, and cross-layout equality can all pass when interpolation weights are stale or when all compared layouts share the same face-sampling defect. The independent oracles prove the numerical producer behavior rather than only the wire format. |
| Alternatives | Retain metadata-only spherical-slice round trips; compare shared/rank/node output only; inspect adaptive output at cycle one before redistribution; use the spherical-slice implementation itself to compute expected values. |
| Why not | Those checks miss ghost-zone face interpolation, stale adaptive ownership, or common-mode numerical defects. The selected shock-tube and cycle-two adaptive cases isolate both risks with deterministic local evidence. |
| Reversal path | Replace the fixtures only with equally independent numerical oracles that still force cross-face and post-redistribution interpolation. |
| Evidence | `tst/test_suite/io/test_output_formats_cpu.py`; `tst/test_suite/io/test_output_formats_mpicpu.py`; focused serial spherical-slice subset returned `12 passed`; strengthened MPI output-format module returned `1 passed`. |
| Follow-up | Retain both numerical regressions in the canonical committed-tree matrix and request final numerical audit acceptance. |

### D-113: Treat Public Node-Restart Manifest Rename As The Commit Point

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-10` restart/MPI correction re-audit |
| Decision | Keep published node-restart payloads rollback-owned until the public manifest rename is attempted. Restore rollback ownership if manifest publication fails. Once the public manifest rename succeeds, treat the checkpoint as committed: later injected failures and generation-reservation removal failures may report an error but must preserve the manifest and declared payloads as a resumable checkpoint. |
| Reason | A visible public manifest is the reader-visible commit record. Removing payloads or the manifest after that point can turn a successfully published checkpoint into a broken public artifact. Precommit rollback and postcommit preservation give the transaction one unambiguous boundary. |
| Alternatives | Roll back the public manifest and payloads after any later cleanup failure; ignore reservation-removal failure; publish payloads without a public commit boundary. |
| Why not | Postcommit rollback can invalidate a checkpoint already visible to readers. Ignoring cleanup failures hides operational defects. Payload-only visibility is not a supported restart entry point. |
| Reversal path | Change the transaction protocol only through a versioned design with crash-recovery and production-filesystem qualification. |
| Evidence | `src/outputs/restart.cpp`; injected `after_manifest_publication` and `reservation_removal` regressions plus one-rank multi-payload span resume in `tst/test_suite/io/test_node_sharding_mpicpu.py`; focused restart subset returned `7 passed`. |
| Follow-up | Retain rank-separated multi-node precommit and postcommit failure-injection evidence in external `RCP-09`. |

### D-114: Tighten Spherical Headers And Document Native-Endian New Payloads

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-10` file-format correction re-audit |
| Decision | Require explicit modern spherical-slice layout metadata, positive radius, at least two angles on each axis, bounded `npoints`, and strict header-only admission. Preserve only the narrow synthetic version-1 shared/rank fallback. Document that new PDF V2 and spherical-slice binary scalars are host-native endian rather than silently changing bytes late in this branch. |
| Reason | Lightweight header APIs should fail closed on intrinsically invalid metadata before payload reads. A byte-order migration is a wire-format change and needs an explicit version rather than an implicit late correction. Historical fallback claims must match available provenance evidence. |
| Alternatives | Infer layout for all modern shards; accept degenerate angular axes or non-positive radii; switch existing new payload bytes to a fixed order without versioning; claim a historical fixture that has not been frozen. |
| Why not | Those choices weaken new-format admission, risk compatibility ambiguity, or overstate evidence. |
| Reversal path | Add a versioned fixed-byte-order format with same-platform and cross-endian fixtures. Freeze a predecessor artifact before upgrading synthetic fallback evidence to provenance-qualified compatibility evidence. |
| Evidence | `vis/python/read_sphslice.py`; strict malformed-header tests in `tst/test_suite/io/test_python_io_readers_cpu.py`; `IO_FORMAT_COMPATIBILITY.md`; deferred Pages fragments. |
| Follow-up | Retain reader regressions and require final file-format audit acceptance. |

### D-115: Require Attributable Filesystem Measurements And Reproducible External Packets

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-08B` and `RCP-09` plan correction |
| Decision | Preserve logical restart-validation pressure as mandatory structured telemetry, but require a separate attributable observed filesystem-read measurement on the target platform. Mark the filesystem-amplification row incomplete when reliable attribution is unavailable. Require rank-separated scheduler logs, exact launcher templates, per-row timeouts, and a two-level checksum contract: `artifacts.sha256` excludes itself and the packet index; the index records the inner digest; the ledger or archive record holds an outer index digest. |
| Reason | Logical pressure is useful protocol telemetry but cannot prove physical read amplification under caches, read-ahead, and parallel-filesystem behavior. Reproducible packet structure is required for an independent external audit. |
| Alternatives | Treat logical pressure as observed IO; accept summary-only scheduler evidence; allow self-referential checksum claims; drop the filesystem row when platform telemetry is inconvenient. |
| Why not | Those choices overclaim performance evidence or make the packet impossible to reproduce and verify. |
| Reversal path | Replace the observed-read row only with a more direct platform-supported filesystem measurement and an explicit auditor-approved definition. |
| Evidence | `IO_EXTERNAL_IO_QUALIFICATION_PLAN.md`; `IO_RESTART_MANIFEST_SCALING_QUALIFICATION.md`; final test-evidence audit. |
| Follow-up | External operators must complete the templates and retain raw accounting logs before recording `RCP-08B` or `RCP-09` acceptance. |

### D-117: Make The External Scheduler Deck Executable

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-08B` and `RCP-09` plan re-audit correction |
| Decision | Add `scripts/run_external_io_qualification_slurm.sh` as the normative Slurm launcher for ED-1, MR-1, MR-2, and representative scaling resumes. Require unique row/sample logs, exact ED-1 and MR-1 placement, a three-line arbitrary-distribution hostfile for asymmetric MR-2 placement, one unmeasured warm-up plus five measured samples, measured-only restart telemetry, monotonic launcher intervals, exit-code and timeout records, and a site-filled attributable-read accounting hook. |
| Reason | Prose templates did not execute the preregistered sampling contract and could overwrite evidence or silently use the wrong asymmetric topology. |
| Alternatives | Leave launcher assembly to each external operator; retain uniform `--ntasks-per-node` examples for MR-2; permit unindexed repeated samples. |
| Why not | Those choices allow irreproducible evidence and cannot prove the frozen topology or measurement protocol. |
| Reversal path | Replace the Slurm script with a scheduler-specific equivalent only if it preserves the same row, topology, sample, timing, filesystem-accounting, and packet-index contract. |
| Evidence | `scripts/run_external_io_qualification_slurm.sh`; `IO_EXTERNAL_IO_QUALIFICATION_PLAN.md`; `IO_RESTART_MANIFEST_SCALING_QUALIFICATION.md`; final process-evidence re-audit. |
| Follow-up | Syntax-check the script locally and require an external operator plus independent auditor to review the site-filled accounting hook before consuming an allocation. |

### D-116: Exercise Mixed-Level Spherical Interpolation And Serialized Narrowing

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-10` numerical correction re-audit |
| Decision | Stage shared spherical-slice values as serialized floats before publication and reject any finite in-memory value that becomes non-finite during narrowing. Add a mixed-level AMR regression with coarse and refined blocks, synchronized ghost-zone binary output, an independent owner-block ghost-snapshot interpolation oracle, and an assertion that sampled stencils cross AMR levels. |
| Reason | A finite `double` can overflow to a non-finite float on disk. The earlier adaptive oracle refined every block uniformly, so it did not exercise coarse-fine ghost stencils even though it validated post-redistribution rebuilding. |
| Alternatives | Let the reader reject serialized infinity after publication; retain the all-refined adaptive oracle; compare only spherical-slice layouts against each other. |
| Why not | Those choices publish invalid files or miss a numerically distinct coarse-fine path. |
| Reversal path | Replace the oracle only with an equally independent coarse-fine numerical check and preserve pre-publication serialized-value admission. |
| Evidence | `src/outputs/spherical_slice.cpp`; `tst/test_suite/io/test_output_formats_cpu.py`; focused spherical-slice subset returned `14 passed`. |
| Follow-up | Run this producer lane on CUDA and any deployment HIP backend during external qualification. |

### D-118: Keep Header-Only Admission Strict And Claims Narrow

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-10` file-format re-audit correction |
| Decision | Make public `read_pdf_header()` reject AKPDFV2 sparse metadata unless it declares a rank or node distribution, required same-family inventory fields, a path-bound shard identity, an in-range shard ID, no opposite-family inventory fields, and node-only `payload_rank` where required. Preserve internal full-reader support for a root metadata copy because full assembly validates the actual payload path separately. Document header-only APIs as single-header declaration checks; reserve sibling completeness and consistency claims for full readers. |
| Reason | The public header API is advertised for bounded preflight and must fail closed without a payload read. Full assembly still needs compatibility with metadata copies stored above shard directories. Documentation must describe the boundary the code actually enforces. |
| Alternatives | Let malformed sparse distributions pass until payload loading; require all internal metadata copies to reside inside shard directories; claim that header-only calls discover sibling families. |
| Why not | Those choices weaken preflight, remove a supported full-reader metadata layout, or overstate API behavior. |
| Reversal path | Add an explicit family-header API if sibling-only inventory inspection becomes a user requirement. Keep `read_pdf_header()` bounded to one metadata artifact. |
| Evidence | `vis/python/read_pdf.py`; strict sparse-header regressions in `tst/test_suite/io/test_python_io_readers_cpu.py`; staged Visualization and File Reference fragments; helper smoke assertions. |
| Follow-up | Retain reader and staged-doc regressions in committed-tree qualification and require fresh file-format audit acceptance. |

### D-119: Index Every External Launch Terminally And Isolate Retries

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-09` scheduler-deck re-audit correction |
| Decision | Use packet-local `ATTEMPT_ID` tokens, Slurm `%J` job-step log IDs, duplicate launch-key rejection, and attempt-scoped output, inventory, accounting, stdout, and stderr paths. Index rank-map launches and every attempted Athena launch only after accounting and forbidden-`.assembled` validation determine the terminal disposition. Preserve Bash 3 compatibility with a portable MR-2 hostfile loop. |
| Reason | Repeated Slurm steps and retries must not overwrite evidence. Accounting-hook failures and forbidden staging artifacts must not disappear from the packet index or remain mislabeled as passed. |
| Alternatives | Retain `%j` job-only log names; rely on operators to avoid retries; require Bash 4 for `mapfile`; append success before artifact validation. |
| Why not | Those choices allow evidence collisions, incomplete failure records, avoidable platform incompatibility, or false pass rows. |
| Reversal path | Replace this runner only with a scheduler-specific equivalent that preserves unique launch keys and one truthful terminal index row per attempt. |
| Evidence | `scripts/run_external_io_qualification_slurm.sh`; local mock lifecycle matrix under Bash `3.2.57`; `IO_EXTERNAL_IO_QUALIFICATION_PLAN.md`; `IO_RESTART_MANIFEST_SCALING_QUALIFICATION.md`; fresh process-evidence re-audit. |
| Follow-up | Require an external operator and independent auditor to inspect the site-filled filesystem-accounting hook before consuming a production allocation. |

### D-120: Time Launches In One Process And Scan Restart Sidecars

| Field | Value |
| --- | --- |
| Status | Accepted during final `RCP-09` scheduler-deck re-audit correction |
| Decision | Measure every launcher interval inside one Python timer process around the subprocess, retain its timing TSV, and reject malformed or negative interval records terminally. After every Athena launch, retain a forbidden-staging scan over the packet, output directory, and restart-manifest directory where applicable; reject both `*.assembled` and `*.assembled.tmp` before appending the terminal packet-index row. |
| Reason | Separate macOS Python processes can expose non-comparable monotonic epochs. Historical restart assembly staging lives beside the manifest and includes a temporary sidecar suffix, outside the output directory scanned by the earlier runner. |
| Alternatives | Assume cross-process monotonic comparability; scan only the `-d` output directory; reject only final `.assembled` files. |
| Why not | Those choices allow negative elapsed intervals or miss the precise unsupported restart-staging artifacts the qualification lane must disprove. |
| Reversal path | Replace the timer or scan only with a scheduler-specific equivalent that retains a nonnegative launcher interval and proves absence of both restart-sidecar suffixes. |
| Evidence | `scripts/run_external_io_qualification_slurm.sh`; `tst/test_suite/io/test_external_io_slurm_runner_cpu.py`; local macOS Python `3.9.6` and Bash `3.2.57` mock lifecycle matrix; fresh process-evidence re-audit. |
| Follow-up | Preserve the checked-in mock suite and require real scheduler evidence before closing external qualification. |

## Pending Decision Queue

Resolve these before merge readiness:

No pending local design decision remains. `RCP-08B` still requires external
scheduler-backed measurements before a keep, refactor, or separate-branch
disposition can be recorded.

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
