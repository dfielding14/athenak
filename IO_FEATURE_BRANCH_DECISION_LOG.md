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

## Pending Decision Queue

Resolve these before merge readiness:

| ID | Question | Required evidence | Owner checkpoint |
| --- | --- | --- | --- |
| D-072 | Select common MPI error-reporting helper and teardown policy. | Diff-driven MPI inventory and helper design. | RCP-02 |
| D-074 | Define checked directory creation, stale temporary handling, and restart generation collision policy. | Filesystem helper design and failure regressions. | RCP-02 |
| D-076 | Define diagnostic clamp, singularity, and zero-density semantics. | Numerical audit and analytic tests. | RCP-04 |
| D-077 | Reject or implement derived ghost-zone output. | Derived-variable audit and focused regression. | RCP-04 |
| D-078 | Select stable spherical-slice radius naming. | Collision analysis and compatibility review. | RCP-04 |
| D-079 | Bound the C++ parser and interface cleanup. | Scope audit after correctness checkpoints. | RCP-05 |
| D-080 | Decide whether to extract shared private Python reader utilities. | Duplication inventory and API review. | RCP-06 |
| D-082 | Select preservation-aware deferred Pages staging implementation. | Live `origin/gh-pages` blob inventory and idempotence design. | RCP-07 |
| D-084 | Select final process-artifact packaging. | Whole-branch review and PR usability audit. | RCP-10 |
| D-086 | Define widened output sequence rendering and counter domain. | Writer/reader inventory and compatibility tests. | RCP-02 |
| D-087 | Define deterministic output namespace collision rejection. | Target map inventory and compatibility tests. | RCP-02 |
| D-088 | Select checked `cbin` Kokkos range representation. | Kernel-range audit and boundary harness. | RCP-03 |
| D-089 | Select writer-side PDF and `sphslice` allocation limits. | Allocation map and reduced-cap tests. | RCP-04 |
| D-090 | Preregister production scaling topology and optimization threshold. | Scheduler qualification plan. | RCP-08 |

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
