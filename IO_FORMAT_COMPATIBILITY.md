# IO Format Compatibility Contract

## Purpose

This file records compatibility decisions and evidence for
`feature/io-output-formats-and-sharding`. It is intentionally separate from the
implementation guide: this contract must be updated as behavior is implemented and
verified.

The implementation is based on `origin/main` commit `886dd2a1`. The Gotham branch is
historical evidence only. No implementation is sourced from
`origin/feature/single-file-per-node-outputs`.

## Compatibility Levels

| Level | Meaning |
| --- | --- |
| Byte-identical | Legacy-mode files remain byte-for-byte equal to baseline fixtures. |
| Header/schema-compatible | Reader-visible schema and payload meaning are unchanged; parameter dumps or nonsemantic bytes may differ. |
| New-reader compatible | The canonical new reader consumes baseline and new layouts; legacy readers are not promised to consume the new layout. |
| Deliberately incompatible | A documented new schema/version or new file family is introduced. |

## Format Matrix

This matrix starts with intended contracts. Each row must be updated with fixture
paths and verification results before the feature is considered complete.

| Output | Layout | Compatibility Decision | Evidence Required | Status |
| --- | --- | --- | --- | --- |
| `.bin` | Existing shared | Header/schema-compatible | Baseline fixture; canonical `bin_convert.py` readback; header/version assertions | Frozen fixture and canonical-reader readback passed |
| `.bin` | Existing per-rank | Header/schema-compatible | Baseline partitioned fixture and shard reconstruction | Frozen shard reconstruction equality passed |
| `.bin` | New per-node | New-reader compatible | Writer/readback equality against shared output; additive shard-inventory metadata and explicit empty-shard tests | One-node MPI aggregation/reconstruction equality passed; canonical reader validates complete dense node IDs and valid empty shards |
| `.bin` | Sliced per-node | New-reader compatible | Empty-owner and selected-slice reconstruction tests | Two-rank single-node test with one non-owning rank passed; explicit empty node-shard reader contract passed; multi-node empty-node qualification remains required |
| `.cbin` | Uniform 3D active-zone full-volume shared/per-rank | Header/schema-compatible | Baseline fixture; direct read; ATHDF/XDMF conversion; moments mode | Frozen readback/reassembly and conversion tests passed; generated moments output round-trips and converts; producer and reader reject invalid metadata before unsafe arithmetic |
| `.cbin` | Uniform 3D active-zone full-volume per-node | New-reader compatible | Shared/per-node numerical comparison; additive shard-inventory metadata; explicit empty-shard tests; moments mode | One-node MPI aggregation/equality and moments-mode conversion passed; canonical reader validates complete dense node IDs and valid empty shards |
| `.cbin` | Lower-dimensional, ghost-zone-expanded, static-refinement, or AMR output | Deliberately excluded | Construction-time rejection before directory creation | Unsupported configurations are rejected explicitly before publication |
| `.cbin` | Sliced output in any shard mode | Deliberately excluded | A sliced emitted extent is incompatible with supported coarsening | Sliced construction is rejected explicitly before data loading |
| `.pdf` | Existing one/two-dimensional text-form data and `.bins.pdf` companion | Header/schema-compatible | Frozen baseline fixture and reader round trip | One- and two-dimensional output byte-identical; reader tests passed |
| `.pdf` | New N-dimensional dense | New-reader compatible | Mixed-scale and weighting readback tests | Three-dimensional mass-weighted and four-dimensional scalar-/volume-weighted writer-readback passed; declared V2 preamble, cycle, finite-metadata, transformed-axis, retained-allocation, and invalid-scale checks passed; histogram accumulation is backend-portable and MPI reduction stages through host memory |
| `.pdf` | New sparse rank/node shards | New-reader compatible | Dense reconstruction equality, strict sibling-inventory metadata, checked per-file publication, practical read/allocation bounds, normalized aggregate metadata, and malformed-shard tests | Rank/shared and node/shared equality, canonical path/identifier/inventory validation, sibling V2-declaration agreement, required declared V2 preambles and matching cycles, bounded-read negatives, aggregate metadata normalization, malformed-shard tests, and temporary-file cleanup passed |
| `.sph.bin` | New shared/rank/node representations | Deliberately incompatible new file family; native state-backed scalar fields only | Writer/reader contract, cross-mode tests, derived-field rejection, explicit empty-shard inventory, layout metadata, bounded reads, normalized aggregate metadata, and atomic publication | Shared/rank/node writer-reader equality, malformed-shard tests, valid empty-shard coverage, layout/distribution rejection, whole-file read plus coordinate/offset bounds, aggregate metadata normalization, derived-field rejection, and temporary-file cleanup passed |
| Existing `file_type=sph` VTK surface output | Legacy `radius`/VTK behavior | Header/schema-compatible and unchanged | Regression output and parameter/filename check | Filename/parameter coexistence regression passed |
| `.rst` | Existing shared/per-rank | Header/schema-compatible | Restart round trips from baseline-compatible files | Generated and frozen `origin/main` shared/per-rank resume regressions passed; per-rank reader flag defect repaired |
| `.rst` | New per-node manifest/payload | Deliberately incompatible new manifest family; payload publication is transactional and restart loading is native and direct with no production `.assembled` staging | Restart/resume, malformed-manifest, terminal-checkpoint numbering, strict path/symlink/content-marker/inventory/replicated-header validation, checked publication, CP-02 forced-small-chunk direct-read coverage, and changed-rank-count resume | Manifest-only resume/numbering, exact/lexical/symlink/hard-link/copied payload-entry rejection, corrupt-marker rejection, path and symlink-escape rejection, bounded ordered inventory and segment negatives, replicated-header mismatch rejection, native no-`.assembled` loading, forced-small-chunk direct reads, and local-node rank-count change passed on two ranks on one node. Real multi-node routing remains a production qualification gate |

## Baseline Fixture Requirements

Frozen legacy fixtures must be generated from an untouched worktree at
`origin/main` commit `886dd2a1`, not from this branch after implementation changes.
A fixture manifest must record:

- producing commit SHA;
- build configuration and command;
- input deck and invocation;
- output filenames, sizes, and SHA-256 checksums; and
- the reader/check command used to establish the baseline.

Generated outputs from executable examples remain separate from immutable
compatibility fixtures.

The frozen fixture set and checksums are stored in
`tst/fixtures/io/origin_main_886dd2a1/README.md` and `SHA256SUMS`. It includes
serial shared and two-rank per-rank binary, coarsened-binary, and restart output,
plus the legacy text-form PDF family.

No baseline fixture will be claimed for `sphslice` or partitioned PDF data from
`origin/main`, because those layouts do not exist on that baseline. They are new
formats measured against their documented writer-reader contract and cross-mode
comparisons.

## Runtime Policy Compatibility

| Behavior | Decision | Rationale | Evidence |
| --- | --- | --- | --- |
| Final output default | `<time>/final_output_policy = all` | Preserve ordinary baseline final outputs unless opted out. | Serial regression passed. |
| Production final output | `restart_only` | Provide the valuable Gotham operational policy explicitly. | Serial regression passed. |
| No final output | `none` | Explicit expert option only. | Serial regression passed. |
| Final checkpoint counter advancement | Advance normally for every file actually written. | Prevent a resumed run from overwriting its terminal checkpoint. | Serial terminal-restart resume regression passed. |
| Output timing default | `<time>/output_timing = false` | Avoid added fences, reductions, and log records by default. | Serial regression passed. |
| Enabled output timing | Fence and report rank-maximum elapsed time per written output with event, block, type, and distribution labels. | Measures the slowest participant and distinguishes shared/rank/node output. | Serial and two-rank shared/rank/node MPI regressions passed. |

## Python API Compatibility

The final public binary conversion module is `vis/python/bin_convert.py`.
`bin_convert_new.py` is not a supported final API.

The canonical module must retain or deliberately migrate all in-repository uses,
including `vis/python/make_athdf.py` consumers of `read_binary`, `write_athdf`, and
`write_xdmf_for`.

| Capability | Compatibility Decision | Reason | Status |
| --- | --- | --- | --- |
| `read_binary` and existing assembly entry points | Preserve base public behavior while adding partitioned/node discovery and optional node inventory metadata additively. | Avoid breaking current analysis scripts while making new node-sharded output complete and auditable. | Frozen shared/rank fixtures, strict canonical shard directories, node metadata checks, explicit empty-shard assembly, duplicate-name and complete grid/logical/exact-logical-geometry rejection with zero-relative storage-aware bounded absolute tolerance, uniform emitted MeshBlock extents within each file, strict coarsening factors, bounded cumulative metadata and parameter reads, retained MeshBlock metadata and transient aggregate-copy budgets, incremental aggregate reconstruction preflight, omitted/mismatched ghost rejection including singleton axes, interior-width ghost placement with extended coordinates, intersecting-cell crop selection, deterministic uncovered levels, restored restriction, cumulative athdf-like allocation and NumPy coordinate/prolongation/restriction temporary preflight, and conversion tested. |
| `read_single_rank_binary_as_athdf` | Preserve the base signature unless a separately documented additive entry point is needed. | Gotham changed this API incompatibly. | Preserved in canonical module; destination-array and dtype contract plus direct signature audit complete; ghost-bearing or sliced emitted extents are rejected because the preserved signature has no ghost argument. |
| `write_athdf` | Retain in canonical `bin_convert.py`. | Base `make_athdf.py` and CLI require it; Gotham variants omit the definition incorrectly. | Preserved and conversion smoke-tested. |
| `write_xdmf_for` and `convert_file` | Retain supported conversion workflow and CLI; strip only the final `.bin`/`.cbin` extension when naming outputs. | Keep established ATHDF/XDMF user path usable without corrupting identifiers containing `bin`. | Preserved; API and CLI conversion tests passed after fixing output naming. |
| `read_pdf.py` | New public reader with dense/sparse node/rank reconstruction and strict validation. | New output format needs supported readback. | Legacy, transitional unversioned dense/sparse, V2 shared/rank/node, sibling V2-declaration agreement, required declared V2 preambles and matching cycles, mandatory V2 sparse inventories, finite metadata, bounded cumulative metadata/payload reads, explicit/generated/legacy edge tokenization and NumPy preflight, strict ASCII-only legacy numeric-row parsing with cumulative row and final-stack bounds, replacement-header retained-state accounting including legacy candidates, duplicate-validation peaks, per-sibling local-array release, inventory-summary retention, aggregate metadata normalization, and malformed shard inputs tested. |
| `read_sphslice.py` | New public reader with complete ownership and metadata checks. | New output format needs supported readback. | Shared/rank/node, strict sibling inventories, bounded whole-file reads plus cumulative header bytes, individual metadata lines, pre-split variable tokens, retained reference variable-metadata summary through coordinate generation, retained arrays, coordinate arrays, payload-copy and duplicate-validation peaks, embedded-header offsets including header-only API, incremental sparse combination with per-sibling array release, aggregate metadata normalization, and malformed shard inputs tested. |
| `bin_convert_new.py` | Do not add as public module or compatibility alias. | The branch removes redundant conversion APIs. | Accepted. |

## Diagnostic Semantics Boundaries

| Diagnostic surface | Supported contract | Evidence |
| --- | --- | --- |
| Generic PDF diagnostics `coord_*`, `mdot_*`, `edot_*`, `vel_*` | Available for single-fluid Hydro/MHD as documented. | Modern N-D PDF writer/readback tests pass. |
| Generic `mdot_*`, `edot_*`, `vel_*` in `<ion-neutral>` runs | Rejected until module-qualified two-fluid semantics exist. | Negative two-fluid PDF construction regression passes. |
| `sphslice` variables | Native state-backed scalar fields only; derived arrays are rejected until ghost-zone-safe interpolation exists. | Negative derived-slice construction regression passes. |

## Deferred Documentation Evidence

Candidate documentation is staged, not applied to the live Pages branch, at
`deferred_docs/gh-pages/io-output-formats-and-sharding/`. The bundle records
the inspected `origin/gh-pages` baseline, exact replacement/insertion target
paths, stale converter/PDF statements to remove, and the post-merge
application procedure.

Native direct node restart loading has replaced the earlier transient rank-0
`<manifest>.assembled` staging design. Do not publish or open a Pages
integration review until the code branch is merged and the staged pages are
reapplied, rebuilt, and re-audited against the then-current `origin/gh-pages`.

A detached temporary worktree at `origin/gh-pages`
`4833aa9341e19861297e330ff02aabfd8001935c` was used to apply the focused
overlays and validate the reference insertions as MyST includes. The command
`make clean html SPHINXOPTS="-W --keep-going"` completed successfully. No
Pages publication or live branch edit is part of this code feature branch.
