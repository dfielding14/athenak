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
| `.cbin` | Uniform 3D active-zone full-volume shared/per-rank | Header/schema-compatible | Baseline fixture; direct read; ATHDF/XDMF conversion; moments mode | Frozen readback/reassembly and conversion tests passed; generated moments output round-trips and converts; producer and reader reject invalid metadata before unsafe arithmetic; every emitted axis extent must be divisible by `coarsen_factor` |
| `.cbin` | Uniform 3D active-zone full-volume per-node | New-reader compatible | Shared/per-node numerical comparison; additive shard-inventory metadata; explicit empty-shard tests; moments mode | One-node MPI aggregation/equality and moments-mode conversion passed; canonical reader validates complete dense node IDs and valid empty shards; every emitted axis extent must be divisible by `coarsen_factor` |
| `.cbin` | Lower-dimensional, ghost-zone-expanded, static-refinement, or AMR output | Deliberately excluded | Construction-time rejection before directory creation | Unsupported configurations are rejected explicitly before publication |
| `.cbin` | Sliced output in any shard mode | Deliberately excluded | A sliced emitted extent is incompatible with supported coarsening | Sliced construction is rejected explicitly before data loading |
| `.pdf` | Existing one/two-dimensional text-form data and `.bins.pdf` companion | Header/schema-compatible | Frozen baseline fixture and reader round trip | One- and two-dimensional output byte-identical; reader tests passed |
| `.pdf` | New N-dimensional dense | New-reader compatible | Mixed-scale and weighting readback tests | Three-dimensional mass-weighted and four-dimensional scalar-/volume-weighted writer-readback passed; mass weighting rejects non-finite and non-positive conserved density while finite variable weighting may be signed; readers reject non-finite histogram payloads while preserving finite signed values; declared V2 preamble and header declaration, cycle, finite-metadata, transformed-axis, retained-allocation, and invalid-scale checks passed; histogram accumulation is backend-portable and MPI reduction stages through host memory |
| `.pdf` | New sparse rank/node shards | New-reader compatible | Dense reconstruction equality, strict sibling-inventory metadata, checked per-file publication, practical read/allocation bounds, normalized aggregate metadata, and malformed-shard tests | Rank/shared and node/shared equality, canonical path/identifier/inventory validation, sibling V2-declaration agreement, required declared V2 preambles and matching cycles, non-finite payload rejection after sparse aggregation, bounded-read negatives, aggregate metadata normalization, malformed-shard tests, and temporary-file cleanup passed |
| `.sph.bin` | New shared/rank/node representations | Deliberately incompatible new file family; native state-backed scalar fields and native multi-field groups. Preserve narrow historical version-1 shared/rank fallback through `single_file_per_rank`; require complete metadata for every new explicit rank/node layout. | Writer/reader contract, cross-mode tests, synthetic historical-fallback regression, native-group readback, derived-field rejection, non-finite producer and forged-payload rejection, explicit empty-shard inventory, layout metadata, bounded reads, normalized aggregate metadata, and atomic publication | Shared/rank/node writer-reader equality, synthetic historical shared/rank fallback, native-group readback, malformed-shard tests, non-finite writer/readback rejection including finite-double-to-float narrowing overflow, valid empty-shard coverage, strict intrinsic dimensions/radius/point-count checks, mandatory explicit layout metadata, header-only admission checks, mixed-level AMR interpolation against a binary ghost-snapshot oracle, aggregate metadata normalization, derived-field rejection, and temporary-file cleanup passed |
| Existing `file_type=sph` VTK surface output | Legacy `radius`/VTK behavior | Header/schema-compatible and unchanged | Regression output and parameter/filename check | Filename/parameter coexistence regression passed |
| `.rst` | Existing shared/per-rank | Header/schema-compatible | Restart round trips from baseline-compatible files | Generated and frozen `origin/main` shared/per-rank resume regressions passed; per-rank reader flag defect repaired |
| `.rst` | New per-node manifest/payload | Deliberately incompatible new manifest family; payload publication is transactional, public-manifest rename is the commit point, and restart loading is native and direct with no production `.assembled` staging | Restart/resume, malformed-manifest, terminal-checkpoint numbering, strict path/symlink/content-marker/inventory/replicated-header validation, checked publication, precommit rollback, postcommit preservation, CP-02 forced-small-chunk direct-read coverage, and changed-rank-count resume | Manifest-only resume/numbering, exact/lexical/symlink/hard-link/copied payload-entry rejection, corrupt-marker rejection, path and symlink-escape rejection, bounded ordered inventory, non-positive payload-total and segment rejection in both writer and loader, replicated-header mismatch rejection, native no-`.assembled` loading, forced-small-chunk direct reads, one-rank multi-payload span routing, post-manifest-publication and reservation-removal failure preservation, and local-node rank-count change passed on two ranks on one node. Real multi-node routing remains a production qualification gate |

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

The version-1 spherical-slice fallback is covered by synthetic parser
regressions only. No provenance-qualified predecessor spherical-slice artifact
has been frozen, so the branch must not describe that fallback as historical
fixture evidence.

## Runtime Policy Compatibility

| Behavior | Decision | Rationale | Evidence |
| --- | --- | --- | --- |
| Final output default | `<time>/final_output_policy = all` | Preserve ordinary baseline final outputs unless opted out. | Serial regression passed. |
| Production final output | `restart_only` | Provide the valuable Gotham operational policy explicitly. | Serial regression passed. |
| No final output | `none` | Explicit expert option only. | Serial regression passed. |
| Final checkpoint counter advancement | Advance normally for every file actually written. | Prevent a resumed run from overwriting its terminal checkpoint. | Serial terminal-restart resume regression passed. |
| Performance timing default | `<time>/performance_timing = false` | Avoid added fences, reductions, and log records by default. | Serial regression passed. |
| Enabled performance timing | Report synchronized max-rank and mean-rank totals for the existing performance regions, including outputs. | Reuse one timing mechanism instead of maintaining a second output-only path. | Serial and two-rank MPI regressions passed. |

## Python API Compatibility

The final public binary conversion module is `vis/python/bin_convert.py`.
`bin_convert_new.py` is not a supported final API.

The canonical module must retain or deliberately migrate all in-repository uses,
including `vis/python/make_athdf.py` consumers of `read_binary`, `write_athdf`, and
`write_xdmf_for`.

| Capability | Compatibility Decision | Reason | Status |
| --- | --- | --- | --- |
| `read_binary` and existing assembly entry points | Preserve base public behavior while adding partitioned/node discovery. Require explicit node inventory metadata for new node shards while preserving inventory-free historical rank fixtures. Keep `mb_data[variable]` as a NumPy array with shape `[n_mbs, nx3, nx2, nx1]`, matching the canonical module contract. | Avoid breaking current analysis scripts while making new node-sharded output complete and auditable. Preallocated arrays avoid an unbounded Python row-list peak and retain ordinary indexing semantics. | Frozen shared/rank fixtures, strict canonical shard directories, mandatory node metadata checks, duplicate-preheader rejection, explicit empty-shard assembly, duplicate-name and complete grid/logical/exact-logical-geometry rejection with zero-relative storage-aware bounded absolute tolerance, uniform emitted MeshBlock extents within each file, strict coarsening factors, bounded cumulative decoded-text metadata before decode and parameter reads, fixed-grammar scalar and preheader parsing without unconstrained token lists, retained MeshBlock metadata, transient duplicate-ownership accounting, direct preallocated shard reconstruction without Python row-list amplification, explicit ndarray contract assertions, internal-only accounting fields, omitted/mismatched ghost rejection including singleton axes, interior-width ghost placement with extended coordinates, intersecting-cell crop selection, deterministic uncovered levels, restored restriction, cumulative athdf-like allocation and NumPy coordinate/prolongation/restriction temporary preflight, and conversion tested. |
| `read_single_rank_binary_as_athdf` | Preserve the base signature unless a separately documented additive entry point is needed. | Gotham changed this API incompatibly. | Preserved in canonical module; destination-array and dtype contract plus direct signature audit complete; ghost-bearing or sliced emitted extents are rejected because the preserved signature has no ghost argument. |
| `write_athdf` | Retain in canonical `bin_convert.py`. | Base `make_athdf.py` and CLI require it; Gotham variants omit the definition incorrectly. | Preserved and conversion smoke-tested. |
| `write_xdmf_for`, `convert_file`, and `make_athdf.py` | Retain supported conversion workflow and CLI; strip only the final `.bin`/`.cbin` extension when naming outputs. Route the batch wrapper through canonical `convert_file`, recognize both conversion suffixes while excluding the separate `.sph.bin` family, expose shard assembly plus reader-budget flags, and preserve direct `main(file_stem=..., verbose=...)` calls by defaulting additive options. | Keep established ATHDF/XDMF user paths usable without corrupting identifiers containing `bin`, misrouting spherical-slice products, duplicating conversion logic, narrowing the canonical format surface, or breaking Python callers of the historical wrapper. | Preserved; package import, API, CLI, direct legacy-wrapper invocation, `.bin`/`.cbin` batch conversion, `.sph.bin` exclusion, shard-option exposure, and budget-propagation tests passed after consolidating output naming and conversion. |
| `read_pdf.py` | New public reader with dense/sparse node/rank reconstruction and strict validation. | New output format needs supported readback. | Legacy, transitional unversioned shared/rank sparse, V2 shared/rank/node, explicit rejection of unversioned node shards, sibling V2-declaration agreement, intrinsic-payload header-declaration enforcement, required declared V2 preambles and matching cycles, mandatory V2 sparse inventories, node-header `payload_rank` binding to the binary preamble, duplicate node payload-rank rejection, duplicate scalar and dimension-alias rejection, finite metadata and histogram payloads, bounded cumulative decoded-text plus per-record object metadata and payload reads, explicit/generated/legacy edge tokenization and NumPy preflight, pre-materialization and pre-decode strict ASCII-only legacy numeric-row admission with upfront dense preallocation, replacement-header retained-state accounting including legacy candidates, duplicate-validation peaks, per-sibling local-array release, inventory-summary retention, aggregate metadata normalization, and malformed shard inputs tested. |
| `read_sphslice.py` | New public reader with complete ownership and metadata checks plus narrow historical version-1 shared/rank fallback. | New output format needs supported readback without dropping already-readable historical layouts. | Shared/rank/node, strict explicit inventories for new rank/node layouts, synthetic historical `single_file_per_rank` inference regression for shared/rank layouts, duplicate scalar and variable-metadata rejection, positive radius, at least two angles per axis, bounded point count, mandatory explicit modern layout, bounded whole-file reads plus cumulative decoded-text and per-record object metadata, individual metadata lines, pre-split variable tokens, retained reference variable-metadata summary through coordinate generation, retained arrays, coordinate arrays, payload-copy and duplicate-validation peaks, embedded-header offsets including strict header-only API, incremental sparse combination with per-sibling array release, aggregate metadata normalization, and malformed shard inputs tested. |
| `bin_convert_new.py` | Do not add as public module or compatibility alias. | The branch removes redundant conversion APIs. | Accepted. |

## Diagnostic Semantics Boundaries

| Diagnostic surface | Supported contract | Evidence |
| --- | --- | --- |
| Generic PDF diagnostics `coord_*`, `mdot_*`, `edot_*`, `vel_*` | Available for single-fluid Hydro/MHD as documented. | Modern N-D PDF writer/readback tests pass. |
| Generic `mdot_*`, `edot_*`, `vel_*` in `<ion-neutral>` runs | Rejected until module-qualified two-fluid semantics exist. | Negative two-fluid PDF construction regression passes. |
| `edot_sph_out`, `edot_sph_in`, `edot_vert_out`, `edot_vert_in` | Positive and negative partitions of signed total energy flux. MHD partitions include the Poynting contribution before classification, so their sign can differ from gas motion. | Adversarial radial and vertical MHD regressions pass. |
| PDF `weight = mass` and `weight = variable` | Mass weighting requires positive finite conserved density times cell volume. Variable weighting accepts finite, possibly signed values times cell volume. | Automated zero, negative, NaN, and infinity mass-density rejection rows pass. Automated finite signed-variable acceptance and NaN/infinity variable-weight rejection rows pass. |
| `sphslice` variables | Native state-backed scalar fields and native multi-field groups are supported; derived arrays are rejected until ghost-zone-safe interpolation exists. The origin-centered spherical surface must fit inside a 3D domain: `slice_r` must be positive and strictly interior to every domain face. Writers reject non-finite interpolated values and reject finite values that become non-finite after narrowing to the serialized float payload. Readers reject non-finite payloads and validate intrinsic header dimensions, radius, and point count, but cannot independently prove the writer-side domain-interior condition from the lightweight file header alone. | Native-group readback, non-finite producer and forged-payload rejection, finite-double-to-float narrowing rejection, mixed-level AMR ghost-snapshot oracle, and negative construction regressions pass. |

## New Binary Payload Portability Boundary

Modern PDF V2 and spherical-slice payload scalars currently use host-native
byte order. Reader/writer tests cover same-platform round trips; this branch
does not claim cross-endian portability for those new binary payloads. A
portable fixed-byte-order migration must introduce an explicit versioned
format decision and reader compatibility tests rather than silently changing
existing bytes.

## Deferred Documentation Evidence

Candidate documentation is staged, not applied to the live Pages branch, at
`deferred_docs/gh-pages/io-output-formats-and-sharding/`. The bundle records
the inspected `origin/gh-pages` baseline, an exact nine-file allowlist,
protected and target Git blobs, bounded section anchors, source hashes,
contradiction checks, and the post-merge application procedure. Apply the
bundle only through `scripts/stage_gh_pages_io_docs.py` in a detached Pages
worktree. The helper fails closed on drift and can emit an external,
write-free reconciliation packet for deliberate review.

Native direct node restart loading has replaced the earlier transient rank-0
`<manifest>.assembled` staging design. Do not publish or open a Pages
integration review until the code branch is merged and the staged pages are
reapplied, rebuilt, and re-audited against the then-current `origin/gh-pages`.

A detached temporary worktree at `origin/gh-pages`
`4833aa9341e19861297e330ff02aabfd8001935c` was used to apply the focused
bounded transformations, verify the exact allowlist, and run the warnings-as-
errors HTML and link-check builds. The command
`make clean html SPHINXOPTS="-W --keep-going"` completed successfully. No
Pages publication or live branch edit is part of this code feature branch.
