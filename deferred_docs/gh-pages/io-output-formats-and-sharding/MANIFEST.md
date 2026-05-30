# Deferred GitHub Pages Integration Manifest: IO Output Formats And Sharding

## Status And Publication Rule

This directory is a deferred documentation integration bundle for the
`gh-pages` branch. It describes code implemented on the IO feature branch and
must not be published before that code is merged.

CP-03 native direct node restart loading has replaced the earlier temporary
rank-0 `<manifest>.assembled` staging design. Do not publish this bundle until
the IO code branch is merged and the staged pages are reconciled, rebuilt, and
re-audited against the then-current Pages baseline.

Target documentation baseline inspected while preparing this bundle:

```text
origin/gh-pages 4833aa9341e19861297e330ff02aabfd8001935c
```

The branch carrying this bundle intentionally does not copy the live
`docs/source/` tree into the code branch. Apply and review this material in a
temporary worktree based on the then-current `origin/gh-pages` after merging
the IO feature branch.

## Code Surface Documented

The documentation in this bundle is tied to these implemented public surfaces:

| Feature | Implemented source or tool | Executable evidence |
| --- | --- | --- |
| Shared, per-rank, and per-node binary output | `src/outputs/binary.cpp`, `src/file_sharding.hpp` | `tst/test_suite/io/test_node_sharding_mpicpu.py` |
| Shared, per-rank, and per-node coarsened binary output | `src/outputs/coarsened_binary.cpp` | `tst/test_suite/io/test_node_sharding_mpicpu.py` |
| N-dimensional PDF output and V2 format | `src/outputs/pdf.cpp`, `vis/python/read_pdf.py` | `tst/test_suite/io/test_output_formats_cpu.py`, `test_output_formats_mpicpu.py` |
| Spherical slices | `src/outputs/spherical_slice.cpp`, `vis/python/read_sphslice.py` | `tst/test_suite/io/test_output_formats_cpu.py`, `test_output_formats_mpicpu.py` |
| Per-node restart manifests and native direct reads | `src/outputs/restart.cpp`, `src/restart_manifest.cpp`, `src/main.cpp`, `src/pgen/pgen.cpp`, `src/outputs/io_wrapper.cpp` | `tst/test_suite/io/test_node_sharding_mpicpu.py`, `test_chunked_io_mpicpu.py` |
| Hardened binary/coarsened shard inventory and atomic spherical-slice publication | `src/outputs/binary.cpp`, `src/outputs/coarsened_binary.cpp`, `src/outputs/spherical_slice.cpp`, `vis/python/bin_convert.py`, `vis/python/read_sphslice.py` | `tst/test_suite/io/test_writer_hardening_cpu.py`, `test_writer_hardening_mpicpu.py` |
| Output timing and final-output policy | `src/driver/driver.cpp` | `tst/test_suite/io/test_io_finalization_timing_cpu.py`, `test_io_finalization_timing_mpicpu.py` |
| Single supported binary converter | `vis/python/bin_convert.py` | `tst/test_suite/io/test_python_io_readers_cpu.py` |
| Runnable examples | `inputs/io/*.athinput`, `vis/python/examples/read_io_outputs.py` | `tst/test_suite/io/test_io_examples_cpu.py`, `test_node_sharding_mpicpu.py` |

## Files To Apply To `gh-pages`

The `overlay/` files are complete page bodies to copy into their listed
destinations after reconciling any intervening Pages changes.

| Staged file | Target path | Operation | Reason |
| --- | --- | --- | --- |
| `overlay/docs/source/modules/outputs.md` | `docs/source/modules/outputs.md` | Replace after review | Correct PDF V2, slice, sharding, restart, timing, and policy guidance. |
| `overlay/docs/source/tools/visualization.md` | `docs/source/tools/visualization.md` | Replace after review | Establish `bin_convert.py` as the sole supported binary converter and add new readers. |
| `overlay/docs/source/configuration.md` | `docs/source/configuration.md` | Replace after review | Expose final parser names and safe output configuration patterns. |
| `overlay/docs/source/running.md` | `docs/source/running.md` | Replace after review | Document run, readback, timing, and per-node restart invocation. |
| `overlay/docs/source/examples/index.md` | `docs/source/examples/index.md` | Replace after review | Register the new worked IO example page. |
| `overlay/docs/source/examples/io_outputs_and_sharding.md` | `docs/source/examples/io_outputs_and_sharding.md` | Add | Provide executable input/readback workflows. |

The `insertions/` files are reviewed MyST sections for existing broad
catalogue pages. Merge them into the live files at the identified subject
headings; do not replace the unrelated reference catalogue.

| Staged file | Target path | Merge location | Reason |
| --- | --- | --- | --- |
| `insertions/docs/source/reference/input_parameters.io_outputs.md` | `docs/source/reference/input_parameters.md` | Output and time-parameter reference sections | Record finalized keys, defaults, accepted values, and failures. |
| `insertions/docs/source/reference/file_reference.io_outputs.md` | `docs/source/reference/file_reference.md` | Output-file format reference section | Record file families, shard paths, manifest/payload rules, and readers. |

## Navigation Decision

`docs/source/modules/index.md` already routes to `modules/outputs.md`, and
`docs/source/index.md` already routes to the Examples and Visualization
sections on the inspected Pages baseline. No edits to those two files are
required for this feature. The only new navigation target is
`examples/io_outputs_and_sharding`, registered by the staged replacement for
`docs/source/examples/index.md`.

If the live Pages routing changes before application, re-check this decision
against the then-current toctrees rather than copying it mechanically.

## Stale Guidance To Remove

When applying this package, remove or replace guidance that:

- promotes `bin_convert_new.py` or describes `bin_convert.py` as legacy;
- describes modern PDF payloads as unversioned doubles or uses `logscaleN`
  instead of `scaleN`;
- states that `mass_weighted` is not accepted for legacy PDF inputs;
- omits `single_file_per_node`, per-node restart manifests, `sphslice`,
  `output_timing`, or `final_output_policy`;
- implies sliced `cbin` output is qualified by this feature.

## Scope Boundaries

- Existing `file_type=sph` output remains supported and distinct from the new
  `file_type=sphslice`.
- `sphslice` currently accepts native state-backed scalar fields only;
  derived-array fields are rejected until ghost-zone-safe interpolation is
  implemented.
- Legacy unsharded PDF input syntax remains supported and is tested against
  frozen `origin/main` one- and two-dimensional output bytes.
- Modern PDF files use the versioned V2 representation documented here.
- Generic `mdot_*`, `edot_*`, and `vel_*` diagnostics reject
  `<ion-neutral>` two-fluid inputs until module-qualified semantics are
  introduced.
- Full-volume node-sharded `cbin` is covered; sliced `cbin` is not promoted and
  is rejected explicitly because its emitted extent is incompatible with
  supported coarsening. The full-volume producer validates `coarsen_factor` as
  a power of two between `2` and the shortest MeshBlock dimension before
  writer construction, then requires every emitted extent, including optional
  ghost zones, to be divisible by the factor during construction.
- Node-sharded binary and full-volume coarsened-binary files add inventory
  metadata without changing legacy shared or rank files. Readers accept
  explicit empty shards, reject incomplete or duplicate node inventories, and
  bound cumulative metadata, fixed MeshBlock metadata, payload, and
  incremental aggregate reconstruction sizes, including transient copies,
  while rejecting duplicate variables, malformed grid/logical metadata,
  nonuniform emitted MeshBlock extents within one file,
  geometry outside each MeshBlock's exact logical physical interval using zero
  relative tolerance and a storage-aware absolute tolerance capped at one
  eighth of the logical block width, and invalid
  variable/coarsening counts. Athdf-like helpers preflight cumulative conversion
  allocations and NumPy coordinate/prolongation/restriction-generation
  temporaries, reject omitted or mismatched ghost counts including impossible
  singleton-axis widths, place requested ghost zones by interior MeshBlock width
  with extended coordinates, and retain cells
  intersecting a selected lower bound, initialize uncovered partial-shard
  levels deterministically, crop prolongation before materialization, and
  restore bounded fine-to-coarse restriction. The preserved legacy
  single-MeshBlock ATHDF helper has no ghost argument and rejects ghost-bearing
  or sliced emitted extents instead of silently truncating them.
- Spherical-slice shards declare dense or sparse-angular layout metadata,
  preserve explicit empty shards, and publish atomically through a temporary
  file followed by rename. Readers combine sparse siblings incrementally,
  bound whole-file reads, cumulative header bytes, individual metadata lines,
  variable-token expansion before splitting, retained reference variable-metadata summary
  through final coordinate generation, cumulative retained arrays, embedded input dumps,
  payload-copy, duplicate-validation, and ownership-diagnostic peaks, release
  incorporated sparse arrays before reading each sibling, bound coordinate
  allocations and embedded-header offsets, and normalize reconstructed
  aggregate metadata.
- Modern PDF sparse shards declare rank/node sibling inventory metadata and
  publish each checked header or payload file atomically through a temporary
  file followed by rename. The reader rejects malformed shard aliases,
  identifier mismatches, inventory gaps, and unreasonable metadata or payload
  sizes. Sibling headers agree on V2 declaration, declared V2 headers require
  V2 payload preambles with matching cycles and finite metadata, require
  writer-mandatory inventory fields for V2 sparse families, reject non-shared
  modern dense distribution metadata, account for retained reference state
  while parsing each replacement shard header including legacy candidates,
  preflight explicit/generated and legacy bin-edge tokenization and NumPy
  materialization, require ASCII-only legacy numeric rows and parse them
  strictly with cumulative retained and final stacking bounds,
  bound payload-copy and duplicate-validation peaks before materialization,
  release incorporated local headers and sparse arrays before reading each
  sibling, and remove shard-local identifiers from reconstructed aggregates. Readers retain
  historical compatibility for transitional unversioned dense/sparse binary
  payloads that may omit additive inventory metadata.
- Node restart loading accepts the public manifest only. It validates generated
  relative paths, symlink containment, ordered inventory, exact segment
  coverage, payload sizes, replicated headers, bounded payload and segment
  inventories, positive segment counts, and an on-disk node-payload marker
  before reading rank-local MeshBlock spans directly through chunked positioned
  reads. Marked hard-link and copied-payload aliases are rejected outside
  manifest loading. Production resume does not create `.assembled` files.
- Multi-rank tests exercised per-node paths on one physical node. A real
  multi-node MPI qualification run remains required before production rollout.

## Application And Validation

Follow [VALIDATION.md](VALIDATION.md) after the code merge. Record:

1. the final code branch merge commit;
2. the `gh-pages` baseline actually used;
3. every applied replacement/insertion and any conflict resolution;
4. stale guidance removed;
5. warnings-as-errors Sphinx build output;
6. the independent code-to-doc audit disposition.
