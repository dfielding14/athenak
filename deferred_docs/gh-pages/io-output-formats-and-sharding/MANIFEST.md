# Deferred GitHub Pages Integration Manifest: IO Output Formats And Sharding

## Status And Publication Rule

This directory is a deferred documentation integration bundle for the
`gh-pages` branch. It describes code implemented on the IO feature branch and
must not be published before that code is merged.

Do not publish this bundle until CP-03 replaces temporary rank-0
`<manifest>.assembled` restart staging with direct distributed payload reads.
The staged user-facing pages intentionally avoid claiming that direct loading
has landed.

Target documentation baseline inspected while preparing this bundle:

```text
origin/gh-pages 8eb329959244d68bffc3b2347432a9c85a622394
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
| Per-node restart manifests | `src/outputs/restart.cpp`, `src/main.cpp` | `tst/test_suite/io/test_node_sharding_mpicpu.py` |
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
- implies sliced `cbin` node-sharded output is qualified by this feature.

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
- Full-volume node-sharded `cbin` is covered; sliced node-sharded `cbin` is not
  promoted because existing shared sliced `cbin` readback has an unresolved
  zero-width meshblock-extent defect.
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
