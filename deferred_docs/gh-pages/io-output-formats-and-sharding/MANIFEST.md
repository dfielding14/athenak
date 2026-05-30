# Deferred GitHub Pages Integration Manifest: IO Output Formats And Sharding

## Status And Publication Rule

This directory is a deferred documentation integration bundle for the
`gh-pages` branch. It describes code implemented on the IO feature branch and
must not be published before that code is merged.

CP-03 native direct node restart loading has replaced the earlier temporary
rank-0 `<manifest>.assembled` staging design. Do not publish this bundle until
the IO code branch is merged and the staged pages are reconciled, rebuilt, and
re-audited against the then-current Pages baseline.

The RCP-07 preview baseline is:

```text
origin/gh-pages 4833aa9341e19861297e330ff02aabfd8001935c
```

The branch carrying this bundle intentionally does not copy the live
`docs/source/` tree into the code branch. Apply it only with
`scripts/stage_gh_pages_io_docs.py` in a detached Pages worktree. The helper
never fetches, commits, pushes, switches branches, or edits a checked-out
`gh-pages` branch.

## Machine-Readable Contract

[`manifest.json`](manifest.json) is the authoritative staging contract. It
records:

- the expected `origin/gh-pages` baseline;
- the exact nine-file public allowlist;
- protected live blobs that must remain unchanged;
- baseline blobs for existing targets and expected absence for the new page;
- source-fragment SHA256 values;
- unique anchors and baseline section SHA256 values for bounded edits; and
- contradiction-search rules.

The helper verifies the complete contract before writing any target file,
builds every transformed document in memory, proves a second pure
transformation is idempotent, and then uses atomic replacement writes.

## Code Surface Documented

| Feature | Implemented source or tool | Executable evidence |
| --- | --- | --- |
| Shared, per-rank, and per-node binary output | `src/outputs/binary.cpp`, `src/file_sharding.hpp` | `tst/test_suite/io/test_node_sharding_mpicpu.py` |
| Uniform 3D active-zone full-volume coarsened binary output | `src/outputs/coarsened_binary.cpp` | `tst/test_suite/io/test_node_sharding_mpicpu.py` |
| N-dimensional PDF output and V2 format | `src/outputs/pdf.cpp`, `vis/python/read_pdf.py` | `tst/test_suite/io/test_output_formats_cpu.py`, `test_output_formats_mpicpu.py` |
| Spherical slices | `src/outputs/spherical_slice.cpp`, `vis/python/read_sphslice.py` | `tst/test_suite/io/test_output_formats_cpu.py`, `test_output_formats_mpicpu.py` |
| Per-node restart manifests and native direct reads | `src/outputs/restart.cpp`, `src/restart_manifest.cpp`, `src/main.cpp`, `src/pgen/pgen.cpp`, `src/outputs/io_wrapper.cpp` | `tst/test_suite/io/test_node_sharding_mpicpu.py`, `test_chunked_io_mpicpu.py` |
| Output timing and final-output policy | `src/driver/driver.cpp` | `tst/test_suite/io/test_io_finalization_timing_cpu.py`, `test_io_finalization_timing_mpicpu.py` |
| Single supported binary converter | `vis/python/bin_convert.py` | `tst/test_suite/io/test_python_io_readers_cpu.py` |
| Runnable examples | `inputs/io/*.athinput`, `vis/python/examples/read_io_outputs.py` | `tst/test_suite/io/test_io_examples_cpu.py`, `test_node_sharding_mpicpu.py` |

## Exact Public Allowlist

| Target path | Operation |
| --- | --- |
| `docs/source/configuration.md` | Replace the unique stale `### Output Blocks` section with a marker-wrapped reviewed fragment. |
| `docs/source/examples/index.md` | Replace the whole examples index only after its baseline blob matches. |
| `docs/source/examples/io_outputs_and_sharding.md` | Add the worked example only when the target is absent. |
| `docs/source/modules/index.md` | Replace the unique Support Systems section with a marker-wrapped reviewed fragment that updates the Outputs count while preserving every neighboring row. |
| `docs/source/modules/outputs.md` | Replace the unique central output-module body and stale PDF implementation-entry row with marker-wrapped reviewed fragments. |
| `docs/source/running.md` | Replace the unique stale `## Output Files` section with a marker-wrapped reviewed fragment. |
| `docs/source/tools/visualization.md` | Replace two unique IO-owned sections with marker-wrapped reviewed fragments. |
| `docs/source/reference/input_parameters.md` | Replace the unique stale `## Output Blocks` section through the line before `## Public Source Terms`. |
| `docs/source/reference/file_reference.md` | Insert a marker-wrapped IO artifact reference before the unique `## Shipped Input Families` anchor. |

The former complete overlays for `configuration.md`, `running.md`,
`tools/visualization.md`, and `modules/outputs.md` were removed. Their bounded
fragments preserve unrelated live Pages content.

## RCP-07 Design Decisions

| Topic | Decision |
| --- | --- |
| Script language | Use Python 3 standard library only. |
| Anchor strategy | Require exact unique heading anchors and baseline section SHA256 values for replacements. Require an exact unique anchor for the insertion. |
| Idempotence | Wrap every bounded edit in an RCP-07 marker pair and require a pure second transform to reproduce identical bytes. |
| Baseline handling | In strict mode, require detached target `HEAD` and local `origin/gh-pages` to equal the recorded baseline. Verify protected and target Git blobs before staging. |
| Drift handling | Strict mode rejects drift. `--reviewed-drift <packet-dir>` writes a three-way reconciliation packet outside the target without target writes. |
| Build execution | Print strict HTML and link-check commands by default. Run them only with explicit `--run-builds`. |
| Unrelated-page proof | Validate the planned file set before writes and require the final Git status path set to equal the exact nine-file allowlist. |

## Navigation Decision

`docs/source/modules/index.md` already routes to `modules/outputs.md`; the
bounded Support Systems replacement updates its registered-format count from
12 to 13 while keeping the Markdown table continuous. Keeping marker comments
outside the table rows is required because comments inserted between rows split
the rendered table even when Sphinx and linkcheck succeed.
`docs/source/index.md` already routes to the Examples and Visualization
sections on the inspected Pages baseline and remains protected. The
examples-index replacement registers `examples/io_outputs_and_sharding`.

## Scope Boundaries

- Existing `file_type = sph` output remains supported and distinct from
  `file_type = sphslice`.
- `sphslice` samples an origin-centered spherical surface on a 3D domain. Its
  positive radius must be strictly interior to every domain face. It
  accepts native state-backed scalar fields and native multi-field groups, but
  not derived arrays.
- Legacy unsharded PDF syntax remains supported; modern PDF files use V2.
- Modern PDF V2 and spherical-slice binary payload scalars currently use
  host-native byte order; cross-endian portability is not claimed.
- Node-sharded uniform 3D active-zone full-volume `cbin` is covered.
  Every emitted axis extent must be divisible by `coarsen_factor`.
  Lower-dimensional, ghost-zone-expanded, static-refinement, AMR, and sliced
  `cbin` producer workflows are not advertised.
- Multi-rank tests exercised per-node paths on one physical node. Real
  multi-node MPI qualification remains required before production rollout.
- This local workflow does not claim an external CUDA-capable qualification
  run.

## Application And Validation

Follow [VALIDATION.md](VALIDATION.md). Preserve this deferred bundle until the
separate Pages change merges.
