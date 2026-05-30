# IO Feature Audit Ledger

## Branch Identity And Rules

| Item | Value |
| --- | --- |
| Feature branch | `feature/io-output-formats-and-sharding` |
| Clean base | `origin/main` at `886dd2a1437e45a3a30b3eeebf2adfa838328f73` |
| Historical evidence | `origin/gotham-1.0` and explicitly listed source commits |
| Forbidden implementation source | `origin/feature/single-file-per-node-outputs` |

This ledger records independent audits, blocking findings, implementation
resolutions, and re-audit results. `IO_FEATURE_BRANCH_GUIDE.md` is the
historical feature-isolation plan. Active execution is governed by the frozen
`IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md`.

## Audit Gates

| Gate | Scope | Status | Evidence Or Required Resolution |
| --- | --- | --- | --- |
| Historical-scope separation | Gotham IO changes versus unrelated source/physics/mesh/test drift | Complete | Independent audit accepted; only listed generic IO behaviors were reconstructed. |
| File-format compatibility | Headers, versions, names, frozen fixtures and reader boundaries | Corrective local matrix passed; fresh re-audit pending | Frozen compatibility, modern writer/reader, node-PDF payload binding, and one-node node-layout evidence pass; real multi-node qualification remains recorded. |
| Python API consolidation | One canonical `bin_convert.py`, consumers, CLI and reader behavior | Locally closed; committed-tree rerun active | Canonical module/readers/examples pass focused tests after fixed-grammar parsing, direct preallocated shard assembly, and legacy wrapper-default correction; external qualification remains recorded separately. |
| MPI and restart correctness | Node sharding, empty shards, manifests and resume numbering | Locally complete; external topology qualification remains | Native direct manifest loading, forced chunks, alias rejection, and one-node MPI tests pass. True multi-node routing with an empty or non-owning node remains required. |
| Tests and examples | Harness placement, fixtures, negative cases and executable usage | Corrective local matrix passed; final audit pending | Corrected serial/MPI suites, style gate, fixture checksums, and executable examples pass; GPU and multi-node execution remain external gates. |
| Deferred Pages integration | Candidate docs content, final API alignment and later Sphinx validation | Locally accepted; publication deferred | The exact nine-file detached stage, strict builds, rendered table checks, and fresh independent structure audit pass. Live Pages remains unchanged and publication waits for code merge. |

## Evidence Log

### 2026-05-29: Clean Base And Branch Creation

| Field | Record |
| --- | --- |
| Snapshot reviewed | `origin/main` / initial branch HEAD `886dd2a1437e45a3a30b3eeebf2adfa838328f73` |
| Commands | `git status --short --branch`; `git rev-parse HEAD origin/main origin/gotham-1.0`; initial branch creation from the clean base; later `git branch -m feature/io-output-formats-and-sharding` |
| Finding | The detached starting worktree was exactly `origin/main`; only the guide was untracked. |
| Resolution | Created the feature branch directly at the clean base and later renamed it to `feature/io-output-formats-and-sharding` at the user's request. |
| Status | Complete. |

### 2026-05-29: Historical Diff Triage

| Field | Record |
| --- | --- |
| Evidence reviewed | Listed IO-related commits from `6774d871` through `4b01b047`; scoped diffs in `src/outputs`, `src/driver`, `src/main.cpp`, `src/mesh`, and `vis/python`. |
| Commands | `git show --stat` over listed commits; `git diff --name-status origin/main..origin/gotham-1.0 -- ...`; scoped `git show` and `git diff` inspection. |
| Finding | The final Gotham branch includes large unrelated mesh, physics, test deletion, and tooling drift. N-D PDF and spherical-slice feature patches are comparatively scoped; node restart/sharding requires a dedicated MPI correctness audit. |
| Resolution | Reconstruct only identified IO behaviors and explicitly exclude coupled files. |
| Status | In progress pending independent scope audit report. |

### 2026-05-29: Timing And Final-Output Runtime Policy

| Field | Record |
| --- | --- |
| Evidence reviewed | Gotham `59f4b3fe`, `e4d5fc2c`, `b313c5a1`, and the correction visible in `4b01b047`; current `src/driver/driver.cpp` baseline. |
| Finding | Gotham timing was unconditional and reported rank-0 time only. Gotham final counter suppression risks overwriting a terminal restart after resume. |
| Implemented resolution | Added opt-in `<time>/output_timing`; MPI maximum-time reporting; explicit `<time>/final_output_policy = all|restart_only|none` with default `all`; final writes use ordinary counter semantics. |
| Files changed | `src/driver/driver.hpp`, `src/driver/driver.cpp` |
| Verification | `cmake -S . -B /tmp/athenak-io-build -DCMAKE_BUILD_TYPE=Debug -DAthena_ENABLE_MPI=OFF && cmake --build /tmp/athenak-io-build -j 4` passed after initializing the required `kokkos` submodule. |
| Status | Implementation, serial runtime regression, and two-rank MPI timing reduction regression complete; independent review required. |

### 2026-05-29: Python Converter And Reader Audit

| Field | Record |
| --- | --- |
| Auditor scope | Read-only subagent audit of `vis/python/bin_convert.py`, Gotham `bin_convert_new.py`, `read_pdf.py`, `read_sphslice.py`, and base `make_athdf.py`; expressly excluded rejected extraction branch. |
| Evidence reviewed | Base `bin_convert.py` and `make_athdf.py`; Gotham Python file history and versions at commits `d8080a40`, `139c094b`, `edbcd4f3`, and `8fae12b1`. |
| Blocking finding | Gotham drops `write_athdf` while both converter paths continue to call it. Promoting either Gotham converter verbatim would break the base ATHDF/XDMF workflow used by `make_athdf.py`. |
| Additional risks | Gotham changes a public single-rank conversion signature; PDF reader under-validates zero-byte/truncated shards and metadata/time consistency; spherical-slice reader does not reject duplicate angle ownership or time/cycle/version mismatch. |
| Resolution direction | Keep one public `vis/python/bin_convert.py`, preserve base conversion writers and public compatibility, add only needed modern partition discovery/reassembly, and port new readers with strict validation. |
| Status | Blocking finding accepted; Python implementation assigned with disjoint file ownership. |

### 2026-05-29: Test, Example, And Deferred Documentation Audit

| Field | Record |
| --- | --- |
| Auditor scope | Read-only audit of current test harness and `origin/gh-pages` routing/content; rejected extraction ref expressly not inspected. |
| Evidence reviewed | `origin/main` test driver/utilities and current output implementations; `origin/gh-pages` at `8eb329959244d68bffc3b2347432a9c85a622394` including Outputs, Visualization, Examples, and routing pages. |
| Blocking findings | Legacy PDF on `origin/main` is text-form 1D/2D and must be frozen before schema changes; existing `file_type=sph` must coexist with new `sphslice`; Pages pages already contradict current code and must be reconciled. |
| Major requirements | Use `inputs/io/` for promoted user decks, `tst/test_suite/io/` for tests, isolated run directories, and a scheduler-backed or injectable multi-node qualification for true node sharding. |
| Resolution direction | Update guide and compatibility contract; produce baseline fixtures from clean main; stage deferred Pages overlay only after implemented interfaces stabilize. |
| Status | Findings accepted; fixture generation and feature implementation in progress. |

### 2026-05-29: C++ Output, MPI, And Restart Historical Audit

| Field | Record |
| --- | --- |
| Auditor scope | Read-only analysis of `origin/main` and listed Gotham IO commits; rejected extraction ref expressly not inspected. |
| Evidence reviewed | Generic derived/PDF series `6774d871`, `82d924d0`, `dcd7c2d1`, `abefea9f`, `d0625a33`, `965cf7ce`, `97fd9e85`; sphslice `d6887323`, `edbcd4f3`, `868c588e`, `95f5bf56`; timing/final policy `59f4b3fe`, `e4d5fc2c`, `b313c5a1`, `4b01b047`. |
| Accepted findings | `abefea9f` fixes are mandatory for PDF/derived variables; later gather/merge sphslice design supersedes initial node aggregation; final counter suppression must not be ported; per-node restart needs the repaired manifest model plus an atomic completion protocol; node setup must be opt-in/lazy; empty coarsened shards need explicit handling. |
| Resolution direction | Core PDF/derived/sphslice implementation assigned separately; node sharding/restart deferred until these blockers are designed and tested; driver policy already implements accepted final/timing design. |
| Status | Historical scope audit complete; MPI/restart implementation gate remains open. |

### 2026-05-29: Frozen Baseline Fixtures And Runtime Regression

| Field | Record |
| --- | --- |
| Fixture source | Detached untouched worktree at `origin/main` / `886dd2a1437e45a3a30b3eeebf2adfa838328f73`. |
| Producer coverage | Serial shared `.bin`, `.cbin`, one- and two-dimensional legacy text-form `.pdf`/`.bins.pdf`, `.rst`; two-rank `single_file_per_rank` `.bin`, `.cbin`, `.rst`. |
| Evidence stored | `tst/fixtures/io/origin_main_886dd2a1/README.md`, producer input decks, artifacts, and `SHA256SUMS`. |
| Runtime policy verification | Serial debug build: `test_io_finalization_timing_cpu.py -q` returned `6 passed` covering default, `restart_only`, `none`, timing opt-in, invalid policy, and terminal-restart resume behavior. Feature-branch MPI debug build: `test_io_finalization_timing_mpicpu.py -q` returned `1 passed`, confirming one rank-maximum timing line per output event. |
| Resolution | Baseline artifacts are frozen before new writer schemas land; runtime policy input parameters are explicitly present in the test deck for command-line overrides. |
| Status | Fixture capture and serial policy regression complete; new reader/readback and MPI assertions pending. |

### 2026-05-29: Core Writer And Python Integration Correction

| Field | Record |
| --- | --- |
| Initial finding | Independent integration review found that the first C++ V2 PDF writer emitted `AKPDFV2` payload preambles and interleaved sparse records while the initial reader understood an unversioned test layout; it also routed pure legacy PDF input through the modern writer, violating the declared compatibility contract. |
| Implemented resolution | Added a byte-compatible legacy text PDF writer path for pure legacy configuration blocks; updated `read_pdf.py` to parse the emitted V2 header/preamble, dense payloads, and interleaved sparse records while retaining frozen legacy readback. |
| Writer-reader verification | Feature builds produced a three-axis mixed `log`/`linear`/`symlog`, mass-weighted V2 PDF and `sphslice`; readers loaded shared files and reconstructed two-rank per-rank output identically to MPI shared output. |
| Legacy verification | Feature build runs from the clean one- and two-dimensional legacy producers compared byte-identically against all frozen `.pdf` and `.bins.pdf` artifacts. |
| Automated verification | Serial IO modules: `32 passed`; MPI IO modules: `2 passed`, including legacy one-/two-dimensional compatibility, three-/four-dimensional writer readback, scalar/volume/mass weights, invalid PDF axes, canonical converter CLI, example readback, writer-reader comparisons, and timing reduction. |
| Status | Core shared/rank PDF, spherical-slice, reader, and legacy compatibility gate passed; node sharding remains pending. |

### 2026-05-29: Node Sharding And Restart Integration

| Field | Record |
| --- | --- |
| Independent implementation review | Bounded worker implemented `FileShardMode`, lazy node communicators, node writers, transient rank-0 `.assembled` restart staging, atomic payload/manifest publication, and strict inventory/path validation; main agent reviewed the returned diff and required path/inventory hardening before integration. |
| Integration findings fixed | Node `.bin`/`.cbin` metadata was initially emitted after `number of variables` despite extending `size of preheader`, breaking canonical readback; moved node metadata into the declared preheader. Existing per-rank resume omitted the mode flag in `Mesh::BuildTreeFromRestart::GetPosition`, invoking MPI-IO on a standard file; fixed and regression-tested. |
| Tests | Two-rank MPI tests cover full and sliced `.bin` node reconstruction, full `.cbin` node reconstruction, PDF and sphslice node/shared equality, canonical node conversion, node restart terminal resume, traversal/absolute-path/incomplete-marker/byte-count manifest rejection, rank restart resume, sharding exclusivity, node timing labels, and promoted node example execution. Additional manifest inventory negatives remain for CP-03. |
| Environment limit | Two local MPI ranks share one physical node; true multi-node file counts and a genuinely empty node shard require scheduler-backed qualification. |
| Status | One-node MPI integration passed; multi-node qualification explicitly remains before production scale claims. |

### 2026-05-29: Deferred Pages Package And Final Local Qualification

| Field | Record |
| --- | --- |
| Documentation bundle | Added `deferred_docs/gh-pages/io-output-formats-and-sharding/` with target manifest, publication/validation procedure, complete focused overlays for Outputs, Visualization, Configuration, Running and Examples, and insertion fragments for the broad input/file reference catalogues. |
| Target baseline | Validated against detached temporary worktree of `origin/gh-pages` at `8eb329959244d68bffc3b2347432a9c85a622394`; no live Pages branch modification or publication performed. |
| Documentation validation | Copied overlays into `/tmp/athenak-gh-pages-io-docs-e948`, included the reference fragments as validation-only `.inc` sections, and ran `make clean html SPHINXOPTS="-W --keep-going"` successfully. |
| Final serial qualification | Rebuilt `/tmp/athenak-io-build/src/athena`; `/Users/dbf75/.uv/envs/interactive/.venv/bin/python -m pytest` over the five serial IO modules returned `38 passed`, including derived-`sphslice` and two-fluid generic-diagnostic rejection regressions. |
| Final MPI qualification | `/Users/dbf75/.uv/envs/interactive/.venv/bin/python -m pytest` over the three MPI IO modules against `/tmp/athenak-io-build-mpi/src/athena` returned `11 passed`. |
| Fixture integrity | `shasum -a 256 -c tst/fixtures/io/origin_main_886dd2a1/SHA256SUMS` passed for every immutable baseline artifact. |
| Remaining qualification | Run the GPU-selectable regression with an actual GPU build and qualify per-node output/restart on more than one physical node before merge-to-production claims. |
| Status | Local implementation and deferred documentation validation complete; final independent findings corrected and re-audited. Superseded by the CP-01 re-audit record below. |

### 2026-05-29: Independent Final Audit Findings And Corrections

| Field | Record |
| --- | --- |
| Documentation/tooling audit findings | Independent audit found invalid staged output variable/cadence examples, an unenabled runtime-policy command, a PDF reader call targeting a header rather than data, an unimplemented promoted `sph` coexistence claim, malformed-PDF tests that did not reach the V2 parser, and overstated binary shard validation language. |
| Documentation/tooling correction | Replaced examples with parser-valid `hydro_w_d` and `dt`/`dcycle`; applied runtime-policy overrides and correct output block label; corrected PDF payload read invocation; added and tested `file_type=sph` in the promoted deck; made malformed sparse tests emit `AKPDFV2`; narrowed binary-assembly claims to implemented checks. |
| Runtime audit findings | Independent audit identified unsafe derived-array interpolation in `sphslice` and ambiguous generic `mdot_*`/`edot_*`/`vel_*` semantics in `<ion-neutral>` two-fluid runs. |
| Runtime correction | `sphslice` now rejects derived-array variables until ghost-zone-safe interpolation is implemented; generic transformed flow/flux diagnostics now reject ion-neutral configurations until module-qualified semantics are defined. Both boundaries are documented and covered by negative executable regressions. |
| Runtime re-audit | The original runtime auditor re-checked both guards and reported no actionable issue in either remediation. |
| Remaining qualification | A real multi-node MPI run and an actual GPU build execution remain required before production/merge-readiness claims; the corrected existing vorticity/current normalization remains an externally visible derived-variable behavior to watch with future analytic coverage. |
| Status | Runtime and documentation/tooling findings resolved; both focused re-audits reported no remaining actionable issues. Superseded by the CP-01 re-audit record below. |

### 2026-05-29: CP-00 Reconstructed Baseline Preservation

| Field | Record |
| --- | --- |
| Snapshot reviewed | `feature/io-output-formats-and-sharding` at clean-base `HEAD=886dd2a1437e45a3a30b3eeebf2adfa838328f73`; comparison ref `origin/feature/single-file-per-node-outputs=47462c5da45d3c763b37fb21505fac5fd3498805`; merge base `886dd2a1437e45a3a30b3eeebf2adfa838328f73`. |
| Inventory | Reviewed tracked modifications under `src/` and `vis/python/bin_convert.py`, plus untracked C++ source, Python readers, examples, pytest modules, frozen fixtures, deferred documentation, and planning records. Generated Python caches were removed before staging. |
| Serial verification | From `/tmp/athenak-io-build/src`: `/Users/dbf75/.uv/envs/interactive/.venv/bin/python -m pytest -q` over `test_io_examples_cpu.py`, `test_io_finalization_timing_cpu.py`, `test_output_formats_cpu.py`, and `test_python_io_readers_cpu.py` returned `37 passed`. |
| MPI verification | From `/tmp/athenak-io-build-mpi/src`: `/Users/dbf75/.uv/envs/interactive/.venv/bin/python -m pytest -q` over `test_io_finalization_timing_mpicpu.py`, `test_output_formats_mpicpu.py`, and `test_node_sharding_mpicpu.py` returned `11 passed`. |
| Fixture integrity | From `tst/fixtures/io/origin_main_886dd2a1`: `shasum -a 256 -c SHA256SUMS` passed for all 27 frozen artifacts. |
| Whitespace handling | Added `.gitattributes` rule `tst/fixtures/io/origin_main_886dd2a1/pdf/** -whitespace` so intentional frozen legacy PDF bytes do not fail `git diff --check`; fixture checksums remain unchanged. |
| Preserved commits | `ca00581b outputs: add reconstructed IO formats and node sharding`; `3049aa92 tools: consolidate IO readers and examples`; `fcc534fe tests: add IO compatibility fixtures and regressions`. |
| Remaining baseline gap | Production node-restart input still reconstructs a temporary `<manifest>.assembled` shared file in `src/main.cpp`; CP-03 must replace staging with direct distributed reads. |
| Status | CP-00 local preservation complete. |

### 2026-05-29: CP-00 Baseline Audit Corrections

| Field | Record |
| --- | --- |
| Independent finding | The first CP-00 auditor observed concurrent preservation commits and identified that twelve frozen `.bin` and `.rst` payloads were present locally but omitted from the test commit because repository-wide ignore rules cover those suffixes. |
| Resolution | Force-added only `tst/fixtures/io/origin_main_886dd2a1/bin/**` and `tst/fixtures/io/origin_main_886dd2a1/rst/**`, then amended the fixture commit to `fcc534fe`. |
| Verification | `git ls-tree -r --name-only HEAD tst/fixtures/io/origin_main_886dd2a1` now includes twelve committed binary/restart payloads; checksum verification and `git diff --check origin/main` pass. |
| Re-audit | Frozen-snapshot auditor confirmed stable `HEAD=fcc534fe`, all 27 checksum-listed fixtures committed and byte-valid, `git diff --check origin/main` clean, 70 committed IO-scoped paths, and no unrelated implementation files. |
| Status | Resolved; CP-00 closed. |

### 2026-05-29: CP-01 Planning And Documentation Consistency Audit

| Field | Record |
| --- | --- |
| Independent findings | Compatibility text overstated the current restart loader as native; guide examples used pre-transactional node restart shard names; the current public restart entry point is the manifest path only; several branch-name references were stale; malformed-manifest coverage was summarized too broadly; `origin/gh-pages` advanced after the earlier deferred-doc validation. |
| Resolution | Corrected branch identity, current-state records, guide shard names, restart-entry-point guidance, compatibility status, and audit-ledger claims. Added a deferred-bundle publication hold until CP-03 removes `.assembled` staging. |
| Decision update | Added D-031 to describe the current staged loader truthfully until CP-03 supersedes it. |
| Remaining work | Reapply and rebuild the deferred Pages bundle against the then-current `origin/gh-pages` tip during CP-07. |
| Re-audit | Independent post-edit consistency review found no blocking or medium documentation defects. Branch identity, staged-loader wording, payload naming, manifest-path-only entry point, narrowed negative-test claims, and the deferred Pages hold are consistent. |
| Status | CP-01 complete. |

### 2026-05-29: CP-02 Preimplementation MPI IO Audit

| Field | Record |
| --- | --- |
| Independent findings | Current `IOWrapper` narrows 64-bit sizes into MPI `int` counts; MPI-file truncation is racy; the reference collective chunk loops correctly schedule zero-byte ranks but can deadlock on rank-local early return; reference offset arithmetic needs stronger checked addition and signed `MPI_Offset` range checks; mesh restart metadata broadcast remains oversized independently of file IO. |
| Resolution direction | Port chunking by intent with overflow-safe math, communicator-agreed or fail-fast collective errors, dummy-buffer participation, one-rank synchronized truncation, and a public chunked `BroadcastBytes()` helper. |
| Explicit exclusions | Do not couple `ParameterInput` to `FileShardMode`; do not change unrelated parameter-header limits. |
| Decision update | Added D-032 through D-036. |
| Status | CP-02 implementation assigned with bounded file ownership. |

### 2026-05-29: CP-02 Chunked MPI IO Implementation Verification

| Field | Record |
| --- | --- |
| Implementation | Added overflow-safe chunked sequential, positioned, and collective MPI byte IO; checked `MPI_Offset` conversion; communicator-agreed collective schedules; dummy-buffer participation for zero-byte ranks; synchronized one-rank MPI-file truncation; and chunked `BroadcastBytes()`. Migrated restart metadata broadcasts without coupling `ParameterInput` to sharding mode. |
| Build verification | Serial `/tmp/athenak-io-build` and MPI `/tmp/athenak-io-build-mpi` builds completed successfully. Existing Clang variable-length-array warnings in `src/mesh/mesh.cpp` remain unchanged. |
| Serial verification | From `/tmp/athenak-io-build/src`, the four focused serial IO modules returned `37 passed`. |
| MPI verification | From `/tmp/athenak-io-build-mpi/src`, the existing MPI IO modules plus `test_chunked_io_mpicpu.py` returned `17 passed`, including shared restart round trip, file truncation reuse, restart resume, and five malformed `ATHENAK_TEST_MAX_MPI_BYTES` values. |
| Fixture integrity | From `tst/fixtures/io/origin_main_886dd2a1`, `shasum -a 256 -c SHA256SUMS` passed for all 27 frozen artifacts. |
| Style checks | `git diff --check origin/main` passes. The bounded worker also reported targeted `cpplint` on owned C++ files and targeted `flake8` on the new pytest passing. |
| Remaining gate | Independent post-integration MPI audit is still required before CP-02 closes and CP-03 source edits begin. |
| Status | Implementation verified locally; awaiting independent audit. |

### 2026-05-29: CP-02 Independent Audit Rejection And Correction

| Field | Record |
| --- | --- |
| Independent blocker | Positioned MPI helpers validated chunk starting offsets but not the inclusive end of the full non-empty byte range. A representable start near `MPI_Offset::max()` could reach MPI even when the range end was not representable. |
| Test-coverage finding | Application-level forced-chunk tests covered shared restart round trip, truncation reuse, migrated broadcasts, malformed override values, and trailing zero-byte collective writes. They did not directly cover zero-byte collective reads, positioned range rejection, multiplication overflow, or communicator-rank disagreement for `ATHENAK_TEST_MAX_MPI_BYTES`. |
| Required correction | Add a shared positioned-range preflight before MPI file operations and a direct MPI wrapper harness for asymmetric zero-byte read/write participation, range overflow, multiplication overflow, and mismatched chunk limits. |
| Decision update | Added D-041. |
| Status | CP-02 stop gate rejected; correction assigned before CP-03 source edits. Superseded by the successful re-audit below. |

### 2026-05-29: CP-02 Correction Re-Audit And Stop-Gate Closure

| Field | Record |
| --- | --- |
| Implemented correction | Added `PreflightPositionedMpiRange()` to validate checked inclusive end offsets before every non-empty positioned MPI read or write. Added a direct MPI wrapper harness compiled against the active MPI build flags. |
| Direct harness verification | `test_io_wrapper_harness_mpicpu.py` returned `5 passed`: asymmetric `7/2/0` byte collective writes and reads, positioned write-range rejection, positioned read-range rejection, `size * count` overflow rejection, and communicator-rank chunk-limit disagreement. |
| Full verification | Corrected serial IO plus GPU-definition smoke returned `38 passed`; corrected MPI IO including the harness and forced chunks returned `22 passed`; all 27 frozen artifacts remained byte-valid; `git diff --check` remained clean. |
| Independent re-audit | The original MPI auditor confirmed the P1 range bug and D-036 coverage gap are resolved, and found no regression in collective symmetry, zero-byte handling, synchronized truncation/open, broadcast chunking, or serial fallback routing. |
| Deferred risks | Caller-side restart arithmetic hardening, serial `fseek`/`ftell` width handling, and forced-small-chunk native node restart remain follow-ups for later checkpoints. |
| Status | CP-02 stop gate closed. |

### 2026-05-29: CP-03 Native Restart Design Audit

| Field | Record |
| --- | --- |
| Independent findings | Production node restart still writes a full rank-0 `<manifest>.assembled` file; manifest parsing is embedded in `src/main.cpp`; reference direct reads provide useful span-routing evidence but its binary manifest, eager communicator setup, ambient inference, and payload-path behavior conflict with the local transactional design. |
| Accepted architecture | Extract a structured text-manifest module; validate once; retain the manifest as the public transactional restart entry point; open a canonical payload header; route MeshBlock reads directly through validated segments; coalesce source spans; use node-collective CP-02 reads with zero-byte participants. |
| Explicit exclusions | Reject payload-path restart, binary manifest replacement, ambient stale-shard discovery, eager node communicator setup, original-rank metadata tables, and unrelated statistics expansion. |
| Decision update | Added D-037. |
| Status | Design accepted; implementation must wait for the CP-02 audit stop gate. |

### 2026-05-29: CP-04 Writer-Hardening Preimplementation Audit

| Field | Record |
| --- | --- |
| Independent findings | Binary and coarsened-binary writers retain legacy large-count branches that should collapse onto the CP-02 wrapper; sliced binary bypasses zero-byte collectives; spherical-slice publication ignores write/close failures and lacks node-leader ownership validation. |
| Accepted architecture | Consume CP-02 collective writes including zero-byte ranks; preserve counter advancement; standardize explicit valid empty node shards; make spherical-slice publication checked and transactional; validate sorted local and world angular ownership before publication; extend strict readers rather than replacing them. |
| Explicit exclusions | Do not skip empty shards, port `ATHENAK_OUTPUT_IO_STATS`, publish spherical-slice files directly to final paths, promote sliced node-sharded `.cbin`, or import unrelated sharding refactors. |
| Decision update | Added D-038 and D-039. |
| Status | Preimplementation direction accepted; execute after CP-03. |

### 2026-05-29: CP-05 Python API Preimplementation Audit

| Field | Record |
| --- | --- |
| Independent findings | The reference `read_rank_binary_as_athdf()` scatter implementation misplaces later logical MeshBlocks at the root-grid origin; its indexed single-block overload changes established positional-call meaning; target-tree readers lack focused athdf-like and malformed-binary coverage. |
| Accepted architecture | Add a thin canonical rank-reader delegate, add only a keyword-only `meshblock_index_in_file=...` selector, omit unneeded `athinput()`, and harden malformed-reader behavior with focused regressions. |
| Decision update | Added D-040, resolving D-017 and D-019. |
| Status | Preimplementation direction accepted; superseded by the successful implementation audit below. |

### 2026-05-29: CP-05 Canonical Python Reader Hardening

| Field | Record |
| --- | --- |
| Implementation | Added `read_rank_binary_as_athdf()` as a thin canonical delegate; added keyword-only `meshblock_index_in_file=...`; retained all existing positional forms; omitted `athinput()`; converted binary readers to context-managed file IO; rejected duplicate logical MeshBlocks across shards; and added explicit empty-shard athdf-like conversion errors. |
| Focused verification | `test_python_io_readers_cpu.py` returned `39 passed`. Python byte-compilation, targeted `flake8`, and `git diff --check` passed. |
| Independent re-audit | A separate Python auditor confirmed the canonical delegate, rank-1 logical placement test, positional compatibility, keyword-only selector, descriptor closure, duplicate-block rejection, empty-shard errors, intended `__all__`, absent `athinput()`, and no unrelated API drift. |
| Deferred optional coverage | Add explicit node-shard delegate readback, parameterize duplicate-block rejection for `.cbin`, and optionally add malformed-header closure and empty-`convert_file()` tests. |
| Status | CP-05 stop gate closed. |

### 2026-05-29: CP-06 Local Qualification Resource Probe

| Field | Record |
| --- | --- |
| MPI launcher helper | Fixed `tst/scripts/utils/athena.py` so MPI execution uses separate `["mpiexec", "-n", ...]` argv entries rather than the invalid single element `"mpiexec -n"`. Python byte-compilation and `git diff --check` passed. |
| Local MPI topology | `mpirun -np 2 hostname | sort -u` reported one physical hostname, `Tin-Drum`. No scheduler launcher command was found locally. |
| Local GPU topology | The host has an Apple M4 Max GPU, but no CUDA or HIP compiler/runtime tool was found. The repository GPU suite configures `Kokkos_ENABLE_CUDA=On`, so this host cannot perform the required GPU qualification. |
| GPU-selectable regression | Strengthened `tst/test_suite/io/test_output_formats_gpu.py` to assert four-dimensional derived PDF axes and scalar weighting metadata. The test logic passed once against the CPU build; that smoke run is not GPU qualification. |
| Remaining external gates | Execute the GPU-selectable regression on a CUDA-capable build and qualify node output/restart on multiple physical nodes, including a genuinely empty or non-owning node. |
| Status | Local resource probe complete; external GPU and multi-node qualification remain required before merge readiness. |

### 2026-05-29: CP-06 Initial Repository Style Gate

| Field | Record |
| --- | --- |
| Command | From `tst/`, `/Users/dbf75/.uv/envs/interactive/.venv/bin/python run_test_suite.py --style`. |
| Result | Failed both C++ and Python style tests. The run generated temporary `tst/test_suite/style/cpplint.py` and Python caches; those generated artifacts were removed after capture. |
| C++ findings | Namespace terminator comments in `src/globals.*`; line wrapping in reconstructed IO files; `long` in `src/outputs/outputs.cpp`; missing `<utility>` include in `src/outputs/pdf.cpp`; and an unapproved `<chrono>` include in `src/outputs/restart.cpp`. |
| Python findings | Line wrapping and slice-spacing issues in reconstructed IO pytest modules, `vis/python/read_pdf.py`, and `vis/python/read_sphslice.py`. |
| Resolution plan | Run one coordinated style sweep after CP-03 and CP-04 stabilize so active-worker files are not edited concurrently. Repeat the full repository style gate before CP-08 closes. |
| Status | Style gate open. |

### 2026-05-29: CP-03 Native Node-Restart Implementation

| Field | Record |
| --- | --- |
| Implementation | Extracted `src/restart_manifest.cpp` and `.hpp`; removed production `.assembled` staging; kept the public text manifest as the only supported node-restart entry point; opened payload 0 for canonical header consumption; routed local MeshBlock data directly through validated segment spans; coalesced adjacent spans; and used node-collective CP-02 reads with zero-byte participants. |
| Validation hardening | Narrowed payload-path rejection to generated `node_XXXXXXXX/<leaf>.g<digits>.payload.rst` artifacts so unrelated shared restarts remain compatible; canonicalized the manifest directory and payload paths; rejected symlink escapes; checked payload sizes, segment coverage, node-local offsets, and close returns; and compared every replicated payload header byte-for-byte against payload 0. |
| Focused verification | Serial and MPI builds passed. From `/tmp/athenak-io-build-mpi/src`, `test_node_sharding_mpicpu.py` returned `32 passed`; `test_chunked_io_mpicpu.py` plus `test_io_wrapper_harness_mpicpu.py` returned `11 passed`. `git diff --check` passed and `rg -n "\\.assembled|StageNodeRestart|CopyFileRange" src` returned no matches. |
| Remaining qualification | Local MPI runs exercise one physical node, forced tiny chunks, zero-byte rank participation, and a two-rank checkpoint resumed on one rank. Real multi-node direct routing with an empty or non-owning node remains required externally. |
| Status | Local implementation checkpoint complete; first independent audit found an alias bypass, corrected and re-tested below. |

### 2026-05-29: CP-04 Writer And Reader Hardening

| Field | Record |
| --- | --- |
| Binary and coarsened binary | Replaced split fallback loops with one CP-02 collective positioned write path including zero-byte ranks; added checked arithmetic, checked metadata and payload writes, checked close returns, additive `number of nodes` metadata, explicit empty-shard publication, and optional sibling-inventory validation in canonical Python readers. Full-volume node `.cbin` remains promoted; sliced node `.cbin` is explicitly rejected. |
| Spherical slice | Validated local and node-merged angular ownership, sorted sparse records with values, checked MPI collectives and file operations, added dense/sparse layout metadata and sibling inventory validation, and published through checked temporary files plus atomic rename. |
| Focused verification | Serial and MPI builds passed. CPU writer-hardening plus canonical-reader matrix returned `51 passed`; `test_writer_hardening_mpicpu.py` returned `2 passed`; the complete local serial IO surface plus CPU smoke of the GPU-selectable regression returned `71 passed`; and the complete local MPI IO surface returned `47 passed`. Fixture checksum verification returned `27` `OK` records. |
| Remaining qualification | Reader-level empty-node fixtures and one-node MPI writers are covered locally. Scheduler-backed multi-node output with a genuinely empty node remains required externally. |
| Status | Local implementation checkpoint complete; first independent audit found three hardening gaps, corrected and re-tested below. |

### 2026-05-29: CP-06 Repository Style Gate Closure

| Field | Record |
| --- | --- |
| Command | From `tst/`, `/Users/dbf75/.uv/envs/interactive/.venv/bin/python run_test_suite.py --style`. |
| Result | Passed: `2 passed in 26.25s`. |
| Resolution | Applied a minimal formatting/include/comment sweep across reconstructed IO source and tests; removed generated `tst/test_suite/style/cpplint.py` and Python `__pycache__` artifacts; and confirmed `git diff --check` has no output. |
| Status | Local style gate closed. |

### 2026-05-29: CP-07 Detached GitHub Pages Validation

| Field | Record |
| --- | --- |
| Pages baseline | Refreshed `origin/gh-pages` and validated against detached commit `4833aa9341e19861297e330ff02aabfd8001935c` in `/private/tmp/athenak-gh-pages-io-docs-e948-final3`. The live `gh-pages` worktree remained clean and untouched. |
| Application | Copied the reviewed overlay page bodies into the detached worktree and appended validation-only MyST include directives for the two insertion fragments beside their live reference targets. |
| Contradiction audit | `rg` over the detached docs found only intended references: canonical `bin_convert.py` guidance, explicit rejection of `bin_convert_new.py`, and explicit statements that production node restart does not create `.assembled` staging files. |
| Build | From the detached `docs/` directory, `make clean html SPHINXOPTS="-W --keep-going"` succeeded. |
| Status | Deferred bundle fits the current Pages framework. Publication remains intentionally deferred until after the IO feature branch is merged. |

### 2026-05-29: CP-03 And CP-04 Post-Audit Corrections

| Field | Record |
| --- | --- |
| Restart-audit blocker | Exact generated payload restart paths were rejected, but lexical `./` aliases reached the same payload as ordinary shared restart input and bypassed the public transactional manifest. |
| Restart correction | Classify generated payload artifacts after lexical normalization and canonical symlink resolution when the path exists. Reject exact, `./`, repeated-separator, symlink-alias, and unpublished `.tmp` payload entry points while preserving unrelated shared `ordinary.payload.rst` compatibility. |
| Restart scaling disposition | Manifest inventory and replicated-header validation currently run on every MPI rank. This is correct but scales as ranks times payloads times header bytes. Central validation plus structured broadcast is recorded as a deferred scaling optimization in D-046 rather than added late without real multi-node qualification. |
| Writer-audit blockers | Sliced node `.cbin` rejection occurred after `LoadOutputData()`, allowing a different one-rank coarsening failure first; Python readers materialized file-controlled expected-ID ranges before checking discovered sibling counts; assembled node `.bin`/`.cbin` objects retained one shard's MeshBlock count. |
| Writer and reader correction | Reject sliced node `.cbin` during writer construction; compare declared and discovered shard counts before bounded expected-ID construction; aggregate assembled `number_of_meshblocks`; and add one-rank, oversized-count, and aggregate-metadata regressions. |
| Focused verification | Serial and MPI builds passed. `test_writer_hardening_cpu.py` returned `17 passed`; `test_writer_hardening_mpicpu.py` returned `3 passed`; and `test_node_sharding_mpicpu.py` returned `36 passed`. |
| Full local verification | Complete serial IO plus CPU smoke of the GPU-selectable test returned `76 passed`; complete MPI IO returned `52 passed`; repository style returned `2 passed`; and `git diff --check` returned no output. |
| Remaining external qualification | CUDA-capable GPU execution and scheduler-backed multi-node output/restart qualification remain required before merge-readiness or production-readiness claims. |
| Status | All locally actionable first-audit findings corrected; fresh final re-audits in progress. |

### 2026-05-29: CP-03 And CP-04 Second-Audit Corrections

| Field | Record |
| --- | --- |
| Restart-audit blockers | Path-based payload rejection still allowed direct restart from a hard link or byte-for-byte copy of a generated node payload. The node writer also did not check replicated-header writes or the payload close return before rename, and the manifest parser accepted an unbounded declared payload count. |
| Restart correction | Add an on-disk node-payload marker after the replicated parameter dump. Require and consume it when loading through a validated manifest; reject it when ordinary restart input opens node-payload bytes directly. Check node-payload header writes and close returns, bound declared payload inventory at `1048576`, and add hard-link, copied-payload, corrupt-marker, and oversized-count regressions. |
| Writer and reader-audit blockers | Modern PDF files were written directly to final names, sparse PDF siblings lacked declared inventory metadata, and PDF, spherical-slice, and binary Python readers accepted file-controlled allocation sizes without practical bounds. Canonical shard-directory spelling and path/header identifiers also needed stricter validation. |
| Writer and reader correction | Publish each modern PDF header and payload through a checked temporary file plus atomic rename; record sparse `rank`/`number_of_ranks` or `node`/`number_of_nodes`; validate complete dense sparse-shard inventories and canonical eight-digit shard directories; reject identifier mismatches; and enforce explicit 512 MiB practical allocation limits in shipped readers. PDF atomicity is per file, not a header-plus-payload family transaction. |
| Focused verification | Serial and MPI builds passed. `test_node_sharding_mpicpu.py` returned `40 passed`; `test_python_io_readers_cpu.py` plus `test_writer_hardening_cpu.py` returned `70 passed`; `test_output_formats_cpu.py` returned `10 passed`; and the focused MPI output/chunk matrix returned `15 passed`. Targeted flake8 and `git diff --check` returned no output. |
| Remaining external qualification | CUDA-capable GPU execution and scheduler-backed multi-node output/restart qualification remain required before merge-readiness or production-readiness claims. |
| Status | All locally actionable second-audit findings corrected; fresh final re-audits and full local verification required before CP-08 closes. |

### 2026-05-29: CP-08 Final-Audit Corrections

| Field | Record |
| --- | --- |
| Test-coverage blockers | Frozen `origin/main` shared and per-rank restart fixtures were checksummed but never resumed, and strict PDF node inventory lacked a locally testable positive synthetic two-shard reconstruction with one explicit empty shard. |
| Python-audit blockers | A header declaring `AKPDFV2` could silently accept a transitional payload without the V2 preamble; V2 payload cycles were not checked against headers; reconstructed PDF and spherical-slice results retained shard-local identifiers; and spherical-slice aggregate `npoints` described shard zero rather than the returned dense surface. |
| Docs-audit blocker | Practical-cap wording exceeded the implementation: selected dense allocations were bounded, but PDF bulk reads, spherical-slice file reads, binary parameter-header reads, multi-MeshBlock accumulation, and shard-aggregate construction still needed preflight limits. |
| Resolution | Require declared V2 preambles and header/payload cycle agreement; normalize reconstructed metadata; bound PDF metadata/payload reads, spherical-slice whole-file reads, binary metadata records and parameter dumps, accumulated binary payloads, and binary shard aggregates; resume frozen restart fixtures; verify frozen fixture manifest/checksums automatically; add positive two-node synthetic PDF inventory reconstruction; and execute the shipped `cbin` summary helper. |
| Focused verification | Python reader/writer hardening returned `83 passed`; serial fixture/examples/policy matrix returned `87 passed`; frozen per-rank restart fixture resume returned `1 passed`; targeted flake8 and `git diff --check` returned no output. |
| Status | Locally actionable final-audit corrections implemented; full local matrix, detached Pages rebuild after final docs refresh, and independent re-audits required before CP-08 closes. |

### 2026-05-29: CP-08 Restart And Aggregate Re-Audit Corrections

| Field | Record |
| --- | --- |
| Restart re-audit blocker | Payload inventory was bounded, but malformed manifests could still append arbitrarily many zero-length segment records before validation, consuming memory and validation time on every rank. |
| Format re-audit blocker | The MPI PDF reconstruction regression still asserted shard-local `rank` and `node` keys after aggregate normalization intentionally removed them. |
| Resolution | Bound segment records before storage by the declared MeshBlock count and a practical `1048576` ceiling; reject non-positive segment counts; add malformed segment-inventory regressions; inspect raw PDF shard headers for local identity while asserting reconstructed aggregates retain total shard counts but no local shard ID. |
| Focused verification | Serial and MPI binaries rebuilt; affected MPI restart/output/chunk/writer matrix returned `53 passed`; targeted flake8 and `git diff --check` returned no output. |
| Status | Locally actionable re-audit corrections implemented; full local matrix, detached Pages rebuild, and independent re-audits remain required before CP-08 closes. |

### 2026-05-29: CP-08 Reader-Preflight And Deferred-Docs Re-Audit Corrections

| Field | Record |
| --- | --- |
| Python re-audit blockers | Sparse PDF sibling headers could disagree on `binary_magic`, and payload parsing reused the initial family header; binary metadata-record bytes were bounded per line but not cumulatively; metadata-only MeshBlock growth and late shard-aggregate checks remained possible; spherical-slice coordinate arrays and input-header read offsets needed preflight. |
| Test and scope re-audit blockers | Frozen fixture automation excluded nested artifacts by basename rather than root-relative path; generated Python and pytest caches repeatedly dirtied the audit workspace. |
| Docs re-audit blockers | The deferred examples index included an unrelated CGM navigation entry; three PDF descriptions overstated the preamble field as a generic shard ID; binary aggregate wording overstated when the check occurred; and Pages publication instructions omitted frozen shared/per-rank restart resume status. |
| Resolution | Compare PDF sibling `binary_magic`, parse each sparse payload against its local header, bound cumulative binary metadata plus fixed MeshBlock and incremental aggregate state, require positive binary variable counts, preflight spherical coordinates and input-header offsets, tighten fixture inventory matching, ignore generated Python caches, remove the CGM navigation entry, correct writer-rank wording, narrow incremental aggregate wording, and add frozen restart resume status to the publication record. |
| Focused verification | Python reader/writer hardening returned `95 passed`; full serial IO plus CPU smoke of the GPU-selectable regression returned `117 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; targeted py_compile, flake8, and `git diff --check` returned no output. |
| Deferred Pages verification | Reapplied the corrected package to detached `/private/tmp/athenak-gh-pages-io-docs-e948-final3` at `origin/gh-pages` commit `4833aa9341e19861297e330ff02aabfd8001935c`; `make clean html SPHINXOPTS="-W --keep-going"` succeeded; live `gh-pages` remained clean. |
| Status | Locally actionable corrections implemented; final independent re-audit dispositions and coherent commits remain required before CP-08 closes. |

### 2026-05-29: CP-08 Retained-Aggregate And Converter Re-Audit Corrections

| Field | Record |
| --- | --- |
| Python re-audit blockers | Spherical-slice reconstruction retained every sparse sibling array before combining; binary serialized metadata accounting undercounted widened retained arrays and aggregate-copy peaks; athdf-like conversion helpers allocated coordinates and dense arrays from parsed grid metadata without preflight; PDF sparse reconstruction retained full sibling headers; spherical header-only readback accepted impossible embedded-header offsets; and coarsened binaries accepted non-positive coarsening factors until division. |
| Docs re-audit blocker | Deferred user pages described V2 writers but omitted intentional historical reader compatibility for transitional unversioned dense/sparse PDF payloads and stated complete inventory requirements too broadly for those additive legacy-compatibility paths. |
| Resolution | Combine spherical shards incrementally, retain only spherical/PDF inventory summaries, account for retained binary metadata and transient aggregate copies, validate coarsening factors and spherical header-only offsets, preflight every athdf-like coordinate/output/level/restriction allocation, and document that writers emit V2 while transitional unversioned reader compatibility remains supported without mandatory additive inventory metadata. |
| Focused verification | Python reader/writer hardening returned `95 passed`; full serial IO plus CPU smoke of the GPU-selectable regression returned `117 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; targeted pycompile, flake8, and `git diff --check` returned no output. |
| Deferred Pages verification | Reapplied the corrected deferred pages into detached `/private/tmp/athenak-gh-pages-io-docs-e948-final3`; `make clean html SPHINXOPTS="-W --keep-going"` succeeded against `origin/gh-pages` `4833aa9341e19861297e330ff02aabfd8001935c`; live `gh-pages` remained clean. |
| Status | Locally actionable corrections implemented; final independent re-audit dispositions and coherent commits remain required before CP-08 closes. |

### 2026-05-29: CP-08 Reconstruction And Producer-Contract Re-Audit Corrections

| Field | Record |
| --- | --- |
| Python re-audit blockers | Athdf-like reconstruction retained placeholder fine-to-coarse branches, full-block prolongation preceded selection slicing, malformed logical levels reached exponent arithmetic, and the single-rank helper replaced caller-provided arrays instead of assigning into them. PDF and spherical readers also enforced several practical bounds per array rather than cumulatively and accepted non-finite PDF metadata. |
| Format re-audit blockers | Passive-scalar group labels wrapped after scalar `99`, allowing distinct fields to collide in serialized metadata; PDF symlog output accepted non-finite thresholds; `cbin` factor validation occurred after unsafe arithmetic; and V2 sparse PDF readers accepted missing writer-mandatory inventory declarations. |
| Resolution | Centralized bounded athdf-like MeshBlock placement with cropped prolongation and restored restriction; validated logical metadata before arithmetic; preserved single-rank destinations and dtype conversion; emitted non-wrapping passive-scalar labels and rejected duplicate binary labels; required finite PDF metadata and mandatory V2 sparse inventory; made PDF/spherical retained-memory checks cumulative; and validated the full `cbin` factor contract before writer construction and reader arithmetic. |
| Focused verification | Python reader/writer hardening returned `115 passed`; full serial IO plus CPU smoke of the GPU-selectable regression returned `146 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Documentation re-audit correction | The first factor wording omitted emitted-extent divisibility, including optional ghost-zone expansion, and the PDF retained-peak wording overstated when sparse shard-local arrays are checked. |
| Resolution | Reject sliced `cbin` consistently in every shard mode; validate full-volume emitted-extent divisibility during writer construction; document the complete factor contract; and state that PDF cumulative retained peaks are checked after each bounded shard load and before aggregate accumulation. |
| Deferred Pages verification | Added the complete `cbin` factor contract to the deferred user-facing insertion and overlays, narrowed PDF retained-peak wording, reapplied the exact files to detached `/private/tmp/athenak-gh-pages-io-docs-e948-final3`, and rebuilt successfully with `make clean html SPHINXOPTS="-W --keep-going"` against `origin/gh-pages` `4833aa9341e19861297e330ff02aabfd8001935c`. |
| Final producer re-audit blocker | Ordered infinite PDF axis bounds passed the range check and could publish malformed metadata with non-finite edges. |
| Resolution | Require finite PDF lower and upper bounds before ordering, log-domain validation, edge construction, or publication; add a producer regression for `output1/bin4_max=inf`. |
| Final verification | Focused producer module returned `19 passed`; full serial IO plus CPU smoke of the GPU-selectable regression returned `146 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; final format-auditor probe confirmed rejection publishes zero PDF files. |
| Status | Locally actionable corrections implemented; final independent re-audit dispositions and coherent commits remain required before CP-08 closes. |

### 2026-05-29: CP-08 Final Python Reconstruction Re-Audit Corrections

| Field | Record |
| --- | --- |
| Python re-audit blockers | Ghost-zone ATHDF placement advanced by emitted rather than interior MeshBlock width and generated coordinates inside root bounds; lower crop bounds dropped intersecting cells; malformed grid, geometry, and oversized logical-location metadata could reach reconstruction; partial-shard level arrays retained uninitialized cells; and PDF, spherical-slice, and ATHDF allocation guards still missed cumulative live peaks. Modern dense PDF headers also accepted arbitrary distribution metadata. |
| Resolution | Validate grid topology, ghost-zone counts, ordered geometry, and bounded logical locations before reconstruction; place ghost-zone blocks by interior width with extended coordinates; centralize intersecting-cell coordinate selection; initialize uncovered level cells to `-1`; preflight retained ATHDF coordinates, outputs, exact-restriction temporaries, PDF centers/payload copies/sparse validation, and spherical-slice dumps/diagnostics; reject non-shared modern dense PDF distributions. |
| Focused verification | Python reader module returned `105 passed`; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Full local verification | Full serial IO plus CPU smoke of the GPU-selectable regression returned `161 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts. |
| Status | Locally actionable Python corrections implemented; final Python re-audit disposition and coherent commits remain required before CP-08 closes. |

### 2026-05-29: CP-08 Final Python Temporary-Budget Re-Audit Corrections

| Field | Record |
| --- | --- |
| Python re-audit blockers | Default `num_ghost=0` could silently accept ghost-bearing products; direct binary reads did not reject duplicate logical MeshBlocks or geometry outside the ghost-extended root domain; prolongation index arrays and center-list construction remained outside ATHDF preflight; PDF symlog vectorization and spherical-slice final coordinate generation retained unbudgeted temporaries. |
| Resolution | Detect emitted ghost zones when the caller omits `num_ghost`; reject duplicate and out-of-domain direct MeshBlocks; preflight prolongation index arrays; generate root centers into preallocated arrays; conservatively preflight symlog transformation peaks and spherical final-coordinate peaks; add focused regressions. |
| Focused verification | Python reader module returned `111 passed`; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Full local verification | Full serial IO plus CPU smoke of the GPU-selectable regression returned `167 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts. |
| Deferred Pages verification | Refreshed the deferred visualization page and manifest with the final reconstruction and live-peak contract, reapplied `tools/visualization.md` to detached `/private/tmp/athenak-gh-pages-io-docs-e948-final3`, and rebuilt successfully with `make clean html SPHINXOPTS="-W --keep-going"` against `origin/gh-pages` `4833aa9341e19861297e330ff02aabfd8001935c`. |
| Status | Locally actionable Python corrections implemented; final Python and documentation re-audit dispositions plus coherent commits remain required before CP-08 closes. |

### 2026-05-29: CP-08 Sparse-Lifecycle, Metadata-Envelope, And API-Coverage Corrections

| Field | Record |
| --- | --- |
| Python re-audit blockers | Prior local PDF shard headers remained live while replacement headers were parsed; binary MeshBlock geometry was checked globally but not against logical-location envelopes; coordinate generation, prolongation indexes, exact-restriction repeats, and sparse duplicate validation retained unbudgeted NumPy temporaries; malformed singleton-axis ghost widths bypassed validation. |
| Test/example re-audit gaps | Three retained ATHDF orchestration APIs lacked positive frozen-fixture automation, and the documented node-sharded `cbin --assemble-shards` helper command was not executed end to end. |
| Resolution | Release incorporated PDF local headers and PDF/spherical sparse arrays before reading each sibling; account for retained PDF reference state during private replacement-header parsing while preserving the public one-argument API; validate exact logical physical intervals and singleton-axis ghost widths; conservatively preflight NumPy source and duplicate-validation arrays; make the preserved single-MeshBlock ATHDF helper reject unsupported ghost-bearing or sliced emitted extents without changing its signature; add weak-reference, reduced-cap, malformed-metadata, frozen-fixture ATHDF, and promoted node-`cbin` helper regressions. |
| Focused verification | Python reader module returned `125 passed`; promoted MPI node example returned `1 passed`; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Full local verification | Full serial IO plus CPU smoke of the GPU-selectable regression returned `181 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts. |
| Deferred Pages verification | Refreshed the deferred visualization page, manifest, and compatibility contract; reapplied `tools/visualization.md` to detached `/private/tmp/athenak-gh-pages-io-docs-e948-final3`; rebuilt successfully with `make clean html SPHINXOPTS="-W --keep-going"` against `origin/gh-pages` `4833aa9341e19861297e330ff02aabfd8001935c`. |
| Status | Locally actionable corrections implemented; final independent Python, documentation, test/example, and scope re-audit dispositions plus coherent commits remain required before CP-08 closes. |

### 2026-05-29: CP-08 Exact-Geometry And PDF Edge-Preflight Corrections

| Field | Record |
| --- | --- |
| Python re-audit blockers | Per-logical geometry containment still admitted shifted or shrunken physical intervals. Modern explicit/generated and legacy PDF edge paths could construct Python or NumPy temporaries before their cumulative retained-memory peak was rejected. |
| Resolution | Require each MeshBlock geometry record to match its exact logical physical interval within floating-point tolerance. Preflight explicit-edge parsing, generated edges, monotonicity validation, and bin centers before NumPy materialization; retain cumulative accounting for other explicit dimensions and reference headers; add reduced-cap tests that prove rejection occurs before `np.array`, `np.fromstring`, or `np.diff`. |
| Focused verification | Python reader module returned `130 passed`; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Full local verification | Full serial IO plus CPU smoke of the GPU-selectable regression returned `186 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts. |
| Deferred Pages verification | Refreshed the deferred visualization page, manifest, and compatibility contract; reapplied `tools/visualization.md` to detached `/private/tmp/athenak-gh-pages-io-docs-e948-final3`; rebuilt successfully with `make clean html SPHINXOPTS="-W --keep-going"` against `origin/gh-pages` `4833aa9341e19861297e330ff02aabfd8001935c`. |
| Status | Locally actionable corrections implemented and local qualification passed; final independent re-audits, coherent commits, and scope closure remain required before CP-08 closes. |

### 2026-05-29: CP-08 Absolute-Tolerance And Legacy-ASCII Corrections

| Field | Record |
| --- | --- |
| Python re-audit blockers | Default `np.isclose` relative tolerance accepted shifted or neighboring MeshBlock intervals on large-offset domains. PDF ASCII splitting occurred before token-list preflight; legacy numeric rows accepted malformed trailing suffixes; legacy payload rows and final `vstack` copies were not bounded cumulatively before materialization. |
| Resolution | Use `rtol=0` and a storage-aware absolute tolerance capped below the logical block width; count ASCII tokens without splitting; preflight token strings, parsed floats, and NumPy arrays before conversion; parse numeric rows strictly; bound cumulative legacy payload rows and the final stacking peak; add large-offset, strict-suffix, reduced-cap, retained-reference, future-explicit-dimension, and `vstack` regressions. |
| Focused verification | Python reader module returned `139 passed`; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Full local verification | Full serial IO plus CPU smoke of the GPU-selectable regression returned `195 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts. |
| Deferred Pages verification | Refreshed the deferred visualization page, manifest, and compatibility contract; reapplied `tools/visualization.md` to detached `/private/tmp/athenak-gh-pages-io-docs-e948-final3`; rebuilt successfully with `make clean html SPHINXOPTS="-W --keep-going"` against `origin/gh-pages` `4833aa9341e19861297e330ff02aabfd8001935c`. |
| Status | Locally actionable corrections implemented and local qualification passed; final independent re-audits, coherent commits, and scope closure remain required before CP-08 closes. |

### 2026-05-29: CP-08 Legacy-Candidate And Spherical-Metadata Corrections

| Field | Record |
| --- | --- |
| Python re-audit blockers | Private legacy `.bins.pdf` header delegation dropped externally retained shard-reference bytes. Spherical-slice headers used unrestricted metadata-line reads and split variable tokens before applying a token-expansion preflight. |
| Resolution | Forward external bytes into legacy header parsing; cap cumulative spherical metadata bytes and individual header lines; count variable tokens without splitting and preflight token strings first; preserve downstream reduced-cap tests by isolating the newly earlier guard; add focused regressions. |
| Focused verification | Python reader module returned `142 passed`; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Full local verification | Full serial IO plus CPU smoke of the GPU-selectable regression returned `198 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts. |
| Deferred Pages verification | Refreshed the deferred visualization page, manifest, and compatibility contract; reapplied `tools/visualization.md` to detached `/private/tmp/athenak-gh-pages-io-docs-e948-final3`; rebuilt successfully with `make clean html SPHINXOPTS="-W --keep-going"` against `origin/gh-pages` `4833aa9341e19861297e330ff02aabfd8001935c`. |
| Status | Locally actionable corrections implemented and local qualification passed; final independent re-audits, coherent commits, and scope closure remain required before CP-08 closes. |

### 2026-05-29: CP-08 ASCII-Model And Spherical-Sibling Corrections

| Field | Record |
| --- | --- |
| Python re-audit blockers | PDF numeric rows accepted Unicode digits even though the peak model assumed ASCII-sized strings. Spherical variable-token accounting did not retain metadata summaries across sibling reads or reject duplicate `variables:` declarations. |
| Resolution | Reject non-ASCII PDF numeric rows before token accounting; reject duplicate spherical variable declarations; retain bounded private variable-metadata summaries; thread reference metadata and dense reconstruction budgets through sibling header and payload reads; account for duplicate-ownership peaks; release non-reference candidates promptly; add focused regressions. |
| Focused verification | Python reader module returned `145 passed`; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Full local verification | Full serial IO plus CPU smoke of the GPU-selectable regression returned `201 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts. |
| Deferred Pages verification | Refreshed the deferred visualization page, manifest, and compatibility contract; reapplied `tools/visualization.md` to detached `/private/tmp/athenak-gh-pages-io-docs-e948-final3`; rebuilt successfully with `make clean html SPHINXOPTS="-W --keep-going"` against `origin/gh-pages` `4833aa9341e19861297e330ff02aabfd8001935c`. |
| Status | Locally actionable corrections implemented and local qualification passed; final independent re-audits, coherent commits, and scope closure remain required before CP-08 closes. |

### 2026-05-29: CP-08 Mixed-MeshBlock-Extent Correction

| Field | Record |
| --- | --- |
| Python re-audit blocker | A malformed `.bin` or `.cbin` file could contain differently shaped MeshBlock records. Direct readers accepted the file, while athdf-like reconstruction reused the first block's output extent and could silently leave trailing cells initialized to zero. |
| Resolution | Require one uniform emitted MeshBlock extent inside every nonempty binary or coarsened-binary file in the shared decoder; reject a later mismatch before exposing file data; add direct malformed two-record regressions for both public readers; document the invariant in the compatibility contract and deferred Pages material. |
| Focused verification | Python reader module returned `148 passed`; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Full local verification | Full serial IO plus CPU smoke of the GPU-selectable regression returned `204 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts. |
| Deferred Pages verification | Refreshed both `tools/visualization.md` and the inserted file-reference section in detached `/private/tmp/athenak-gh-pages-io-docs-e948-final3`; rebuilt successfully with `make clean html SPHINXOPTS="-W --keep-going"` against `origin/gh-pages` `4833aa9341e19861297e330ff02aabfd8001935c`. |
| Status | Locally actionable correction implemented and local qualification passed; final independent re-audits, coherent commits, and scope closure remain required before CP-08 closes. |

### 2026-05-29: CP-08 Local Closure And Final Independent Audit Disposition

| Field | Record |
| --- | --- |
| Independent audit disposition | C++ format/output registration, MPI/restart, Python tooling, tests/examples/fixtures, deferred Pages, and final three-branch scope auditors reported no remaining locally actionable findings after corrections and rechecks. The scope auditor confirmed that implementation commit `685d04f1` contains the final mixed-extent decoder guard and its `.bin`/`.cbin` regressions together. |
| Final local verification | Python reader module returned `148 passed`; full serial IO plus CPU smoke of the GPU-selectable regression returned `204 passed`; full MPI IO returned `59 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; targeted py_compile, flake8, and `git diff --check origin/main` returned no output. |
| Production contradiction searches | `rg -n "\\.assembled|StageNodeRestart|CopyFileRange|bin_convert_new" src vis inputs` returned no matches. The broad documentation/test search contains only intentional historical records, explicit no-staging guidance, rejected-API guidance, and negative assertions. |
| Deferred Pages disposition | Reapplied the reviewed visualization overlay and file-reference insertion to detached `/private/tmp/athenak-gh-pages-io-docs-e948-final3`; warnings-as-errors build passed against `origin/gh-pages` `4833aa9341e19861297e330ff02aabfd8001935c`; live `/Users/dbf75/Work/Research/AthenaK/athenak-DF` remained clean and untouched. |
| Remaining external gates | Execute the GPU-selectable IO regression on a CUDA-capable build and qualify per-node output plus direct restart routing under a scheduler on multiple physical nodes, including a genuinely empty or non-owning node. |
| Status | CP-08 is complete for all locally actionable work. Historical rows that mention pending final re-audits are superseded by this closing disposition. The branch is locally complete but must not be described as merge-ready until both external qualification gates pass. |

## Blocking Findings And Resolutions

| ID | Finding | Severity | Resolution | Re-Audit |
| --- | --- | --- | --- | --- |
| PY-001 | Gotham converter variants call missing `write_athdf`, breaking supported base HDF5 conversion if promoted. | Blocking | Preserve base `write_athdf`/`write_xdmf_for`/CLI and merge modern behavior into canonical `bin_convert.py`; implementation delegated. | Resolved; public symbols and conversion tested. |
| PY-002 | Gotham PDF and spherical-slice readers insufficiently validate malformed or inconsistent shards. | Blocking | Require enhanced new readers to validate empty/truncated payloads, metadata/time consistency, and ownership/index correctness. | Resolved for shared/rank/node tested paths. |
| PY-003 | Canonical converter formed output names using global `.replace(".bin", "")`, corrupting identifiers such as `bin_shared` and failing to normalize `.cbin`. | Major | Strip only the final input suffix with `os.path.splitext`; exercise shared/rank conversion through the CLI. | Resolved; serial CLI regression passes. |
| TEST-001 | Current legacy PDF format is text-form and incompatible with assuming a modern binary baseline. | Blocking | Freeze legacy `.pdf`/`.bins.pdf` outputs from clean `origin/main`; declare N-D storage as new-reader compatible. | Resolved for 1-D and 2-D; byte-identical producer regressions pass. |
| CPP-001 | Existing `file_type=sph` VTK behavior could be accidentally replaced by new `sphslice`. | Blocking | Retain `sph` unchanged; add `sphslice` as separate format and test both compatibility and new readback. | Resolved; `sphslice` shared/rank/node readback and legacy `sph` filename/parameter regression passed. |
| CPP-002 | Derived-variable/PDF history before `abefea9f` contains indexing/allocation and device-memory correctness defects. | Blocking | Treat `abefea9f` corrections as mandatory in the clean core output implementation. | Implemented; mixed-derived PDF builds and executes in serial/MPI. |
| CPP-003 | `sphslice` trilinear interpolation can access ghost cells while derived arrays are populated only over active cells. | Blocking | Reject derived-array `sphslice` variables until ghost-zone-safe angular sampling is implemented; document native-only contract. | Resolved; negative construction regression passes and independent re-audit accepted guard. |
| CPP-004 | Generic transformed flow/flux names silently select inconsistent fluids when both Hydro and MHD are active under `<ion-neutral>`. | Blocking | Reject generic `mdot_*`/`edot_*`/`vel_*` for ion-neutral runs until module-qualified semantics exist. | Resolved; negative two-fluid PDF regression passes and independent re-audit accepted guard. |
| GPU-001 | The local qualification environment does not demonstrate device-memory correctness for new derived-variable paths. | Merge qualification | Add an IO `_gpu` regression using scalar and spherical derived PDF axes; run it in a GPU-capable build before merge readiness. | Test target added; GPU execution remains required. |
| MPI-001 | Gotham per-node restart may expose a manifest before every payload write is complete. | Blocking | Add atomic payload publication, inventory validation, then atomic public-manifest publication. | Resolved locally: native direct manifest-only resume, inventory/path/header negatives, alias rejection, and forced-small-chunk restart pass. |
| MPI-002 | Gotham performs global node communicator setup even when node sharding is unused and has ambiguous empty-shard behavior. | Blocking | Implement opt-in/lazy node-shard infrastructure and explicit valid-empty or skipped-shard reader contract. | Resolved for shared-node execution; scheduler-backed multi-node empty-shard run remains a production qualification item. |
| MPI-003 | Existing MPI-file truncation is racy and wrapper transfers narrow 64-bit sizes into MPI `int` counts. | Blocking | Implement synchronized one-rank truncation plus overflow-safe chunked byte IO, broadcast chunking, and forced-small-chunk tests during CP-02. | Resolved; independent re-audit accepted range preflight and wrapper harness correction. |
| TEST-002 | Existing sliced `cbin` producer emits a zero-width meshblock extent that canonical readback rejects even in shared mode. | Scope boundary | Reject sliced `cbin` consistently in every shard mode; keep full-volume `.cbin` equality testing and sliced `.bin` empty-owner testing. | Resolved as an explicit construction-time exclusion. |
| TEST-003 | Repository ignore rules omitted twelve frozen `.bin` and `.rst` fixture payloads from the first preservation commit. | Blocking | Force-add only the immutable fixture payloads and amend the test-fixture commit. | Resolved in `fcc534fe`; frozen-snapshot re-audit passed. |
| RST-001 | Existing per-rank restart resume calls `IOWrapper::GetPosition()` without the per-rank mode flag, invoking MPI-IO on a standard file handle. | Blocking | Forward `single_file_per_rank` in `Mesh::BuildTreeFromRestart` and add an MPI per-rank resume regression. | Resolved; two-rank resume regression passes. |
| DOC-001 | Live Pages pages advertise converter/PDF behavior inconsistent with clean code baseline. | Blocking | Stage a deferred overlay that reconciles existing pages after code stabilizes; do not publish before code merge. | Resolved locally: detached CP-07 application and warnings-as-errors build passed against refreshed Pages baseline; publication remains deferred. |
| DOC-002 | Planning records described node restart loading as native and used stale pre-transactional payload names. | Blocking | Describe current strict manifest validation plus transient `.assembled` staging truthfully; keep manifest path as the only supported restart entry point until CP-03. | Resolved: CP-01 corrected the intermediate state and CP-07 revalidated native direct-loading docs after CP-03. |
| RST-002 | Generated payload-path aliases could bypass manifest-only restart entry policy. | Blocking | Normalize lexically, resolve canonical symlink aliases when possible, reject temporary generated payload artifacts, and add alias regressions. | Resolved; focused restart suite returned `36 passed`. |
| PY-004 | File-controlled declared shard counts could cause unbounded expected-ID allocation; assembled node binary metadata retained a shard-local MeshBlock count. | Blocking | Compare declared and discovered counts before bounded ID construction and aggregate assembled MeshBlock metadata. | Resolved; CPU writer-reader hardening returned `17 passed`. |
| CPP-005 | Sliced node `.cbin` rejection occurred after data loading, so one-rank execution could fail for the wrong reason. | Blocking | Reject sliced `.cbin` consistently during construction and parameterize the MPI regression for one and two ranks. | Resolved; MPI writer hardening plus serial construction negatives pass. |
| RST-003 | Generated node-payload bytes could be opened directly through hard links or copies, bypassing the manifest-only transaction boundary; node-payload header writes and close returns were not checked; and declared payload inventory was unbounded. | Blocking | Mark node payload bytes explicitly, consume the marker only under validated manifest loading, reject marked bytes under ordinary restart input, check publication writes and closes, bound manifest inventory, and add focused negatives. | Resolved locally; native restart suite returned `40 passed`; final independent re-audit pending. |
| PY-005 | Sparse PDF shards lacked strict sibling inventory metadata, and shipped readers accepted malformed aliases, identifier mismatches, and file-controlled dense allocations without practical bounds. | Blocking | Emit PDF shard IDs/counts, publish modern PDF files atomically per file, enforce canonical sibling directories and identifier matching, and cap practical allocations before construction. | Resolved locally; CPU reader/writer matrix returned `70 passed`; final independent re-audit pending. |
| TEST-004 | Frozen legacy restart fixtures were checksummed but never resumed, strict node-PDF inventory lacked a positive synthetic two-shard reconstruction, fixture integrity was manual-only, and the shipped `cbin` helper branch lacked executable coverage. | Blocking | Add frozen shared/per-rank resume tests, synthetic two-node sparse-PDF reconstruction with an explicit empty shard, automatic manifest/checksum verification, and a `cbin` summary-helper regression. | Resolved locally; focused matrices pass; final independent re-audit pending. |
| PY-006 | Declared PDF V2 files could fall back to transitional parsing, V2 cycles were not checked against headers, reconstructed PDF/spherical-slice metadata retained shard-local state, and bulk-read limits were incomplete. | Blocking | Require V2 preambles and matching cycles when declared; normalize aggregate metadata; cap PDF header/payload reads, spherical-slice whole-file reads, and binary metadata, accumulated payload, and aggregate reconstruction sizes. | Resolved locally; Python reader/writer hardening returned `83 passed`; final independent re-audit pending. |
| RST-004 | Malformed node manifests could append unbounded zero-length segment records even after payload inventory was capped. | Blocking | Bound segment records before storage by declared MeshBlock count and a practical ceiling; reject non-positive segments; add focused negatives. | Resolved locally; affected MPI matrix returned `53 passed`; final independent re-audit pending. |
| TEST-005 | MPI PDF reconstruction regression still expected shard-local identifiers after aggregate metadata normalization. | Blocking | Inspect raw shard headers for local identity and assert reconstructed aggregate metadata retains totals without a local shard ID. | Resolved locally; affected MPI matrix returned `53 passed`; final independent re-audit pending. |
| PY-007 | Reader preflight remained incomplete for mixed PDF sibling V2 declarations, cumulative binary metadata, metadata-only MeshBlocks, late aggregate checks, and spherical coordinate/offset inputs. | Blocking | Validate sibling `binary_magic` and local headers; bound cumulative binary metadata, positive variables, fixed MeshBlock metadata, incremental aggregate state, spherical coordinate arrays, and input-header offsets. | Resolved locally; reader/writer matrix returned `92 passed`; final independent re-audit pending. |
| TEST-006 | Fixture automation ignored nested artifacts named like root metadata and generated Python caches repeatedly dirtied audit status. | Major | Compare fixture inventory exclusions by root-relative path; ignore Python bytecode and pytest caches; remove generated artifacts before staging. | Resolved locally; focused fixture regression, exact 27-artifact comparison, and cleanup pass. |
| DOC-003 | Deferred Pages overlay carried unrelated CGM navigation and misstated PDF writer-rank plus binary incremental-bound details; publication record omitted frozen restart resumes. | Blocking | Remove unrelated navigation, correct payload writer-rank semantics, narrow incremental aggregate wording, require frozen restart resume status, rebuild detached Pages candidate. | Resolved locally; warnings-as-errors detached Pages rebuild passed; final docs re-audit pending. |
| PY-008 | Retained sparse sibling arrays, widened binary metadata, aggregate-copy peaks, unbounded athdf-like allocations, spherical header-only offsets, and non-positive coarsening factors remained outside preflight. | Blocking | Incrementally combine spherical shards; retain inventory summaries only; budget retained/copy bytes; validate offsets and coarsening; add shared athdf-like allocation preflights. | Resolved locally; reader/writer matrix returned `95 passed`; final independent re-audit pending. |
| DOC-004 | Deferred pages omitted historical read compatibility for transitional unversioned binary PDFs and overstated additive V2 inventory metadata as mandatory for those artifacts. | Major | State that writers emit V2 while readers retain transitional dense/sparse compatibility and only new V2 sparse shards declare complete inventory. | Resolved locally; detached Pages rebuild passed; final docs re-audit pending. |
| PY-009 | Athdf-like helpers retained placeholder restriction branches, materialized full prolonged blocks before selection, accepted malformed logical levels, and discarded single-rank destination arrays; PDF/spherical retained-memory checks were incomplete and PDF metadata accepted non-finite values. | Blocking | Centralize bounded MeshBlock placement, crop prolongation, restore restriction, validate logical metadata, assign into supplied destinations, require finite PDF metadata, and budget cumulative retained arrays. | Resolved locally; reader/writer matrix returned `115 passed`; final independent re-audit pending. |
| CPP-006 | Passive-scalar group labels wrapped after scalar `99`, and PDF symlog thresholds accepted non-finite values. | Blocking | Emit non-wrapping scalar labels, reject duplicate binary labels defensively, require finite symlog thresholds, and add producer/reader regressions. | Resolved locally; serial producer surface returned `146 passed`; final independent re-audit pending. |
| CPP-007 | `cbin` accepted factors outside its documented power-of-two range and could reach division before rejecting malformed input; range-valid factors could still produce incompatible ghost-expanded extents. | Blocking | Validate the factor range before writer construction, validate emitted-extent divisibility during writer construction, reject sliced output consistently, and reject malformed factors before reader division. | Resolved locally; producer and reader negatives pass; final independent re-audit pending. |
| PY-010 | V2 sparse PDF readers accepted missing writer-mandatory distribution, local-ID, or sibling-count inventory metadata. | Blocking | Require complete sparse V2 inventory declarations while retaining permissive transitional unversioned compatibility. | Resolved locally; V2 inventory negatives and transitional compatibility pass; final independent re-audit pending. |
| CPP-008 | Ordered infinite PDF axis bounds passed the producer range check and could publish malformed metadata with non-finite bin edges. | Blocking | Require finite lower and upper bounds for every PDF dimension before edge construction and publication; add an `inf` producer negative. | Resolved locally; focused producer module returned `19 passed`, full serial IO returned `146 passed`, and final independent producer re-audit accepted the fix. |
| PY-011 | Ghost-zone ATHDF offsets and coordinates, intersecting-cell lower crops, malformed-grid validation, partial-shard level initialization, cumulative ATHDF/PDF/spherical live peaks, and dense PDF distribution validation remained incomplete. | Blocking | Validate grid topology before reconstruction; centralize correct ghost/crop placement; initialize uncovered levels to `-1`; add cumulative allocation preflight at each materialization peak; validate modern dense PDF distribution metadata; add focused regressions. | Resolved locally; Python reader module returned `105 passed`, full serial IO returned `161 passed`, full MPI IO returned `59 passed`, style returned `2 passed`; final independent Python re-audit pending. |
| PY-012 | Omitted ghost counts, direct duplicate/out-of-domain MeshBlocks, prolongation indexes, center-list construction, symlog transformation arrays, and final spherical-coordinate temporaries remained outside strict reconstruction validation or live-peak budgets. | Blocking | Reject omitted ghost declarations and malformed direct inventories; preallocate centers; bound prolongation indexes, symlog transforms, and spherical coordinate generation; add focused regressions. | Resolved locally; Python reader module returned `111 passed`, full serial IO returned `167 passed`, full MPI IO returned `59 passed`, style returned `2 passed`; final independent Python re-audit pending. |
| PY-013 | PDF local shard headers, per-logical binary geometry, singleton-axis ghost widths, ATHDF NumPy source arrays, sparse duplicate-validation arrays, and the preserved single-MeshBlock ATHDF helper remained outside strict lifecycle, live-peak, or unsupported-extent validation. | Blocking | Release incorporated sibling state, account for retained replacement-header state privately, validate logical geometry envelopes and singleton widths, budget source and duplicate-validation temporaries, reject unsupported single-MeshBlock emitted extents without changing its signature, and add focused regressions. | Resolved locally; Python reader module returned `125 passed`, full serial IO returned `181 passed`, full MPI IO returned `59 passed`, style returned `2 passed`; final independent Python and documentation re-audits pending. |
| TEST-007 | Retained ATHDF orchestration APIs and the documented node-sharded `cbin --assemble-shards` helper path lacked positive end-to-end automation. | Major | Add frozen-fixture ATHDF equality regressions and invoke the node `cbin` helper from the promoted MPI example. | Resolved locally; reader module returned `125 passed`, promoted MPI node example returned `1 passed`, full serial IO returned `181 passed`, and full MPI IO returned `59 passed`; final independent test/example re-audit pending. |
| PY-014 | Per-logical geometry containment accepted shifted or shrunken records, and PDF edge parsing/generation/validation could materialize Python or NumPy temporaries before cumulative retained-memory rejection. | Blocking | Require exact logical physical intervals; preflight modern explicit/generated and legacy edge parsing, validation, and center peaks before materialization; include later explicit dimensions and retained reference headers in live budgets; add pre-materialization regressions. | Resolved locally; Python reader module returned `130 passed`, full serial IO returned `186 passed`, full MPI IO returned `59 passed`, style returned `2 passed`, fixture checksum verification passed for all `27` artifacts, and detached Pages warnings-as-errors build passed; final independent Python re-audit pending. |
| PY-015 | Default relative geometry tolerance admitted shifted large-offset blocks; PDF token lists formed before preflight; legacy ASCII rows accepted malformed suffixes; cumulative payload rows and final stacks materialized without retained-peak bounds. | Blocking | Use zero relative tolerance with a storage-aware bounded absolute tolerance; pre-count and budget tokens before splitting; parse ASCII rows strictly; account for cumulative payload rows and final `vstack`; add boundary regressions. | Resolved locally; Python reader module returned `139 passed`, full serial IO returned `195 passed`, full MPI IO returned `59 passed`, style returned `2 passed`, fixture checksum verification passed for all `27` artifacts, and detached Pages warnings-as-errors build passed; final independent Python re-audit pending. |
| TEST-008 | Generated-edge peak rejection was implemented before NumPy edge construction but lacked a fail-fast ordering regression. | Major | Replace edge generation with a failing stub under a reduced cap and assert the retained-peak guard rejects first. | Resolved locally; Python reader module returned `139 passed`, full serial IO returned `195 passed`, and style returned `2 passed`; final independent test/example re-audit pending. |
| PY-016 | Legacy PDF candidate parsing dropped retained reference bytes, and spherical-slice header lines plus variable-token expansion were not bounded before materialization. | Blocking | Forward private legacy-header retained budgets; cap cumulative spherical metadata and line reads; pre-count and budget variable tokens before splitting; add focused regressions. | Resolved locally; Python reader module returned `142 passed`, full serial IO returned `198 passed`, full MPI IO returned `59 passed`, style returned `2 passed`, fixture checksum verification passed for all `27` artifacts, and detached Pages warnings-as-errors build passed; final independent Python and documentation re-audits pending. |
| PY-017 | PDF token parsing accepted Unicode digits outside its ASCII peak model, and spherical sibling variable metadata was not retained or released consistently across parsing, payload, and ownership peaks. | Blocking | Reject non-ASCII PDF numeric rows; reject duplicate spherical variable declarations; retain bounded metadata summaries; budget sibling coexistence through payload and ownership checks; release non-reference candidates promptly; add focused regressions. | Resolved locally; Python reader module returned `145 passed`, full serial IO returned `201 passed`, full MPI IO returned `59 passed`, style returned `2 passed`, fixture checksum verification passed for all `27` artifacts, and detached Pages warnings-as-errors build passed; final independent Python re-audit pending. |
| TEST-009 | Private spherical retained-budget regression did not prove that public multi-shard reconstruction forwards retained state into later sibling reads. | Major | Extend the public two-shard lifecycle regression to assert a nonzero retained-byte argument after the first sibling. | Resolved locally; Python reader module returned `145 passed`, full serial IO returned `201 passed`, and style returned `2 passed`; final independent test/example re-audit pending. |
| PY-018 | Final spherical coordinate-generation preflight removed the retained variable-metadata summary before constructing theta/phi arrays even though the returned header still owned those strings. | Blocking | Preserve and charge the private metadata summary until after coordinate construction; add a focused reduced-cap regression. | Resolved locally; Python reader module returned `146 passed`, full serial IO returned `202 passed`, full MPI IO returned `59 passed`, style returned `2 passed`, fixture checksum verification passed for all `27` artifacts, and detached Pages warnings-as-errors build passed; final independent Python and documentation re-audits pending. |
| PY-019 | Binary and coarsened-binary direct readers accepted mixed emitted MeshBlock extents within one file even though athdf-like reconstruction derived one block size from the first record. | Blocking | Reject later MeshBlock records whose emitted extents differ from the first record in the shared decoder; add malformed two-record `.bin` and `.cbin` regressions. | Resolved and independently accepted; Python reader module returned `148 passed`, full serial IO returned `204 passed`, full MPI IO returned `59 passed`, style returned `2 passed`, fixture checksum verification passed for all `27` artifacts, detached Pages warnings-as-errors build passed, and final Python, docs, tests/examples, and scope rechecks reported clean. |

## Robustification Pass

### Approved Guide Freeze

| Field | Record |
| --- | --- |
| Approval instruction | User requested execution of `IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md` in its entirety and explicitly requested subagent review. |
| Frozen path | `IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md` |
| Process-only commit | `e40621d81e1948567415a1435f531b275fe14c2e` |
| SHA-256 | `cdf53351104c135f2d79f9ee2f36c2908001b00a31339db18eb1cf66bcae11ff` |
| Line count | `2222` |
| Byte count | `89840` |
| Approval date | `2026-05-29 EDT` |
| Rule | Recompute the checksum, line count, and byte count at every work session and checkpoint. Stop on drift until a changed guide is explicitly approved, frozen, and re-audited. |

### Live Robustification Checkpoint Board

| Checkpoint | Status | Decisions | Pre-edit auditors | Implementation commit(s) | Focused tests | Post-edit auditors | Reflection | Remaining risk |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RCP-00 | Locally closed | D-069, D-070, D-071, D-073, D-075, D-081, D-083, D-085 | Runtime baseline; format/tooling baseline; documentation/process baseline | `e40621d8` guide freeze | Fresh floors registered below | Three baseline reports accepted | R-00 recorded below | External CUDA and multi-node topology remain unavailable locally |
| RCP-01 | Locally closed | D-070, D-071, D-085, D-091, D-092, D-093 | Focused layout auditor accepted | `3b96e2d1` | Full serial `238 passed`; full MPI `67 passed`; style `2 passed`; fixtures `27` verified; exact pre-edit byte comparisons pass | Arithmetic, compatibility, topology, admission, and single-issue re-auditors accepted after corrections | R-01 recorded below | Signed output-sequence domain remains queued for RCP-02 |
| RCP-02 | Locally closed | D-072, D-073, D-074, D-086, D-087 accepted | MPI call inventory; filesystem-publication inventory; sequence and namespace inventory accepted | `e6e6ee7c` runtime hardening | Focused publication `58 passed`; injected cleanup `2 passed`; full local MPI `76 passed`; style `2 passed` | Three publication audits accepted after corrections | R-02 recorded below | Real cross-node cleanup ordering and production-filesystem behavior remain RCP-09 gates |
| RCP-03 | Locally closed | D-094 and D-095 accepted; supersede D-075 and resolve D-088 | Kernel-range, producer-contract, and reader-compatibility audits accepted direction | `e6e6ee7c` runtime hardening | Layout, node64, reader, conversion, explicit adaptive-AMR rejection, and rank/node sliced rejection pass | Fresh final acceptance audit accepted | R-03 recorded below | External multi-node writer qualification remains open |
| RCP-04 | Locally closed | D-076, D-077, D-078, D-089, D-096, D-097, D-112, D-116 accepted | Numerical semantics, analytic-test, backend portability, and late adaptive-shape/staging audits completed | `e6e6ee7c` runtime hardening | Serial format matrix and mixed-level oracle pass; sharded sparse-overflow MPI `3 passed`; repeated shared/rank/node smoke passed | Fresh numerical/MPI re-auditor accepted; mixed-level MPI composition retained as nonblocking enhancement | R-04 refreshed below | Representative CUDA execution remains an RCP-09 gate |
| RCP-05 | Locally closed | D-079 accepted | Bounded registration cleanup and tracked-particle repair implemented | `e6e6ee7c` runtime hardening | Parser-sensitive serial `88 passed`; MPI format hardening `8 passed`; tracked serial `5 passed`; tracked MPI `3 passed`; source audit `1 passed`; style `2 passed`; post-correction focused serial `113 passed`; post-correction MPI `63 passed` | Behavior-preservation auditor accepted; scope auditor accepted after narrow correction | R-05 recorded below | Tracked-particle bytes remain intentionally native-endian legacy format |
| RCP-06 | Locally closed; committed-tree rerun passed | D-080, D-081, D-108, D-110, D-111, D-114, D-118 accepted | Python API and memory-budget audits completed | `e6e6ee7c` runtime hardening | Focused readers and writer hardening `298 passed`; committed canonical serial `639 passed, 4 skipped`; committed canonical MPI `85 passed, 3 skipped`; committed collection `731 tests` | Fresh Python/file-format re-auditor accepted checked-in fail-closed boundary coverage | Closure refreshed below | Public limits remain keyword-only overrides and CLI flags |
| RCP-07 | Locally closed again; post-merge restaging required | D-082, D-098, D-099, D-104, D-107 accepted | Live Pages structure and protected-blob inventory established | `6d64ab74` deferred Pages integration | Helper `25 passed`; nine-file strict detached stage, build, linkcheck, drift packet, and rendered browser QA passed | Fresh file-format/Pages re-auditor accepted exact allowlist, hashes, prose, and helper assertions | R-07 refreshed below | Publish only after code merge through a separate Pages review |
| RCP-08 | RCP-08A locally closed; RCP-08B external | D-083 direction retained; D-090 accepted | Benchmark-design auditor accepted corrected instrumentation and preregistration | `e6e6ee7c` instrumentation; process snapshot pending | Focused timing `3 passed`; full node-sharding MPI `58 passed`; style passed | Fresh benchmark-design correction auditor accepted | R-08A closure recorded below | RCP-08B requires scheduler-backed measurements |
| RCP-09 | Blocked externally; local scheduler-deck tooling accepted | D-109, D-115, D-117, D-119 through D-134, D-136, D-138, D-140, and D-141 accepted; D-135, D-137, and D-139 rejected | Local environment and evidence-plan auditors completed | `f5c29fc5` process baseline; `5ee873e2` archive-admission and descriptor-hygiene correction | Local environment probe recorded; Slurm runner syntax-check and expanded macOS Python/Bash-3 mock lifecycle matrix passed; committed strict runner plus packet-finalizer suites `116 passed`; immutable local evidence packet finalized idempotently | Fresh narrow and broad process-evidence acceptance audits accepted the settled tree with no P1/P2/P3 findings | Plan corrected after durable-publication, crash-consistency, archive-admission, and descriptor-hygiene audits | No CUDA/HIP toolchain, scheduler launcher, second physical host, or attributable deployment-filesystem measurements locally |
| RCP-10 | Committed-tree local closeout passed; external gates remain | D-084, D-100 through D-134, D-136, D-138, D-140, and D-141 accepted; D-135, D-137, and D-139 rejected | Final lane auditors found bounded reader, diagnostic, timeout, package-summary, restart-payload, rendered-doc, transient-admission, evidence-boundary, token-expansion, shard-bookkeeping, wrapper-default, commit-point, intrinsic-header, coarse-fine-oracle, serialized-narrowing, external-packet, header-only-contract, retained-timer-grammar, MR-2-hostfile, scan-terminality, rank-map, packet-publication, canonical-inventory, permission-scan, structured-index, packet-tree-alias, rank-map-alias, lifecycle-prefix, dangling-archive-symlink, outer-identity, post-launch-alias, accounting-history, canonical-archive-identity, timing-publication, accounting-hook-window, timing-parent-alias, failing-hook-index-alias, packet-root-alias, packet-replacement, rejected-append-mutation, late-read-only-artifact, late-archive-sink-alias, child-object-replacement, cross-filesystem-metadata-publication, short-append, process-termination, parent-sync-retry, ambiguous-reserved-temp, archive-adjacent-temp, and imported-descriptor-hygiene blockers | `e6e6ee7c` runtime; `6d64ab74` Pages; `f5c29fc5` process baseline; `5ee873e2` archive-admission and descriptor-hygiene correction | Readers `298 passed`; spherical slices `14 passed`; restart MPI `7 passed`; sharded sparse-overflow MPI `3 passed`; scheduler/finalizer mock `116 passed`; committed canonical serial `639 passed, 4 skipped`; committed canonical MPI `85 passed, 3 skipped`; committed collection `731 tests`; tracked particles `5 passed` and `3 passed`; deferred Pages helper `25 passed`; style `2 passed`; frozen fixtures `27` verified; static gates passed; detached Pages strict verification and rendered browser QA passed | Fresh narrow and broad whole-branch auditors accepted settled tree `1e474976fddd019d7d237b9b81828b99844621b3` with no P1/P2/P3 findings | R-10 committed-tree evidence packet recorded below | Cannot close merge-readiness gate before external CUDA, optional deployment-specific HIP, physical multi-node, filesystem, and RCP-08B evidence |

### Expanded Robustification Finding Register

The frozen guide defines `ROB-001` through `ROB-022`. `RCP-00` baseline auditors
added the following durable rows. Do not close or merge them silently into broad
cleanup claims.

| ID | Severity | Finding | Primary checkpoint |
| --- | --- | --- | --- |
| ROB-023 | P2 | Restart parameter-header scanning ignores a failed `Seek()` reposition after reading the parameter dump | RCP-01 |
| ROB-024 | P2 | Signed `file_number` accepts negative values and can overflow at `INT_MAX`; widened rendering alone does not define a safe sequence-counter domain | RCP-02 |
| ROB-025 | P2 | `cbin` `gid` filtering compares a pack-local index against a global MeshBlock ID and can select the wrong block when a pack offset is nonzero | RCP-03 |
| ROB-026 | P3 | Public Python converter paths still need explicit allocation preflight and deliberate budget override behavior | RCP-06 |
| ROB-027 | P1 docs | Literal deferred Pages whole-page replacement can delete unrelated live `origin/gh-pages` guidance even when Sphinx builds successfully | RCP-07 |
| ROB-028 | P1 | Restart metadata reconstruction allocates and broadcasts `listsize * nmb_total` bytes with unchecked arithmetic and advances through the buffer with a signed cursor | RCP-01 |
| ROB-029 | P1 | Restart startup consumes serialized mesh dimensions, ghost counts, refinement state, and coordinate bounds before rejecting structurally invalid values | RCP-01 |
| ROB-030 | P2 | Serial positioned IO validates the starting offset but not the terminal offset, and serial byte operations narrow counts to `size_t` without representability checks | RCP-01 |
| ROB-031 | P2 | Node-restart writer can publish a manifest above the strict reader's 64 MiB limit and can generate payload paths above the parser's accepted path bound | RCP-01 |
| ROB-032 | P2 | Restart allocation and MHD face-stride calculations need checked `size_t` preflight and checked subtraction as well as checked add and multiply operations | RCP-01 |
| ROB-033 | P1 | Restart leaf metadata can reach tree shifts, refinement recursion, and load balancing without validating topology completeness, canonical persisted order, signed-safe logical levels, positive finite costs, or rank-map completion | RCP-01 |
| ROB-034 | P1 | Multilevel restart headers can accept odd active MeshBlock extents even though scratch construction rejects them before coarse-grid refinement arithmetic | RCP-01 |
| ROB-035 | P2 | Adaptive `num_levels + root_level - 1` arithmetic can overflow a signed `int` before restart or scratch admission checks reject the requested refinement depth | RCP-01 |
| ROB-036 | P1 | Serialized restart dimensions can disagree with dimension flags derived from the persisted parameter dump, allowing self-consistent but incompatible header metadata into tree and physics setup | RCP-01 |
| ROB-037 | P2 | Adaptive `num_levels` still narrows through `atoi()` inside `GetOrAddInteger()` before wide logical-level validation unless mesh construction parses the textual value directly | RCP-01 |
| ROB-038 | P1 | Inactive coarse MeshBlock axes can retain noncanonical serialized bounds even though scratch construction initializes them to `0/0` | RCP-01 |
| ROB-039 | P1 | Branch-added and inherited-but-touched MPI calls still have unchecked returns or fragmented low-context reporting across node communicator setup, output timing, reductions, and restart publication | RCP-02 |
| ROB-040 | P1 | Rank-local fatal exits in touched output and restart-read paths can strand peers in later world or node collectives, including output-timing reduction and restart metadata broadcasts | RCP-02 |
| ROB-041 | P1 | `.bin`, `.cbin`, and shared or per-rank restart writers publish directly to public paths while modern PDF, `sphslice`, and node restart use temporary-file rename publication | RCP-02 |
| ROB-042 | P2 | Node-restart generation names are clock-derived without exclusive reservation; collision and rename-failure paths can reuse or strand owned temporary artifacts | RCP-02 |
| ROB-043 | P2 | Output path components are not validated consistently before writers create directories and public namespaces | RCP-02 |

The `ROB-017` inventory is expanded to include inherited table and VTK writers in
addition to `.bin`, `.cbin`, PDF, `sphslice`, and restart output.

### 2026-05-29: RCP-00 Reopen Baseline And Freeze Facts

| Field | Record |
| --- | --- |
| Session-start refresh | `git fetch --prune origin` returned success before local inspection. |
| Snapshot | `HEAD=e40621d81e1948567415a1435f531b275fe14c2e`; `origin/main=886dd2a1437e45a3a30b3eeebf2adfa838328f73`; `origin/gh-pages=4833aa9341e19861297e330ff02aabfd8001935c`; merge base equals `origin/main`. |
| Worktree classification | Clean after the guide-only freeze commit. Existing worktrees were enumerated. Temporary Debug builds live under `/tmp/athenak-io-robust-build` and `/tmp/athenak-io-robust-build-mpi`, outside the repository. |
| Guide verification | Checksum, line count, and byte count match the approved frozen-guide record above. |
| Independent auditors | Runtime baseline auditor; format/tooling baseline auditor; documentation/process baseline auditor. All were read-only and reported clean repository status. |
| Accepted recommendations | Retain all seeded `ROB-*` findings; add `ROB-023` through `ROB-027`; include `src/parameter_input.cpp` in `RCP-01`; expand sequence-format inventory to table and VTK writers; keep Pages publication deferred until preservation-aware staging exists; register external-build-safe baseline commands. |
| Rejected recommendations | None. Recommendations that belong to later checkpoints remain queued rather than implemented during `RCP-00`. |
| Immutable fixtures | `(cd tst/fixtures/io/origin_main_886dd2a1 && shasum -a 256 -c SHA256SUMS)` returned `OK` for all `27` artifacts. |
| Python syntax | `python3 -m py_compile vis/python/bin_convert.py vis/python/read_pdf.py vis/python/read_sphslice.py vis/python/examples/read_io_outputs.py tst/test_suite/io/*.py` passed. |
| Pytest collection | `python3 -m pytest --collect-only -q tst/test_suite/io` reported `263 tests collected`. |
| Fresh serial Debug build | `cmake -S . -B /tmp/athenak-io-robust-build -DCMAKE_BUILD_TYPE=Debug -DAthena_ENABLE_MPI=OFF`; `cmake --build /tmp/athenak-io-robust-build -j 4`; both passed. |
| Fresh MPI Debug build | `cmake -S . -B /tmp/athenak-io-robust-build-mpi -DCMAKE_BUILD_TYPE=Debug -DAthena_ENABLE_MPI=ON`; `cmake --build /tmp/athenak-io-robust-build-mpi -j 4`; both passed. |
| Fresh serial floor | From `/tmp/athenak-io-robust-build/src` after linking `inputs -> <repo>/tst/inputs`: `python -m pytest -q <repo>/tst/test_suite/io/*_cpu.py <repo>/tst/test_suite/io/test_output_formats_gpu.py` returned `204 passed`. The `_gpu` module is a CPU smoke only on this host. |
| Fresh MPI floor | From `/tmp/athenak-io-robust-build-mpi/src` after linking `inputs -> <repo>/tst/inputs`: `python -m pytest -q <repo>/tst/test_suite/io/*_mpicpu.py` returned `59 passed`. |
| Fresh style floor | From `<repo>/tst`: `python -m pytest -q test_suite/style` returned `2 passed`. The style wrapper transiently downloads and removes `tst/test_suite/style/cpplint.py`. |
| Initial invocation correction | An initial isolated-build run linked `<repo>/inputs`, which lacks the flat regression decks expected by the pytest modules. Those failures were classified as invocation errors. Replacing the temporary link with `<repo>/tst/inputs` restored the registered floors without repository edits. |
| Local topology | `mpirun -np 2 hostname | sort -u` reported one physical host: `Tin-Drum`. No `srun`, `jsrun`, `aprun`, `qsub`, `sbatch`, `bsub`, or `flux` launcher was found. No `nvcc`, `nvidia-smi`, `hipcc`, or `rocm-smi` tool was found. |
| Status | `RCP-00` locally closed. Implementation may begin at `RCP-01`. `RCP-09` remains blocked locally and cannot be described as externally qualified. |

### Reflection R-00: Reassess The Robustification Direction

| Question | Record |
| --- | --- |
| Did auditors find defects not present in the guide? | Yes. Added `ROB-023` through `ROB-027`. |
| Are planned refactors too broad? | The risk is highest in restart and Python tooling. Use a small shared restart-layout descriptor, narrow checked helpers, and keyword-only Python overrides. Avoid broad format rewrites. |
| Should checkpoint order change? | No. Restart arithmetic remains the first implementation checkpoint. MPI/publication helpers remain a prerequisite for namespace work. Conservative `cbin` enforcement follows. |
| Did prior resolved findings become invalid? | Prior local evidence remains valid, but the new register reopens stronger robustness claims. Prior Pages HTML success is insufficient evidence because a literal overlay can remove unrelated live content. |
| What changed direction? | `RCP-07` must build a preservation-aware staging helper and review whole-page overlays against live `origin/gh-pages`; it must not mechanically replace live pages. |

### Canonical External-Build-Safe Local Matrix

Use these commands after linking each temporary build tree's `src/inputs` to
`<repo>/tst/inputs`. Record exact logs and exit codes for checkpoint closure.

| Matrix | Command | RCP-00 floor |
| --- | --- | --- |
| Python readers | From `<repo>`: `PYTHONDONTWRITEBYTECODE=1 python -m pytest -p no:cacheprovider -q tst/test_suite/io/test_python_io_readers_cpu.py` | `148 passed` |
| Serial IO and CPU smoke of GPU-selectable regression | From `/tmp/athenak-io-robust-build/src`: `python -m pytest -q <repo>/tst/test_suite/io/*_cpu.py <repo>/tst/test_suite/io/test_output_formats_gpu.py` | `204 passed` |
| MPI IO | From `/tmp/athenak-io-robust-build-mpi/src`: `python -m pytest -q <repo>/tst/test_suite/io/*_mpicpu.py` | `59 passed` |
| Style | From `<repo>/tst`: `python -m pytest -q test_suite/style` | `2 passed` |
| Fixtures | From `<repo>/tst/fixtures/io/origin_main_886dd2a1`: `shasum -a 256 -c SHA256SUMS` | `27` artifacts verified |

### Canonical External Qualification Matrix

The local matrix does not satisfy these production gates:

| Gate | Requirement | Local disposition |
| --- | --- | --- |
| CUDA | Configure a real CUDA-capable build and execute `tst/test_suite/io/test_output_formats_gpu.py` plus representative neighboring regressions | Blocked: no CUDA compiler or runtime tool on `Tin-Drum` |
| Multi-node | Execute per-node output equality, direct manifest restart, changed-rank restart, empty or non-owning node, timing labels, and scaling measurements on at least two physical scheduler-backed nodes | Blocked: one physical host and no scheduler launcher locally |

### 2026-05-29: RCP-01 Pre-Edit Restart-Layout Audit

| Field | Record |
| --- | --- |
| Auditor | Focused read-only restart-layout auditor. Repository status remained clean and no runtime tests or edits were performed. |
| Accepted design | Use a small shared restart-layout descriptor plus checked unsigned add and multiply helpers. Keep serialization loops, MPI routing, file opening, and module-specific state reads and writes local. This confirms D-070 and expands the boundary through D-091. |
| Persisted layout | Preserve `H = P + M + F + L + C + S + sizeof(IOWrapperSizeT)`, where `P` is the parameter dump, `M` is the node-only marker, `F` is the fixed mesh header, `L` is logical-location bytes, `C` is MeshBlock-cost bytes, `S` is optional Z4c tracker and turbulence state, and the final field stores per-MeshBlock payload bytes. Do not add serialized per-field descriptors. |
| Payload descriptor fields | Hydro cell-centered; MHD cell-centered; MHD `x1f`, `x2f`, and `x3f`; radiation; forcing; Z4c; and mutually exclusive ADM fallback. Calculate every field element count, byte count, cumulative payload size, rank offset, node offset, payload size, and local subview count with checked arithmetic. |
| New scope finding | Add `ROB-028`: `Mesh::BuildTreeFromRestart()` reconstructs replicated metadata with unchecked `listsize * nmb_total` allocation and broadcast arithmetic plus a signed cursor. Include `src/mesh/build_tree.cpp` in `RCP-01`. |
| Existing findings confirmed | `ROB-004`: unchecked serial `fseek()` and negative `ftell()` conversion. `ROB-016`: missing restart path logs and continues. `ROB-019`: unbounded manifest probe and parser lines. `ROB-023`: parameter scanning discards the final reposition result. |
| Compatibility constraints | Preserve legacy shared and per-rank restart bytes, node-only marker placement, manifest-only node entry semantics, and direct local-span routing. Do not broaden the descriptor into a restart rewrite. |
| Accepted recommendations | Include metadata reconstruction; preflight signed counts before unsigned conversion; replace `int` Kokkos subview-count narrowing; use descriptor field offsets for MHD face strides; add artificial-count, serial-wrapper, bounded-manifest, and immediate missing-path tests. |
| Rejected recommendations | None. Local Kokkos allocation preflight is accepted as part of the descriptor use sites rather than a separate abstraction. |
| Stop-gate status | Pre-edit layout map and abstraction decision are recorded. Runtime editing may begin. |

### 2026-05-29: RCP-01 Restart-Layout Implementation And Red-Team Corrections

| Field | Record |
| --- | --- |
| Implementation scope | Added the header-only `src/restart_layout.hpp` arithmetic facility; used it from restart write, restart read, and replicated metadata reconstruction; added checked POSIX serial positioning; made parameter-header reposition failures fatal; made missing restart input fail immediately; bounded node-manifest materialization; and added focused harnesses and regressions. |
| Arithmetic auditor findings | Added `ROB-029` through `ROB-032`: validate serialized mesh structure before division and allocation; reject serial terminal-offset and `size_t` narrowing overflow; enforce reader-writer manifest-limit symmetry; preflight memory-size representability and use checked subtraction for MHD face-stride remainders. |
| Corrections | Added structural restart-header validation in `Mesh::BuildTreeFromRestart()` before root-grid division; expanded serial wrapper checks; added `CheckedSizeT`, `CheckedMemorySize`, `CheckedSubtract`, and `ManifestBudget`; enforced generated payload-path and final manifest-byte budgets in the writer; stored checked MHD face-stride remainder bytes in the descriptor; and extended harness coverage. |
| Deferred finding | The auditor reconfirmed signed `file_number` truncation and overflow exposure across output writers. That is already registered as `ROB-024` and remains intentionally queued for the cross-format inventory in `RCP-02`; it is not silently closed by restart work. |
| Focused verification | Restart-layout and serial-wrapper CPU harnesses plus restart finalization tests returned `21 passed`; MPI restart and node-sharding tests returned `56 passed`; repository style returned `2 passed`; `git diff --check` returned no output. |
| Legacy byte preservation | Built detached pre-edit `9fe385b8b900f99beecdc9bd4440819ebd67aac1` serial and MPI executables. Deterministic generated shared restart files matched post-edit bytes exactly for file numbers `00000` and `00001` (`56183` bytes each). Deterministic generated per-rank restart files matched post-edit bytes exactly for both ranks and both file numbers (`34886` bytes each). |
| Compatibility auditor | Accepted preservation of shared and per-rank wire order, node marker placement, manifest-only node entry, direct span routing, and MHD field strides. The detached pre-edit comparison closes the auditor's request for producer-to-producer byte evidence. |
| Remaining RCP-01 action | Wait for the correction re-auditor and run the full local matrix before closing and committing the checkpoint. |

### Reflection R-01: Keep The Restart Facility Narrow

| Question | Record |
| --- | --- |
| Did the shared layout work reveal duplicated logic elsewhere? | Yes. Replicated metadata reconstruction and MHD face-stride arithmetic were part of the same restart safety boundary. They are covered by the shared checked helpers without moving serialization control flow. |
| Is the abstraction narrow enough? | Yes. `src/restart_layout.hpp` owns checked arithmetic and payload offsets only. Restart format order, MPI routing, file opening, Kokkos view construction, and module-specific loops remain in their existing owners. |
| Did any legacy fixture change unexpectedly? | No. Frozen resume coverage still passes, and detached pre-edit producer comparisons are byte-identical for generated shared and per-rank restart files. |
| Are wrapper-level chunk tests still representative after caller changes? | Yes. Caller arithmetic is checked before entry, and the existing forced-small-chunk node restart tests still exercise wrapper chunking. |
| Should any helper move into a general checked-IO utility for RCP-02? | Not yet. Keep restart-layout helpers scoped to restart until the MPI/publication inventory demonstrates an actually shared boundary. |

### 2026-05-29: RCP-01 Topology Re-Audit And Second Correction Cycle

| Field | Record |
| --- | --- |
| Correction re-audit blockers | The first correction re-auditor rejected payload-path validation after mutation, noncanonical positive mesh spacing, and unvalidated leaf metadata before tree shifts. The metadata-domain sidecar then identified the 3D x2-axis flag, signed-safe level cap, zero-cost, full-child topology, and aggregate-load risks. A second correction auditor added multilevel MeshBlock parity and wide adaptive-level arithmetic. |
| Durable findings | Added `ROB-033` through `ROB-035` for restart topology admission, multilevel parity, and adaptive level-count overflow. |
| Resolution | Validate node payload paths before opening temporary payloads; validate canonical spacing, root level, active MeshBlock parity, and signed-safe level `<=30`; validate decoded leaf axes, duplicates, ancestors, complete child sets, canonical persisted Z order, positive finite costs, finite aggregate cost, neighbor balance, and post-load-balance rank maps before payload routing; compute adaptive maximum level in wide arithmetic. |
| Focused coverage | Added corruption regressions for generated payload path before open, spacing, MeshBlock extent, multilevel parity, adaptive level count, negative and excessive levels, active and inactive axis bounds, duplicate locations, permuted records, ancestor overlap, incomplete refinement, 2:1 imbalance, zero cost, infinite cost, finite-value aggregate overflow, and valid 3D restart with nonzero x2 location. Serial restart/layout harness set returned `39 passed`; MPI restart set returned `58 passed`. |
| Full local verification | Full serial IO plus CPU smoke of the GPU-selectable regression returned `236 passed`; full MPI IO returned `67 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; targeted py_compile and `git diff --check` returned no output. |
| Compatibility note | Restart wire bytes are unchanged. The additional checks reject malformed persisted metadata before unsafe arithmetic or silent payload relabeling. Previously recorded detached pre-edit producer comparisons remain byte-identical for valid shared and per-rank output. |
| Remaining RCP-01 action | Wait for the fresh final topology auditor, then commit the coherent checkpoint if accepted. |

### 2026-05-29: RCP-01 Final Admission Corrections

| Field | Record |
| --- | --- |
| Final topology-audit blockers | Serialized restart dimensions were validated internally but were not required to match the dimensional flags derived from the persisted parameter dump. Adaptive `num_levels` validation performed wide addition only after `GetOrAddInteger()` had already narrowed text through `atoi()`. |
| Durable findings | Added `ROB-036` and `ROB-037`. |
| Resolution | Require serialized mesh and MeshBlock active dimensions to agree with retained `multi_d` and `three_d` flags before index arithmetic. Read existing adaptive `num_levels` as text, parse with checked `std::stoll`, reject malformed text and out-of-range values before addition, and use `GetOrAddInteger()` only for the safe missing-default case. |
| Focused verification | Restart/layout harness set returned `40 passed`; MPI restart set returned `58 passed`; dimensional mismatch and oversized adaptive-level regressions pass. |
| Full local verification | Full serial IO plus CPU smoke of the GPU-selectable regression returned `237 passed`; full MPI IO returned `67 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; targeted py_compile and `git diff --check` returned no output. |
| Remaining RCP-01 action | Wait for a fresh admission-focused auditor, then commit the coherent checkpoint if accepted. |

### 2026-05-29: RCP-01 Inactive Coarse-Axis Correction

| Field | Record |
| --- | --- |
| Admission re-audit blocker | Inactive coarse MeshBlock axes were not required to retain canonical `0/0` bounds. |
| Durable finding | Added `ROB-038`. |
| Resolution | Reject nonzero inactive `cjs/cje` and `cks/cke` bounds while preserving the existing active-axis coarse-grid checks. |
| Focused verification | Restart/layout harness set returned `41 passed`; MPI restart set returned `58 passed`; inactive coarse-axis regression passes. |
| Full local verification | Full serial IO plus CPU smoke of the GPU-selectable regression returned `238 passed`; full MPI IO returned `67 passed`; repository style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; targeted py_compile and `git diff --check` returned no output. |
| Independent acceptance | Fresh single-issue re-auditor returned `ACCEPTED` after inspecting inactive coarse-axis bounds, active-axis checks, dimensional agreement, and checked adaptive-level parsing. |
| Status | `RCP-01` is locally closed. Commit the coherent checkpoint before starting `RCP-02`. |

### 2026-05-29: RCP-02 Pre-Edit MPI And Publication Audit

| Field | Record |
| --- | --- |
| Independent auditors | A diff-driven MPI auditor, a filesystem-publication auditor, and a sequence or namespace auditor completed read-only passes before implementation. Their recommendations agree on the narrow helper boundaries and compatibility constraints recorded below. |
| MPI inventory | The diff-driven auditor classified `83` direct `MPI_*` call sites across `13` touched source files: `45` branch-added and in scope, `9` inherited-but-touched requiring disposition, `18` inherited and outside the narrow checkpoint scope, and `11` fatal-shutdown calls that must remain simple. |
| MPI blockers | Add `ROB-039` and `ROB-040`. Node communicator metadata and teardown, timing reductions, PDF reductions, and node-restart publication contain unchecked calls. Rank-local exits in touched output and restart-read paths can strand peers before later collectives. |
| Publication inventory | `.bin`, `.cbin`, and shared or per-rank restart write public targets directly. Modern PDF and `sphslice` write deterministic temporary paths and rename. Node restart publishes payload temporaries before a manifest temporary. No audited path uses payload sync plus containing-directory sync, so the truthful guarantee is namespace atomicity rather than crash durability. |
| Publication blockers | Add `ROB-041` through `ROB-043`. Directory creation is unchecked, node-restart generations lack exclusive reservation, deterministic temporary files need explicit stale-file semantics, and generated path components are not validated consistently. |
| Sequence inventory | Fixed `%05d` buffers truncate at `100000` in `.bin`, `.cbin`, PDF, `sphslice`, restart, table, VTK, particle VTK, Cartesian-grid, and spherical-surface writers. Existing `ROB-017` and `ROB-024` remain open. Python PDF header inference additionally hardcodes exactly five digits. |
| Namespace inventory | Existing `ROB-021` remains open. Construction has no deterministic target registry. Keys must cover `.bin`, `.cbin`, modern and legacy PDF overlap, `sphslice`, and restart. PDF directories omit the first axis variable, so blocks with the same explicit `id` can collide even when histogram semantics differ. |
| Accepted helper boundary | Add a narrow shared MPI helper for best-effort error rendering, checked returns, and simple world abort. Add a narrow output-filesystem helper for checked directories, minimum-five-digit unbounded sequence rendering, checked counter advancement, temporary naming, rename publication, and safe generated path components. Keep communicator participation and publication ordering visible at call sites. |
| Compatibility boundary | Preserve legacy PDF bytes and direct append behavior. Preserve names below `100000`. Widen tokens at and above `100000`. Do not promise crash durability. Preserve inherited non-feature writer formats while migrating their sequence rendering away from fixed buffers. |
| Stop-gate status | Reconnaissance recorded. Implementation may begin, but `RCP-02` remains open until focused tests and two independent post-edit audits accept the result. |

### 2026-05-29: RCP-02 Publication Re-Audit Corrections

| Field | Record |
| --- | --- |
| Re-audit blockers | The first post-edit publication auditor found that node-restart failure cleanup could abort before remote node leaders finished unlinking owned payloads, `.bin` and `.cbin` header failures bypassed temporary cleanup, restart wrapper-level fatal exits bypassed attempt cleanup, unlink failures were silent, and producer-surface failure coverage was incomplete. |
| Resolution | Added process-local fatal cleanup hooks through `src/mpi_utils.hpp`; made node-restart failures that are already world-coordinated clean every rank, wait at a world barrier, and only then abort; retained manifest ownership until reservation removal completes; cleaned `.bin` and `.cbin` temporaries on header and wrapper-level failures; and report failed cleanup attempts to stderr. |
| Focused verification | Shared helper, namespace, and serial publication coverage returned `35 passed`; targeted node-restart MPI regressions returned `4 passed`; the diff-driven MPI source gate passed after dispositioning the coordinated cleanup barrier. Producer-level publication coverage then returned `29 passed`, including stale temporary replacement, failed rename cleanup for `.bin`, `.cbin`, modern PDF, `sphslice`, and restart, plus portable unwritable-directory rejection. |
| Namespace expansion | Construction-time namespace coverage returned `17 passed`, including cadence-independent `.bin`, `.cbin`, `sphslice`, negative-PVTK-`gid`, and PDF histogram-definition collision cases. Legacy PDF widened-name inference now includes `100001`. |
| External limit | Real cross-node failure-injection qualification remains an `RCP-09` scheduler-backed gate. Local evidence uses one physical host only. |
| Remaining action | Wait for the third independent filesystem-publication re-audit before locally closing `RCP-02`. |

### Reflection R-02: Keep Publication Guarantees Narrow And Observable

| Question | Record |
| --- | --- |
| Did re-audit change implementation direction? | Yes. Rename publication alone was insufficient for node restart. Expected multi-rank transaction failures now clean before abort; unexpected wrapper or MPI fatals still perform best-effort local cleanup. |
| Is crash durability promised? | No. The branch promises namespace-atomic rename publication and observable best-effort cleanup, not payload or directory sync durability. |
| Did the helper boundary remain narrow? | Yes. Shared hooks, checked MPI rendering, directories, sequence formatting, temporary naming, rename, and cleanup reporting are centralized. Communicator ordering remains visible in writers. |
| What remains externally unproven? | Cleanup ordering across at least two physical nodes and production filesystems. Preserve that as an explicit qualification gate rather than inferring it from one-host MPI. |

### 2026-05-29: RCP-03 Coarsened-Binary Contract And Evidence

| Field | Record |
| --- | --- |
| Contract correction | Superseded the earlier ghost-expanded producer allowance. The writer now supports uniform three-dimensional full-volume active-zone output only and rejects ghost zones, AMR, slices, and lower-dimensional meshes during construction. The Python reader deliberately remains broader for validated historical files. |
| Arithmetic hardening | Added checked `std::int64_t` Kokkos plane, factor, and iteration ranges; checked `std::size_t` allocation elements and bytes; 64-bit node sum and prefix helpers; and a direct two-rank MPI harness proving node-local totals and prefixes above `INT_MAX`. |
| Reader hardening | Validate moment count, suffix groups, shared roots, and variable-count divisibility before reading payload bytes. |
| Conversion evidence | The moment-enabled shared, rank-sharded, and node-sharded regression now inspects generated ATHDF variable names, `uov` shape, `uov` values, XDMF attributes, and HDF hyperslab references rather than only checking that files exist. |
| Focused verification | Checked-layout harness returned `8 passed`; direct node-64 MPI harness returned `1 passed`; shared/rank/node moment ATHDF/XDMF regression returned `1 passed`; widened PDF reader subset returned `6 passed`; `git diff --check` returned no output. |
| Remaining action | Complete rebuild, run the focused reader and MPI writer matrix, then ask the RCP-03 auditor for a final read-only acceptance pass. |

### Reflection R-03: Prefer A Narrow Producer And A Compatible Reader

| Question | Record |
| --- | --- |
| Should the writer expand to match every reader capability? | No. Historical reader breadth is useful for analysis but is not evidence for newly supported producer rows. |
| Did the Kokkos refactor stay scoped? | Yes. Checked arithmetic is isolated in `coarsened_binary_layout.hpp`; serialization and format ownership remain in `coarsened_binary.cpp`. |
| Did testing become more representative? | Yes. Boundary tests now distinguish valid values above `INT_MAX` from true overflow, exercise node-local 64-bit collectives, and inspect generated conversion content. |
| Should lower-dimensional, AMR, sliced, or ghost-zone output be promoted now? | No. Add rows later only with explicit format semantics and end-to-end evidence. |

### 2026-05-29: RCP-02 Secondary Cleanup Reporting Correction

| Field | Record |
| --- | --- |
| Fresh re-audit blockers | The third independent publication auditor accepted the coordinated node-restart transaction but found secondary raw `remove()` calls after reservation-close, PDF write or close, and `sphslice` write or close failures. Those branches could silently fail to discard an attempt-owned path. |
| Resolution | Route each secondary discard through `output_file_utils::DiscardOwnedPath()` so unlink failures are reported to stderr. Keep the normal node-restart reservation removal explicitly world-coordinated because it participates in transaction completion. |
| Focused verification | Shared helper, namespace, serial publication, and MPI-source coverage returned `58 passed`; direct helper coverage now includes a failed owned-path discard and returned `10 passed`; injected node-restart payload-write and post-publication cleanup returned `2 passed`. |
| External limit | Cross-node cleanup ordering, remote-node completion before abort, and production-filesystem rename or unlink behavior remain scheduler-backed `RCP-09` requirements. Namespace atomicity and observable best-effort cleanup are claimed locally; crash durability is not. |
| Remaining action | Wait for a fresh correction-only acceptance auditor before marking `RCP-02` locally closed. |

### 2026-05-29: RCP-04 Diagnostic And Writer-Budget Implementation

| Field | Record |
| --- | --- |
| Numerical policy | Clamp bounded coordinate projections. Define zero radial or cylindrical projections on singular geometric axes. Reject non-positive-density and non-finite generic fluid diagnostics at runtime. Preserve construction-time rejection for ambiguous generic ion-neutral diagnostics. |
| Ghost-zone boundary | Reject derived-array `ghost_zones=true` output during construction. Reject PDF ghost sampling explicitly. Keep `sphslice` native-state-backed until ghost-zone-safe derived interpolation exists. |
| Spherical-slice naming | Replace low-precision radius names with deterministic round-trip scientific tokens and permit distinct radii with the same `id`. |
| Writer allocation caps | Add positive per-output `max_writer_allocation_bytes` with a default of `536870912` bytes. Preflight PDF and spherical-slice retained and staging allocations before materialization. |
| Backend portability correction | Remove PDF `ScatterView` accumulation, use `Kokkos::atomic_add()`, and stage MPI reductions through host mirrors so PDF output does not assume GPU-aware MPI. Account for host mirrors in budget preflight. |
| Focused verification | Direct semantic harness passed; focused PDF and `sphslice` subset returned `45 passed`; MPI PDF and restart subset returned `3 passed`; serial canonical matrix returned `361 passed`; MPI canonical matrix returned `74 passed`; Python reader subset returned `159 passed`; style returned `2 passed`; fixture checks verified `27` artifacts. |
| External limit | Representative CUDA execution is still required under `RCP-09`; the local host has no CUDA toolchain or runtime. |
| Remaining action | Wait for a fresh documentation-to-code and backend-portability auditor before marking `RCP-04` locally closed. |

### Reflection R-04: Keep Diagnostics Explicit And Backend-Neutral

| Question | Record |
| --- | --- |
| Are these diagnostics general output fields or PDF-focused fields? | They are general derived output fields. PDFs are a major consumer but do not own the semantics. |
| Should any advertised diagnostic be removed? | No. Every advertised name has analytic coverage or an explicit construction-time rejection boundary. |
| Are ghost-zone restrictions consistent? | Yes. Derived ghost zones are rejected before loading because current kernels do not populate them. |
| Does precise radius naming remain usable? | Yes. The token is longer but deterministic, readable, and collision-resistant for adjacent representable radii. |
| Did GPU readiness change implementation direction? | Yes. PDF accumulation now uses backend-portable atomics and MPI reductions use host staging. Real CUDA execution remains a required external gate. |

### 2026-05-29: RCP-02 Through RCP-04 Local Closure

| Field | Record |
| --- | --- |
| RCP-02 acceptance | A fresh correction-only auditor accepted centralized reporting for attempt-owned cleanup, intentional normal-path reservation removal with coordinated fallback, checked directories, widened sequences, namespace preflight, MPI dispositions, wrapper path context, and local transaction cleanup. |
| RCP-03 acceptance | A fresh final auditor accepted the narrow uniform-three-dimensional active-zone `.cbin` producer, broader validated legacy reader, checked Kokkos arithmetic, 64-bit node accounting, docs agreement, and focused evidence. Added explicit adaptive-AMR and rank-sharded sliced negative rows after the audit. |
| RCP-04 acceptance | A fresh docs-to-code and backend-portability auditor accepted diagnostic singularity semantics, density failure policy, derived ghost-zone rejection, radius naming, writer caps, PDF atomics, host-staged reductions, and deferred docs agreement. |
| Canonical verification | Python readers returned `159 passed`; serial IO plus GPU-selectable CPU smoke returned `363 passed`; MPI IO returned `76 passed`; style returned `2 passed`; fixture checksum verification passed for all `27` artifacts; `py_compile`, collection (`439 tests`), and `git diff --check origin/main --` passed. The serial count includes the explicit adaptive-AMR negative row. The MPI count includes shared rank/node sliced `.cbin` rejection rows. |
| Frozen-guide check | `IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md` remains SHA-256 `cdf53351104c135f2d79f9ee2f36c2908001b00a31339db18eb1cf66bcae11ff`, `2222` lines, and `89840` bytes. |
| External boundary | Local closure is not production qualification. Real CUDA execution and scheduler-backed multi-node cleanup, empty-node, routing, restart, timing, scaling, and production-filesystem behavior remain open under `RCP-09`. |
| Status | `RCP-02`, `RCP-03`, and `RCP-04` are locally closed. Commit the coherent correctness snapshot before starting `RCP-05`. |

### 2026-05-29: Reopened RCP-04 Adaptive-Shape And Spherical-Slice Staging Correction

| Field | Record |
| --- | --- |
| Late independent findings | Added `ROB-044` and `ROB-045`. A retained derived-variable view was reallocated only on first use, so a derived PDF emitted after adaptive local-pack growth could write beyond its original meshblock extent. The `sphslice` constructor cap did not explicitly preflight later ownership, node-gather, sorting, metadata, or dense-serialization peaks. |
| Resolution | Added full-shape derived-view admission before every shared derived-field kernel. Split `sphslice` retained geometry from transient staging estimates; preflight global ownership scratch, shared dense storage, local sparse buffers, sort copies, node metadata and gather buffers, publication sorting, and dense serialization. After a correction audit rejected the first text model, replaced the parameter-dump counting invocation with structural sizing over the parsed hierarchy, streamed the actual serialized header directly to the temporary file, admitted path and token staging, and moved local and node sort checks ahead of retained-vector mutation. Metadata still uses an allocation-free counting stream. |
| Focused regressions | Added a GPU-selectable adaptive derived-PDF regression that emits through an AMR pack expansion and a serial `sphslice` serialization-cap regression with an oversized input dump. Both targeted regressions pass locally on the serial CPU build. |
| Scope boundary | The adaptive regression is CPU-executed locally and remains part of the CUDA qualification row under `RCP-09`. `sphslice` still rejects derived interpolation until ghost-zone-safe sampling is implemented. |
| Status | `RCP-04` is reopened pending focused matrices and a fresh independent acceptance audit. |

### 2026-05-29: RCP-04 Direct-Streaming Re-Audit Correction

| Field | Record |
| --- | --- |
| Fresh audit blockers | The first late-correction auditor found that `ParameterDump()` itself allocated padded strings during the sizing pass, materialized stringstream buffers and `.str()` copies overlapped without being counted twice, and local/node sorting admission occurred after retained vectors had already been populated or swapped. |
| Resolution | Size parameter dumps structurally without calling the serializer, stream metadata and the actual parameter dump directly through a checked file-backed stream buffer, include path/token staging in the budget, and perform local/node sort admission before retained shard mutation. |
| Focused verification | All four serial/MPI ordinary and particle builds completed. Serial output, GPU-selectable CPU smoke, writer hardening, namespace, and MPI-source coverage returned `113 passed`; MPI formats, writer hardening, and node sharding returned `63 passed`; tracked-particle serial returned `5 passed`; tracked-particle MPI returned `3 passed`; `git diff --check` returned no output. |
| Status | Awaiting the correction auditor's fresh acceptance report. |

### 2026-05-29: RCP-08A Manifest-Scaling Preparation And RCP-09 External Plan

| Field | Record |
| --- | --- |
| Current protocol | Every rank validates the strict public manifest and replicated payload headers before direct local-span loading. Retain that correct path until production evidence exists. |
| Instrumentation | Added default-off `ATHENAK_RESTART_MANIFEST_TIMING=1` structured diagnostics for per-rank validation, startup parse, and routed local-block loading. Validation records include payload count, manifest bytes, replicated header bytes, and the per-rank replicated-header read estimate. |
| Preregistration | `IO_RESTART_MANIFEST_SCALING_QUALIFICATION.md` fixes required topology tiers, payload targets, one warm-up plus five measured resumes, percentile reporting, environment capture, empty/non-owning-node arrangement, and keep/follow-up/separate-branch thresholds. |
| External plan | `IO_EXTERNAL_IO_QUALIFICATION_PLAN.md` records CUDA and scheduler-backed multi-node evidence requirements, including rank-to-host mapping, reader equality, restart rollback, no `.assembled` files, production-filesystem behavior, and independent audits. |
| Local environment | `Tin-Drum` exposes `/opt/homebrew/bin/mpirun`, but no `nvidia-smi`, `nvcc`, `hipcc`, `sbatch`, `srun`, `qsub`, or `bsub`. External CUDA and physical multi-node execution remain unavailable locally. |
| Status | `RCP-08A` implementation awaits rebuild, focused regression, and benchmark-design audit. `RCP-08B` and `RCP-09` remain externally blocked until scheduler-backed evidence is supplied. |

### 2026-05-29: RCP-10 Local Process-Packaging Preparation

| Field | Record |
| --- | --- |
| Decision | Added `IO_FEATURE_BRANCH_PROCESS_ARTIFACT_DISPOSITION.md` and accepted `D-084` for the local packaging stage. Preserve the active frozen guide, decision log, audit ledger, compatibility contract, qualification plans, and deferred Pages bundle while external gates remain open. |
| Durable branch records | Retain `IO_FORMAT_COMPATIBILITY.md`, the scaling preregistration, the external qualification plan, the Pages bundle, and its staging helper. |
| Review-only records | Preserve large guides and implementation diaries for review now; migrate durable content into the compatibility contract, Pages documentation, and pull-request record before any later cleanup. |
| Scope rule | No process artifact may be removed, moved, or summarized silently. The frozen robustification guide stays at its recorded root path until checksum-based execution closes. |
| Status | Local packaging preparation complete. Full matrices, final lane audits, detached Pages validation, and external qualification evidence remain required before final closure. |

### 2026-05-29: RCP-05 Bounded Output-Registration Refactor And Tracked-Particle Repair

| Field | Record |
| --- | --- |
| Refactor scope | Extracted common variable and ID parsing, MeshBlock selection, slice parsing, and PDF parsing into private helpers while preserving the existing constructor chain. Added safe in-class defaults for touched output parameters, converted the branch-added `sphslice` owner to `std::unique_ptr`, and documented byte units, collective participation, node-communicator lifetime, and manifest post-load invariants. |
| Adjacent inherited defect | The duplicated parser incorrectly required `variable` for `file_type=trk`. The touched tracked writer also used a host-stack counter inside a device kernel, computed payload offsets without `sizeof(float)`, and did not prove dense unique tag ownership. |
| Tracked resolution | Exclude `trk` from variable and ID parsing; validate particle-module presence and requested count; use a device counter with explicit synchronization; validate global dense unique tags; compute checked byte offsets; preserve the native-endian legacy payload; and retain MPI ordering with checked collectives. |
| Focused verification | Ordinary parser-sensitive serial matrix returned `88 passed`; ordinary MPI format-hardening matrix returned `8 passed`; dedicated `PROBLEM=part_random` serial tracked suite returned `5 passed`; dedicated tracked MPI suite returned `3 passed`; the checked-MPI source disposition test returned `1 passed`; `git diff --check` returned no output. |
| Scope boundary | Retain the established format-construction chain and raw list ownership outside the branch-added spherical-slice owner. Do not redesign inherited output formats or tracked-particle wire bytes in this checkpoint. |
| Status | Implementation complete locally; independent behavior-preservation and scope audits required before `RCP-05` closure. |

### Reflection R-05: Keep The Refactor Bounded

| Question | Record |
| --- | --- |
| Did parser extraction preserve the established construction chain? | Yes. Common parsing moved into focused helpers while writer dispatch and inherited output ownership stayed in place. |
| Did adjacent tracked-particle repair change compatibility? | No. It restores correct dense-tag ownership and byte offsets while preserving six native-endian floats per tag slot. |
| Did scope audit find maintainability defects? | Yes. It caught an allocation-domain cast, a world-vs-node collective comment ambiguity, stale vector wording, and board drift. All were corrected and re-audited. |
| Should a broader output-factory rewrite happen here? | No. It remains a separate refactor candidate after every inherited format has behavior snapshots. |
| Status | `RCP-05` is locally closed after independent behavior acceptance and correction-only scope re-audit. |

### 2026-05-30: RCP-07 Detached Pages Staging Helper Implementation

| Field | Record |
| --- | --- |
| Preservation correction | Removed destructive complete overlays for `configuration.md`, `running.md`, `tools/visualization.md`, and `modules/outputs.md`. Replaced them with reconciliation-aware marker-wrapped fragments that edit only live IO-owned sections. Retained a hash-checked examples-index replacement and an absence-checked worked-example add. |
| Machine contract | Added `deferred_docs/gh-pages/io-output-formats-and-sharding/manifest.json` with the frozen Pages baseline, exact eight-file public allowlist, five protected blobs, seven target baseline blobs, expected add absence, source SHA256 values, unique anchors, replacement-section SHA256 values, and contradiction rules. |
| Helper behavior | Added Python-standard-library `scripts/stage_gh_pages_io_docs.py`. Strict mode requires a clean detached linked worktree, frozen target `HEAD`, matching local `origin/gh-pages`, protected and target blobs, source hashes, unique anchors, idempotent in-memory transforms, contradiction acceptance, atomic writes, and exact allowlist dirtiness. `--verify-staged` permits only intended dirtiness and reproves byte identity and idempotence. `--reviewed-drift` emits an external three-way packet and fingerprints the target before and after generation. |
| Focused verification | `PYTHONDONTWRITEBYTECODE=1 python3 -m flake8 scripts/stage_gh_pages_io_docs.py tst/test_suite/io/test_stage_gh_pages_io_docs_cpu.py` passed. `PYTHONDONTWRITEBYTECODE=1 python3 -m pytest -p no:cacheprovider -q tst/test_suite/io/test_stage_gh_pages_io_docs_cpu.py` returned `11 passed`. |
| Detached Pages preview | Created `/tmp/athenak-gh-pages-io-docs` detached at `4833aa9341e19861297e330ff02aabfd8001935c`. Strict stage and `--verify-staged` passed with exactly eight intended paths. Contradiction search found zero unsupported alternate-converter, stale modern-PDF-scale, or stale legacy-weight rejection matches. |
| Documentation gates | From detached `docs/`, `make clean html SPHINXOPTS="-W --keep-going"` passed and `make linkcheck SPHINXOPTS="-W --keep-going"` passed with an empty `build/linkcheck/output.txt`. Rendered browser inspection loaded all eight staged pages, found every expected staged section anchor, and found the new example link in the examples index and navigation. |
| Drift-packet exercise | Ran real `--reviewed-drift /tmp/athenak-gh-pages-io-docs-drift /tmp/athenak-gh-pages-io-docs`. It wrote base/current/proposed trees, diffs, source copies, and `report.json` outside the target while leaving target status unchanged. |
| Baseline protection | Before removing the preview, protected working blobs matched the recorded IDs for `docs/source/index.md`, `docs/source/modules/index.md`, `docs/Makefile`, `docs/requirements.txt`, and `.github/workflows/docs.yml`. Removed the generated site, packet, logs, and detached worktree. Local `origin/gh-pages` remained `4833aa9341e19861297e330ff02aabfd8001935c`. |
| Status | Local implementation and validation passed. Independent helper-safety, docs-to-code, and navigation/build auditors remain required before `RCP-07` checkpoint closure. Repeat staging against a refreshed detached Pages baseline after code merge before opening a separate Pages review. |

### Reflection R-07: Keep Pages Staging Bounded And Drift-Aware

| Question | Record |
| --- | --- |
| Does the helper fail safely when Pages moves? | Yes locally. Strict mode rejects baseline, target-blob, source-hash, section-hash, anchor, protected-blob, and dirty-target drift before writes. Reviewed-drift mode emits an external reconciliation packet without target writes. |
| Are insertion anchors stable enough? | Yes for the inspected baseline. Replacements use unique live section boundaries plus section SHA256 values. The one insertion uses the unique `## Shipped Input Families` anchor. Any ambiguity is a hard failure. |
| Does staged documentation describe only qualified behavior? | The public payload preserves explicit boundaries: no external CUDA qualification claim, no real multi-node qualification claim, native state-backed `sphslice` scalar fields and native groups, and narrow uniform 3D active-zone full-volume `cbin` production support. Independent docs-to-code audit remains pending. |
| Should any page remain manual because automation would be brittle? | No page in the eight-file allowlist needs manual staging on the inspected baseline. Drift is deliberately escalated to packet review rather than guessed reconciliation. |
| Is the later publication procedure clear? | Yes. `VALIDATION.md` separates detached preview, write-free drift reconciliation, refreshed post-code-merge restaging, strict Sphinx gates, independent audit, and a separate Pages review. |

### 2026-05-30: RCP-02 Fresh Current-Tree Publication Reconciliation

| Field | Record |
| --- | --- |
| Trigger | A late report from an older snapshot required a fresh read-only audit of the current dirty tree before relying on prior publication-cleanup closure. |
| Fresh acceptance | The correction auditor accepted predictable node-restart header, payload-write, close, manifest-publication, and reservation-removal failures as coordinated cleanup paths. Unexpected wrapper or MPI fatal exits retain explicitly best-effort process-local cleanup rather than being mislabeled as coordinated recovery. |
| Source-audit correction | The current MPI source inventory unions tracked diffs with untracked `src/**` files. The accepted inventory contained `69` `(path, symbol)` rows and `112` direct MPI calls, including the branch-added `mpi_utils.hpp`. |
| Focused evidence | Serial and MPI Debug rebuilds passed. Shared-helper plus MPI-source audit returned `11 passed`; serial publication, modern PDF/`sphslice`, and reduced-cap coverage returned `31 passed`; node-restart reservation collision and injected cleanup stages returned `3 passed`; scoped `git diff --check` returned no output. |
| External boundary | Cross-node cleanup ordering, remote-node unlink completion before abort, and production parallel-filesystem rename/unlink behavior remain scheduler-backed `RCP-09` requirements. |
| Status | `RCP-02` remains locally closed. |

### 2026-05-30: RCP-04 Repeated Spherical-Slice Retention Correction

| Field | Record |
| --- | --- |
| Fresh audit blocker | A correction auditor reproduced a repeated shared-`sphslice` case where the prior dense `outarray` remained live during the next ownership-validation scratch allocation but was omitted from the admission check. With `ntheta=128`, `nphi=128`, and a `1900000`-byte cap, the unmodeled live peak was `1966080` bytes. |
| Resolution | Reserve the persistent ownership vector at construction. Before adaptive rebuild or retained sparse-vector mutation, include any retained dense and sparse writer buffers in the ownership-validation admission. After publication, release write-only dense and sparse staging buffers on every participating rank so the next cycle begins without stale staging. |
| Regression | Added a repeated shared-output reduced-cap test with `dcycle=1`, `ntheta=128`, `nphi=128`, and `max_writer_allocation_bytes=1900000`. It emits two spherical-slice files and leaves no temporary artifact. The focused cap subset returned `4 passed`. |
| Status | `RCP-04` remains reopened pending a fresh correction-only auditor and the broader format matrix. |

### 2026-05-30: RCP-06 Python Reader Audit Reopen

| Field | Record |
| --- | --- |
| Independent blockers | The API auditor found package-form imports such as `import vis.python.bin_convert` broken by top-level-only helper imports. The same audit demonstrated that `.bin` and `.cbin` variable-list tokenization could materialize unbudgeted metadata under a one-byte live limit. |
| Required correction | Preserve package-form and standalone/top-level imports. Preflight binary variable token expansion and retain metadata accounting before split and dictionary materialization. Add empty-shard reduced-budget regressions. |
| Status | `RCP-06` is reopened. A bounded correction worker owns the Python-reader files and tests; independent re-audit is required. |

### 2026-05-30: RCP-08A Benchmark-Contract Audit Reopen

| Field | Record |
| --- | --- |
| Independent blockers | The benchmark-design auditor found that a logical comparison-byte estimate was mislabeled as observed filesystem amplification; the p95 decision statistic and external startup endpoint were underspecified; the empty-node arrangement was not reproducible from the named decks; the timing regression checked labels rather than the record contract; and default-off wording overstated unconditional timer construction and potentially fatal report-only arithmetic. |
| Required correction | Relabel the scoped estimate, freeze aggregation and endpoint rules, freeze reproducible deck overrides and placement, strengthen parser-level and synthetic multi-payload regressions, gate disabled probes, and make diagnostic reporting nonfatal or explicitly bounded. |
| Status | `RCP-08A` is reopened. A bounded correction worker owns instrumentation, tests, and the preregistration document; fresh benchmark-design audit is required before scheduler use. |

### 2026-05-30: RCP-10 Comment-Hygiene And Publication-Prose Classification

| Field | Record |
| --- | --- |
| Contradiction scan | Production code contains no `.assembled`, `StageNodeRestart`, `CopyFileRange`, or `bin_convert_new` implementation path. Remaining matches are explicit negative tests, compatibility boundaries, deferred documentation, and qualification-plan checks. |
| Marker classification | Remaining `TODO`, `FIXME`, and `DBF` matches in `src/outputs/vtk_mesh.cpp`, `src/outputs/basetype_output.cpp`, and the top-level instructional comment in `src/outputs/outputs.hpp` are inherited from `origin/main`. The touched `src/main.cpp` restart flag comment was normalized during the RCP-08A correction. |
| Compatibility prose | Updated `IO_FORMAT_COMPATIBILITY.md` to describe the deterministic detached-worktree helper, exact allowlist, bounded anchors, hashes, contradiction checks, and write-free drift packet instead of the superseded manual overlay/include workflow. |
| Status | Local hygiene classification is recorded. Final whole-branch audits remain required. |

### 2026-05-30: RCP-04 Repeated-Retention Correction Closure

| Field | Record |
| --- | --- |
| Independent acceptance | A fresh correction-only auditor accepted ownership admission before adaptive rebuild and sparse mutation, retained dense/sparse accounting, bounded persistent ownership storage, admission-before-allocation in shared/rank/node and serialization phases, direct header streaming, and release of write-only buffers on publishing and non-publishing ranks. |
| Focused evidence | Serial format matrix returned `71 passed`; MPI format and writer-hardening matrix returned `8 passed`; a temporary two-rank repeated-output smoke emitted shared `2`, rank `4`, and node `2` spherical-slice files with no temporary artifacts; scoped and whole-worktree `git diff --check` returned no output. |
| External boundary | The local node smoke uses one physical machine. Representative CUDA and true multi-node qualification remain `RCP-09` gates. |
| Status | `RCP-04` is locally closed again. |

### 2026-05-30: RCP-06 Shared Python Reader Utility Closure

| Field | Record |
| --- | --- |
| Decision | Accepted `D-080`: extract private shared reader-limit, checked-arithmetic, token, and CLI helpers while retaining `bin_convert.py` as the one supported binary converter. |
| Corrected blockers | Added package-relative imports with direct-script fallbacks. Added `.bin` and `.cbin` variable-list tokenization preflight before split, retained variable-name and dictionary accounting before construction, and propagation through shard assembly and ATHDF conversion. |
| Evidence | Focused reader module returned `181 passed`. Package imports, top-level imports, direct scripts, all `16` keyword-only public signatures, `py_compile`, targeted `flake8`, and `git diff --check` passed. An independent API and memory-budget auditor accepted the corrected current tree. |
| Scope boundary | `bin_convert_new.py` remains absent and unsupported. Reader limits remain keyword-only overrides and CLI flags rather than a new public configuration module. |
| Status | `RCP-06` is locally closed. |

### 2026-05-30: RCP-07 Documentation-To-Code Inflection

| Field | Record |
| --- | --- |
| Independent blockers | The docs-to-code auditor found scalar-only `sphslice` prose inconsistent with the implemented multi-variable payload, omitted public `sphslice` 3D/interior-radius construction constraints, and omitted per-axis `cbin` divisibility requirements. |
| Chosen disposition | Accepted `D-098`: document and test native state-backed spherical-slice scalar fields and native multi-field groups rather than artificially reject groups. Document that the origin-centered spherical surface must fit inside a 3D domain, with a positive radius strictly interior to every domain face. Add per-axis `cbin` divisibility wherever the supported row is advertised. |
| Evidence | Added producer-reader coverage for `variable=hydro_w`, which returns `dens`, `velx`, `vely`, `velz`, and `eint`; added negative 2D and boundary-radius construction rows; focused new regressions returned `3 passed`; broader serial format plus cbin-layout matrix returned `79 passed`. |
| Status | RCP-07 documentation alignment remains open pending helper-safety correction, detached restaging, and fresh independent audits. |

### 2026-05-30: RCP-07 Fail-Closed Pages Helper Correction

| Field | Record |
| --- | --- |
| Independent blockers | Two helper-safety auditors reproduced hidden `assume-unchanged` and `skip-worktree` mutations, accepted intent-to-add state, accepted allowlisted symlink or mode substitutions, and omission of target-root metadata from reviewed-drift fingerprints. |
| Decision | Accepted `D-099`: compare the full index to `HEAD`, reject noncanonical index flags before trusting porcelain status, require regular-file type and expected mode for every payload target, and fingerprint the worktree root as well as descendants. |
| Regressions | Added intent-to-add, both hidden-index-flag variants, allowlisted symlink, allowlisted mode, and root-metadata mutation rows. The focused no-cache helper suite returned `22 passed`; `py_compile`, targeted `flake8`, and `git diff --check` passed. |
| Detached verification | Recreated `/tmp/athenak-gh-pages-io-docs` detached at `4833aa9341e19861297e330ff02aabfd8001935c`. Strict stage, `--verify-staged --run-builds`, warnings-as-errors HTML, linkcheck, and external reviewed-drift packet generation passed with exactly eight intended payload paths. |
| Status | Helper-safety correction is implemented. Fresh independent helper-safety and navigation/build acceptance remain required before local `RCP-07` closure. |

### 2026-05-30: RCP-07 Local Closure

| Field | Record |
| --- | --- |
| Independent helper-safety acceptance | A fresh adversarial auditor accepted full-index equality, intent-to-add rejection, hidden-index-flag rejection, exact allowlist proof, regular-file and mode enforcement, reviewed-drift destination exclusion, and root-plus-descendant fingerprinting. Disposable probes and the focused suite exercised the exploit classes. |
| Independent docs acceptance | A fresh docs-to-code auditor accepted the implemented `sphslice` surface contract, native group support, derived rejection, supported `cbin` divisibility wording, source hashes, and exact staging allowlist. An asymmetric 3D spherical-slice smoke accepted a shifted containing domain and rejected inadequate face clearance. |
| Independent navigation acceptance | A separate detached-worktree auditor staged the exact eight-file payload, ran warnings-as-errors HTML and linkcheck builds, verified no tracked deletions or protected-blob drift, and followed rendered navigation from the site index through the examples index to the IO example. |
| Residual boundary | The deferred bundle remains intentionally unpublished on this code branch. Refresh `origin/gh-pages`, repeat detached staging, and open a separate Pages review after code merge. External CUDA and physical multi-node qualification remain open. |
| Status | `RCP-07` is locally closed. |

### 2026-05-30: RCP-08A Benchmark-Contract Closure

| Field | Record |
| --- | --- |
| Correction | `NodeRestartManifest::Load()` now captures validation metrics without emitting records. `main` snapshots `startup_parse` before emitting either validation or startup records. Diagnostic estimates remain default-off and saturating. The preregistration document pairs `external_minimal_resume_wall_s` with every measured representative no-step launcher invocation. |
| Evidence | Focused timing regressions returned `3 passed`; full node-sharding MPI returned `58 passed`; style and `git diff --check` passed. A fresh benchmark-design correction auditor independently accepted the implementation, formula, endpoint, and external protocol. |
| External boundary | `RCP-08B` remains intentionally open until scheduler-backed representative measurements are collected. |
| Status | `RCP-08A` is locally closed. |

### 2026-05-30: RCP-10 Full Local Validation Matrix

| Field | Record |
| --- | --- |
| Builds | Rebuilt serial Debug `/tmp/athenak-io-robust-build`, MPI Debug `/tmp/athenak-io-robust-build-mpi`, tracked-particle serial `/tmp/athenak-io-robust-build-part`, and tracked-particle MPI `/tmp/athenak-io-robust-build-part-mpi`. |
| Full local IO matrices | Serial CPU IO returned `412 passed, 4 skipped`. MPI IO returned `79 passed, 3 skipped`. Collection returned `500 tests`. The skips are environment-selective rows, not unexpected failures. |
| Focused compatibility evidence | Frozen shared restart plus promoted serial examples returned `5 passed`. Frozen per-rank restart plus generated node-manifest resume returned `2 passed`. The promoted node-sharded reader/manifest example returned `1 passed`. Tracked-particle serial returned `5 passed`; tracked-particle MPI returned `3 passed`; GPU-selectable CPU smoke returned `2 passed`. |
| Fixture integrity | `shasum -a 256 -c SHA256SUMS` passed for all `27` frozen `origin/main` artifacts. |
| Static gates | Python `py_compile`, targeted `flake8`, repository C++ style, and `git diff --check origin/main --` passed. IO collection returned `500 tests`. |
| Contradiction and hygiene scan | No production implementation path matched `.assembled`, `StageNodeRestart`, `CopyFileRange`, or `bin_convert_new`. Remaining `TODO`, `FIXME`, and `DBF` markers are inherited baseline comments already classified above. |
| Deferred Pages | Fresh detached strict stage, exact `--verify-staged`, warnings-as-errors HTML, linkcheck, reviewed-drift packet generation, rendered inspection, helper-safety audit, docs-to-code audit, and navigation/build audit passed against local `origin/gh-pages=4833aa9341e19861297e330ff02aabfd8001935c`. |
| External probe | Host `Tin-Drum` provides `/opt/homebrew/bin/mpirun` but no `nvidia-smi`, `nvcc`, `hipcc`, `sbatch`, `srun`, `qsub`, or `bsub`. CUDA execution, physical multi-node MPI, scheduler-backed restart scaling, and deployment-filesystem qualification remain external gates. |
| Status | All locally executable matrix rows pass. `RCP-10` remains open pending final independent red-team lanes, deliberate commits, and the explicit external `RCP-08B` and `RCP-09` block. |

### 2026-05-30: RCP-10 Red-Team Reopen

| Field | Record |
| --- | --- |
| File-format blockers | The file-format lane reproduced inventory-free node `.bin` and `.cbin` relocation, intrinsic sparse `AKPDFV2` payload downgrade after header-declaration removal, inventory-free spherical-slice shards, and duplicate scalar metadata accepted with last-write-wins behavior. |
| Python API blockers | The Python lane reproduced scalable embedded athinput, PDF-header, and spherical-header text accepted under a low `max_live_bytes` budget; duplicate `.cbin` preheaders; and `make_athdf.py` corruption of interior `.bin` basename components. |
| Numerical blockers | The numerical lane showed that MHD Poynting flux can oppose gas motion while `edot_*_{out,in}` used gas velocity for classification. The same lane identified missing fail-closed handling for non-positive PDF mass density. |
| Packaging blockers | The test-evidence lane identified two MPI rejection subprocesses without explicit timeouts and a summary-only validation record below the frozen exact-command evidence standard. The scope lane identified a retained historical integration plan presented as governing, a decision-summary entry point that stopped at `D-067`, and a stale live checkpoint board. |
| External boundary | The test-evidence lane reconfirmed the intended external block: no local CUDA execution, scheduler-backed physical multi-node qualification, restart-manifest scaling measurement, or deployment-filesystem qualification. |
| Status | Local `RCP-10` closure was reopened. No earlier successful run is treated as final evidence until the corrected full matrix and fresh independent audits pass. |

### 2026-05-30: RCP-10 Bounded Red-Team Corrections

| Field | Record |
| --- | --- |
| Decisions | Accepted `D-100`, `D-101`, and `D-102`. Preserve historical inventory-free rank binary compatibility and bounded transitional unversioned sparse PDF support, but require fail-closed metadata for new node and spherical-slice layouts. Count scalable text materialization against live budgets. Partition energy channels by signed total flux. |
| Reader correction | Require node binary inventory metadata, require spherical-slice shard inventory metadata, require intrinsic dense and sparse V2 payloads to carry V2 header declarations, reject duplicate metadata, account retained decoded text, and replace only the final batch-converter suffix. |
| Runtime correction | Classify radial and vertical energy `_out` and `_in` channels by signed total energy flux after MHD Poynting contributions. Reject non-finite or non-positive conserved density for PDF `weight = mass`; retain finite signed `weight = variable`. |
| Test and package correction | Added explicit timeouts to the corrupt-manifest and conflicting rank/node MPI rejection rows. Marked `IO_FEATURE_BRANCH_INTEGRATION_PLAN.md` historical and superseded. Redirected the decision-summary entry point to controlling post-`D-067` decisions and the live ledger. |
| Focused evidence | Corrected reader module returned `197 passed`. Serial and MPI Debug rebuilds passed. Diagnostic harness, adversarial radial and vertical MHD regressions, and non-positive mass-density regressions returned `5 passed`. Broad Python `flake8`, `py_compile`, and `git diff --check origin/main --` passed before the final record update. Deferred Pages helper suite returned `22 passed`. |
| Status | Bounded corrections are implemented. Full local matrices, detached Pages restaging, fresh correction audits, deliberate commits, and clean-tree inspection remain required. |

### 2026-05-30: RCP-10 Corrected Local Matrix And Nine-File Pages Verification

| Field | Record |
| --- | --- |
| Decisions | Accepted `D-103`, `D-104`, and `D-105`. Bind node-PDF headers to binary payload leader ranks, reject non-positive restart segments before publication, supersede only the eight-file portion of `D-082`, and require immutable external qualification packets. |
| Builds | From the repository root, `cmake --build /tmp/athenak-io-robust-build -j 4`, `cmake --build /tmp/athenak-io-robust-build-mpi -j 4`, `cmake --build /tmp/athenak-io-robust-build-part -j 4`, and `cmake --build /tmp/athenak-io-robust-build-part-mpi -j 4` returned exit `0`. Retained logs: `/tmp/athenak-io-robust-logs/build-serial.log`, `build-mpi.log`, `build-part-serial.log`, and `build-part-mpi.log`. |
| Canonical serial matrix | From `/tmp/athenak-io-robust-build/src`, `PYTHONDONTWRITEBYTECODE=1 /Users/dbf75/.uv/envs/interactive/.venv/bin/python -m pytest -p no:cacheprovider -q /Users/dbf75/.codex/worktrees/e948/athenak-DF/tst/test_suite/io/*_cpu.py /Users/dbf75/.codex/worktrees/e948/athenak-DF/tst/test_suite/io/test_output_formats_gpu.py` returned exit `0`: `441 passed, 4 skipped`. Retained log: `/tmp/athenak-io-robust-logs/pytest-serial-canonical.log`. |
| Canonical MPI matrix | From `/tmp/athenak-io-robust-build-mpi/src`, `PYTHONDONTWRITEBYTECODE=1 /Users/dbf75/.uv/envs/interactive/.venv/bin/python -m pytest -p no:cacheprovider -q /Users/dbf75/.codex/worktrees/e948/athenak-DF/tst/test_suite/io/*_mpicpu.py` returned exit `0`: `79 passed, 3 skipped`. Retained log: `/tmp/athenak-io-robust-logs/pytest-mpi-canonical.log`. |
| Dedicated tracked-particle matrix | From `/tmp/athenak-io-robust-build-part/src`, the explicit-interpreter `test_tracked_particle_output_cpu.py` command returned exit `0`: `5 passed`. From `/tmp/athenak-io-robust-build-part-mpi/src`, the matching `test_tracked_particle_output_mpicpu.py` command returned exit `0`: `3 passed`. Retained logs: `/tmp/athenak-io-robust-logs/pytest-part-serial.log` and `pytest-part-mpi.log`. |
| Reader correction subset | From the repository root, `PYTHONDONTWRITEBYTECODE=1 /Users/dbf75/.uv/envs/interactive/.venv/bin/python -m pytest -p no:cacheprovider -q tst/test_suite/io/test_python_io_readers_cpu.py tst/test_suite/io/test_writer_hardening_cpu.py` returned exit `0`: `230 passed`. Retained log: `/tmp/athenak-io-robust-logs/python-readers.log`. |
| Collection | From the repository root, `PYTHONDONTWRITEBYTECODE=1 /Users/dbf75/.uv/envs/interactive/.venv/bin/python -m pytest -p no:cacheprovider --collect-only -q tst/test_suite/io` returned exit `0`: `527 tests collected`. Retained log: `/tmp/athenak-io-robust-logs/pytest-collect.log`. |
| Fixtures | From `tst/fixtures/io/origin_main_886dd2a1`, `shasum -a 256 -c SHA256SUMS` returned exit `0` for all `27` immutable artifacts. Retained log: `/tmp/athenak-io-robust-logs/fixture-checksums.log`. |
| Static gates | From `tst`, the explicit-interpreter `pytest -p no:cacheprovider -q test_suite/style` command returned exit `0`: `2 passed`. Broad explicit-interpreter `flake8`, `py_compile`, `git diff --check`, and the production contradiction search returned exit `0`. Retained logs: `/tmp/athenak-io-robust-logs/style.log`, `flake8.log`, `static-gates.log`, and `contradiction.log`. |
| Nine-file detached Pages stage | Recreated `/tmp/athenak-gh-pages-io-docs` detached at local `origin/gh-pages=4833aa9341e19861297e330ff02aabfd8001935c`. From the repository root, `python scripts/stage_gh_pages_io_docs.py /tmp/athenak-gh-pages-io-docs`, `python scripts/stage_gh_pages_io_docs.py --verify-staged --run-builds /tmp/athenak-gh-pages-io-docs`, and `python scripts/stage_gh_pages_io_docs.py --reviewed-drift /tmp/athenak-gh-pages-io-docs-drift /tmp/athenak-gh-pages-io-docs` returned exit `0`. The exact dirty path set is nine files, linkcheck output is empty, and the reviewed-drift packet remained outside the target. Retained logs: `/tmp/athenak-io-robust-logs-pages/stage.log`, `verify-build.log`, and `drift.log`. |
| Rendered Pages QA | Served the detached HTML preview locally and inspected rendered navigation with the in-app browser. Confirmed `Data I/O (13 registered formats)`, followed the rendered IO example link, verified promoted shard-reader and manifest-resume commands, verified `payload_rank` fail-closed prose, checked Outputs-module mass-weight wording, and preserved the home-page iframe. Retained record: `/tmp/athenak-io-robust-logs-pages/browser.log`. |
| External boundary | Host `Tin-Drum` still has `/opt/homebrew/bin/mpirun` but no CUDA/HIP toolchain, scheduler launcher, or second physical host. `RCP-08B` and `RCP-09` remain external gates governed by immutable evidence packets. |
| Status | Locally executable validation passed. Fresh correction audits, deliberate commits, final whole-branch audit, clean-tree inspection, and the explicit external gates remain required. |

### 2026-05-30: RCP-10 Late Admission, Restart, And Rendered-Pages Correction

| Field | Record |
| --- | --- |
| Independent late blockers | The Python API lane found binary metadata records and legacy PDF text rows that still decoded before low-memory rejection, public binary-reader dictionaries that leaked internal accounting, and a batch wrapper that bypassed canonical `convert_file`. The restart lane demonstrated an unused zero-block header-only payload accepted by the loader and a weak unversioned node-PDF relocation path. The numerical lane requested permanent signed-variable and non-finite-weight producer locks. The deferred-Pages lane found twice that marker comments inside Markdown tables split rendered HTML even though Sphinx and linkcheck passed: first in the Support Systems table, then in the Implementation Entry Points table. It also found that the retained PDF implementation row still said one or two dimensions. |
| Decisions | Accepted `D-106` and `D-107`. Reject non-positive restart payload counts; reject unversioned node PDFs while retaining shared/rank transitional compatibility; enforce text admission before decode; keep parser accounting private; route `make_athdf.py` through canonical conversion; wrap the full bounded Support Systems and Implementation Entry Points sections so markers remain outside table rows; refresh the PDF implementation entry. |
| Focused correction evidence | Reader and writer-hardening subset returned `236 passed`. Serial output-format producer module returned `79 passed`. Corrupted node-manifest subset returned `29 passed`. Deferred-Pages helper suite returned `23 passed`. Incremental serial, MPI, particle serial, and particle MPI builds returned exit `0`. |
| Corrected full local matrix | Serial canonical IO including the GPU-selectable CPU smoke returned `453 passed, 4 skipped`. MPI IO returned `80 passed, 3 skipped`. Collection returned `540 tests`. Reader plus writer-hardening subset returned `236 passed`. Tracked-particle serial and MPI rows returned `5 passed` and `3 passed`. Repository style returned `2 passed`; targeted `flake8`, `py_compile`, fixture checksums, and `git diff --check` returned exit `0`. The preceding `441`/`79`/`527` matrix row is retained as historical evidence and superseded by this corrected rerun. |
| Corrected detached Pages evidence | Recreated the detached preview from local `origin/gh-pages=4833aa9341e19861297e330ff02aabfd8001935c`. Strict stage, `--verify-staged --run-builds`, empty linkcheck, and reviewed-drift packet generation passed for the exact nine-file payload. Browser QA recorded all checks true, including Outputs and Boundary Values in one Support Systems table, both links, the 13-format count, one-through-four PDF implementation wording, mass-weight wording, example route and commands, `payload_rank` fail-closed prose, and the preserved home iframe. |
| External boundary | This workstation still cannot execute representative CUDA/HIP or scheduler-backed physical multi-node qualification. `RCP-08B` and `RCP-09` remain explicit external gates. |
| Status | Local correction matrix passed. Fresh independent correction audits, deliberate commits, final committed-tree audit, and clean-tree inspection remain required. |

### 2026-05-30: RCP-10 Superseding Admission And Rendered-Structure Closure

| Field | Record |
| --- | --- |
| Decisions | Accepted `D-108` and `D-109`. Enforce reader live-memory caps before textual materialization, stream embedded athinput records with retained-object accounting, reject non-finite PDF histogram payloads while preserving finite signed values, extend the batch wrapper to canonical `.bin`/`.cbin` options, enforce root-side positive restart payload totals, and state the HIP qualification boundary explicitly. |
| Builds | Incremental serial Debug `/tmp/athenak-io-robust-build`, MPI Debug `/tmp/athenak-io-robust-build-mpi`, tracked-particle serial `/tmp/athenak-io-robust-build-part`, and tracked-particle MPI `/tmp/athenak-io-robust-build-part-mpi` builds returned exit `0`. Retained logs: `/tmp/athenak-io-robust-logs-final/build-serial.log`, `build-mpi.log`, `build-part-serial.log`, and `build-part-mpi.log`. |
| Canonical local matrix | Serial canonical IO including GPU-selectable CPU smoke returned `467 passed, 4 skipped`. MPI IO returned `80 passed, 3 skipped`. Dedicated tracked-particle serial and MPI rows returned `5 passed` and `3 passed`. Collection returned `554 tests`. Reader plus writer-hardening subset returned `250 passed`. Repository style returned `2 passed`. |
| Static and fixtures | Broad `py_compile`, `flake8`, `git diff --check origin/main --`, production contradiction search, and all `27` frozen fixture SHA-256 checks passed. Retained logs live under `/tmp/athenak-io-robust-logs-final/`. |
| Corrected detached Pages evidence | Recreated `/tmp/athenak-gh-pages-io-docs` detached from local `origin/gh-pages=4833aa9341e19861297e330ff02aabfd8001935c`. Strict stage, exact nine-file verification, warnings-as-errors HTML, empty linkcheck, and reviewed-drift packet generation passed. In-app browser QA proved that Support Systems is one HTML table containing Outputs and Boundary Values with working links, Implementation Entry Points is one HTML table containing `src/outputs/pdf.cpp`, no stray PDF paragraph exists, the updated `.bin`/`.cbin` batch-wrapper row renders, the example route resolves, `payload_rank` prose renders, and the home iframe remains present. Retained log: `/tmp/athenak-io-robust-logs-final/browser.log`. |
| Independent Pages acceptance | A fresh auditor independently recreated the preview and accepted both repaired HTML tables, marker placement outside tables, exact allowlist, source hash, protected blobs, empty linkcheck, preserved iframe, and deferred-publication boundary. |
| Restart re-audit | A fresh restart/MPI auditor accepted current writer-loader symmetry, cleanup, routing, and one-node MPI evidence and requested one P3 source-test strengthening: assert both the per-rank zero-segment and aggregate per-node zero-payload guards before public manifest creation. The regression was strengthened accordingly. |
| External boundary | CUDA remains the frozen required device lane. HIP readiness requires a separate packet when a HIP deployment is intended. Scheduler-backed physical multi-node routing, deployment-filesystem cleanup behavior, and `RCP-08B` scaling measurements remain external. |
| Status | Superseding local matrix passed. Fresh reader, file-format, numerical, test-evidence, scope-packaging, and final whole-branch audit acceptance, deliberate commits, and committed-tree inspection remain required. |

### 2026-05-30: RCP-10 Python Container And Spherical Numerical Closure Correction

| Field | Record |
| --- | --- |
| Independent Python blockers | A fresh Python audit found that historical spherical-slice version-1 shared/rank fallback had been removed accidentally, post-discovery validation still rebuilt unaccounted shard identifier containers, direct binary validation and athdf-like slice classification created unaccounted tuple sets, legacy PDF reconstruction retained row arrays before stacking, `make_athdf.py` materialized an unbounded glob inventory, and high-cardinality metadata dictionaries remained undercharged. |
| Decisions | Accepted `D-111` and `D-112`. Preserve narrow historical spherical-slice fallback while keeping new explicit rank/node layouts strict. Stream sibling validation, charge sort peaks and transient tuple ownership, bound batch discovery, preallocate legacy PDF rows, charge per-record header objects, keep the documented direct-reader ndarray contract, and add independent numerical producer oracles. |
| Reader correction evidence | `PYTHONDONTWRITEBYTECODE=1 /Users/dbf75/.uv/envs/interactive/.venv/bin/python -m pytest -p no:cacheprovider -q tst/test_suite/io/test_python_io_readers_cpu.py tst/test_suite/io/test_writer_hardening_cpu.py` returned exit `0`: `270 passed`. The suite now includes actual high-cardinality modern PDF, legacy PDF, and spherical-slice headers; historical spherical-slice shared/rank fixtures; direct logical-owner and slice-classification preflight; bounded batch discovery; and public ndarray assertions. |
| Numerical correction evidence | From `/tmp/athenak-io-robust-build/src`, the focused spherical-slice serial subset returned `12 passed, 71 deselected`. From `/tmp/athenak-io-robust-build-mpi/src`, the strengthened MPI output-format module returned `1 passed`. The serial regressions force MeshBlock-face interpolation and compare post-redistribution spherical slices to an independent binary-snapshot oracle. |
| Test-failure triage | Stricter metadata accounting initially moved several low-budget tests to earlier rejection phases. Recalibrated each test against parser-entry, retained-header, or initial-inventory state as appropriate rather than weakening validation. A later serial matrix intentionally failed because reviewed documentation fragments had stale machine-manifest hashes; refreshed only the five reviewed source hashes, then updated the visualization-utility hash again after documenting the restored historical reader fallback. |
| Frozen guide | `IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md` remains SHA-256 `cdf53351104c135f2d79f9ee2f36c2908001b00a31339db18eb1cf66bcae11ff`, `2222` lines, and `89840` bytes. |
| Status | Local corrections are implemented. Canonical serial rerun, detached Pages restaging, fresh Python audit disposition, final lane audits, deliberate commits, committed-tree verification, and the external `RCP-08B`/`RCP-09` gates remain open. |

### Reflection R-10.1: Keep Compatibility Narrow And Make Admission Claims Literal

| Question | Record |
| --- | --- |
| Did strictness remove useful compatibility? | Yes once: historical spherical-slice version-1 fallback disappeared while the format version remained unchanged. Restored only the shared/rank inference path; new explicit rank/node layouts remain strict. |
| Is a text-byte multiplier sufficient for Python readers? | No. It bounds decoded text expansion but not many short retained dictionary and set records. Per-record object admission now supplements the text charge. |
| Should direct binary readers expose lists again? | No. The documented `mb_data[variable]` array contract is clearer, ordinary indexing remains intact, and preallocation avoids scalable row-list overhead. |
| Did numerical test direction change? | Yes. Shape and cross-layout comparisons were insufficient for interpolation. Added analytical cross-face and independent post-redistribution binary-backed oracles. |
| Is the branch ready for merge? | No. Local correction qualification and committed-tree audit remain active, and external CUDA plus scheduler-backed multi-node evidence is still unavailable on this workstation. |

### 2026-05-30: RCP-10 Commit-Point, Intrinsic-Header, And Coarse-Fine Correction

| Field | Record |
| --- | --- |
| Independent restart/MPI blockers | The restart lane found that cleanup still owned published manifests and payloads after public-manifest rename. A later reservation-removal or injected post-publication failure could invalidate an already visible checkpoint. It also requested a local one-rank resume that combines spans from more than one payload under forced-small-chunk reads. |
| Independent file-format blockers | The format lane found that explicit spherical-slice layouts were not mandatory, degenerate dimensions and non-positive radii reached payload processing, header-only PDF V2 sparse APIs accepted incomplete inventories, and historical spherical fallback evidence was described more strongly than the synthetic regression justified. |
| Independent numerical blockers | The numerical lane found finite `double` spherical samples that overflowed to non-finite serialized floats after the writer's finite-value check. It also found that the adaptive binary-backed oracle refined every block uniformly and therefore did not force coarse-fine interpolation stencils. |
| Independent evidence blockers | The evidence lane found that the external packet checksum rule was self-referential, scheduler instructions were not directly executable, per-rank logs and timeout policy were underspecified, and logical restart-validation pressure was standing in for the guide-required observed filesystem-read-amplification measurement. |
| Decisions | Accepted `D-113` through `D-116`. Treat public-manifest rename as the transaction commit point; require postcommit checkpoint preservation; tighten intrinsic header and header-only reader admission; document host-native new-payload byte order; require attributable observed filesystem-read measurements; make packet manifests non-self-referential; and add serialized-float plus mixed-level AMR spherical-slice regressions. |
| Focused correction evidence | Reader and writer-hardening suite returned `288 passed`. Focused spherical-slice producer subset returned `14 passed, 71 deselected`. Focused node-restart MPI subset returned `7 passed, 55 deselected`. Broad static gates and the canonical full local matrix remain to be rerun after deferred Pages hash refresh. |
| External boundary | The external plan now requires rank-separated scheduler logs, explicit launcher substitutions, timeout dispositions, CUDA spherical-slice execution, observed filesystem reads distinct from logical pressure, and packet checksum manifests that exclude themselves. This workstation cannot execute those external rows. |
| Status | Bounded corrections are implemented. Deferred Pages hashes, detached Pages restaging, canonical matrices, fresh settled-tree audits, deliberate commits, immutable local evidence, and committed-tree inspection remain required. |

### Reflection R-10.2: Keep Commit And Evidence Boundaries Explicit

| Question | Record |
| --- | --- |
| What is the node-restart transaction boundary? | Successful rename of the public manifest. Before that point rollback removes owned publication artifacts. After that point cleanup failures preserve the resumable checkpoint and return an error. |
| Did existing adaptive coverage force coarse-fine interpolation? | No. The earlier slope-driven case refined all blocks. The new location-driven case retains both level-0 and level-1 blocks and asserts that sampled owner-block ghost stencils cross levels. |
| Is logical validation pressure an observed filesystem measurement? | No. Keep it as required protocol telemetry and collect attributable observed filesystem-read bytes separately on the deployment platform. |
| Can historical spherical fallback be claimed as frozen compatibility evidence? | No. The regression is synthetic until a provenance-qualified predecessor artifact is frozen. |

### 2026-05-30: RCP-10 Settled Correction Matrix And Detached Pages Restage

| Field | Record |
| --- | --- |
| Frozen guide | `sha256sum IO_FEATURE_BRANCH_ROBUSTIFICATION_GUIDE.md` returned `cdf53351104c135f2d79f9ee2f36c2908001b00a31339db18eb1cf66bcae11ff`; `wc -lc` returned `2222` lines and `89840` bytes. |
| Builds | `cmake --build /tmp/athenak-io-robust-build -j 4`, `cmake --build /tmp/athenak-io-robust-build-mpi -j 4`, `cmake --build /tmp/athenak-io-robust-build-part -j 4`, and `cmake --build /tmp/athenak-io-robust-build-part-mpi -j 4` returned exit `0`. Retained logs: `/tmp/athenak-io-robust-logs-settled/build-*.log`. |
| Canonical serial matrix | From `/tmp/athenak-io-robust-build/src`, explicit-interpreter pytest over `tst/test_suite/io/*_cpu.py` plus `test_output_formats_gpu.py` returned exit `0`: `511 passed, 4 skipped`. Retained log: `/tmp/athenak-io-robust-logs-settled/pytest-serial-canonical.log`. |
| Canonical MPI matrix | From `/tmp/athenak-io-robust-build-mpi/src`, explicit-interpreter pytest over `tst/test_suite/io/*_mpicpu.py` returned exit `0`: `83 passed, 3 skipped`. Retained log: `/tmp/athenak-io-robust-logs-settled/pytest-mpi-canonical.log`. |
| Dedicated tracked-particle matrix | Explicit-interpreter serial and MPI tracked-particle modules returned exit `0`: `5 passed` and `3 passed`. Retained logs: `/tmp/athenak-io-robust-logs-settled/pytest-part-serial.log` and `pytest-part-mpi.log`. |
| Reader correction subset | Explicit-interpreter pytest over `test_python_io_readers_cpu.py` and `test_writer_hardening_cpu.py` returned exit `0`: `288 passed`. Retained log: `/tmp/athenak-io-robust-logs-settled/python-readers.log`. |
| Collection, style, fixtures, and static gates | Collection returned `601 tests`. Repository style returned `2 passed`. All `27` frozen fixture SHA-256 checks passed. Broad `py_compile`, `flake8`, `git diff --check`, and contradiction search returned exit `0`. Retained logs: `/tmp/athenak-io-robust-logs-settled/`. |
| Detached Pages stage | Recreated `/tmp/athenak-gh-pages-io-docs` detached at local `origin/gh-pages=4833aa9341e19861297e330ff02aabfd8001935c`. Strict stage, `--verify-staged --run-builds`, reviewed-drift generation, exact nine-file status, and empty `docs/build/linkcheck/output.txt` passed. Retained logs: `/tmp/athenak-io-robust-logs-settled-pages/`. |
| Rendered Pages QA | In-app browser QA confirmed 13 registered formats, one Support Systems table containing Outputs and Boundary Values, added binary/coarsened-binary/spherical-slice implementation entries, mass-weight prose, serialized-float overflow prose, native-endian limitation, ReaderLimits and both header-only APIs, promoted example readback plus manifest-resume commands, `payload_rank` prose, public restart commit-point prose, and one preserved home-page iframe. Retained record: `/tmp/athenak-io-robust-logs-settled-pages/browser-final.log`. |
| External boundary | CUDA, optional deployment-specific HIP, scheduler-backed physical multi-node routing, target-filesystem rename/unlink behavior, attributable observed filesystem-read amplification, and `RCP-08B` scaling remain external. |
| Status | Superseding local correction matrix passed. Fresh settled-tree audits, deliberate commits, immutable local packet assembly, committed-tree rerun, and external qualification remain open. |

### 2026-05-30: RCP-09 External Packet And Scheduler-Deck Reopen

| Field | Record |
| --- | --- |
| Independent finding | Fresh process-evidence audit rejected the packet checksum cycle, prose-only scheduler deck, stale live-board wording, and opening `RCP-08A` label in the scaling preregistration. |
| Checksum correction | `artifacts.sha256` now excludes itself and the packet index. The index records the inner digest; the ledger or archive record retains an outer index digest. |
| Scheduler correction | Added `scripts/run_external_io_qualification_slurm.sh` with fixed ED-1 and MR-1 mappings, hostfile-forced asymmetric MR-2 mapping, row/sample-separated logs, warm-up plus five measured resumes, measured-only restart telemetry, monotonic launcher timing, filesystem-accounting hooks, exit codes, and timeout dispositions. |
| Board correction | Refreshed RCP-06 and RCP-10 to retain the already-passed settled canonical matrix while leaving fresh audits, committed-tree rerun, and external gates open. |
| Status | Local syntax-check passed. Fresh mock lifecycle checks and process-evidence re-audit remain required before the scheduler deck is accepted for external use. |

### 2026-05-30: RCP-10 Header-Only Workflow And Sparse-Overflow Correction

| Field | Record |
| --- | --- |
| Independent finding | Fresh file-format audit found that public `read_pdf_header()` accepted malformed sparse AKPDFV2 distributions and path-ID mismatches, the staged Visualization example passed a payload path to the header API, and the public prose overstated sibling-family validation. Fresh numerical audit requested direct sparse sharded serialized-overflow coverage. |
| Decision | Accepted `D-118`. Keep public PDF header admission path-bound and fail-closed while preserving internal full-reader root-metadata compatibility. Narrow public prose to single-header declaration checks and use the companion `.header.pdf` path in the example. Add rank- and node-sharded MPI sparse serialized-overflow regressions. |
| Focused reader evidence | Reader plus writer-hardening subset returned exit `0`: `298 passed`, including public unbound-root rejection, internal full-reader root-metadata compatibility, opposite-family ID rejection, and rank-header `payload_rank` rejection. Staging-helper suite returned exit `0`: `25 passed`. Targeted `flake8`, `bash -n scripts/run_external_io_qualification_slurm.sh`, and `git diff --check` returned exit `0`. |
| Detached Pages evidence | Recreated the detached `origin/gh-pages` preview after fragment changes. Strict staging with builds, reviewed-drift packet generation, exact nine-file status, and empty linkcheck output passed. |
| Status | Sharded MPI correction, rendered Pages QA refresh, fresh targeted re-audits, canonical rerun, deliberate commits, immutable local packet assembly, committed-tree rerun, and external qualification remain open. |

### 2026-05-30: RCP-09 Terminal Scheduler-Deck Indexing Correction

| Field | Record |
| --- | --- |
| Independent finding | Fresh process-evidence re-audit found Slurm retry collisions, missing or false-pass terminal rows after accounting and `.assembled` failures, Bash-4-only MR-2 parsing, and stale syntax-check wording. |
| Decision | Accepted `D-119`. Add attempt-scoped packet paths, Slurm `%J` job-step logs, duplicate-key rejection, rank-map timeout indexing, terminal post-validation Athena indexing, and portable MR-2 hostfile parsing. |
| Local mock evidence | A `/tmp/athenak-slurm-runner-mock` lifecycle matrix passed under `/bin/bash` `3.2.57`: ordinary ED-1 generation, duplicate-attempt rejection, MR-2 hostfile parsing, unset-hook incompleteness, post-hook failure indexing, forbidden-`.assembled` indexing, and rank-map timeout indexing. |
| Static evidence | `bash -n scripts/run_external_io_qualification_slurm.sh`, scoped `flake8`, and `git diff --check` returned exit `0`. |
| Status | Fresh process-evidence re-audit remains required before external scheduler use. Physical external execution remains open. |

### 2026-05-30: RCP-09 Single-Process Timing And Restart-Sidecar Scan Correction

| Field | Record |
| --- | --- |
| Independent finding | Fresh process-evidence re-audit found that separate Python timer processes can produce non-comparable macOS monotonic epochs and that output-directory-only staging scans miss `<manifest>.assembled` and `<manifest>.assembled.tmp`. |
| Decision | Accepted `D-120`. Time each launcher within one Python process, retain its timing TSV, and scan packet, output, and manifest directories for both forbidden sidecar suffixes before terminal indexing. Add a checked-in pytest mock suite. |
| Local mock evidence | Expanded `/tmp/athenak-slurm-runner-mock-final` matrix passed under macOS `/usr/bin/python3` `3.9.6` and `/bin/bash` `3.2.57`: nonnegative timers, one warm-up plus five measured samples, duplicate rejection, MR-2 parsing, accounting failures, output staging, both manifest-sidecar suffixes, and rank-map timeout. |
| Checked-in regression evidence | `PYTHONDONTWRITEBYTECODE=1 /Users/dbf75/.uv/envs/interactive/.venv/bin/python -m pytest -p no:cacheprovider -q tst/test_suite/io/test_external_io_slurm_runner_cpu.py` returned exit `0`: `9 passed`. The suite includes direct negative-timer-artifact rejection for rank-map and Athena launches. |
| Canonical refresh | Serial IO plus GPU-selectable CPU smoke returned `532 passed, 4 skipped`; MPI IO remained `85 passed, 3 skipped`; collection returned `624 tests`; repository style returned `2 passed`; static gates passed. |
| Status | Fresh process-evidence re-audit remains required before external scheduler use. |

### 2026-05-30: RCP-09 Strict Retained-Timer Grammar Correction

| Field | Record |
| --- | --- |
| Independent finding | Fresh process-evidence re-audit rejected the prior scheduler deck because retained timing TSVs with extra rows still indexed as passed. It also found stale `.assembled`-only preregistration summaries and a premature re-audit citation. |
| Decision | Accepted `D-121`. Require exactly one retained two-column `monotonic_elapsed_ns` timing row with a nonnegative decimal integer; reject every other grammar terminally. Align all preregistration summaries with both forbidden sidecar suffixes and record the rejected audit truthfully. |
| Checked-in regression evidence | Expanded runner suite returned exit `0`: `19 passed`. It covers extra-row, missing-row, extra-column, wrong-key, malformed-value, and negative-value timing artifacts at both rank-map and Athena-launch boundaries. |
| Canonical refresh | Serial IO plus GPU-selectable CPU smoke returned `542 passed, 4 skipped`; MPI IO remained `85 passed, 3 skipped`; collection returned `634 tests`; repository style returned `2 passed`; fixture checksums and static gates passed. |
| Status | Local correction validation passed. A fresh independent process-evidence re-audit remains required before external scheduler use. |

### 2026-05-30: RCP-09 MR-2 Admission And Executable Packet-Finalization Correction

| Field | Record |
| --- | --- |
| Independent finding | Fresh process-evidence re-audit accepted `D-121` but rejected the deck because blank or whitespace-only MR-2 host-B lines still indexed as passing scheduler rows. It also identified operator-driven packet finalization as an avoidable packaging weakness. |
| Decisions | Accepted `D-122` and `D-123`. Reject blank MR-2 host lines before launch. Add a checked-in packet finalizer that constructs the acyclic inner manifest and outer archive record and removes packet write permissions. |
| Checked-in regression evidence | Runner plus packet-finalizer suites returned exit `0`: `26 passed`. They cover empty and whitespace-only MR-2 hostfile rejection; acyclic manifests; root-only metadata exclusion; nested retained filename handling; outer records outside the packet; removed write permissions; repeat-finalization rejection; and packet/archive symlink rejection. |
| Canonical refresh | Serial IO plus GPU-selectable CPU smoke returned `549 passed, 4 skipped`; MPI IO remained `85 passed, 3 skipped`; collection returned `641 tests`; repository style returned `2 passed`; fixture checksums and static gates passed. |
| Status | Local correction validation passed. A fresh independent process-evidence re-audit remains required before external scheduler use. |

### 2026-05-30: RCP-09 Fail-Closed Scheduler Evidence And Transactional Finalization Correction

| Field | Record |
| --- | --- |
| Independent finding | Fresh process-evidence re-audit rejected the prior deck because recursive-output scan failures could escape without terminal Athena rows, forbidden-staging scan failures could false-pass, requested MR-2 placement was not proven by retained observed rank maps, archive hard-link aliases could corrupt covered packet artifacts, and interrupted finalization was not retryable. |
| Decisions | Accepted `D-124` and `D-125`. Retain and validate observed rank maps, require canonical hostname tokens, terminally index scan failures with retained stderr, reject aliases and control-character paths, and make packet finalization checksum-verified, permission-verified, idempotent, and retryable before the outer archive claim. |
| Checked-in regression evidence | Runner plus packet-finalizer suites returned exit `0`: `45 passed`. They cover malformed hostfiles; missing and incorrect observed rank maps; recursive-inventory and staging-scan failures; strict timing TSV grammar; symlink and hard-link aliases; control-character root and artifact paths; Markdown-breaking packet roots; acyclic manifests; writable partial metadata recovery; permission-removal retry; archive-append retry; ordinary filename spaces; and idempotent repeats. |
| Canonical refresh | Serial IO plus GPU-selectable CPU smoke returned `568 passed, 4 skipped`; MPI IO returned `85 passed, 3 skipped`; reader subset returned `298 passed`; tracked-particle serial and MPI returned `5 passed` and `3 passed`; deferred-doc helper returned `25 passed`; collection returned `660 tests`; repository style returned `2 passed`; fixture checksums and static gates passed. |
| Status | Local correction validation passed. A fresh independent process-evidence re-audit remains required before external scheduler use. |

### 2026-05-30: RCP-09 Canonical Inventory And Structured Index Correction

| Field | Record |
| --- | --- |
| Independent finding | Fresh process-evidence re-audit accepted the D-124 topology and scan fixes but rejected the deck because retries could archive incomplete checksum inventories, failed permission scans could false-pass as empty, and preseeded or control-character runner indexes could be archived as checksum-valid malformed evidence. |
| Decision | Accepted `D-126`. Add one shared strict validator for packet-index and archive TSV records; admit packet aliases and serialized controls before scheduler work; regenerate and compare canonical inventories on every retry; reject stale helper-temp artifacts; and require successful permission scans before outer publication. |
| Checked-in regression evidence | Runner plus packet-finalizer suites returned exit `0`: `57 passed`. Added malformed preseeded index, index-symlink write-through, hook-control, stale helper-temp, incomplete retained-manifest, late retry mutation, failed permission-scan, malformed runner-index, header-only runner-index, mutated retained Markdown index, and malformed outer-archive rows to the existing fail-closed matrix. |
| Canonical refresh | Serial IO plus GPU-selectable CPU smoke returned `580 passed, 4 skipped`; MPI IO remained `85 passed, 3 skipped`; collection returned `672 tests`; repository style returned `2 passed`; fixture checksums and static gates passed. |
| Status | Local correction validation passed. A fresh independent process-evidence re-audit remains required before external scheduler use. |

### 2026-05-30: RCP-10 Superseding Precommit Matrix And Rendered Pages QA

| Field | Record |
| --- | --- |
| Canonical serial matrix | Explicit-interpreter serial IO plus GPU-selectable CPU smoke returned exit `0`: `532 passed, 4 skipped`. Retained log: `/tmp/athenak-io-robust-logs-precommit-final/pytest-serial-canonical.log`. |
| Canonical MPI matrix | Explicit-interpreter MPI IO returned exit `0`: `85 passed, 3 skipped`. Retained log: `/tmp/athenak-io-robust-logs-precommit-final/pytest-mpi-canonical.log`. |
| Readers, collection, and style | Reader plus writer-hardening subset returned `298 passed`; collection returned `624 tests`; checked-in scheduler-runner mock suite returned `9 passed`; repository style returned `2 passed`. |
| Static gates | Scoped `flake8`, broad reader/tool `py_compile`, `bash -n scripts/run_external_io_qualification_slurm.sh`, `git diff --check`, frozen fixture checksums, and frozen-guide verification passed. |
| Detached Pages | Fresh strict staging with warnings-as-errors builds, `--verify-staged --run-builds`, reviewed-drift generation, exact nine-file status, and empty linkcheck output passed. |
| Rendered Pages QA | In-app browser QA confirmed output commit-point, native-endian, serialized-overflow, spherical-slice, implementation-entry, and mass-weighting prose; Configuration and Running commit-point/native-endian prose; the 13-format module-index row and Boundary Values neighbor; Visualization header-path, declaration-scope, `ReaderLimits`, and both header APIs; promoted example readback, `--assemble-shards`, manifest resume, and spherical-slice prose; File Reference `payload_rank`, native-endian, commit-point, and declaration-scope prose; and one preserved home iframe. Retained record: `/tmp/athenak-io-robust-logs-settled-pages-final/browser-final.log`. |
| Status | Precommit local matrix passed. Fresh targeted re-audits, deliberate commits, immutable local packet assembly, committed-tree rerun, final whole-branch audit, and external qualification remain open. |

### 2026-05-30: RCP-09 Packet-Lifecycle And Alias-Admission Correction

| Field | Record |
| --- | --- |
| Independent finding | A fresh process-evidence re-audit rejected the `D-126` deck. It reproduced packet deck-copy write-through through a hard-link alias before scheduler work, observed rank-map symlink false passes, checksum-valid contradictory and impossible TSV histories, dangling archive-symlink write-through before rejection, and conflicting durable claims for one packet path. |
| Decision | Accepted `D-127`. Fail-closed scan packet symlinks, hard-link aliases, control-character paths, and scanner failures before deck copying and before every scheduler launch; require regular single-link observed rank-map artifacts; bind exit codes to dispositions; validate attempt-local lifecycle prefixes; reject dangling archive symlinks; and key outer archive claims by absolute packet path. |
| Checked-in regression evidence | Expanded runner plus packet-finalizer suites returned exit `0`: `70 passed`. Added preexisting deck hard-link, nested control-character path, packet-tree scanner-failure, rank-map symlink, rank-map hard-link, dangling archive-symlink, conflicting archive-claim, contradictory exit-status, negative exit-status, missing-prefix, and post-terminal continuation regressions to the existing fail-closed matrix. |
| Static evidence | `bash -n` over the runner and finalizer, Python bytecode compilation of the shared validator, and `git diff --check` returned exit `0`. |
| Status | Local correction validation passed. A fresh independent process-evidence re-audit remains required before external scheduler use. |

### 2026-05-30: RCP-09 Postlaunch Publication And Canonical-Identity Correction

| Field | Record |
| --- | --- |
| Independent finding | A fresh process-evidence re-audit rejected the `D-127` deck. It reproduced scheduler-created rank-map directory and inventory-target symlinks that reached topology publication before the next prelaunch scan, checksum-valid impossible accounting histories, and conflicting outer claims for one packet through lexically distinct absolute spellings. |
| Decision | Accepted `D-128`. Repeat packet-tree admission immediately after each scheduler return; publish topology inventories exclusively after rank-map directory admission; retain validation stderr through an external temporary followed by packet-local publication; constrain accounting dispositions to measured rows with truthful hook state; and require canonical outer packet identities. |
| Checked-in regression evidence | Expanded runner plus packet-finalizer suites returned exit `0`: `76 passed`. Added scheduler-created rank-map directory-symlink, scheduler-created topology-inventory symlink, noncanonical outer-packet path, unmeasured accounting-before, warm-up accounting-after, and measured pass with unset-hook regressions to the preceding matrix. |
| Static evidence | `bash -n` over the runner and finalizer, Python bytecode compilation of the shared validator, and `git diff --check` returned exit `0`. |
| Status | Local correction validation passed. A fresh independent process-evidence re-audit remains required before external scheduler use. |

### 2026-05-30: RCP-09 Timer And Accounting-Hook Admission Correction

| Field | Record |
| --- | --- |
| Reflection finding | The `D-128` rank-map publication fix exposed adjacent retained-evidence windows: ordinary timing TSV publication after scheduler return and packet mutation by site-filled accounting hooks between launcher-boundary scans. |
| Decision | Accepted `D-129`. Publish timing TSV files exclusively after subprocess return and repeat packet-tree admission after successful pre-launch hooks and after post-launch hooks before retaining inventories. |
| Checked-in regression evidence | Expanded runner plus packet-finalizer suites returned exit `0`: `78 passed`. Added scheduler-created timing-symlink and accounting-hook-created alias regressions to the preceding matrix. |
| Static evidence | `bash -n` over the runner and finalizer plus `git diff --check` returned exit `0`. |
| Status | Local correction validation passed. The active fresh process-evidence re-auditor was instructed to assess the latest live tree before external scheduler use. |

### 2026-05-30: RCP-09 Descriptor-Bound Timing And Index-Append Correction

| Field | Record |
| --- | --- |
| Independent finding | A fresh process-evidence re-audit rejected the `D-129` deck. It reproduced timing-file write-through through a scheduler-created `logs/` parent symlink despite final-component `O_EXCL`, and packet-index write-through when a failing pre-launch accounting hook replaced `packet-index.tsv` with an external symlink before returning. |
| Decision | Accepted `D-130`. Bind timing publication to preopened and postreturn-revalidated packet-root and log-directory descriptors. Move packet-index append into the shared validator with packet scan, index validation, parent descriptor/path identity checks, final-component no-follow open, regular single-link inode matching, append sync, and complete postappend validation. Scan after every hook return, including failures, before appending. |
| Checked-in regression evidence | Expanded runner plus packet-finalizer suites returned exit `0`: `81 passed`. Added scheduler-created timing-parent symlink, scheduler-created packet-root symlink, and failing-accounting-hook packet-index symlink regressions to the preceding matrix. |
| Static evidence | `bash -n` over the runner and finalizer, Python bytecode compilation of the shared validator, and `git diff --check` returned exit `0`. |
| Status | Local correction validation passed. A new independent process-evidence acceptance audit remains required before external scheduler use. |

### 2026-05-30: RCP-09 Stable Packet And Archive-Publication Correction

| Field | Record |
| --- | --- |
| Independent finding | A fresh process-evidence re-audit rejected the `D-130` deck. It reproduced successful-hook clean packet-root replacement, packet-index mutation before rejected-row validation, late read-only artifact archival outside `artifacts.sha256`, and outer-record write-through through a symlink substituted after shell preflight. |
| Decision | Accepted `D-131`. Pin and recheck packet-root identity throughout runner admission; validate complete candidate TSV histories before append and require the expected packet parent; pin packet-root and archive-directory identity during finalization; regenerate the canonical packet inventory after recursive write-permission removal; and publish outer archive rows through the shared descriptor-bound validator primitive with the expected archive parent and packet identities. |
| Checked-in regression evidence | Expanded runner plus packet-finalizer suites returned exit `0`: `85 passed`. Added clean packet-root replacement, rejected-candidate immutability, read-only late-artifact, and late archive-sink substitution regressions to the preceding matrix. |
| Static evidence | `bash -n` over the runner and finalizer, Python bytecode compilation of the shared validator, and `git diff --check` returned exit `0`. |
| Status | Local correction validation passed. A new independent process-evidence acceptance audit remains required before external scheduler use. |

### 2026-05-30: RCP-09 Post-Permission Alias-Admission Reflection

| Field | Record |
| --- | --- |
| Reflection finding | The `D-131` stable-publication correction exposed one adjacent boundary worth making explicit: recursive write-permission removal precedes the durable archive claim and can be wrapped or observed as a packet-topology mutation window. |
| Decision | Accepted `D-132`. Repeat the complete packet-tree safety scan after `chmod -R a-w` and at each later durable-publication boundary. Keep the pinned packet-root check and reject symlinks, hard-link aliases, control-character paths, and reserved helper-temporary paths before archival. |
| Checked-in regression evidence | Added direct fault injection that introduces a symlink during the permission-removal command and requires terminal rejection before any archive row is retained. |
| Status | Narrow correction validation and an independent latest-tree acceptance audit remain required before external scheduler use. |

### 2026-05-30: RCP-09 Root-Only Identity Admission Rejection

| Field | Record |
| --- | --- |
| Independent finding | A fresh read-only adversarial acceptance audit rejected the `D-131` staged tree. It reproduced clean replacement of packet-local index, logs, rank-map, and launch-output objects within a stable packet root; symlink insertion during canonical-manifest construction; a late artifact inserted through the postpublication `grep` window; and clean replacement of the archive sink inode. |
| Decision | Rejected root-only admission as `D-133`. The post-`chmod` alias scan in `D-132` remains necessary but is not sufficient. |
| Bound snapshot | Rejected staged tree `673c9a70c1d8de24c99d7e435f28f9410e2e0dd7`; `HEAD=f5c29fc564f17343ed594c1aa7d2442c99baf38c`; frozen guide SHA-256 `cdf53351104c135f2d79f9ee2f36c2908001b00a31339db18eb1cf66bcae11ff`. |
| Status | Rejected for external scheduler-deck use. Superseding correction required. |

### 2026-05-30: RCP-09 Child-Object And Archive-Sink Continuity Correction

| Field | Record |
| --- | --- |
| Decision | Accepted `D-134`. Pin and revalidate runner child evidence objects; require packet-index inode continuity before and after append; bracket manifest construction with full packet-tree safety admission; create and pin the archive sink through the shared descriptor-safe helper; remove the postpublication shell `grep`; and make descriptor-bound archive append the final substantive operation. |
| Trust boundary | Invoke finalization from a trusted environment without a concurrent same-owner packet or archive mutator. Owner write-bit removal provides an operator-visible immutable packet transition; it cannot prevent an owner from deliberately restoring permissions after finalization returns. |
| Checked-in regression evidence | Added direct clean-replacement probes for index, logs, inventory, rank-map directory, launch-output directory, and archive sink; manifest-scan symlink injection; and proof that no postpublication `grep` mutation window remains. |
| Status | Focused validation passed, but a fresh independent audit rejected best-effort publication. Superseded by `D-136`. |

### 2026-05-30: RCP-09 Best-Effort Publication Rejection

| Field | Record |
| --- | --- |
| Independent finding | A fresh read-only adversarial acceptance audit rejected the `D-134` staged tree. External-temporary `mv` publication could become cross-filesystem copy-and-remove and leave a non-retryable partial metadata pair. One unchecked `os.write()` in either append path could leave packet-local or shared archive TSV state malformed after rejection. |
| Decision | Rejected best-effort publication as `D-135`. Invocation-local child pins remain necessary but do not establish durable metadata or append retryability. |
| Bound snapshot | Rejected staged tree `e15fcd89b2ce5e6439ba6b37d72441b0d84f0b06`; `HEAD=f5c29fc564f17343ed594c1aa7d2442c99baf38c`; frozen guide SHA-256 `cdf53351104c135f2d79f9ee2f36c2908001b00a31339db18eb1cf66bcae11ff`. |
| Status | Rejected for external scheduler-deck use. Superseding durable-publication correction required. |

### 2026-05-30: RCP-09 Durable Metadata And Append Retry Correction

| Field | Record |
| --- | --- |
| Decision | Accepted `D-136`. Publish finalizer metadata through packet-local exclusive temporaries, complete writes, file sync, same-directory atomic rename, and packet-directory sync. Descriptor-safely reset exact root-level owned helper temporaries and writable inconsistent pairs. Append complete rows with checked write loops and truncate-plus-sync rollback on any write, sync, or postvalidation failure. Sync the archive parent after sink creation. |
| Operator lifecycle | From the first runner invocation through finalization, do not replace packet-local children between invocations and do not permit a concurrent same-owner packet or archive mutator. |
| Checked-in regression evidence | Added owned stale-helper-temporary recovery, writable inconsistent-pair recovery, packet-index short-write rollback with successful retry, archive-sink short-write rollback with successful retry, and packet-local metadata short-write completion regressions. |
| Status | Focused validation and a fresh independent latest-tree acceptance audit remain required before external scheduler use. |

### 2026-05-30: RCP-10 Durable-Publication Canonical Refresh

| Field | Record |
| --- | --- |
| Focused packet-tool evidence | Runner plus packet-finalizer fault-injection suites returned exit `0`: `97 passed`. The matrix now covers owned stale-helper-temporary cleanup, writable inconsistent metadata-pair recovery, checked completion under partial metadata writes, and truncate-plus-sync rollback with successful retry after packet-index and archive-row write failures. |
| Canonical local matrix | Explicit-interpreter serial IO plus GPU-selectable CPU smoke returned exit `0`: `620 passed, 4 skipped`. Explicit-interpreter MPI IO returned exit `0`: `85 passed, 3 skipped`. IO collection returned `712 tests`. |
| Supporting gates | Repository style returned `2 passed`; Bash syntax, targeted `py_compile`, targeted `flake8`, `git diff --check`, and frozen-guide verification returned exit `0`. The preceding four-build D-134 refresh, reader `298 passed`, tracked-particle `5 passed` and `3 passed`, deferred-doc helper `25 passed`, frozen-fixture checksum pass, fresh detached Pages strict stage and browser QA remain applicable because D-136 changes only process tooling and process documentation. |
| Status | Local correction validation passed. A fresh independent process-evidence acceptance audit remains required before external scheduler use. |

### 2026-05-30: RCP-09 Exception-Only Crash-Recovery Rejection

| Field | Record |
| --- | --- |
| Independent finding | A fresh read-only adversarial acceptance audit rejected the `D-136` staged tree. Process termination could retain poisoned packet-index or shared-archive append suffixes without invoking exception rollback. Directory-sync failure could leave visible metadata or an archive sink that retry accepted without re-syncing. Exact reserved-looking root files were silently deleted without ownership proof. |
| Decision | Rejected exception-only recovery and filename ownership as `D-137`. |
| Bound snapshot | Rejected staged tree `e276eb9855f729fd3c83d4f31d70a133a8e742a6`; `HEAD=f5c29fc564f17343ed594c1aa7d2442c99baf38c`; frozen guide SHA-256 `cdf53351104c135f2d79f9ee2f36c2908001b00a31339db18eb1cf66bcae11ff`. |
| Status | Rejected for external scheduler-deck use. Superseding crash-consistency correction required. |

### 2026-05-30: RCP-09 Whole-Target Atomic-Append Correction

| Field | Record |
| --- | --- |
| Decision | Accepted `D-138`. Publish complete packet-index and archive-row replacements through sibling temporaries, checked complete writes, file sync, same-directory atomic rename, and parent-directory sync. Refresh the runner's pinned packet-index identity after each append. Re-sync admitted retained metadata and archive parents on retry. Reject ambiguous retained packet-local reserved temporaries for explicit operator adjudication rather than silent deletion. |
| Checked-in regression evidence | Updated direct short-write probes to require preservation of the preceding complete packet-index and prior shared archive history; expanded finalizer reserved-temporary negatives; added runner rejection of ambiguous packet-index atomic temporaries before scheduler work; and injected packet-local metadata, archive-sink admission, packet-index replacement, and archive replacement parent-directory sync failures followed by explicit retry or re-sync. Runner plus packet-finalizer fault-injection suites returned exit `0`: `103 passed`. IO collection returned `718 tests`. Bash syntax, targeted `py_compile`, targeted `flake8`, `git diff --check`, frozen-guide verification, and serial-build example execution passed. |
| Status | Focused validation passed, but a fresh independent audit rejected archive-adjacent replacement-temporary admission. Superseded by `D-140`. |

### 2026-05-30: RCP-10 Crash-Consistency Correction Refresh

| Field | Record |
| --- | --- |
| Focused packet-tool evidence | Runner plus packet-finalizer fault-injection suites returned exit `0`: `103 passed`. The matrix now covers whole-target atomic replacement after short writes, parent-directory sync failure after visible replacement, retry admission of an existing archive sink after creation-sync failure, retained packet-local metadata re-sync, and fail-closed handling of ambiguous reserved temporaries. |
| Collection | Explicit-interpreter IO collection returned exit `0`: `718 tests collected`. |
| Supporting gates | Serial-build examples returned exit `0`: `4 passed`. Bash syntax, targeted `py_compile`, targeted `flake8`, `git diff --check`, and frozen-guide verification returned exit `0`. |
| Status | Focused local correction validation passed. Canonical D-138 matrix refresh and fresh independent acceptance audits remain required before commit. |

### 2026-05-30: RCP-09 Archive-Adjacent Temporary Admission Rejection

| Field | Record |
| --- | --- |
| Independent finding | A fresh read-only adversarial acceptance audit rejected the `D-138` staged tree. Process termination before archive replacement rename could retain `.archive.tsv.tmp.*` beside the shared archive, outside the packet tree. Retry admitted the canonical archive while the ambiguous sibling remained. |
| Decision | Rejected archive-adjacent replacement-temporary blind spots as `D-139`. |
| Bound snapshot | Rejected staged tree `ccc80754849e59cb6f45f6428d4bb371f8784b65`; `HEAD=f5c29fc564f17343ed594c1aa7d2442c99baf38c`; frozen guide SHA-256 `cdf53351104c135f2d79f9ee2f36c2908001b00a31339db18eb1cf66bcae11ff`. |
| Status | Rejected for external scheduler-deck use. Superseding archive-directory admission required. |

### 2026-05-30: RCP-09 Archive-Adjacent Temporary Admission Correction

| Field | Record |
| --- | --- |
| Decision | Accepted `D-140`. Before every archive-record open, scan the verified archive-parent descriptor for destination-scoped `.<archive-name>.tmp.*` siblings. Reject and preserve every match for explicit operator adjudication. |
| Checked-in regression evidence | Added exact and near archive-adjacent finalizer rejection, operator-cleanup recovery, actual subprocess kill before archive rename with retained-candidate adjudication, and actual subprocess kill after archive rename with idempotent retry regressions. Runner plus packet-finalizer fault-injection suites returned exit `0`: `108 passed`. IO collection returned `723 tests`. Bash syntax, targeted `py_compile`, targeted `flake8`, and `git diff --check` passed. |
| Status | Focused validation passed. A fresh independent latest-tree acceptance audit remains required before external scheduler use. |

### 2026-05-30: RCP-10 Archive-Admission Correction Refresh

| Field | Record |
| --- | --- |
| Focused packet-tool evidence | Runner plus packet-finalizer fault-injection suites returned exit `0`: `108 passed`. The matrix now includes exact and near archive-adjacent candidate rejection, operator-cleanup recovery, actual subprocess termination before archive rename, and actual subprocess termination after archive rename with idempotent retry. |
| Collection | Explicit-interpreter IO collection returned exit `0`: `723 tests collected`. |
| Canonical local matrix | Explicit-interpreter serial IO plus GPU-selectable CPU smoke returned exit `0`: `631 passed, 4 skipped`. Explicit-interpreter MPI IO returned exit `0`: `85 passed, 3 skipped`. Reader hardening returned `298 passed`; tracked-particle modules returned `5 passed` and `3 passed`; deferred Pages helper returned `25 passed`; repository style returned `2 passed`; all `27` frozen fixtures verified. |
| Supporting gates | Bash syntax, targeted `py_compile`, targeted `flake8`, `git diff --check`, frozen-guide verification, and detached Pages strict stage verification with warnings-as-errors HTML plus empty linkcheck returned exit `0`. |
| Status | Full local correction validation passed. Fresh independent acceptance audits remain required before commit. |

### 2026-05-30: RCP-10 Imported-Validator Descriptor-Hygiene Refinement

| Field | Record |
| --- | --- |
| Independent finding | The fresh `D-140` adversarial audit accepted external scheduler-deck use under the documented trust boundary and found no P1/P2 issues. It reported a P3 imported-module hygiene gap: malformed archive validation and expected-identity rejection closed descriptors only on success. |
| Decision | Accepted `D-141`. Close verified parent-directory, packet-index, archive-record, and archive-parent descriptors through `finally` throughout the shared validator. |
| Checked-in regression evidence | Added repeated parent-identity, directory-sync, existing-publication-target, publication-rollback-sync, aliased-metadata-reset, packet-index-identity, malformed-archive, and archive-identity rejection loops that preserve the baseline `/dev/fd` count. Focused descriptor-hygiene subset returned exit `0`: `8 passed`; runner plus packet-finalizer fault-injection suites returned exit `0`: `116 passed`. IO collection returned `731 tests`. |
| Canonical serial refresh | Explicit-interpreter serial IO plus GPU-selectable CPU smoke returned exit `0`: `639 passed, 4 skipped`. |
| Supporting gates | Explicit-interpreter MPI IO returned exit `0`: `85 passed, 3 skipped`. Repository style returned exit `0`: `2 passed`. Bash syntax, targeted `py_compile`, targeted `flake8`, `git diff --check`, frozen-guide verification, and detached Pages strict verification with warnings-as-errors HTML plus empty linkcheck returned exit `0`. |
| Status | Full local D-141 refinement validation passed. A fresh independent latest-tree acceptance audit remains required before commit. |

### 2026-05-30: RCP-10 Committed-Tree Verification And Immutable Local Packet

| Field | Record |
| --- | --- |
| Tested runtime/tooling snapshot | Commit `5ee873e2b099800a9e86f64f33b27901c0bb1510`; tree `1e474976fddd019d7d237b9b81828b99844621b3`. This ledger closeout is intentionally a later documentation-only commit. |
| Independent audit disposition | Fresh narrow and broad read-only auditors accepted settled tree `1e474976fddd019d7d237b9b81828b99844621b3` with no P1/P2/P3 findings. The narrow auditor independently reproduced the eight descriptor-hygiene regressions, the complete `116`-test scheduler/finalizer suite, and direct archive-adjacent and process-termination fault injection. The broad auditor independently reproduced the complete local qualification matrix in an isolated mirror. |
| Four builds | Serial, MPI, particle, and particle-plus-MPI builds all returned exit `0`: `[100%] Built target athena`. |
| Committed local matrix | Canonical serial IO plus GPU-selectable CPU smoke returned exit `0`: `639 passed, 4 skipped`. Canonical MPI IO returned exit `0`: `85 passed, 3 skipped`. Scheduler/finalizer tools returned exit `0`: `116 passed`. Reader hardening returned exit `0`: `298 passed`. Tracked-particle modules returned exit `0`: `5 passed` and `3 passed`. Deferred Pages helper returned exit `0`: `25 passed`. Repository style returned exit `0`: `2 passed`. IO collection returned `731 tests`. All `27` frozen fixtures verified. |
| Static gates | Bash syntax, targeted Python bytecode compilation, targeted `flake8`, `git diff --check`, and frozen-guide checksum, line-count, and byte-count verification returned exit `0`. |
| Deferred Pages verification | Fresh nine-file detached `origin/gh-pages` staging passed warnings-as-errors HTML and linkcheck. Rendered browser QA verified the preserved home-page iframe, module navigation placement, outputs semantics, reader API references, runnable-example references, and file-format contract text. No Pages worktree was published. |
| Immutable local packet | `/Users/dbf75/Work/Research/AthenaK/evidence/io-output-formats-and-sharding/2026-05-30-local-5ee873e2` |
| Inner artifact-manifest SHA-256 | `f310b319b2b8cfe2a5b7d7a99c63b9fb4603547d06d09badb1d708450236c3bf` |
| Outer packet-index SHA-256 | `058160efda3e7c29a0ad16f58a9a963bc9b8ee42ef6b493ac7eaecd6f683f119` |
| Archive record | Shared archive `/Users/dbf75/Work/Research/AthenaK/evidence/io-output-formats-and-sharding/io-qualification-archive.tsv` retains exactly one row after repeat finalization. Every retained artifact verifies against the inner manifest. Packet paths are read-only. |
| External boundary | This local packet demonstrates the committed tooling and packet lifecycle on one workstation only. It does not close `RCP-08B` or `RCP-09`: representative CUDA, optional deployment-specific HIP, scheduler-backed physical multi-node, attributable deployment-filesystem read amplification, and production-filesystem behavior remain external gates. |
| Status | Local robustification closeout passed for committed runtime/tooling snapshot `5ee873e2b099800a9e86f64f33b27901c0bb1510`. The branch must not be described as merge-ready or production-qualified until the explicit external gates pass and independent auditors accept their immutable evidence packets. |
