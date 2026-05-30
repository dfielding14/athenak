# IO Feature Audit Ledger

## Branch Identity And Rules

| Item | Value |
| --- | --- |
| Feature branch | `feature/io-output-formats-and-sharding` |
| Clean base | `origin/main` at `886dd2a1437e45a3a30b3eeebf2adfa838328f73` |
| Historical evidence | `origin/gotham-1.0` and explicitly listed source commits |
| Forbidden implementation source | `origin/feature/single-file-per-node-outputs` |

This ledger records independent audits, blocking findings, implementation
resolutions, and re-audit results required by `IO_FEATURE_BRANCH_GUIDE.md`.

## Audit Gates

| Gate | Scope | Status | Evidence Or Required Resolution |
| --- | --- | --- | --- |
| Historical-scope separation | Gotham IO changes versus unrelated source/physics/mesh/test drift | Complete | Independent audit accepted; only listed generic IO behaviors were reconstructed. |
| File-format compatibility | Headers, versions, names, frozen fixtures and reader boundaries | Locally complete | Frozen compatibility, modern writer/reader, and one-node node-layout evidence pass; real multi-node qualification remains recorded. |
| Python API consolidation | One canonical `bin_convert.py`, consumers, CLI and reader behavior | Complete locally | Canonical module/readers/examples pass tests; final documentation/tooling re-audit found no remaining actionable issue. |
| MPI and restart correctness | Node sharding, empty shards, manifests and resume numbering | Locally complete; external topology qualification remains | Native direct manifest loading, forced chunks, alias rejection, and one-node MPI tests pass. True multi-node routing with an empty or non-owning node remains required. |
| Tests and examples | Harness placement, fixtures, negative cases and executable usage | Complete locally | Final serial/MPI suites, style gate, fixture checksums, and executable examples pass; GPU and multi-node execution remain external gates. |
| Deferred Pages integration | Candidate docs content, final API alignment and later Sphinx validation | Validated in detached worktree; publication deferred | Deferred package builds against refreshed `origin/gh-pages` at `4833aa9341e19861297e330ff02aabfd8001935c`; live branch remains unchanged and publication waits for code merge. |

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
