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
| MPI and restart correctness | Node sharding, empty shards, manifests and resume numbering | Staged loader gap remains | Transactional publication and one-node MPI tests pass, but the node-restart loader still creates a temporary `.assembled` file; CP-03 must replace it with direct distributed reads. True multi-node empty-node qualification also remains required. |
| Tests and examples | Harness placement, fixtures, negative cases and executable usage | Complete locally | Final serial/MPI suites and fixture checksums pass; independent audits addressed and accepted corrections. |
| Deferred Pages integration | Candidate docs content, final API alignment and later Sphinx validation | Staged; revalidation required | Earlier deferred package built under a temporary `origin/gh-pages` worktree, but CP-01 found documentation corrections and the remote Pages tip advanced. Reapply and rebuild during CP-07; live branch remains unchanged. |

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
| Status | Preimplementation direction accepted; execute during CP-05. |

### 2026-05-29: CP-06 Local Qualification Resource Probe

| Field | Record |
| --- | --- |
| MPI launcher helper | Fixed `tst/scripts/utils/athena.py` so MPI execution uses separate `["mpiexec", "-n", ...]` argv entries rather than the invalid single element `"mpiexec -n"`. Python byte-compilation and `git diff --check` passed. |
| Local MPI topology | `mpirun -np 2 hostname | sort -u` reported one physical hostname, `Tin-Drum`. No scheduler launcher command was found locally. |
| Local GPU topology | The host has an Apple M4 Max GPU, but no CUDA or HIP compiler/runtime tool was found. The repository GPU suite configures `Kokkos_ENABLE_CUDA=On`, so this host cannot perform the required GPU qualification. |
| GPU-selectable regression | Strengthened `tst/test_suite/io/test_output_formats_gpu.py` to assert four-dimensional derived PDF axes and scalar weighting metadata. The test logic passed once against the CPU build; that smoke run is not GPU qualification. |
| Remaining external gates | Execute the GPU-selectable regression on a CUDA-capable build and qualify node output/restart on multiple physical nodes, including a genuinely empty or non-owning node. |
| Status | Local resource probe complete; external GPU and multi-node qualification remain required before merge readiness. |

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
| MPI-001 | Gotham per-node restart may expose a manifest before every payload write is complete. | Blocking | Add atomic payload publication, inventory validation, then atomic public-manifest publication. | Resolved on one node for resume plus traversal, absolute-path, incomplete-marker, and byte-count rejection. Additional inventory negatives remain for CP-03. |
| MPI-002 | Gotham performs global node communicator setup even when node sharding is unused and has ambiguous empty-shard behavior. | Blocking | Implement opt-in/lazy node-shard infrastructure and explicit valid-empty or skipped-shard reader contract. | Resolved for shared-node execution; scheduler-backed multi-node empty-shard run remains a production qualification item. |
| MPI-003 | Existing MPI-file truncation is racy and wrapper transfers narrow 64-bit sizes into MPI `int` counts. | Blocking | Implement synchronized one-rank truncation plus overflow-safe chunked byte IO, broadcast chunking, and forced-small-chunk tests during CP-02. | Implementation assigned; independent post-integration audit required. |
| TEST-002 | Existing sliced `cbin` producer emits a zero-width meshblock extent that canonical readback rejects even in shared mode. | Scope boundary | Do not claim sliced `cbin` node support in this branch; keep full-volume `.cbin` node equality testing and sliced `.bin` empty-owner testing. | Recorded for separate baseline repair. |
| TEST-003 | Repository ignore rules omitted twelve frozen `.bin` and `.rst` fixture payloads from the first preservation commit. | Blocking | Force-add only the immutable fixture payloads and amend the test-fixture commit. | Resolved in `fcc534fe`; frozen-snapshot re-audit passed. |
| RST-001 | Existing per-rank restart resume calls `IOWrapper::GetPosition()` without the per-rank mode flag, invoking MPI-IO on a standard file handle. | Blocking | Forward `single_file_per_rank` in `Mesh::BuildTreeFromRestart` and add an MPI per-rank resume regression. | Resolved; two-rank resume regression passes. |
| DOC-001 | Live Pages pages advertise converter/PDF behavior inconsistent with clean code baseline. | Blocking | Stage a deferred overlay that reconciles existing pages after code stabilizes; do not publish before code merge. | Earlier staged package passed warnings-as-errors validation; superseded by the CP-01 re-audit and requires CP-07 reapplication, rebuild, and re-audit. |
| DOC-002 | Planning records described node restart loading as native and used stale pre-transactional payload names. | Blocking | Describe current strict manifest validation plus transient `.assembled` staging truthfully; keep manifest path as the only supported restart entry point until CP-03. | Corrected during CP-01; final docs revalidation required after CP-03. |
