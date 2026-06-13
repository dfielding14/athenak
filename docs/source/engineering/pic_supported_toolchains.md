# MHD-PIC Supported Toolchains

## Scope

This declaration defines the bounded toolchain scope for clean-launch work on
the active `paper_mhd_pic_vl2_tsc` and `extended_mhd_pic` identities. The
historical `paper_mhd_pic` identity remains available for restart and
source-history compatibility only; do not use it for new publication work. A
successful configure, build, or startup check is not a scientific qualification
record and does not qualify Frontier portability, decomposition parity, GPU
parity, or paper reproduction.

| Profile | Supported scope now | Qualification limit |
| --- | --- | --- |
| Host serial | Debug and Release configure, build, `-c`, parser suite, and bounded serial smokes | Host-only development evidence |
| Host OpenMP | Separate configure, build, and `-c` startup check | Not a Frontier CPU or threaded-runtime qualification |
| Host MPI | Separate configure, build, and `-c` startup check | Rank launch requires a suitable launcher; qualifying MPI runtime is deferred to controlled Frontier work |
| Frontier HIP/MPI | Declared AMD `gfx90a` candidate profile | Build/runtime/profile-selection evidence requires Frontier allocation and the submission unlock |

Keep the runtime identity separate from the build profile. Both
`paper_mhd_pic_vl2_tsc` and `extended_mhd_pic` use mass-normalized momentum
state, but extension controls must never be presented as paper-mode evidence.

## Host Configure Commands

Serial:

```bash
cmake -S . -B /tmp/athenak-pic-host-debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DAthena_ENABLE_MPI=OFF \
  -DAthena_ENABLE_OPENMP=OFF
cmake --build /tmp/athenak-pic-host-debug -j 8
/tmp/athenak-pic-host-debug/src/athena -c
```

OpenMP:

```bash
cmake -S . -B /tmp/athenak-pic-host-openmp \
  -DCMAKE_BUILD_TYPE=Debug \
  -DAthena_ENABLE_MPI=OFF \
  -DAthena_ENABLE_OPENMP=ON \
  -DKokkos_ENABLE_OPENMP=ON
cmake --build /tmp/athenak-pic-host-openmp -j 8
OMP_PROC_BIND=false /tmp/athenak-pic-host-openmp/src/athena -c
```

MPI compile and startup identity:

```bash
cmake -S . -B /tmp/athenak-pic-host-mpi \
  -DCMAKE_BUILD_TYPE=Debug \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_MPI=ON
cmake --build /tmp/athenak-pic-host-mpi -j 8
/tmp/athenak-pic-host-mpi/src/athena -c
```

Do not infer MPI runtime qualification from the MPI `-c` output.

## Frontier HIP Candidate

The Frontier candidate profile uses the Cray `CC` wrapper, AMD `gfx90a`, HIP,
and MPI. Keep every bulk artifact under
`/lustre/orion/ast207/proj-shared/dfielding/PIC`; Kronos is not a build,
execution, continuation, or publication path. Before this workflow, install the
paired production control plane from reviewed clean tracked source files. The
production installer rejects modified or untracked control-plane source and
publishes its immutable tree relative to a pinned authorized-parent
descriptor. Source that installed environment profile, then invoke only the
installed build-profile writer:

```bash
set -euo pipefail

export PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
export CONTROL_PLANE_DIR="${PIC_ROOT}/control_plane/<reviewed-control-plane-digest>"
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
CONTROL_PLANE=("$PYTHON" -I "${CONTROL_PLANE_DIR}/run_control_plane.py")
export SRC_DIR=/ccs/home/dfielding/athenak-pic
source "${CONTROL_PLANE_DIR}/frontier_pic_environment.sh" || exit $?
export GIT_COMMIT_FULL="$(git -C "$SRC_DIR" rev-parse HEAD)"

"${CONTROL_PLANE[@]}" write_orion_build_profile.py \
  --source-root "$SRC_DIR" \
  --expected-git-commit "$GIT_COMMIT_FULL" \
  --profile-id hip-mpi-release-paper-pic

export GIT_COMMIT="${GIT_COMMIT_FULL:0:12}"
export CONFIG=hip-mpi-release-paper-pic
export BIN_DIR="${PIC_ROOT}/bin/${GIT_COMMIT}/${CONFIG}"
"${CONTROL_PLANE[@]}" create_clean_candidate_freeze.py \
  --source-root "$SRC_DIR" \
  --executable "${BIN_DIR}/athena" \
  --build-profile "${BIN_DIR}/build_profile.json" \
  --build-profile-id hip-mpi-release-paper-pic \
  --prepared-artifact-inventory \
    tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json
```

Before committing a candidate, regenerate
`tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json`
with `generate_prepared_pic_artifact_inventory.py` and review its diff. The
freeze re-derives the committed inventory from `source.tar`, so a stale deck or
analyzer checksum fails closed.

`write_orion_build_profile.py` accepts no operator-selected executable, output,
build-directory, log, command-file, or provenance-input paths. For
`hip-mpi-release-paper-pic`, it derives the Orion layout, requires fresh build,
bin, and log paths, materializes a detached local checkout plus recursive local
submodule checkouts from the authorized clean source closure, and invokes the
closed direct `/usr/bin/cmake` configure and build argv with
`/opt/cray/pe/craype/2.7.33/bin/CC`, `ROCM_PATH=/opt/rocm-6.2.4`, explicit
double precision through `Athena_SINGLE_PRECISION=OFF`, and a minimal subprocess
environment that excludes caller Git, Python, CMake, and loader overrides.
Single-precision shock runs are outside this release profile and require
separate qualification. The installed profile explicitly unloads the inactive
default `darshan-runtime` module so its site-Spack pkg-config path cannot drift into
the closed PIC compiler-wrapper environment. It captures the exact
argv as `build-invocations.json`, binds empty `git_status.preconfigure.txt` and
post-build `git_status.txt` captures from the fresh checkout, and records the
cache, module list, fixed toolchain description, recursive submodule status,
redacted runtime environment allowlist, retained `build-environment.json`, and
configure/build logs. It publishes the
schema-v3 `build_profile.json` and adjacent `profile_receipt.json` exclusively
under Orion.

The Frontier production environment sets
`MPICH_GPU_EAGER_REGISTER_HOST_MEM=0`, `MPICH_GPU_NO_ASYNC_COPY=1`, and
`MPICH_GPU_IPC_ENABLED=0`. Q043 teardown reproductions on the pinned CPE 24.11
stack showed intermittent HSA memory faults during process-exit cleanup with
asynchronous GPU copies and IPC handle caching enabled.

The clean-candidate freeze revalidates the profile and receipt and retains
read-only copies of all eleven provenance inputs. The closed argv, hashes, and
fresh detached checkout reduce stale or mixed-build mistakes; they do not
cryptographically prove that the compiler honored the recorded command or
establish compiler semantics. Reproducibility and qualification remain
separate gates. The freeze and later pre-submit snapshots use
descriptor-relative staging and publication below pinned authorized parents,
followed by a lexical-parent identity check. This closes ancestor-swap pathname
redirection during publication; it does not protect writable artifacts against
a malicious same-UID process.

This workflow creates a candidate build profile, not permission to submit. Do not invoke
`sbatch`, initialize the production ledger, or bypass the installed immutable
submission wrapper before the control-plane and storage unlock described in
[MHD-PIC Clean-Launch Runbook](pic_clean_launch_runbook.md).

## Sanitizers

The current Cray/AMD login toolchain exhibits an ASan startup issue: an
AddressSanitizer-instrumented executable aborts before useful runtime checks
with an ODR-violation diagnostic for duplicated `.str` globals. The same class
of startup failure remains visible when the build adds
`-fno-sanitize-address-use-odr-indicator`.

Do not disable ASan ODR detection and call the result an ASan pass. Until the
toolchain startup issue is resolved or baselined against an approved compiler,
use a UBSan-only host fallback:

```bash
cmake -S . -B /tmp/athenak-pic-host-ubsan \
  -DCMAKE_BUILD_TYPE=Debug \
  -DAthena_ENABLE_MPI=OFF \
  -DCMAKE_CXX_FLAGS="-fsanitize=undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=undefined"
cmake --build /tmp/athenak-pic-host-ubsan -j 8
UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1" \
  /tmp/athenak-pic-host-ubsan/src/athena -c
```

Record UBSan as UBSan-only evidence. It does not replace ASan memory-safety
coverage.

## Verified Against

- `CMakeLists.txt:35`
- `CMakeLists.txt:36`
- `CMakeLists.txt:50`
- `CMakeLists.txt:66`
- `configure.sh:21`
- `tst/publication/frontier_control_plane/frontier_pic_environment.sh:1`
- `tst/publication/frontier_control_plane/submit_frontier_job.sh:2`
- `tst/publication/frontier_control_plane/control_plane_common.py:187`
