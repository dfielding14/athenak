# MHD-PIC Supported Toolchains

## Scope

This declaration defines the bounded toolchain scope for clean-launch work on
`paper_mhd_pic` and `extended_mhd_pic`. It is not a scientific qualification
record. A successful configure, build, or startup check does not qualify
Frontier portability, decomposition parity, GPU parity, or paper reproduction.

| Profile | Supported scope now | Qualification limit |
| --- | --- | --- |
| Host serial | Debug and Release configure, build, `-c`, parser suite, and bounded serial smokes | Host-only development evidence |
| Host OpenMP | Separate configure, build, and `-c` startup check | Not a Frontier CPU or threaded-runtime qualification |
| Host MPI | Separate configure, build, and `-c` startup check | Rank launch requires a suitable launcher; qualifying MPI runtime is deferred to controlled Frontier work |
| Frontier HIP/MPI | Declared AMD `gfx90a` candidate profile | Build/runtime/profile-selection evidence requires Frontier allocation and the submission unlock |

Keep the runtime identity separate from the build profile. Both
`paper_mhd_pic` and `extended_mhd_pic` use mass-normalized momentum state, but
extension controls must never be presented as paper-mode evidence.

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
and MPI:

```bash
module reset
module load PrgEnv-amd
module load amd
module load rocm
module load craype-accel-amd-gfx90a
export MPICH_GPU_SUPPORT_ENABLED=1

cmake -S . -B /tmp/athenak-pic-frontier-hip \
  -DCMAKE_BUILD_TYPE=Release \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_MPI=ON \
  -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_VEGA90A=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include" \
  -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64"
cmake --build /tmp/athenak-pic-frontier-hip -j 8
```

This is a candidate build recipe, not permission to submit. Do not invoke
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
