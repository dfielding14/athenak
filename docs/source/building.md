# Building AthenaK

## Prerequisites

- C++17 compiler (GCC 10+, Clang 12+, or NVCC 11+ for CUDA)
- CMake 3.18 or newer
- Git (for submodules) and `make`/`ninja`
- Optional: MPI distribution, ROCm, or CUDA toolkit depending on your target

Before configuring, clone the repository and pull bundled dependencies:

```bash
git clone https://github.com/IAS-Astrophysics/athenak.git
cd athenak
git submodule update --init --recursive
```

## CPU Reference Build (Out-of-Source)

```bash
cmake -S . -B build-cpu -DCMAKE_BUILD_TYPE=Release
cmake --build build-cpu -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)
```

Artifacts live in `build-cpu/src/athena`. Re-run `cmake --build …` whenever you change code.

## MPI Build

MPI support requires both Athena and Kokkos to be configured with MPI enabled. Use a dedicated build directory:

```bash
cmake -S . -B build-mpi \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_MPI=ON
cmake --build build-mpi -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)
```

Run with `mpirun -np <ranks> ./build-mpi/src/athena -i inputs/hydro/sod.athinput`.

## GPU Builds

### NVIDIA (CUDA)

```bash
cmake -S . -B build-cuda \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ARCH_AMPERE80=ON \
  -DCMAKE_CXX_COMPILER=${PWD}/kokkos/bin/nvcc_wrapper
cmake --build build-cuda -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)
```

Choose the matching `Kokkos_ARCH_*` value (e.g. `VOLTA70`, `AMPERE86`) for your GPU.

### AMD (HIP/ROCm)

```bash
cmake -S . -B build-hip \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_VEGA90A=ON
cmake --build build-hip -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)
```

Ensure `HIP_PATH`/`ROCM_PATH` is set if the toolchain is installed outside the default location.

## Common CMake Options

| Option | Default | Notes |
|--------|---------|-------|
| `Athena_ENABLE_MPI` | `OFF` | Requires `Kokkos_ENABLE_MPI=ON` and an MPI toolchain |
| `Athena_ENABLE_OPENMP` | `OFF` | Also set `Kokkos_ENABLE_OPENMP=ON` |
| `Athena_SINGLE_PRECISION` | `OFF` | Enable single precision build |
| `PROBLEM` | `built_in_pgens` | Select a problem generator (reconfigure after changing) |
| `CMAKE_BUILD_TYPE` | `Release` | Use `Debug`/`RelWithDebInfo` for development |

Configure custom problem generators with:

```bash
cmake -S . -B build-problem -DPROBLEM=my_problem
cmake --build build-problem -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)
```

## Recommended Tools

- Use `cmake --build <dir> --target install` to stage binaries (optional).
- Add `-G Ninja` to `cmake -S . -B …` for faster incremental builds if Ninja is available.

## Troubleshooting

### Missing MPI

```bash
cmake -S . -B build-mpi -DAthena_ENABLE_MPI=ON -DMPI_CXX_COMPILER=mpicxx
```

### Missing CUDA Toolkit

```bash
export CUDA_HOME=/usr/local/cuda
cmake -S . -B build-cuda -DKokkos_ENABLE_CUDA=ON
```

### Reconfiguring an Existing Build

When switching options (e.g. enabling MPI), run `cmake -S . -B <build-dir> ...` again to regenerate cache entries before building.
