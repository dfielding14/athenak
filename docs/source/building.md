# Building AthenaK

## Prerequisites

- A compiler capable of C++17
- CMake 3.16 or newer; the bundled Kokkos dependency requires 3.16 even
  though AthenaK's top-level project declaration is older
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

Enable MPI through AthenaK's top-level option. The top-level CMake file finds
MPI and sets `Kokkos_ENABLE_MPI=ON` in the cache automatically:

```bash
cmake -S . -B build-mpi \
  -DAthena_ENABLE_MPI=ON
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
  -DKokkos_ARCH_VEGA90A=ON \
  -DCMAKE_CXX_COMPILER=${ROCM_PATH:-/opt/rocm}/bin/hipcc
cmake --build build-hip -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)
```

Set `ROCM_PATH` if ROCm is installed outside `/opt/rocm`, and point
`CMAKE_CXX_COMPILER` at its `hipcc` (or a ROCm clang compiler configured for
HIP). The bundled Kokkos configuration inspects the selected compiler when
enabling HIP support.

## Common CMake Options

| Option | Default | Notes |
|--------|---------|-------|
| `Athena_ENABLE_MPI` | `OFF` | Finds MPI and enables Kokkos MPI support when set |
| `Athena_ENABLE_OPENMP` | `OFF` | Finds OpenMP and enables Kokkos OpenMP support when set |
| `Athena_SINGLE_PRECISION` | `OFF` | Enable single precision build |
| `PROBLEM` | `built_in_pgens` | Select a problem generator (reconfigure after changing) |
| `CMAKE_BUILD_TYPE` | `Release` | Use `Debug`/`RelWithDebInfo` for development |

For a custom problem generator source at `src/pgen/my_problem.cpp`, configure with:

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
export CUDAToolkit_ROOT=/usr/local/cuda
export PATH="$CUDAToolkit_ROOT/bin:$PATH"
cmake -S . -B build-cuda \
  -DKokkos_ENABLE_CUDA=ON \
  -DCMAKE_CXX_COMPILER=${PWD}/kokkos/bin/nvcc_wrapper
```

### Reconfiguring an Existing Build

When switching options (e.g. enabling MPI), run `cmake -S . -B <build-dir> ...` again to regenerate cache entries before building.
