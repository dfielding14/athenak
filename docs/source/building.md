# Building AthenaK

## Requirements

- C++17 compiler (GCC 8+, Clang 10+, NVCC 11+)
- CMake 3.13+
- MPI (optional)

## Basic CPU Build

```bash
mkdir build
cd build
cmake ..
make -j8
```

## MPI Build

```bash
cmake -B build -DAthena_ENABLE_MPI=ON
cd build && make -j8
```

## GPU Build (NVIDIA)

```bash
cmake -B build \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ARCH_VOLTA70=ON \
  -DCMAKE_CXX_COMPILER=$KOKKOS_PATH/bin/nvcc_wrapper
cd build && make -j8
```

## GPU Build (AMD)

```bash
cmake -B build \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_VEGA906=ON
cd build && make -j8
```

## Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `Athena_ENABLE_MPI` | OFF | Enable MPI parallelization |
| `Athena_ENABLE_OPENMP` | OFF | Enable OpenMP |
| `Athena_SINGLE_PRECISION` | OFF | Use single precision |
| `CMAKE_BUILD_TYPE` | Release | Debug/Release/RelWithDebInfo |
| `PROBLEM` | built_in_pgens | Problem generator name |

## Custom Problem Generator

```bash
cmake -B build -DPROBLEM=my_problem
```

## Troubleshooting

### MPI Not Found
```bash
cmake -DCMAKE_CXX_COMPILER=mpicxx ..
```

### CUDA Not Found
```bash
export CUDA_HOME=/usr/local/cuda
cmake -DKokkos_ENABLE_CUDA=ON ..
```