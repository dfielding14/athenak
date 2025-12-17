#!/bin/bash

module purge
module load modules/2.3-20240529 slurm cuda/12.3 openmpi/cuda-4.0.7
export LD_PRELOAD=/mnt/sw/fi/cephtweaks/lib/libcephtweaks.so
export CEPHTWEAKS_LAZYIO=1

athenak=/mnt/home/dfielding/athenak-df
build=$athenak/build

cmake \
  -D CMAKE_CXX_COMPILER=$athenak/kokkos/bin/nvcc_wrapper \
  -D Kokkos_ENABLE_CUDA=On \
  -D Kokkos_ARCH_HOPPER90=On \
  -D Athena_ENABLE_MPI=On \
  -D PROBLEM=TRML \
  -B $build \
  -S $athenak

cd $build
make -j
