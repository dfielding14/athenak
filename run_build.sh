#!/bin/bash

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/23.12
module load cce/17.0.0
module load rocm/5.7.1
module load cmake cray-python emacs
module unload darshan-runtime
export MPICH_GPU_SUPPORT_ENABLED=1
export FI_MR_CACHE_MONITOR=disabled #memhooks
export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_CXI_RX_MATCH_MODE=software


athenak='/ccs/home/pkempski/athenak-PK'
build='/ccs/home/pkempski/athenak-PK/build'

cd ${athenak}
cmake -Bbuild -DAthena_ENABLE_MPI=ON -DKokkos_ARCH_ZEN3=ON -DKokkos_ARCH_VEGA90A=ON \
      -DKokkos_ENABLE_HIP=ON -DCMAKE_CXX_COMPILER=CC \
      -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64" \
      -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include  -munsafe-fp-atomics" \
      -DPROBLEM=part_static_turb
cd ${build}
make -j
