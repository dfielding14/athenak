#!/bin/bash
module restore
module load cpe/24.07 PrgEnv-amd cray-mpich/8.1.30 craype-accel-amd-gfx90a amd/6.2.0 rocm/6.2.0
module load cmake cray-python emacs
module unload darshan-runtime
module load perftools-base perftools

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export FI_MR_CACHE_MONITOR=kdreg2
export MPICH_SMP_SINGLE_COPY_MODE=NONE

export FI_CXI_DEFAULT_CQ_SIZE=131072
export FI_CXI_DEFAULT_TX_SIZE=2048
export FI_CXI_RX_MATCH_MODE=hybrid

build="build"

cmake -B $build -DAthena_ENABLE_MPI=ON -DKokkos_ARCH_ZEN3=ON -DKokkos_ARCH_VEGA90A=ON \
   -DKokkos_ENABLE_HIP=ON -DCMAKE_CXX_COMPILER=CC \
   -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64" \
   -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include" \
   -DPROBLEM=cgm_cooling_flow_amr_metals

cd $build
#make clean
make -j 64

rm src/athena_instrumented
pat_build -w -o src/athena_instrumented src/athena
