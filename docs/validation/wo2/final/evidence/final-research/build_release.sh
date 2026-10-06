#!/bin/bash
set -euo pipefail
cd /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2
module restore
module load PrgEnv-cray cpe/25.09 cray-mpich/9.0.1 cce/20.0.0
module load craype-accel-amd-gfx90a rocm/6.4.2
module unload darshan-runtime
unset CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH
export INCLUDE_PATH_X86_64="$(CC --print-resource-dir)/include:${CC_X86_64}/include/craylibs"
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
backend="${1:?cpu or hip}"
export TMPDIR="$PWD/final-research/release-build-${backend}/tmp"
mkdir -p "$TMPDIR"
if [[ "$backend" == hip ]]; then
  cmake -S final-research/source-release -B final-research/release-build-hip \
    -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=ON -DKokkos_ENABLE_HIP=ON \
    -DKokkos_ARCH_ZEN3=ON -DKokkos_ARCH_VEGA90A=ON -DCMAKE_CXX_COMPILER=CC \
    -DCMAKE_CXX_FLAGS="-fno-cray -mno-daz-ftz" -DCMAKE_EXE_LINKER_FLAGS="-no-pie"
elif [[ "$backend" == cpu ]]; then
  module unload craype-accel-amd-gfx90a
  cmake -S final-research/source-release -B final-research/release-build-cpu \
    -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=ON -DKokkos_ENABLE_SERIAL=ON \
    -DCMAKE_CXX_COMPILER=CC -DCMAKE_EXE_LINKER_FLAGS="-ffp-model=precise -mno-daz-ftz"
else
  exit 2
fi
cmake --build "final-research/release-build-${backend}" --parallel 8
sha256sum "final-research/release-build-${backend}/src/athena"
