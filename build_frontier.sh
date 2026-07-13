#!/bin/bash

# Simple Frontier build template for AthenaK.
# Example: PROBLEM=turb BUILD_DIR=build_frontier_turb ./build_frontier.sh

set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
build_dir=${BUILD_DIR:-"${repo_root}/build_frontier"}
problem=${PROBLEM:-built_in_pgens}
build_jobs=${BUILD_JOBS:-16}

source /opt/cray/pe/lmod/lmod/init/bash
module reset
module load cpe/24.11 PrgEnv-amd/8.6.0 amd/6.2.4 rocm/6.2.4
module load cray-mpich/8.1.31 craype-accel-amd-gfx90a
if module is-loaded darshan-runtime; then
  module unload darshan-runtime
fi

printf 'Configuring AthenaK for Frontier\n'
printf '  problem: %s\n' "${problem}"
printf '  build:   %s\n' "${build_dir}"

/usr/bin/cmake -S "${repo_root}" -B "${build_dir}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DAthena_SINGLE_PRECISION=OFF \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_AMD_GFX90A=ON \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include -munsafe-fp-atomics" \
  -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64" \
  -DPROBLEM="${problem}"

/usr/bin/cmake --build "${build_dir}" --parallel "${build_jobs}"

printf 'Build complete: %s/src/athena\n' "${build_dir}"
