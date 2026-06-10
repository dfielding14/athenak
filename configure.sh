#!/usr/bin/env bash

# Frontier HIP build.
# Optional overrides:
#   BUILD_DIR=build-frontier BUILD_JOBS=16 BUILD_TYPE=Release ./configure.sh
#   CONFIGURE_ONLY=1 ./configure.sh
#   CPE_VERSION=26.03 MPICH_VERSION=9.1.0 ROCM_VERSION=7.2.0 \
#     CCE_VERSION=21.0.0 BUILD_DIR=build-frontier-cpe2603-rocm720 ./configure.sh
# Additional arguments are forwarded to CMake after the defaults below.

set -euo pipefail

source_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${BUILD_DIR:-build-frontier}"
build_jobs="${BUILD_JOBS:-16}"
build_type="${BUILD_TYPE:-Release}"
configure_only="${CONFIGURE_ONLY:-0}"
cpe_version="${CPE_VERSION:-25.09}"
mpich_version="${MPICH_VERSION:-9.0.1}"
rocm_version="${ROCM_VERSION:-6.4.2}"
cce_version="${CCE_VERSION:-20.0.0}"
cmake_version="${CMAKE_VERSION:-3.30.5}"

if [[ "${build_dir}" != /* ]]; then
  build_dir="${source_dir}/${build_dir}"
fi

if ! [[ "${build_jobs}" =~ ^[1-9][0-9]*$ ]]; then
  echo "BUILD_JOBS must be a positive integer; got '${build_jobs}'." >&2
  exit 2
fi

if ! type module >/dev/null 2>&1; then
  echo "The environment-modules command is unavailable." >&2
  exit 1
fi

module reset
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load "cpe/${cpe_version}" "cray-mpich/${mpich_version}" "rocm/${rocm_version}"
module load "cce/${cce_version}"
module load "cmake/${cmake_version}"
if module is-loaded darshan-runtime; then
  module unload darshan-runtime
fi

if [[ ! -f "${source_dir}/kokkos/CMakeLists.txt" ]]; then
  echo "Kokkos submodule is not initialized. Run:" >&2
  echo "  git -C '${source_dir}' submodule update --init --recursive" >&2
  exit 1
fi

export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export MPICH_GPU_SUPPORT_ENABLED=1
export FI_MR_CACHE_MONITOR=kdreg2

cmake \
  -S "${source_dir}" \
  -B "${build_dir}" \
  -DCMAKE_BUILD_TYPE="${build_type}" \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_VEGA90A=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64" \
  -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include -munsafe-fp-atomics" \
  -DPROBLEM=uniform_hydro \
  "$@"

if [[ "${configure_only}" == "1" ]]; then
  echo "Configured AthenaK in ${build_dir}"
else
  cmake --build "${build_dir}" --parallel "${build_jobs}"
fi
