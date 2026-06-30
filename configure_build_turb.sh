#!/usr/bin/env bash
# Configure AthenaK for the turb problem with the Frontier toolchain and
# CMake options used by the 2026-06-26 CGL weak-scaling campaign.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${BUILD_DIR:-${script_dir}/build-frontier-hip-cpe25.09-cce20-rocm6.4.2-turb}"
build_jobs="${BUILD_JOBS:-16}"

if git -C "$script_dir" submodule status --recursive | grep -Eq '^[+-U]'; then
  echo "The recorded Kokkos submodule is not initialized at the committed revision." >&2
  echo "Run: git -C $script_dir submodule update --init --recursive" >&2
  exit 2
fi

module reset
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
if module is-loaded darshan-runtime; then
  module unload darshan-runtime
fi

# Prevent inherited interactive CCE 18 include paths from contaminating the
# selected CCE 20 wrapper environment. These are the exact environment fixes
# used for the scaling executable.
export INCLUDE_PATH_X86_64=/opt/cray/pe/cce/20.0.0/cce-clang/x86_64/lib/clang/20/include:/opt/cray/pe/cce/20.0.0/cce/x86_64/include/craylibs
unset CPLUS_INCLUDE_PATH C_INCLUDE_PATH
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"

printf 'Source directory: %s\nBuild directory:  %s\nLoaded modules:\n' \
  "$script_dir" "$build_dir"
module -t list 2>&1

cmake --fresh -S "$script_dir" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DPROBLEM=turb \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_VEGA90A=ON \
  -DCMAKE_CXX_COMPILER=CC

printf '\nConfigured successfully. Building with %s parallel jobs...\n' "$build_jobs"
cmake --build "$build_dir" --parallel "$build_jobs"

printf '\nBuild complete. Executable: %s\n' "$build_dir/src/athena"
