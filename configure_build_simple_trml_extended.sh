#!/usr/bin/env bash
# Configure and build AthenaK's simple_TRML_extended problem on Frontier.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${BUILD_DIR:-${script_dir}/build-frontier-simple-trml-extended}"

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

export INCLUDE_PATH_X86_64=/opt/cray/pe/cce/20.0.0/cce-clang/x86_64/lib/clang/20/include:/opt/cray/pe/cce/20.0.0/cce/x86_64/include/craylibs
unset CPLUS_INCLUDE_PATH C_INCLUDE_PATH
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"

printf 'Source directory: %s\nBuild directory:  %s\nLoaded modules:\n' \
  "$script_dir" "$build_dir"
module -t list 2>&1

cmake --fresh -S "$script_dir" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DPROBLEM=simple_TRML_extended \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_VEGA90A=ON \
  -DCMAKE_CXX_FLAGS=-munsafe-fp-atomics \
  -DCMAKE_CXX_COMPILER=CC

printf '\nConfigured successfully. Building with 16 parallel jobs...\n'
cmake --build "$build_dir" --parallel 16
