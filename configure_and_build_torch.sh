#!/usr/bin/env bash
# Configure and build AthenaK for the Tiegan SGS turbulence runs on Torch H200s.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${BUILD_DIR:-${script_dir}/build-torch-sgs-cuda12.9-h200}"
build_jobs="${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-8}}"

kokkos_status="$(git -C "$script_dir" submodule status --recursive)"
if grep -Eq '^[+-U]' <<< "$kokkos_status"; then
  echo 'Initialize Kokkos at the revision recorded by AthenaK.' >&2
  printf 'Run git -C %q submodule update --init --recursive\n' "$script_dir" >&2
  exit 2
fi

# Torch's container wrapper loads nvhpc-hpcx-cuda12-openmpi5/26.5.
# Build on a Torch compute node and use the same wrapper to run AthenaK on H200s.
toolchain=/share/apps/images/run-nvhpc-26.5-cuda-12.9.bash

printf 'Source directory %s\nBuild directory  %s\n' "$script_dir" "$build_dir"

"$toolchain" cmake --fresh -S "$script_dir" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DPROBLEM=turb \
  -DAthena_ENABLE_MPI=ON \
  -DAthena_SINGLE_PRECISION=OFF \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ARCH_HOPPER90=ON \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_CXX_COMPILER="$script_dir/kokkos/bin/nvcc_wrapper"

printf '\nBuilding with %s parallel jobs\n' "$build_jobs"
"$toolchain" cmake --build "$build_dir" --parallel "$build_jobs"

printf '\nBuilt executable at %s\n' "$build_dir/src/athena"
