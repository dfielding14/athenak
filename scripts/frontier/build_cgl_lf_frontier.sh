#!/bin/bash
set -euo pipefail

# Produce an immutable baseline HIP/MPI executable for CGL-LF qualification.
# Run this script on Frontier after verifying the live module stack.

CGL_ROOT="${CGL_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/CGL}"
SRC_DIR="${SRC_DIR:-${CGL_ROOT}/repo/athenak-DF}"
BUILD_JOBS="${BUILD_JOBS:-16}"
STACK_TAG="${STACK_TAG:-cpe25.09-cce20-rocm6.4.2}"

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime

test -d "${SRC_DIR}/.git"
GIT_SHA="$(git -C "${SRC_DIR}" rev-parse HEAD)"
GIT_SHORT="$(git -C "${SRC_DIR}" rev-parse --short=12 HEAD)"
BUILD_DIR="${BUILD_DIR:-${CGL_ROOT}/build/frontier-hip-${GIT_SHORT}-${STACK_TAG}}"
BUILD_LOG_DIR="${CGL_ROOT}/logs/build"
MANIFEST_DIR="${CGL_ROOT}/runs/build-manifests/${GIT_SHORT}-${STACK_TAG}"
mkdir -p "${BUILD_DIR}" "${BUILD_LOG_DIR}" "${MANIFEST_DIR}"

if ! git -C "${SRC_DIR}" diff --quiet --ignore-submodules -- ||
   ! git -C "${SRC_DIR}" diff --cached --quiet --ignore-submodules --; then
  echo "Refusing to build Frontier evidence from tracked source changes." >&2
  exit 1
fi

{
  date -u +"created_utc=%Y-%m-%dT%H:%M:%SZ"
  hostname
  echo "source_dir=${SRC_DIR}"
  echo "git_revision=${GIT_SHA}"
  echo "build_dir=${BUILD_DIR}"
  module -t list 2>&1
  CC --version 2>&1 | head -n 4
  cmake --version | head -n 1
} | tee "${MANIFEST_DIR}/environment.txt"

cmake -S "${SRC_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DPROBLEM=built_in_pgens \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_VEGA90A=ON \
  -DCMAKE_CXX_COMPILER=CC \
  2>&1 | tee "${BUILD_LOG_DIR}/configure-${GIT_SHORT}-${STACK_TAG}.log"

cmake --build "${BUILD_DIR}" --parallel "${BUILD_JOBS}" \
  2>&1 | tee "${BUILD_LOG_DIR}/build-${GIT_SHORT}-${STACK_TAG}.log"

ATHENA="${BUILD_DIR}/src/athena"
test -x "${ATHENA}"
sha256sum "${ATHENA}" | tee "${MANIFEST_DIR}/athena.sha256"
"${ATHENA}" -c > "${MANIFEST_DIR}/athena-config.txt"
cp "${BUILD_DIR}/CMakeCache.txt" "${MANIFEST_DIR}/CMakeCache.txt"

printf 'Built executable: %s\nGit revision: %s\n' "${ATHENA}" "${GIT_SHA}"
