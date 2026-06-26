#!/bin/bash
# Configure, build, and install an exact-commit Q011 Frontier executable.
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=develop
#SBATCH --nodes=1
#SBATCH --time=00:20:00
#SBATCH --job-name=pic-s54-build
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

repo_root="${PIC_SOURCE_ROOT:-/autofs/nccs-svm1_home2/dfielding/athenak-pic}"
shared_root="${PIC_SHARED_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/PIC}"
source_commit="$(git -C "${repo_root}" rev-parse HEAD)"
requested_commit="${PIC_BUILD_COMMIT:-${source_commit}}"

if [[ "${requested_commit}" != "${source_commit}" ]]; then
  printf 'PIC_BUILD_COMMIT=%s does not match source HEAD=%s\n' \
    "${requested_commit}" "${source_commit}" >&2
  exit 2
fi
if [[ -n "$(git -C "${repo_root}" status --porcelain)" ]]; then
  printf 'Refusing to build a dirty or untracked source tree:\n' >&2
  git -C "${repo_root}" status --short >&2
  exit 2
fi

commit_key="${source_commit:0:12}"
profile="hip-mpi-release-paper-pic"
build_root="${shared_root}/build/${commit_key}/${profile}/cmake"
install_root="${shared_root}/bin/${commit_key}/${profile}"
receipt_root="${shared_root}/build/${commit_key}/${profile}/receipts/${SLURM_JOB_ID}"

export PIC_FRONTIER_PROFILE=frontier_minimum_supported
source "${repo_root}/tst/publication/frontier_control_plane/frontier_pic_environment.sh"

mkdir -p "${build_root}" "${install_root}" "${receipt_root}"
cmake -S "${repo_root}" -B "${build_root}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ARCH_AMD_GFX90A=ON \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64" \
  -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include" \
  -DPROBLEM=built_in_pgens
cmake --build "${build_root}" --parallel 16

install -m 0755 "${build_root}/src/athena" "${install_root}/athena"
printf '%s\n' "${source_commit}" >"${receipt_root}/source_commit.txt"
git -C "${repo_root}" status --porcelain=v1 >"${receipt_root}/source_status.txt"
sha256sum \
  "${install_root}/athena" \
  "${repo_root}/tst/publication/frontier_build_q011_stencil_v1.sh" \
  >"${receipt_root}/bindings.sha256"
cp "${build_root}/CMakeCache.txt" "${receipt_root}/CMakeCache.txt"
printf 'build_complete=true\nexecutable=%s\n' "${install_root}/athena" \
  >"${receipt_root}/build.status"

printf 'Q011 executable: %s\n' "${install_root}/athena"
sha256sum "${install_root}/athena"
