#!/usr/bin/env bash

# Configure and build the full GOTHAM AthenaK problem on Frontier.  The script
# deliberately re-executes in a clean environment: after Frontier's July 2026
# module transition, `module restore` can retain CCE 18 include paths and mix
# them with CCE 20.

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  printf 'configure_and_build_gotham.sh: execute this script; do not source it\n' >&2
  return 1
fi

set -euo pipefail

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}")
REPO_ROOT=$(dirname "$SCRIPT_PATH")

die() {
  printf 'configure_and_build_gotham.sh: %s\n' "$*" >&2
  exit 1
}

usage() {
  printf '%s\n' \
    'Usage: ./configure_and_build_gotham.sh [OPTIONS]' \
    '' \
    'Configure and build GOTHAM with Frontier HIP and MPI.' \
    '' \
    'Options:' \
    '  --problem NAME      AthenaK problem generator (default: gotham)' \
    '  --build-dir DIR     Build directory (default: build)' \
    '  --build-type TYPE   Release, Debug, RelWithDebInfo, or MinSizeRel' \
    '                      (default: Release)' \
    '  --jobs N            Parallel build jobs (default: 16)' \
    '  --reuse-cache       Reuse the existing CMake cache instead of --fresh' \
    '  --configure-only    Configure without compiling' \
    '  -h, --help          Show this help' \
    '' \
    'Equivalent defaults may be set with GOTHAM_BUILD_DIR,' \
    'GOTHAM_BUILD_TYPE, GOTHAM_BUILD_JOBS, and GOTHAM_PROBLEM.'
}

# Keep only the small set of caller values needed by the build.  This prevents
# stale compiler variables from surviving a CPE/module-family transition.
if [[ "${GOTHAM_FRONTIER_CLEAN_ENV:-0}" != "1" ]]; then
  caller_user=${USER:-$(id -un)}
  caller_home=${HOME:-}
  [[ -n "$caller_home" ]] || die 'HOME is not set'
  exec /usr/bin/env -i \
    HOME="$caller_home" \
    USER="$caller_user" \
    LOGNAME="${LOGNAME:-$caller_user}" \
    SHELL=/bin/bash \
    PATH=/usr/local/bin:/usr/bin:/bin:/opt/cray/pe/lmod/lmod/libexec \
    GOTHAM_FRONTIER_CLEAN_ENV=1 \
    GOTHAM_BUILD_DIR="${GOTHAM_BUILD_DIR:-}" \
    GOTHAM_BUILD_TYPE="${GOTHAM_BUILD_TYPE:-}" \
    GOTHAM_BUILD_JOBS="${GOTHAM_BUILD_JOBS:-}" \
    GOTHAM_PROBLEM="${GOTHAM_PROBLEM:-}" \
    /bin/bash "$SCRIPT_PATH" "$@"
fi
[[ -z "${LOADEDMODULES:-}" && -z "${INCLUDE_PATH_X86_64:-}" ]] || {
  die 'the internal clean-environment marker was set in a contaminated shell'
}

build_dir=${GOTHAM_BUILD_DIR:-build}
build_type=${GOTHAM_BUILD_TYPE:-Release}
build_jobs=${GOTHAM_BUILD_JOBS:-16}
problem=${GOTHAM_PROBLEM:-gotham}
fresh_configure=1
configure_only=0

while (( $# > 0 )); do
  case "$1" in
    --problem)
      (( $# >= 2 )) || die '--problem requires a value'
      problem=$2
      shift 2
      ;;
    --build-dir)
      (( $# >= 2 )) || die '--build-dir requires a value'
      build_dir=$2
      shift 2
      ;;
    --build-type)
      (( $# >= 2 )) || die '--build-type requires a value'
      build_type=$2
      shift 2
      ;;
    --jobs)
      (( $# >= 2 )) || die '--jobs requires a value'
      build_jobs=$2
      shift 2
      ;;
    --reuse-cache)
      fresh_configure=0
      shift
      ;;
    --configure-only)
      configure_only=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown option: $1"
      ;;
  esac
done

[[ "$build_jobs" =~ ^[1-9][0-9]*$ ]] || die '--jobs must be a positive integer'
[[ "$problem" =~ ^[A-Za-z0-9_]+$ ]] || die '--problem contains invalid characters'
case "$build_type" in
  Release|Debug|RelWithDebInfo|MinSizeRel) ;;
  *) die "unsupported build type: $build_type" ;;
esac

if [[ "$build_dir" != /* ]]; then
  build_dir="${REPO_ROOT}/${build_dir}"
fi
build_dir=$(realpath -m "$build_dir")
[[ "$build_dir" != "$REPO_ROOT" ]] || die 'the build directory cannot be the source root'
if [[ "$problem" != built_in_pgens ]]; then
  [[ -f "${REPO_ROOT}/src/pgen/${problem}.cpp" ]] || {
    die "problem generator not found: src/pgen/${problem}.cpp"
  }
fi

lmod_init=/opt/cray/pe/lmod/lmod/init/bash
[[ -r "$lmod_init" ]] || die "Lmod initializer not found: $lmod_init"
# shellcheck source=/dev/null
source "$lmod_init"

module_roots=(
  /opt/cray/pe/lmod/modulefiles/core
  /opt/cray/pe/lmod/lmod/modulefiles/Core
  /opt/cray/pe/lmod/modulefiles/craype-targets/default
  /sw/frontier/spack-envs/modules/Core/25.03
)
for module_root in "${module_roots[@]}"; do
  [[ -d "$module_root" ]] || die "required module path not found: $module_root"
done
module use "${module_roots[@]}"
module load cpe/25.09
module load craype-x86-trento craype-network-ofi craype-accel-amd-gfx90a
module load PrgEnv-cray cray-mpich/9.0.1 rocm/6.4.2 cce/20.0.0 cmake/3.30.5

required_modules=(
  cpe/25.09
  PrgEnv-cray/8.6.0
  craype-x86-trento
  craype-network-ofi
  craype-accel-amd-gfx90a
  cray-mpich/9.0.1
  rocm/6.4.2
  cce/20.0.0
  cmake/3.30.5
)
loaded_modules=$(module -t list 2>&1)
for required_module in "${required_modules[@]}"; do
  grep -Fxq "$required_module" <<<"$loaded_modules" || {
    die "required module was not loaded: $required_module"
  }
done

command -v CC >/dev/null 2>&1 || die 'Cray C++ wrapper CC is unavailable'
command -v cc >/dev/null 2>&1 || die 'Cray C wrapper cc is unavailable'
command -v cmake >/dev/null 2>&1 || die 'cmake is unavailable'
: "${ROCM_PATH:?ROCM_PATH was not set by the ROCm module}"
[[ -d "${ROCM_PATH}/include" && -d "${ROCM_PATH}/lib" ]] || {
  die "incomplete ROCm installation: ${ROCM_PATH}"
}

export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH:-}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

cmake_args=(
  -S "$REPO_ROOT"
  -B "$build_dir"
  -DCMAKE_BUILD_TYPE="$build_type"
  -DCMAKE_C_COMPILER=cc
  -DCMAKE_CXX_COMPILER=CC
  -DMPI_CXX_COMPILER=CC
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
  -DAthena_ENABLE_MPI=ON
  -DKokkos_ENABLE_HIP=ON
  -DKokkos_ARCH_ZEN3=ON
  -DKokkos_ARCH_VEGA90A=ON
  -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64"
  -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include -munsafe-fp-atomics"
  -DPROBLEM="$problem"
)

printf 'Configuring GOTHAM AthenaK\n'
printf '  source:     %s\n' "$REPO_ROOT"
printf '  build:      %s\n' "$build_dir"
printf '  build type: %s\n' "$build_type"
printf '  problem:    %s\n' "$problem"
printf '  jobs:       %s\n' "$build_jobs"

if (( fresh_configure )); then
  cmake --fresh "${cmake_args[@]}"
else
  cmake "${cmake_args[@]}"
fi

if (( configure_only )); then
  printf 'Configuration complete: %s\n' "$build_dir"
  exit 0
fi

cmake --build "$build_dir" --parallel "$build_jobs"

athena_executable="${build_dir}/src/athena"
[[ -x "$athena_executable" ]] || die "build did not produce ${athena_executable}"
printf 'Build complete: %s\n' "$athena_executable"
