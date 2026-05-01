#!/bin/bash

# Exit on any error
set -Eeuo pipefail

# Create log file with timestamp
timestamp="$(date +%Y%m%d_%H%M%S)"
log_file="build_frontier_scalar_mixing_perfect_powerlaw_${timestamp}.log"
echo "Build log will be saved to: ${log_file}"

# Function to echo to both stdout and log file
log_echo() {
    echo "$@" | tee -a "${log_file}"
}

# Redirect all output to both terminal and log file
exec > >(tee -a "${log_file}")
exec 2>&1

# Frontier-optimized build script for AthenaK with scalar_mixing_perfect_powerlaw problem
log_echo "=== Building AthenaK on Frontier with scalar_mixing_perfect_powerlaw ==="
log_echo "=== Build started at $(date) ==="

# Load the same Frontier module stack used by configure_divb_testing.sh.
module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module load cmake
module unload darshan-runtime


echo "=== Loaded modules ==="
module -t list

# Set environment variables
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export FI_MR_CACHE_MONITOR=kdreg2

# Define paths
athenak_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
build_dir="${athenak_dir}/build_scalar_mixing_perfect_powerlaw"
active_build_dir="${athenak_dir}/.build_scalar_mixing_perfect_powerlaw_${timestamp}_$$"
lock_dir="${build_dir}.lock"
jobs="${BUILD_JOBS:-16}"

release_lock() {
    if [ -r "${lock_dir}/pid" ] && [ "$(cat "${lock_dir}/pid" 2>/dev/null)" = "$$" ]; then
        rm -rf "${lock_dir}"
    fi
}

acquire_lock() {
    if mkdir "${lock_dir}" 2>/dev/null; then
        printf '%s\n' "$$" > "${lock_dir}/pid"
        trap release_lock EXIT
        return
    fi

    if [ -r "${lock_dir}/pid" ]; then
        lock_pid="$(cat "${lock_dir}/pid" 2>/dev/null || true)"
        if [ -n "${lock_pid}" ] && kill -0 "${lock_pid}" 2>/dev/null; then
            echo "ERROR: another scalar_mixing_perfect_powerlaw build is already running."
            echo "       build directory: ${build_dir}"
            echo "       lock held by PID: ${lock_pid}"
            echo "       Wait for it to finish, or remove ${lock_dir} if that PID is stale."
            exit 1
        fi
    fi

    echo "Removing stale build lock: ${lock_dir}"
    rm -rf "${lock_dir}"
    if ! mkdir "${lock_dir}" 2>/dev/null; then
        echo "ERROR: failed to acquire build lock after removing stale lock."
        echo "       another build may have started concurrently: ${lock_dir}"
        exit 1
    fi
    printf '%s\n' "$$" > "${lock_dir}/pid"
    trap release_lock EXIT
}

# Clean and create build directory
echo "=== Setting up build directory ==="
acquire_lock

echo "Building in private directory: ${active_build_dir}"
rm -rf "${active_build_dir}"
mkdir -p "${active_build_dir}"

# Configure with CMake
echo "=== Configuring with CMake ==="
cd "${athenak_dir}"

cmake -B"${active_build_dir}" \
      -DAthena_ENABLE_MPI=ON \
      -DKokkos_ARCH_ZEN3=ON \
      -DKokkos_ARCH_VEGA90A=ON \
      -DKokkos_ENABLE_HIP=ON \
      -DCMAKE_CXX_COMPILER=CC \
      -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include -munsafe-fp-atomics" \
      -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64" \
      -DPROBLEM=scalar_mixing_perfect_powerlaw

# Build
echo "=== Building AthenaK ==="
cmake --build "${active_build_dir}" --parallel "${jobs}"

echo "=== Publishing successful build ==="
previous_build_dir="${build_dir}.previous"
rm -rf "${previous_build_dir}"
if [ -e "${build_dir}" ]; then
    mv "${build_dir}" "${previous_build_dir}"
fi
mv "${active_build_dir}" "${build_dir}"
rm -rf "${previous_build_dir}"

echo "=== Build complete! ==="
echo "Executable location: ${build_dir}/src/athena"
echo "=== Build finished at $(date) ==="
echo "Build log saved to: ${log_file}"
