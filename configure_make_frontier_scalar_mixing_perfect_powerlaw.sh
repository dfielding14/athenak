#!/bin/bash

# Exit on any error
set -e

# Create log file with timestamp
log_file="build_frontier_scalar_mixing_perfect_powerlaw_$(date +%Y%m%d_%H%M%S).log"
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
athenak_dir='/ccs/home/dfielding/athenak-df-fractal'
build_dir="${athenak_dir}/build_scalar_mixing_perfect_powerlaw"

# Clean and create build directory
echo "=== Setting up build directory ==="

if [ -d "${build_dir}" ]; then
    echo "Removing existing build directory..."
    rm -rf "${build_dir}"
fi
mkdir -p "${build_dir}"

# Configure with CMake
echo "=== Configuring with CMake ==="
cd "${athenak_dir}"

cmake -B"${build_dir}" \
      -DAthena_ENABLE_MPI=ON \
      -DKokkos_ARCH_ZEN3=ON \
      -DKokkos_ARCH_VEGA90A=ON \
      -DKokkos_ENABLE_HIP=ON \
      -DCMAKE_CXX_COMPILER=CC \
      -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include" \
      -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64" \
      -DPROBLEM=scalar_mixing_perfect_powerlaw

# Build
echo "=== Building AthenaK ==="
cd "${build_dir}"
make -j16

echo "=== Build complete! ==="
echo "Executable location: ${build_dir}/src/athena"
echo "=== Build finished at $(date) ==="
echo "Build log saved to: ${log_file}"
