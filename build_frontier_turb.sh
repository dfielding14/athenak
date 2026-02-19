#!/bin/bash

# Exit on any error
set -e

# Create log file with timestamp
log_file="build_frontier_turb_$(date +%Y%m%d_%H%M%S).log"
echo "Build log will be saved to: ${log_file}"

# Function to echo to both stdout and log file
log_echo() {
    echo "$@" | tee -a "${log_file}"
}

# Redirect all output to both terminal and log file
exec > >(tee -a "${log_file}")
exec 2>&1

# Frontier-optimized build script for AthenaK with AMR turbulence
log_echo "=== Building AthenaK on Frontier with AMR turbulence ==="
log_echo "=== Build started at $(date) ==="

# Load optimal modules for Frontier
module restore
module load cpe/24.07 PrgEnv-amd cray-mpich/8.1.30 craype-accel-amd-gfx90a amd/6.2.0 rocm/6.2.0
module load cmake cray-python emacs
module unload darshan-runtime


echo "=== Loaded modules ==="
module -t list

# Set environment variables
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

# Define paths
athenak_dir='/ccs/home/dfielding/athenak-df'
build_dir="${athenak_dir}/build_turb"

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
      -DCMAKE_CXX_FLAGS="-I${ROCM_PATH}/include -munsafe-fp-atomics" \
      -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/lib -lamdhip64" \
      -DPROBLEM=turb_timed_amr

# Build
echo "=== Building AthenaK ==="
cd "${build_dir}"
make -j16

echo "=== Build complete! ==="
echo "Executable location: ${build_dir}/src/athena"
echo "=== Build finished at $(date) ==="
echo ""
echo "Build log saved to: ${log_file}"
