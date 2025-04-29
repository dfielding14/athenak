module load cmake

build="build"

#cmake \
#    -D CMAKE_POLICY_DEFAULT_CMP0074=NEW \
#    -D Athena_ENABLE_MPI=ON \
#    -D CMAKE_CXX_COMPILER=icpx \
#    -D Kokkos_ENABLE_SERIAL=ON \
#    -D Kokkos_ENABLE_SYCL=ON \
#    -D Kokkos_ARCH_INTEL_PVC=ON \
#    -D Kokkos_ENABLE_SYCL_RELOCATABLE_DEVICE_CODE=ON \
#    -D Kokkos_ENABLE_OPENMP=ON \
#    -D CMAKE_INSTALL_PREFIX=$build \
#    -D PROBLEM=cgm_cooling_flow \
#    -B $build

cmake \
    -D Athena_ENABLE_MPI=ON \
    -D CMAKE_CXX_COMPILER=icpx \
    -D Kokkos_ENABLE_SERIAL=ON \
    -D Kokkos_ENABLE_SYCL=ON \
    -D Kokkos_ARCH_INTEL_PVC=ON \
    -D Kokkos_ENABLE_SYCL_RELOCATABLE_DEVICE_CODE=ON \
    -D PROBLEM=turb \
    -B $build

cd $build
#make clean
make -j 104 all
