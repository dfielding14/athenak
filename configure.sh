module --force purge
module load modules/2.3-20240529 cuda/12.3.2 openmpi/cuda-4.0.7
export LD_PRELOAD=/mnt/sw/fi/cephtweaks/lib/libcephtweaks.so
export CEPHTWEAKS_LAZYIO=1

athenak="/mnt/home/btan1/Work/athenak"
build="build"

cmake \
  -D CMAKE_CXX_COMPILER=$athenak/kokkos/bin/nvcc_wrapper \
  -D Kokkos_ENABLE_CUDA=On \
  -D Kokkos_ARCH_AMPERE80=On \
  -D Athena_ENABLE_MPI=On \
  -D PROBLEM=cgm_cooling_flow\
  -B $build

cd $build
#make clean
make -j 16
