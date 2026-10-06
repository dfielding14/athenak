#!/bin/bash
set -euo pipefail
cd /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2
backend=${1:?cpu or hip}
shift
module restore
module load PrgEnv-cray cpe/25.09 cray-mpich/9.0.1 cce/20.0.0 cray-python/3.11.7
module unload darshan-runtime
case "$backend" in
  hip)
    module load craype-accel-amd-gfx90a rocm/6.4.2
    export MPICH_GPU_SUPPORT_ENABLED=1 MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
    export MPICH_OFI_NIC_POLICY=GPU MPICH_GPU_IPC_CACHE_MAX_SIZE=1000 HSA_XNACK=1
    export MPICH_OFI_NUM_CQ_ENTRIES=131072 FI_MR_CACHE_MONITOR=kdreg2
    export FI_CXI_RX_MATCH_MODE=software
    ;;
  cpu)
    module unload craype-accel-amd-gfx90a
    export MPICH_GPU_SUPPORT_ENABLED=0
    ;;
  *) echo 'backend must be cpu or hip' >&2; exit 2 ;;
esac
unset CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH
export INCLUDE_PATH_X86_64="$(CC --print-resource-dir)/include:${CC_X86_64}/include/craylibs"
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export MPICH_MPIIO_HINTS='*:romio_cb_write=disable' OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$PWD/final-research/release-equivalence-runtime/tmp"
export MPLCONFIGDIR="$PWD/final-research/release-equivalence-runtime/mpl-cache"
export XDG_CACHE_HOME="$PWD/final-research/release-equivalence-runtime/cache"
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
python3 final-research/run_release_equivalence.py --backend "$backend" "$@"
