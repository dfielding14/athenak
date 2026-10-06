#!/bin/bash
set -euo pipefail
cd /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/task2-research
module restore
module load PrgEnv-cray cpe/25.09 cray-mpich/9.0.1 cce/20.0.0 cray-python/3.11.7
module load craype-accel-amd-gfx90a rocm/6.4.2
module unload darshan-runtime
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export MPICH_GPU_SUPPORT_ENABLED=1 MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export MPICH_OFI_NIC_POLICY=GPU MPICH_GPU_IPC_CACHE_MAX_SIZE=1000 HSA_XNACK=1
export MPICH_MPIIO_HINTS='*:romio_cb_write=disable'
export MPICH_OFI_NUM_CQ_ENTRIES=131072 FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software OMP_NUM_THREADS=1
export MERGE_BACKEND=hip
export TMPDIR=$PWD/tmp PYTHONDONTWRITEBYTECODE=1
python run_merge_timing.py --label task1-refill-hip --binary integrated-task1-refill/build/athena-hip-merge --manifest integrated-task1-refill/build/manifest-hip.json
