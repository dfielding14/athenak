#!/bin/bash
set -uo pipefail
cd /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2
module restore
module load PrgEnv-cray craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2 cce/20.0.0
module load cray-python/3.11.7
module unload darshan-runtime
unset CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH
export INCLUDE_PATH_X86_64="$(CC --print-resource-dir)/include:${CC_X86_64}/include/craylibs"
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export MPICH_GPU_SUPPORT_ENABLED=1 MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export MPICH_OFI_NIC_POLICY=GPU MPICH_GPU_IPC_CACHE_MAX_SIZE=1000 HSA_XNACK=1
export MPICH_MPIIO_HINTS='*:romio_cb_write=disable'
export MPICH_OFI_NUM_CQ_ENTRIES=131072 FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$PWD/baseline/tmp" MPLCONFIGDIR="$PWD/baseline/mpl-cache"
export XDG_CACHE_HOME="$PWD/baseline/cache"
export PATH="$PWD/scripts:$PATH"
export WO2_BACKEND=hip WO2_BINARY="$PWD/baseline/bin/athena-hip-nocray"
export WO2_TEST_BINARY="$WO2_BINARY" WO2_TEST_LABEL=hip-nocray
export CXX="$PWD/compiler-research/wrappers/wo2-cxx-noftz"
export WO2_PYTEST_K="_cpu and not test_cgl_lf_stage_i"
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$XDG_CACHE_HOME"


backend="${1:?cpu or hip}"
suite="${2:?cpu/gpu/mpi/mpi-gpu}"
label="${3:?unique label}"
export WO2_TEST_SOURCE="${WO2_TEST_SOURCE:-/autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2}"
export WO2_TEST_BINARY="$PWD/final-research/bin/athena-${backend}"
export WO2_TEST_LABEL="$label"
export WO2_PYTEST_K="${4:-_cpu and not test_cgl_lf_stage_i}"
export ATHENAK_RUN_MPI_GPU=1
export ATHENAK_MPI_GPU_LAUNCHER="srun --exact --threads-per-core=1 --cpu-bind=threads -c7 --gpus-per-task=1 --gpu-bind=closest"
if [[ "$backend" == cpu ]]; then export MPICH_GPU_SUPPORT_ENABLED=0; fi
python3 scripts/pytest_baseline.py "$backend" "$suite" > "final-research/${label}-${suite}.log" 2>&1
