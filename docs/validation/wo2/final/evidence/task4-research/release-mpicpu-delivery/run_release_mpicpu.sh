#!/bin/bash
set -euo pipefail
cd /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/task4-research
: "${SLURM_JOB_ID:?Coordinated allocation required}"
label=${1:?Fresh evidence label}
module restore
module load PrgEnv-cray cpe/25.09 cray-mpich/9.0.1 cce/20.0.0 cray-python/3.11.7
module unload darshan-runtime rocm craype-accel-amd-gfx90a
unset CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH
export INCLUDE_PATH_X86_64="$(CC --print-resource-dir)/include:${CC_X86_64}/include/craylibs"
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export MPICH_GPU_SUPPORT_ENABLED=0 MPICH_MPIIO_HINTS='*:romio_cb_write=disable'
export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
final=/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research
export ATHENAK_CGL_PASSIVE_BINARY="$final/release-bin/athena-cpu"
export ATHENAK_CGL_PASSIVE_MPI_LAUNCHER='srun --exact -N1 -n4 --ntasks-per-node=4 --threads-per-core=1 --cpu-bind=cores -c1 --gpus-per-node=0 --overlap'
if [[ -n ${TASK4_NODE:-} ]]; then
  export ATHENAK_CGL_PASSIVE_MPI_LAUNCHER="$ATHENAK_CGL_PASSIVE_MPI_LAUNCHER -w $TASK4_NODE"
fi
export PYTHONPATH="$final/source-release/tst${PYTHONPATH:+:$PYTHONPATH}"
export TMPDIR="$PWD/release-mpicpu-runs/$label/tmp"
export MPLCONFIGDIR="$PWD/release-mpicpu-runs/$label/mpl-cache"
export XDG_CACHE_HOME="$PWD/release-mpicpu-runs/$label/cache"
[[ ! -e "release-mpicpu-runs/$label" ]] || { echo 'Fresh label required' >&2; exit 2; }
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
module -t list > "release-mpicpu-runs/$label/modules.txt" 2>&1
python3 run_release_mpicpu.py "$label"
