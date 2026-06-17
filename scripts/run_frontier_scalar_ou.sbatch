#!/bin/bash
# Submit with:
#   sbatch -A PROJECT_ID [-N NODES] scripts/run_frontier_scalar_ou.sbatch \
#     [block/parameter=value ...]

#SBATCH -J scalar-ou
#SBATCH -o %x-%j.out
#SBATCH -t 01:00:00
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-node=8

set -Eeuo pipefail

# Keep this tuple synchronized with the Frontier build. Override all four
# module variables together when using another OLCF-compatible toolchain.
cpe_module="${FRONTIER_CPE_MODULE:-cpe/25.09}"
mpich_module="${FRONTIER_MPICH_MODULE:-cray-mpich/9.0.1}"
rocm_module="${FRONTIER_ROCM_MODULE:-rocm/6.4.2}"
cce_module="${FRONTIER_CCE_MODULE:-cce/20.0.0}"

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load "${cpe_module}" "${mpich_module}" "${rocm_module}"
module load "${cce_module}"
module unload darshan-runtime

: "${CRAY_LD_LIBRARY_PATH:?Frontier modules did not define CRAY_LD_LIBRARY_PATH}"
: "${SLURM_JOB_ID:?This script must run inside a Slurm allocation}"
: "${SLURM_JOB_ACCOUNT:?Submit with sbatch -A PROJECT_ID}"
: "${SLURM_NNODES:?Slurm did not define SLURM_NNODES}"

export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_VERSION_DISPLAY="${MPICH_VERSION_DISPLAY:-1}"
export FI_MR_CACHE_MONITOR=kdreg2
export OMP_NUM_THREADS=1

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
athena_exe="${ATHENA_EXE:-${repo_dir}/build_scalar_mixing/src/athena}"
athena_input="${ATHENA_INPUT:-${repo_dir}/inputs/hydro/scalar_mixing_stochastic_forcing.athinput}"

if [[ -z "${RUN_DIR:-}" ]]; then
  : "${MEMBERWORK:?Set RUN_DIR or load the Frontier environment defining MEMBERWORK}"
  run_dir="${MEMBERWORK}/${SLURM_JOB_ACCOUNT}/scalar_ou/${SLURM_JOB_ID}"
else
  run_dir="${RUN_DIR}"
fi

if [[ ! -x "${athena_exe}" ]]; then
  echo "ERROR: AthenaK executable not found or not executable: ${athena_exe}" >&2
  exit 1
fi
if [[ ! -r "${athena_input}" ]]; then
  echo "ERROR: AthenaK input file not readable: ${athena_input}" >&2
  exit 1
fi

mkdir -p "${run_dir}"
cp "${athena_input}" "${run_dir}/submitted.athinput"

echo "=== Frontier scalar OU run ==="
echo "nodes:       ${SLURM_NNODES}"
echo "MPI ranks:   $((8 * SLURM_NNODES))"
echo "executable:  ${athena_exe}"
echo "input:       ${athena_input}"
echo "run dir:     ${run_dir}"
echo "modules:"
module -t list

srun -N "${SLURM_NNODES}" \
     -n "$((8 * SLURM_NNODES))" \
     --ntasks-per-node=8 \
     --cpus-per-task=7 \
     --threads-per-core=1 \
     --gpus-per-task=1 \
     --gpu-bind=closest \
     --unbuffered \
     "${athena_exe}" -d "${run_dir}" -i "${athena_input}" "$@"
