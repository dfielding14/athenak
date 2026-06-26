#!/bin/bash
# Compact nonlinear-window Bell engineering run.
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --job-name=pic-nonlinear-bell-engineering
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

repo_root="${PIC_SOURCE_ROOT:-/autofs/nccs-svm1_home2/dfielding/athenak-pic}"
executable="${PIC_EXECUTABLE:-/lustre/orion/ast207/proj-shared/dfielding/PIC/bin/a8eed16a6dd3/hip-mpi-release-paper-pic/athena}"
deck="${repo_root}/inputs/publication/q019_q023_carrier_nonlinear_bell_redesign_v1/q019-q023-carrier-s2-resolution-coarse-s0.athinput"
output_root="${PIC_OUTPUT_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/PIC/engineering_readiness/nonlinear-bell-${SLURM_JOB_ID}}"

export PIC_FRONTIER_PROFILE=frontier_minimum_supported
source "${repo_root}/tst/publication/frontier_control_plane/frontier_pic_environment.sh"
export OMP_NUM_THREADS=1
export MPICH_GPU_EAGER_REGISTER_HOST_MEM=0
export MPICH_GPU_NO_ASYNC_COPY=1
export MPICH_GPU_IPC_ENABLED=0

mkdir -p "${output_root}"
sha256sum "${executable}" "${deck}" >"${output_root}/bindings.sha256"
/usr/bin/time -p -o "${output_root}/run.time" \
  srun -N1 -n8 --ntasks-per-node=8 --cpus-per-task=1 \
    --gpus-per-task=1 --gpu-bind=closest \
    "${executable}" -i "${deck}" -d "${output_root}" \
    >"${output_root}/athena.stdout" \
    2>"${output_root}/athena.stderr"

cd "${repo_root}"
/opt/cray/pe/python/3.11.7/bin/python3 -B -m \
  tst.publication.bell_engineering_quicklook_v1 \
  nonlinear "${output_root}" --output "${output_root}/quicklook.json" --require-pass

printf '%s\n' 'execution_complete=true' 'physics_analysis_complete=true' \
  >"${output_root}/engineering_readiness.status"
