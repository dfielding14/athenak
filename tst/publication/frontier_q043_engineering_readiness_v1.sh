#!/bin/bash
# One-allocation engineering gate. This is not registered publication evidence.
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=pic-q043-engineering-readiness
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

repo_root="${PIC_SOURCE_ROOT:-/autofs/nccs-svm1_home2/dfielding/athenak-pic}"
executable="${PIC_EXECUTABLE:-/lustre/orion/ast207/proj-shared/dfielding/PIC/bin/a8eed16a6dd3/hip-mpi-release-paper-pic/athena}"
output_root="${PIC_OUTPUT_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/PIC/engineering_readiness/q043-${SLURM_JOB_ID}}"

export PIC_FRONTIER_PROFILE=frontier_minimum_supported
source "${repo_root}/tst/publication/frontier_control_plane/frontier_pic_environment.sh"

# These settings eliminate the intermittent post-run HSA/GTL teardown fault
# observed with CPE 24.11 while leaving the numerical kernels unchanged.
export MPICH_GPU_EAGER_REGISTER_HOST_MEM=0
export MPICH_GPU_NO_ASYNC_COPY=1
export MPICH_GPU_IPC_ENABLED=0

cd "${repo_root}"
/opt/cray/pe/python/3.11.7/bin/python3 -B -m \
  tst.publication.q043_engineering_readiness_v1 \
  --run \
  --executable "${executable}" \
  --output-root "${output_root}"
