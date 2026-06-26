#!/bin/bash
# Compact shock engineering gate. This is not publication qualification.
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --job-name=pic-shock-engineering-readiness
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

repo_root="${PIC_SOURCE_ROOT:-/autofs/nccs-svm1_home2/dfielding/athenak-pic}"
executable="${PIC_EXECUTABLE:-/lustre/orion/ast207/proj-shared/dfielding/PIC/bin/a8eed16a6dd3/hip-mpi-release-paper-pic/athena}"
run_id="shock-a8eed16a-${SLURM_JOB_ID}"
output_root="${PIC_OUTPUT_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/PIC/engineering_readiness/${run_id}}"

export PIC_FRONTIER_PROFILE=frontier_minimum_supported
source "${repo_root}/tst/publication/frontier_control_plane/frontier_pic_environment.sh"
export OMP_NUM_THREADS=1
export MPICH_GPU_EAGER_REGISTER_HOST_MEM=0
export MPICH_GPU_NO_ASYNC_COPY=1
export MPICH_GPU_IPC_ENABLED=0

mkdir -p "${output_root}"
sha256sum "${executable}" >"${output_root}/executable.sha256"

run_case() {
  local label="$1"
  local deck="$2"
  local basename="$3"
  shift 3
  /usr/bin/time -p -o "${output_root}/${label}.time" \
    srun -N1 -n1 --ntasks-per-node=1 --cpus-per-task=1 \
      --gpus-per-task=1 --gpu-bind=closest \
      "${executable}" -i "${repo_root}/${deck}" \
      -d "${output_root}" "job/basename=${basename}" "$@" \
      >"${output_root}/${label}.stdout" \
      2>"${output_root}/${label}.stderr"
}

run_case coarse \
  inputs/tests/pic_parallel_shock_coarse_uniform.athinput \
  "${run_id}_coarse" time/nlim=12 time/tlim=0.02
run_case amr \
  inputs/tests/pic_parallel_shock_amr_fiducial.athinput \
  "${run_id}_amr" time/nlim=12 time/tlim=0.02
run_case fine \
  inputs/tests/pic_parallel_shock_fine_uniform.athinput \
  "${run_id}_fine" time/nlim=12 time/tlim=0.02

# This case exercises live CR deposition, injection, momentum/energy feedback,
# and gas subtraction. It is deliberately one cycle and makes no acceleration claim.
run_case coupled_feedback \
  inputs/tests/pic_parallel_shock_section54_stage_timing_acceptance_vl2_tsc.athinput \
  "${run_id}_coupled_feedback"

cd "${repo_root}"
/opt/cray/pe/python/3.11.7/bin/python3 -B \
  tst/publication/run_pic_parallel_shock_benchmark.py \
  --athena-cwd "${output_root}" \
  --output-root "${output_root}/analysis" \
  --run-id "${run_id}"

printf '%s\n' \
  'engineering_ready=true' \
  'transport_triad=complete' \
  'coupled_feedback_smoke=complete' \
  'nonrelativistic_acceleration_claim=false' \
  >"${output_root}/engineering_readiness.status"
