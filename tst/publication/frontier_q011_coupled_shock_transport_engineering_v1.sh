#!/bin/bash
# Compact coupled non-relativistic shock transport engineering run.
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --job-name=pic-shock-transport-engineering
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

repo_root="${PIC_SOURCE_ROOT:-/autofs/nccs-svm1_home2/dfielding/athenak-pic}"
executable="${PIC_EXECUTABLE:-/lustre/orion/ast207/proj-shared/dfielding/PIC/bin/a8eed16a6dd3/hip-mpi-release-paper-pic/athena}"
deck="${repo_root}/inputs/publication/pic_parallel_shock_section54_production_science_successor_v1_vl2_tsc.athinput"
basename="q011_shock_transport_engineering_${SLURM_JOB_ID}"
output_root="${PIC_OUTPUT_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/PIC/engineering_readiness/shock-transport-${SLURM_JOB_ID}}"

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
    "job/basename=${basename}" \
    mesh/nx1=256 mesh/x1max=768 mesh/nx2=64 mesh/x2max=192 \
    meshblock/nx1=32 meshblock/nx2=16 \
    mesh_refinement/refinement=none mesh_refinement/num_levels=1 \
    time/tlim=20 time/nlim=5000 time/ndiag=50 \
    particles/deposit_qscale=9.0e-4 \
    problem/ps_remove_birth_time_before=-1 \
    problem/ps_enable_curvature_amr=false \
    problem/ps_feedback_diag_dcycle=100 \
    output1/dt=1 output2/dt=1 output3/dt=1 output4/dt=1 \
    output5/dt=1 output6/dt=1 output7/dt=1 output8/dt=0.1 output9/dt=10 \
    >"${output_root}/athena.stdout" \
    2>"${output_root}/athena.stderr"

cd "${repo_root}"
/opt/cray/pe/python/3.11.7/bin/python3 -B -m \
  tst.publication.shock_transport_engineering_quicklook_v1 \
  "${output_root}" "${basename}" --output "${output_root}/quicklook.json" \
  --require-transport

printf '%s\n' \
  'execution_complete=true' \
  'coupled_shock_transport_analysis_complete=true' \
  >"${output_root}/engineering_readiness.status"
