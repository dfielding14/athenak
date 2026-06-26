#!/bin/bash
# Restartable Section 5.4 shock/DSA production run for Frontier.
#
# Fresh launch:
#   sbatch --export=ALL,Q011_FRESH=1,Q011_TLIM=1200,Q011_GRID_MODE=static_dx3 <this-file>
# Continuation:
#   sbatch --export=ALL,Q011_RUN_ROOT=<root>,Q011_TLIM=1200,Q011_GRID_MODE=static_dx3 <this-file>
#
# A continuation selects only restart files carrying AthenaK's complete marker
# and removes history rows newer than that checkpoint before execution.
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=normal
#SBATCH --nodes=92
#SBATCH --time=06:00:00
#SBATCH --job-name=pic-s54-dsa
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

repo_root="${PIC_SOURCE_ROOT:-/autofs/nccs-svm1_home2/dfielding/athenak-pic}"
shared_root="${PIC_SHARED_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/PIC}"
source_commit="$(git -C "${repo_root}" rev-parse HEAD)"
commit_key="${source_commit:0:12}"
executable="${PIC_EXECUTABLE:-${shared_root}/bin/${commit_key}/hip-mpi-release-paper-pic/athena}"
deck="${Q011_DECK:-${repo_root}/inputs/q011_section54_static_dx3_final_v1_vl2_tsc.athinput}"
target_time="${Q011_TLIM:-1200}"
grid_mode="${Q011_GRID_MODE:-static_dx3}"
basename="${Q011_BASENAME:-q011_section54_dsa}"

if [[ -n "$(git -C "${repo_root}" status --porcelain)" ]]; then
  printf 'Refusing to run from a dirty or untracked source tree:\n' >&2
  git -C "${repo_root}" status --short >&2
  exit 2
fi
if [[ ! -x "${executable}" ]]; then
  printf 'Q011 executable is absent or not executable: %s\n' "${executable}" >&2
  exit 2
fi
if [[ ! -f "${deck}" ]]; then
  printf 'Q011 deck is absent: %s\n' "${deck}" >&2
  exit 2
fi

if [[ "${Q011_FRESH:-0}" == "1" ]]; then
  run_root="${Q011_RUN_ROOT:-${shared_root}/production_science/shock-section54-${SLURM_JOB_ID}}"
  restart_file=""
else
  : "${Q011_RUN_ROOT:?Q011_RUN_ROOT is required for a continuation job}"
  run_root="${Q011_RUN_ROOT}"
  restart_file="$({
    find "${run_root}/rst" -maxdepth 1 -type f -name '*.rst.complete' -print0
  } | sort -zV | tail -z -n 1 | tr -d '\0')"
  if [[ -z "${restart_file}" ]]; then
    printf 'No complete restart is available below %s/rst\n' "${run_root}" >&2
    exit 1
  fi
  restart_file="${restart_file%.complete}"
fi

export PIC_FRONTIER_PROFILE=frontier_minimum_supported
source "${repo_root}/tst/publication/frontier_control_plane/frontier_pic_environment.sh"
export OMP_NUM_THREADS=1
export MPICH_GPU_EAGER_REGISTER_HOST_MEM=0
export MPICH_GPU_NO_ASYNC_COPY=1
export MPICH_GPU_IPC_ENABLED=0

mkdir -p "${run_root}" "${run_root}/segments"
segment_root="${run_root}/segments/${SLURM_JOB_ID}"
mkdir -p "${segment_root}"

segment_state=false
finish_segment() {
  rc=$?
  trap - EXIT
  printf 'execution_complete=%s\nexit_code=%s\ntarget_time=%s\n' \
    "${segment_state}" "${rc}" "${target_time}" >"${segment_root}/segment.status"
  exit "${rc}"
}
trap finish_segment EXIT

sha256sum \
  "${executable}" \
  "${deck}" \
  "${repo_root}/src/pgen/tests/pic_parallel_shock.cpp" \
  "${repo_root}/tst/publication/frontier_q011_section54_long_dsa_v1.sh" \
  "${repo_root}/tst/publication/prepare_q011_restart_history_v1.py" \
  "${repo_root}/tst/publication/make_q011_section54_figure8_v1.py" \
  "${repo_root}/tst/publication/make_q011_section54_evolution_v1.py" \
  "${repo_root}/tst/publication/make_q011_section54_dsa_spectrum_v1.py" \
  >"${segment_root}/bindings.sha256"
printf '%s\n' "${source_commit}" >"${segment_root}/source_commit.txt"
git -C "${repo_root}" status --porcelain=v1 >"${segment_root}/source_status.txt"
printf '%s\n' \
  "job_id=${SLURM_JOB_ID}" \
  "run_root=${run_root}" \
  "target_time=${target_time}" \
  "grid_mode=${grid_mode}" \
  "deck=${deck}" \
  "restart_file=${restart_file:-fresh}" \
  "nodes=${SLURM_JOB_NUM_NODES}" \
  "ranks=$((SLURM_JOB_NUM_NODES * 8))" \
  >"${segment_root}/invocation.txt"

if [[ -n "${restart_file}" ]]; then
  history_file="${run_root}/${basename}.mhd.hst"
  if [[ ! -f "${history_file}" ]]; then
    printf 'Continuation history is absent: %s\n' "${history_file}" >&2
    exit 1
  fi
  python3 "${repo_root}/tst/publication/prepare_q011_restart_history_v1.py" \
    --restart "${restart_file}" \
    --history "${history_file}" \
    --archive "${segment_root}/${basename}.mhd.hst.pre-continuation" \
    --receipt "${segment_root}/history_preparation.json" \
    >"${segment_root}/history_preparation.stdout"
fi

athena_args=(
  -d "${run_root}"
  "job/basename=${basename}"
  "time/tlim=${target_time}"
  "time/nlim=200000"
  "time/ndiag=500"
  "particles/pic_load_balance_cost_per_particle=0.0"
  "problem/ps_feedback_diag_dcycle=5000"
)
case "${grid_mode}" in
  static_dx3|static_dx6)
    athena_args+=("problem/ps_enable_curvature_amr=false")
    ;;
  uniform_dx6)
    athena_args+=(
      "mesh/nx1=8000" "mesh/nx2=520"
      "mesh_refinement/refinement=none" "mesh_refinement/num_levels=1"
      "problem/ps_enable_curvature_amr=false"
    )
    ;;
  uniform_dx3)
    athena_args+=(
      "mesh/nx1=16000" "mesh/nx2=1040"
      "mesh_refinement/refinement=none" "mesh_refinement/num_levels=1"
      "problem/ps_enable_curvature_amr=false"
    )
    ;;
  uniform_dx12_stress)
    athena_args+=(
      "mesh_refinement/refinement=none" "mesh_refinement/num_levels=1"
      "problem/ps_enable_curvature_amr=false"
      "particles/deposit_qscale=0.0144"
    )
    ;;
  *)
    printf 'Unsupported Q011_GRID_MODE: %s\n' "${grid_mode}" >&2
    exit 2
    ;;
esac
if [[ -n "${restart_file}" ]]; then
  athena_args=(-r "${restart_file}" "${athena_args[@]}")
else
  athena_args=(-i "${deck}" "${athena_args[@]}")
fi

printf '%q ' "${executable}" "${athena_args[@]}" >"${segment_root}/athena_command.sh"
printf '\n' >>"${segment_root}/athena_command.sh"

# The pgen's optional raw escape streams use relative paths. Isolate them by
# segment/run instead of allowing concurrent runs to write into the source tree.
cd "${run_root}"
/usr/bin/time -p -o "${segment_root}/run.time" \
  srun -N"${SLURM_JOB_NUM_NODES}" -n"$((SLURM_JOB_NUM_NODES * 8))" \
    --ntasks-per-node=8 --cpus-per-task=1 --gpus-per-task=1 --gpu-bind=closest \
    "${executable}" "${athena_args[@]}" \
    >"${segment_root}/athena.stdout" \
    2>"${segment_root}/athena.stderr"

segment_state=true
