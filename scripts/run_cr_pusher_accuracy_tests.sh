#!/usr/bin/env bash
# Run the CR pusher accuracy-test matrix on Frontier GPUs.
#
# This script is both sbatch-able and runnable inside an existing allocation.
#
#SBATCH -A AST207
#SBATCH -J CR_PUSH_ACC
#SBATCH -o /lustre/orion/ast207/proj-shared/dfielding/AMR/particles/testing/%x.%j.out
#SBATCH -t 02:00:00
#SBATCH -p batch
#SBATCH --qos=debug
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1

set -euo pipefail

REPO=${AMR_CR_ACCURACY_REPO:-/ccs/home/dfielding/athenak-cr-tracers-followup-architecture}
EXE=${AMR_CR_ACCURACY_EXE:-${REPO}/build-frontier/src/athena}
INPUT_DIR=${REPO}/inputs/particles
ANALYZER=${REPO}/scripts/analyze_cr_pusher_accuracy.py
TESTING_ROOT=${AMR_CR_ACCURACY_TESTING_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/AMR/particles/testing}
STAMP=${AMR_CR_ACCURACY_STAMP:-$(date +%Y%m%d_%H%M%S)}
RUN_ROOT=${AMR_CR_ACCURACY_RUN_ROOT:-${TESTING_ROOT}/cr_pusher_accuracy_${STAMP}}
MANIFEST=${RUN_ROOT}/cases.jsonl

mkdir -p "${RUN_ROOT}"
: >"${MANIFEST}"

if command -v module >/dev/null 2>&1; then
  module restore
  module load cpe/24.07 PrgEnv-amd cray-mpich/8.1.30 craype-accel-amd-gfx90a amd/6.2.0 rocm/6.2.0
  module unload darshan-runtime || true
fi

export MPICH_GPU_SUPPORT_ENABLED=${MPICH_GPU_SUPPORT_ENABLED:-1}
export MPICH_SMP_SINGLE_COPY_MODE=${MPICH_SMP_SINGLE_COPY_MODE:-NONE}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export HSA_ENABLE_DEBUG=0
ulimit -c 0 || true

test -x "${EXE}"
test -f "${ANALYZER}"

git -C "${REPO}" rev-parse HEAD >"${RUN_ROOT}/source.commit"
git -C "${REPO}" status --short >"${RUN_ROOT}/source.status"
sha256sum "${EXE}" "${ANALYZER}" \
  "${INPUT_DIR}/cr_pusher_uniform_b.athinput" \
  "${INPUT_DIR}/cr_pusher_smooth_divb0.athinput" \
  "${INPUT_DIR}/cr_pusher_boundary_crossing.athinput" \
  >"${RUN_ROOT}/bindings.sha256"
env | sort >"${RUN_ROOT}/environment.txt"

append_manifest() {
  local case_name=$1
  local test_name=$2
  local mode=$3
  local gyro_fraction=$4
  local ranks=$5
  local input_file=$6
  local run_dir=$7
  local ntrack=$8
  local nspecies=$9
  local min_mass=${10}
  local mass_log_spacing=${11}
  local b0x=${12}
  local b0y=${13}
  local b0z=${14}
  local b_profile=${15}
  local bgrad=${16}
  local bamp=${17}
  local bwave=${18}

  python3 - "${MANIFEST}" \
    "${case_name}" "${test_name}" "${mode}" "${gyro_fraction}" "${ranks}" \
    "${input_file}" "${run_dir}" "${ntrack}" "${nspecies}" "${min_mass}" \
    "${mass_log_spacing}" "${b0x}" "${b0y}" "${b0z}" "${b_profile}" \
    "${bgrad}" "${bamp}" "${bwave}" <<'PY'
import json
import sys
from pathlib import Path

(
    manifest, case_name, test_name, mode, gyro_fraction, ranks, input_file,
    run_dir, ntrack, nspecies, min_mass, mass_log_spacing, b0x, b0y, b0z,
    b_profile, bgrad, bamp, bwave,
) = sys.argv[1:]
record = {
    "case": case_name,
    "test": test_name,
    "mode": mode,
    "gyro_fraction": float(gyro_fraction),
    "ranks": int(ranks),
    "input": input_file,
    "run_dir": run_dir,
    "ntrack": int(ntrack),
    "nspecies": int(nspecies),
    "min_mass": float(min_mass),
    "mass_log_spacing": float(mass_log_spacing),
    "B0x": float(b0x),
    "B0y": float(b0y),
    "B0z": float(b0z),
    "B_profile": b_profile,
    "Bgrad": float(bgrad),
    "Bamp": float(bamp),
    "Bwave_number": float(bwave),
    "x1min": -0.5,
    "x1max": 0.5,
    "x2min": -0.5,
    "x2max": 0.5,
    "x3min": -0.5,
    "x3max": 0.5,
}
with Path(manifest).open("a") as handle:
    handle.write(json.dumps(record, sort_keys=True) + "\n")
Path(run_dir, "case.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
PY
}

run_case() {
  local case_name=$1
  local test_name=$2
  local mode=$3
  local gyro_fraction=$4
  local ranks=$5
  local input_file=$6
  local ntrack=$7
  local nspecies=$8
  local min_mass=$9
  local mass_log_spacing=${10}
  local b0x=${11}
  local b0y=${12}
  local b0z=${13}
  local b_profile=${14}
  local bgrad=${15}
  local bamp=${16}
  local bwave=${17}
  local per_particle=${18}
  shift 18
  local overrides=("$@")

  local run_dir=${RUN_ROOT}/${case_name}
  if [[ -e "${run_dir}/.complete" ]]; then
    printf 'Skipping complete case: %s\n' "${case_name}"
    append_manifest "${case_name}" "${test_name}" "${mode}" "${gyro_fraction}" \
      "${ranks}" "${input_file}" "${run_dir}" "${ntrack}" "${nspecies}" \
      "${min_mass}" "${mass_log_spacing}" "${b0x}" "${b0y}" "${b0z}" \
      "${b_profile}" "${bgrad}" "${bamp}" "${bwave}"
    return
  fi
  if [[ -e "${run_dir}" && "${AMR_CR_ACCURACY_OVERWRITE:-0}" != "1" ]]; then
    printf 'Refusing to overwrite existing incomplete case: %s\n' "${run_dir}" >&2
    exit 1
  fi
  rm -rf "${run_dir}"
  mkdir -p "${run_dir}"
  append_manifest "${case_name}" "${test_name}" "${mode}" "${gyro_fraction}" \
    "${ranks}" "${input_file}" "${run_dir}" "${ntrack}" "${nspecies}" \
    "${min_mass}" "${mass_log_spacing}" "${b0x}" "${b0y}" "${b0z}" \
    "${b_profile}" "${bgrad}" "${bamp}" "${bwave}"

  local command=(
    "${EXE}" -i "${input_file}" -d "${run_dir}"
    "job/basename=${case_name}"
    "particles/subcycle=true"
    "particles/subcycle_gyro_fraction=${gyro_fraction}"
    "particles/subcycle_per_particle_gyro=${per_particle}"
    "particles/exchange_gyro_only_substeps=false"
    "particles/update_global_counts_each_exchange=false"
  )
  command+=("${overrides[@]}")
  printf '%q ' "${command[@]}" >"${run_dir}/command.txt"
  printf '\n' >>"${run_dir}/command.txt"

  local launcher=()
  if [[ -n "${SLURM_JOB_ID:-}" || "${AMR_CR_ACCURACY_USE_SRUN:-0}" == "1" ]]; then
    launcher=(
      srun -N1 -n"${ranks}" --ntasks-per-node="${ranks}"
      --cpus-per-task=1 --gpus-per-task=1 --gpu-bind=closest
      --kill-on-bad-exit=1
    )
  fi

  printf 'Running %-34s ranks=%s mode=%s gyro=%s\n' \
    "${case_name}" "${ranks}" "${mode}" "${gyro_fraction}"
  /usr/bin/time -p -o "${run_dir}/run.time" \
    "${launcher[@]}" "${command[@]}" \
    >"${run_dir}/stdout.txt" 2>"${run_dir}/stderr.txt"
  touch "${run_dir}/.complete"
}

UNIFORM_INPUT=${INPUT_DIR}/cr_pusher_uniform_b.athinput
SMOOTH_INPUT=${INPUT_DIR}/cr_pusher_smooth_divb0.athinput
BOUNDARY_INPUT=${INPUT_DIR}/cr_pusher_boundary_crossing.athinput

# Test 1: uniform B, analytic helix.
run_case uniform_global_g005 uniform_b global 0.05 1 "${UNIFORM_INPUT}" \
  8 4 0.25 2.0 0.0 0.0 1.0 uniform 1.0 0.05 1.0 false
run_case uniform_pp_g005 uniform_b per_particle 0.05 1 "${UNIFORM_INPUT}" \
  8 4 0.25 2.0 0.0 0.0 1.0 uniform 1.0 0.05 1.0 true
run_case uniform_pp_g010 uniform_b per_particle 0.10 1 "${UNIFORM_INPUT}" \
  8 4 0.25 2.0 0.0 0.0 1.0 uniform 1.0 0.05 1.0 true
run_case uniform_pp_g020 uniform_b per_particle 0.20 1 "${UNIFORM_INPUT}" \
  8 4 0.25 2.0 0.0 0.0 1.0 uniform 1.0 0.05 1.0 true
run_case uniform_pp_g030 uniform_b per_particle 0.30 1 "${UNIFORM_INPUT}" \
  8 4 0.25 2.0 0.0 0.0 1.0 uniform 1.0 0.05 1.0 true

# Test 2: smooth divergence-free B, strict numerical reference plus sweep.
run_case smooth_ref_g001 smooth_divb0 reference 0.01 1 "${SMOOTH_INPUT}" \
  8 4 0.25 2.0 0.1 -0.2 1.0 sinusoidal_divb_free 1.0 0.015 1.0 true
run_case smooth_pp_g002 smooth_divb0 per_particle 0.02 1 "${SMOOTH_INPUT}" \
  8 4 0.25 2.0 0.1 -0.2 1.0 sinusoidal_divb_free 1.0 0.015 1.0 true
run_case smooth_pp_g005 smooth_divb0 per_particle 0.05 1 "${SMOOTH_INPUT}" \
  8 4 0.25 2.0 0.1 -0.2 1.0 sinusoidal_divb_free 1.0 0.015 1.0 true
run_case smooth_pp_g010 smooth_divb0 per_particle 0.10 1 "${SMOOTH_INPUT}" \
  8 4 0.25 2.0 0.1 -0.2 1.0 sinusoidal_divb_free 1.0 0.015 1.0 true
run_case smooth_pp_g020 smooth_divb0 per_particle 0.20 1 "${SMOOTH_INPUT}" \
  8 4 0.25 2.0 0.1 -0.2 1.0 sinusoidal_divb_free 1.0 0.015 1.0 true
run_case smooth_pp_g030 smooth_divb0 per_particle 0.30 1 "${SMOOTH_INPUT}" \
  8 4 0.25 2.0 0.1 -0.2 1.0 sinusoidal_divb_free 1.0 0.015 1.0 true
run_case smooth_global_g005 smooth_divb0 global 0.05 1 "${SMOOTH_INPUT}" \
  8 4 0.25 2.0 0.1 -0.2 1.0 sinusoidal_divb_free 1.0 0.015 1.0 false

# Test 3: boundary/exchange stress.  A single-block reference supplies the
# low-crossing trajectory, while the 8-rank cases force frequent exchanges.
run_case boundary_ref_singleblock boundary_reference reference 0.05 1 "${BOUNDARY_INPUT}" \
  16 2 1.0 1.0 0.0 0.0 1.0 uniform 1.0 0.05 1.0 false \
  "meshblock/nx1=16" "meshblock/nx2=16" "meshblock/nx3=16" \
  "particles/check_motion_bounds=false"
run_case boundary_global_g005 boundary_crossing global 0.05 8 "${BOUNDARY_INPUT}" \
  16 2 1.0 1.0 0.0 0.0 1.0 uniform 1.0 0.05 1.0 false
run_case boundary_pp_g005 boundary_crossing per_particle 0.05 8 "${BOUNDARY_INPUT}" \
  16 2 1.0 1.0 0.0 0.0 1.0 uniform 1.0 0.05 1.0 true
run_case boundary_pp_g010 boundary_crossing per_particle 0.10 8 "${BOUNDARY_INPUT}" \
  16 2 1.0 1.0 0.0 0.0 1.0 uniform 1.0 0.05 1.0 true

python3 "${ANALYZER}" \
  --manifest "${MANIFEST}" \
  --write-json "${RUN_ROOT}/analysis/summary.json" \
  --write-csv "${RUN_ROOT}/analysis/error_table.csv" \
  --write-md "${RUN_ROOT}/analysis/summary.md" \
  --fail-on-error

ln -sfn "${RUN_ROOT}" "${TESTING_ROOT}/cr_pusher_accuracy_latest"
printf 'CR pusher accuracy run complete: %s\n' "${RUN_ROOT}"

