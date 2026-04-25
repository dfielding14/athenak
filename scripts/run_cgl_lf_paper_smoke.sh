#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

BUILD_DIR="${BUILD_DIR:-build-cgl-lf-paper}"
ATHENA_BIN="${ATHENA_BIN:-${BUILD_DIR}/src/athena}"
RUN_DIR="${RUN_DIR:-/tmp/athenak_cgl_lf_paper_smoke}"
JOBS="${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"

if [[ ! -x "${ATHENA_BIN}" ]]; then
  cmake -S . -B "${BUILD_DIR}" -DPROBLEM=cgl_lf_paper -DCMAKE_BUILD_TYPE=Release
  cmake --build "${BUILD_DIR}" -j "${JOBS}"
fi

rm -rf "${RUN_DIR}"
mkdir -p "${RUN_DIR}"

run_case() {
  local input="$1"
  shift
  local name
  name="$(basename "${input}" .athinput)"
  echo "Running ${name}"
  (
    cd "${RUN_DIR}"
    "${ROOT_DIR}/${ATHENA_BIN}" -i "${ROOT_DIR}/${input}" \
      mesh/nx1=8 mesh/nx2=8 mesh/nx3=16 \
      meshblock/nx1=8 meshblock/nx2=8 meshblock/nx3=16 \
      time/nlim=2 time/tlim=0.02 \
      output1/dt=0.01 output2/dt=10.0 "$@"
  )
}

run_case inputs/cgl_lf_paper/cgl_lf_paper_turb_active.athinput
run_case inputs/cgl_lf_paper/cgl_lf_paper_turb_active.athinput \
  job/basename=cgl_lf_paper_turb_random \
  problem/forcing_mode=random \
  turb_driving/driving_type=2
run_case inputs/cgl_lf_paper/cgl_lf_paper_turb_passive.athinput
run_case inputs/cgl_lf_paper/cgl_lf_paper_turb_limiter_off.athinput
run_case inputs/cgl_lf_paper/cgl_lf_paper_np_mode.athinput
run_case inputs/cgl_lf_paper/cgl_lf_paper_fast_wave.athinput
run_case inputs/cgl_lf_paper/cgl_lf_paper_oblique_iaw.athinput
run_case inputs/cgl_lf_paper/cgl_lf_paper_linear_wave_scan.athinput

python3 "${ROOT_DIR}/scripts/analyze_cgl_lf_paper.py" "${RUN_DIR}" \
  --output-dir "${RUN_DIR}/analysis" \
  --synthetic-wave-test

echo "CGL-LF paper smoke outputs: ${RUN_DIR}"
