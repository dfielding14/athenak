#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

BUILD_DIR="${BUILD_DIR:-build-cgl-implementation}"
ATHENA_BIN="${ATHENA_BIN:-${BUILD_DIR}/src/athena}"
OUTPUT_DIR="${OUTPUT_DIR:-${BUILD_DIR}/cgl_lf_validation}"
DATA_DIR="${DATA_DIR:-${OUTPUT_DIR}/data}"
LOG_DIR="${LOG_DIR:-${OUTPUT_DIR}/logs}"
FIGURE_DIR="${FIGURE_DIR:-${OUTPUT_DIR}/figures}"

if [[ ! -x "${ATHENA_BIN}" ]]; then
  echo "AthenaK binary not found at ${ATHENA_BIN}; configuring ${BUILD_DIR}" >&2
  cmake -S . -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
  cmake --build "${BUILD_DIR}" -j "${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"
fi

mkdir -p "${DATA_DIR}" "${LOG_DIR}" "${FIGURE_DIR}"

inputs=(
  inputs/unit_tests/cgl_lf_quant_parallel.athinput
  inputs/unit_tests/cgl_lf_quant_parallel_collisional.athinput
  inputs/unit_tests/cgl_lf_quant_perp.athinput
  inputs/unit_tests/cgl_lf_quant_perp_collisional.athinput
  inputs/unit_tests/cgl_lf_quant_grad_b.athinput
  inputs/unit_tests/cgl_lf_flux_limiter.athinput
  inputs/unit_tests/cgl_lf_limiter_heat_flux_suppression.athinput
  inputs/unit_tests/cgl_lf_limiter_mirror.athinput
  inputs/unit_tests/cgl_lf_limiter_firehose.athinput
  inputs/unit_tests/cgl_lf_field_wave.athinput
  inputs/unit_tests/cgl_lf_paper_oblique_wave.athinput
  inputs/unit_tests/cgl_pure_paper_oblique_wave.athinput
  inputs/unit_tests/cgl_pure_paper_eigen_alfven.athinput
  inputs/unit_tests/cgl_pure_paper_eigen_slow.athinput
  inputs/unit_tests/cgl_pure_paper_eigen_fast.athinput
  inputs/unit_tests/cgl_lf_paper_eigen_alfven.athinput
  inputs/unit_tests/cgl_lf_paper_eigen_slow.athinput
  inputs/unit_tests/cgl_lf_paper_eigen_fast.athinput
)

echo "Running CGL-LF validation tests with ${ATHENA_BIN}"
for input in "${inputs[@]}"; do
  name="$(basename "${input}" .athinput)"
  log="${LOG_DIR}/${name}.log"
  echo "  ${name}"
  "${ATHENA_BIN}" -i "${input}" \
    problem/validation_output=true \
    problem/validation_output_dir="${DATA_DIR}" | tee "${log}"
done

python3 scripts/plot_cgl_lf_validation.py \
  --data-dir "${DATA_DIR}" \
  --figure-dir "${FIGURE_DIR}"

echo "CGL-LF validation complete"
