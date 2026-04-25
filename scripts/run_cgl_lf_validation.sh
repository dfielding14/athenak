#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

BUILD_DIR="${BUILD_DIR:-build-codex-cgl-lf-quant}"
ATHENA_BIN="${ATHENA_BIN:-${BUILD_DIR}/src/athena}"
DATA_DIR="${DATA_DIR:-docs/figures/data}"
LOG_DIR="${LOG_DIR:-docs/figures/logs}"
FIGURE_DIR="${FIGURE_DIR:-docs/figures}"

if [[ ! -x "${ATHENA_BIN}" ]]; then
  echo "AthenaK binary not found at ${ATHENA_BIN}; building ${BUILD_DIR}" >&2
  cmake --build "${BUILD_DIR}" -j "${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"
fi

mkdir -p "${DATA_DIR}" "${LOG_DIR}" "${FIGURE_DIR}"

inputs=(
  inputs/unit_tests/cgl_lf_quant_parallel.athinput
  inputs/unit_tests/cgl_lf_quant_perp.athinput
  inputs/unit_tests/cgl_lf_quant_grad_b.athinput
  inputs/unit_tests/cgl_lf_flux_limiter.athinput
  inputs/unit_tests/cgl_lf_limiter_mirror.athinput
  inputs/unit_tests/cgl_lf_limiter_firehose.athinput
  inputs/unit_tests/cgl_lf_field_wave.athinput
  inputs/unit_tests/cgl_lf_paper_oblique_wave.athinput
  inputs/unit_tests/cgl_pure_paper_oblique_wave.athinput
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

if command -v pdflatex >/dev/null 2>&1; then
  (cd docs && pdflatex -interaction=nonstopmode cgl_lf_validation.tex >/dev/null)
  (cd docs && pdflatex -interaction=nonstopmode cgl_lf_validation.tex >/dev/null)
  echo "Updated docs/cgl_lf_validation.pdf"
else
  echo "pdflatex not found; skipped PDF rebuild" >&2
fi

echo "CGL-LF validation complete"
