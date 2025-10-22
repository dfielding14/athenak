#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOCS_DIR="${ROOT_DIR}/docs"
BUILD_DIR="${DOCS_DIR}/build/html"
PORT="${PORT:-8000}"

echo "Building documentation (make -C ${DOCS_DIR} html)…"
make -C "${DOCS_DIR}" html

if [[ ! -d "${BUILD_DIR}" ]]; then
  echo "error: expected build output at ${BUILD_DIR}" >&2
  exit 1
fi

cd "${BUILD_DIR}"
echo "Serving ${BUILD_DIR} at http://127.0.0.1:${PORT} (Ctrl-C to stop)"
python3 -m http.server "${PORT}"
