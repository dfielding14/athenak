#!/usr/bin/env bash

# Shared fail-closed helpers for GOTHAM PDF reconstruction jobs.

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-"/ccs/home/dfielding/athenak-gotham-pdf-rebuild"}
ARCHIVE_ROOT=${ARCHIVE_ROOT:-"/lustre/orion/ast207/proj-shared/brent/gotham"}
OUTPUT_ROOT=${OUTPUT_ROOT:-"/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild"}
RELEASE_ARTIFACT_DIR=${RELEASE_ARTIFACT_DIR:-"${OUTPUT_ROOT}/release"}
RUN_MODE=${RUN_MODE:-"plan"}
ALLOW_EXISTING_OUTPUT=${ALLOW_EXISTING_OUTPUT:-"NO"}
if [[ -z "${PYTHON:-}" ]]; then
    if [[ -x /opt/cray/pe/python/3.11.7/bin/python3 ]]; then
        PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
    else
        PYTHON=python3
    fi
fi

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

fail() {
    log "ERROR: $*" >&2
    exit 1
}

require_file() {
    local path="$1"
    local description="$2"
    [[ -f "$path" ]] || fail "${description} not found: ${path}"
}

require_directory() {
    local path="$1"
    local description="$2"
    [[ -d "$path" ]] || fail "${description} not found: ${path}"
}

require_executable() {
    local path="$1"
    local description="$2"
    [[ -x "$path" ]] || fail "${description} is not executable: ${path}"
}

verify_release_identity() {
    local identity="$1"
    local executable="$2"
    require_file "$identity" "frozen Frontier release identity"
    require_file "${identity}.sha256" "frozen Frontier release identity checksum"
    "$PYTHON" "${JOBS_DIR}/release_gate.py" verify-identity \
        --identity "$identity" \
        --repo-root "$REPO_ROOT" \
        --executable "$executable"
}

verify_release_pass() {
    local identity="$1"
    local executable="$2"
    local report="$3"
    local expected_kind="$4"
    require_file "$report" "${expected_kind} pass report"
    require_file "${report}.sha256" "${expected_kind} pass report checksum"
    "$PYTHON" "${JOBS_DIR}/release_gate.py" verify-pass \
        --identity "$identity" \
        --repo-root "$REPO_ROOT" \
        --executable "$executable" \
        --pass-report "$report" \
        --expected-kind "$expected_kind"
}

require_unsigned_integer() {
    local value="$1"
    local description="$2"
    [[ "$value" =~ ^[0-9]+$ ]] || fail "${description} must be an unsigned integer: ${value}"
}

require_positive_integer() {
    local value="$1"
    local description="$2"
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || fail "${description} must be a positive integer: ${value}"
}

require_run_mode() {
    case "$RUN_MODE" in
        plan|execute)
            ;;
        *)
            fail "RUN_MODE must be plan or execute, not: ${RUN_MODE}"
            ;;
    esac
}

require_products() {
    case "$1" in
        original|science|all)
            ;;
        *)
            fail "PRODUCTS must be original, science, or all, not: $1"
            ;;
    esac
}

slurm_cluster_name() {
    scontrol show config 2>/dev/null |
        awk -F= '/^ClusterName/ && !found {gsub(/[[:space:]]/, "", $2); print $2; found=1}'
}

require_slurm_cluster() {
    local expected="$1"
    local actual
    actual=$(slurm_cluster_name)
    [[ -n "$actual" ]] || fail "unable to determine Slurm cluster name"
    [[ "$actual" == "$expected" ]] || fail "this operation requires Slurm cluster ${expected}, found ${actual}"
}

print_command() {
    printf 'COMMAND:'
    printf ' %q' "$@"
    printf '\n'
}

count_source_shards() {
    local input_dir="$1"
    local sequence="$2"
    find -L "$input_dir" -mindepth 1 -maxdepth 2 -type f \
        -name "gotham.hydro_w.${sequence}.bin" -print | wc -l
}

require_source_shards() {
    local input_dir="$1"
    local sequence="$2"
    local expected="$3"
    local found

    require_directory "$input_dir" "source input directory"
    require_positive_integer "$expected" "EXPECTED_SOURCE_SHARDS"
    found=$(count_source_shards "$input_dir" "$sequence")
    [[ "$found" == "$expected" ]] || {
        fail "source shard count mismatch for ${input_dir} sequence ${sequence}: expected ${expected}, found ${found}"
    }
    log "Verified ${found} source shards for ${input_dir} sequence ${sequence}"
}

prepare_output_directory() {
    local output_dir="$1"
    if [[ "$RUN_MODE" == "plan" ]]; then
        log "Plan mode: output directory will not be created: ${output_dir}"
        return
    fi

    if [[ -e "$output_dir" ]]; then
        [[ -d "$output_dir" ]] || fail "output target exists and is not a directory: ${output_dir}"
        [[ "$ALLOW_EXISTING_OUTPUT" == "YES" ]] || {
            fail "refusing existing output directory; set ALLOW_EXISTING_OUTPUT=YES only for an intentional rerun: ${output_dir}"
        }
        return
    fi

    mkdir -p "$(dirname "$output_dir")"
    mkdir "$output_dir" || fail "unable to claim output directory; another job may have won the race: ${output_dir}"
}

validate_rebuild_manifest() {
    local output_dir="$1"
    local expected_available="$2"
    local expected_processed="$3"

    "$PYTHON" - "$output_dir" "$expected_available" "$expected_processed" <<'PY'
import json
import math
import pathlib
import sys

output_dir = pathlib.Path(sys.argv[1])
expected_available = int(sys.argv[2])
expected_processed = int(sys.argv[3])
manifest_path = output_dir / "rebuild_manifest.json"
if not manifest_path.is_file():
    raise SystemExit(f"missing rebuild manifest: {manifest_path}")

manifest = json.loads(manifest_path.read_text())
checks = {
    "shards_available": expected_available,
    "shards_processed": expected_processed,
}
for key, expected in checks.items():
    actual = manifest.get(key)
    if actual != expected:
        raise SystemExit(f"{key}: expected {expected}, found {actual}")

for key in ("meshblocks_processed", "cells_processed", "payload_bytes_read"):
    value = manifest.get(key)
    if not isinstance(value, int) or value <= 0:
        raise SystemExit(f"{key} is not positive: {value!r}")

products = manifest.get("products")
if not isinstance(products, list) or not products:
    raise SystemExit("manifest has no products")

output_number = manifest.get("output_number")
for product in products:
    product_id = product.get("id")
    total = product.get("sum")
    if not isinstance(total, (int, float)) or not math.isfinite(total):
        raise SystemExit(f"non-finite product sum for {product_id}: {total!r}")
    product_dir = output_dir / product_id
    header = product_dir / "gotham.header.pdf"
    payload = product_dir / f"gotham.{output_number}.pdf"
    if not header.is_file() or not payload.is_file():
        raise SystemExit(f"missing output files for {product_id}")

print(
    f"validated {manifest_path}: "
    f"{manifest['shards_processed']} shards, "
    f"{manifest['cells_processed']} cells, {len(products)} products"
)
PY
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    fail "common.sh is a library and must be sourced"
fi
