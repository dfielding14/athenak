#!/usr/bin/env bash

# Group one simulation's explicit snapshot manifest by required node count and
# submit sparse Slurm arrays. Default behavior only prints the sbatch commands.

set -euo pipefail

JOBS_DIR=${JOBS_DIR:-"/ccs/home/dfielding/athenak-gotham-pdf-rebuild/tools/gotham_pdf_rebuild/jobs"}
source "${JOBS_DIR}/common.sh"

SIMULATION=${SIMULATION:-""}
SUBMIT_MODE=${SUBMIT_MODE:-"print"}
ARRAY_CONCURRENCY=${ARRAY_CONCURRENCY:-"1"}
PRODUCTS=${PRODUCTS:-"science"}
CHUNK_BLOCKS=${CHUNK_BLOCKS:-"16"}
CONFIRM_FRONTIER_PRODUCTION=${CONFIRM_FRONTIER_PRODUCTION:-""}
WORKER_SCRIPT=${WORKER_SCRIPT:-"${JOBS_DIR}/frontier_production_array.sbatch"}
EXECUTABLE=${EXECUTABLE:-"${REPO_ROOT}/build-gotham-pdf-frontier/gotham_pdf_rebuild"}
RELEASE_IDENTITY=${RELEASE_IDENTITY:-""}
GOLDEN_ORIGINAL_PASS=${GOLDEN_ORIGINAL_PASS:-""}
GOLDEN_SCIENCE_PASS=${GOLDEN_SCIENCE_PASS:-""}
GOLDEN_ALL_PASS=${GOLDEN_ALL_PASS:-""}
CANARY_256_PASS=${CANARY_256_PASS:-""}
CANARY_NO_CONTROL_PASS=${CANARY_NO_CONTROL_PASS:-""}

case "$SIMULATION" in
    res_4pc|res_4pc_highmdot|res_8pc|res_8pc_lowmdot|res_8pc_highmdot)
        ;;
    *)
        fail "SIMULATION must name one of the five production simulations"
        ;;
esac
case "$SUBMIT_MODE" in
    print|submit)
        ;;
    *)
        fail "SUBMIT_MODE must be print or submit"
        ;;
esac
require_run_mode
require_products "$PRODUCTS"
require_positive_integer "$ARRAY_CONCURRENCY" "ARRAY_CONCURRENCY"
require_positive_integer "$CHUNK_BLOCKS" "CHUNK_BLOCKS"
require_file "$WORKER_SCRIPT" "Frontier production worker"

MANIFEST_PATH=${MANIFEST_PATH:-"${JOBS_DIR}/manifests/${SIMULATION}.tsv"}
require_file "$MANIFEST_PATH" "simulation manifest"

expected_header=$'simulation\tphase\tsequence\tinput_relative\texpected_shards\texpected_nodes\toutput_number\toutput_time\tcontrol_pdf_sequence'
actual_header=$(head -n 1 "$MANIFEST_PATH")
[[ "$actual_header" == "$expected_header" ]] || fail "unexpected manifest header: ${MANIFEST_PATH}"

if [[ "$RUN_MODE" == "execute" ]]; then
    [[ "$CONFIRM_FRONTIER_PRODUCTION" == "GOTHAM_PDF_REBUILD_PRODUCTION" ]] || {
        fail "production execution requires CONFIRM_FRONTIER_PRODUCTION=GOTHAM_PDF_REBUILD_PRODUCTION"
    }
    verify_release_pass "$RELEASE_IDENTITY" "$EXECUTABLE" "$GOLDEN_ORIGINAL_PASS" "golden-original"
    verify_release_pass "$RELEASE_IDENTITY" "$EXECUTABLE" "$GOLDEN_SCIENCE_PASS" "golden-science"
    verify_release_pass "$RELEASE_IDENTITY" "$EXECUTABLE" "$GOLDEN_ALL_PASS" "golden-all"
    verify_release_pass "$RELEASE_IDENTITY" "$EXECUTABLE" "$CANARY_256_PASS" "canary-256"
    verify_release_pass "$RELEASE_IDENTITY" "$EXECUTABLE" "$CANARY_NO_CONTROL_PASS" "canary-no-control"
fi
if [[ "$SUBMIT_MODE" == "submit" ]]; then
    require_slurm_cluster frontier
fi

declare -A indices_by_nodes=()
index=0
while IFS=$'\t' read -r simulation phase sequence input_relative expected_shards expected_nodes \
    output_number output_time control_pdf_sequence; do
    [[ "$simulation" == "$SIMULATION" ]] || fail "manifest row ${index} names ${simulation}, expected ${SIMULATION}"
    require_positive_integer "$expected_nodes" "manifest expected_nodes at row ${index}"
    if [[ -n "${indices_by_nodes[$expected_nodes]:-}" ]]; then
        indices_by_nodes[$expected_nodes]+=",${index}"
    else
        indices_by_nodes[$expected_nodes]="${index}"
    fi
    index=$((index + 1))
done < <(tail -n +2 "$MANIFEST_PATH")
(( index > 0 )) || fail "manifest has no snapshots: ${MANIFEST_PATH}"

log "Prepared ${index} ${SIMULATION} snapshots from ${MANIFEST_PATH}"
for nodes in $(printf '%s\n' "${!indices_by_nodes[@]}" | sort -n); do
    array_spec="${indices_by_nodes[$nodes]}%${ARRAY_CONCURRENCY}"
    export_spec="ALL,JOBS_DIR=${JOBS_DIR},MANIFEST_PATH=${MANIFEST_PATH},RUN_MODE=${RUN_MODE},PRODUCTS=${PRODUCTS},CHUNK_BLOCKS=${CHUNK_BLOCKS},OUTPUT_ROOT=${OUTPUT_ROOT},ARCHIVE_ROOT=${ARCHIVE_ROOT},REPO_ROOT=${REPO_ROOT},ALLOW_EXISTING_OUTPUT=${ALLOW_EXISTING_OUTPUT},CONFIRM_FRONTIER_PRODUCTION=${CONFIRM_FRONTIER_PRODUCTION},EXECUTABLE=${EXECUTABLE},RELEASE_IDENTITY=${RELEASE_IDENTITY},GOLDEN_ORIGINAL_PASS=${GOLDEN_ORIGINAL_PASS},GOLDEN_SCIENCE_PASS=${GOLDEN_SCIENCE_PASS},GOLDEN_ALL_PASS=${GOLDEN_ALL_PASS},CANARY_256_PASS=${CANARY_256_PASS},CANARY_NO_CONTROL_PASS=${CANARY_NO_CONTROL_PASS}"
    command=(
        sbatch
        --nodes="$nodes"
        --array="$array_spec"
        --export="$export_spec"
        "$WORKER_SCRIPT"
    )
    log "Node group ${nodes}: array=${array_spec}"
    print_command "${command[@]}"
    if [[ "$SUBMIT_MODE" == "submit" ]]; then
        "${command[@]}"
    fi
done
