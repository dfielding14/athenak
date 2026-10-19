#!/usr/bin/env bash
#SBATCH -A AST207
#SBATCH -J gotham_pdf_frontier
#SBATCH -o /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_frontier_%A_%a.out
#SBATCH -e /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_frontier_%A_%a.err
#SBATCH -t 06:00:00
#SBATCH -p batch
#SBATCH -q normal
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7

set -euo pipefail

JOBS_DIR=${JOBS_DIR:-"/ccs/home/dfielding/athenak-gotham-pdf-rebuild/tools/gotham_pdf_rebuild/jobs"}
source "${JOBS_DIR}/common.sh"

EXECUTABLE=${EXECUTABLE:-"${REPO_ROOT}/build-gotham-pdf-frontier/gotham_pdf_rebuild"}
MANIFEST_PATH=${MANIFEST_PATH:-""}
PRODUCTS=${PRODUCTS:-"science"}
CHUNK_BLOCKS=${CHUNK_BLOCKS:-"16"}
CONFIRM_FRONTIER_PRODUCTION=${CONFIRM_FRONTIER_PRODUCTION:-""}
RELEASE_IDENTITY=${RELEASE_IDENTITY:-""}
GOLDEN_ORIGINAL_PASS=${GOLDEN_ORIGINAL_PASS:-""}
GOLDEN_SCIENCE_PASS=${GOLDEN_SCIENCE_PASS:-""}
GOLDEN_ALL_PASS=${GOLDEN_ALL_PASS:-""}
CANARY_256_PASS=${CANARY_256_PASS:-""}
CANARY_NO_CONTROL_PASS=${CANARY_NO_CONTROL_PASS:-""}

require_run_mode
require_products "$PRODUCTS"
require_positive_integer "$CHUNK_BLOCKS" "CHUNK_BLOCKS"
require_file "$MANIFEST_PATH" "Frontier snapshot manifest"
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || fail "submit this worker as a Slurm array"
require_unsigned_integer "$SLURM_ARRAY_TASK_ID" "SLURM_ARRAY_TASK_ID"
require_positive_integer "${SLURM_JOB_NUM_NODES:-0}" "SLURM_JOB_NUM_NODES"
require_slurm_cluster frontier

manifest_line=$(awk -F $'\t' -v row="$((SLURM_ARRAY_TASK_ID + 2))" 'NR == row {print; exit}' "$MANIFEST_PATH")
[[ -n "$manifest_line" ]] || fail "manifest has no data row for array index ${SLURM_ARRAY_TASK_ID}"
IFS=$'\t' read -r simulation phase sequence input_relative expected_shards expected_nodes \
    output_number output_time control_pdf_sequence <<< "$manifest_line"

[[ "$simulation" =~ ^res_(4|8)pc(_highmdot|_lowmdot)?$ ]] || fail "invalid simulation in manifest: ${simulation}"
[[ "$phase" =~ ^phase[0-9a-z_]+$ ]] || fail "invalid phase in manifest: ${phase}"
[[ "$sequence" =~ ^[0-9]{5}$ ]] || fail "invalid sequence in manifest: ${sequence}"
[[ "$input_relative" != /* ]] || fail "manifest input path must be relative: ${input_relative}"
[[ "$input_relative" != ../* && "$input_relative" != *"/../"* && "$input_relative" != */.. ]] || {
    fail "manifest input path escapes its root: ${input_relative}"
}
require_positive_integer "$expected_shards" "manifest expected_shards"
require_positive_integer "$expected_nodes" "manifest expected_nodes"
[[ "$output_number" =~ ^[0-9]{5}$ ]] || fail "invalid output_number in manifest: ${output_number}"
[[ "$output_time" =~ ^-?[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$ ]] || {
    fail "invalid output_time in manifest: ${output_time}"
}
[[ "$SLURM_JOB_NUM_NODES" == "$expected_nodes" ]] || {
    fail "allocation mismatch: manifest requires ${expected_nodes} nodes, job has ${SLURM_JOB_NUM_NODES}"
}
if [[ "$RUN_MODE" == "execute" ]]; then
    [[ "$CONFIRM_FRONTIER_PRODUCTION" == "GOTHAM_PDF_REBUILD_PRODUCTION" ]] || {
        fail "production execution requires CONFIRM_FRONTIER_PRODUCTION=GOTHAM_PDF_REBUILD_PRODUCTION"
    }
fi

case "$input_relative" in
    "${simulation}/"*)
        INPUT_DIR="${ARCHIVE_ROOT}/${input_relative}"
        ;;
    source_views/*)
        INPUT_DIR="${OUTPUT_ROOT}/${input_relative}"
        ;;
    *)
        fail "manifest input path is not an archive or curated source view: ${input_relative}"
        ;;
esac
OUTPUT_DIR="${OUTPUT_ROOT}/production/${simulation}/${phase}/${sequence}"
NTASKS=$((expected_nodes * 8))

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module load cmake
module unload darshan-runtime

export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH:-}:${LD_LIBRARY_PATH:-}"
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export FI_CXI_DEFAULT_CQ_SIZE=131072
export FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software
export HSA_XNACK=1

require_executable "$EXECUTABLE" "Frontier HIP reducer"
if [[ "$RUN_MODE" == "execute" ]]; then
    verify_release_pass "$RELEASE_IDENTITY" "$EXECUTABLE" "$GOLDEN_ORIGINAL_PASS" "golden-original"
    verify_release_pass "$RELEASE_IDENTITY" "$EXECUTABLE" "$GOLDEN_SCIENCE_PASS" "golden-science"
    verify_release_pass "$RELEASE_IDENTITY" "$EXECUTABLE" "$GOLDEN_ALL_PASS" "golden-all"
    verify_release_pass "$RELEASE_IDENTITY" "$EXECUTABLE" "$CANARY_256_PASS" "canary-256"
    verify_release_pass "$RELEASE_IDENTITY" "$EXECUTABLE" "$CANARY_NO_CONTROL_PASS" "canary-no-control"
fi
require_source_shards "$INPUT_DIR" "$sequence" "$expected_shards"
prepare_output_directory "$OUTPUT_DIR"

command=(
    srun --nodes="$expected_nodes" --ntasks="$NTASKS" --ntasks-per-node=8
    --cpus-per-task=7 --gpus-per-node=8 --gpu-bind=closest
    "$EXECUTABLE"
    --input-dir "$INPUT_DIR"
    --sequence "$sequence"
    --output-number "$output_number"
    --output-time "$output_time"
    --expected-shards "$expected_shards"
    --output-dir "$OUTPUT_DIR"
    --products "$PRODUCTS"
    --chunk-blocks "$CHUNK_BLOCKS"
)
[[ "$RUN_MODE" == "plan" ]] && command+=(--dry-run)

log "Frontier array index=${SLURM_ARRAY_TASK_ID}, RUN_MODE=${RUN_MODE}"
log "Snapshot=${simulation}/${phase}/${sequence}, output=${output_number}@${output_time}, control=${control_pdf_sequence}, shards=${expected_shards}, nodes=${expected_nodes}"
print_command "${command[@]}"
"${command[@]}"

if [[ "$RUN_MODE" == "execute" ]]; then
    validate_rebuild_manifest "$OUTPUT_DIR" "$expected_shards" "$expected_shards"
fi
