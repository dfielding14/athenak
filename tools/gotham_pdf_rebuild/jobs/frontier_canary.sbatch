#!/usr/bin/env bash
#SBATCH -A AST207
#SBATCH -J gotham_pdf_canary
#SBATCH -o /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_canary_%j.out
#SBATCH -e /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_canary_%j.err
#SBATCH -t 06:00:00
#SBATCH -p batch
#SBATCH -q normal
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7

# Explicit post-golden release canaries. Submit with --nodes=256 or --nodes=64.

set -euo pipefail

JOBS_DIR=${JOBS_DIR:-"/ccs/home/dfielding/athenak-gotham-pdf-rebuild/tools/gotham_pdf_rebuild/jobs"}
source "${JOBS_DIR}/common.sh"

EXECUTABLE=${EXECUTABLE:-"${REPO_ROOT}/build-gotham-pdf-frontier/gotham_pdf_rebuild"}
PRODUCTS=${PRODUCTS:-"all"}
CHUNK_BLOCKS=${CHUNK_BLOCKS:-"16"}
CANARY_KIND=${CANARY_KIND:-""}
CONFIRM_FRONTIER_CANARY=${CONFIRM_FRONTIER_CANARY:-""}
RELEASE_IDENTITY=${RELEASE_IDENTITY:-""}
BUDGET_FILE=${BUDGET_FILE:-"${JOBS_DIR}/performance_budgets.json"}
VALIDATOR=${VALIDATOR:-"${REPO_ROOT}/tools/gotham_pdf_rebuild/validate_real_8pc.py"}

case "$CANARY_KIND" in
    canary-256)
        SIMULATION="res_4pc"
        PHASE="phase2a"
        SEQUENCE="00003"
        OUTPUT_NUMBER="00012"
        OUTPUT_TIME="10.020026243958727"
        EXPECTED_SOURCE_SHARDS="2048"
        EXPECTED_NODES="256"
        CONTROL_MODE="archived"
        ;;
    canary-no-control)
        SIMULATION="res_8pc"
        PHASE="phase2_uniform"
        SEQUENCE="00003"
        OUTPUT_NUMBER="00003"
        OUTPUT_TIME="10.0201"
        EXPECTED_SOURCE_SHARDS="512"
        EXPECTED_NODES="64"
        CONTROL_MODE="none"
        ;;
    *)
        fail "CANARY_KIND must be canary-256 or canary-no-control"
        ;;
esac

INPUT_DIR="${ARCHIVE_ROOT}/${SIMULATION}/${PHASE}/bin"
OUTPUT_DIR=${OUTPUT_DIR:-"${OUTPUT_ROOT}/validation/frontier_canary/${CANARY_KIND}/${SIMULATION}/${PHASE}/${SEQUENCE}/${SLURM_JOB_ID:-manual}"}
VALIDATOR_REPORT=${VALIDATOR_REPORT:-"${OUTPUT_DIR}/validator_report.json"}

require_run_mode
[[ "$PRODUCTS" == "all" ]] || fail "Frontier canaries require PRODUCTS=all"
require_positive_integer "$CHUNK_BLOCKS" "CHUNK_BLOCKS"
require_slurm_cluster frontier
if [[ "$RUN_MODE" == "execute" ]]; then
    [[ "$CONFIRM_FRONTIER_CANARY" == "$CANARY_KIND" ]] || {
        fail "execution requires CONFIRM_FRONTIER_CANARY=${CANARY_KIND}"
    }
    [[ "${SLURM_JOB_NUM_NODES:-0}" == "$EXPECTED_NODES" ]] || {
        fail "${CANARY_KIND} execution requires ${EXPECTED_NODES} nodes"
    }
    require_file "$RELEASE_IDENTITY" "frozen Frontier release identity"
fi

module restore
module load PrgEnv-cray craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2 cce/20.0.0 cmake
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
require_file "$BUDGET_FILE" "Frontier performance budgets"
require_file "$VALIDATOR" "full snapshot validator"
if [[ "$RUN_MODE" == "execute" ]]; then
    verify_release_identity "$RELEASE_IDENTITY" "$EXECUTABLE"
fi
require_source_shards "$INPUT_DIR" "$SEQUENCE" "$EXPECTED_SOURCE_SHARDS"
prepare_output_directory "$OUTPUT_DIR"

if [[ "$RUN_MODE" == "execute" ]]; then
    launch=(
        srun --nodes="$EXPECTED_NODES" --ntasks="$EXPECTED_SOURCE_SHARDS"
        --ntasks-per-node=8 --cpus-per-task=7 --gpus-per-node=8 --gpu-bind=closest
    )
else
    launch=(srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=7 --gpus-per-node=8)
fi
command=(
    "${launch[@]}"
    "$EXECUTABLE"
    --input-dir "$INPUT_DIR"
    --sequence "$SEQUENCE"
    --output-number "$OUTPUT_NUMBER"
    --output-time "$OUTPUT_TIME"
    --expected-shards "$EXPECTED_SOURCE_SHARDS"
    --output-dir "$OUTPUT_DIR"
    --products "$PRODUCTS"
    --chunk-blocks "$CHUNK_BLOCKS"
)
[[ "$RUN_MODE" == "plan" ]] && command+=(--dry-run)

log "Frontier ${CANARY_KIND}: RUN_MODE=${RUN_MODE}, source=${SIMULATION}/${PHASE}/${SEQUENCE}"
print_command "${command[@]}"
"${command[@]}"

if [[ "$RUN_MODE" == "execute" ]]; then
    validate_rebuild_manifest "$OUTPUT_DIR" "$EXPECTED_SOURCE_SHARDS" "$EXPECTED_SOURCE_SHARDS"
    validator_command=(
        "$PYTHON" "$VALIDATOR"
        "$OUTPUT_DIR"
        "${ARCHIVE_ROOT}/${SIMULATION}/${PHASE}"
        --expected-source-sequence "$SEQUENCE"
        --skip-source-cycle-check
        --expected-shards "$EXPECTED_SOURCE_SHARDS"
        --expected-domain-volume 64000000
        --expected-output-number "$OUTPUT_NUMBER"
        --expected-output-time "$OUTPUT_TIME"
        --archived-sequence "$OUTPUT_NUMBER"
        --report-path "$VALIDATOR_REPORT"
    )
    [[ "$CONTROL_MODE" == "none" ]] && validator_command+=(--no-controls)
    print_command "${validator_command[@]}"
    "${validator_command[@]}"
    pass_report=$("$PYTHON" "${JOBS_DIR}/release_gate.py" attest-run \
        --identity "$RELEASE_IDENTITY" \
        --repo-root "$REPO_ROOT" \
        --executable "$EXECUTABLE" \
        --output-dir "$OUTPUT_DIR" \
        --validator-report "$VALIDATOR_REPORT" \
        --budget-file "$BUDGET_FILE" \
        --budget-key "$CANARY_KIND" \
        --kind "$CANARY_KIND" \
        --nodes "$EXPECTED_NODES" \
        --artifact-dir "$RELEASE_ARTIFACT_DIR")
    log "Published immutable canary pass report: ${pass_report}"
fi
