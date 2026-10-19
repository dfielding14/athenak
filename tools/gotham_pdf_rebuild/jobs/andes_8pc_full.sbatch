#!/usr/bin/env bash
#SBATCH -A AST207
#SBATCH -J gotham_pdf_8pc_full
#SBATCH -o /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_8pc_full_%j.out
#SBATCH -e /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_8pc_full_%j.err
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -N 6
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=4

set -euo pipefail

# Exact archived-control snapshot; this is a correctness job, not a generic runner.

JOBS_DIR=${JOBS_DIR:-"/ccs/home/dfielding/athenak-gotham-pdf-rebuild/tools/gotham_pdf_rebuild/jobs"}
source "${JOBS_DIR}/common.sh"

EXECUTABLE=${EXECUTABLE:-"${REPO_ROOT}/build-gotham-pdf-andes-cuda/gotham_pdf_rebuild"}
SIMULATION="res_8pc"
PHASE="phase2"
SEQUENCE="00028"
OUTPUT_NUMBER="00062"
OUTPUT_TIME="10.520002963096713"
INPUT_DIR="${ARCHIVE_ROOT}/${SIMULATION}/${PHASE}/bin"
EXPECTED_SOURCE_SHARDS="1024"
PRODUCTS=${PRODUCTS:-"original"}
CHUNK_BLOCKS=${CHUNK_BLOCKS:-"32"}
CONFIRM_FULL_8PC=${CONFIRM_FULL_8PC:-""}
OUTPUT_DIR=${OUTPUT_DIR:-"${OUTPUT_ROOT}/validation/andes_8pc_full/${SIMULATION}/${PHASE}/${SEQUENCE}/${SLURM_JOB_ID:-manual}"}
GPU_NODES=${GPU_NODES:-"${SLURM_JOB_NUM_NODES:-6}"}

require_run_mode
require_products "$PRODUCTS"
require_positive_integer "$CHUNK_BLOCKS" "CHUNK_BLOCKS"
require_positive_integer "$GPU_NODES" "GPU_NODES"
require_slurm_cluster andes
if [[ "$RUN_MODE" == "execute" ]]; then
    expected_confirmation="${SIMULATION}/${PHASE}/${SEQUENCE}"
    [[ "$CONFIRM_FULL_8PC" == "$expected_confirmation" ]] || {
        fail "full 8 pc execution requires CONFIRM_FULL_8PC=${expected_confirmation}"
    }
fi
[[ "${SLURM_JOB_NUM_NODES:-$GPU_NODES}" == "$GPU_NODES" ]] || {
    fail "GPU_NODES=${GPU_NODES} does not match the ${SLURM_JOB_NUM_NODES} allocated nodes"
}
NTASKS=$((GPU_NODES * 4))

module reset
module load cuda/11.2.2

require_executable "$EXECUTABLE" "Andes CUDA reducer"
require_source_shards "$INPUT_DIR" "$SEQUENCE" "$EXPECTED_SOURCE_SHARDS"
prepare_output_directory "$OUTPUT_DIR"

command=(
    srun --nodes="$GPU_NODES" --ntasks="$NTASKS" --ntasks-per-node=4 --cpus-per-task=8
    --gpus-per-task=1 --gpu-bind=closest
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

log "Andes full 8 pc validation: RUN_MODE=${RUN_MODE}, source=${SIMULATION}/${PHASE}/${SEQUENCE}"
print_command "${command[@]}"
"${command[@]}"

if [[ "$RUN_MODE" == "execute" ]]; then
    validate_rebuild_manifest "$OUTPUT_DIR" "$EXPECTED_SOURCE_SHARDS" "$EXPECTED_SOURCE_SHARDS"
fi
