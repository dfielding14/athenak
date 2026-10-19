#!/usr/bin/env bash
#SBATCH -A AST207
#SBATCH -J gotham_pdf_k80_check
#SBATCH -o /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_k80_check_%j.out
#SBATCH -e /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_k80_check_%j.err
#SBATCH -t 02:00:00
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=4

set -euo pipefail

JOBS_DIR=${JOBS_DIR:-"/ccs/home/dfielding/athenak-gotham-pdf-rebuild/tools/gotham_pdf_rebuild/jobs"}
source "${JOBS_DIR}/common.sh"

EXECUTABLE=${EXECUTABLE:-"${REPO_ROOT}/build-gotham-pdf-andes-cuda/gotham_pdf_rebuild"}
COMPARATOR=${COMPARATOR:-"${JOBS_DIR}/compare_rebuild_outputs.py"}
INPUT_DIR=${INPUT_DIR:-"${ARCHIVE_ROOT}/res_8pc_highmdot/phase2/bin"}
SEQUENCE=${SEQUENCE:-"00020"}
EXPECTED_SOURCE_SHARDS=${EXPECTED_SOURCE_SHARDS:-"1024"}
PRODUCTS=${PRODUCTS:-"all"}
CHUNK_BLOCKS=${CHUNK_BLOCKS:-"8"}
MAX_SHARDS=${MAX_SHARDS:-"4"}
MAX_BLOCKS_PER_SHARD=${MAX_BLOCKS_PER_SHARD:-"8"}
RTOL=${RTOL:-"1e-9"}
GLOBAL_RTOL=${GLOBAL_RTOL:-"1e-12"}
ATOL=${ATOL:-"1e-12"}
OUTPUT_BASE=${OUTPUT_BASE:-"${OUTPUT_ROOT}/validation/andes_k80_correctness/${SLURM_JOB_ID:-manual}"}
SINGLE_OUTPUT="${OUTPUT_BASE}/single_rank"
FOUR_OUTPUT="${OUTPUT_BASE}/four_rank"

require_run_mode
require_products "$PRODUCTS"
require_positive_integer "$CHUNK_BLOCKS" "CHUNK_BLOCKS"
require_positive_integer "$MAX_SHARDS" "MAX_SHARDS"
require_positive_integer "$MAX_BLOCKS_PER_SHARD" "MAX_BLOCKS_PER_SHARD"
require_slurm_cluster andes

module reset
module load cuda/11.2.2

require_executable "$EXECUTABLE" "Andes CUDA reducer"
require_file "$COMPARATOR" "streaming output comparator"
require_source_shards "$INPUT_DIR" "$SEQUENCE" "$EXPECTED_SOURCE_SHARDS"

base_arguments=(
    --input-dir "$INPUT_DIR"
    --sequence "$SEQUENCE"
    --products "$PRODUCTS"
    --chunk-blocks "$CHUNK_BLOCKS"
    --max-shards "$MAX_SHARDS"
    --max-blocks-per-shard "$MAX_BLOCKS_PER_SHARD"
)

if [[ "$RUN_MODE" == "plan" ]]; then
    command=(
        srun --nodes=1 --ntasks=4 --ntasks-per-node=4 --cpus-per-task=8
        --gpus-per-task=1 --gpu-bind=closest
        "$EXECUTABLE" "${base_arguments[@]}" --output-dir "$FOUR_OUTPUT" --dry-run
    )
    log "Andes K80 correctness plan"
    print_command "${command[@]}"
    "${command[@]}"
    exit 0
fi

prepare_output_directory "$SINGLE_OUTPUT"
prepare_output_directory "$FOUR_OUTPUT"

single_command=(
    srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=8
    --gpus-per-task=1 --gpu-bind=closest
    "$EXECUTABLE" "${base_arguments[@]}" --output-dir "$SINGLE_OUTPUT"
)
four_command=(
    srun --nodes=1 --ntasks=4 --ntasks-per-node=4 --cpus-per-task=8
    --gpus-per-task=1 --gpu-bind=closest
    "$EXECUTABLE" "${base_arguments[@]}" --output-dir "$FOUR_OUTPUT"
)

log "Running one-rank K80 reference"
print_command "${single_command[@]}"
"${single_command[@]}"
validate_rebuild_manifest "$SINGLE_OUTPUT" "$MAX_SHARDS" "$MAX_SHARDS"

log "Running four-rank K80 candidate"
print_command "${four_command[@]}"
"${four_command[@]}"
validate_rebuild_manifest "$FOUR_OUTPUT" "$MAX_SHARDS" "$MAX_SHARDS"

log "Comparing one-rank and four-rank K80 outputs"
python3 "$COMPARATOR" "$SINGLE_OUTPUT" "$FOUR_OUTPUT" \
    --rtol "$RTOL" --global-rtol "$GLOBAL_RTOL" --atol "$ATOL"
