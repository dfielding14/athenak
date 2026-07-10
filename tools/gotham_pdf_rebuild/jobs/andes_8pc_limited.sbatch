#!/usr/bin/env bash
#SBATCH -A AST207
#SBATCH -J gotham_pdf_8pc_limited
#SBATCH -o /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_8pc_limited_%j.out
#SBATCH -e /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_8pc_limited_%j.err
#SBATCH -t 04:00:00
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=4

set -euo pipefail

JOBS_DIR=${JOBS_DIR:-"/ccs/home/dfielding/athenak-gotham-pdf-rebuild/tools/gotham_pdf_rebuild/jobs"}
source "${JOBS_DIR}/common.sh"

EXECUTABLE=${EXECUTABLE:-"${REPO_ROOT}/build-gotham-pdf-andes-cuda/gotham_pdf_rebuild"}
SIMULATION=${SIMULATION:-"res_8pc_highmdot"}
PHASE=${PHASE:-"phase2"}
SEQUENCE=${SEQUENCE:-"00020"}
INPUT_DIR=${INPUT_DIR:-"${ARCHIVE_ROOT}/${SIMULATION}/${PHASE}/bin"}
EXPECTED_SOURCE_SHARDS=${EXPECTED_SOURCE_SHARDS:-"1024"}
PRODUCTS=${PRODUCTS:-"all"}
CHUNK_BLOCKS=${CHUNK_BLOCKS:-"32"}
MAX_SHARDS=${MAX_SHARDS:-"16"}
MAX_BLOCKS_PER_SHARD=${MAX_BLOCKS_PER_SHARD:-"0"}
OUTPUT_DIR=${OUTPUT_DIR:-"${OUTPUT_ROOT}/validation/andes_8pc_limited/${SIMULATION}/${PHASE}/${SEQUENCE}/${SLURM_JOB_ID:-manual}"}

require_run_mode
require_products "$PRODUCTS"
require_positive_integer "$CHUNK_BLOCKS" "CHUNK_BLOCKS"
require_positive_integer "$MAX_SHARDS" "MAX_SHARDS"
require_unsigned_integer "$MAX_BLOCKS_PER_SHARD" "MAX_BLOCKS_PER_SHARD"
require_slurm_cluster andes

module reset
module load cuda/11.2.2

require_executable "$EXECUTABLE" "Andes CUDA reducer"
require_source_shards "$INPUT_DIR" "$SEQUENCE" "$EXPECTED_SOURCE_SHARDS"
prepare_output_directory "$OUTPUT_DIR"

command=(
    srun --nodes=1 --ntasks=4 --ntasks-per-node=4 --cpus-per-task=8
    --gpus-per-task=1 --gpu-bind=closest
    "$EXECUTABLE"
    --input-dir "$INPUT_DIR"
    --sequence "$SEQUENCE"
    --output-dir "$OUTPUT_DIR"
    --products "$PRODUCTS"
    --chunk-blocks "$CHUNK_BLOCKS"
    --max-shards "$MAX_SHARDS"
)
if (( MAX_BLOCKS_PER_SHARD > 0 )); then
    command+=(--max-blocks-per-shard "$MAX_BLOCKS_PER_SHARD")
fi
[[ "$RUN_MODE" == "plan" ]] && command+=(--dry-run)

log "Andes limited 8 pc validation: RUN_MODE=${RUN_MODE}, source=${SIMULATION}/${PHASE}/${SEQUENCE}"
print_command "${command[@]}"
"${command[@]}"

if [[ "$RUN_MODE" == "execute" ]]; then
    validate_rebuild_manifest "$OUTPUT_DIR" "$MAX_SHARDS" "$MAX_SHARDS"
fi
