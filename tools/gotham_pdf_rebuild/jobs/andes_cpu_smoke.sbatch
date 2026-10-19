#!/usr/bin/env bash
#SBATCH -A AST207
#SBATCH -J gotham_pdf_cpu_smoke
#SBATCH -o /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_cpu_smoke_%j.out
#SBATCH -e /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_cpu_smoke_%j.err
#SBATCH -t 00:30:00
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

JOBS_DIR=${JOBS_DIR:-"/ccs/home/dfielding/athenak-gotham-pdf-rebuild/tools/gotham_pdf_rebuild/jobs"}
source "${JOBS_DIR}/common.sh"

EXECUTABLE=${EXECUTABLE:-"${REPO_ROOT}/build-gotham-pdf-cpu/gotham_pdf_rebuild"}
INPUT_DIR=${INPUT_DIR:-"${ARCHIVE_ROOT}/res_8pc_highmdot/phase2/bin"}
SEQUENCE=${SEQUENCE:-"00020"}
EXPECTED_SOURCE_SHARDS=${EXPECTED_SOURCE_SHARDS:-"1024"}
PRODUCTS=${PRODUCTS:-"all"}
CHUNK_BLOCKS=${CHUNK_BLOCKS:-"1"}
MAX_SHARDS=${MAX_SHARDS:-"1"}
MAX_BLOCKS_PER_SHARD=${MAX_BLOCKS_PER_SHARD:-"1"}
OUTPUT_DIR=${OUTPUT_DIR:-"${OUTPUT_ROOT}/validation/andes_cpu_smoke/${SLURM_JOB_ID:-manual}"}

require_run_mode
require_products "$PRODUCTS"
require_positive_integer "$CHUNK_BLOCKS" "CHUNK_BLOCKS"
require_positive_integer "$MAX_SHARDS" "MAX_SHARDS"
require_positive_integer "$MAX_BLOCKS_PER_SHARD" "MAX_BLOCKS_PER_SHARD"
require_slurm_cluster andes

module reset
module load gcc/9.3.0 openmpi/4.0.4

require_executable "$EXECUTABLE" "CPU reducer"
require_source_shards "$INPUT_DIR" "$SEQUENCE" "$EXPECTED_SOURCE_SHARDS"
prepare_output_directory "$OUTPUT_DIR"

command=(
    srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=2
    "$EXECUTABLE"
    --input-dir "$INPUT_DIR"
    --sequence "$SEQUENCE"
    --output-dir "$OUTPUT_DIR"
    --products "$PRODUCTS"
    --chunk-blocks "$CHUNK_BLOCKS"
    --max-shards "$MAX_SHARDS"
    --max-blocks-per-shard "$MAX_BLOCKS_PER_SHARD"
)
[[ "$RUN_MODE" == "plan" ]] && command+=(--dry-run)

log "Andes CPU smoke: RUN_MODE=${RUN_MODE}"
print_command "${command[@]}"
"${command[@]}"

if [[ "$RUN_MODE" == "execute" ]]; then
    validate_rebuild_manifest "$OUTPUT_DIR" "$MAX_SHARDS" "$MAX_SHARDS"
fi
