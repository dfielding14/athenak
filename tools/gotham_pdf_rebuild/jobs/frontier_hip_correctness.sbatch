#!/usr/bin/env bash
#SBATCH -A AST207
#SBATCH -J gotham_pdf_hip_check
#SBATCH -o /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_hip_check_%j.out
#SBATCH -e /lustre/orion/ast207/proj-shared/gotham/logs/gotham_pdf_hip_check_%j.err
#SBATCH -t 00:30:00
#SBATCH -p batch
#SBATCH -q debug
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7

# Pre-freeze Frontier HIP and MPI correctness check on limited real input.

set -euo pipefail

JOBS_DIR=${JOBS_DIR:-"/ccs/home/dfielding/athenak-gotham-pdf-rebuild/tools/gotham_pdf_rebuild/jobs"}
source "${JOBS_DIR}/common.sh"

EXECUTABLE=${EXECUTABLE:-"${REPO_ROOT}/build-gotham-pdf-frontier/gotham_pdf_rebuild"}
COMPARATOR=${COMPARATOR:-"${JOBS_DIR}/compare_rebuild_outputs.py"}
INPUT_DIR=${INPUT_DIR:-"${ARCHIVE_ROOT}/res_8pc/phase2/bin"}
SEQUENCE=${SEQUENCE:-"00028"}
EXPECTED_SOURCE_SHARDS=${EXPECTED_SOURCE_SHARDS:-"1024"}
PRODUCTS=${PRODUCTS:-"all"}
CHUNK_BLOCKS=${CHUNK_BLOCKS:-"2"}
MAX_SHARDS=${MAX_SHARDS:-"8"}
MAX_BLOCKS_PER_SHARD=${MAX_BLOCKS_PER_SHARD:-"2"}
OUTPUT_BASE=${OUTPUT_BASE:-"${OUTPUT_ROOT}/validation/frontier_hip_correctness/${SLURM_JOB_ID:-manual}"}
SINGLE_OUTPUT="${OUTPUT_BASE}/single_rank"
EIGHT_OUTPUT="${OUTPUT_BASE}/eight_rank"

require_run_mode
require_products "$PRODUCTS"
require_positive_integer "$CHUNK_BLOCKS" "CHUNK_BLOCKS"
require_positive_integer "$MAX_SHARDS" "MAX_SHARDS"
require_positive_integer "$MAX_BLOCKS_PER_SHARD" "MAX_BLOCKS_PER_SHARD"
require_slurm_cluster frontier

module restore
module load PrgEnv-cray craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2 cce/20.0.0 cmake
module unload darshan-runtime

export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH:-}:${LD_LIBRARY_PATH:-}"
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export HSA_XNACK=1

require_executable "$EXECUTABLE" "Frontier HIP reducer"
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
        srun --nodes=1 --ntasks=8 --ntasks-per-node=8 --cpus-per-task=7
        --gpus-per-node=8 --gpu-bind=closest
        "$EXECUTABLE" "${base_arguments[@]}" --output-dir "$EIGHT_OUTPUT" --dry-run
    )
    print_command "${command[@]}"
    "${command[@]}"
    exit 0
fi

prepare_output_directory "$SINGLE_OUTPUT"
prepare_output_directory "$EIGHT_OUTPUT"

single_command=(
    srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=7
    --gpus-per-node=8 --gpu-bind=closest
    "$EXECUTABLE" "${base_arguments[@]}" --output-dir "$SINGLE_OUTPUT"
)
eight_command=(
    srun --nodes=1 --ntasks=8 --ntasks-per-node=8 --cpus-per-task=7
    --gpus-per-node=8 --gpu-bind=closest
    "$EXECUTABLE" "${base_arguments[@]}" --output-dir "$EIGHT_OUTPUT"
)

print_command "${single_command[@]}"
"${single_command[@]}"
validate_rebuild_manifest "$SINGLE_OUTPUT" "$MAX_SHARDS" "$MAX_SHARDS"

print_command "${eight_command[@]}"
"${eight_command[@]}"
validate_rebuild_manifest "$EIGHT_OUTPUT" "$MAX_SHARDS" "$MAX_SHARDS"

"$PYTHON" "$COMPARATOR" "$SINGLE_OUTPUT" "$EIGHT_OUTPUT" \
    --rtol 1e-9 --global-rtol 1e-12 --atol 1e-12
