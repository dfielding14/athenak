#!/bin/bash
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:15:00
#SBATCH --job-name=pic-f0-hipmpi-smoke
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.out

set -euo pipefail

PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
export PATH=/usr/bin:/bin
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
SUBMISSION_DIR="${PIC_ROOT}/manifests/f0_hipmpi_smoke/${PIC_SUBMISSION_ID}"
SNAPSHOT_DIR="${SUBMISSION_DIR}/snapshot"
ARTIFACT_DIR="${PIC_ROOT}/runs/f0_hipmpi_smoke/${PIC_SUBMISSION_ID}"

mkdir -p "$ARTIFACT_DIR"
exec > >(tee "${ARTIFACT_DIR}/job.log") 2>&1

"$PYTHON" "${SNAPSHOT_DIR}/verify_compute_node_snapshot.py" \
  --manifest "${SUBMISSION_DIR}/pre_submit_manifest.json" \
  --submission-id "$PIC_SUBMISSION_ID" \
  --reservation-id "$PIC_RESERVATION_ID" \
  --manifest-sha256 "$PIC_MANIFEST_SHA256"
printf 'immutable snapshot verified\n' > "${ARTIFACT_DIR}/snapshot_verification.txt"

if ! type module >/dev/null 2>&1; then
  source /etc/profile
fi
source "${SNAPSHOT_DIR}/frontier_pic_environment.sh"
export SLURM_EXPORT_ENV=ALL
ROCM_SMI="${ROCM_PATH}/bin/rocm-smi"
test -x "$ROCM_SMI"

record_pic_environment > "${ARTIFACT_DIR}/environment.allowlist.txt"
module -t list 2> "${ARTIFACT_DIR}/modules.txt"

srun -N1 -n1 -c1 --gpus-per-task=1 --gpu-bind=closest \
  "$ROCM_SMI" --showproductname > "${ARTIFACT_DIR}/rocm_smi.txt"
srun -N1 -n8 -c7 --gpus-per-task=1 --gpu-bind=closest \
  /bin/bash -c \
  'printf "rank=%s host=%s ROCR_VISIBLE_DEVICES=%s GPU_DEVICE_ORDINAL=%s HIP_VISIBLE_DEVICES=%s\n" \
    "${SLURM_PROCID}" "$(hostname)" "${ROCR_VISIBLE_DEVICES:-unset}" \
    "${GPU_DEVICE_ORDINAL:-unset}" "${HIP_VISIBLE_DEVICES:-unset}"' \
  | sort -n > "${ARTIFACT_DIR}/gpu_mapping.txt"

srun -N1 -n1 -c7 --gpus-per-task=1 --gpu-bind=closest \
  "${SNAPSHOT_DIR}/athena" \
  -i "${SNAPSHOT_DIR}/pic_parser_contract_guards.athinput" \
  -n \
  > "${ARTIFACT_DIR}/athena_parser_stdout.txt" \
  2> "${ARTIFACT_DIR}/athena_parser_stderr.txt"

"$PYTHON" "${SNAPSHOT_DIR}/analysis/000-frontier_f0_smoke_analysis.py" \
  --artifact-dir "$ARTIFACT_DIR"
(
  cd "$ARTIFACT_DIR"
  sha256sum \
    snapshot_verification.txt \
    environment.allowlist.txt \
    modules.txt \
    rocm_smi.txt \
    gpu_mapping.txt \
    athena_parser_stdout.txt \
    athena_parser_stderr.txt \
    analysis.json \
    > SHA256SUMS
)
