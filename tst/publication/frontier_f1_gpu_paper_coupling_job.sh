#!/bin/bash
# Run the bounded one-rank Frontier GPU paper-mode coupling oracle.
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:15:00
#SBATCH --job-name=pic-f1-gpu-paper-coupling
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.out

set -euo pipefail

PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
export PATH=/usr/bin:/bin
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
CAMPAIGN=f1_gpu_paper_coupling
SUBMISSION_DIR="${PIC_ROOT}/manifests/${CAMPAIGN}/${PIC_SUBMISSION_ID}"
SNAPSHOT_DIR="${SUBMISSION_DIR}/snapshot"
ARTIFACT_DIR="${PIC_ROOT}/runs/${CAMPAIGN}/${PIC_SUBMISSION_ID}"
OUTPUT_ROOT="${ARTIFACT_DIR}/output"

mkdir -p "${OUTPUT_ROOT}/coeff0" "${OUTPUT_ROOT}/coeff7"
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
/usr/bin/ldd "${SNAPSHOT_DIR}/athena" > "${ARTIFACT_DIR}/athena_ldd.txt"

srun -N1 -n1 -c1 --gpus-per-task=1 --gpu-bind=closest \
  "$ROCM_SMI" --showproductname > "${ARTIFACT_DIR}/rocm_smi.txt"
srun -N1 -n1 -c7 --gpus-per-task=1 --gpu-bind=closest \
  /bin/bash -c \
  'printf "rank=%s host=%s ROCR_VISIBLE_DEVICES=%s GPU_DEVICE_ORDINAL=%s HIP_VISIBLE_DEVICES=%s\n" \
    "${SLURM_PROCID}" "$(hostname)" "${ROCR_VISIBLE_DEVICES:-unset}" \
    "${GPU_DEVICE_ORDINAL:-unset}" "${HIP_VISIBLE_DEVICES:-unset}"' \
  > "${ARTIFACT_DIR}/gpu_mapping.txt"

srun -N1 -n1 -c7 --gpus-per-task=1 --gpu-bind=closest \
  "${SNAPSHOT_DIR}/athena" \
  -i "${SNAPSHOT_DIR}/pic_paper_coupling_conservation.athinput" \
  -d "${OUTPUT_ROOT}/coeff0" \
  job/basename=f1_gpu_paper_coupling_coeff0 \
  time/nlim=2 \
  particles/couple_j_to_efield_coeff=0.0 \
  > "${ARTIFACT_DIR}/athena_coeff0_stdout.txt" \
  2> "${ARTIFACT_DIR}/athena_coeff0_stderr.txt"

srun -N1 -n1 -c7 --gpus-per-task=1 --gpu-bind=closest \
  "${SNAPSHOT_DIR}/athena" \
  -i "${SNAPSHOT_DIR}/pic_paper_coupling_conservation.athinput" \
  -d "${OUTPUT_ROOT}/coeff7" \
  job/basename=f1_gpu_paper_coupling_coeff7 \
  time/nlim=2 \
  particles/couple_j_to_efield_coeff=7.0 \
  > "${ARTIFACT_DIR}/athena_coeff7_stdout.txt" \
  2> "${ARTIFACT_DIR}/athena_coeff7_stderr.txt"

"$PYTHON" "${SNAPSHOT_DIR}/analysis/000-frontier_f1_gpu_paper_coupling_analysis.py" \
  --artifact-dir "$ARTIFACT_DIR"
(
  cd "$ARTIFACT_DIR"
  find . -type f ! -name SHA256SUMS -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > SHA256SUMS
)
