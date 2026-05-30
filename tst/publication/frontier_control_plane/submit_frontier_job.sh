#!/bin/bash
# Submit only immutable, validated Frontier PIC manifests through this wrapper.
set -euo pipefail

PIC_ROOT="${PIC_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/PIC}"
PROJECT_HOME_MIRROR_ROOT="${PROJECT_HOME_MIRROR_ROOT:-/ccs/proj/ast207/proj-shared/PIC}"
MANIFEST="${1:?usage: submit_frontier_job.sh PRE_SUBMIT_MANIFEST}"
CONTROL_PLANE_DIR="$(cd "$(/usr/bin/dirname "$0")" && pwd)"
VALIDATOR="${CONTROL_PLANE_DIR}/validate_and_reserve_frontier_job.py"
TRAMPOLINE="${CONTROL_PLANE_DIR}/launch_trampoline.py"
PYTHON="/opt/cray/pe/python/3.11.7/bin/python3"
SBATCH="/usr/bin/sbatch"
SCANCEL="/usr/bin/scancel"
LEDGER_JSONL="${PIC_ROOT}/ledger/node_hours.jsonl"
LEDGER_CSV="${PIC_ROOT}/ledger/node_hours.csv"
RECEIPTS_JSONL="${PIC_ROOT}/ledger/mirror_receipts.jsonl"
MIRROR_JSONL="${PROJECT_HOME_MIRROR_ROOT}/ledger/node_hours.jsonl"
PENDING_FILE="${PIC_ROOT}/ledger/pending_submission.json"

[[ "$PIC_ROOT" == "/lustre/orion/ast207/proj-shared/dfielding/PIC" ]] || {
  printf "Unauthorized PIC_ROOT: %s\n" "$PIC_ROOT" >&2
  exit 1
}
[[ "$PROJECT_HOME_MIRROR_ROOT" == "/ccs/proj/ast207/proj-shared/PIC" ]] || {
  printf "Unauthorized PROJECT_HOME_MIRROR_ROOT: %s\n" "$PROJECT_HOME_MIRROR_ROOT" >&2
  exit 1
}
[[ "$(/usr/bin/dirname "$CONTROL_PLANE_DIR")" == "${PIC_ROOT}/control_plane" ]] || {
  printf "Run the submission wrapper from an installed control-plane version.\n" >&2
  exit 1
}
[[ ! -L "$0" && "$CONTROL_PLANE_DIR" == "$(/usr/bin/readlink -f -- "$CONTROL_PLANE_DIR")" ]] || {
  printf "Refusing a symlink alias for the installed control-plane version.\n" >&2
  exit 1
}
"$PYTHON" "$VALIDATOR" verify-control-plane >/dev/null

reservation_id=""
job_id=""
attached=0

cancel_unsubmitted_reservation() {
  status="$?"
  set +e
  if [[ -n "$job_id" && "$attached" -eq 0 ]]; then
    "$SCANCEL" "$job_id"
    printf "Job %s requires manual attachment and reconciliation; see %s\n" \
      "$job_id" "$PENDING_FILE" >&2
  elif [[ -n "$reservation_id" ]]; then
    "$PYTHON" "$VALIDATOR" cancel-reservation \
      --ledger-jsonl "$LEDGER_JSONL" \
      --ledger-csv "$LEDGER_CSV" \
      --receipts-jsonl "$RECEIPTS_JSONL" \
      --mirror-jsonl "$MIRROR_JSONL" \
      --reservation-id "$reservation_id" \
      --notes "sbatch failed before scheduler job ID assignment"
  fi
  exit "$status"
}
trap cancel_unsubmitted_reservation ERR INT TERM

test -r "$LEDGER_JSONL" || {
  printf "Missing initialized mirrored ledger: %s\n" "$LEDGER_JSONL" >&2
  exit 1
}

reservation_id="$(
  "$PYTHON" "$VALIDATOR" reserve \
    --manifest "$MANIFEST" \
    --ledger-jsonl "$LEDGER_JSONL" \
    --ledger-csv "$LEDGER_CSV" \
    --receipts-jsonl "$RECEIPTS_JSONL" \
    --mirror-jsonl "$MIRROR_JSONL" \
    --node-hour-cap 10000 \
    --pending-marker "$PENDING_FILE"
)"
submission_id="$("$PYTHON" "$VALIDATOR" submission-id \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
manifest_sha256="$("$PYTHON" "$VALIDATOR" manifest-sha256 \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
job_script_sha256="$("$PYTHON" "$VALIDATOR" snapshot-sha256 \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --role job-script \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
executable_sha256="$("$PYTHON" "$VALIDATOR" snapshot-sha256 \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --role executable \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
account="$("$PYTHON" "$VALIDATOR" directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key account \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
partition="$("$PYTHON" "$VALIDATOR" directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key partition \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
qos="$("$PYTHON" "$VALIDATOR" directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key qos \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
nodes="$("$PYTHON" "$VALIDATOR" directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key nodes \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
walltime="$("$PYTHON" "$VALIDATOR" directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key time \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
output="$("$PYTHON" "$VALIDATOR" directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key output \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"

job_id="$("$SBATCH" --parsable \
  --comment "pic-reservation=${reservation_id}" \
  --account "$account" \
  --partition "$partition" \
  --qos "$qos" \
  --nodes "$nodes" \
  --time "$walltime" \
  --output "$output" \
  --export "PIC_RESERVATION_ID=${reservation_id},PIC_SUBMISSION_ID=${submission_id},PIC_MANIFEST_SHA256=${manifest_sha256}" \
  "$TRAMPOLINE" \
  --manifest "$MANIFEST" \
  --manifest-sha256 "$manifest_sha256" \
  --job-script-sha256 "$job_script_sha256" \
  --executable-sha256 "$executable_sha256" \
  --reservation-id "$reservation_id" \
  --submission-id "$submission_id" \
  --ledger-jsonl "$LEDGER_JSONL" \
  --receipts-jsonl "$RECEIPTS_JSONL" \
  --mirror-jsonl "$MIRROR_JSONL")"
job_id="${job_id%%;*}"

"$PYTHON" "$VALIDATOR" mark-submitted \
  --ledger-jsonl "$LEDGER_JSONL" \
  --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" \
  --mirror-jsonl "$MIRROR_JSONL" \
  --reservation-id "$reservation_id" \
  --job-id "$job_id"

"$PYTHON" "$VALIDATOR" attach-job-id \
  --ledger-jsonl "$LEDGER_JSONL" \
  --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" \
  --mirror-jsonl "$MIRROR_JSONL" \
  --reservation-id "$reservation_id" \
  --job-id "$job_id"
attached=1
trap - ERR INT TERM
printf "Submitted %s with reservation %s\n" "$job_id" "$reservation_id"
