#!/bin/bash
# Submit only immutable, validated Frontier PIC manifests through this wrapper.
set -euo pipefail

PIC_ROOT="${PIC_ROOT:-/lustre/orion/ast207/proj-shared/dfielding/PIC}"
PROJECT_HOME_MIRROR_ROOT="${PROJECT_HOME_MIRROR_ROOT:-/ccs/proj/ast207/proj-shared/PIC}"
MANIFEST="${1:?usage: submit_frontier_job.sh PRE_SUBMIT_MANIFEST}"
CONTROL_PLANE_DIR="$(cd "$(/usr/bin/dirname "$0")" && pwd)"
VALIDATOR="${CONTROL_PLANE_DIR}/validate_and_reserve_frontier_job.py"
TRAMPOLINE="${CONTROL_PLANE_DIR}/launch_trampoline.py"
RUNNER="${CONTROL_PLANE_DIR}/run_control_plane.py"
BATCH_DELEGATE="${CONTROL_PLANE_DIR}/run_installed_control_plane_job.sh"
PYTHON=(/opt/cray/pe/python/3.11.7/bin/python3 -I)
CONTROL_PLANE=("${PYTHON[@]}" "$RUNNER")
SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin SLURM_CLUSTERS=frontier)
SBATCH="/usr/bin/sbatch"
SCANCEL="/usr/bin/scancel"
SCONTROL="/usr/bin/scontrol"
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
"${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py verify-control-plane >/dev/null

reservation_id=""
job_id=""
attached=0
dispatch_started=0

cancel_unsubmitted_reservation() {
  status="$?"
  set +e
  if [[ -n "$job_id" ]]; then
    "${SLURM_ENV[@]}" "$SCANCEL" "$job_id"
    printf "Job %s requires reconciliation after submission failure; see %s\n" \
      "$job_id" "$PENDING_FILE" >&2
  elif [[ -n "$reservation_id" && "$dispatch_started" -eq 0 ]]; then
    if ! "${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py cancel-reservation \
        --ledger-jsonl "$LEDGER_JSONL" \
        --ledger-csv "$LEDGER_CSV" \
        --receipts-jsonl "$RECEIPTS_JSONL" \
        --mirror-jsonl "$MIRROR_JSONL" \
        --reservation-id "$reservation_id" \
        --notes "submission wrapper failed before scheduler dispatch"; then
      printf "Reservation %s may have entered scheduler dispatch; inspect Slurm and %s before reconciliation.\n" \
        "$reservation_id" "$PENDING_FILE" >&2
    fi
  elif [[ -n "$reservation_id" ]]; then
    printf "Reservation %s may have entered scheduler dispatch; inspect Slurm and %s before reconciliation.\n" \
      "$reservation_id" "$PENDING_FILE" >&2
  fi
  exit "$status"
}
trap cancel_unsubmitted_reservation ERR INT TERM

test -r "$LEDGER_JSONL" || {
  printf "Missing initialized mirrored ledger: %s\n" "$LEDGER_JSONL" >&2
  exit 1
}

reservation_id="$(
  "${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py reserve \
    --manifest "$MANIFEST" \
    --ledger-jsonl "$LEDGER_JSONL" \
    --ledger-csv "$LEDGER_CSV" \
    --receipts-jsonl "$RECEIPTS_JSONL" \
    --mirror-jsonl "$MIRROR_JSONL" \
    --node-hour-cap 10000 \
    --pending-marker "$PENDING_FILE"
)"
submission_id="$("${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py submission-id \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
manifest_sha256="$("${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py manifest-sha256 \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
job_script_sha256="$("${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py snapshot-sha256 \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --role job-script \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
executable_sha256="$("${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py snapshot-sha256 \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --role executable \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
account="$("${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key account \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
partition="$("${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key partition \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
qos="$("${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key qos \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
nodes="$("${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key nodes \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
walltime="$("${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key time \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"
output="$("${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py directive \
  --manifest "$MANIFEST" --reservation-id "$reservation_id" --key output \
  --ledger-jsonl "$LEDGER_JSONL" --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" --mirror-jsonl "$MIRROR_JSONL")"

"${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py mark-dispatch-started \
  --ledger-jsonl "$LEDGER_JSONL" \
  --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" \
  --mirror-jsonl "$MIRROR_JSONL" \
  --reservation-id "$reservation_id"
dispatch_started=1

job_id="$("${SLURM_ENV[@]}" "$SBATCH" --parsable --hold \
  --comment "pic-reservation=${reservation_id}" \
  --account "$account" \
  --partition "$partition" \
  --qos "$qos" \
  --nodes "$nodes" \
  --time "$walltime" \
  --output "$output" \
  --export=NIL \
  "$BATCH_DELEGATE" "$RUNNER" launch_trampoline.py \
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

"${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py mark-submitted \
  --ledger-jsonl "$LEDGER_JSONL" \
  --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" \
  --mirror-jsonl "$MIRROR_JSONL" \
  --reservation-id "$reservation_id" \
  --job-id "$job_id"

"${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py attach-job-id \
  --ledger-jsonl "$LEDGER_JSONL" \
  --ledger-csv "$LEDGER_CSV" \
  --receipts-jsonl "$RECEIPTS_JSONL" \
  --mirror-jsonl "$MIRROR_JSONL" \
  --reservation-id "$reservation_id" \
  --job-id "$job_id"
attached=1
"${SLURM_ENV[@]}" "$SCONTROL" release "$job_id"
trap - ERR INT TERM
printf "Submitted %s with reservation %s\n" "$job_id" "$reservation_id"
