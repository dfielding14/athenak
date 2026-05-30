#!/usr/bin/env bash
#
# Packet-local Slurm launcher for the external IO qualification rows.
#
# Required environment:
#   PACKET REPO BUILD PYTHON
#
# Optional environment:
#   ROW_TIMEOUT          Per-launch timeout in seconds. Default: 300.
#   ATTEMPT_ID           Unique packet-local retry or topology token.
#                       Default: attempt-1. Reusing one launch key is rejected.
#   FS_ACCOUNTING_HOOK   Site-provided executable invoked as:
#                         hook before ROW SAMPLE OUTPUT_DIR
#                         hook after  ROW SAMPLE OUTPUT_DIR
#                       Its stdout/stderr are retained. An unset hook leaves the
#                       observed-filesystem-read row explicitly incomplete.
#   MR2_HOSTFILE         Three-line Slurm arbitrary-distribution hostfile for
#                       MR-2: host A once, then host B twice.
#
# Usage:
#   run_external_io_qualification_slurm.sh ED-1 generate
#   run_external_io_qualification_slurm.sh MR-1 generate
#   run_external_io_qualification_slurm.sh MR-1 resume MANIFEST
#   run_external_io_qualification_slurm.sh MR-2 generate
#   run_external_io_qualification_slurm.sh MR-2 resume MANIFEST
#   run_external_io_qualification_slurm.sh scaling resume MANIFEST
#
# For "scaling", set NODES, RANKS, and RANKS_PER_NODE explicitly.

set -euo pipefail

die() {
  printf 'error: %s\n' "$*" >&2
  exit 2
}

require_env() {
  local name="$1"
  [[ -n "${!name:-}" ]] || die "required environment variable $name is unset"
}

run_timed() {
  local timing_log="$1"
  shift
  "$PYTHON" - "$timing_log" "$@" <<'PY'
from pathlib import Path
import subprocess
import sys
import time

timing_log = Path(sys.argv[1])
command = sys.argv[2:]
start_ns = time.monotonic_ns()
completed = subprocess.run(command, check=False)
elapsed_ns = time.monotonic_ns() - start_ns
timing_log.write_text(f"monotonic_elapsed_ns\t{elapsed_ns}\n")
if elapsed_ns < 0:
    sys.exit(2)
if completed.returncode < 0:
    sys.exit(128 - completed.returncode)
sys.exit(completed.returncode)
PY
}

read_elapsed_ns() {
  local timing_log="$1"
  local elapsed_ns
  elapsed_ns="$(awk -F '\t' '$1 == "monotonic_elapsed_ns" { print $2 }' \
    "$timing_log")"
  [[ "$elapsed_ns" =~ ^[0-9]+$ ]] || return 1
  printf '%s\n' "$elapsed_ns"
}

append_index() {
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$ROW" "$ACTION" "$ATTEMPT_ID" "$1" "$2" "$3" "$4" "$5" "$6" "$7" \
    >> "$INDEX"
}

reject_duplicate_launch() {
  local sample="$1"
  if awk -F '\t' -v row="$ROW" -v action="$ACTION" -v attempt="$ATTEMPT_ID" \
      -v sample="$sample" \
      'NR > 1 && $1 == row && $2 == action && $3 == attempt && $4 == sample {
         found = 1
       }
       END { exit(found ? 0 : 1) }' "$INDEX"; then
    die "packet index already contains $ROW $ACTION $ATTEMPT_ID $sample"
  fi
}

run_accounting_hook() {
  local phase="$1"
  local sample="$2"
  local output_dir="$3"
  local log="$PACKET/logs/$ROW.$ACTION.$ATTEMPT_ID.$sample.fs-accounting.$phase.log"
  if [[ -z "${FS_ACCOUNTING_HOOK:-}" ]]; then
    printf 'INCOMPLETE: FS_ACCOUNTING_HOOK is unset\n' > "$log"
    return 3
  fi
  if [[ ! -x "$FS_ACCOUNTING_HOOK" ]]; then
    printf 'INCOMPLETE: FS_ACCOUNTING_HOOK is not executable: %s\n' \
      "$FS_ACCOUNTING_HOOK" > "$log"
    return 3
  fi
  "$FS_ACCOUNTING_HOOK" "$phase" "$ROW" "$sample" "$output_dir" \
    > "$log" 2>&1
}

scan_forbidden_staging() {
  local sample="$1"
  local output_dir="$2"
  local inventory="$PACKET/inventory/$ROW.$ACTION.$ATTEMPT_ID.$sample.staging-scan.txt"
  local -a roots=("$PACKET" "$output_dir")
  if [[ "$ACTION" == resume && -n "$MANIFEST" ]]; then
    roots+=("$(dirname "$MANIFEST")")
  fi
  find "${roots[@]}" -type f \
    \( -name '*.assembled' -o -name '*.assembled.tmp' \) -print \
    | sort -u > "$inventory"
  [[ ! -s "$inventory" ]]
}

set_topology() {
  SRUN_TOPOLOGY=()
  case "$ROW" in
    ED-1|MR-1)
      NODES=2
      RANKS=4
      RANKS_PER_NODE=2
      SRUN_TOPOLOGY=(--nodes=2 --ntasks=4 --ntasks-per-node=2
                     --distribution=block:block)
      ;;
    MR-2)
      NODES=2
      RANKS=3
      RANKS_PER_NODE=asymmetric
      require_env MR2_HOSTFILE
      [[ -f "$MR2_HOSTFILE" ]] || die "MR2_HOSTFILE does not exist: $MR2_HOSTFILE"
      mr2_hosts=()
      while IFS= read -r host || [[ -n "$host" ]]; do
        mr2_hosts+=("$host")
      done < "$MR2_HOSTFILE"
      [[ "${#mr2_hosts[@]}" -eq 3 ]] ||
        die "MR2_HOSTFILE must contain exactly three host lines"
      [[ "${mr2_hosts[0]}" != "${mr2_hosts[1]}" ]] ||
        die "MR2_HOSTFILE must place rank 0 on host A and rank 1 on host B"
      [[ "${mr2_hosts[1]}" == "${mr2_hosts[2]}" ]] ||
        die "MR2_HOSTFILE must place ranks 1 and 2 on host B"
      SRUN_TOPOLOGY=(--nodes=2 --ntasks=3 --distribution=arbitrary)
      ;;
    scaling)
      require_env NODES
      require_env RANKS
      require_env RANKS_PER_NODE
      SRUN_TOPOLOGY=(--nodes="$NODES" --ntasks="$RANKS"
                     --ntasks-per-node="$RANKS_PER_NODE"
                     --distribution=block:block)
      ;;
    *)
      die "unsupported row: $ROW"
      ;;
  esac
}

run_srun() {
  local sample="$1"
  local measured="$2"
  local output_dir="$3"
  shift 3
  local stdout="$PACKET/logs/$ROW.$ACTION.$ATTEMPT_ID.$sample.%J.%t.out"
  local stderr="$PACKET/logs/$ROW.$ACTION.$ATTEMPT_ID.$sample.%J.%t.err"
  local timing_log="$PACKET/logs/$ROW.$ACTION.$ATTEMPT_ID.$sample.launch-timing.tsv"
  local elapsed_ns=0 exit_code=0 disposition=passed
  local before_status=0 after_status=0
  local -a launcher=(srun "${SRUN_TOPOLOGY[@]}" --label
                     --output="$stdout" --error="$stderr" "$@")

  reject_duplicate_launch "$sample"
  mkdir -p "$output_dir"
  if [[ "$measured" == 1 ]]; then
    set +e
    run_accounting_hook before "$sample" "$output_dir"
    before_status=$?
    set -e
    if [[ "$before_status" -ne 0 ]]; then
      disposition=incomplete-accounting-before
      append_index "$sample" "$measured" "$before_status" "$disposition" \
        "$elapsed_ns" "$output_dir" "${FS_ACCOUNTING_HOOK:-INCOMPLETE-unset}"
      return "$before_status"
    fi
  fi
  set +e
  if [[ "$ROW" == MR-2 ]]; then
    run_timed "$timing_log" env SLURM_HOSTFILE="$MR2_HOSTFILE" \
      ATHENAK_RESTART_MANIFEST_TIMING="$measured" \
      timeout "$ROW_TIMEOUT" "${launcher[@]}"
  else
    run_timed "$timing_log" env ATHENAK_RESTART_MANIFEST_TIMING="$measured" \
      timeout "$ROW_TIMEOUT" "${launcher[@]}"
  fi
  exit_code=$?
  set -e
  if ! elapsed_ns="$(read_elapsed_ns "$timing_log")"; then
    elapsed_ns=0
    exit_code=2
    disposition=failed-invalid-launch-timing
  elif [[ "$exit_code" -eq 124 ]]; then
    disposition=timeout-incomplete
  elif [[ "$exit_code" -ne 0 ]]; then
    disposition=failed
  else
    disposition=passed
  fi
  if [[ "$measured" == 1 ]]; then
    set +e
    run_accounting_hook after "$sample" "$output_dir"
    after_status=$?
    set -e
    if [[ "$after_status" -ne 0 ]]; then
      disposition=incomplete-accounting-after
      [[ "$exit_code" -ne 0 ]] && disposition=failed-and-incomplete-accounting-after
    fi
  fi
  find "$output_dir" -print | sort \
    > "$PACKET/inventory/$ROW.$ACTION.$ATTEMPT_ID.$sample.files.txt"
  if ! scan_forbidden_staging "$sample" "$output_dir"; then
    disposition=failed-forbidden-assembled
    [[ "$exit_code" -eq 0 ]] && exit_code=2
  elif [[ "$after_status" -ne 0 && "$exit_code" -eq 0 ]]; then
    exit_code="$after_status"
  fi
  append_index "$sample" "$measured" "$exit_code" "$disposition" \
    "$elapsed_ns" "$output_dir" "${FS_ACCOUNTING_HOOK:-INCOMPLETE-unset}"
  [[ "$exit_code" -eq 0 ]] || exit "$exit_code"
}

run_rank_map() {
  local sample=rank-map
  local stdout="$PACKET/logs/$ROW.$ACTION.$ATTEMPT_ID.$sample.%J.%t.out"
  local stderr="$PACKET/logs/$ROW.$ACTION.$ATTEMPT_ID.$sample.%J.%t.err"
  local timing_log="$PACKET/logs/$ROW.$ACTION.$ATTEMPT_ID.$sample.launch-timing.tsv"
  local elapsed_ns=0 exit_code=0 disposition
  local -a launcher=(srun "${SRUN_TOPOLOGY[@]}" --label
                     --output="$stdout" --error="$stderr" hostname)

  reject_duplicate_launch "$sample"
  set +e
  if [[ "$ROW" == MR-2 ]]; then
    run_timed "$timing_log" env SLURM_HOSTFILE="$MR2_HOSTFILE" \
      timeout "$ROW_TIMEOUT" "${launcher[@]}"
  else
    run_timed "$timing_log" timeout "$ROW_TIMEOUT" "${launcher[@]}"
  fi
  exit_code=$?
  set -e
  if ! elapsed_ns="$(read_elapsed_ns "$timing_log")"; then
    elapsed_ns=0
    exit_code=2
    disposition=failed-invalid-rank-map-timing
  elif [[ "$exit_code" -eq 124 ]]; then
    disposition=timeout-rank-map
  elif [[ "$exit_code" -ne 0 ]]; then
    disposition=failed-rank-map
  else
    disposition=passed-rank-map
  fi
  append_index "$sample" 0 "$exit_code" "$disposition" "$elapsed_ns" \
    "N/A" "N/A"
  [[ "$exit_code" -eq 0 ]] || exit "$exit_code"
}

require_env PACKET
require_env REPO
require_env BUILD
require_env PYTHON
command -v srun >/dev/null || die "srun is unavailable"
command -v timeout >/dev/null || die "timeout is unavailable"
[[ -x "$BUILD/src/athena" ]] || die "AthenaK executable is unavailable: $BUILD/src/athena"

ROW="${1:-}"
ACTION="${2:-}"
MANIFEST="${3:-}"
ROW_TIMEOUT="${ROW_TIMEOUT:-300}"
ATTEMPT_ID="${ATTEMPT_ID:-attempt-1}"
INDEX="$PACKET/packet-index.tsv"
DECK="$REPO/tst/inputs/io_node_sharding.athinput"

[[ "$ACTION" == generate || "$ACTION" == resume ]] ||
  die "action must be generate or resume"
[[ "$ATTEMPT_ID" =~ ^[A-Za-z0-9._-]+$ ]] ||
  die "ATTEMPT_ID must contain only letters, digits, dot, underscore, or hyphen"
[[ -f "$DECK" ]] || die "frozen correctness deck is unavailable: $DECK"
mkdir -p "$PACKET"/{decks,logs,outputs,inventory}
cp -p "$DECK" "$PACKET/decks/io_node_sharding.athinput"
if [[ ! -e "$INDEX" ]]; then
  printf 'row\taction\tattempt_id\tsample\tmeasured\texit_code\tdisposition\tmonotonic_elapsed_ns\toutput_dir\tfs_accounting_hook\n' \
    > "$INDEX"
fi

set_topology
run_rank_map

case "$ROW:$ACTION" in
  ED-1:generate)
    run_srun run 0 "$PACKET/outputs/$ROW/$ATTEMPT_ID/generate" \
      "$BUILD/src/athena" -i "$DECK" -d "$PACKET/outputs/$ROW/$ATTEMPT_ID/generate" \
      mesh/nx1=32 output1/single_file_per_node=true \
      output2/single_file_per_node=true output3/single_file_per_node=true \
      output4/single_file_per_node=true output5/single_file_per_node=true \
      output6/single_file_per_node=true time/final_output_policy=none
    ;;
  MR-1:generate|MR-2:generate)
    nx1=32
    [[ "$ROW" == MR-2 ]] && nx1=24
    run_srun run 0 "$PACKET/outputs/$ROW/$ATTEMPT_ID/generate" \
      "$BUILD/src/athena" -i "$DECK" -d "$PACKET/outputs/$ROW/$ATTEMPT_ID/generate" \
      mesh/nx1="$nx1" output1/dt=-1 output2/dt=-1 output3/dt=-1 \
      output4/dt=-1 output5/dt=-1 output6/single_file_per_node=true \
      time/final_output_policy=none
    ;;
  MR-1:resume|MR-2:resume|scaling:resume)
    [[ -f "$MANIFEST" ]] || die "resume manifest does not exist: $MANIFEST"
    run_srun warmup 0 "$PACKET/outputs/$ROW/$ATTEMPT_ID/resume/warmup" \
      "$BUILD/src/athena" -r "$MANIFEST" \
      -d "$PACKET/outputs/$ROW/$ATTEMPT_ID/resume/warmup" \
      time/tlim=0 time/nlim=0 time/final_output_policy=none \
      time/output_timing=false
    for sample in 1 2 3 4 5; do
      run_srun "sample-$sample" 1 \
        "$PACKET/outputs/$ROW/$ATTEMPT_ID/resume/sample-$sample" \
        "$BUILD/src/athena" -r "$MANIFEST" \
        -d "$PACKET/outputs/$ROW/$ATTEMPT_ID/resume/sample-$sample" \
        time/tlim=0 time/nlim=0 time/final_output_policy=none \
        time/output_timing=false
    done
    ;;
  *)
    die "unsupported row/action combination: $ROW:$ACTION"
    ;;
esac

printf 'completed row=%s action=%s index=%s\n' "$ROW" "$ACTION" "$INDEX"
