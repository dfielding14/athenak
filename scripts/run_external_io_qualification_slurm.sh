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
import os
import stat
import subprocess
import sys
import time

timing_log = Path(sys.argv[1])
command = sys.argv[2:]
directory_flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
directory_flags |= getattr(os, "O_NOFOLLOW", 0)
packet_fd = os.open(timing_log.parent.parent, directory_flags)
packet_metadata = os.fstat(packet_fd)
directory_fd = os.open(timing_log.parent, directory_flags)
directory_metadata = os.fstat(directory_fd)
start_ns = time.monotonic_ns()
completed = subprocess.run(command, check=False)
elapsed_ns = time.monotonic_ns() - start_ns
try:
    current_packet_metadata = timing_log.parent.parent.lstat()
    current_directory_metadata = timing_log.parent.lstat()
except OSError:
    sys.exit(2)
if (
    not stat.S_ISDIR(current_packet_metadata.st_mode)
    or not os.path.samestat(packet_metadata, current_packet_metadata)
    or
    not stat.S_ISDIR(current_directory_metadata.st_mode)
    or not os.path.samestat(directory_metadata, current_directory_metadata)
):
    sys.exit(2)
file_flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
file_flags |= getattr(os, "O_NOFOLLOW", 0)
descriptor = os.open(timing_log.name, file_flags, 0o644, dir_fd=directory_fd)
with os.fdopen(descriptor, "w") as stream:
    stream.write(f"monotonic_elapsed_ns\t{elapsed_ns}\n")
os.close(directory_fd)
os.close(packet_fd)
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
  elapsed_ns="$(awk -F '\t' '
    NR == 1 && NF == 2 && $1 == "monotonic_elapsed_ns" &&
        $2 ~ /^[0-9]+$/ {
      elapsed_ns = $2
      valid = 1
      next
    }
    {
      valid = 0
      exit 1
    }
    END {
      if (NR != 1 || !valid) {
        exit 1
      }
      print elapsed_ns
    }
  ' "$timing_log")" || return 1
  [[ "$elapsed_ns" =~ ^[0-9]+$ ]] || return 1
  printf '%s\n' "$elapsed_ns"
}

append_index() {
  validate_serialized_fields "$@"
  scan_packet_tree_safety
  validate_packet_index_file
  INDEX_ID="$("$INDEX_VALIDATOR" --append-row "$INDEX" "$PACKET_ID" "$INDEX_ID" \
    "$ROW" "$ACTION" "$ATTEMPT_ID" "$1" "$2" "$3" "$4" "$5" "$6" "$7")"
}

validate_serialized_fields() {
  local field
  for field in "$ROW" "$ACTION" "$ATTEMPT_ID" "$@"; do
    [[ -n "$field" && ! "$field" =~ [[:cntrl:]] ]] ||
      die "packet index fields must be nonempty and contain no control characters"
  done
}

reject_control_characters() {
  local field
  for field in "$@"; do
    [[ ! "$field" =~ [[:cntrl:]] ]] ||
      die "qualification paths and hooks must contain no control characters"
  done
}

stat_links() {
  if stat -f '%l' "$1" >/dev/null 2>&1; then
    stat -f '%l' "$1"
  else
    stat -c '%h' "$1"
  fi
}

stat_identity() {
  if stat -f '%d:%i' "$1" >/dev/null 2>&1; then
    stat -f '%d:%i' "$1"
  else
    stat -c '%d:%i' "$1"
  fi
}

validate_packet_root_identity() {
  [[ -d "$PACKET" && ! -L "$PACKET" ]] ||
    die "packet directory must remain a regular non-symlink directory"
  [[ "$(stat_identity "$PACKET")" == "$PACKET_ID" ]] ||
    die "packet directory identity changed during qualification"
}

validate_directory_identity() {
  local label="$1"
  local path="$2"
  local expected_identity="$3"
  [[ -d "$path" && ! -L "$path" ]] ||
    die "$label must remain a regular non-symlink directory: $path"
  [[ "$(stat_identity "$path")" == "$expected_identity" ]] ||
    die "$label identity changed during qualification: $path"
}

validate_fixed_packet_children() {
  validate_directory_identity "packet decks directory" "$PACKET/decks" "$DECKS_ID"
  validate_directory_identity "packet logs directory" "$PACKET/logs" "$LOGS_ID"
  validate_directory_identity "packet outputs directory" "$PACKET/outputs" "$OUTPUTS_ID"
  validate_directory_identity "packet inventory directory" "$PACKET/inventory" \
    "$INVENTORY_ID"
}

validate_packet_index_file() {
  [[ -f "$INDEX" && ! -L "$INDEX" ]] ||
    die "packet index must be a regular non-symlink file: $INDEX"
  [[ "$(stat_links "$INDEX")" -eq 1 ]] ||
    die "packet index must not have hard-link aliases: $INDEX"
  if [[ -n "${INDEX_ID:-}" ]]; then
    [[ "$(stat_identity "$INDEX")" == "$INDEX_ID" ]] ||
      die "packet index identity changed during qualification: $INDEX"
  fi
  "$INDEX_VALIDATOR" "$INDEX"
}

scan_packet_tree_safety() {
  local packet_symlinks packet_hardlinks packet_control_paths packet_temp_paths
  validate_packet_root_identity
  validate_fixed_packet_children
  if ! packet_symlinks="$(find "$PACKET" -type l -print)"; then
    die "packet symlink scan failed"
  fi
  [[ -z "$packet_symlinks" ]] ||
    die "packet directory tree must not contain symlinks"
  if ! packet_hardlinks="$(find "$PACKET" -type f -links +1 -print)"; then
    die "packet hard-link scan failed"
  fi
  [[ -z "$packet_hardlinks" ]] ||
    die "packet directory tree must not contain hard-linked files"
  if ! packet_control_paths="$(find "$PACKET" -name '*[[:cntrl:]]*' -print)"; then
    die "packet control-character path scan failed"
  fi
  [[ -z "$packet_control_paths" ]] ||
    die "packet artifact paths must not contain control characters"
  if ! packet_temp_paths="$(find "$PACKET" -name '.packet-index.tsv.tmp.*' -print)"; then
    die "packet helper-temporary scan failed"
  fi
  [[ -z "$packet_temp_paths" ]] ||
    die "packet directory tree contains reserved helper-temporary paths"
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
  local errors="$inventory.stderr"
  local -a roots=("$PACKET" "$output_dir")
  if [[ "$ACTION" == resume && -n "$MANIFEST" ]]; then
    roots+=("$(dirname "$MANIFEST")")
  fi
  if ! find "${roots[@]}" -type f \
      \( -name '*.assembled' -o -name '*.assembled.tmp' \) -print \
      2> "$errors" | LC_ALL=C sort -u > "$inventory"; then
    return 2
  fi
  [[ ! -s "$inventory" ]]
}

retain_output_inventory() {
  local sample="$1"
  local output_dir="$2"
  local inventory="$PACKET/inventory/$ROW.$ACTION.$ATTEMPT_ID.$sample.files.txt"
  local errors="$inventory.stderr"
  find "$output_dir" -print 2> "$errors" | LC_ALL=C sort > "$inventory"
}

validate_rank_map() {
  local rank_map_dir="$1"
  local inventory="$PACKET/inventory/$ROW.$ACTION.$ATTEMPT_ID.rank-map.tsv"
  local errors="$inventory.stderr"
  local validation_errors validation_status
  [[ ! -e "$inventory" && ! -L "$inventory" ]] ||
    die "rank-map inventory publication target already exists: $inventory"
  [[ ! -e "$errors" && ! -L "$errors" ]] ||
    die "rank-map error publication target already exists: $errors"
  validation_errors="$(mktemp "${TMPDIR:-/tmp}/athenak-rank-map-errors.XXXXXX")"
  if [[ "$ROW" == MR-2 ]]; then
    set -- "$rank_map_dir" "$inventory" "$ROW" "$RANKS" "$NODES" \
      "$RANKS_PER_NODE" "${mr2_hosts[@]}"
  else
    set -- "$rank_map_dir" "$inventory" "$ROW" "$RANKS" "$NODES" \
      "$RANKS_PER_NODE"
  fi
  set +e
  "$PYTHON" - "$@" 2> "$validation_errors" <<'PY'
from collections import Counter
import os
from pathlib import Path
import re
import stat
import sys

rank_map_dir = Path(sys.argv[1])
inventory = Path(sys.argv[2])
row = sys.argv[3]
ranks = int(sys.argv[4])
nodes = int(sys.argv[5])
ranks_per_node = sys.argv[6]
expected_hosts = sys.argv[7:]
hostname_pattern = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")
rank_pattern = re.compile(r"rank-(0|[1-9][0-9]*)\.tsv\Z")
records = {}

if not stat.S_ISDIR(rank_map_dir.lstat().st_mode):
    raise SystemExit(f"rank-map directory must be a regular directory: {rank_map_dir}")
for path in rank_map_dir.iterdir():
    match = rank_pattern.fullmatch(path.name)
    metadata = path.lstat()
    if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1 or match is None:
        raise SystemExit(f"unexpected rank-map artifact: {path.name}")
    lines = path.read_text().splitlines()
    if len(lines) != 1:
        raise SystemExit(f"rank-map artifact must contain one line: {path.name}")
    fields = lines[0].split("\t")
    if len(fields) != 2 or not fields[0].isdigit():
        raise SystemExit(f"malformed rank-map artifact: {path.name}")
    rank = int(fields[0])
    host = fields[1]
    if rank != int(match.group(1)) or rank in records:
        raise SystemExit(f"rank-map rank mismatch or duplicate: {path.name}")
    if hostname_pattern.fullmatch(host) is None:
        raise SystemExit(f"noncanonical rank-map hostname: {host!r}")
    records[rank] = host

if sorted(records) != list(range(ranks)):
    raise SystemExit(f"rank-map ranks mismatch: {sorted(records)}")
hosts = [records[rank] for rank in range(ranks)]
if row in ("ED-1", "MR-1"):
    if not (hosts[0] == hosts[1] and hosts[2] == hosts[3] and hosts[0] != hosts[2]):
        raise SystemExit(f"rank-map topology mismatch: {hosts}")
elif row == "MR-2":
    if hosts != expected_hosts:
        raise SystemExit(f"rank-map MR-2 placement mismatch: {hosts} != {expected_hosts}")
elif row == "scaling":
    if not ranks_per_node.isdigit():
        raise SystemExit(f"non-numeric scaling ranks-per-node: {ranks_per_node}")
    counts = Counter(hosts)
    if len(counts) != nodes or set(counts.values()) != {int(ranks_per_node)}:
        raise SystemExit(f"rank-map scaling topology mismatch: {counts}")
else:
    raise SystemExit(f"unsupported rank-map row: {row}")

payload = "rank\thost\n" + "".join(
    f"{rank}\t{records[rank]}\n" for rank in range(ranks)
)
descriptor = os.open(inventory, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
with os.fdopen(descriptor, "w") as stream:
    stream.write(payload)
PY
  validation_status=$?
  set -e
  [[ ! -e "$errors" && ! -L "$errors" ]] || {
    rm -f "$validation_errors"
    die "rank-map error publication target already exists: $errors"
  }
  mv "$validation_errors" "$errors"
  return "$validation_status"
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
      for host in "${mr2_hosts[@]}"; do
        [[ "$host" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] ||
          die "MR2_HOSTFILE host lines must be canonical hostname tokens"
      done
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
  local elapsed_ns=0 exit_code=0 disposition=passed output_dir_id
  local before_status=0 after_status=0 inventory_status=0 scan_status=0
  local -a launcher=(srun "${SRUN_TOPOLOGY[@]}" --label
                     --output="$stdout" --error="$stderr" "$@")

  scan_packet_tree_safety
  reject_duplicate_launch "$sample"
  mkdir -p "$output_dir"
  output_dir_id="$(stat_identity "$output_dir")"
  validate_directory_identity "launch output directory" "$output_dir" "$output_dir_id"
  if [[ "$measured" == 1 ]]; then
    set +e
    run_accounting_hook before "$sample" "$output_dir"
    before_status=$?
    set -e
    scan_packet_tree_safety
    validate_directory_identity "launch output directory" "$output_dir" "$output_dir_id"
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
  scan_packet_tree_safety
  validate_directory_identity "launch output directory" "$output_dir" "$output_dir_id"
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
    scan_packet_tree_safety
    validate_directory_identity "launch output directory" "$output_dir" "$output_dir_id"
  fi
  validate_directory_identity "launch output directory" "$output_dir" "$output_dir_id"
  set +e
  retain_output_inventory "$sample" "$output_dir"
  inventory_status=$?
  scan_forbidden_staging "$sample" "$output_dir"
  scan_status=$?
  set -e
  scan_packet_tree_safety
  validate_directory_identity "launch output directory" "$output_dir" "$output_dir_id"
  if [[ "$inventory_status" -ne 0 ]]; then
    disposition=failed-inventory-scan
    [[ "$exit_code" -eq 0 ]] && exit_code=2
  elif [[ "$scan_status" -eq 2 ]]; then
    disposition=failed-staging-scan
    [[ "$exit_code" -eq 0 ]] && exit_code=2
  elif [[ "$scan_status" -ne 0 ]]; then
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
  local rank_map_dir="$PACKET/inventory/$ROW.$ACTION.$ATTEMPT_ID.rank-map"
  local elapsed_ns=0 exit_code=0 disposition rank_map_dir_id
  local rank_map_script='
set -eu
rank="${SLURM_PROCID:?}"
case "$rank" in ""|*[!0-9]*) exit 2;; esac
host="$(hostname)"
case "$host" in ""|*[!A-Za-z0-9._-]*) exit 2;; esac
printf "%s\t%s\n" "$rank" "$host" > "$1/rank-$rank.tsv"
'
  local -a launcher=(srun "${SRUN_TOPOLOGY[@]}" --label
                     --output="$stdout" --error="$stderr"
                     bash -c "$rank_map_script" rank-map-task "$rank_map_dir")

  scan_packet_tree_safety
  reject_duplicate_launch "$sample"
  mkdir -p "$rank_map_dir"
  rank_map_dir_id="$(stat_identity "$rank_map_dir")"
  validate_directory_identity "rank-map directory" "$rank_map_dir" "$rank_map_dir_id"
  set +e
  if [[ "$ROW" == MR-2 ]]; then
    run_timed "$timing_log" env SLURM_HOSTFILE="$MR2_HOSTFILE" \
      timeout "$ROW_TIMEOUT" "${launcher[@]}"
  else
    run_timed "$timing_log" timeout "$ROW_TIMEOUT" "${launcher[@]}"
  fi
  exit_code=$?
  set -e
  scan_packet_tree_safety
  validate_directory_identity "rank-map directory" "$rank_map_dir" "$rank_map_dir_id"
  if ! elapsed_ns="$(read_elapsed_ns "$timing_log")"; then
    elapsed_ns=0
    exit_code=2
    disposition=failed-invalid-rank-map-timing
  elif [[ "$exit_code" -eq 124 ]]; then
    disposition=timeout-rank-map
  elif [[ "$exit_code" -ne 0 ]]; then
    disposition=failed-rank-map
  elif ! validate_rank_map "$rank_map_dir"; then
    exit_code=2
    disposition=failed-invalid-rank-map
  else
    disposition=passed-rank-map
  fi
  validate_directory_identity "rank-map directory" "$rank_map_dir" "$rank_map_dir_id"
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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
INDEX_VALIDATOR="$SCRIPT_DIR/validate_external_io_packet_index.py"

[[ "$ACTION" == generate || "$ACTION" == resume ]] ||
  die "action must be generate or resume"
[[ "$ATTEMPT_ID" =~ ^[A-Za-z0-9._-]+$ ]] ||
  die "ATTEMPT_ID must contain only letters, digits, dot, underscore, or hyphen"
[[ -f "$DECK" ]] || die "frozen correctness deck is unavailable: $DECK"
[[ -f "$INDEX_VALIDATOR" ]] || die "packet-index validator is unavailable: $INDEX_VALIDATOR"
[[ ! -L "$PACKET" ]] || die "packet directory must not be a symlink"
reject_control_characters "$PACKET" "${FS_ACCOUNTING_HOOK:-INCOMPLETE-unset}"
mkdir -p "$PACKET"/{decks,logs,outputs,inventory}
PACKET_ID="$(stat_identity "$PACKET")"
DECKS_ID="$(stat_identity "$PACKET/decks")"
LOGS_ID="$(stat_identity "$PACKET/logs")"
OUTPUTS_ID="$(stat_identity "$PACKET/outputs")"
INVENTORY_ID="$(stat_identity "$PACKET/inventory")"
INDEX_ID=
scan_packet_tree_safety
cp -p "$DECK" "$PACKET/decks/io_node_sharding.athinput"
if [[ ! -e "$INDEX" ]]; then
  printf 'row\taction\tattempt_id\tsample\tmeasured\texit_code\tdisposition\tmonotonic_elapsed_ns\toutput_dir\tfs_accounting_hook\n' \
    > "$INDEX"
fi
INDEX_ID="$(stat_identity "$INDEX")"
validate_packet_index_file

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
