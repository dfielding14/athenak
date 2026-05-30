#!/usr/bin/env bash
# Finalize one external IO qualification evidence packet without checksum cycles.

set -euo pipefail

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

[[ "$#" -eq 2 ]] ||
  die "usage: finalize_external_io_qualification_packet.sh PACKET ARCHIVE_RECORD"

PACKET="$1"
ARCHIVE_RECORD="$2"
MANIFEST="$PACKET/artifacts.sha256"
INDEX="$PACKET/packet-index.md"
RUNNER_INDEX="$PACKET/packet-index.tsv"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
INDEX_VALIDATOR="$SCRIPT_DIR/validate_external_io_packet_index.py"
tmp_manifest=
tmp_index=

[[ -d "$PACKET" ]] || die "packet directory does not exist: $PACKET"
[[ ! -L "$PACKET" ]] || die "packet directory must not be a symlink"
[[ -f "$RUNNER_INDEX" && ! -L "$RUNNER_INDEX" ]] ||
  die "runner packet index must be a regular non-symlink file: $RUNNER_INDEX"
[[ -f "$INDEX_VALIDATOR" ]] || die "packet-index validator is unavailable: $INDEX_VALIDATOR"
packet_symlinks="$(find "$PACKET" -type l -print)"
[[ -z "$packet_symlinks" ]] || die "packet must not contain symlinks"
packet_hardlinks="$(find "$PACKET" -type f -links +1 -print)"
[[ -z "$packet_hardlinks" ]] || die "packet must not contain hard-linked files"
packet_control_paths="$(find "$PACKET" -name '*[[:cntrl:]]*' -print)"
[[ -z "$packet_control_paths" ]] ||
  die "packet artifact paths must not contain control characters"
packet_abs="$(cd "$PACKET" && pwd -P)"
archive_dir="$(dirname "$ARCHIVE_RECORD")"
[[ -d "$archive_dir" ]] || die "archive-record directory does not exist: $archive_dir"
archive_dir_abs="$(cd "$archive_dir" && pwd -P)"
archive_abs="$archive_dir_abs/$(basename "$ARCHIVE_RECORD")"
[[ ! "$packet_abs" =~ [[:cntrl:]] ]] || die "packet path must not contain control characters"
[[ "$packet_abs" != *'|'* && "$packet_abs" != *'`'* ]] ||
  die "packet path must not contain Markdown table delimiters"
[[ ! "$archive_abs" =~ [[:cntrl:]] ]] ||
  die "archive-record path must not contain control characters"
case "$archive_abs" in
  "$packet_abs"|"$packet_abs"/*)
    die "archive record must live outside the immutable packet"
    ;;
esac

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
    die "packet directory identity changed during finalization"
}

validate_archive_dir_identity() {
  [[ -d "$archive_dir_abs" && ! -L "$archive_dir_abs" ]] ||
    die "archive-record directory must remain a regular non-symlink directory"
  [[ "$(stat_identity "$archive_dir_abs")" == "$ARCHIVE_DIR_ID" ]] ||
    die "archive-record directory identity changed during finalization"
}

scan_packet_tree_safety() {
  local unsafe
  validate_packet_root_identity
  unsafe="$(find "$PACKET" -type l -print -o -type f -links +1 -print -o \
    -name '*[[:cntrl:]]*' -print -o \
    \( -name '.artifacts.sha256.tmp.*' -o -name '.packet-index.md.tmp.*' -o \
      -name '.packet-index.tsv.tmp.*' \) \
    -print | head -n 1)"
  [[ -z "$unsafe" ]] ||
    die "packet contains an unsafe alias or reserved helper path: $unsafe"
}

preflight_archive_record() {
  local observed_identity
  validate_archive_dir_identity
  [[ ! -L "$ARCHIVE_RECORD" ]] || die "archive record must not be a symlink"
  if [[ -e "$ARCHIVE_RECORD" ]]; then
    [[ -f "$ARCHIVE_RECORD" ]] || die "archive record must be a regular file"
    [[ "$(stat_links "$ARCHIVE_RECORD")" -eq 1 ]] ||
      die "archive record must not have hard-link aliases"
  fi
  observed_identity="$(
    "$INDEX_VALIDATOR" --ensure-archive-record "$ARCHIVE_RECORD" "$ARCHIVE_DIR_ID"
  )" || die "archive record descriptor admission failed"
  if [[ -z "${ARCHIVE_RECORD_ID:-}" ]]; then
    ARCHIVE_RECORD_ID="$observed_identity"
  fi
  [[ "$observed_identity" == "$ARCHIVE_RECORD_ID" ]] ||
    die "archive record identity changed during finalization"
}

validate_runner_index() {
  [[ "$(stat_links "$RUNNER_INDEX")" -eq 1 ]] ||
    die "runner packet index must not have hard-link aliases"
  "$INDEX_VALIDATOR" --require-rows "$RUNNER_INDEX"
}

scan_packet_write_bits() {
  local output="$1"
  local errors="$2"
  validate_packet_root_identity
  find "$PACKET" \( -perm -0200 -o -perm -0020 -o -perm -0002 \) -print \
    > "$output" 2> "$errors"
}

packet_has_write_bits() {
  local output errors status
  output="$(mktemp "${TMPDIR:-/tmp}/athenak-io-write-bits.XXXXXX")"
  errors="$output.stderr"
  scan_packet_write_bits "$output" "$errors" ||
    die "packet write-permission scan failed during metadata recovery"
  status=1
  [[ -s "$output" ]] && status=0
  rm -f "$output" "$errors"
  return "$status"
}

reset_writable_metadata_pair() {
  packet_has_write_bits ||
    die "inconsistent finalization metadata exists in a read-only packet"
  "$INDEX_VALIDATOR" --reset-metadata-pair "$MANIFEST" "$INDEX" "$PACKET_ID"
}

publish_metadata_pair() {
  "$INDEX_VALIDATOR" --publish-file "$tmp_manifest" "$MANIFEST" "$PACKET_ID"
  "$INDEX_VALIDATOR" --publish-file "$tmp_index" "$INDEX" "$PACKET_ID"
  scan_packet_tree_safety
}

write_canonical_manifest() {
  local output="$1"
  scan_packet_tree_safety
  (
    cd "$PACKET"
    find . -type f ! -path ./artifacts.sha256 ! -path ./packet-index.md -print0 |
      LC_ALL=C sort -z |
      xargs -0 shasum -a 256
  ) > "$output"
  scan_packet_tree_safety
}

PACKET_ID="$(stat_identity "$PACKET")"
ARCHIVE_DIR_ID="$(stat_identity "$archive_dir_abs")"
ARCHIVE_RECORD_ID=

cleanup() {
  rm -f "$tmp_manifest" "$tmp_index"
}
trap cleanup EXIT

scan_packet_tree_safety
preflight_archive_record
validate_runner_index
if [[ -e "$MANIFEST" && ! -e "$INDEX" ]] ||
    [[ ! -e "$MANIFEST" && -e "$INDEX" ]]; then
  reset_writable_metadata_pair
fi

tmp_manifest="$(mktemp "${TMPDIR:-/tmp}/athenak-io-artifacts.XXXXXX")"
write_canonical_manifest "$tmp_manifest"
manifest_digest="$(shasum -a 256 "$tmp_manifest" | awk '{print $1}')"
tmp_index="$(mktemp "${TMPDIR:-/tmp}/athenak-io-packet-index.XXXXXX")"
cat > "$tmp_index" <<EOF
# AthenaK External IO Qualification Packet

| Field | Value |
| --- | --- |
| Packet root | \`$packet_abs\` |
| Inner artifact manifest | \`artifacts.sha256\` |
| Inner artifact-manifest SHA-256 | \`$manifest_digest\` |
| Checksum rule | The inner manifest excludes itself and this packet index. The archive record outside the packet retains the outer index digest. |
EOF
if [[ ! -e "$MANIFEST" && ! -e "$INDEX" ]]; then
  publish_metadata_pair
elif [[ ! -f "$MANIFEST" || ! -f "$INDEX" ]]; then
  die "finalization metadata must be regular files"
elif ! cmp -s "$tmp_manifest" "$MANIFEST" ||
    ! cmp -s "$tmp_index" "$INDEX"; then
  reset_writable_metadata_pair
  publish_metadata_pair
fi
"$INDEX_VALIDATOR" --sync-directory "$PACKET" "$PACKET_ID"
rm -f "$tmp_manifest"
tmp_manifest=
rm -f "$tmp_index"
tmp_index=

manifest_digest="$(shasum -a 256 "$MANIFEST" | awk '{print $1}')"
(cd "$PACKET" && shasum -a 256 -c artifacts.sha256 >/dev/null) ||
  die "packet artifact checksum verification failed"
index_digest="$(shasum -a 256 "$INDEX" | awk '{print $1}')"

chmod -R a-w "$PACKET"
scan_packet_tree_safety
write_bits="$(mktemp "${TMPDIR:-/tmp}/athenak-io-write-bits.XXXXXX")"
write_errors="$write_bits.stderr"
scan_packet_write_bits "$write_bits" "$write_errors" ||
  die "packet write-permission scan failed after chmod"
[[ ! -s "$write_bits" ]] || die "packet remains writable after chmod"
rm -f "$write_bits" "$write_errors"
(cd "$PACKET" && shasum -a 256 -c artifacts.sha256 >/dev/null) ||
  die "packet artifact checksum verification failed after chmod"
tmp_manifest="$(mktemp "${TMPDIR:-/tmp}/athenak-io-artifacts.XXXXXX")"
write_canonical_manifest "$tmp_manifest"
cmp -s "$tmp_manifest" "$MANIFEST" ||
  die "packet inventory changed after write-permission removal"
rm -f "$tmp_manifest"
tmp_manifest=
scan_packet_tree_safety
preflight_archive_record
scan_packet_tree_safety
"$INDEX_VALIDATOR" --append-archive-row "$ARCHIVE_RECORD" "$ARCHIVE_DIR_ID" \
  "$ARCHIVE_RECORD_ID" "$PACKET_ID" "$packet_abs" artifacts.sha256 "$manifest_digest" \
  packet-index.md "$index_digest"

printf 'packet=%s\n' "$packet_abs"
printf 'artifacts_sha256=%s\n' "$manifest_digest"
printf 'packet_index_sha256=%s\n' "$index_digest"
