#!/bin/bash
# Apply the inventory-bound Frontier profile, then exec the trusted argv.
set -euo pipefail

export HOME=/
export USER="$(/usr/bin/id -un)"
source /opt/cray/pe/lmod/lmod/init/profile
unset BASH_ENV ENV
: "${PIC_FRONTIER_PROFILE_FD:?PIC_FRONTIER_PROFILE_FD is required}"
: "${PIC_RUNTIME_ALLOWLIST_FD:?PIC_RUNTIME_ALLOWLIST_FD is required}"
: "${PIC_RUNTIME_ALLOWLIST_DIR_FD:?PIC_RUNTIME_ALLOWLIST_DIR_FD is required}"
case "${PIC_FRONTIER_PROFILE_FD}:${PIC_RUNTIME_ALLOWLIST_FD}:${PIC_RUNTIME_ALLOWLIST_DIR_FD}" in
  *[!0-9:]*|:*|*:)
    printf 'PIC runtime bindings must be file descriptors\n' >&2
    exit 1
    ;;
esac
if ! source "/proc/self/fd/${PIC_FRONTIER_PROFILE_FD}"; then
  exit 1
fi
record_pic_environment >&"$PIC_RUNTIME_ALLOWLIST_FD"
/opt/cray/pe/python/3.11.7/bin/python3 -E -s - \
  "$PIC_RUNTIME_ALLOWLIST_FD" "$PIC_RUNTIME_ALLOWLIST_DIR_FD" <<'PY'
import os
import stat
import sys

allowlist_fd = int(sys.argv[1])
directory_fd = int(sys.argv[2])
if not stat.S_ISREG(os.fstat(allowlist_fd).st_mode):
    raise SystemExit("PIC runtime allowlist is not a regular file")
if not stat.S_ISDIR(os.fstat(directory_fd).st_mode):
    raise SystemExit("PIC runtime allowlist parent is not a directory")
os.fsync(allowlist_fd)
os.fchmod(allowlist_fd, 0o400)
os.fsync(allowlist_fd)
os.fsync(directory_fd)
PY
eval "exec ${PIC_FRONTIER_PROFILE_FD}>&-"
eval "exec ${PIC_RUNTIME_ALLOWLIST_FD}>&-"
eval "exec ${PIC_RUNTIME_ALLOWLIST_DIR_FD}>&-"
unset PIC_FRONTIER_PROFILE_FD PIC_RUNTIME_ALLOWLIST_FD PIC_RUNTIME_ALLOWLIST_DIR_FD
unset -f module ml clearMT clearLmod xSetTitleLmod 2>/dev/null || true
exec "$@"
