#!/bin/bash
# Slurm spools this script, so receive the installed runner as an absolute argv.
set -euo pipefail

RUNNER="${1:?usage: run_installed_control_plane_job.sh INSTALLED_RUNNER [ARG ...]}"
shift
RUNNER_PREFIX="/lustre/orion/ast207/proj-shared/dfielding/PIC/control_plane/"
case "$RUNNER" in
  "${RUNNER_PREFIX}"*/run_control_plane.py)
    version="${RUNNER#"${RUNNER_PREFIX}"}"
    version="${version%/run_control_plane.py}"
    ;;
  *)
    printf 'Refusing non-installed control-plane runner: %s\n' "$RUNNER" >&2
    exit 1
    ;;
esac
[[ "$version" =~ ^[0-9a-f]{64}$ ]] || {
  printf 'Refusing non-installed control-plane runner: %s\n' "$RUNNER" >&2
  exit 1
}
exec /opt/cray/pe/python/3.11.7/bin/python3 -I "$RUNNER" "$@"
