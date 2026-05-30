#!/bin/bash
# Apply the inventory-bound Frontier profile, then exec the trusted argv.
set -euo pipefail

CONTROL_PLANE_DIR="$(cd "$(/usr/bin/dirname "$0")" && pwd)"
if ! type module >/dev/null 2>&1; then
  source /etc/profile
fi
source "${CONTROL_PLANE_DIR}/frontier_pic_environment.sh"
exec "$@"
