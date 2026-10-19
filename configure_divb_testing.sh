#!/usr/bin/env bash

set -euo pipefail

script_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
exec "${script_dir}/configure_and_build_gotham.sh" \
  --problem built_in_pgens \
  --build-dir build_divb_testing \
  "$@"
