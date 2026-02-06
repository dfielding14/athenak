#!/usr/bin/env bash

set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

if [ "$#" -gt 1 ]; then
  echo "Usage: $0 [base_ref]"
  exit 2
fi

base_ref="${1:-}"

changed_py_files=()
while IFS= read -r file; do
  if [ -n "$file" ]; then
    changed_py_files+=("$file")
  fi
done < <(
  {
    if [ -n "$base_ref" ]; then
      git diff --name-only --diff-filter=ACMRTUXB "${base_ref}...HEAD"
    fi
    git diff --name-only --diff-filter=ACMRTUXB
    git diff --cached --name-only --diff-filter=ACMRTUXB
  } | awk 'NF' | sort -u | grep -E '^(tst|vis)/.*\.py$' || true
)

if [ "${#changed_py_files[@]}" -eq 0 ]; then
  echo "No changed Python files under tst/ or vis/ to lint."
  exit 0
fi

echo "Linting changed Python files:"
printf '  %s\n' "${changed_py_files[@]}"

flake8 "${changed_py_files[@]}"

echo "Changed-file Python style checks passed."
