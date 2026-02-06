#!/usr/bin/env bash

set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

if [ "$#" -gt 1 ]; then
  echo "Usage: $0 [base_ref]"
  exit 2
fi

base_ref="${1:-}"

changed_cpp_files=()
while IFS= read -r file; do
  if [ -n "$file" ]; then
    changed_cpp_files+=("$file")
  fi
done < <(
  {
    if [ -n "$base_ref" ]; then
      git diff --name-only --diff-filter=ACMRTUXB "${base_ref}...HEAD"
    fi
    git diff --name-only --diff-filter=ACMRTUXB
    git diff --cached --name-only --diff-filter=ACMRTUXB
  } | awk 'NF' | sort -u | grep -E '^src/.*\.(cpp|hpp)$' || true
)

if [ "${#changed_cpp_files[@]}" -eq 0 ]; then
  echo "No changed C++ files under src/ to lint."
  exit 0
fi

echo "Linting changed C++ files:"
printf '  %s\n' "${changed_cpp_files[@]}"

cpplint_py="${repo_root}/tst/scripts/style/cpplint.py"
python3 -u "$cpplint_py" \
  --filter=-build/include_subdir \
  --counting=detailed \
  "${changed_cpp_files[@]}"

for file in "${changed_cpp_files[@]}"; do
  if grep -n $'\t' "$file"; then
    echo "ERROR: Do not use tab characters in $file"
    exit 1
  fi

  if grep -nri "}}" "$file" | grep -v "//" >/dev/null; then
    echo "ERROR: Use one closing brace per line in $file"
    exit 1
  fi

  if grep -nrEi '^\s+#pragma' "$file" >/dev/null; then
    echo "ERROR: Left justify #pragma statements in $file"
    exit 1
  fi

  if grep -n -E ' +$' "$file" >/dev/null; then
    echo "ERROR: Found trailing whitespace in $file"
    exit 1
  fi
done

echo "Changed-file C++ style checks passed."
