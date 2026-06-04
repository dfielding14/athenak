#!/bin/bash
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --job-name=pic-q011-pressure-review
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

export PATH=/usr/bin:/bin
export PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR=/tmp
unset BASH_ENV ENV GIT_ALTERNATE_OBJECT_DIRECTORIES GIT_CONFIG_COUNT
unset GIT_CONFIG_GLOBAL GIT_CONFIG_KEY_0 GIT_CONFIG_NOSYSTEM GIT_CONFIG_SYSTEM
unset GIT_CONFIG_VALUE_0 GIT_DIR GIT_INDEX_FILE GIT_OBJECT_DIRECTORY GIT_WORK_TREE
unset PYTHONHOME PYTHONINSPECT PYTHONPATH PYTHONSTARTUP PYTHONUSERBASE

REPO_ROOT=/autofs/nccs-svm1_home2/dfielding/athenak-pic
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
EXPECTED_GIT_COMMIT="${1:?usage: sbatch $0 FULL_GIT_COMMIT}"
[[ "$EXPECTED_GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]]
GIT=(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null
)
SOURCE_CLOSURE=(
  tst/publication/frontier_q011_section54_pressure_pilot_publish_job.sh
  tst/publication/publish_q011_section54_pressure_pilot_bundle.py
  tst/publication/analyze_q011_section54_pressure_pilot.py
  tst/publication/analyze_q011_section54_pressure_pilot_case.py
  tst/publication/analyze_q011_section54_outputs.py
  tst/publication/frontier_f1_structured_artifacts.py
  tst/publication/pvtk_particles.py
  tst/publication/render_q011_section54_pressure_pilot_review_packet.py
  tst/publication/frontier_q011_section54_pressure_pilot_review_packet_job.sh
  tst/publication/readiness/q011_section54_pressure_pilot_preregistration_2026-06-01.json
  tst/publication/readiness/q011_section54_pressure_pilot_aggregate_analysis_compatibility_successor_2026-06-02.json
  tst/publication/readiness/q011_section54_pressure_pilot_registered_execution_preregistration_2026-06-02.json
  tst/publication/readiness/q011_section54_pressure_pilot_registered_execution_retry_successor_v2_2026-06-02.json
  tst/publication/readiness/q011_section54_pressure_pilot_postrun_aggregate_source_authorization_successor_2026-06-02.json
  tst/publication/readiness/q011_section54_pressure_pilot_postrun_aggregate_source_authorization_successor_v2_2026-06-03.json
  tst/publication/readiness/q011_section54_pressure_pilot_postrun_aggregate_source_authorization_successor_v3_2026-06-04.json
  tst/publication/readiness/plotting_environment_lock_candidate_2026-05-30.json
)

cd "$REPO_ROOT"
"${GIT[@]}" diff --quiet "$EXPECTED_GIT_COMMIT" -- "${SOURCE_CLOSURE[@]}"
UNTRACKED_SOURCE=$(
  "${GIT[@]}" ls-files --others --exclude-standard -- "${SOURCE_CLOSURE[@]}"
)
test -z "$UNTRACKED_SOURCE"
SOURCE_COMMIT=$("${GIT[@]}" rev-parse HEAD)
test "$SOURCE_COMMIT" = "$EXPECTED_GIT_COMMIT"
test "$("${GIT[@]}" rev-parse origin/PIC)" = "$EXPECTED_GIT_COMMIT"
echo "source_commit=$SOURCE_COMMIT"
/usr/bin/sha256sum "${SOURCE_CLOSURE[@]}"

SNAPSHOT_ROOT=$(/usr/bin/mktemp -d "${TMPDIR}/pic-q011-pressure-review.XXXXXXXX")
SOURCE_ARCHIVE=$(/usr/bin/mktemp "${TMPDIR}/pic-q011-pressure-review-source.XXXXXXXX.tar")
trap '/usr/bin/chmod -R u+w "$SNAPSHOT_ROOT" 2>/dev/null || true; /usr/bin/rm -rf "$SNAPSHOT_ROOT"; /usr/bin/rm -f "$SOURCE_ARCHIVE"' EXIT
"${GIT[@]}" archive "$EXPECTED_GIT_COMMIT" > "$SOURCE_ARCHIVE"
/usr/bin/chmod 0444 "$SOURCE_ARCHIVE"
export PIC_PRESSURE_PUBLICATION_SOURCE_ARCHIVE_PATH="$SOURCE_ARCHIVE"
export PIC_PRESSURE_PUBLICATION_SOURCE_SNAPSHOT_ROOT="$SNAPSHOT_ROOT"
SOURCE_ARCHIVE_SHA256=$(/usr/bin/sha256sum "$SOURCE_ARCHIVE" | /usr/bin/awk '{print $1}')
echo "source_archive_sha256=$SOURCE_ARCHIVE_SHA256"
test "$("${GIT[@]}" get-tar-commit-id < "$SOURCE_ARCHIVE")" = "$SOURCE_COMMIT"
/usr/bin/tar -xf "$SOURCE_ARCHIVE" -C "$SNAPSHOT_ROOT"
/usr/bin/chmod -R a-w "$SNAPSHOT_ROOT"
cd "$SNAPSHOT_ROOT/tst/publication"

run_snapshot_python() {
  local script="${1:?usage: run_snapshot_python SCRIPT [ARG ...]}"
  shift
  "$PYTHON" -I -B -c \
    'import runpy, sys; root, script, *args = sys.argv[1:]; sys.path.insert(0, root); sys.argv = [script, *args]; runpy.run_path(script, run_name="__main__")' \
    "$PWD" "${PWD}/${script}" "$@"
}

run_snapshot_plot_python() {
  local script="${1:?usage: run_snapshot_plot_python SCRIPT [ARG ...]}"
  shift
  "$PYTHON" -I -B -c \
    'import importlib.metadata, json, runpy, sys; root, script, package_root, lock_path, *args = sys.argv[1:]; sys.path.insert(0, package_root); lock = json.load(open(lock_path, encoding="utf-8")); actual_python = sys.version.split()[0]; expected_python = lock["python"]; actual_python == expected_python or sys.exit(f"plot Python drifted: {actual_python} != {expected_python}"); [(importlib.metadata.version(name) == version) or sys.exit(f"plot dependency drifted: {name}") for name, version in lock["dependencies"].items()]; sys.path.insert(0, root); sys.argv = [script, *args]; runpy.run_path(script, run_name="__main__")' \
    "$PWD" "${PWD}/${script}" \
    /autofs/nccs-svm1_home2/dfielding/.local/lib/python3.11/site-packages \
    "${PWD}/readiness/plotting_environment_lock_candidate_2026-05-30.json" "$@"
}

run_snapshot_plot_python render_q011_section54_pressure_pilot_review_packet.py
