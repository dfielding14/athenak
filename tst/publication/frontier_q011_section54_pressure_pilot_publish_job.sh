#!/bin/bash
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=pic-q011-pressure-publish
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

SNAPSHOT_ROOT=$(/usr/bin/mktemp -d "${TMPDIR}/pic-q011-pressure-publish.XXXXXXXX")
SOURCE_ARCHIVE=$(/usr/bin/mktemp "${TMPDIR}/pic-q011-pressure-publish-source.XXXXXXXX.tar")
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

run_snapshot_python publish_q011_section54_pressure_pilot_bundle.py \
  "$PIC_ROOT/publication/q011_section54_pressure_pilot_bundle" \
  --receipt-path \
  "$PIC_ROOT/publication/q011_section54_pressure_pilot_bundle_receipt.json" \
  --analysis-result-path \
  "$PIC_ROOT/publication/q011_section54_pressure_pilot_analysis.json" \
  --case-artifact-dir \
  "ps_p0_1p00=$PIC_ROOT/runs/q011_section54_pressure_ps_p0_1p00/a9db403f-a174-45ae-956d-f96185a94db1" \
  --case-descriptor-sha256 \
  "ps_p0_1p00=da322586da4831c7747e3a2f5774adf6b999835dc8636b9d971411e1dc9471b2" \
  --case-artifact-dir \
  "ps_p0_0p05=$PIC_ROOT/runs/q011_section54_pressure_ps_p0_0p05/3eb84dfa-0be7-4a62-b77b-64131c08dc6f" \
  --case-descriptor-sha256 \
  "ps_p0_0p05=978cfc63efc9503b9dd0a043793d759e930f390340b843a39d621c2b859ff85c" \
  --case-artifact-dir \
  "ps_p0_0p10=$PIC_ROOT/runs/q011_section54_pressure_ps_p0_0p10/07bd4d96-ffc8-4446-98ef-562d4d83da0c" \
  --case-descriptor-sha256 \
  "ps_p0_0p10=dba35e7b157c1646444792cfda68963dbfd2062b4887961bfb380e51a113f5d4" \
  --case-artifact-dir \
  "ps_p0_0p20=$PIC_ROOT/runs/q011_section54_pressure_ps_p0_0p20/8b4065a3-b02a-4312-8238-20d72d958bbb" \
  --case-descriptor-sha256 \
  "ps_p0_0p20=ef842e1440cbe1fd51e3ee79af5c838704e782708e2a827ec0bfcc89ca1ae312"

run_snapshot_python publish_q011_section54_pressure_pilot_bundle.py \
  --verify-published-receipt \
  "$PIC_ROOT/publication/q011_section54_pressure_pilot_bundle_receipt.json"
