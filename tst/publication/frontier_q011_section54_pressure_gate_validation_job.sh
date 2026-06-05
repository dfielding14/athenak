#!/bin/bash
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=pic-q011-pressure-gate-validate
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

export PATH=/opt/cray/pe/python/3.11.7/bin:/usr/bin:/bin
export PYTHONDONTWRITEBYTECODE=1
unset BASH_ENV ENV GIT_ALTERNATE_OBJECT_DIRECTORIES GIT_CONFIG_COUNT
unset GIT_CONFIG_GLOBAL GIT_CONFIG_KEY_0 GIT_CONFIG_NOSYSTEM GIT_CONFIG_SYSTEM
unset GIT_CONFIG_VALUE_0 GIT_DIR GIT_INDEX_FILE GIT_OBJECT_DIRECTORY GIT_WORK_TREE
unset PYTHONHOME PYTHONINSPECT PYTHONPATH PYTHONSTARTUP PYTHONUSERBASE

REPO_ROOT=/autofs/nccs-svm1_home2/dfielding/athenak-pic
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
EXPECTED_GIT_COMMIT="${1:?usage: sbatch $0 FULL_GIT_COMMIT}"
[[ "$EXPECTED_GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]]
cd "$REPO_ROOT"
SOURCE_STATUS=$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    status --porcelain --untracked-files=all
)
test -z "$SOURCE_STATUS"
SOURCE_COMMIT=$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null rev-parse HEAD
)
test "$SOURCE_COMMIT" = "$EXPECTED_GIT_COMMIT"
test "$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    rev-parse origin/PIC
)" = "$EXPECTED_GIT_COMMIT"
echo "source_commit=$SOURCE_COMMIT"

SNAPSHOT_ROOT=$(/usr/bin/mktemp -d "${TMPDIR:-/tmp}/pic-q011-pressure-gate-validate.XXXXXXXX")
trap '/usr/bin/rm -rf "$SNAPSHOT_ROOT"' EXIT
/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
  /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null \
  archive "$SOURCE_COMMIT" | /usr/bin/tar -xf - -C "$SNAPSHOT_ROOT"
cd "$SNAPSHOT_ROOT"
export PYTHONPATH="$PWD:$PWD/tst/publication/frontier_control_plane"

publication_python_output=$(/usr/bin/find tst/publication -type f -name '*.py' | /usr/bin/sort)
publication_shell_output=$(/usr/bin/find tst/publication -type f -name '*.sh' | /usr/bin/sort)
publication_json_output=$(/usr/bin/find tst/publication -type f -name '*.json' | /usr/bin/sort)
test -n "$publication_python_output"
test -n "$publication_shell_output"
test -n "$publication_json_output"
mapfile -t publication_python <<< "$publication_python_output"
mapfile -t publication_shell <<< "$publication_shell_output"
mapfile -t publication_json <<< "$publication_json_output"
test "${#publication_python[@]}" -eq 145
test "${#publication_shell[@]}" -eq 17
test "${#publication_json[@]}" -eq 290

"$PYTHON" -B - "${publication_python[@]}" <<'PY'
import sys
from pathlib import Path

for relative in sys.argv[1:]:
    path = Path(relative)
    compile(path.read_bytes(), str(path), "exec")
PY

/usr/bin/bash -n "${publication_shell[@]}"

"$PYTHON" -B - "${publication_json[@]}" <<'PY'
import json
import sys
from pathlib import Path

for relative in sys.argv[1:]:
    path = Path(relative)
    with path.open(encoding="utf-8") as stream:
        json.load(stream)
PY

"$PYTHON" -B - <<'PY'
import json
from pathlib import Path

from tst.publication.frontier_control_plane import generate_prepared_pic_artifact_inventory

path = Path("tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json")
with path.open(encoding="utf-8") as stream:
    committed = json.load(stream)
generated = generate_prepared_pic_artifact_inventory.prepared_artifact_inventory()
if committed != generated:
    raise SystemExit("Committed prepared PIC artifact inventory is stale")
PY

"$PYTHON" -B - <<'PY'
import fcntl
import os

from tst.publication import publish_q011_section54_pressure_pilot_bundle as publisher

pic_root, _publication_root = publisher._publication_root(publisher.AUTHORIZED_PIC_ROOT)
acceptance_root = publisher._publication_acceptance_root(pic_root)
anchor = publisher._publication_transaction_anchor(pic_root)
anchor_descriptor = publisher._open_absolute_directory(anchor)
acceptance_descriptor = publisher._open_absolute_directory(acceptance_root)
try:
    publisher._lock_publication_transaction(anchor_descriptor, acceptance_descriptor)
    publisher._require_same_directory(
        anchor, anchor_descriptor, "stable publication transaction anchor"
    )
    publisher._require_same_directory(
        acceptance_root,
        acceptance_descriptor,
        "authorized PIC publication acceptance root",
    )
    for path in (anchor, acceptance_root):
        competitor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        try:
            try:
                fcntl.flock(competitor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                pass
            else:
                raise SystemExit(f"publication transaction lock did not serialize: {path}")
        finally:
            os.close(competitor)
finally:
    os.close(acceptance_descriptor)
    os.close(anchor_descriptor)
PY

"$PYTHON" -B -m unittest \
  tst.publication.frontier_control_plane.test_q011_pressure_review_packet_verifier \
  tst.publication.test_q011_section54_pressure_selection \
  tst.publication.test_q011_section54_qualifying_campaign_execution \
  tst.publication.test_q011_section54_attempt_manifest_materializer \
  tst.publication.test_analyze_q011_section54_campaign \
  tst.publication.test_analyze_q011_section54_numerical_qualification \
  tst.publication.test_pic_qualification_manifest \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_q011_helper_source_order_matches_execution_and_analyzer \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_missing_pressure_packet_binding \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_pressure_packet_hash_drift \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_cross_bound_pressure_packet_aggregate \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_pressure_bundle_manifest_hash_drift \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_pressure_aggregate_analysis_hash_drift \
  tst.publication.frontier_control_plane.test_control_plane.SnapshotTests.test_planner_retention_rejects_pressure_aggregate_descriptor_drift

cd "$REPO_ROOT"
test "$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null rev-parse HEAD
)" = "$SOURCE_COMMIT"
SOURCE_STATUS=$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    status --porcelain --untracked-files=all
)
test -z "$SOURCE_STATUS"
TAXONOMY_DRY_RUN_SEED="$REPO_ROOT/tst/build/src/bin/pic_bell_pub_taxonomy.taxonomy.bin"
test ! -e "$TAXONOMY_DRY_RUN_SEED"
export PYTHONPATH="$PWD:$PWD/tst/publication/frontier_control_plane"

modules_output=$(
  /usr/bin/find tst/publication -type f -name 'test_*.py' |
    /usr/bin/sed -e 's#/#.#g' -e 's#\.py$##' |
    /usr/bin/sort
)
test -n "$modules_output"
mapfile -t modules <<< "$modules_output"
test "${#modules[@]}" -eq 65

"$PYTHON" -B -m unittest "${modules[@]}"
test "$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null rev-parse HEAD
)" = "$SOURCE_COMMIT"
SOURCE_STATUS=$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    status --porcelain --untracked-files=all
)
test -z "$SOURCE_STATUS"
test ! -e "$TAXONOMY_DRY_RUN_SEED"
