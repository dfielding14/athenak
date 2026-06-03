#!/bin/bash
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=pic-q011-repair-validate
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

export PYTHONDONTWRITEBYTECODE=1

REPO_ROOT=/ccs/home/dfielding/athenak-pic
cd "$REPO_ROOT"
test -z "$(git status --porcelain --untracked-files=all)"
echo "source_commit=$(git rev-parse HEAD)"

SNAPSHOT_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/pic-q011-repair-validate.XXXXXXXX")
trap 'rm -rf "$SNAPSHOT_ROOT"' EXIT
git archive HEAD | tar -xf - -C "$SNAPSHOT_ROOT"
cd "$SNAPSHOT_ROOT"
export PYTHONPATH="$PWD:$PWD/tst/publication/frontier_control_plane"

mapfile -t publication_python < <(find tst/publication -type f -name '*.py' | sort)
mapfile -t publication_shell < <(find tst/publication -type f -name '*.sh' | sort)
mapfile -t publication_json < <(find tst/publication -type f -name '*.json' | sort)

python3 -B - "${publication_python[@]}" <<'PY'
import sys
from pathlib import Path

for relative in sys.argv[1:]:
    path = Path(relative)
    compile(path.read_bytes(), str(path), "exec")
PY

bash -n "${publication_shell[@]}"

python3 -B - "${publication_json[@]}" <<'PY'
import json
import sys
from pathlib import Path

for relative in sys.argv[1:]:
    path = Path(relative)
    with path.open(encoding="utf-8") as stream:
        json.load(stream)
PY

python3 -B - <<'PY'
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

python3 -B -m unittest \
  tst.publication.test_analyze_q011_section54_campaign \
  tst.publication.test_q011_section54_model \
  tst.publication.test_q011_section54_particles \
  tst.publication.test_q011_section54_spatial \
  tst.publication.test_q011_section54_restart \
  tst.publication.test_q011_section54_artifacts \
  tst.publication.test_q011_section54_attempt_manifest_materializer \
  tst.publication.test_q011_section54_qualifying_campaign_preregistration \
  tst.publication.test_q011_section54_pressure_selection \
  tst.publication.test_q011_section54_qualifying_campaign_execution \
  tst.publication.test_publish_q011_section54_campaign_attempt \
  tst.publication.test_analyze_q011_section54_numerical_qualification

mapfile -t modules < <(
  rg --files tst/publication -g 'test_*.py' |
    sed -e 's#/#.#g' -e 's#\.py$##' |
    sort
)

python3 -B -m unittest "${modules[@]}"
