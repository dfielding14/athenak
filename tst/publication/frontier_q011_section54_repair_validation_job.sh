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
export PYTHONPATH="$PWD:$PWD/tst/publication/frontier_control_plane"

cd /ccs/home/dfielding/athenak-pic
export PYTHONPATH="$PWD:$PWD/tst/publication/frontier_control_plane"

git diff --check

mapfile -t changed_python < <(
  {
    git diff --name-only -- '*.py'
    git ls-files --others --exclude-standard -- '*.py'
  } | sort -u
)

python3 -B - "${changed_python[@]}" <<'PY'
import sys
from pathlib import Path

for relative in sys.argv[1:]:
    path = Path(relative)
    compile(path.read_bytes(), str(path), "exec")
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
