#!/bin/bash
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --job-name=pic-q011-pressure-review
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

export PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
export PYTHONDONTWRITEBYTECODE=1

REPO_ROOT=/ccs/home/dfielding/athenak-pic
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
)

cd "$REPO_ROOT"
git diff --quiet HEAD -- "${SOURCE_CLOSURE[@]}"
test -z "$(git ls-files --others --exclude-standard -- "${SOURCE_CLOSURE[@]}")"
SOURCE_COMMIT=$(git rev-parse HEAD)
echo "source_commit=$SOURCE_COMMIT"
sha256sum "${SOURCE_CLOSURE[@]}"

SNAPSHOT_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/pic-q011-pressure-review.XXXXXXXX")
SOURCE_ARCHIVE=$(mktemp "${TMPDIR:-/tmp}/pic-q011-pressure-review-source.XXXXXXXX.tar")
trap 'chmod -R u+w "$SNAPSHOT_ROOT" 2>/dev/null || true; rm -rf "$SNAPSHOT_ROOT"; rm -f "$SOURCE_ARCHIVE"' EXIT
git archive HEAD > "$SOURCE_ARCHIVE"
chmod 0444 "$SOURCE_ARCHIVE"
export PIC_PRESSURE_PUBLICATION_SOURCE_ARCHIVE_PATH="$SOURCE_ARCHIVE"
export PIC_PRESSURE_PUBLICATION_SOURCE_SNAPSHOT_ROOT="$SNAPSHOT_ROOT"
echo "source_archive_sha256=$(sha256sum "$SOURCE_ARCHIVE" | awk '{print $1}')"
test "$(git get-tar-commit-id < "$SOURCE_ARCHIVE")" = "$SOURCE_COMMIT"
tar -xf "$SOURCE_ARCHIVE" -C "$SNAPSHOT_ROOT"
chmod -R a-w "$SNAPSHOT_ROOT"
cd "$SNAPSHOT_ROOT/tst/publication"
python3 -B render_q011_section54_pressure_pilot_review_packet.py
