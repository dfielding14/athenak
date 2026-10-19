#!/bin/bash
#SBATCH -A AST207
#SBATCH -p service
#SBATCH --qos=develop
#SBATCH -J gotham_pdf_plots
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH -t 02:00:00

set -euo pipefail

: "${SNAPSHOT_LIST:?SNAPSHOT_LIST must name the ordered production snapshot list}"
: "${STAGING_ROOT:?STAGING_ROOT must name the isolated plot staging root}"
: "${PLOTTER:?PLOTTER must name plot_rebuilt_pdfs.py}"
: "${CAMPAIGN_LOG_ROOT:?CAMPAIGN_LOG_ROOT must name the per-snapshot log root}"

mkdir -p "$STAGING_ROOT" "$STAGING_ROOT/mpl" "$CAMPAIGN_LOG_ROOT"

render_one() {
  local index="$1"
  local input_dir output_dir
  input_dir="$(sed -n "$((index + 1))p" "$SNAPSHOT_LIST")"
  output_dir="$STAGING_ROOT/$index"
  if [[ -f "$output_dir/metadata.json" ]]; then
    echo "SKIP completed snapshot $index"
    return 0
  fi
  rm -rf "$output_dir"
  mkdir -p "$output_dir" "$STAGING_ROOT/mpl/$index"
  (
    export MPLCONFIGDIR="$STAGING_ROOT/mpl/$index"
    /opt/cray/pe/python/3.11.7/bin/python3 "$PLOTTER" \
      "$input_dir" \
      --output-dir "$output_dir" \
      --render-profile production \
      --workers 8 \
      --no-gallery
  ) >"$CAMPAIGN_LOG_ROOT/$index.out" 2>"$CAMPAIGN_LOG_ROOT/$index.err"
  echo "DONE snapshot $index"
}
export -f render_one
export SNAPSHOT_LIST STAGING_ROOT PLOTTER CAMPAIGN_LOG_ROOT

snapshot_count="$(wc -l < "$SNAPSHOT_LIST")"
seq 0 "$((snapshot_count - 1))" | xargs -n 1 -P 7 bash -c 'render_one "$1"' _

completed="$(find "$STAGING_ROOT" -mindepth 2 -maxdepth 2 -name metadata.json | wc -l)"
if [[ "$completed" -ne "$snapshot_count" ]]; then
  echo "Only $completed of $snapshot_count snapshots completed" >&2
  exit 1
fi
echo "Completed all $completed staged snapshots"
