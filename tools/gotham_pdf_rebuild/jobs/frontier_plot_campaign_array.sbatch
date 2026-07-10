#!/bin/bash
#SBATCH -A AST207
#SBATCH -J gotham_pdf_plots
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -t 01:30:00

set -euo pipefail

: "${SNAPSHOT_LIST:?SNAPSHOT_LIST must name the ordered production snapshot list}"
: "${STAGING_ROOT:?STAGING_ROOT must name the isolated plot staging root}"
: "${PLOTTER:?PLOTTER must name plot_rebuilt_pdfs.py}"

index="${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array}"
input_dir="$(sed -n "$((index + 1))p" "$SNAPSHOT_LIST")"
if [[ -z "$input_dir" || ! -d "$input_dir" ]]; then
  echo "Invalid snapshot list entry ${index}: ${input_dir}" >&2
  exit 2
fi

output_dir="$STAGING_ROOT/$index"
mkdir -p "$output_dir" "$STAGING_ROOT/mpl/$index"
export MPLCONFIGDIR="$STAGING_ROOT/mpl/$index"

exec /opt/cray/pe/python/3.11.7/bin/python3 "$PLOTTER" \
  "$input_dir" \
  --output-dir "$output_dir" \
  --render-profile production \
  --workers "${SLURM_CPUS_PER_TASK:-16}" \
  --no-gallery
