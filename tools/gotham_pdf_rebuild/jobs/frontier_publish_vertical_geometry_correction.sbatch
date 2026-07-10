#!/bin/bash
#SBATCH -A AST207
#SBATCH -p service
#SBATCH --qos=develop
#SBATCH -J gotham_vertical_publish
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH -t 02:00:00

set -euo pipefail

: "${STAGING_ROOT:?STAGING_ROOT must name completed correction staging}"
: "${FINAL_ROOT:?FINAL_ROOT must name production plot root}"
: "${PUBLISHER:?PUBLISHER must name publish_vertical_geometry_correction.py}"
: "${REPATHER:?REPATHER must name repath_rebuilt_pdf_plot_metadata.py}"
: "${MOVIE_SCRIPT:?MOVIE_SCRIPT must name make_rebuilt_pdf_movies.py}"
: "${FFMPEG:?FFMPEG must name static ffmpeg}"
: "${FFPROBE:?FFPROBE must name static ffprobe}"

/opt/cray/pe/python/3.11.7/bin/python3 "$PUBLISHER" \
  "$STAGING_ROOT" \
  --final-root "$FINAL_ROOT" \
  --snapshot-count 161

/opt/cray/pe/python/3.11.7/bin/python3 "$REPATHER" \
  "$FINAL_ROOT/metadata" \
  --campaign-root "$FINAL_ROOT/../.campaign_v2/staging" \
  --final-root "$FINAL_ROOT" \
  --workers 28

rm -rf "$FINAL_ROOT/movies"
/opt/cray/pe/python/3.11.7/bin/python3 "$MOVIE_SCRIPT" \
  "$FINAL_ROOT/png" \
  --output-dir "$FINAL_ROOT/movies" \
  --ffmpeg "$FFMPEG" \
  --ffprobe "$FFPROBE" \
  --fps 5 \
  --minimum-duration 3 \
  --max-dimension 1440 \
  --workers 14 \
  --ffmpeg-threads 4 \
  --preset medium \
  --crf 20 \
  --overwrite
