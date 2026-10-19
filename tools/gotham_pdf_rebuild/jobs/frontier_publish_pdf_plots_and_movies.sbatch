#!/bin/bash
#SBATCH -A AST207
#SBATCH -p service
#SBATCH --qos=develop
#SBATCH -J gotham_pdf_publish
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH -t 02:00:00

set -euo pipefail

: "${STAGING_ROOT:?STAGING_ROOT must name the completed plot staging root}"
: "${FINAL_ROOT:?FINAL_ROOT must name the empty final production plot root}"
: "${ASSEMBLER:?ASSEMBLER must name assemble_rebuilt_pdf_plots.py}"
: "${MOVIE_SCRIPT:?MOVIE_SCRIPT must name make_rebuilt_pdf_movies.py}"
: "${FFMPEG:?FFMPEG must name the static ffmpeg executable}"
: "${FFPROBE:?FFPROBE must name the static ffprobe executable}"

/opt/cray/pe/python/3.11.7/bin/python3 "$ASSEMBLER" \
  "$STAGING_ROOT" \
  --output-dir "$FINAL_ROOT" \
  --snapshot-count 161 \
  --expected-products 64 \
  --expected-plots 239

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
