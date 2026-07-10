#!/bin/bash
#SBATCH -A AST207
#SBATCH -J gotham_pdf_movies
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH -t 06:00:00

set -euo pipefail

: "${PNG_DIR:?PNG_DIR must name the published PNG directory}"
: "${MOVIE_DIR:?MOVIE_DIR must name the movie output directory}"
: "${MOVIE_SCRIPT:?MOVIE_SCRIPT must name make_rebuilt_pdf_movies.py}"
: "${FFMPEG:?FFMPEG must name the static ffmpeg executable}"
: "${FFPROBE:?FFPROBE must name the static ffprobe executable}"

exec /opt/cray/pe/python/3.11.7/bin/python3 "$MOVIE_SCRIPT" \
  "$PNG_DIR" \
  --output-dir "$MOVIE_DIR" \
  --ffmpeg "$FFMPEG" \
  --ffprobe "$FFPROBE" \
  --fps 5 \
  --minimum-duration 3 \
  --max-dimension 1440 \
  --workers 16 \
  --ffmpeg-threads 4 \
  --preset medium \
  --crf 20 \
  --overwrite
