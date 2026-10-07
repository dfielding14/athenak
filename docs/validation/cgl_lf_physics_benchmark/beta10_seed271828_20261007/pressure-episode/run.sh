#!/bin/bash
set -euo pipefail
cd /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$PWD/analysis-checks/pressure-episode/tmp"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=1
mkdir -p "$TMPDIR"
srun --jobid=5629677 --nodes=1 --ntasks=1 --cpus-per-task=8 --gpus=0 --exact --overlap \
  python3 analysis-checks/pressure-episode/analyze_pressure_episode.py \
  --source /autofs/nccs-svm1_home2/dfielding/athenak-cgl-wo2 \
  --metrics analysis-fixed-6-14/metrics.json \
  --output analysis-checks/pressure-episode
