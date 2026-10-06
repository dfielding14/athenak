#!/bin/bash
set -euo pipefail
cd /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research
module restore
module load PrgEnv-cray cpe/25.09 cray-mpich/9.0.1 cce/20.0.0 cray-python/3.11.7
module unload darshan-runtime
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$PWD/tmp"
python3 amr_cap_probe.py
