#!/bin/bash
set -euo pipefail
cd /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/passive-cast-proof
module restore
module load PrgEnv-cray cpe/25.09 cray-mpich/9.0.1 cce/20.0.0 cray-python/3.11.7
module load craype-accel-amd-gfx90a rocm/6.4.2
module unload darshan-runtime
unset CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH
export INCLUDE_PATH_X86_64="$(CC --print-resource-dir)/include:${CC_X86_64}/include/craylibs"
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH:-}"
export TMPDIR="$PWD/tmp" PYTHONDONTWRITEBYTECODE=1
backend="${1:?cpu or hip}"
if [[ "$backend" == cpu ]]; then module unload craype-accel-amd-gfx90a; fi
python3 compile_proof.py "$backend"
