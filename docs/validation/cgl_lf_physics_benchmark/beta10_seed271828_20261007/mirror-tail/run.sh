#!/bin/bash
set -euo pipefail
cd /lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/physics-benchmark
source scripts/runtime_cpu.sh
export SLURM_JOB_ID=5629677
set -x
srun --exact --overlap -N1 -n1 --ntasks-per-node=1 --threads-per-core=1 --cpu-bind=threads -c4 --gpus-per-node=0 python3 analysis-checks/mirror-tail-t14p25/audit_snapshot.py
