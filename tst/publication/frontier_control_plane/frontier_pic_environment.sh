#!/bin/bash
# Baseline Frontier HIP/MPI profile. Promote alternatives only with A/B evidence.
module reset
module load PrgEnv-amd
module load amd
module load rocm
module load craype-accel-amd-gfx90a
export MPICH_GPU_SUPPORT_ENABLED=1

record_pic_environment() {
  printf 'HSA_XNACK=%s\n' "${HSA_XNACK:-0}"
  printf 'MPICH_GPU_SUPPORT_ENABLED=%s\n' "$MPICH_GPU_SUPPORT_ENABLED"
  printf 'OMP_NUM_THREADS=%s\n' "${OMP_NUM_THREADS:-unset}"
  printf 'SLURM_EXPORT_ENV=%s\n' "${SLURM_EXPORT_ENV:-unset}"
}
