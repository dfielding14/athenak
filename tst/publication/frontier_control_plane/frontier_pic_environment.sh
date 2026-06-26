#!/bin/bash
# Frontier HIP/MPI profiles. Promote alternatives only with matched A/B evidence.
PIC_FRONTIER_PROFILE="${PIC_FRONTIER_PROFILE:-frontier_minimum_supported}"
export PIC_FRONTIER_PROFILE

case "$PIC_FRONTIER_PROFILE" in
  frontier_minimum_supported|frontier_xnack1_experimental|frontier_ofi_tuned_experimental)
    ;;
  *)
    printf 'Unsupported PIC_FRONTIER_PROFILE=%s\n' "$PIC_FRONTIER_PROFILE" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac

if ! module --force purge; then
  printf 'Failed to purge inherited Frontier module profile\n' >&2
  return 1 2>/dev/null || exit 1
fi

export MODULEPATH=/opt/cray/pe/modulefiles/Linux:/opt/cray/pe/modulefiles/Core:/opt/cray/pe/lmod/lmod/modulefiles/Core:/opt/cray/pe/lmod/modulefiles/craype-targets/default:/sw/frontier/modulefiles

if ! module use /opt/cray/pe/lmod/modulefiles/core \
    || ! module use /opt/cray/pe/lmod/modulefiles/craype-targets/1.15.0 \
    || ! module use /opt/cray/modulefiles \
    || ! module load cpe/24.11 \
    || ! module load craype-x86-trento \
    || ! module load libfabric/2.3.1 \
    || ! module load craype-network-ofi \
    || ! module load xpmem/1.0.1-1.5_1_gfb6998056825 \
    || ! module load perftools-base/24.11.0 \
    || ! module load cray-pmi/6.1.15 \
    || ! module load cray-dsmml/0.3.0 \
    || ! module load PrgEnv-amd/8.6.0 \
    || ! module load amd/6.2.4 \
    || ! module load rocm/6.2.4 \
    || ! module load craype/2.7.33 \
    || ! module load cray-mpich/8.1.31 \
    || ! module load cray-libsci/24.11.0 \
    || ! module load craype-accel-amd-gfx90a; then
  printf 'Failed to load Frontier PIC module profile\n' >&2
  return 1 2>/dev/null || exit 1
fi

if module is-loaded darshan-runtime && ! module unload darshan-runtime; then
  printf 'Failed to unload inactive darshan-runtime module\n' >&2
  return 1 2>/dev/null || exit 1
fi

# Site Lmod hooks are intentionally absent from the stripped compute-node
# launcher.  Publish the reviewed final value explicitly after all module work.
export MODULEPATH=/sw/frontier/spack-envs/modules/rocmcc/6.2.4/cray-mpich-8.1.31/rocm-6.2.4/rocmcc-6.2.4:/sw/frontier/spack-envs/modules/rocmcc/6.2.4/rocm-6.2.4/rocmcc-6.2.4:/sw/frontier/spack-envs/modules/rocmcc/6.2.4/cray-mpich-8.1.31/rocmcc-6.2.4:/sw/frontier/spack-envs/modules/rocmcc/6.2.4/rocmcc-6.2.4:/opt/cray/pe/lmod/modulefiles/mpi/amd/4.0/ofi/1.0/cray-mpich/8.0:/opt/cray/pe/lmod/modulefiles/comnet/amd/4.0/ofi/1.0:/opt/cray/pe/lmod/modulefiles/compiler/amd/4.0:/opt/cray/pe/lmod/modulefiles/mix_compilers:/opt/cray/pe/lmod/modulefiles/perftools/24.11.0:/opt/cray/pe/lmod/modulefiles/net/ofi/1.0:/opt/cray/pe/lmod/modulefiles/cpu/x86-trento/1.0:/opt/cray/modulefiles:/opt/cray/pe/lmod/modulefiles/craype-targets/1.15.0:/opt/cray/pe/lmod/modulefiles/core:/opt/cray/pe/modulefiles/Linux:/opt/cray/pe/modulefiles/Core:/opt/cray/pe/lmod/lmod/modulefiles/Core:/opt/cray/pe/lmod/modulefiles/craype-targets/default:/sw/frontier/modulefiles

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_EAGER_REGISTER_HOST_MEM=0
export MPICH_GPU_NO_ASYNC_COPY=1
export MPICH_GPU_IPC_ENABLED=0
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export SLURM_EXPORT_ENV=ALL

unset HSA_XNACK
unset MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED
unset MPICH_OFI_NIC_POLICY
unset MPICH_GPU_IPC_CACHE_MAX_SIZE
unset MPICH_MPIIO_HINTS
unset MPICH_OFI_NUM_CQ_ENTRIES
unset FI_MR_CACHE_MONITOR
unset FI_CXI_RX_MATCH_MODE

case "$PIC_FRONTIER_PROFILE" in
  frontier_minimum_supported)
    ;;
  frontier_xnack1_experimental)
    export HSA_XNACK=1
    export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=1
    ;;
  frontier_ofi_tuned_experimental)
    export MPICH_OFI_NIC_POLICY=GPU
    export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
    export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
    export MPICH_OFI_NUM_CQ_ENTRIES=131072
    export FI_MR_CACHE_MONITOR=kdreg2
    export FI_CXI_RX_MATCH_MODE=software
    ;;
  *)
    printf 'Unsupported PIC_FRONTIER_PROFILE=%s\n' "$PIC_FRONTIER_PROFILE" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac

record_pic_environment() {
  local name value
  printf 'PIC_FRONTIER_PROFILE=%s\n' "$PIC_FRONTIER_PROFILE"
  printf 'HSA_XNACK=%s\n' "${HSA_XNACK:-0}"
  for name in MPICH_ENV_DISPLAY MPICH_VERSION_DISPLAY \
      MPICH_GPU_SUPPORT_ENABLED MPICH_GPU_EAGER_REGISTER_HOST_MEM \
      MPICH_GPU_NO_ASYNC_COPY MPICH_GPU_IPC_ENABLED \
      MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED \
      MPICH_OFI_NIC_POLICY MPICH_GPU_IPC_CACHE_MAX_SIZE MPICH_MPIIO_HINTS \
      MPICH_OFI_NUM_CQ_ENTRIES FI_MR_CACHE_MONITOR FI_CXI_RX_MATCH_MODE \
      OMP_NUM_THREADS SLURM_EXPORT_ENV ROCM_PATH LOADEDMODULES _LMFILES_ \
      MODULEPATH; do
    if value="$(/usr/bin/printenv "$name")"; then
      printf '%s=%s\n' "$name" "$value"
    else
      printf '%s=<unset>\n' "$name"
    fi
  done
}
