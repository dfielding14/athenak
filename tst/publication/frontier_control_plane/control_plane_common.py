#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Shared helpers for immutable Frontier PIC submission snapshots."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import tarfile
from typing import Callable, Iterable
import uuid


AUTHORIZED_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
AUTHORIZED_PROJECT_HOME_ROOT = Path("/ccs/proj/ast207/proj-shared/PIC")
AUTHORIZED_CLEAN_CANDIDATE_SOURCE_ROOT = Path("/ccs/home/dfielding/athenak-pic")
AUTHORIZED_ACCOUNT = "AST207"
BUILD_PROVENANCE_FILENAMES = {
    "configure_log": "configure.log",
    "build_log": "build.log",
    "cmake_cache": "CMakeCache.txt",
    "module_list": "modules.txt",
    "toolchain": "toolchain.txt",
    "build_invocations": "build-invocations.json",
    "git_status_preconfigure": "git_status.preconfigure.txt",
    "git_status": "git_status.txt",
    "submodule_status": "submodule_status.txt",
    "environment_allowlist": "environment.allowlist.txt",
    "build_environment": "build-environment.json",
}
AUTHORIZED_PARTITION = "batch"
AUTHORIZED_NODE_HOUR_CAP = 10000.0
AUTHORIZED_LEDGER_MIRROR_TRANSPORT = "filesystem_copy"
AUTHORIZED_PROJECT_HOME_USAGE = [
    "small_append_only_ledger_and_control_plane_mirror",
]
AUTHORIZED_PROJECT_HOME_RETENTION_ROLE = "operational_ledger_mirror_only"
AUTHORIZED_ORION_BULK_EVIDENCE_USAGE = [
    "simulation_outputs",
    "immutable_bulk_signoff_bundles",
    "private_reference_artifact_staging",
    "restore_drill_evidence",
]
AUTHORIZED_ORION_RETENTION_ROLE = (
    "user_selected_sole_bulk_evidence_root_with_documented_durability_risk"
)
AUTHORIZED_LONG_TERM_STORAGE_STATUS = (
    "user_selected_orion_only_with_documented_durability_risk"
)
AUTHORIZED_LONG_TERM_STORAGE_RISK = (
    "Orion-only retention is user-directed and does not provide an institutional "
    "or approved off-site durable archive."
)
AUTHORIZED_LONG_TERM_STORAGE_BLOCKS = [
    "terminal_durable_retention_signoff_pending_external_review",
]
AUTHORIZED_OLCF_SIDE_STORAGE_STATUS = "passed_user_authorized_orion_only_storage"
AUTHORIZED_HISTORICAL_PROJECT_HOME_BULK_ARTIFACTS = (
    "chronology_only_superseded_by_orion_policy_copies_do_not_add_new_bulk_artifacts"
)
AUTHORIZED_LEDGER_GENESIS_AUTHORIZATION = (
    "user_removed_kronos_dependency_and_selected_orion_only_bulk_evidence_root"
)
AUTHORIZED_STORAGE_PREFLIGHT_METHOD = "local_create_write_sync_remove_probe"
PREPARED_ARTIFACT_INVENTORY_PATH = (
    "tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json"
)


def scheduler_account_matches_authorized(value: object) -> bool:
    """Accept the configured account or Slurm's canonical lowercase spelling."""
    return isinstance(value, str) and value in {
        AUTHORIZED_ACCOUNT,
        AUTHORIZED_ACCOUNT.lower(),
    }


def _is_lowercase_sha256(value: object) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


TRUSTED_GIT = "/usr/bin/git"
TRUSTED_GIT_OPTIONS = [
    "-c",
    "core.fsmonitor=false",
    "-c",
    "core.hooksPath=/dev/null",
]
TRUSTED_PYTHON = "/opt/cray/pe/python/3.11.7/bin/python3"
TRUSTED_CMAKE = "/usr/bin/cmake"
TRUSTED_CXX_COMPILER = "/opt/cray/pe/craype/2.7.33/bin/CC"
TRUSTED_ROCM_PATH = "/opt/rocm-6.2.4"
PRODUCTION_BUILD_PROFILE = "hip-mpi-release-paper-pic"
PRODUCTION_TOOLCHAIN_DESCRIPTION = (
    "Frontier PrgEnv-amd/8.6.0 CC wrapper /opt/cray/pe/craype/2.7.33/bin/CC, "
    "ROCm/6.2.4 HIP /opt/rocm-6.2.4, GPU-aware Cray MPICH/8.1.31, gfx90a"
)
PRODUCTION_ENVIRONMENT_ALLOWLIST = [
    ("PIC_FRONTIER_PROFILE", "frontier_minimum_supported"),
    ("HSA_XNACK", "<unset>"),
    ("MPICH_ENV_DISPLAY", "1"),
    ("MPICH_VERSION_DISPLAY", "1"),
    ("MPICH_GPU_SUPPORT_ENABLED", "1"),
    ("MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED", "<unset>"),
    ("MPICH_OFI_NIC_POLICY", "<unset>"),
    ("MPICH_GPU_IPC_CACHE_MAX_SIZE", "<unset>"),
    ("MPICH_MPIIO_HINTS", "<unset>"),
    ("MPICH_OFI_NUM_CQ_ENTRIES", "<unset>"),
    ("FI_MR_CACHE_MONITOR", "<unset>"),
    ("FI_CXI_RX_MATCH_MODE", "<unset>"),
    ("OMP_NUM_THREADS", "<unset>"),
    ("SLURM_EXPORT_ENV", "ALL"),
    ("ROCM_PATH", TRUSTED_ROCM_PATH),
]
PRODUCTION_RUNTIME_LOADED_MODULES = (
    "cpe/24.11",
    "craype-x86-trento",
    "libfabric/2.3.1",
    "craype-network-ofi",
    "xpmem/1.0.1-1.5_1_gfb6998056825",
    "perftools-base/24.11.0",
    "cray-pmi/6.1.15",
    "cray-dsmml/0.3.0",
    "PrgEnv-amd/8.6.0",
    "amd/6.2.4",
    "rocm/6.2.4",
    "craype/2.7.33",
    "cray-mpich/8.1.31",
    "cray-libsci/24.11.0",
    "craype-accel-amd-gfx90a",
)
PRODUCTION_RUNTIME_MODULEFILES = (
    "/opt/cray/pe/lmod/modulefiles/core/cpe/24.11.lua",
    "/opt/cray/pe/lmod/modulefiles/craype-targets/1.15.0/craype-x86-trento.lua",
    "/opt/cray/modulefiles/libfabric/2.3.1",
    "/opt/cray/pe/lmod/modulefiles/craype-targets/1.15.0/craype-network-ofi.lua",
    "/opt/cray/modulefiles/xpmem/1.0.1-1.5_1_gfb6998056825",
    "/opt/cray/pe/lmod/modulefiles/core/perftools-base/24.11.0.lua",
    "/opt/cray/pe/lmod/modulefiles/core/cray-pmi/6.1.15.lua",
    "/opt/cray/pe/lmod/modulefiles/core/cray-dsmml/0.3.0.lua",
    "/opt/cray/pe/lmod/modulefiles/core/PrgEnv-amd/8.6.0.lua",
    "/opt/cray/pe/lmod/modulefiles/core/amd/6.2.4.lua",
    "/opt/cray/pe/lmod/modulefiles/core/rocm/6.2.4.lua",
    "/opt/cray/pe/lmod/modulefiles/core/craype/2.7.33.lua",
    "/opt/cray/pe/lmod/modulefiles/comnet/amd/4.0/ofi/1.0/cray-mpich/8.1.31.lua",
    "/opt/cray/pe/lmod/modulefiles/compiler/amd/4.0/cray-libsci/24.11.0.lua",
    "/opt/cray/pe/lmod/modulefiles/craype-targets/1.15.0/craype-accel-amd-gfx90a.lua",
)
PRODUCTION_RUNTIME_MODULEPATH = (
    "/sw/frontier/spack-envs/modules/rocmcc/6.2.4/cray-mpich-8.1.31/rocm-6.2.4/rocmcc-6.2.4:"
    "/sw/frontier/spack-envs/modules/rocmcc/6.2.4/rocm-6.2.4/rocmcc-6.2.4:"
    "/sw/frontier/spack-envs/modules/rocmcc/6.2.4/cray-mpich-8.1.31/rocmcc-6.2.4:"
    "/sw/frontier/spack-envs/modules/rocmcc/6.2.4/rocmcc-6.2.4:"
    "/opt/cray/pe/lmod/modulefiles/mpi/amd/4.0/ofi/1.0/cray-mpich/8.0:"
    "/opt/cray/pe/lmod/modulefiles/comnet/amd/4.0/ofi/1.0:"
    "/opt/cray/pe/lmod/modulefiles/compiler/amd/4.0:"
    "/opt/cray/pe/lmod/modulefiles/mix_compilers:"
    "/opt/cray/pe/lmod/modulefiles/perftools/24.11.0:"
    "/opt/cray/pe/lmod/modulefiles/net/ofi/1.0:"
    "/opt/cray/pe/lmod/modulefiles/cpu/x86-trento/1.0:"
    "/opt/cray/modulefiles:"
    "/opt/cray/pe/lmod/modulefiles/craype-targets/1.15.0:"
    "/opt/cray/pe/lmod/modulefiles/core:"
    "/opt/cray/pe/modulefiles/Linux:"
    "/opt/cray/pe/modulefiles/Core:"
    "/opt/cray/pe/lmod/lmod/modulefiles/Core:"
    "/opt/cray/pe/lmod/modulefiles/craype-targets/default:"
    "/sw/frontier/modulefiles"
)
PRODUCTION_REQUIRED_MODULES = set(PRODUCTION_RUNTIME_LOADED_MODULES)
PRODUCTION_BUILD_ENVIRONMENT = {
    "CMAKE_PREFIX_PATH": f"{TRUSTED_ROCM_PATH}/lib/cmake/hip:{TRUSTED_ROCM_PATH}",
    "CRAYPAT_LD_LIBRARY_PATH": "/opt/cray/pe/perftools/24.11.0/lib64",
    "CRAYPAT_OPTS_EXECUTABLE": "libexec64/opts",
    "CRAYPAT_ROOT": "/opt/cray/pe/perftools/24.11.0",
    "CRAYPE_DIR": "/opt/cray/pe/craype/2.7.33",
    "CRAYPE_LINK_TYPE": "dynamic",
    "CRAYPE_NETWORK_TARGET": "ofi",
    "CRAYPE_VERSION": "2.7.33",
    "CRAY_ACCEL_TARGET": "amd_gfx90a",
    "CRAY_ACCEL_VENDOR": "amd",
    "CRAY_AMD_COMPILER_PREFIX": TRUSTED_ROCM_PATH,
    "CRAY_AMD_COMPILER_VERSION": "6.2.4",
    "CRAY_CPU_TARGET": "x86-trento",
    "CRAY_DSMML_BASEDIR": "/opt/cray/pe/dsmml/0.3.0",
    "CRAY_DSMML_DIR": "/opt/cray/pe/dsmml/0.3.0/dsmml",
    "CRAY_DSMML_PREFIX": "/opt/cray/pe/dsmml/0.3.0/dsmml",
    "CRAY_DSMML_ROOTDIR": "/opt/cray/pe/dsmml/0.3.0",
    "CRAY_DSMML_VER": "0.3.0",
    "CRAY_DSMML_VERSION": "0.3.0",
    "CRAY_LD_LIBRARY_PATH": (
        "/opt/cray/pe/libsci/24.11.0/AMD/6.0/x86_64/lib:"
        "/opt/cray/pe/pmi/6.1.15/lib:"
        "/opt/cray/pe/mpich/8.1.31/ofi/amd/6.0/lib:"
        "/opt/cray/pe/mpich/8.1.31/gtl/lib:"
        "/opt/cray/pe/perftools/24.11.0/lib64:"
        "/opt/cray/pe/dsmml/0.3.0/dsmml/lib"
    ),
    "CRAY_LIBSCI_BASE_DIR": "/opt/cray/pe/libsci/24.11.0",
    "CRAY_LIBSCI_PREFIX": "/opt/cray/pe/libsci/24.11.0/AMD/6.0/x86_64",
    "CRAY_LIBSCI_PREFIX_DIR": "/opt/cray/pe/libsci/24.11.0/AMD/6.0/x86_64",
    "CRAY_LIBSCI_VERSION": "24.11.0",
    "CRAY_LMOD_COMPILER": "amd/4.0",
    "CRAY_LMOD_CPU": "x86-trento/1.0",
    "CRAY_LMOD_MPI": "cray-mpich/8.0",
    "CRAY_LMOD_NET": "ofi/1.0",
    "CRAY_MPICH_BASEDIR": "/opt/cray/pe/mpich/8.1.31/ofi",
    "CRAY_MPICH_DIR": "/opt/cray/pe/mpich/8.1.31/ofi/amd/6.0",
    "CRAY_MPICH_PREFIX": "/opt/cray/pe/mpich/8.1.31/ofi/amd/6.0",
    "CRAY_MPICH_ROOTDIR": "/opt/cray/pe/mpich/8.1.31",
    "CRAY_MPICH_VER": "8.1.31",
    "CRAY_MPICH_VERSION": "8.1.31",
    "CRAY_PERFTOOLS_PREFIX": "/opt/cray/pe/perftools/24.11.0",
    "CRAY_PERFTOOLS_VERSION": "24.11.0",
    "CRAY_PMI_INCLUDE_OPTS": "-I/opt/cray/pe/pmi/6.1.15/include",
    "CRAY_PMI_POST_LINK_OPTS": "-L/opt/cray/pe/pmi/6.1.15/lib",
    "CRAY_PMI_PREFIX": "/opt/cray/pe/pmi/6.1.15",
    "CRAY_PMI_VERSION": "6.1.15",
    "CRAY_ROCM_DIR": TRUSTED_ROCM_PATH,
    "CRAY_ROCM_INCLUDE_OPTS": (
        "-I/opt/rocm-6.2.4/include -I/opt/rocm-6.2.4/include/rocprofiler "
        "-I/opt/rocm-6.2.4/include/roctracer -I/opt/rocm-6.2.4/include/hip "
        "-D__HIP_PLATFORM_AMD__"
    ),
    "CRAY_ROCM_POST_LINK_OPTS": (
        " -L/opt/rocm-6.2.4/lib -L/opt/rocm-6.2.4/lib/rocprofiler "
        "-L/opt/rocm-6.2.4/lib/roctracer -lamdhip64"
    ),
    "CRAY_ROCM_PREFIX": TRUSTED_ROCM_PATH,
    "CRAY_ROCM_VERSION": "6.2.4",
    "CRAY_TCMALLOC_MEMFS_FORCE": "1",
    "CRAY_XPMEM_INCLUDE_OPTS": "-I/opt/xpmem/include",
    "CRAY_XPMEM_POST_LINK_OPTS": "-L/opt/xpmem/lib64",
    "FI_CXI_ATS": "0",
    "HIP_LIB_PATH": f"{TRUSTED_ROCM_PATH}/lib",
    "HOME": "/",
    "LANG": "C",
    "LC_ALL": "C",
    "MPICH_DIR": "/opt/cray/pe/mpich/8.1.31/ofi/amd/6.0",
    "PATH": f"{TRUSTED_ROCM_PATH}/bin:/opt/cray/pe/craype/2.7.33/bin:/usr/bin:/bin",
    "PE_AMD_FIXED_PKGCONFIG_PATH": "/opt/cray/pe/mpich/8.1.31/ofi/amd/6.0/lib/pkgconfig",
    "PE_DSMML_MODULE_NAME": "cray-dsmml",
    "PE_DSMML_PKGCONFIG_LIBS": "dsmml",
    "PE_ENV": "AMD",
    "PE_FORTRAN_PKGCONFIG_LIBS": "mpichf90",
    "PE_LIBSCI_GENCOMPILERS_AMD_x86_64": "6.0",
    "PE_LIBSCI_GENCOMPS_AMD_x86_64": "60",
    "PE_LIBSCI_MODULE_NAME": "cray-libsci/24.11.0",
    "PE_LIBSCI_OMP_REQUIRES": " ",
    "PE_LIBSCI_OMP_REQUIRES_openmp": "_mp",
    "PE_LIBSCI_PKGCONFIG_LIBS": "libsci_mpi:libsci",
    "PE_LIBSCI_PKGCONFIG_VARIABLES": (
        "PE_LIBSCI_OMP_REQUIRES_@openmp@:PE_SCI_EXT_LIBPATH:PE_SCI_EXT_LIBNAME"
    ),
    "PE_LIBSCI_REQUIRED_PRODUCTS": "PE_MPICH",
    "PE_LIBSCI_VOLATILE_PKGCONFIG_PATH": (
        "/opt/cray/pe/libsci/24.11.0/@PRGENV@/@PE_LIBSCI_GENCOMPS@/"
        "@PE_LIBSCI_TARGET@/lib/pkgconfig"
    ),
    "PE_LIBSCI_VOLATILE_PRGENV": "AMD",
    "PE_MPICH_FIXED_PRGENV": "AMD",
    "PE_MPICH_FORTRAN_PKGCONFIG_LIBS": "mpichf90",
    "PE_MPICH_GENCOMPILERS_AMD": "6.0",
    "PE_MPICH_GTL_DIR_amd_gfx906": "-L/opt/cray/pe/mpich/8.1.31/gtl/lib",
    "PE_MPICH_GTL_DIR_amd_gfx908": "-L/opt/cray/pe/mpich/8.1.31/gtl/lib",
    "PE_MPICH_GTL_DIR_amd_gfx90a": "-L/opt/cray/pe/mpich/8.1.31/gtl/lib",
    "PE_MPICH_GTL_DIR_amd_gfx940": "-L/opt/cray/pe/mpich/8.1.31/gtl/lib",
    "PE_MPICH_GTL_DIR_amd_gfx942": "-L/opt/cray/pe/mpich/8.1.31/gtl/lib",
    "PE_MPICH_GTL_DIR_nvidia70": "-L/opt/cray/pe/mpich/8.1.31/gtl/lib",
    "PE_MPICH_GTL_DIR_nvidia80": "-L/opt/cray/pe/mpich/8.1.31/gtl/lib",
    "PE_MPICH_GTL_DIR_nvidia90": "-L/opt/cray/pe/mpich/8.1.31/gtl/lib",
    "PE_MPICH_GTL_DIR_ponteVecchio": "-L/opt/cray/pe/mpich/8.1.31/gtl/lib",
    "PE_MPICH_GTL_LIBS_amd_gfx906": "-lmpi_gtl_hsa",
    "PE_MPICH_GTL_LIBS_amd_gfx908": "-lmpi_gtl_hsa",
    "PE_MPICH_GTL_LIBS_amd_gfx90a": "-lmpi_gtl_hsa",
    "PE_MPICH_GTL_LIBS_amd_gfx940": "-lmpi_gtl_hsa",
    "PE_MPICH_GTL_LIBS_amd_gfx942": "-lmpi_gtl_hsa",
    "PE_MPICH_GTL_LIBS_nvidia70": "-lmpi_gtl_cuda",
    "PE_MPICH_GTL_LIBS_nvidia80": "-lmpi_gtl_cuda",
    "PE_MPICH_GTL_LIBS_nvidia90": "-lmpi_gtl_cuda",
    "PE_MPICH_GTL_LIBS_ponteVecchio": "-lmpi_gtl_ze",
    "PE_MPICH_MODULE_NAME": "cray-mpich",
    "PE_MPICH_PKGCONFIG_LIBS": "mpich",
    "PE_MPICH_PKGCONFIG_VARIABLES": (
        "PE_MPICH_GTL_DIR_@accelerator@:PE_MPICH_GTL_LIBS_@accelerator@"
    ),
    "PE_PERFTOOLS_MPICH_LIBDIR": "/opt/cray/pe/mpich/8.1.31/ofi/amd/6.0/lib",
    "PE_PKGCONFIG_LIBS": "libsci_mpi:libsci:mpich:rocm-6.2.4:dsmml",
    "PE_PKGCONFIG_PRODUCTS": "PE_LIBSCI:PE_PMI:PE_MPICH:PE_DSMML:PE_XPMEM",
    "PE_PMI_PKGCONFIG_LIBS": "cray-pmi",
    "PE_PRODUCT_LIST": "CRAY_PMI:CRAYPE:CRAYPE_X86_TRENTO:PERFTOOLS:CRAYPAT:CRAY_ROCM:CRAY_ACCEL",
    "PE_XPMEM_PKGCONFIG_LIBS": "xpmem",
    "PKG_CONFIG_PATH": (
        "/opt/cray/pe/pmi/6.1.15/lib/pkgconfig:"
        "/opt/cray/pe/craype/2.7.33/pkg-config:/usr/lib64/pkgconfig:"
        "/opt/cray/pe/dsmml/0.3.0/dsmml/lib/pkgconfig:/opt/cray/libfabric/2.3.1/lib64/pkgconfig"
    ),
    "ROCM_PATH": TRUSTED_ROCM_PATH,
}
TRUSTED_SQUEUE = "/usr/bin/squeue"
TRUSTED_SCONTROL = "/usr/bin/scontrol"
TRUSTED_SACCT = "/usr/bin/sacct"
TRUSTED_SBATCH = "/usr/bin/sbatch"
TRUSTED_SCANCEL = "/usr/bin/scancel"
AUTHORIZED_SLURM_CLUSTER = "frontier"
REGISTERED_SCIENCE_SCOPE = "registered_science"
FRONTIER_ADMISSION_SMOKE_SCOPE = "frontier_admission_smoke"
SUBMISSION_SCOPES = {
    REGISTERED_SCIENCE_SCOPE,
    FRONTIER_ADMISSION_SMOKE_SCOPE,
}
PENDING_CLEAN_CANDIDATE_FREEZE = "pending_clean_candidate_freeze"
AUTHORIZED_CLEAN_CANDIDATE_FREEZE = "authorized"
AUTHORIZED_ADMISSION_SMOKE_STATUS = "authorized_f0_parser_contract_only"
PENDING_ADMISSION_SMOKE_STATUS = "pending_exact_executable_binding"
CLOSED_ADMISSION_SMOKE_STATUS = "closed_after_pass"
AUTHORIZED_REGISTERED_SCIENCE_SLICE_STATUS = "authorized"
TRUSTED_LAUNCH_EXECUTOR = "trusted_trampoline_athena_argv_v1"
CANONICAL_POLICY_RELATIVE = Path("policy/storage_policy.json")
ACTIVE_PROMOTION_RELATIVE = Path("policy/active_promotion.json")
SITE_POLICY_MAX_AGE_SECONDS = 24 * 60 * 60
CONTROL_PLANE_FILES = [
    "clean_candidate.schema.json",
    "control_plane.schema.json",
    "control_plane_common.py",
    "create_clean_candidate_freeze.py",
    "create_pre_submit_manifest.py",
    "frontier_pic_environment.sh",
    "initialize_frontier_ledger.py",
    "launch_trampoline.py",
    "launch_with_frontier_profile.sh",
    "ledger.py",
    "promote_active_policy.py",
    "reconcile_frontier_job.py",
    "reconcile_manual_frontier_allocations.py",
    "run_installed_control_plane_job.sh",
    "run_control_plane.py",
    "submit_frontier_job.sh",
    "terminal_recovery_handoff.py",
    "validate_and_reserve_frontier_job.py",
    "verify_compute_node_snapshot.py",
    "write_orion_build_profile.py",
]
PLACEHOLDER_PATTERN = re.compile(rb"REPLACE_[A-Z0-9_]+")
SENSITIVE_PATTERN = re.compile(
    rb"(TOKEN|PASSWORD|SECRET|PRIVATE[_-]?KEY|BEGIN [A-Z ]*PRIVATE KEY)",
    re.IGNORECASE,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def trusted_git_environment() -> dict[str, str]:
    """Return a Git environment that does not load caller-controlled config."""
    return {
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1",
        "HOME": "/",
        "LANG": "C",
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
    }


def trusted_git_command(*arguments: str) -> list[str]:
    """Return a Git argv that disables repository-local command execution hooks."""
    return [TRUSTED_GIT, *TRUSTED_GIT_OPTIONS, *arguments]


def trusted_slurm_environment() -> dict[str, str]:
    """Return a scheduler environment that cannot inherit caller routing controls."""
    return {
        "HOME": "/",
        "LANG": "C",
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
        "SLURM_CLUSTERS": AUTHORIZED_SLURM_CLUSTER,
    }


def production_environment_allowlist_bytes() -> bytes:
    return "".join(f"{name}={value}\n" for name, value in PRODUCTION_ENVIRONMENT_ALLOWLIST).encode(
        "utf-8"
    )


def production_module_list_bytes() -> bytes:
    """Return the exact reviewed module selection recorded for production builds."""
    return "".join(
        f"{module}\t{modulefile}\n"
        for module, modulefile in zip(
            PRODUCTION_RUNTIME_LOADED_MODULES, PRODUCTION_RUNTIME_MODULEFILES
        )
    ).encode("utf-8")


def require_production_module_environment(environment: dict[str, str]) -> None:
    """Require the exact clean Frontier module stack and modulefile provenance."""
    if tuple(environment.get("LOADEDMODULES", "").split(":")) != (
        PRODUCTION_RUNTIME_LOADED_MODULES
    ):
        raise ValueError("Loaded Frontier modules differ from the reviewed exact stack")
    if tuple(environment.get("_LMFILES_", "").split(":")) != (
        PRODUCTION_RUNTIME_MODULEFILES
    ):
        raise ValueError("Loaded Frontier modulefiles differ from the reviewed exact stack")
    if environment.get("MODULEPATH") != PRODUCTION_RUNTIME_MODULEPATH:
        raise ValueError("Frontier MODULEPATH differs from the reviewed exact value")


def measured_production_module_list_bytes(
    environment: dict[str, str] | None = None,
) -> bytes:
    """Measure and serialize the reviewed production module stack."""
    require_production_module_environment(dict(os.environ if environment is None else environment))
    return production_module_list_bytes()


def require_production_build_environment(environment: dict[str, object]) -> None:
    """Reject caller-controlled process state from the reviewed build subprocesses."""
    if any(not isinstance(key, str) or not isinstance(value, str) for key, value in environment.items()):
        raise ValueError("Production build environment must contain only text keys and values")
    if environment != PRODUCTION_BUILD_ENVIRONMENT:
        raise ValueError("Production build environment differs from the reviewed exact values")


def production_build_invocations(
    *,
    authorized_pic_root: Path,
    git_commit: str,
    profile_id: str,
) -> dict[str, list[str]]:
    """Return the one reviewed production build argv pair."""
    if re.fullmatch(r"[0-9a-f]{40}", git_commit) is None:
        raise ValueError("Build-profile Git commit must be a full lowercase hexadecimal commit")
    if profile_id != PRODUCTION_BUILD_PROFILE:
        raise ValueError(f"Unsupported installed Frontier build profile: {profile_id}")
    root = Path(os.path.abspath(authorized_pic_root))
    cmake_dir = root / "build" / git_commit[:12] / profile_id / "cmake"
    fresh_source = root / "build" / git_commit[:12] / profile_id / "source"
    return {
        "configure": [
            TRUSTED_CMAKE,
            "-S",
            str(fresh_source),
            "-B",
            str(cmake_dir),
            "-DCMAKE_BUILD_TYPE=Release",
            "-DAthena_ENABLE_MPI=ON",
            "-DKokkos_ENABLE_HIP=ON",
            "-DKokkos_ARCH_ZEN3=ON",
            "-DKokkos_ARCH_AMD_GFX90A=ON",
            f"-DCMAKE_CXX_COMPILER={TRUSTED_CXX_COMPILER}",
            f"-DCMAKE_CXX_FLAGS=-I{TRUSTED_ROCM_PATH}/include",
            f"-DCMAKE_EXE_LINKER_FLAGS=-L{TRUSTED_ROCM_PATH}/lib -lamdhip64",
            "-DPROBLEM=built_in_pgens",
        ],
        "build": [TRUSTED_CMAKE, "--build", str(cmake_dir), "--parallel", "32"],
    }


def require_production_build_provenance(
    *,
    authorized_pic_root: Path,
    git_commit: str,
    profile_id: str,
    toolchain: str,
    invocations: dict[str, object],
    module_list: bytes,
    environment_allowlist: bytes,
    build_environment: bytes,
) -> None:
    """Require exact production semantics for artifacts rooted in the Orion PIC tree."""
    if Path(os.path.abspath(authorized_pic_root)) != Path(os.path.abspath(AUTHORIZED_PIC_ROOT)):
        return
    if profile_id != PRODUCTION_BUILD_PROFILE:
        raise ValueError(f"Unsupported installed Frontier build profile: {profile_id}")
    if toolchain != PRODUCTION_TOOLCHAIN_DESCRIPTION:
        raise ValueError("Production build toolchain description differs from the reviewed value")
    if invocations != production_build_invocations(
        authorized_pic_root=authorized_pic_root,
        git_commit=git_commit,
        profile_id=profile_id,
    ):
        raise ValueError("Production build invocations differ from the reviewed argv")
    if environment_allowlist != production_environment_allowlist_bytes():
        raise ValueError("Production build environment differs from the reviewed allowlist")
    try:
        environment = read_json_bytes(build_environment, label="production build environment")
    except ValueError as error:
        raise ValueError("Production build environment is not valid JSON") from error
    require_production_build_environment(environment)
    if module_list != production_module_list_bytes():
        raise ValueError("Production module list differs from the reviewed exact selection")


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"Non-finite JSON number is not allowed: {value}")


def _reject_duplicate_json_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON object key is not allowed: {key}")
        value[key] = item
    return value


def read_json_bytes(data: bytes, *, label: str) -> dict[str, object]:
    try:
        value = json.loads(
            data.decode("utf-8"),
            parse_constant=_reject_json_constant,
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"Expected valid UTF-8 JSON in {label}") from error
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {label}")
    return value


def read_stable_regular_file(path: Path, *, require_read_only_mode: bool = False) -> bytes:
    """Read one non-symlink regular file once so later checks use identical bytes."""
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"Artifact is not a regular file: {path}")
        if require_read_only_mode and metadata.st_mode & 0o222:
            raise ValueError(f"Artifact is not read-only: {path}")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            return stream.read()
    finally:
        os.close(descriptor)


def read_stable_regular_file_below(
    path: Path,
    root: Path,
    *,
    require_read_only_mode: bool = False,
) -> bytes:
    """Read a file through anchored non-symlink directories below one trusted root."""
    lexical_root = Path(os.path.abspath(root))
    lexical_root.resolve(strict=True)
    lexical_path = Path(os.path.abspath(path))
    try:
        relative = lexical_path.relative_to(lexical_root)
    except ValueError as error:
        raise ValueError(f"Path is outside authorized lexical root: {lexical_path}") from error
    if not relative.parts:
        raise ValueError(f"Expected a file below authorized root: {lexical_path}")
    directory_flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    file_flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    root_descriptor = os.open(
        lexical_root,
        os.O_RDONLY | os.O_DIRECTORY,
    )
    directory_descriptor = root_descriptor
    descriptor: int | None = None
    try:
        for part in relative.parts[:-1]:
            next_descriptor = os.open(part, directory_flags, dir_fd=directory_descriptor)
            if directory_descriptor != root_descriptor:
                os.close(directory_descriptor)
            directory_descriptor = next_descriptor
        descriptor = os.open(relative.parts[-1], file_flags, dir_fd=directory_descriptor)
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"Artifact is not a regular file: {path}")
        if require_read_only_mode and metadata.st_mode & 0o222:
            raise ValueError(f"Artifact is not read-only: {path}")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            return stream.read()
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if directory_descriptor != root_descriptor:
            os.close(directory_descriptor)
        os.close(root_descriptor)


def git_archive_commit_from_bytes(data: bytes) -> str:
    return subprocess.check_output(
        trusted_git_command("get-tar-commit-id"),
        input=data,
        env=trusted_git_environment(),
    ).decode("ascii").strip()


def git_commit_tree_from_bytes(data: bytes, *, expected_commit: str) -> str:
    """Verify one raw Git commit object and return its referenced tree."""
    if not re.fullmatch(r"[0-9a-f]{40}", expected_commit):
        raise ValueError("Git commit object ID is malformed")
    try:
        actual_commit = subprocess.check_output(
            trusted_git_command("hash-object", "--stdin", "-t", "commit"),
            input=data,
            stderr=subprocess.PIPE,
            env=trusted_git_environment(),
        ).decode("ascii").strip()
    except subprocess.CalledProcessError as error:
        raise ValueError("Raw Git commit object is malformed") from error
    if actual_commit != expected_commit:
        raise ValueError("Raw Git commit object does not match its declared ID")
    headers, separator, _ = data.partition(b"\n\n")
    if not separator:
        raise ValueError("Raw Git commit object has no message separator")
    names: list[bytes] = []
    values: list[bytes] = []
    for line in headers.splitlines():
        if line.startswith(b" "):
            if not names:
                raise ValueError("Raw Git commit object starts with a continuation header")
            continue
        match = re.fullmatch(rb"([a-z][a-z0-9-]*) (.+)", line)
        if match is None:
            raise ValueError("Raw Git commit object has a malformed header")
        names.append(match.group(1))
        values.append(match.group(2))
    if not names:
        raise ValueError("Raw Git commit object has no headers")
    match = re.fullmatch(rb"([0-9a-f]{40})", values[0]) if names[0] == b"tree" else None
    if match is None:
        raise ValueError("Raw Git commit object has no canonical tree header")
    index = 1
    while index < len(names) and names[index] == b"parent":
        if re.fullmatch(rb"[0-9a-f]{40}", values[index]) is None:
            raise ValueError("Raw Git commit object has a malformed parent header")
        index += 1
    if names[index:index + 2] != [b"author", b"committer"]:
        raise ValueError("Raw Git commit object has noncanonical identity headers")
    reserved = {b"tree", b"parent", b"author", b"committer"}
    if any(name in reserved for name in names[index + 2:]):
        raise ValueError("Raw Git commit object repeats a reserved header")
    return match.group(1).decode("ascii")


def canonical_relative_posix_path(value: object, *, field: str) -> PurePosixPath:
    if not isinstance(value, str):
        raise ValueError(f"{field} must be a string")
    path = PurePosixPath(value)
    if (
        not value
        or path.is_absolute()
        or not path.parts
        or value != path.as_posix()
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise ValueError(f"{field} is not a canonical relative path: {value!r}")
    return path


def _source_archive_regular_files(data: bytes) -> dict[str, bytes]:
    """Return exact regular-file members from one validated source archive."""
    records: dict[str, bytes] = {}
    names: set[str] = set()
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:*") as stream:
        for member in stream.getmembers():
            path = PurePosixPath(member.name)
            canonical_name = path.as_posix()
            if path.is_absolute() or not path.parts or any(
                part in {"", ".", ".."} for part in path.parts
            ) or member.name != canonical_name:
                raise ValueError(f"Unsafe path in source archive: {member.name!r}")
            if canonical_name in names:
                raise ValueError(f"Duplicate path in source archive: {member.name!r}")
            names.add(canonical_name)
            if not member.isfile():
                continue
            extracted = stream.extractfile(member)
            if extracted is None:
                raise ValueError(f"Cannot read source-archive member: {member.name!r}")
            records[canonical_name] = extracted.read()
    return records


def _prepared_artifact_records(
    value: object, *, field: str, source_files: dict[str, bytes]
) -> list[dict[str, str]]:
    if not isinstance(value, list) or not 1 <= len(value) <= 4096:
        raise ValueError(f"{field} must contain between 1 and 4096 records")
    records: list[dict[str, str]] = []
    for raw in value:
        if not isinstance(raw, dict) or set(raw) != {"path", "sha256"}:
            raise ValueError(f"{field} record is malformed")
        path = canonical_relative_posix_path(
            raw.get("path"), field=f"{field} path"
        ).as_posix()
        digest = str(raw.get("sha256", ""))
        if re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            raise ValueError(f"{field} SHA-256 is malformed: {path}")
        data = source_files.get(path)
        if data is None:
            raise ValueError(
                f"{field} path is not an exact regular member of the source archive: {path}"
            )
        if sha256_bytes(data) != digest:
            raise ValueError(f"{field} archived-byte checksum mismatch: {path}")
        records.append({"path": path, "sha256": digest})
    paths = [record["path"] for record in records]
    if paths != sorted(set(paths)):
        raise ValueError(f"{field} must use unique canonical path order")
    return records


def prepared_artifact_manifest_from_source_archive(
    source_archive: bytes, *, inventory_path: object
) -> dict[str, object]:
    """Derive one prepared-artifact manifest from a committed source inventory."""
    source_files = _source_archive_regular_files(source_archive)
    normalized_inventory_path = canonical_relative_posix_path(
        inventory_path, field="Prepared-artifact inventory path"
    ).as_posix()
    if normalized_inventory_path != PREPARED_ARTIFACT_INVENTORY_PATH:
        raise ValueError(
            "Prepared-artifact inventory path must use the canonical source-relative path: "
            f"{PREPARED_ARTIFACT_INVENTORY_PATH}"
        )
    inventory_bytes = source_files.get(normalized_inventory_path)
    if inventory_bytes is None:
        raise ValueError(
            "Prepared-artifact inventory path is not an exact regular member "
            f"of the source archive: {normalized_inventory_path}"
        )
    inventory = read_json_bytes(
        inventory_bytes, label=f"prepared-artifact inventory {normalized_inventory_path}"
    )
    if set(inventory) != {"schema_version", "paper_decks", "analyzers"}:
        raise ValueError("Prepared-artifact inventory has unexpected fields")
    if type(inventory.get("schema_version")) is not int or inventory.get("schema_version") != 1:
        raise ValueError("Unsupported prepared-artifact inventory schema")
    paper_decks = _prepared_artifact_records(
        inventory.get("paper_decks"),
        field="Prepared paper-deck inventory",
        source_files=source_files,
    )
    analyzers = _prepared_artifact_records(
        inventory.get("analyzers"),
        field="Prepared analyzer inventory",
        source_files=source_files,
    )
    expected_paper_decks = sorted(
        path
        for path in source_files
        if path.startswith("inputs/tests/pic")
        and path.endswith(".athinput")
        and "/" not in path[len("inputs/tests/") :]
    )
    expected_analyzers = sorted(
        path
        for path in source_files
        if path.startswith("tst/publication/analyze_")
        and path.endswith(".py")
        and "/" not in path[len("tst/publication/") :]
    )
    if [record["path"] for record in paper_decks] != expected_paper_decks:
        raise ValueError(
            "Prepared paper-deck inventory must exactly cover archived inputs/tests/pic*.athinput"
        )
    if [record["path"] for record in analyzers] != expected_analyzers:
        raise ValueError(
            "Prepared analyzer inventory must exactly cover archived tst/publication/analyze_*.py"
        )
    paths = [record["path"] for record in paper_decks + analyzers]
    if normalized_inventory_path in paths or len(paths) != len(set(paths)):
        raise ValueError(
            "Prepared-artifact paths must be distinct from each other and the inventory"
        )
    return {
        "inventory_path": normalized_inventory_path,
        "inventory_sha256": sha256_bytes(inventory_bytes),
        "paper_decks": paper_decks,
        "analyzers": analyzers,
    }


def validate_prepared_artifact_closure(
    value: object, *, source_archive: bytes
) -> dict[str, object]:
    """Revalidate candidate-level prepared artifacts against archived source bytes."""
    if not isinstance(value, dict) or set(value) != {
        "inventory_path",
        "inventory_sha256",
        "paper_decks",
        "analyzers",
    }:
        raise ValueError("Clean-candidate prepared-artifact attestation has unexpected fields")
    inventory_sha256 = str(value.get("inventory_sha256", ""))
    if re.fullmatch(r"[0-9a-f]{64}", inventory_sha256) is None:
        raise ValueError("Prepared-artifact inventory SHA-256 is malformed")
    expected = prepared_artifact_manifest_from_source_archive(
        source_archive, inventory_path=value.get("inventory_path")
    )
    if value != expected:
        raise ValueError(
            "Clean-candidate prepared-artifact manifest differs from archived source inventory"
        )
    return expected


def _git_object_sha1(kind: str, data: bytes) -> bytes:
    header = f"{kind} {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).digest()


def direct_submodule_gitlinks(
    records: Iterable[dict[str, object]], *, parent_path: str | None = None
) -> dict[str, str]:
    """Return the immediate Git links represented in one archived repository."""
    parsed: list[tuple[PurePosixPath, str]] = []
    seen: set[str] = set()
    for record in records:
        value = str(record["path"])
        path = canonical_relative_posix_path(value, field="Git-link path")
        commit = str(record["git_commit"])
        if (
            not value
            or path.is_absolute()
            or not path.parts
            or not re.fullmatch(r"[0-9a-f]{40}", commit)
        ):
            raise ValueError(f"Unsafe Git-link attestation: {value!r}")
        if value in seen:
            raise ValueError(f"Duplicate Git-link attestation: {value!r}")
        seen.add(value)
        parsed.append((path, commit))
    parent = (
        canonical_relative_posix_path(parent_path, field="Git-link parent path")
        if parent_path is not None
        else None
    )
    result: dict[str, str] = {}
    for path, commit in parsed:
        ancestors = [
            candidate
            for candidate, _ in parsed
            if len(candidate.parts) < len(path.parts)
            and path.parts[: len(candidate.parts)] == candidate.parts
        ]
        nearest = max(ancestors, key=lambda candidate: len(candidate.parts), default=None)
        if nearest != parent:
            continue
        relative = path.relative_to(parent).as_posix() if parent else path.as_posix()
        result[relative] = commit
    return result


def source_bundle_sha256(
    source_archive_sha256: str,
    source_commit_sha256: str,
    submodules: Iterable[dict[str, object]],
) -> str:
    """Bind parent and recursive-submodule archives plus raw Git commit objects."""
    value = {
        "source_archive_sha256": source_archive_sha256,
        "source_commit_sha256": source_commit_sha256,
        "submodules": list(submodules),
    }
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def git_tree_sha1_from_archive_bytes(
    data: bytes,
    *,
    gitlinks: dict[str, str] | None = None,
    reject_symlinks: bool = False,
) -> str:
    """Reconstruct the Git tree object ID represented by a git-archive tarball."""
    root: dict[str, object] = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:*") as stream:
        for member in stream.getmembers():
            path = PurePosixPath(member.name)
            canonical_name = path.as_posix()
            if path.is_absolute() or not path.parts or any(
                part in {"", ".", ".."} for part in path.parts
            ) or member.name != canonical_name:
                raise ValueError(f"Unsafe path in source archive: {member.name!r}")
            node = root
            for part in path.parts[:-1]:
                existing = node.setdefault(part, {})
                if not isinstance(existing, dict):
                    raise ValueError(f"Archive path collides with a file: {member.name!r}")
                node = existing
            name = path.parts[-1]
            if member.isdir():
                existing = node.setdefault(name, {})
                if not isinstance(existing, dict):
                    raise ValueError(f"Archive directory collides with a file: {member.name!r}")
                continue
            if name in node:
                raise ValueError(f"Duplicate path in source archive: {member.name!r}")
            if member.isfile():
                extracted = stream.extractfile(member)
                if extracted is None:
                    raise ValueError(f"Cannot read source-archive member: {member.name!r}")
                payload = extracted.read()
                mode = "100755" if member.mode & 0o111 else "100644"
            elif member.issym():
                if reject_symlinks:
                    raise ValueError(f"Symlink is not allowed in source archive: {member.name!r}")
                payload = member.linkname.encode("utf-8", "surrogateescape")
                mode = "120000"
            else:
                raise ValueError(f"Unsupported source-archive member: {member.name!r}")
            node[name] = (mode, _git_object_sha1("blob", payload))

    for value, commit in (gitlinks or {}).items():
        path = canonical_relative_posix_path(value, field="Git-link path")
        if (
            not re.fullmatch(r"[0-9a-f]{40}", commit)
        ):
            raise ValueError(f"Unsafe Git-link attestation: {value!r}")
        node = root
        for part in path.parts[:-1]:
            existing = node.setdefault(part, {})
            if not isinstance(existing, dict):
                raise ValueError(f"Git-link path collides with a file: {value!r}")
            node = existing
        name = path.parts[-1]
        if name not in node or node[name] != {}:
            raise ValueError(f"Git-link path is not an empty archive directory: {value!r}")
        node[name] = ("160000", bytes.fromhex(commit))

    def tree_sha1(node: dict[str, object]) -> bytes:
        entries = []
        for name, value in node.items():
            encoded_name = name.encode("utf-8", "surrogateescape")
            if isinstance(value, dict):
                mode = "40000"
                digest = tree_sha1(value)
                sort_key = encoded_name + b"/"
            else:
                mode, digest = value
                sort_key = encoded_name
            entry = mode.encode("ascii") + b" " + encoded_name + b"\0" + digest
            entries.append((sort_key, entry))
        payload = b"".join(entry for _, entry in sorted(entries))
        return _git_object_sha1("tree", payload)

    return tree_sha1(root).hex()


def git_tree_sha1_from_archive(
    archive: Path,
    *,
    gitlinks: dict[str, str] | None = None,
    reject_symlinks: bool = False,
) -> str:
    return git_tree_sha1_from_archive_bytes(
        archive.read_bytes(), gitlinks=gitlinks, reject_symlinks=reject_symlinks
    )


def _documented_build_provenance_paths(
    *,
    authorized_pic_root: Path,
    git_commit: str,
    profile_id: str,
) -> dict[str, Path]:
    """Return the reviewed Orion provenance layout for one build profile."""
    if re.fullmatch(r"[0-9a-f]{40}", git_commit) is None:
        raise ValueError("Build-profile Git commit must be a full lowercase hexadecimal commit")
    profile_id = profile_id.strip()
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", profile_id) is None:
        raise ValueError("Build-profile ID must be a safe path component")
    authorized_pic_root = Path(os.path.abspath(authorized_pic_root))
    stem = f"{git_commit[:12]}.{profile_id}"
    artifact_dir = authorized_pic_root / "bin" / git_commit[:12] / profile_id
    log_dir = authorized_pic_root / "logs" / "build"
    return {
        "configure_log": log_dir / f"{stem}.configure.log",
        "build_log": log_dir / f"{stem}.build.log",
        "cmake_cache": artifact_dir / "CMakeCache.txt",
        "module_list": artifact_dir / "modules.txt",
        "toolchain": artifact_dir / "toolchain.txt",
        "build_invocations": artifact_dir / "build-invocations.json",
        "git_status_preconfigure": artifact_dir / "git_status.preconfigure.txt",
        "git_status": artifact_dir / "git_status.txt",
        "submodule_status": artifact_dir / "submodule_status.txt",
        "environment_allowlist": artifact_dir / "environment.allowlist.txt",
        "build_environment": artifact_dir / "build-environment.json",
    }


def validate_clean_candidate_bundle(
    candidate: dict[str, object],
    *,
    source_archive: bytes,
    source_commit: bytes,
    submodule_archives: list[bytes],
    submodule_commits: list[bytes],
    build_profile: bytes,
    build_profile_receipt: bytes,
    build_provenance: dict[str, bytes],
    executable_sha256: str,
    expected_control_plane_version: str | None = None,
    authorized_pic_root: Path | None = None,
    authorized_source_root: Path | None = None,
) -> list[dict[str, str]]:
    """Validate the cryptographic closure and production-authorized build layout."""
    authorized_pic_root = Path(
        os.path.abspath(AUTHORIZED_PIC_ROOT if authorized_pic_root is None else authorized_pic_root)
    )
    authorized_source_root = Path(
        os.path.abspath(
            AUTHORIZED_CLEAN_CANDIDATE_SOURCE_ROOT
            if authorized_source_root is None
            else authorized_source_root
        )
    )
    if set(candidate) != {
        "schema_version",
        "freeze_id",
        "created_utc",
        "prepared_artifacts",
        "source",
        "build",
    }:
        raise ValueError("Clean-candidate manifest has unexpected top-level fields")
    if type(candidate.get("schema_version")) is not int or candidate.get("schema_version") != 4:
        raise ValueError("Unsupported clean-candidate manifest schema")
    source = candidate.get("source")
    build = candidate.get("build")
    if not isinstance(source, dict) or set(source) != {
        "archive_path",
        "archive_sha256",
        "commit_path",
        "commit_sha256",
        "source_bundle_sha256",
        "git_commit",
        "git_tree",
        "worktree_status",
        "submodule_status",
        "submodules",
    }:
        raise ValueError("Clean-candidate source attestation has unexpected fields")
    if not isinstance(build, dict) or set(build) != {
        "profile_id",
        "profile_path",
        "profile_sha256",
        "profile_receipt_path",
        "profile_receipt_sha256",
        "source_archive_sha256",
        "source_commit_sha256",
        "source_bundle_sha256",
        "toolchain",
        "build_invocations_sha256",
        "executable_path",
        "executable_sha256",
    }:
        raise ValueError("Clean-candidate build attestation has unexpected fields")
    source_digest = str(source.get("archive_sha256", ""))
    if sha256_bytes(source_archive) != source_digest:
        raise ValueError("Clean-candidate source archive checksum mismatch")
    source_commit_digest = str(source.get("commit_sha256", ""))
    if sha256_bytes(source_commit) != source_commit_digest:
        raise ValueError("Clean-candidate source commit-object checksum mismatch")
    git_commit = str(source.get("git_commit", ""))
    git_tree = str(source.get("git_tree", ""))
    if not re.fullmatch(r"[0-9a-f]{40}", git_commit):
        raise ValueError("Clean-candidate Git commit is malformed")
    if not re.fullmatch(r"[0-9a-f]{40}", git_tree):
        raise ValueError("Clean-candidate Git tree is malformed")
    if source.get("worktree_status") != "clean":
        raise ValueError("Clean-candidate source worktree is not attested clean")
    freeze_id = str(candidate.get("freeze_id", ""))
    try:
        uuid.UUID(freeze_id)
    except ValueError as error:
        raise ValueError("Clean-candidate freeze ID is malformed") from error
    utc_datetime(candidate.get("created_utc"), field="clean_candidate.created_utc")
    records = source.get("submodules")
    if (
        not isinstance(records, list)
        or len(records) != len(submodule_archives)
        or len(records) != len(submodule_commits)
    ):
        raise ValueError("Clean-candidate submodule archive count mismatch")
    profile_records: list[dict[str, str]] = []
    for index, (record, archive, commit_object) in enumerate(
        zip(records, submodule_archives, submodule_commits)
    ):
        if not isinstance(record, dict) or set(record) != {
            "path",
            "archive_path",
            "archive_sha256",
            "commit_path",
            "commit_sha256",
            "git_commit",
            "git_tree",
            "worktree_status",
        }:
            raise ValueError("Clean-candidate submodule attestation has unexpected fields")
        path = canonical_relative_posix_path(
            record.get("path"), field="Clean-candidate submodule path"
        ).as_posix()
        commit = str(record.get("git_commit", ""))
        tree = str(record.get("git_tree", ""))
        archive_digest = str(record.get("archive_sha256", ""))
        commit_digest = str(record.get("commit_sha256", ""))
        if not re.fullmatch(r"[0-9a-f]{40}", commit):
            raise ValueError("Clean-candidate submodule Git commit is malformed")
        if not re.fullmatch(r"[0-9a-f]{40}", tree):
            raise ValueError("Clean-candidate submodule Git tree is malformed")
        if not re.fullmatch(r"[0-9a-f]{64}", archive_digest):
            raise ValueError("Clean-candidate submodule archive checksum is malformed")
        if not re.fullmatch(r"[0-9a-f]{64}", commit_digest):
            raise ValueError("Clean-candidate submodule commit-object checksum is malformed")
        if record.get("worktree_status") != "clean":
            raise ValueError("Clean-candidate submodule worktree is not attested clean")
        if sha256_bytes(archive) != archive_digest:
            raise ValueError("Clean-candidate submodule archive checksum mismatch")
        if sha256_bytes(commit_object) != commit_digest:
            raise ValueError("Clean-candidate submodule commit-object checksum mismatch")
        if git_archive_commit_from_bytes(archive) != commit:
            raise ValueError("Clean-candidate submodule archive does not identify its commit")
        if git_commit_tree_from_bytes(commit_object, expected_commit=commit) != tree:
            raise ValueError("Clean-candidate submodule commit object does not match its Git tree")
        profile_records.append(
            {
                "path": path,
                "archive_sha256": archive_digest,
                "commit_sha256": commit_digest,
                "git_commit": commit,
                "git_tree": tree,
            }
        )
    paths = [record["path"] for record in profile_records]
    if paths != sorted(set(paths)):
        raise ValueError("Clean-candidate submodules must use unique canonical path order")
    expected_status = "clean_pinned_archived" if records else "absent"
    if source.get("submodule_status") != expected_status:
        raise ValueError("Clean-candidate submodule status does not match its archives")
    for archive, record in zip(submodule_archives, profile_records):
        if (
            git_tree_sha1_from_archive_bytes(
                archive,
                gitlinks=direct_submodule_gitlinks(
                    profile_records, parent_path=record["path"]
                ),
                reject_symlinks=True,
            )
            != record["git_tree"]
        ):
            raise ValueError("Clean-candidate submodule archive does not match its Git tree")
    if git_archive_commit_from_bytes(source_archive) != git_commit:
        raise ValueError("Clean-candidate source archive does not identify its Git commit")
    if git_commit_tree_from_bytes(source_commit, expected_commit=git_commit) != git_tree:
        raise ValueError("Clean-candidate source commit object does not match its Git tree")
    if (
        git_tree_sha1_from_archive_bytes(
            source_archive,
            gitlinks=direct_submodule_gitlinks(profile_records),
            reject_symlinks=True,
        )
        != git_tree
    ):
        raise ValueError("Clean-candidate source archive does not match its Git tree")
    validate_prepared_artifact_closure(
        candidate.get("prepared_artifacts"), source_archive=source_archive
    )
    bundle_digest = source_bundle_sha256(
        source_digest, source_commit_digest, profile_records
    )
    if (
        source.get("source_bundle_sha256") != bundle_digest
        or build.get("source_archive_sha256") != source_digest
        or build.get("source_commit_sha256") != source_commit_digest
        or build.get("source_bundle_sha256") != bundle_digest
        or build.get("executable_sha256") != executable_sha256
    ):
        raise ValueError("Clean-candidate build is not bound to its frozen source bundle")
    if sha256_bytes(build_profile) != build.get("profile_sha256"):
        raise ValueError("Clean-candidate build-profile checksum mismatch")
    if sha256_bytes(build_profile_receipt) != build.get("profile_receipt_sha256"):
        raise ValueError("Clean-candidate build-profile receipt checksum mismatch")
    try:
        profile = read_json_bytes(build_profile, label="clean-candidate build profile")
    except ValueError as error:
        raise ValueError("Clean-candidate build profile is not valid JSON") from error
    for field in ("profile_id", "toolchain"):
        if not isinstance(build.get(field), str) or not str(build[field]).strip():
            raise ValueError(f"Clean-candidate build {field} must not be blank")
    if not re.fullmatch(r"[0-9a-f]{64}", str(build.get("build_invocations_sha256", ""))):
        raise ValueError("Clean-candidate build invocation checksum is malformed")
    expected_provenance_paths = _documented_build_provenance_paths(
        authorized_pic_root=authorized_pic_root,
        git_commit=git_commit,
        profile_id=str(build["profile_id"]),
    )
    if set(build_provenance) != set(BUILD_PROVENANCE_FILENAMES):
        raise ValueError("Frozen build provenance input set is incomplete")
    provenance_inputs = profile.get("provenance_inputs")
    if (
        not isinstance(provenance_inputs, dict)
        or set(provenance_inputs) != set(BUILD_PROVENANCE_FILENAMES)
    ):
        raise ValueError("Frozen build profile provenance input set is incomplete")
    provenance_records: dict[str, dict[str, str]] = {}
    for label in BUILD_PROVENANCE_FILENAMES:
        record = provenance_inputs[label]
        if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
            raise ValueError(f"Malformed frozen build provenance input: {label}")
        path = str(record.get("path", ""))
        digest = str(record.get("sha256", ""))
        if not path or not Path(path).is_absolute():
            raise ValueError(f"Frozen build provenance path is not absolute: {label}")
        if path != str(expected_provenance_paths[label]):
            raise ValueError(
                f"Frozen build provenance path does not use documented Orion layout: {label}"
            )
        if not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError(f"Malformed frozen build provenance checksum: {label}")
        if sha256_bytes(build_provenance[label]) != digest:
            raise ValueError(f"Frozen build provenance checksum mismatch: {label}")
        provenance_records[label] = {"path": path, "sha256": digest}
    if build_provenance["git_status_preconfigure"]:
        raise ValueError("Frozen preconfigure Git status must be empty")
    if build_provenance["git_status"]:
        raise ValueError("Frozen build Git status must be empty")
    for label in set(BUILD_PROVENANCE_FILENAMES) - {
        "git_status_preconfigure",
        "git_status",
        "submodule_status",
    }:
        if not build_provenance[label]:
            raise ValueError(f"Frozen build provenance input must not be empty: {label}")
    try:
        toolchain = build_provenance["toolchain"].decode("utf-8").strip()
        invocations = read_json_bytes(
            build_provenance["build_invocations"], label="frozen build invocations"
        )
    except UnicodeDecodeError as error:
        raise ValueError("Frozen build reviewed provenance must be UTF-8 text") from error
    profile_authorized_source_root = profile.get("authorized_source_root")
    if profile_authorized_source_root != str(authorized_source_root):
        raise ValueError("Frozen build profile does not use the authorized source root")
    if (
        set(invocations) != {"configure", "build"}
        or any(
            not isinstance(invocations[key], list)
            or not invocations[key]
            or any(not isinstance(token, str) or not token for token in invocations[key])
            for key in ["configure", "build"]
        )
    ):
        raise ValueError("Frozen build invocations are malformed")
    require_production_build_provenance(
        authorized_pic_root=authorized_pic_root,
        git_commit=git_commit,
        profile_id=str(build["profile_id"]),
        toolchain=toolchain,
        invocations=invocations,
        module_list=build_provenance["module_list"],
        environment_allowlist=build_provenance["environment_allowlist"],
        build_environment=build_provenance["build_environment"],
    )
    invocation_sha256 = sha256_bytes(build_provenance["build_invocations"])
    fresh_source_root = (
        authorized_pic_root
        / "build"
        / git_commit[:12]
        / str(build["profile_id"])
        / "source"
    )
    if (
        type(profile.get("schema_version")) is not int
        or profile != {
        "schema_version": 3,
        "profile_id": build.get("profile_id"),
        "authorized_source_root": str(authorized_source_root),
        "fresh_source_root": str(fresh_source_root),
        "git_commit": git_commit,
        "git_tree": git_tree,
        "source_archive_sha256": source_digest,
        "source_commit_sha256": source_commit_digest,
        "source_bundle_sha256": bundle_digest,
        "toolchain": toolchain,
        "build_invocations_sha256": invocation_sha256,
        "executable_sha256": executable_sha256,
        "provenance_inputs": provenance_records,
        "submodules": profile_records,
        }
    ):
        raise ValueError("Frozen build profile does not match clean-candidate attestation")
    if (
        build.get("toolchain") != toolchain
        or build.get("build_invocations_sha256") != invocation_sha256
    ):
        raise ValueError("Clean-candidate build metadata differs from frozen provenance")
    try:
        receipt = read_json_bytes(
            build_profile_receipt, label="clean-candidate build-profile receipt"
        )
    except ValueError as error:
        raise ValueError("Clean-candidate build-profile receipt is not valid JSON") from error
    artifact_dir = (
        authorized_pic_root / "bin" / git_commit[:12] / str(build["profile_id"])
    )
    control_plane_version = str(receipt.get("control_plane_version", ""))
    if not re.fullmatch(r"[0-9a-f]{64}", control_plane_version):
        raise ValueError("Build-profile receipt control-plane version is malformed")
    if (
        expected_control_plane_version is not None
        and control_plane_version != expected_control_plane_version
    ):
        raise ValueError("Build-profile receipt belongs to another control-plane version")
    expected_receipt = {
        "schema_version": 1,
        "control_plane_version": control_plane_version,
        "profile_path": str(artifact_dir / "build_profile.json"),
        "profile_sha256": sha256_bytes(build_profile),
        "source_bundle_sha256": bundle_digest,
        "fresh_source_root": str(fresh_source_root),
        "build_invocations_sha256": invocation_sha256,
        "git_status_preconfigure_sha256": provenance_records[
            "git_status_preconfigure"
        ]["sha256"],
        "git_status_sha256": provenance_records["git_status"]["sha256"],
        "configure_log_sha256": provenance_records["configure_log"]["sha256"],
        "build_log_sha256": provenance_records["build_log"]["sha256"],
        "executable_path": str(artifact_dir / "athena"),
        "executable_sha256": executable_sha256,
    }
    if type(receipt.get("schema_version")) is not int or receipt != expected_receipt:
        raise ValueError("Frozen build-profile receipt does not match clean candidate")
    return profile_records


def read_json(path: Path) -> dict[str, object]:
    return read_json_bytes(read_stable_regular_file(path), label=str(path))


def utc_datetime(value: object, *, field: str) -> datetime:
    if not isinstance(value, str):
        raise ValueError(f"{field} must be a canonical RFC-3339 UTC timestamp")
    text = value
    if re.fullmatch(
        r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}"
        r"(?:\.[0-9]{1,6})?Z",
        text,
    ) is None:
        raise ValueError(f"{field} must be a canonical RFC-3339 UTC timestamp")
    try:
        result = datetime.fromisoformat(text[:-1] + "+00:00")
    except ValueError as error:
        raise ValueError(f"Invalid {field}: {text}") from error
    if result.tzinfo != timezone.utc:
        raise ValueError(f"{field} must use UTC")
    return result


def write_json_exclusive(path: Path, value: dict[str, object]) -> None:
    durable_mkdir_parents(path.parent)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def fsync_directory(path: Path) -> None:
    """Sync one existing directory without following a symlink alias."""
    descriptor = os.open(
        path,
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def durable_mkdir_parents(
    path: Path, *, mode: int = 0o755, root: Path | None = None
) -> None:
    """Create missing descendants below an optional trusted lexical root."""
    path = Path(os.path.abspath(path))
    if root is None:
        anchor = path
        while not anchor.exists():
            if anchor.parent == anchor:
                raise FileNotFoundError(
                    f"No existing ancestor for directory creation: {path}"
                )
            anchor = anchor.parent
    else:
        lexical_root = Path(os.path.abspath(root))
        try:
            path.relative_to(lexical_root)
        except ValueError as error:
            raise ValueError(
                f"Directory path is outside trusted lexical root: {path}"
            ) from error
        anchor = lexical_root
        while not anchor.exists():
            if anchor.parent == anchor:
                raise FileNotFoundError(
                    f"No existing trusted-root ancestor for directory creation: {path}"
                )
            anchor = anchor.parent
    try:
        relative = path.relative_to(anchor)
    except ValueError as error:
        raise ValueError(f"Directory path is outside trusted lexical root: {path}") from error
    flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    anchor_flags = flags
    if root is not None and anchor == Path(os.path.abspath(root)):
        anchor_flags = os.O_RDONLY | os.O_DIRECTORY
    parent_descriptor = os.open(anchor, anchor_flags)
    try:
        for part in relative.parts:
            created = False
            try:
                child_descriptor = os.open(part, flags, dir_fd=parent_descriptor)
            except FileNotFoundError:
                os.mkdir(part, mode=mode, dir_fd=parent_descriptor)
                child_descriptor = os.open(part, flags, dir_fd=parent_descriptor)
                created = True
            if created:
                os.fsync(child_descriptor)
                os.fsync(parent_descriptor)
            os.close(parent_descriptor)
            parent_descriptor = child_descriptor
    finally:
        os.close(parent_descriptor)


def _open_parent_directory(
    path: Path, *, root: Path | None = None
) -> int:
    """Open an output parent through no-follow traversal below a trusted root."""
    path = Path(os.path.abspath(path))
    anchor = Path("/") if root is None else Path(os.path.abspath(root))
    try:
        relative = path.relative_to(anchor)
    except ValueError as error:
        raise ValueError(f"Output parent is outside trusted lexical root: {path}") from error
    flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    parent_descriptor = os.open(
        anchor,
        os.O_RDONLY | os.O_DIRECTORY if root is not None else flags,
    )
    try:
        for part in relative.parts:
            child_descriptor = os.open(part, flags, dir_fd=parent_descriptor)
            os.close(parent_descriptor)
            parent_descriptor = child_descriptor
    except BaseException:
        os.close(parent_descriptor)
        raise
    return parent_descriptor


def open_directory_below(path: Path, *, root: Path) -> int:
    """Open one directory through anchored no-follow traversal below root."""
    path = Path(os.path.abspath(path))
    root = Path(os.path.abspath(root))
    if path == root:
        return os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    parent_descriptor = _open_parent_directory(path.parent, root=root)
    try:
        return os.open(
            path.name,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_descriptor,
        )
    finally:
        os.close(parent_descriptor)


class PinnedDirectoryAncestry:
    """Retain and recheck every directory component below one stable anchor."""

    def __init__(self, path: Path, *, root: Path) -> None:
        self.path = Path(os.path.abspath(path))
        self.root = Path(os.path.abspath(root))
        try:
            self.relative = self.path.relative_to(self.root)
        except ValueError as error:
            raise ValueError(
                f"Directory path is outside trusted lexical root: {self.path}"
            ) from error
        self._descriptors: list[tuple[Path, int]] = []
        flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
        descriptor = os.open(self.root, flags)
        self._descriptors.append((self.root, descriptor))
        try:
            current = self.root
            for part in self.relative.parts:
                descriptor = os.open(part, flags, dir_fd=descriptor)
                current /= part
                self._descriptors.append((current, descriptor))
            self.require_same()
        except BaseException:
            self.close()
            raise

    @property
    def descriptor(self) -> int:
        if not self._descriptors:
            raise ValueError("Pinned directory ancestry is closed")
        return self._descriptors[-1][1]

    def __enter__(self) -> "PinnedDirectoryAncestry":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def require_same(self) -> None:
        if not self._descriptors:
            raise ValueError("Pinned directory ancestry is closed")
        flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
        lexical_descriptor = os.open(self.root, flags)
        try:
            for index, (_, expected_descriptor) in enumerate(self._descriptors):
                expected = os.fstat(expected_descriptor)
                actual = os.fstat(lexical_descriptor)
                if (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino):
                    raise ValueError(
                        f"Directory ancestry changed below trusted root: {self.path}"
                    )
                if index + 1 < len(self._descriptors):
                    child_descriptor = os.open(
                        self._descriptors[index + 1][0].name,
                        flags,
                        dir_fd=lexical_descriptor,
                    )
                    os.close(lexical_descriptor)
                    lexical_descriptor = child_descriptor
        finally:
            os.close(lexical_descriptor)

    def close(self) -> None:
        for _, descriptor in reversed(self._descriptors):
            os.close(descriptor)
        self._descriptors.clear()


def atomic_write_bytes(
    path: Path,
    data: bytes,
    *,
    mode: int = 0o444,
    replace: bool = True,
    root: Path | None = None,
) -> None:
    path = Path(os.path.abspath(path))
    if not path.name or path.name in {".", ".."}:
        raise ValueError("Atomic output path must name one file")
    durable_mkdir_parents(path.parent, root=root)
    parent_descriptor = _open_parent_directory(path.parent, root=root)
    try:
        _require_same_directory(path.parent, parent_descriptor, root=root)
        atomic_write_bytes_at(
            parent_descriptor,
            path.name,
            data,
            mode=mode,
            replace=replace,
            post_publish_check=lambda: _require_same_directory(
                path.parent, parent_descriptor, root=root
            ),
        )
    finally:
        os.close(parent_descriptor)


def atomic_write_bytes_at(
    parent_descriptor: int,
    name: str,
    data: bytes,
    *,
    mode: int = 0o444,
    replace: bool = True,
    post_publish_check: Callable[[], None] | None = None,
) -> None:
    """Publish one file relative to a retained trusted parent descriptor."""
    if not name or name in {".", ".."} or Path(name).name != name:
        raise ValueError("Atomic output path must name one file")
    temporary_name = f".{name}.tmp-{uuid.uuid4()}"
    temporary_exists = False
    rollback_name: str | None = None
    try:
        descriptor = os.open(
            temporary_name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL,
            0o600,
            dir_fd=parent_descriptor,
        )
        temporary_exists = True
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
            os.fchmod(stream.fileno(), mode)
            os.fsync(stream.fileno())
        if replace:
            rollback_name = f".{name}.rollback-{uuid.uuid4()}"
            try:
                os.link(
                    name,
                    rollback_name,
                    src_dir_fd=parent_descriptor,
                    dst_dir_fd=parent_descriptor,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                rollback_name = None
            replacement_visible = False
            try:
                os.replace(
                    temporary_name,
                    name,
                    src_dir_fd=parent_descriptor,
                    dst_dir_fd=parent_descriptor,
                )
                temporary_exists = False
                replacement_visible = True
                os.fsync(parent_descriptor)
                if post_publish_check is not None:
                    post_publish_check()
            except BaseException as error:
                if replacement_visible:
                    try:
                        if rollback_name is None:
                            os.unlink(name, dir_fd=parent_descriptor)
                        else:
                            os.replace(
                                rollback_name,
                                name,
                                src_dir_fd=parent_descriptor,
                                dst_dir_fd=parent_descriptor,
                            )
                            rollback_name = None
                        os.fsync(parent_descriptor)
                    except BaseException as rollback_error:
                        raise RuntimeError(
                            f"Failed to roll back replacement publication: {name}"
                        ) from rollback_error
                raise error
            if rollback_name is not None:
                os.unlink(rollback_name, dir_fd=parent_descriptor)
                rollback_name = None
        else:
            os.link(
                temporary_name,
                name,
                src_dir_fd=parent_descriptor,
                dst_dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
            try:
                os.unlink(temporary_name, dir_fd=parent_descriptor)
                temporary_exists = False
                os.fsync(parent_descriptor)
                if post_publish_check is not None:
                    post_publish_check()
            except BaseException as error:
                try:
                    os.unlink(name, dir_fd=parent_descriptor)
                    os.fsync(parent_descriptor)
                except BaseException as rollback_error:
                    raise RuntimeError(
                        f"Failed to roll back exclusive publication: {name}"
                    ) from rollback_error
                raise error
    finally:
        if temporary_exists:
            try:
                os.unlink(temporary_name, dir_fd=parent_descriptor)
            except FileNotFoundError:
                pass
        if rollback_name is not None:
            try:
                os.unlink(rollback_name, dir_fd=parent_descriptor)
            except FileNotFoundError:
                pass


def atomic_write_json(
    path: Path,
    value: dict[str, object],
    *,
    mode: int = 0o444,
    replace: bool = True,
    root: Path | None = None,
) -> None:
    data = (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    atomic_write_bytes(path, data, mode=mode, replace=replace, root=root)


def atomic_write_json_at(
    parent_descriptor: int,
    name: str,
    value: dict[str, object],
    *,
    mode: int = 0o444,
    replace: bool = True,
    post_publish_check: Callable[[], None] | None = None,
) -> None:
    data = (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    atomic_write_bytes_at(
        parent_descriptor,
        name,
        data,
        mode=mode,
        replace=replace,
        post_publish_check=post_publish_check,
    )


def _descriptor_path(parent_descriptor: int, name: str) -> Path:
    return Path("/proc") / str(os.getpid()) / "fd" / str(parent_descriptor) / name


def _require_same_directory(path: Path, descriptor: int, *, root: Path | None) -> None:
    lexical_descriptor = _open_parent_directory(path, root=root)
    try:
        expected = os.fstat(descriptor)
        actual = os.fstat(lexical_descriptor)
        if (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino):
            raise ValueError(f"Writable parent path changed during publication: {path}")
    finally:
        os.close(lexical_descriptor)


def require_same_directory(
    path: Path, descriptor: int, *, root: Path | None = None
) -> None:
    """Reject replacement of a lexical directory retained by descriptor."""
    _require_same_directory(path, descriptor, root=root)


class PinnedStagingDirectory:
    """Stage and publish beneath one retained no-follow parent descriptor."""

    def __init__(self, parent: Path, *, prefix: str, root: Path | None = None) -> None:
        self.parent = Path(os.path.abspath(parent))
        self.root = Path(os.path.abspath(root)) if root is not None else None
        self.name = f"{prefix}{uuid.uuid4()}"
        self.parent_descriptor: int | None = None
        self.path: Path | None = None
        self._staging_exists = False

    def __enter__(self) -> "PinnedStagingDirectory":
        durable_mkdir_parents(self.parent, root=self.root)
        self.parent_descriptor = _open_parent_directory(self.parent, root=self.root)
        try:
            os.mkdir(self.name, mode=0o700, dir_fd=self.parent_descriptor)
            self._staging_exists = True
            self.path = _descriptor_path(self.parent_descriptor, self.name)
            self.require_lexical_parent()
            return self
        except BaseException:
            if self._staging_exists:
                remove_tree(self._path())
                os.fsync(self.parent_descriptor)
                self._staging_exists = False
            os.close(self.parent_descriptor)
            self.parent_descriptor = None
            raise

    def _descriptor(self) -> int:
        if self.parent_descriptor is None:
            raise ValueError("Pinned staging directory is not open")
        return self.parent_descriptor

    def _path(self) -> Path:
        if self.path is None:
            raise ValueError("Pinned staging directory is not open")
        return self.path

    def require_lexical_parent(self) -> None:
        _require_same_directory(self.parent, self._descriptor(), root=self.root)

    def publish_tree(self, destination: Path) -> None:
        destination = Path(os.path.abspath(destination))
        if destination.parent != self.parent:
            raise ValueError("Tree destination does not use the pinned staging parent")
        durable_replace_tree(
            self._path(),
            destination,
            parent_descriptor=self._descriptor(),
            temporary_name=self.name,
            root=self.root,
        )
        self._staging_exists = False

    def write_json(
        self,
        destination: Path,
        value: dict[str, object],
        *,
        mode: int = 0o444,
        replace: bool = True,
    ) -> None:
        destination = Path(os.path.abspath(destination))
        if destination.parent != self.parent:
            raise ValueError("JSON destination does not use the pinned staging parent")
        self.require_lexical_parent()
        atomic_write_json_at(
            self._descriptor(),
            destination.name,
            value,
            mode=mode,
            replace=replace,
            post_publish_check=self.require_lexical_parent,
        )

    def __exit__(self, *_: object) -> None:
        try:
            if self._staging_exists:
                remove_tree(self._path())
                os.fsync(self._descriptor())
        finally:
            if self.parent_descriptor is not None:
                os.close(self.parent_descriptor)
                self.parent_descriptor = None


def durable_replace_tree(
    temporary: Path,
    destination: Path,
    *,
    parent_descriptor: int | None = None,
    temporary_name: str | None = None,
    root: Path | None = None,
) -> None:
    """Publish one staged tree only after its files and directories are durable."""
    temporary = Path(os.path.abspath(temporary))
    destination = Path(os.path.abspath(destination))
    owns_parent_descriptor = parent_descriptor is None
    if parent_descriptor is None:
        if temporary.parent != destination.parent:
            raise ValueError("Staged tree and destination must share one parent")
        parent_descriptor = _open_parent_directory(destination.parent, root=root)
    source_name = temporary.name if temporary_name is None else temporary_name
    temporary = _descriptor_path(parent_descriptor, source_name)
    stable_destination = _descriptor_path(parent_descriptor, destination.name)
    paths = [temporary, *temporary.rglob("*")]
    for path in paths:
        if path.is_symlink():
            raise ValueError(f"Durable tree publication rejects symlinks: {path}")
    for path in paths:
        if not path.is_file():
            continue
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
        try:
            if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                raise ValueError(f"Durable tree entry is not a regular file: {path}")
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    directories = sorted(
        (path for path in paths if path.is_dir()),
        key=lambda path: len(path.parts),
        reverse=True,
    )
    for path in directories:
        descriptor = os.open(
            path,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
        )
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    try:
        _require_same_directory(destination.parent, parent_descriptor, root=root)
        os.replace(
            source_name,
            destination.name,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
        )
        try:
            os.fsync(parent_descriptor)
            _require_same_directory(destination.parent, parent_descriptor, root=root)
        except BaseException as error:
            try:
                remove_tree(stable_destination)
                os.fsync(parent_descriptor)
            except BaseException as rollback_error:
                raise RuntimeError(
                    f"Failed to roll back durable tree publication: {destination}"
                ) from rollback_error
            raise error
    finally:
        if owns_parent_descriptor:
            os.close(parent_descriptor)


def make_tree_read_only(root: Path, *, executable_names: set[str] | None = None) -> None:
    executables = executable_names or set()
    for path in sorted(root.rglob("*"), reverse=True):
        if path.is_dir():
            path.chmod(0o555)
        else:
            path.chmod(0o555 if path.name in executables else 0o444)
    root.chmod(0o555)


def remove_tree(root: Path) -> None:
    """Remove a generated tree even if it was already staged read-only."""
    if not root.exists():
        return
    for path in root.rglob("*"):
        if path.is_dir():
            path.chmod(0o700)
    root.chmod(0o700)
    shutil.rmtree(root)


def require_read_only(path: Path) -> None:
    if path.stat().st_mode & 0o222:
        raise ValueError(f"Artifact is not read-only: {path}")


def require_not_symlink(path: Path) -> None:
    if path.is_symlink():
        raise ValueError(f"Path must not be a symlink alias: {path}")


def canonical_policy_path(authorized_pic_root: Path = AUTHORIZED_PIC_ROOT) -> Path:
    return Path(os.path.abspath(authorized_pic_root)) / CANONICAL_POLICY_RELATIVE


def active_promotion_path(root: Path) -> Path:
    return Path(os.path.abspath(root)) / ACTIVE_PROMOTION_RELATIVE


def require_below(path: Path, root: Path) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError as error:
        raise ValueError(f"Path is outside authorized PIC root: {resolved}") from error
    return resolved


def require_no_symlink_components_below(path: Path, root: Path) -> Path:
    """Reject pre-existing symlink traversal below one trusted lexical root."""
    lexical_root = Path(os.path.abspath(root))
    lexical_path = Path(os.path.abspath(path))
    try:
        relative = lexical_path.relative_to(lexical_root)
    except ValueError as error:
        raise ValueError(f"Path is outside authorized lexical root: {lexical_path}") from error
    current = lexical_root
    for part in relative.parts:
        current /= part
        if current.is_symlink():
            raise ValueError(f"Path traverses a symlink below authorized root: {current}")
    return lexical_path


def require_canonical_path_below(path: Path, root: Path) -> Path:
    """Require one lexical path under root with no symlink aliases."""
    lexical = require_no_symlink_components_below(path, root)
    lexical_root = Path(os.path.abspath(root))
    relative = lexical.relative_to(lexical_root)
    resolved = lexical.resolve()
    expected = root.resolve().joinpath(*relative.parts)
    if resolved != expected:
        raise ValueError(f"Path must use its trusted-root spelling: {lexical}; expected {expected}")
    return lexical


def _relative_artifact_path(value: object, *, field: str) -> str:
    text = str(value)
    path = PurePosixPath(text)
    if (
        not text
        or path.is_absolute()
        or path.as_posix() != text
        or not path.parts
        or path.parts != tuple(part for part in path.parts if part not in {"", ".", ".."})
    ):
        raise ValueError(f"{field} must be a non-empty relative artifact path")
    return text


def validate_launch_contract(value: object) -> dict[str, object]:
    """Validate the closed structured argv contract executed by the trampoline."""
    if not isinstance(value, dict) or set(value) != {
        "schema_version",
        "executor",
        "pre_actions",
        "actions",
        "post_actions",
    }:
        raise ValueError(
            "Launch contract must contain only schema_version, executor, "
            "pre_actions, actions and post_actions"
        )
    if type(value.get("schema_version")) is not int or value.get("schema_version") != 1:
        raise ValueError("Unsupported launch-contract schema")
    if value.get("executor") != TRUSTED_LAUNCH_EXECUTOR:
        raise ValueError("Launch contract does not select the trusted Athena executor")
    actions = value.get("actions")
    if not isinstance(actions, list) or not 1 <= len(actions) <= 16:
        raise ValueError("Launch contract requires between one and sixteen Athena actions")
    identifiers = set()
    for phase in ["pre_actions", "post_actions"]:
        bounded_actions = value.get(phase)
        if not isinstance(bounded_actions, list) or len(bounded_actions) > 16:
            raise ValueError(f"Launch contract {phase} must contain at most sixteen actions")
        for action in bounded_actions:
            if not isinstance(action, dict):
                raise ValueError(f"Launch-contract {phase} action is malformed")
            identifier = str(action.get("action_id", ""))
            if not re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,63}", identifier):
                raise ValueError(f"Launch-contract {phase} action ID is malformed")
            if identifier in identifiers:
                raise ValueError(f"Duplicate launch action ID: {identifier}")
            identifiers.add(identifier)
            kind = action.get("kind")
            if kind == "snapshot_sha256":
                if set(action) != {
                    "action_id",
                    "kind",
                    "snapshot_role",
                    "output_artifact",
                }:
                    raise ValueError(f"Launch-contract {phase} snapshot action is malformed")
                role = str(action.get("snapshot_role", ""))
                if role not in {
                    "job-script",
                    "executable",
                    "input-deck",
                    "environment-profile",
                    "timeout-margin",
                    "queue-snapshot",
                } and not re.fullmatch(r"analysis-script-[0-9]{3}", role):
                    raise ValueError(f"Launch-contract {phase} snapshot role is malformed")
                _relative_artifact_path(
                    action.get("output_artifact"), field="output_artifact"
                )
            elif kind == "artifact_sha256":
                if set(action) != {
                    "action_id",
                    "kind",
                    "artifact",
                    "output_artifact",
                }:
                    raise ValueError(f"Launch-contract {phase} artifact action is malformed")
                artifact = _relative_artifact_path(action.get("artifact"), field="artifact")
                output = _relative_artifact_path(
                    action.get("output_artifact"), field="output_artifact"
                )
                if artifact == output:
                    raise ValueError("Artifact checksum output must differ from its input")
            elif kind == "artifact_nonempty":
                if set(action) != {"action_id", "kind", "artifact"}:
                    raise ValueError(f"Launch-contract {phase} assertion action is malformed")
                _relative_artifact_path(action.get("artifact"), field="artifact")
            else:
                raise ValueError(f"Launch-contract {phase} accepts only bounded built-in actions")
    for action in actions:
        if not isinstance(action, dict) or set(action) != {
            "action_id",
            "kind",
            "resources",
            "arguments",
            "stdout_artifact",
            "stderr_artifact",
        }:
            raise ValueError("Launch action has unexpected or missing fields")
        identifier = str(action.get("action_id", ""))
        if not re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,63}", identifier):
            raise ValueError("Launch action ID is malformed")
        if identifier in identifiers:
            raise ValueError(f"Duplicate launch action ID: {identifier}")
        identifiers.add(identifier)
        if action.get("kind") != "athena":
            raise ValueError("Trusted trampoline accepts only Athena actions")
        resources = action.get("resources")
        if not isinstance(resources, dict) or set(resources) != {
            "nodes",
            "tasks",
            "cpus_per_task",
            "gpus_per_task",
            "gpu_bind",
        }:
            raise ValueError("Launch-action resources are malformed")
        for field in ["nodes", "tasks", "cpus_per_task", "gpus_per_task"]:
            number = resources.get(field)
            if not isinstance(number, int) or isinstance(number, bool) or number <= 0:
                raise ValueError(f"Launch-action resources.{field} must be a positive integer")
        if resources.get("gpu_bind") != "closest":
            raise ValueError("Launch-action resources.gpu_bind must be closest")
        arguments = action.get("arguments")
        if not isinstance(arguments, list):
            raise ValueError("Launch-action arguments must be an array")
        argument_index = 0
        while argument_index < len(arguments):
            argument = arguments[argument_index]
            if not isinstance(argument, dict) or len(argument) != 1:
                raise ValueError("Launch-action argument must be one structured token")
            if "literal" in argument:
                literal = argument["literal"]
                if (
                    not isinstance(literal, str)
                    or not literal
                    or "\0" in literal
                    or len(literal) > 4096
                ):
                    raise ValueError("Launch-action literal is malformed")
                if literal in {"-i", "-d"}:
                    argument_index += 1
                    if argument_index >= len(arguments):
                        raise ValueError(f"Launch-action {literal} requires one value")
                    structured = arguments[argument_index]
                    if literal == "-i" and structured != {"snapshot_role": "input-deck"}:
                        raise ValueError("Launch-action -i requires the frozen input deck")
                    if literal == "-d":
                        if (
                            not isinstance(structured, dict)
                            or set(structured) != {"artifact_directory"}
                        ):
                            raise ValueError(
                                "Launch-action -d requires one artifact directory"
                            )
                        _relative_artifact_path(
                            structured["artifact_directory"],
                            field="artifact_directory",
                        )
                elif literal == "-r" or literal.startswith("-r="):
                    raise ValueError(
                        "Launch-action restart input is not authorized before trusted "
                        "restart snapshots are implemented"
                    )
                elif literal in {"-n", "-c"}:
                    pass
                elif literal.startswith("-"):
                    raise ValueError(f"Launch-action CLI flag is not authorized: {literal}")
                else:
                    key, separator, override = literal.partition("=")
                    if (
                        not separator
                        or not re.fullmatch(r"[A-Za-z0-9_]+(?:/[A-Za-z0-9_]+)+", key)
                        or not override
                        or "/" in override
                        or "\\" in override
                        or ".." in override
                    ):
                        raise ValueError(
                            f"Launch-action Athena override is not authorized: {literal}"
                        )
            else:
                raise ValueError(
                    "Launch-action structured path token must immediately follow -i or -d"
                )
            argument_index += 1
        stdout = _relative_artifact_path(
            action.get("stdout_artifact"), field="stdout_artifact"
        )
        stderr = _relative_artifact_path(
            action.get("stderr_artifact"), field="stderr_artifact"
        )
        if stdout == stderr:
            raise ValueError("Launch-action stdout and stderr artifacts must differ")
    return value


def launch_contract_sha256(value: object) -> str:
    contract = validate_launch_contract(value)
    payload = json.dumps(
        contract, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return sha256_bytes(payload)


def scrub_file(path: Path) -> None:
    data = path.read_bytes()
    if PLACEHOLDER_PATTERN.search(data):
        raise ValueError(f"Unresolved REPLACE_* placeholder in {path}")
    if SENSITIVE_PATTERN.search(data):
        raise ValueError(f"Potential sensitive string in submission artifact: {path}")


def snapshot_file(
    source: Path,
    destination: Path,
    *,
    role: str,
    destination_root: Path,
    scrub: bool = True,
) -> dict[str, str]:
    if not source.is_file():
        raise FileNotFoundError(f"Missing snapshot source for {role}: {source}")
    if scrub:
        scrub_file(source)
    try:
        destination.resolve().relative_to(destination_root.resolve())
    except ValueError as error:
        raise ValueError(
            f"Snapshot destination for {role} escapes staging root: {destination}"
        ) from error
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    if scrub:
        scrub_file(destination)
    return {
        "role": role,
        "path": str(destination),
        "sha256": sha256(destination),
        "source_path": str(source.resolve()),
        "source_sha256": sha256(source),
    }


def verify_snapshot_files(manifest: dict[str, object], *, root: Path) -> None:
    snapshot_files = manifest.get("snapshot_files")
    if not isinstance(snapshot_files, list) or not snapshot_files:
        raise ValueError("Manifest has no snapshot_files")
    for record in snapshot_files:
        if not isinstance(record, dict) or set(record) != {
            "role",
            "path",
            "sha256",
            "source_path",
            "source_sha256",
        }:
            raise ValueError("Malformed snapshot file record")
        if (
            not isinstance(record["role"], str)
            or not record["role"]
            or not isinstance(record["path"], str)
            or not record["path"]
            or not isinstance(record["source_path"], str)
            or not record["source_path"]
            or not _is_lowercase_sha256(record["sha256"])
            or not _is_lowercase_sha256(record["source_sha256"])
        ):
            raise ValueError("Malformed snapshot file record")
        path = Path(record["path"])
        require_canonical_path_below(path, root)
        data = read_stable_regular_file_below(path, root)
        if sha256_bytes(data) != record.get("sha256"):
            raise ValueError(f"Snapshot checksum mismatch: {path}")


def record_for_role(
    manifest: dict[str, object], role: str
) -> dict[str, object]:
    snapshot_files = manifest.get("snapshot_files", [])
    matches = [
        record for record in snapshot_files
        if isinstance(record, dict) and record.get("role") == role
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one snapshot role={role}")
    return matches[0]


def checksum_records(paths: Iterable[Path]) -> list[dict[str, str]]:
    return [{"path": str(path.resolve()), "sha256": sha256(path)} for path in paths]


def inventory_digest(records: list[dict[str, str]]) -> str:
    payload = json.dumps(records, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def validate_control_plane_inventory(inventory: dict[str, object]) -> list[dict[str, str]]:
    """Require the exact closed inventory shape used by installed generations."""
    if (
        set(inventory) != {"schema_version", "version", "files"}
        or type(inventory.get("schema_version")) is not int
        or inventory.get("schema_version") != 1
        or not isinstance(inventory.get("version"), str)
        or re.fullmatch(r"[0-9a-f]{64}", str(inventory["version"])) is None
    ):
        raise ValueError("Unsupported control-plane inventory schema")
    raw_records = inventory["files"]
    if not isinstance(raw_records, list):
        raise ValueError("Malformed control-plane inventory files")
    records: list[dict[str, str]] = []
    for raw in raw_records:
        if (
            not isinstance(raw, dict)
            or set(raw) != {"path", "sha256"}
            or not isinstance(raw["path"], str)
            or not isinstance(raw["sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", raw["sha256"]) is None
        ):
            raise ValueError("Malformed control-plane inventory record")
        records.append({"path": raw["path"], "sha256": raw["sha256"]})
    if [record["path"] for record in records] != CONTROL_PLANE_FILES:
        raise ValueError("Control-plane inventory file list differs from required list")
    if inventory["version"] != inventory_digest(records):
        raise ValueError("Control-plane inventory digest mismatch")
    return records


def validate_historical_control_plane_inventory(
    inventory: dict[str, object],
) -> list[dict[str, str]]:
    """Require a closed self-authenticating inventory for an older generation."""
    if (
        set(inventory) != {"schema_version", "version", "files"}
        or type(inventory.get("schema_version")) is not int
        or inventory.get("schema_version") != 1
        or not isinstance(inventory.get("version"), str)
        or re.fullmatch(r"[0-9a-f]{64}", str(inventory["version"])) is None
    ):
        raise ValueError("Unsupported historical control-plane inventory schema")
    raw_records = inventory["files"]
    if not isinstance(raw_records, list) or not raw_records:
        raise ValueError("Malformed historical control-plane inventory files")
    records: list[dict[str, str]] = []
    for raw in raw_records:
        if (
            not isinstance(raw, dict)
            or set(raw) != {"path", "sha256"}
            or not isinstance(raw["path"], str)
            or not raw["path"]
            or "/" in raw["path"]
            or Path(raw["path"]).name != raw["path"]
            or raw["path"] == "inventory.json"
            or not isinstance(raw["sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", raw["sha256"]) is None
        ):
            raise ValueError("Malformed historical control-plane inventory record")
        records.append({"path": raw["path"], "sha256": raw["sha256"]})
    if len({record["path"] for record in records}) != len(records):
        raise ValueError("Historical control-plane inventory has duplicate files")
    if inventory["version"] != inventory_digest(records):
        raise ValueError("Historical control-plane inventory digest mismatch")
    return records


def _read_read_only_regular_file_at(
    directory_descriptor: int, name: str, *, label: str
) -> bytes:
    if not name or "/" in name or Path(name).name != name:
        raise ValueError(f"{label} has an invalid installed name")
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_descriptor,
    )
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"{label} is not a regular file")
        if metadata.st_mode & 0o222:
            raise ValueError(f"{label} is not read-only")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            return stream.read()
    finally:
        os.close(descriptor)


def verify_installed_control_plane(
    control_plane_dir: Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    lexical_root = Path(os.path.abspath(authorized_pic_root))
    expected_parent = lexical_root / "control_plane"
    resolved = require_canonical_path_below(
        control_plane_dir, lexical_root
    )
    if resolved.parent != expected_parent:
        raise ValueError(
            f"Control plane is not installed under authorized root: {resolved}"
        )
    try:
        directory_descriptor = open_directory_below(resolved, root=lexical_root)
    except FileNotFoundError as error:
        raise ValueError(f"Missing installed control-plane directory: {resolved}") from error
    try:
        if os.fstat(directory_descriptor).st_mode & 0o222:
            raise ValueError(f"Installed control-plane directory is not read-only: {resolved}")
        expected_names = {*CONTROL_PLANE_FILES, "inventory.json"}
        if set(os.listdir(directory_descriptor)) != expected_names:
            raise ValueError("Installed control-plane entries differ from required list")
        inventory = read_json_bytes(
            _read_read_only_regular_file_at(
                directory_descriptor,
                "inventory.json",
                label="Installed control-plane inventory",
            ),
            label=str(resolved / "inventory.json"),
        )
        records = validate_control_plane_inventory(inventory)
        if resolved.name != inventory["version"]:
            raise ValueError("Control-plane inventory digest mismatch")
        for record in records:
            data = _read_read_only_regular_file_at(
                directory_descriptor,
                record["path"],
                label=f"Installed control-plane file {record['path']}",
            )
            if sha256_bytes(data) != record["sha256"]:
                raise ValueError(
                    f"Installed control-plane checksum mismatch: {resolved / record['path']}"
                )
        return inventory
    finally:
        os.close(directory_descriptor)


def verify_historical_installed_control_plane(
    control_plane_dir: Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    """Verify an immutable predecessor using its own closed inventory."""
    lexical_root = Path(os.path.abspath(authorized_pic_root))
    expected_parent = lexical_root / "control_plane"
    resolved = require_canonical_path_below(control_plane_dir, lexical_root)
    if resolved.parent != expected_parent:
        raise ValueError(
            f"Historical control plane is not installed under authorized root: {resolved}"
        )
    try:
        directory_descriptor = open_directory_below(resolved, root=lexical_root)
    except FileNotFoundError as error:
        raise ValueError(
            f"Missing historical installed control-plane directory: {resolved}"
        ) from error
    try:
        if os.fstat(directory_descriptor).st_mode & 0o222:
            raise ValueError(
                f"Historical installed control-plane directory is not read-only: {resolved}"
            )
        inventory = read_json_bytes(
            _read_read_only_regular_file_at(
                directory_descriptor,
                "inventory.json",
                label="Historical installed control-plane inventory",
            ),
            label=str(resolved / "inventory.json"),
        )
        records = validate_historical_control_plane_inventory(inventory)
        expected_names = {*(record["path"] for record in records), "inventory.json"}
        if set(os.listdir(directory_descriptor)) != expected_names:
            raise ValueError("Historical installed control-plane entries differ from inventory")
        if resolved.name != inventory["version"]:
            raise ValueError("Historical control-plane inventory digest mismatch")
        for record in records:
            data = _read_read_only_regular_file_at(
                directory_descriptor,
                record["path"],
                label=f"Historical installed control-plane file {record['path']}",
            )
            if sha256_bytes(data) != record["sha256"]:
                raise ValueError(
                    f"Historical installed control-plane checksum mismatch: "
                    f"{resolved / record['path']}"
                )
        return inventory
    finally:
        os.close(directory_descriptor)


def validate_storage_policy(
    policy: dict[str, object],
    *,
    control_plane_version: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
    ledger_mirror_transport: str = AUTHORIZED_LEDGER_MIRROR_TRANSPORT,
    allow_pending_genesis: bool = False,
) -> dict[str, object]:
    allowed_policy_keys = {
        "schema_version",
        "frontier",
        "science_submission_freeze",
        "registered_science_slices",
        "frontier_admission_smoke",
        "olcf_side_storage",
        "long_term_storage",
        "authorization_date",
        "authorized_by",
        "reviewer",
    }
    required_policy_keys = allowed_policy_keys - {
        "authorization_date",
        "authorized_by",
        "reviewer",
    }
    if (
        type(policy.get("schema_version")) is not int
        or policy.get("schema_version") != 1
        or not required_policy_keys <= set(policy)
        or not set(policy) <= allowed_policy_keys
    ):
        raise ValueError("Storage policy root schema is malformed")
    for key in ["authorization_date", "authorized_by", "reviewer"]:
        if key in policy and (
            not isinstance(policy[key], str) or not policy[key].strip()
        ):
            raise ValueError(f"Storage policy root metadata {key} must be non-empty text")
    frontier = policy.get("frontier")
    storage = policy.get("olcf_side_storage")
    long_term = policy.get("long_term_storage")
    science_freeze = policy.get("science_submission_freeze")
    admission_smoke = policy.get("frontier_admission_smoke")
    registered_science_slices = policy.get("registered_science_slices")
    if (
        not isinstance(frontier, dict)
        or not isinstance(storage, dict)
        or not isinstance(long_term, dict)
        or not isinstance(science_freeze, dict)
        or not isinstance(admission_smoke, dict)
        or not isinstance(registered_science_slices, list)
    ):
        raise ValueError(
            "Storage policy is missing Frontier, OLCF-side storage, science-freeze, "
            "admission-smoke, or registered-science data"
        )
    required_storage_keys = {
        "installed_control_plane_version",
        "staged_control_plane_candidate_version",
        "installed_control_plane_lifecycle",
        "orion_simulation_root_preflight",
        "project_home_mirror_root",
        "project_home_usage",
        "project_home_retention_role",
        "project_home_ledger_mirror_transport",
        "project_home_preflight",
        "orion_bulk_evidence_root",
        "orion_bulk_evidence_usage",
        "orion_retention_role",
        "manual_accounting_authorizations",
        "ledger_genesis_allowed",
    }
    allowed_storage_keys = required_storage_keys | {
        "status",
        "last_preflight_utc",
        "historical_project_home_bulk_artifacts",
        "ledger_genesis_authorization",
        "ledger_genesis",
    }
    if (
        not required_storage_keys <= set(storage)
        or not set(storage) <= allowed_storage_keys
    ):
        raise ValueError("Storage policy OLCF-side storage schema is malformed")
    if set(long_term) != {"status", "selected_destination", "risk", "blocks"}:
        raise ValueError("Storage policy long-term storage schema is malformed")
    if (
        "status" in storage
        and storage["status"] != AUTHORIZED_OLCF_SIDE_STORAGE_STATUS
    ):
        raise ValueError("Storage policy OLCF-side storage status is not authorized")
    if "last_preflight_utc" in storage:
        utc_datetime(
            storage["last_preflight_utc"],
            field="olcf_side_storage.last_preflight_utc",
        )
    if (
        "historical_project_home_bulk_artifacts" in storage
        and storage["historical_project_home_bulk_artifacts"]
        != AUTHORIZED_HISTORICAL_PROJECT_HOME_BULK_ARTIFACTS
    ):
        raise ValueError("Storage policy historical Project Home artifact role is not authorized")
    if (
        "ledger_genesis_authorization" in storage
        and storage["ledger_genesis_authorization"]
        != AUTHORIZED_LEDGER_GENESIS_AUTHORIZATION
    ):
        raise ValueError("Storage policy ledger-genesis authorization is not authorized")
    expected_frontier = {
        "account": authorized_account,
        "partition": AUTHORIZED_PARTITION,
        "simulation_root": str(authorized_pic_root.resolve()),
    }
    for key, expected in expected_frontier.items():
        if frontier.get(key) != expected:
            raise ValueError(f"Storage policy frontier.{key} must be {expected!r}")
    maximum_node_hours = frontier.get("maximum_node_hours")
    if (
        not isinstance(maximum_node_hours, (int, float))
        or isinstance(maximum_node_hours, bool)
        or float(maximum_node_hours) != AUTHORIZED_NODE_HOUR_CAP
    ):
        raise ValueError(
            f"Storage policy frontier.maximum_node_hours must be {AUTHORIZED_NODE_HOUR_CAP!r}"
        )
    if frontier.get("serial_pic_submissions") is not True:
        raise ValueError("Storage policy frontier.serial_pic_submissions must be True")
    allowed_frontier_keys = {
        *expected_frontier,
        "maximum_node_hours",
        "serial_pic_submissions",
        "qos_policy",
    }
    if not set(frontier) <= allowed_frontier_keys:
        raise ValueError("Storage policy frontier schema is malformed")
    if (
        "qos_policy" in frontier
        and frontier["qos_policy"] != "debug_preferred_normal_fallback"
    ):
        raise ValueError("Storage policy Frontier QoS policy is not authorized")
    project_home_mirror_root = storage.get("project_home_mirror_root")
    if (
        not isinstance(project_home_mirror_root, str)
        or Path(project_home_mirror_root).resolve()
        != authorized_project_home_root.resolve()
    ):
        raise ValueError("Storage policy Project Home mirror root is not authorized")
    if storage.get("project_home_usage") != AUTHORIZED_PROJECT_HOME_USAGE:
        raise ValueError("Storage policy Project Home usage is not authorized")
    if (
        storage.get("project_home_retention_role")
        != AUTHORIZED_PROJECT_HOME_RETENTION_ROLE
    ):
        raise ValueError("Storage policy Project Home retention role is not authorized")
    if ledger_mirror_transport != AUTHORIZED_LEDGER_MIRROR_TRANSPORT:
        raise ValueError("Only filesystem_copy Project Home ledger mirroring is authorized")
    if storage.get("project_home_ledger_mirror_transport") != ledger_mirror_transport:
        raise ValueError("Storage policy does not authorize the ledger mirror transport")
    orion_bulk_evidence_root = storage.get("orion_bulk_evidence_root")
    if (
        not isinstance(orion_bulk_evidence_root, str)
        or Path(orion_bulk_evidence_root).resolve() != authorized_pic_root.resolve()
    ):
        raise ValueError("Storage policy Orion bulk-evidence root is not authorized")
    if storage.get("orion_bulk_evidence_usage") != AUTHORIZED_ORION_BULK_EVIDENCE_USAGE:
        raise ValueError("Storage policy Orion bulk-evidence usage is not authorized")
    if storage.get("orion_retention_role") != AUTHORIZED_ORION_RETENTION_ROLE:
        raise ValueError("Storage policy Orion retention role is not authorized")
    if long_term.get("status") != AUTHORIZED_LONG_TERM_STORAGE_STATUS:
        raise ValueError("Storage policy must preserve the Orion-only durability risk")
    selected_destination = long_term.get("selected_destination")
    if (
        not isinstance(selected_destination, str)
        or Path(selected_destination).resolve() != authorized_pic_root.resolve()
    ):
        raise ValueError("Storage policy long-term destination is not authorized")
    if long_term.get("risk") != AUTHORIZED_LONG_TERM_STORAGE_RISK:
        raise ValueError("Storage policy long-term risk statement is not authorized")
    if long_term.get("blocks") != AUTHORIZED_LONG_TERM_STORAGE_BLOCKS:
        raise ValueError("Storage policy long-term terminal block is not authorized")
    for key, expected_path in [
        ("orion_simulation_root_preflight", authorized_pic_root.resolve()),
        ("project_home_preflight", authorized_project_home_root.resolve()),
    ]:
        record = storage.get(key)
        if (
            not isinstance(record, dict)
            or "status" not in record
            or not set(record) <= {"status", "path", "method"}
            or record.get("status") != "passed"
        ):
            raise ValueError(f"Storage policy {key} has not passed")
        if (
            "path" in record
            and (
                not isinstance(record["path"], str)
                or Path(record["path"]).resolve() != expected_path
            )
        ):
            raise ValueError(f"Storage policy {key}.path is not authorized")
        if (
            "method" in record
            and record["method"] != AUTHORIZED_STORAGE_PREFLIGHT_METHOD
        ):
            raise ValueError(f"Storage policy {key}.method is not authorized")
    manual_authorizations = storage.get("manual_accounting_authorizations")
    if not isinstance(manual_authorizations, list) or len(manual_authorizations) > 16:
        raise ValueError("Storage policy manual-accounting authorizations are malformed")
    manual_authorization_ids: set[str] = set()
    manual_authorization_paths: set[Path] = set()
    manual_authorization_root = (
        Path(os.path.abspath(authorized_pic_root))
        / "policy"
        / "manual_accounting_authorizations"
    )
    for authorization in manual_authorizations:
        if not isinstance(authorization, dict) or set(authorization) != {
            "authorization_id",
            "path",
            "project_home_path",
            "sha256",
        }:
            raise ValueError("Storage policy manual-accounting authorization is malformed")
        authorization_id = authorization.get("authorization_id")
        if (
            not isinstance(authorization_id, str)
            or re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,127}", authorization_id) is None
            or authorization_id in manual_authorization_ids
        ):
            raise ValueError("Storage policy manual-accounting authorization ID is invalid")
        if not isinstance(authorization.get("path"), str):
            raise ValueError("Storage policy manual-accounting authorization path is invalid")
        path = require_canonical_path_below(
            Path(authorization["path"]), manual_authorization_root
        )
        project_home_authorization_root = (
            Path(os.path.abspath(authorized_project_home_root))
            / "policy"
            / "manual_accounting_authorizations"
        )
        if not isinstance(authorization.get("project_home_path"), str):
            raise ValueError(
                "Storage policy manual-accounting authorization mirror path is invalid"
            )
        project_home_path = require_canonical_path_below(
            Path(authorization["project_home_path"]), project_home_authorization_root
        )
        if path.name != f"{authorization_id}.json" or path in manual_authorization_paths:
            raise ValueError("Storage policy manual-accounting authorization path is invalid")
        if project_home_path.name != path.name:
            raise ValueError(
                "Storage policy manual-accounting authorization mirror path is invalid"
            )
        digest = authorization.get("sha256")
        if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            raise ValueError("Storage policy manual-accounting authorization SHA-256 is invalid")
        data = read_stable_regular_file_below(
            path, manual_authorization_root, require_read_only_mode=True
        )
        project_home_data = read_stable_regular_file_below(
            project_home_path,
            project_home_authorization_root,
            require_read_only_mode=True,
        )
        if sha256_bytes(data) != digest or project_home_data != data:
            raise ValueError("Storage policy manual-accounting authorization bytes differ")
        manual_authorization_ids.add(authorization_id)
        manual_authorization_paths.add(path)
    genesis_allowed = storage.get("ledger_genesis_allowed")
    genesis = storage.get("ledger_genesis")
    if genesis_allowed is True:
        if not allow_pending_genesis:
            raise ValueError("Storage policy has not closed Frontier PIC ledger genesis")
        if genesis is not None:
            raise ValueError("Pending ledger genesis must not retain initialized fields")
    elif genesis_allowed is False:
        if not isinstance(genesis, dict) or set(genesis) != {
            "status",
            "timestamp",
            "control_plane_version",
            "event_sha256",
            "mirror_ack_sha256",
            "mirror_transport",
        }:
            raise ValueError("Storage policy initialized ledger genesis is malformed")
        if genesis.get("status") != "initialized":
            raise ValueError("Storage policy ledger genesis is not initialized")
        utc_datetime(genesis.get("timestamp"), field="ledger_genesis.timestamp")
        if not _is_lowercase_sha256(genesis.get("control_plane_version")):
            raise ValueError("Storage policy ledger genesis control-plane version is malformed")
        for field in ["event_sha256", "mirror_ack_sha256"]:
            if not _is_lowercase_sha256(genesis.get(field)):
                raise ValueError(f"Storage policy ledger_genesis.{field} is malformed")
        if genesis.get("mirror_transport") != AUTHORIZED_LEDGER_MIRROR_TRANSPORT:
            raise ValueError("Storage policy ledger genesis mirror transport is invalid")
    else:
        raise ValueError("Storage policy ledger genesis authorization is malformed")
    if storage.get("installed_control_plane_version") != control_plane_version:
        raise ValueError("Storage policy does not authorize this control-plane version")
    if storage.get("staged_control_plane_candidate_version") != control_plane_version:
        raise ValueError("Storage policy staged candidate does not match this control-plane version")
    if storage.get("installed_control_plane_lifecycle") != (
        "paired_installed_reviewed_generation"
    ):
        raise ValueError("Storage policy does not attest a paired installed reviewed generation")
    freeze_status = science_freeze.get("status")
    if freeze_status == PENDING_CLEAN_CANDIDATE_FREEZE:
        if set(science_freeze) != {"status"}:
            raise ValueError("Pending science freeze must not retain candidate fields")
    elif freeze_status == AUTHORIZED_CLEAN_CANDIDATE_FREEZE:
        if set(science_freeze) != {
            "status",
            "manifest_path",
            "manifest_sha256",
            "build_profile_control_plane_version",
        }:
            raise ValueError("Authorized science freeze must identify one exact manifest")
        if not isinstance(science_freeze["manifest_path"], str):
            raise ValueError("Authorized science-freeze manifest path is malformed")
        manifest_path = Path(science_freeze["manifest_path"])
        require_canonical_path_below(
            manifest_path, authorized_pic_root / "clean_candidates"
        )
        if manifest_path.name != "clean_candidate_manifest.json":
            raise ValueError("Authorized science-freeze manifest has an invalid path")
        if not _is_lowercase_sha256(science_freeze["manifest_sha256"]):
            raise ValueError("Authorized science-freeze manifest digest is malformed")
        if not _is_lowercase_sha256(
            science_freeze["build_profile_control_plane_version"]
        ):
            raise ValueError(
                "Authorized science-freeze build-profile control-plane version is malformed"
            )
    else:
        raise ValueError("Storage policy does not declare a recognized science freeze state")
    identifiers = set()
    for registered_slice in registered_science_slices:
        if not isinstance(registered_slice, dict) or set(registered_slice) != {
            "authorization_id",
            "status",
            "campaign",
            "test_id",
            "evidence_class",
            "physical_mode",
            "runtime_profile",
            "selected_qos",
            "registered_short_nonproduction",
            "maximum_nodes",
            "maximum_walltime_seconds",
            "maximum_attempts",
            "job_script_sha256",
            "input_deck_sha256",
            "environment_profile_sha256",
            "analysis_script_sha256",
            "executable_sha256",
            "launch_contract_sha256",
            "clean_candidate_manifest_sha256",
        }:
            raise ValueError("Storage policy registered-science authorization is malformed")
        identifier = registered_slice["authorization_id"]
        if (
            not isinstance(identifier, str)
            or not re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,63}", identifier)
            or identifier in identifiers
        ):
            raise ValueError("Storage policy registered-science authorization ID is invalid")
        identifiers.add(identifier)
        if registered_slice["status"] != AUTHORIZED_REGISTERED_SCIENCE_SLICE_STATUS:
            raise ValueError("Storage policy registered-science slice is not authorized")
        for key in ["campaign", "test_id", "evidence_class", "physical_mode"]:
            if (
                not isinstance(registered_slice[key], str)
                or not registered_slice[key].strip()
            ):
                raise ValueError(f"Storage policy registered-science {key} is empty")
        if registered_slice["runtime_profile"] != "frontier_minimum_supported":
            raise ValueError("Storage policy registered-science runtime profile is invalid")
        if registered_slice["selected_qos"] not in {"debug", "normal"}:
            raise ValueError("Storage policy registered-science QoS is invalid")
        if not isinstance(registered_slice["registered_short_nonproduction"], bool):
            raise ValueError("Storage policy registered-science short-job flag is invalid")
        for key in ["maximum_nodes", "maximum_walltime_seconds", "maximum_attempts"]:
            number = registered_slice[key]
            if not isinstance(number, int) or isinstance(number, bool) or number <= 0:
                raise ValueError(f"Storage policy registered-science {key} is invalid")
        for key in [
            "job_script_sha256",
            "input_deck_sha256",
            "environment_profile_sha256",
            "executable_sha256",
            "launch_contract_sha256",
            "clean_candidate_manifest_sha256",
        ]:
            if not _is_lowercase_sha256(registered_slice[key]):
                raise ValueError(f"Storage policy registered-science {key} is malformed")
        analysis_sha256 = registered_slice["analysis_script_sha256"]
        if (
            not isinstance(analysis_sha256, list)
            or not 1 <= len(analysis_sha256) <= 16
            or any(
                not _is_lowercase_sha256(digest)
                for digest in analysis_sha256
            )
        ):
            raise ValueError(
                "Storage policy registered-science analysis_script_sha256 is malformed"
            )
        if freeze_status != AUTHORIZED_CLEAN_CANDIDATE_FREEZE:
            raise ValueError("Registered-science slices require an authorized clean freeze")
        if (
            registered_slice["clean_candidate_manifest_sha256"]
            != science_freeze["manifest_sha256"]
        ):
            raise ValueError("Registered-science slice belongs to another clean freeze")
    if (
        registered_science_slices
        and admission_smoke.get("status") != CLOSED_ADMISSION_SMOKE_STATUS
    ):
        raise ValueError("Registered-science slices require closed admission smoke")
    if admission_smoke.get("status") == CLOSED_ADMISSION_SMOKE_STATUS:
        if set(admission_smoke) != {"status"}:
            raise ValueError("Closed admission smoke must not retain executable fields")
        return policy
    if admission_smoke.get("status") == PENDING_ADMISSION_SMOKE_STATUS:
        if set(admission_smoke) != {"status"}:
            raise ValueError("Pending admission smoke must not retain executable fields")
        return policy
    if set(admission_smoke) != {
        "status",
        "campaign",
        "test_id",
        "evidence_class",
        "physical_mode",
        "selected_qos",
        "registered_short_nonproduction",
        "maximum_nodes",
        "maximum_walltime_seconds",
        "job_script_sha256",
        "input_deck_sha256",
        "environment_profile_sha256",
        "analysis_script_sha256",
        "executable_sha256",
        "launch_contract_sha256",
    }:
        raise ValueError("Storage policy admission-smoke authorization is malformed")
    expected_admission_smoke = {
        "status": AUTHORIZED_ADMISSION_SMOKE_STATUS,
        "campaign": "f0_hipmpi_smoke",
        "test_id": "pic_parser_contract_guards",
        "evidence_class": "frontier_f0_admission_smoke_candidate",
        "physical_mode": "extended_mhd_pic_parser_contract",
        "selected_qos": "debug",
        "registered_short_nonproduction": True,
        "maximum_nodes": 1,
        "maximum_walltime_seconds": 15 * 60,
    }
    for key, expected in expected_admission_smoke.items():
        value = admission_smoke.get(key)
        if type(value) is not type(expected) or value != expected:
            raise ValueError(f"Storage policy frontier_admission_smoke.{key} is invalid")
    for key in [
        "job_script_sha256",
        "input_deck_sha256",
        "environment_profile_sha256",
        "executable_sha256",
        "launch_contract_sha256",
    ]:
        if not _is_lowercase_sha256(admission_smoke.get(key)):
            raise ValueError(f"Storage policy frontier_admission_smoke.{key} is malformed")
    analysis_sha256 = admission_smoke.get("analysis_script_sha256")
    if (
        not isinstance(analysis_sha256, list)
        or len(analysis_sha256) != 1
        or not _is_lowercase_sha256(analysis_sha256[0])
    ):
        raise ValueError(
            "Storage policy frontier_admission_smoke.analysis_script_sha256 is malformed"
        )
    return policy


def require_storage_policy_unlock_snapshot(
    *,
    control_plane_version: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
    ledger_mirror_transport: str = AUTHORIZED_LEDGER_MIRROR_TRANSPORT,
    allow_pending_genesis: bool = False,
) -> tuple[dict[str, object], dict[str, str]]:
    policy_path = canonical_policy_path(authorized_pic_root)
    mirror_policy_path = canonical_policy_path(authorized_project_home_root)
    promotion_path = active_promotion_path(authorized_pic_root)
    mirror_promotion_path = active_promotion_path(authorized_project_home_root)
    artifacts: dict[Path, bytes] = {}
    for path, root in [
        (policy_path, authorized_pic_root),
        (mirror_policy_path, authorized_project_home_root),
        (promotion_path, authorized_pic_root),
        (mirror_promotion_path, authorized_project_home_root),
    ]:
        lexical_root = Path(os.path.abspath(root))
        lexical_root.resolve(strict=True)
        require_canonical_path_below(path, lexical_root)
        artifacts[path] = read_stable_regular_file_below(
            path, lexical_root, require_read_only_mode=True
        )
    promotion_bytes = artifacts[promotion_path]
    mirror_promotion_bytes = artifacts[mirror_promotion_path]
    if promotion_bytes != mirror_promotion_bytes:
        raise ValueError("Orion and Project Home active-policy promotion bytes differ")
    promotion = read_json_bytes(promotion_bytes, label=str(promotion_path))
    mirror_promotion = read_json_bytes(
        mirror_promotion_bytes, label=str(mirror_promotion_path)
    )
    if promotion != mirror_promotion:
        raise ValueError("Orion and Project Home active-policy promotion records differ")
    policy_bytes = artifacts[policy_path]
    mirror_policy_bytes = artifacts[mirror_policy_path]
    if policy_bytes != mirror_policy_bytes:
        raise ValueError("Orion and Project Home active-policy bytes differ")
    policy_sha256 = sha256_bytes(policy_bytes)
    expected = {
        "schema_version": 1,
        "control_plane_version": control_plane_version,
        "policy_path": str(policy_path),
        "project_home_policy_path": str(mirror_policy_path),
        "policy_sha256": policy_sha256,
    }
    if type(promotion.get("schema_version")) is not int or promotion != expected:
        raise ValueError("Active-policy promotion record is not anchored to this control plane")
    if sha256_bytes(mirror_policy_bytes) != promotion["policy_sha256"]:
        raise ValueError("Project Home active-policy mirror checksum differs")
    policy = validate_storage_policy(
        read_json_bytes(policy_bytes, label=str(policy_path)),
        control_plane_version=control_plane_version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        ledger_mirror_transport=ledger_mirror_transport,
        allow_pending_genesis=allow_pending_genesis,
    )
    return policy, {
        "active_policy_sha256": policy_sha256,
        "active_promotion_sha256": sha256_bytes(promotion_bytes),
    }


def require_storage_policy_unlock(
    *,
    control_plane_version: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
    ledger_mirror_transport: str = AUTHORIZED_LEDGER_MIRROR_TRANSPORT,
    allow_pending_genesis: bool = False,
) -> dict[str, object]:
    policy, _ = require_storage_policy_unlock_snapshot(
        control_plane_version=control_plane_version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        ledger_mirror_transport=ledger_mirror_transport,
        allow_pending_genesis=allow_pending_genesis,
    )
    return policy


def require_ledger_paths(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> None:
    expected = [
        (
            ledger_jsonl,
            Path(os.path.abspath(authorized_pic_root)) / "ledger" / "node_hours.jsonl",
            authorized_pic_root,
        ),
        (
            ledger_csv,
            Path(os.path.abspath(authorized_pic_root)) / "ledger" / "node_hours.csv",
            authorized_pic_root,
        ),
        (
            receipts_jsonl,
            Path(os.path.abspath(authorized_pic_root)) / "ledger" / "mirror_receipts.jsonl",
            authorized_pic_root,
        ),
        (
            mirror_jsonl,
            Path(os.path.abspath(authorized_project_home_root)) / "ledger" / "node_hours.jsonl",
            authorized_project_home_root,
        ),
    ]
    for supplied, required, root in expected:
        actual = Path(os.path.abspath(supplied))
        if actual != required:
            raise ValueError(f"Unauthorized ledger path: {actual}; expected {required}")
        require_canonical_path_below(actual, root)


def stable_serialization_anchor(authorized_pic_root: Path) -> Path:
    lexical_root = Path(os.path.abspath(authorized_pic_root))
    production_root = Path(os.path.abspath(AUTHORIZED_PIC_ROOT))
    if lexical_root == production_root:
        return production_root.parents[2]
    return lexical_root.parent
