#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Shared helpers for immutable Frontier PIC submission snapshots."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import ctypes
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

from operator_attestation import validate_sealed_operator_attestation
from q011_pressure_review_packet_verifier import (
    consume_sealed_pressure_reanalysis_attestation,
    consume_sealed_pressure_reviewer_attestation,
    consume_published_pressure_pilot_review_packet,
    validate_pressure_reanalysis_source_snapshot,
)


AUTHORIZED_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
AUTHORIZED_PROJECT_HOME_ROOT = Path(
    "/autofs/nccs-svm1_proj/ast207/proj-shared/PIC"
)
AUTHORIZED_PROJECT_HOME_LEDGER_ROOT = Path("/ccs/proj/ast207/proj-shared/PIC")
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
AUTHORIZED_STORAGE_PREFLIGHT_OPERATIONS = [
    "create_exclusive",
    "write_all",
    "fsync_file",
    "read_back_exact",
    "unlink",
    "fsync_parent",
    "verify_absent",
]
Q043_REGISTERED_MATRIX_RELATIVE = Path(
    "analysis/q043_registered_execution_raw_oracle_successor_v1/"
    "q043_registered_matrix_qualification.json"
)
Q043_REGISTERED_MATRIX_RECORD_TYPE = (
    "q043_registered_execution_raw_oracle_matrix_qualification"
)
Q043_REGISTERED_MATRIX_STATUS = (
    "complete_registered_Q043_matrix_prerequisite_non_authorizing"
)
Q043_REGISTERED_CASE_STATUS = "registered_case_evidence_admitted_non_authorizing"
Q043_REGISTERED_CAMPAIGN = "q043_registered_execution_raw_oracle_successor_v1"
Q043_REGISTERED_CASE_COUNT = 132
Q023_REGISTERED_MATRIX_RELATIVE = Path(
    "analysis/q023_paper_bell_linear_joverc_registered_successor_v1/"
    "q023_registered_matrix_qualification.json"
)
Q023_REGISTERED_MATRIX_RECORD_TYPE = (
    "q023_registered_execution_linear_matrix_qualification"
)
Q023_REGISTERED_MATRIX_STATUS = (
    "complete_registered_Q023_linear_matrix_pass_non_authorizing"
)
Q023_REGISTERED_CASE_STATUS = (
    "registered_execution_case_admitted_non_authorizing"
)
Q023_REGISTERED_CAMPAIGN_ID = "Q023-PAPER-BELL-LINEAR-JOVERC"
Q023_REGISTERED_POLICY_CAMPAIGN = (
    "q023_paper_bell_linear_joverc_registered_successor_v1"
)
Q023_REGISTERED_CASE_COUNT = 55
# The capture helper authenticates a clean tracked HEAD before and after probing.
# Pin the cycle-free reviewed blob closure here; its common-module digest is
# checked against this executing controller below.
AUTHORIZED_STORAGE_PREFLIGHT_CAPTURE_SOURCE_BLOBS = {
    "entrypoint_sha256": (
        "b6dae64b28dbcc7ce82877ad25d53bc4f0637016c4bd274431c1a4ba947ec94b"
    ),
    "runner_sha256": (
        "bab2de7e25d61fd84c1f5a3ee08ef1563b106daa3210517ee717db1ca3e2f5fe"
    ),
    "schema_sha256": (
        "348b80f6b56fa56da57a4d30932b41784c939f5a2b1bf49c82d3a6acd024ee3f"
    ),
}
AUTHORIZED_HISTORICAL_STORAGE_PREFLIGHT_RETIREMENT_POLICY_SHA256 = (
    "647573e109852da0e343cd1c80033dbbe4ad4d2d672588cdea4e852340bad3d8"
)
AUTHORIZED_HISTORICAL_STORAGE_PREFLIGHT_RETIREMENT_PROMOTION_SHA256 = (
    "073da1d4fe7f2eb054da9f3ec2bf2860f2643c3f44ca88d6183f4592eb0681bf"
)
# Exact launch-prohibited Stage-4 predecessor retained after p0=1 publication.
AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_POLICY_SHA256 = (
    "aeab7e4ef92c7cbbd5b84fcd139c046f2a96f21b4f280fa6dd89c80deca981d1"
)
AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_PROMOTION_SHA256 = (
    "ef11cb301ec4917cd32aaca8af56e4f2d753367682c6ba28fc048c004904613e"
)
AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_CONTROL_PLANE_VERSION = (
    "930a04d1d39c873ea49abfcf500069011f6d5759240a8f5c5b3341a6d243b246"
)
AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_PROBE_ID = (
    "d94e275e-ba45-4bbc-9849-125e90103abe"
)
AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_EVIDENCE_SHA256 = (
    "66b06c0b0d3f9af54c63c117aa9d494dbebdcd529518f0b36a0288d2fac8e518"
)
AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_SOURCE_AUTHENTICATION = {
    "common_sha256": (
        "8f2d0e38c3ae79a69caf87c98490d77f7c5e9a540903bd3917049b752aa83b52"
    ),
    "entrypoint_sha256": (
        "b6dae64b28dbcc7ce82877ad25d53bc4f0637016c4bd274431c1a4ba947ec94b"
    ),
    "git_commit": "67a418c432e2d424aa9e6cf5ed16316ea40fc0a4",
    "runner_sha256": (
        "6053f190ed5bea093537ca5e6aef110212d54861f6726a6fa7294a1716eca2d5"
    ),
    "schema_sha256": (
        "348b80f6b56fa56da57a4d30932b41784c939f5a2b1bf49c82d3a6acd024ee3f"
    ),
    "tracked_clean_head_blobs": True,
}
PREPARED_ARTIFACT_INVENTORY_PATH = (
    "tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json"
)
PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS = (
    "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
)
Q011_PLANNER_RETENTION_ROLE = "q011_section54_deterministic_retained_attempt"
Q011_SECTION54_PHYSICAL_MODE = "paper_mhd_pic_vl2_tsc"
Q011_SECTION54_VARIANTS = (
    (
        "coarse_uniform_dx12",
        (
            "mesh_refinement/refinement=none",
            "mesh_refinement/num_levels=1",
            "problem/ps_enable_curvature_amr=false",
        ),
    ),
    ("three_level_amr_root_dx12_finest_dx3", ()),
    (
        "fine_uniform_dx3",
        (
            "mesh/nx1=16000",
            "mesh/nx2=1040",
            "mesh_refinement/refinement=none",
            "mesh_refinement/num_levels=1",
            "problem/ps_enable_curvature_amr=false",
        ),
    ),
)
Q011_SECTION54_SEEDS = (
    23050101,
    23050102,
    23050103,
    23050104,
    23050105,
    23050106,
    23050107,
    23050108,
)
Q011_SECTION54_SOURCE_BINDING_PATHS = {
    "pressure_selection_receipt": "bindings/human_pressure_selection_receipt.json",
    "clean_candidate_manifest": "bindings/clean_candidate_manifest.json",
    "environment_profile": "bindings/environment_profile.sh",
    "qualifying_preregistration": (
        "bindings/q011_section54_qualifying_campaign_preregistration.json"
    ),
    "restart_preregistration": (
        "bindings/q011_section54_restart_continuation_preregistration.json"
    ),
    "paper_deck": "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
}
Q011_SECTION54_ARCHIVE_SOURCE_PATHS = {
    "environment_profile": (
        "tst/publication/frontier_control_plane/frontier_pic_environment.sh"
    ),
    "qualifying_preregistration": (
        "tst/publication/readiness/"
        "q011_section54_qualifying_campaign_preregistration_successor_v3_2026-06-06.json"
    ),
    "restart_preregistration": (
        "tst/publication/readiness/"
        "q011_section54_restart_continuation_preregistration_successor_2026-06-06.json"
    ),
    "paper_deck": "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
}
Q011_SECTION54_HELPER_SOURCES = (
    "tst/publication/q011_section54_model.py",
    "tst/publication/q011_section54_pressure_pilot_execution.py",
    "tst/publication/q011_section54_pressure_selection.py",
    "tst/publication/q011_section54_historical_pressure_pilot_consumer.py",
    "tst/publication/frontier_control_plane/q011_pressure_review_packet_verifier.py",
    "tst/publication/q011_section54_restart.py",
    "tst/publication/analyze_q011_section54_outputs.py",
    "tst/publication/analyze_q011_section54_campaign.py",
    "tst/publication/analyze_q011_section54_numerical_qualification.py",
    "tst/publication/q011_section54_particles.py",
    "tst/publication/q011_section54_spatial.py",
    "tst/publication/q011_section54_artifacts.py",
    "tst/publication/publish_q011_section54_pressure_pilot_bundle.py",
    "tst/publication/analyze_q011_section54_pressure_pilot.py",
    "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
    "tst/publication/frontier_f1_structured_artifacts.py",
    "tst/publication/q011_section54_attempt_manifest_materializer.py",
    "tst/publication/publish_q011_section54_campaign_attempt.py",
    "tst/publication/immutable_orion_tree.py",
    "tst/publication/pvtk_particles.py",
    "tst/publication/q011_parallel_shock_storage_estimator.py",
    "tst/publication/frontier_control_plane/control_plane_common.py",
    "tst/publication/frontier_control_plane/ledger.py",
    "tst/publication/frontier_control_plane/operator_attestation.py",
    "tst/publication/q011_section54_qualifying_campaign_execution.py",
)


def scheduler_account_matches_authorized(value: object) -> bool:
    """Accept the configured account or Slurm's canonical lowercase spelling."""
    return isinstance(value, str) and value in {
        AUTHORIZED_ACCOUNT,
        AUTHORIZED_ACCOUNT.lower(),
    }


def _is_lowercase_sha256(value: object) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _is_canonical_uuid(value: object) -> bool:
    if not isinstance(value, str):
        return False
    try:
        parsed = uuid.UUID(value)
    except ValueError:
        return False
    return str(parsed) == value


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
ACTIVE_PROMOTION_TRANSACTION_RELATIVE = Path(
    "policy/.active_promotion_transaction.json"
)
ACTIVE_PROMOTION_ROLLBACK_ANCHOR_PREFIXES = (
    ".storage_policy.json.transaction-rollback-",
    ".active_promotion.json.transaction-rollback-",
)
SITE_POLICY_MAX_AGE_SECONDS = 24 * 60 * 60
CONTROL_PLANE_FILES = [
    "clean_candidate.schema.json",
    "control_plane.schema.json",
    "control_plane_common.py",
    "create_clean_candidate_freeze.py",
    "create_pre_submit_manifest.py",
    "frontier_job.sh",
    "frontier_pic_environment.sh",
    "initialize_frontier_ledger.py",
    "launch_trampoline.py",
    "launch_with_frontier_profile.sh",
    "ledger.py",
    "operator_attestation.py",
    "promote_active_policy.py",
    "q011_pressure_review_packet_verifier.py",
    "reconcile_frontier_job.py",
    "reconcile_manual_frontier_allocations.py",
    "reconcile_q023_registered_execution.py",
    "reconcile_q019_registered_execution.py",
    "reconcile_q043_registered_execution.py",
    "revalidate_clean_candidate.py",
    "run_installed_control_plane_job.sh",
    "run_control_plane.py",
    "storage_preflight.schema.json",
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


def _exact_json_equal(left: object, right: object) -> bool:
    """Compare JSON values without Python's bool/int/float equality aliases."""
    return json.dumps(
        left,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ) == json.dumps(
        right,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


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


def project_home_ledger_root(authorized_project_home_root: Path) -> Path:
    """Keep the production append-only ledger on its frozen lexical spelling."""
    lexical = Path(os.path.abspath(authorized_project_home_root))
    if lexical == Path(os.path.abspath(AUTHORIZED_PROJECT_HOME_ROOT)):
        return Path(os.path.abspath(AUTHORIZED_PROJECT_HOME_LEDGER_ROOT))
    return lexical


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
            "-DAthena_SINGLE_PRECISION=OFF",
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
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ValueError(f"Artifact is not a regular file: {path}")
        if require_read_only_mode and before.st_mode & 0o222:
            raise ValueError(f"Artifact is not read-only: {path}")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            data = stream.read()
        after = os.fstat(descriptor)
        stable_fields = (
            "st_dev",
            "st_ino",
            "st_mode",
            "st_size",
            "st_mtime_ns",
            "st_ctime_ns",
        )
        if (
            any(
                getattr(before, field) != getattr(after, field)
                for field in stable_fields
            )
            or len(data) != after.st_size
        ):
            raise ValueError(f"Artifact changed while reading: {path}")
        return data
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
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ValueError(f"Artifact is not a regular file: {path}")
        if require_read_only_mode and before.st_mode & 0o222:
            raise ValueError(f"Artifact is not read-only: {path}")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            data = stream.read()
        after = os.fstat(descriptor)
        stable_fields = (
            "st_dev",
            "st_ino",
            "st_mode",
            "st_size",
            "st_mtime_ns",
            "st_ctime_ns",
        )
        if (
            any(
                getattr(before, field) != getattr(after, field)
                for field in stable_fields
            )
            or len(data) != after.st_size
        ):
            raise ValueError(f"Artifact changed while reading: {path}")
        return data
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
        [
            *(
                path
                for path in source_files
                if path.startswith("inputs/tests/pic")
                and path.endswith(".athinput")
                and "/" not in path[len("inputs/tests/") :]
            ),
            *(
                path
                for path in source_files
                if path.startswith(
                    "inputs/tests/"
                    "q043_bell_current_volume_aware_deposited_current_oracle/"
                )
                and path.endswith(".athinput")
                and "/"
                not in path[
                    len(
                        "inputs/tests/"
                        "q043_bell_current_volume_aware_deposited_current_oracle/"
                    ) :
                ]
            ),
            *(
                path
                for path in source_files
                if path.startswith(
                    "inputs/publication/"
                    "q019_physics_first_nonlinear_bell_successor_v2/"
                )
                and path.endswith(".athinput")
                and "/"
                not in path[
                    len(
                        "inputs/publication/"
                        "q019_physics_first_nonlinear_bell_successor_v2/"
                    ) :
                ]
            ),
            *PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS,
        ]
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
            "Prepared paper-deck inventory must exactly cover archived "
            "inputs/tests/pic*.athinput, the nested Q043 and Q019 matrices, "
            "and required publication decks"
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

    @property
    def descriptors(self) -> tuple[int, ...]:
        if not self._descriptors:
            raise ValueError("Pinned directory ancestry is closed")
        return tuple(descriptor for _, descriptor in self._descriptors)

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


def active_promotion_transaction_path(root: Path) -> Path:
    return Path(os.path.abspath(root)) / ACTIVE_PROMOTION_TRANSACTION_RELATIVE


def active_promotion_recovery_entry_names(root: Path) -> list[str]:
    """List names that reserve the active-policy transaction namespace."""
    lexical_root = Path(os.path.abspath(root))
    policy_parent = active_promotion_transaction_path(lexical_root).parent
    try:
        descriptor = open_directory_below(policy_parent, root=lexical_root)
    except FileNotFoundError:
        return []
    try:
        marker_name = active_promotion_transaction_path(lexical_root).name
        names = sorted(
            name
            for name in os.listdir(descriptor)
            if name == marker_name
            or any(
                name.startswith(prefix)
                for prefix in ACTIVE_PROMOTION_ROLLBACK_ANCHOR_PREFIXES
            )
        )
        require_same_directory(policy_parent, descriptor, root=lexical_root)
        return names
    finally:
        os.close(descriptor)


def require_no_active_promotion_transaction(
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> None:
    """Fail closed while a durable active-policy rollback transaction exists."""
    for root in [authorized_pic_root, authorized_project_home_root]:
        if active_promotion_recovery_entry_names(root):
            raise ValueError(
                "Active-policy promotion transaction requires locked recovery"
            )


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


def _clean_candidate_fixed_layout_name(name: str, *, label: str) -> None:
    if (
        not name
        or name in {".", ".."}
        or "/" in name
        or Path(name).name != name
    ):
        raise ValueError(f"{label} has an invalid fixed-layout name")


def _clean_candidate_stable_metadata(metadata: os.stat_result) -> tuple[int, ...]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_mode,
        metadata.st_nlink,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )


def _require_clean_candidate_regular_metadata(
    metadata: os.stat_result, *, label: str
) -> None:
    if not stat.S_ISREG(metadata.st_mode):
        raise ValueError(f"{label} is not a regular file")
    if metadata.st_nlink != 1:
        raise ValueError(f"{label} does not have exactly one hard link")
    if metadata.st_mode & 0o222:
        raise ValueError(f"{label} is not read-only")


def _read_clean_candidate_regular_descriptor(
    descriptor: int, *, label: str
) -> tuple[bytes, tuple[int, ...]]:
    before = os.fstat(descriptor)
    _require_clean_candidate_regular_metadata(before, label=label)
    payload = bytearray()
    while chunk := os.read(descriptor, 1024 * 1024):
        payload.extend(chunk)
    after = os.fstat(descriptor)
    if (
        _clean_candidate_stable_metadata(before)
        != _clean_candidate_stable_metadata(after)
        or len(payload) != after.st_size
    ):
        raise ValueError(f"{label} changed while reading")
    return bytes(payload), _clean_candidate_stable_metadata(after)


def _read_clean_candidate_regular_file_at(
    directory_descriptor: int, name: str, *, label: str
) -> bytes:
    """Read one fixed-layout candidate member through a retained parent descriptor."""
    _clean_candidate_fixed_layout_name(name, label=label)
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_descriptor,
    )
    try:
        payload, _ = _read_clean_candidate_regular_descriptor(
            descriptor, label=label
        )
        return payload
    finally:
        os.close(descriptor)


class _RetainedCleanCandidateClosure:
    """Retain and recheck every descriptor opened for one candidate closure."""

    def __init__(
        self, *, read_regular_file_at: Callable[..., bytes] | None = None
    ) -> None:
        self._read_regular_file_at = read_regular_file_at
        self._directories: list[
            tuple[int, str, int, tuple[int, ...], str]
        ] = []
        self._files: list[tuple[int, str, int, tuple[int, ...], str]] = []
        self._entries: dict[int, tuple[set[str], str]] = {}
        self._watch_descriptor = self._open_watch_descriptor()

    def _open_watch_descriptor(self) -> int:
        libc = ctypes.CDLL(None, use_errno=True)
        try:
            inotify_init1 = libc.inotify_init1
        except AttributeError as error:
            raise OSError("Clean-candidate retained closure requires Linux inotify") from error
        inotify_init1.argtypes = [ctypes.c_int]
        inotify_init1.restype = ctypes.c_int
        descriptor = inotify_init1(os.O_NONBLOCK | os.O_CLOEXEC)
        if descriptor < 0:
            number = ctypes.get_errno()
            raise OSError(number, os.strerror(number))
        return descriptor

    def _watch(self, descriptor: int) -> None:
        libc = ctypes.CDLL(None, use_errno=True)
        inotify_add_watch = libc.inotify_add_watch
        inotify_add_watch.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_uint32]
        inotify_add_watch.restype = ctypes.c_int
        mask = (
            0x00000002  # IN_MODIFY
            | 0x00000004  # IN_ATTRIB
            | 0x00000008  # IN_CLOSE_WRITE
            | 0x00000040  # IN_MOVED_FROM
            | 0x00000080  # IN_MOVED_TO
            | 0x00000100  # IN_CREATE
            | 0x00000200  # IN_DELETE
            | 0x00000400  # IN_DELETE_SELF
            | 0x00000800  # IN_MOVE_SELF
            | 0x00002000  # IN_UNMOUNT
        )
        path = os.fsencode(f"/proc/{os.getpid()}/fd/{descriptor}")
        if inotify_add_watch(self._watch_descriptor, path, mask) < 0:
            number = ctypes.get_errno()
            raise OSError(number, os.strerror(number), os.fsdecode(path))

    def _require_no_watch_events(self) -> None:
        try:
            payload = os.read(self._watch_descriptor, 1024 * 1024)
        except BlockingIOError:
            return
        if payload:
            raise ValueError("Clean-candidate tree changed during retained closure")

    def watch_ancestry(self, descriptors: Iterable[int]) -> None:
        """Watch retained root ancestry so rename-away/restore cannot evade identity checks."""
        self._require_no_watch_events()
        for descriptor in descriptors:
            self._watch(descriptor)
        self._require_no_watch_events()

    def open_directory_at(
        self, directory_descriptor: int, name: str, *, label: str
    ) -> int:
        self._require_no_watch_events()
        _clean_candidate_fixed_layout_name(name, label=label)
        descriptor = os.open(
            name,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_descriptor,
        )
        metadata = os.fstat(descriptor)
        if not stat.S_ISDIR(metadata.st_mode) or metadata.st_mode & 0o222:
            os.close(descriptor)
            raise ValueError(f"{label} is not a read-only directory")
        try:
            self._watch(descriptor)
        except BaseException:
            os.close(descriptor)
            raise
        self._directories.append(
            (
                directory_descriptor,
                name,
                descriptor,
                _clean_candidate_stable_metadata(metadata),
                label,
            )
        )
        return descriptor

    def require_exact_entries(
        self, directory_descriptor: int, expected: set[str], *, label: str
    ) -> None:
        self._require_no_watch_events()
        prior = self._entries.get(directory_descriptor)
        if prior is not None and prior != (expected, label):
            raise ValueError(f"{label} fixed layout changed during capture")
        self._entries[directory_descriptor] = (set(expected), label)
        self._require_exact_entries_now(directory_descriptor, expected, label=label)
        self._require_no_watch_events()

    def _require_exact_entries_now(
        self, directory_descriptor: int, expected: set[str], *, label: str
    ) -> None:
        if set(os.listdir(directory_descriptor)) != expected:
            raise ValueError(f"{label} entries do not match the fixed layout")

    def read_regular_file_at(
        self, directory_descriptor: int, name: str, *, label: str
    ) -> bytes:
        self._require_no_watch_events()
        _clean_candidate_fixed_layout_name(name, label=label)
        descriptor = os.open(
            name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_descriptor,
        )
        try:
            _require_clean_candidate_regular_metadata(
                os.fstat(descriptor), label=label
            )
            self._watch(descriptor)
            injected_payload: bytes | None = None
            if (
                self._read_regular_file_at is not None
                and self._read_regular_file_at
                is not _read_clean_candidate_regular_file_at
            ):
                injected_payload = self._read_regular_file_at(
                    directory_descriptor, name, label=label
                )
                if not isinstance(injected_payload, bytes):
                    raise ValueError(f"{label} injected reader returned malformed bytes")
            payload, metadata = _read_clean_candidate_regular_descriptor(
                descriptor, label=label
            )
            if injected_payload is not None and injected_payload != payload:
                raise ValueError(f"{label} changed while reading")
            self._require_no_watch_events()
            self._files.append(
                (directory_descriptor, name, descriptor, metadata, label)
            )
            descriptor = -1
            return payload
        finally:
            if descriptor >= 0:
                os.close(descriptor)

    def _require_directory_same(
        self,
        parent_descriptor: int,
        name: str,
        descriptor: int,
        expected: tuple[int, ...],
        *,
        label: str,
    ) -> None:
        retained = os.fstat(descriptor)
        if (
            not stat.S_ISDIR(retained.st_mode)
            or retained.st_mode & 0o222
            or _clean_candidate_stable_metadata(retained) != expected
        ):
            raise ValueError(f"{label} changed during retained closure")
        lexical_descriptor = os.open(
            name,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_descriptor,
        )
        try:
            if _clean_candidate_stable_metadata(
                os.fstat(lexical_descriptor)
            ) != expected:
                raise ValueError(f"{label} changed during retained closure")
        finally:
            os.close(lexical_descriptor)

    def _require_file_same(
        self,
        parent_descriptor: int,
        name: str,
        descriptor: int,
        expected: tuple[int, ...],
        *,
        label: str,
    ) -> None:
        retained = os.fstat(descriptor)
        if (
            not stat.S_ISREG(retained.st_mode)
            or retained.st_nlink != 1
            or retained.st_mode & 0o222
            or _clean_candidate_stable_metadata(retained) != expected
        ):
            raise ValueError(f"{label} changed during retained closure")
        lexical_descriptor = os.open(
            name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_descriptor,
        )
        try:
            lexical = os.fstat(lexical_descriptor)
            if (
                not stat.S_ISREG(lexical.st_mode)
                or lexical.st_nlink != 1
                or lexical.st_mode & 0o222
                or _clean_candidate_stable_metadata(lexical) != expected
            ):
                raise ValueError(f"{label} changed during retained closure")
        finally:
            os.close(lexical_descriptor)

    def require_same(self) -> None:
        self._require_no_watch_events()
        for parent, name, descriptor, metadata, label in self._directories:
            self._require_directory_same(
                parent, name, descriptor, metadata, label=label
            )
        for descriptor, (expected, label) in self._entries.items():
            self._require_exact_entries_now(descriptor, expected, label=label)
        for parent, name, descriptor, metadata, label in self._files:
            self._require_file_same(
                parent, name, descriptor, metadata, label=label
            )
        for parent, name, descriptor, metadata, label in self._directories:
            self._require_directory_same(
                parent, name, descriptor, metadata, label=label
            )
        for descriptor, (expected, label) in self._entries.items():
            self._require_exact_entries_now(descriptor, expected, label=label)
        self._require_no_watch_events()

    def close(self, *, close_ancestry: Callable[[], None] | None = None) -> None:
        error: BaseException | None = None
        for _, _, descriptor, _, _ in reversed(self._files):
            try:
                os.close(descriptor)
            except BaseException as caught:
                error = error or caught
        self._files.clear()
        for _, _, descriptor, _, _ in reversed(self._directories):
            try:
                os.close(descriptor)
            except BaseException as caught:
                error = error or caught
        self._directories.clear()
        self._entries.clear()
        if close_ancestry is not None:
            try:
                close_ancestry()
            except BaseException as caught:
                error = error or caught
        try:
            self._require_no_watch_events()
        except BaseException as caught:
            error = error or caught
        try:
            os.close(self._watch_descriptor)
        except BaseException as caught:
            error = error or caught
        if error is not None:
            raise error


def _clean_candidate_mapping(
    record: dict[str, object], key: str
) -> dict[str, object]:
    value = record.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Clean-candidate manifest has no {key} object")
    return value


def _clean_candidate_text(record: dict[str, object], key: str) -> str:
    value = record.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Clean-candidate manifest has no {key}")
    return value.strip()


def _require_exact_clean_candidate_layout_path(
    record: dict[str, object], key: str, expected: Path
) -> None:
    if _clean_candidate_text(record, key) != str(expected):
        raise ValueError(f"Clean-candidate {key} does not match the fixed layout")


def _read_clean_candidate_submodules(
    source: dict[str, object],
    *,
    candidate_dir: Path,
    candidate_descriptor: int,
    closure: _RetainedCleanCandidateClosure,
) -> tuple[list[bytes], list[bytes]]:
    records = source.get("submodules")
    if not isinstance(records, list):
        raise ValueError("Clean-candidate submodules must be a list")
    if not records:
        return [], []
    submodules_descriptor = closure.open_directory_at(
        candidate_descriptor,
        "submodules",
        label="Clean-candidate submodule directory",
    )
    expected: set[str] = set()
    for index in range(len(records)):
        expected.add(f"{index:04d}.tar")
        expected.add(f"{index:04d}.commit")
    closure.require_exact_entries(
        submodules_descriptor,
        expected,
        label="Clean-candidate submodule directory",
    )
    archives: list[bytes] = []
    commits: list[bytes] = []
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            raise ValueError(
                "Clean-candidate submodule attestation must be an object"
            )
        archive_name = f"{index:04d}.tar"
        commit_name = f"{index:04d}.commit"
        _require_exact_clean_candidate_layout_path(
            record,
            "archive_path",
            candidate_dir / "submodules" / archive_name,
        )
        archives.append(
            closure.read_regular_file_at(
                submodules_descriptor,
                archive_name,
                label=f"Clean-candidate submodule archive {index}",
            )
        )
        _require_exact_clean_candidate_layout_path(
            record,
            "commit_path",
            candidate_dir / "submodules" / commit_name,
        )
        commits.append(
            closure.read_regular_file_at(
                submodules_descriptor,
                commit_name,
                label=f"Clean-candidate submodule commit object {index}",
            )
        )
    closure.require_exact_entries(
        submodules_descriptor,
        expected,
        label="Clean-candidate submodule directory",
    )
    return archives, commits


def read_clean_candidate_tree(
    candidate_manifest_path: Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    read_regular_file_at: Callable[..., bytes] | None = None,
) -> dict[str, object]:
    """Capture one fixed one-level clean-candidate closure without path traversal."""
    lexical_pic_root = Path(os.path.abspath(authorized_pic_root))
    candidate_root = lexical_pic_root / "clean_candidates"
    candidate_path = require_canonical_path_below(
        Path(os.path.abspath(candidate_manifest_path)), candidate_root
    )
    if (
        candidate_path.name != "clean_candidate_manifest.json"
        or candidate_path.parent.parent != candidate_root
    ):
        raise ValueError(
            "Clean-candidate manifest must be one fixed level below clean_candidates"
        )
    closure = _RetainedCleanCandidateClosure(
        read_regular_file_at=read_regular_file_at
    )
    try:
        candidate_root_ancestry = PinnedDirectoryAncestry(
            candidate_root, root=lexical_pic_root
        )
    except BaseException:
        closure.close()
        raise
    with candidate_root_ancestry:
        try:
            closure.watch_ancestry(candidate_root_ancestry.descriptors)
            candidate_root_descriptor = candidate_root_ancestry.descriptor
            candidate_descriptor = closure.open_directory_at(
                candidate_root_descriptor,
                candidate_path.parent.name,
                label="Clean-candidate directory",
            )
        except BaseException:
            closure.close(close_ancestry=candidate_root_ancestry.close)
            raise
        try:
            candidate_bytes = closure.read_regular_file_at(
                candidate_descriptor,
                "clean_candidate_manifest.json",
                label="Clean-candidate manifest",
            )
            candidate = read_json_bytes(
                candidate_bytes, label="Clean-candidate manifest"
            )
            freeze_id = _clean_candidate_text(candidate, "freeze_id")
            try:
                parsed_freeze_id = uuid.UUID(freeze_id)
            except ValueError as error:
                raise ValueError("Clean-candidate freeze ID is malformed") from error
            if str(parsed_freeze_id) != freeze_id:
                raise ValueError("Clean-candidate freeze ID is not canonical")
            if candidate_path.parent.name != freeze_id:
                raise ValueError("Clean-candidate path does not match its freeze ID")
            utc_datetime(
                candidate.get("created_utc"), field="clean_candidate.created_utc"
            )
            source = _clean_candidate_mapping(candidate, "source")
            build = _clean_candidate_mapping(candidate, "build")
            _require_exact_clean_candidate_layout_path(
                source, "archive_path", candidate_path.parent / "source.tar"
            )
            _require_exact_clean_candidate_layout_path(
                source, "commit_path", candidate_path.parent / "source.commit"
            )
            _require_exact_clean_candidate_layout_path(
                build, "profile_path", candidate_path.parent / "build_profile.json"
            )
            _require_exact_clean_candidate_layout_path(
                build,
                "profile_receipt_path",
                candidate_path.parent / "profile_receipt.json",
            )
            _require_exact_clean_candidate_layout_path(
                build, "executable_path", candidate_path.parent / "athena"
            )
            expected = {
                "athena",
                "build_provenance",
                "build_profile.json",
                "clean_candidate_manifest.json",
                "profile_receipt.json",
                "source.commit",
                "source.tar",
            }
            if source.get("submodules"):
                expected.add("submodules")
            closure.require_exact_entries(
                candidate_descriptor,
                expected,
                label="Clean-candidate directory",
            )
            source_archive = closure.read_regular_file_at(
                candidate_descriptor,
                "source.tar",
                label="Clean-candidate source archive",
            )
            source_commit = closure.read_regular_file_at(
                candidate_descriptor,
                "source.commit",
                label="Clean-candidate source commit object",
            )
            submodule_archives, submodule_commits = _read_clean_candidate_submodules(
                source,
                candidate_dir=candidate_path.parent,
                candidate_descriptor=candidate_descriptor,
                closure=closure,
            )
            build_profile = closure.read_regular_file_at(
                candidate_descriptor,
                "build_profile.json",
                label="Clean-candidate build profile",
            )
            build_profile_receipt = closure.read_regular_file_at(
                candidate_descriptor,
                "profile_receipt.json",
                label="Clean-candidate build-profile receipt",
            )
            executable = closure.read_regular_file_at(
                candidate_descriptor,
                "athena",
                label="Clean-candidate executable",
            )
            provenance_descriptor = closure.open_directory_at(
                candidate_descriptor,
                "build_provenance",
                label="Frozen build provenance directory",
            )
            provenance_expected = set(BUILD_PROVENANCE_FILENAMES.values())
            closure.require_exact_entries(
                provenance_descriptor,
                provenance_expected,
                label="Frozen build provenance directory",
            )
            build_provenance = {
                label: closure.read_regular_file_at(
                    provenance_descriptor,
                    filename,
                    label=f"Frozen build provenance {label}",
                )
                for label, filename in BUILD_PROVENANCE_FILENAMES.items()
            }
            closure.require_exact_entries(
                provenance_descriptor,
                provenance_expected,
                label="Frozen build provenance directory",
            )
            closure.require_exact_entries(
                candidate_descriptor,
                expected,
                label="Clean-candidate directory",
            )
            closure.require_same()
            candidate_root_ancestry.require_same()
            # The ancestry recheck itself is inside the watched lifetime. Drain
            # again so rename-away-and-restore cannot hide in its return window.
            closure.require_same()
            return {
                "candidate_manifest_path": candidate_path,
                "candidate_manifest_bytes": candidate_bytes,
                "candidate": candidate,
                "source_archive": source_archive,
                "source_commit": source_commit,
                "submodule_archives": submodule_archives,
                "submodule_commits": submodule_commits,
                "build_profile": build_profile,
                "build_profile_receipt": build_profile_receipt,
                "build_provenance": build_provenance,
                "executable": executable,
            }
        finally:
            closure.close(close_ancestry=candidate_root_ancestry.close)


def _planner_relative_path(value: object, *, label: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{label} must be text")
    path = PurePosixPath(value)
    if (
        not value
        or path.is_absolute()
        or value != path.as_posix()
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise ValueError(f"{label} must be one canonical relative path")
    return value


def _planner_binding(value: object, *, label: str) -> dict[str, str]:
    if not isinstance(value, dict) or set(value) != {"path", "sha256"}:
        raise ValueError(f"{label} schema is malformed")
    relative = _planner_relative_path(value["path"], label=f"{label} path")
    digest = value["sha256"]
    if not _is_lowercase_sha256(digest):
        raise ValueError(f"{label} digest is malformed")
    return {"path": relative, "sha256": digest}


def _planner_json(payload: bytes, *, label: str) -> dict[str, object]:
    def reject_duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, item in pairs:
            if key in result:
                raise ValueError(f"{label} repeats JSON key {key!r}")
            result[key] = item
        return result

    try:
        value = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=reject_duplicates,
            parse_constant=lambda item: (_ for _ in ()).throw(
                ValueError(f"{label} contains forbidden JSON constant {item}")
            ),
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not canonical UTF-8 JSON") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value


def _planner_regular_bytes_at(parent_descriptor: int, name: str, *, label: str) -> bytes:
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=parent_descriptor,
    )
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_mode & 0o222
        ):
            raise ValueError(f"{label} is not one read-only regular file")
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        current = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
        identity = lambda item: (  # noqa: E731
            item.st_dev,
            item.st_ino,
            item.st_mode,
            item.st_nlink,
            item.st_size,
            item.st_mtime_ns,
            item.st_ctime_ns,
        )
        if (
            identity(before) != identity(after)
            or (after.st_dev, after.st_ino) != (current.st_dev, current.st_ino)
        ):
            raise ValueError(f"{label} changed while reading")
        return bytes(payload)
    finally:
        os.close(descriptor)


def _scan_read_only_planner_tree(
    descriptor: int, *, prefix: str = ""
) -> tuple[dict[str, bytes], set[str]]:
    metadata = os.fstat(descriptor)
    if not stat.S_ISDIR(metadata.st_mode) or metadata.st_mode & 0o222:
        raise ValueError("Immutable qualifying planner contains a writable directory")
    files: dict[str, bytes] = {}
    directories: set[str] = set()
    for name in sorted(os.listdir(descriptor)):
        if name in {"", ".", ".."} or "/" in name:
            raise ValueError("Immutable qualifying planner contains an unsafe member")
        relative = f"{prefix}/{name}" if prefix else name
        observed = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
        if stat.S_ISDIR(observed.st_mode):
            child = os.open(
                name,
                os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=descriptor,
            )
            try:
                opened = os.fstat(child)
                if (opened.st_dev, opened.st_ino) != (observed.st_dev, observed.st_ino):
                    raise ValueError("Immutable qualifying planner directory changed")
                child_files, child_directories = _scan_read_only_planner_tree(
                    child, prefix=relative
                )
                current = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
                if (current.st_dev, current.st_ino) != (opened.st_dev, opened.st_ino):
                    raise ValueError("Immutable qualifying planner directory changed")
                files.update(child_files)
                directories.update(child_directories)
                directories.add(relative)
            finally:
                os.close(child)
        elif stat.S_ISREG(observed.st_mode):
            files[relative] = _planner_regular_bytes_at(
                descriptor, name, label=f"immutable qualifying planner member {relative}"
            )
        else:
            raise ValueError("Immutable qualifying planner contains a non-file member")
    return files, directories


def _planner_required_directories(paths: Iterable[str]) -> set[str]:
    directories: set[str] = set()
    for relative in paths:
        parent = PurePosixPath(relative).parent
        while parent.as_posix() != ".":
            directories.add(parent.as_posix())
            parent = parent.parent
    return directories


def _planner_inventory(payload: bytes) -> dict[str, str]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("Immutable qualifying planner inventory is not UTF-8") from error
    records: dict[str, str] = {}
    for line in text.splitlines(keepends=True):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)\n", line)
        if match is None:
            raise ValueError("Immutable qualifying planner inventory is malformed")
        relative = _planner_relative_path(match.group(2), label="planner inventory member")
        if relative == "artifact_inventory.sha256" or relative in records:
            raise ValueError("Immutable qualifying planner inventory is malformed")
        records[relative] = match.group(1)
    if text != "".join(f"{records[path]}  {path}\n" for path in sorted(records)):
        raise ValueError("Immutable qualifying planner inventory is noncanonical")
    return records


def _planner_member(
    files: dict[str, bytes], binding: object, *, label: str
) -> tuple[dict[str, str], bytes]:
    normalized = _planner_binding(binding, label=label)
    try:
        payload = files[normalized["path"]]
    except KeyError as error:
        raise ValueError(f"{label} is absent from immutable qualifying planner") from error
    if hashlib.sha256(payload).hexdigest() != normalized["sha256"]:
        raise ValueError(f"{label} checksum drifted")
    return normalized, payload


def _planner_digest_value(value: object) -> str:
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    ).hexdigest()


def _planner_archive_member(archive_payload: bytes, relative: str) -> bytes:
    """Read one regular source-archive member without accepting aliases."""
    try:
        with tarfile.open(fileobj=io.BytesIO(archive_payload), mode="r:") as archive:
            members = [member for member in archive.getmembers() if member.name == relative]
            if len(members) != 1 or not members[0].isfile():
                raise ValueError(
                    f"Clean-candidate source archive omits reviewed member: {relative}"
                )
            stream = archive.extractfile(members[0])
            if stream is None:
                raise ValueError(
                    f"Clean-candidate source archive member is unreadable: {relative}"
                )
            return stream.read()
    except tarfile.TarError as error:
        raise ValueError("Clean-candidate source archive is not a readable tar file") from error


def _planner_exact_primitive_types(actual: object, expected: object, *, label: str) -> None:
    """Reject JSON boolean aliases and other primitive-type substitutions."""
    if type(actual) is not type(expected):
        raise ValueError(f"{label} primitive type drifted")
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise ValueError(f"{label} keys drifted")
        for key in sorted(expected):
            _planner_exact_primitive_types(actual[key], expected[key], label=f"{label}/{key}")
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError(f"{label} length drifted")
        for index, (actual_item, expected_item) in enumerate(zip(actual, expected)):
            _planner_exact_primitive_types(
                actual_item, expected_item, label=f"{label}[{index}]"
            )


def _planner_expected_matrix() -> dict[str, object]:
    return {
        "physical_mode": Q011_SECTION54_PHYSICAL_MODE,
        "grid_variants": [variant for variant, _ in Q011_SECTION54_VARIANTS],
        "qualifying_seeds": list(Q011_SECTION54_SEEDS),
        "expected_baseline_attempts": 24,
        "paired_seed_rule": (
            "Use the same qualifying seed for coarse-uniform, AMR and fine-uniform variants."
        ),
    }


def _planner_attempt_id(index: int, variant: str, seed: int) -> str:
    return f"baseline-{index:03d}-{variant}-seed-{seed}"


def _planner_source_archive_and_candidate(
    files: dict[str, bytes],
    normalized_source_bindings: dict[str, dict[str, str]],
    candidate: object,
    *,
    pic_root: Path,
) -> tuple[bytes, dict[str, object]]:
    """Derive the planner candidate solely from its bound frozen candidate bytes."""
    if not isinstance(candidate, dict):
        raise ValueError("Planner candidate binding must be an object")
    manifest_binding = normalized_source_bindings["clean_candidate_manifest"]
    manifest_payload = files[manifest_binding["path"]]
    manifest = _planner_json(manifest_payload, label="planner clean-candidate manifest")
    if (
        set(manifest) != {
            "schema_version",
            "freeze_id",
            "created_utc",
            "prepared_artifacts",
            "source",
            "build",
        }
        or type(manifest.get("schema_version")) is not int
        or manifest["schema_version"] != 4
    ):
        raise ValueError("Planner clean-candidate manifest schema drifted")
    freeze_id = str(manifest.get("freeze_id", ""))
    try:
        uuid.UUID(freeze_id)
    except ValueError as error:
        raise ValueError("Planner clean-candidate freeze ID is malformed") from error
    candidate_root = pic_root / "clean_candidates" / freeze_id
    external_manifest = candidate_root / "clean_candidate_manifest.json"
    external_manifest_payload = read_stable_regular_file_below(
        external_manifest, pic_root, require_read_only_mode=True
    )
    if external_manifest_payload != manifest_payload:
        raise ValueError("Planner clean-candidate manifest differs from frozen candidate")
    source = manifest.get("source")
    build = manifest.get("build")
    prepared = manifest.get("prepared_artifacts")
    if not isinstance(source, dict) or not isinstance(build, dict) or not isinstance(prepared, dict):
        raise ValueError("Planner clean-candidate manifest objects drifted")
    source_archive_path = candidate_root / "source.tar"
    if source.get("archive_path") != str(source_archive_path):
        raise ValueError("Planner clean-candidate source archive path drifted")
    source_archive = read_stable_regular_file_below(
        source_archive_path, pic_root, require_read_only_mode=True
    )
    if hashlib.sha256(source_archive).hexdigest() != source.get("archive_sha256"):
        raise ValueError("Planner clean-candidate source archive checksum drifted")
    submodules = source.get("submodules")
    if not isinstance(submodules, list):
        raise ValueError("Planner clean-candidate submodules drifted")
    validated_submodules = []
    for record in submodules:
        if not isinstance(record, dict):
            raise ValueError("Planner clean-candidate submodule binding drifted")
        validated_submodules.append(
            {
                "path": record["path"],
                "archive_sha256": record["archive_sha256"],
                "commit_sha256": record["commit_sha256"],
                "git_commit": record["git_commit"],
                "git_tree": record["git_tree"],
            }
        )
    environment_payload = files[normalized_source_bindings["environment_profile"]["path"]]
    environment_digest = hashlib.sha256(environment_payload).hexdigest()
    reviewed_environment = _planner_archive_member(
        source_archive, Q011_SECTION54_ARCHIVE_SOURCE_PATHS["environment_profile"]
    )
    if environment_payload != reviewed_environment:
        raise ValueError("Planner environment profile differs from reviewed source archive")
    expected_candidate = {
        "clean_candidate_manifest": {
            "path": str(external_manifest),
            "sha256": hashlib.sha256(manifest_payload).hexdigest(),
        },
        "freeze_id": freeze_id,
        "git_commit": source["git_commit"],
        "git_tree": source["git_tree"],
        "source_archive_sha256": source["archive_sha256"],
        "source_commit_sha256": source["commit_sha256"],
        "source_bundle_sha256": source["source_bundle_sha256"],
        "prepared_artifact_inventory_sha256": prepared["inventory_sha256"],
        "validated_submodules": validated_submodules,
        "build_profile": {
            "path": build["profile_path"],
            "sha256": build["profile_sha256"],
        },
        "build_profile_receipt": {
            "path": build["profile_receipt_path"],
            "sha256": build["profile_receipt_sha256"],
        },
        "build_invocations_sha256": build["build_invocations_sha256"],
        "executable": {
            "path": build["executable_path"],
            "sha256": build["executable_sha256"],
        },
        "environment_profile": {
            "path": str(
                pic_root
                / "control_plane"
                / candidate["environment_profile"]["control_plane_version"]
                / "frontier_pic_environment.sh"
            ),
            "sha256": environment_digest,
            "control_plane_version": candidate["environment_profile"][
                "control_plane_version"
            ],
            "reviewed_source": {
                "path": Q011_SECTION54_ARCHIVE_SOURCE_PATHS["environment_profile"],
                "sha256": environment_digest,
            },
        },
    }
    _planner_exact_primitive_types(candidate, expected_candidate, label="planner candidate binding")
    if candidate != expected_candidate:
        raise ValueError("Planner candidate binding drifted from frozen candidate")
    return source_archive, expected_candidate


def _planner_expected_baseline_contract(
    *,
    attempt_id: str,
    variant: str,
    model_overrides: tuple[str, ...],
    seed: int,
    selected_ps_p0: float,
    candidate: dict[str, object],
    paper_deck_binding: dict[str, str],
    attempt_root: Path,
) -> dict[str, object]:
    return {
        "record_type": "q011_section54_launch_prohibited_handoff_contract",
        "schema_version": 1,
        "contract_role": "source_local_review_handoff_only",
        "launch_authorized": False,
        "scheduler_submission_authorized": False,
        "live_policy_mutation_authorized": False,
        "attempt_id": attempt_id,
        "variant": variant,
        "qualifying_seed": seed,
        "selected_problem_ps_p0": selected_ps_p0,
        "executable": candidate["executable"],
        "environment_profile": candidate["environment_profile"],
        "paper_deck": paper_deck_binding,
        "authorized_orion_attempt_root": str(attempt_root),
        "argv": [
            "-i",
            Q011_SECTION54_SOURCE_BINDING_PATHS["paper_deck"],
            "-d",
            str(attempt_root / "raw"),
            f"job/basename={attempt_id}",
            f"problem/ps_p0={selected_ps_p0!r}",
            f"particles/pic_random_seed={seed}",
            f"problem/ps_inject_seed={seed}",
            f"problem/ps_seed_noise_seed={seed}",
            *model_overrides,
        ],
        "required_separate_boundary": (
            "review_and_promote_a_registered_frontier_submission_policy_then_use_"
            "the_installed_control_plane_wrapper"
        ),
    }


def _planner_expected_baseline_descriptor(
    *,
    index: int,
    attempt_id: str,
    variant: str,
    seed: int,
    selected_ps_p0: float,
    candidate: dict[str, object],
    attempt_root: Path,
    contract_path: str,
    contract_payload: bytes,
) -> dict[str, object]:
    return {
        "record_type": "q011_section54_baseline_attempt_descriptor",
        "schema_version": 1,
        "attempt_index": index,
        "attempt_id": attempt_id,
        "status": "planned_not_authorized",
        "variant": variant,
        "qualifying_seed": seed,
        "physical_mode": Q011_SECTION54_PHYSICAL_MODE,
        "selected_problem_ps_p0": selected_ps_p0,
        "candidate_binding": {
            "git_commit": candidate["git_commit"],
            "source_bundle_sha256": candidate["source_bundle_sha256"],
            "executable_sha256": candidate["executable"]["sha256"],
            "environment_profile_sha256": candidate["environment_profile"]["sha256"],
        },
        "authorized_orion_attempt_root": str(attempt_root),
        "launch_contract": {
            "path": contract_path,
            "sha256": hashlib.sha256(contract_payload).hexdigest(),
        },
        "artifact_retention": (
            "retain_every_emitted_raw_artifact_for_every_attempt_including_failed_attempts"
        ),
    }


def _planner_restart_carrier_id(seed: int) -> str:
    return f"amr-restart-continuation-seed-{seed}"


def _planner_expected_restart_contract(
    *,
    carrier_id: str,
    source_attempt: dict[str, object],
    restart_preregistration: dict[str, object],
    candidate: dict[str, object],
    paper_deck_binding: dict[str, str],
    attempt_root: Path,
) -> dict[str, object]:
    continuation = restart_preregistration["continuation_contract"]
    checkpoint = continuation["checkpoint_nominal_slot_omega0_inverse"]
    return {
        "record_type": "q011_section54_launch_prohibited_handoff_contract",
        "schema_version": 1,
        "contract_role": "source_local_restart_review_handoff_only",
        "launch_authorized": False,
        "scheduler_submission_authorized": False,
        "live_policy_mutation_authorized": False,
        "carrier_id": carrier_id,
        "source_baseline_attempt_id": source_attempt["attempt_id"],
        "variant": source_attempt["variant"],
        "qualifying_seed": source_attempt["qualifying_seed"],
        "selected_problem_ps_p0": source_attempt["selected_problem_ps_p0"],
        "executable": candidate["executable"],
        "environment_profile": candidate["environment_profile"],
        "paper_deck": paper_deck_binding,
        "authorized_orion_attempt_root": str(attempt_root),
        "checkpoint_nominal_slot_omega0_inverse": checkpoint,
        "checkpoint_observed_commit_binding_required": True,
        "checkpoint_input": (
            f"retain_from_{source_attempt['attempt_id']}_at_nominal_slot_t{int(checkpoint)}_"
            "then_bind_exact_observed_cycle_time_and_checksum_before_any_separately_"
            "authorized_continuation"
        ),
        "argv_template": [
            "-r",
            "<exact-retained-checkpoint-path-bound-by-separate-reviewed-successor>",
            "-d",
            str(attempt_root / "raw"),
        ],
        "required_separate_boundary": (
            "materialize_an_observed_cycle_time_and_checkpoint_checksum_bound_"
            "registered_restart_successor_then_review_and_promote_policy_before_"
            "using_the_installed_control_plane_wrapper"
        ),
    }


def _planner_expected_restart_carrier(
    *,
    carrier_id: str,
    source_attempt: dict[str, object],
    restart_preregistration: dict[str, object],
    restart_preregistration_binding: dict[str, str],
    attempt_root: Path,
    contract_path: str,
    contract_payload: bytes,
) -> dict[str, object]:
    continuation = restart_preregistration["continuation_contract"]
    return {
        "record_type": "q011_section54_amr_restart_continuation_carrier",
        "schema_version": 1,
        "carrier_id": carrier_id,
        "status": "planned_not_authorized",
        "source_baseline_attempt_id": source_attempt["attempt_id"],
        "variant": source_attempt["variant"],
        "qualifying_seed": source_attempt["qualifying_seed"],
        "selected_problem_ps_p0": source_attempt["selected_problem_ps_p0"],
        "authorized_orion_attempt_root": str(attempt_root),
        "restart_preregistration": restart_preregistration_binding,
        "checkpoint_nominal_slot_omega0_inverse": continuation[
            "checkpoint_nominal_slot_omega0_inverse"
        ],
        "checkpoint_observed_commit_binding_required": True,
        "retained_output_nominal_slots_after_checkpoint_omega0_inverse": continuation[
            "retained_output_nominal_slots_after_checkpoint_omega0_inverse"
        ],
        "retained_output_pairing_policy": continuation[
            "retained_output_pairing_policy"
        ],
        "comparison_tolerances_max_absolute_difference": continuation[
            "comparison_tolerances_max_absolute_difference"
        ],
        "launch_contract": {
            "path": contract_path,
            "sha256": hashlib.sha256(contract_payload).hexdigest(),
        },
    }


def _planner_expected_independent_recompute_plan(
    *,
    plan_id: str,
    campaign_root: Path,
    qualifying_preregistration_binding: dict[str, str],
) -> dict[str, object]:
    return {
        "record_type": "q011_section54_independent_raw_artifact_recompute_plan",
        "schema_version": 1,
        "plan_id": plan_id,
        "status": "plan_frozen_independent_implementation_and_review_artifacts_open",
        "qualifying_preregistration": qualifying_preregistration_binding,
        "authorized_orion_campaign_root": str(campaign_root),
        "implementation_rule": (
            "Recompute every primary metric from archived raw artifacts using a "
            "reviewer-owned script or an independently implemented analyzer. The "
            "recompute script must not import production metric-extraction functions "
            "or the local analyzer and helper sources bound by this plan."
        ),
        "raw_input_policy": (
            "consume_root_relative_sha256_inventories_and_archived_raw_bin_pvtk_rst_"
            "bytes_for_every_attempt_including_failed_attempts"
        ),
        "required_records": [
            "independent_script_path_or_archive_locator",
            "independent_script_sha256",
            "independent_environment_lock",
            "attempt_inventory_sha256_values",
            "input_artifact_checksums",
            "metric_comparison_table",
            "reviewer_identity",
            "reviewer_disposition",
        ],
        "required_metric_table_columns": [
            "attempt_id",
            "variant",
            "qualifying_seed",
            "observable",
            "production_metric",
            "independent_metric",
            "absolute_difference",
            "relative_difference",
            "declared_tolerance",
            "disposition",
        ],
        "production_helper_imports_authorized": False,
        "claim_closure_authorized": False,
        "frontier_execution_authorized": False,
    }


def _planner_expected_policy_fragment(
    *,
    plan_id: str,
    pic_root: Path,
    campaign_root: Path,
    candidate: dict[str, object],
    pressure_receipt_binding: dict[str, str],
    contract_bindings: list[dict[str, str]],
    restart_contract_binding: dict[str, str],
) -> dict[str, object]:
    return {
        "record_type": "q011_section54_qualifying_campaign_nonauthorizing_policy_fragment",
        "schema_version": 1,
        "plan_id": plan_id,
        "status": "review_fragment_only_not_live_policy",
        "qualification_effect": "plan_only_no_execution_authorization_no_claim_closure",
        "authorized_orion_root": str(pic_root),
        "authorized_orion_campaign_root": str(campaign_root),
        "selected_pressure_receipt": pressure_receipt_binding,
        "candidate_binding": candidate,
        "baseline_launch_contracts": contract_bindings,
        "restart_continuation_launch_contract": restart_contract_binding,
        "integration_policy": (
            "requires_separate_reviewed_registered_frontier_submission_policy_"
            "successor_and_installed_control_plane_promotion"
        ),
        "mutates_live_policy": False,
        "scheduler_calls_authorized": False,
        "scheduler_submission_authorized": False,
        "frontier_execution_authorized": False,
        "launch_authorized": False,
        "claim_closure_authorized": False,
    }


def _planner_normalized_pressure_receipt(
    value: object, *, pic_root: Path
) -> dict[str, object]:
    """Require the canonical source-local human pressure-selection receipt shape."""
    if (
        not isinstance(value, dict)
        or set(value)
        != {
            "schema_version",
            "record_type",
            "selection_method",
            "published_pressure_pilot_receipt",
            "published_pressure_pilot_review_packet_receipt",
            "pilot_bundle_manifest_sha256",
            "aggregate_pilot_analysis_sha256",
            "case_descriptors",
            "selected_case",
            "authoritative_reanalysis_attestation",
            "reviewer_attestation",
        }
        or type(value.get("schema_version")) is not int
        or value["schema_version"] != 3
        or value.get("record_type") != "q011_section54_pressure_selection_receipt"
        or value.get("selection_method") != "human_review_only"
    ):
        raise ValueError("Planner human pressure-selection receipt schema drifted")
    published = value.get("published_pressure_pilot_receipt")
    packet = value.get("published_pressure_pilot_review_packet_receipt")
    if (
        not isinstance(published, dict)
        or set(published) != {"path", "sha256"}
        or not isinstance(published.get("path"), str)
        or not Path(published["path"]).is_absolute()
        or not _is_lowercase_sha256(published.get("sha256"))
        or not isinstance(packet, dict)
        or set(packet) != {"path", "sha256"}
        or not isinstance(packet.get("path"), str)
        or not Path(packet["path"]).is_absolute()
        or not _is_lowercase_sha256(packet.get("sha256"))
        or not _is_lowercase_sha256(value.get("pilot_bundle_manifest_sha256"))
        or not _is_lowercase_sha256(value.get("aggregate_pilot_analysis_sha256"))
    ):
        raise ValueError("Planner human pressure-selection publication binding drifted")
    published_path = Path(published["path"])
    packet_path = Path(packet["path"])
    if (
        Path(os.path.abspath(published_path)) != published_path
        or published_path.parent != pic_root / "publication"
        or Path(os.path.abspath(packet_path)) != packet_path
        or packet_path.parent != pic_root / "publication"
    ):
        raise ValueError("Planner human pressure-selection publication path drifted")
    try:
        packet_verification = consume_published_pressure_pilot_review_packet(
            packet["path"],
            aggregate_receipt_binding=published,
            authorized_pic_root=pic_root,
        )
    except ValueError as error:
        raise ValueError(
            f"Planner human pressure-selection review packet verification failed: {error}"
        ) from error
    if (
        not isinstance(packet_verification, dict)
        or set(packet_verification)
        != {
            "receipt_binding",
            "aggregate_receipt_binding",
            "packet_receipt",
            "aggregate_receipt",
            "aggregate_bundle",
            "aggregate_analysis",
            "source_bindings",
            "inventory",
        }
        or packet_verification["receipt_binding"] != packet
        or packet_verification["aggregate_receipt_binding"] != published
    ):
        raise ValueError(
            "Planner human pressure-selection review packet verifier result drifted"
        )
    aggregate_receipt = packet_verification["aggregate_receipt"]
    if (
        not isinstance(aggregate_receipt, dict)
        or packet_verification["aggregate_bundle"]
        != aggregate_receipt.get("aggregate_bundle")
        or packet_verification["aggregate_analysis"]
        != aggregate_receipt.get("aggregate_analysis")
    ):
        raise ValueError(
            "Planner human pressure-selection aggregate verifier result drifted"
        )
    registered = (
        ("ps_p0_1p00", 1.0),
        ("ps_p0_0p05", 0.05),
        ("ps_p0_0p10", 0.1),
        ("ps_p0_0p20", 0.2),
    )
    descriptors = value.get("case_descriptors")
    if not isinstance(descriptors, list) or len(descriptors) != len(registered):
        raise ValueError("Planner human pressure-selection descriptor set drifted")
    descriptor_digests = []
    for descriptor, (case_id, problem_ps_p0) in zip(descriptors, registered):
        if (
            not isinstance(descriptor, dict)
            or set(descriptor) != {"case_id", "problem_ps_p0", "descriptor_sha256"}
            or descriptor.get("case_id") != case_id
            or type(descriptor.get("problem_ps_p0")) is not float
            or descriptor["problem_ps_p0"] != problem_ps_p0
            or not _is_lowercase_sha256(descriptor.get("descriptor_sha256"))
        ):
            raise ValueError("Planner human pressure-selection descriptor set drifted")
        descriptor_digests.append(descriptor["descriptor_sha256"])
    if len(set(descriptor_digests)) != len(descriptor_digests):
        raise ValueError("Planner human pressure-selection descriptors are not unique")
    aggregate_bundle = aggregate_receipt.get("aggregate_bundle")
    aggregate_analysis = aggregate_receipt.get("aggregate_analysis")
    aggregate_cases = aggregate_receipt.get("raw_cases")
    if (
        not isinstance(aggregate_bundle, dict)
        or set(aggregate_bundle) != {"path", "manifest_sha256"}
        or not isinstance(aggregate_analysis, dict)
        or set(aggregate_analysis) != {"path", "sha256"}
        or not isinstance(aggregate_cases, list)
        or len(aggregate_cases) != len(descriptors)
        or value["pilot_bundle_manifest_sha256"]
        != aggregate_bundle.get("manifest_sha256")
        or value["aggregate_pilot_analysis_sha256"]
        != aggregate_analysis.get("sha256")
    ):
        raise ValueError(
            "Planner human pressure-selection aggregate evidence binding drifted"
        )
    for descriptor, aggregate_case in zip(descriptors, aggregate_cases):
        if (
            not isinstance(aggregate_case, dict)
            or set(aggregate_case)
            != {
                "case_id",
                "artifact_dir",
                "descriptor_path",
                "descriptor_sha256",
                "artifact_inventory_sha256",
                "runtime_artifacts",
            }
            or descriptor["case_id"] != aggregate_case.get("case_id")
            or descriptor["descriptor_sha256"]
            != aggregate_case.get("descriptor_sha256")
        ):
            raise ValueError(
                "Planner human pressure-selection aggregate descriptor binding drifted"
            )
    selected = value.get("selected_case")
    if (
        not isinstance(selected, dict)
        or set(selected) != {"case_id", "problem_ps_p0"}
        or type(selected.get("problem_ps_p0")) is not float
        or (selected.get("case_id"), selected.get("problem_ps_p0")) not in registered
    ):
        raise ValueError("Planner human pressure-selection review record drifted")
    try:
        reanalysis = consume_sealed_pressure_reanalysis_attestation(
            value["authoritative_reanalysis_attestation"],
            aggregate_receipt_binding=published,
            packet_receipt_binding=packet,
            pilot_bundle_manifest_sha256=value["pilot_bundle_manifest_sha256"],
            aggregate_pilot_analysis_sha256=value["aggregate_pilot_analysis_sha256"],
            authorized_pic_root=pic_root,
        )
        reviewer = consume_sealed_pressure_reviewer_attestation(
            value["reviewer_attestation"],
            aggregate_receipt_binding=published,
            packet_receipt_binding=packet,
            reanalysis_verification=reanalysis,
            selected_case=selected,
            authorized_pic_root=pic_root,
        )
    except ValueError as error:
        raise ValueError(
            f"Planner human pressure-selection attestation verification failed: {error}"
        ) from error
    if (
        reanalysis.get("binding") != value["authoritative_reanalysis_attestation"]
        or reviewer.get("binding") != value["reviewer_attestation"]
    ):
        raise ValueError("Planner human pressure-selection attestation binding drifted")
    return value


def validate_planner_retention_binding(
    value: object,
    *,
    authorized_pic_root: Path,
    expected_clean_candidate_manifest_sha256: str | None = None,
) -> dict[str, object]:
    """Re-derive one Q011 retention overlay from its immutable qualifying plan."""
    expected_keys = {
        "schema_version",
        "retention_role",
        "planner_root",
        "planner_inventory_sha256",
        "planner_plan_id",
        "planner_materialization_receipt",
        "attempt_id",
        "authorized_orion_attempt_root",
        "authorized_orion_raw_root",
        "argv",
    }
    if not isinstance(value, dict) or set(value) != expected_keys:
        raise ValueError("Planner-retention binding schema is malformed")
    attempt_id = value.get("attempt_id")
    plan_id = value.get("planner_plan_id")
    inventory_sha256 = value.get("planner_inventory_sha256")
    if (
        type(value.get("schema_version")) is not int
        or value["schema_version"] != 1
        or value.get("retention_role") != Q011_PLANNER_RETENTION_ROLE
        or not isinstance(attempt_id, str)
        or re.fullmatch(r"[a-z0-9][a-z0-9._-]{0,127}", attempt_id) is None
        or not _is_lowercase_sha256(plan_id)
        or not _is_lowercase_sha256(inventory_sha256)
    ):
        raise ValueError("Planner-retention binding identity is malformed")
    pic_root = Path(os.path.abspath(authorized_pic_root))
    planner_root = Path(str(value.get("planner_root", "")))
    if (
        not planner_root.is_absolute()
        or Path(os.path.abspath(planner_root)) != planner_root
        or planner_root.parent != pic_root / "plans"
        or planner_root.name != f"q011-section54-qualifying-campaign-plan-{plan_id}"
    ):
        raise ValueError("Planner-retention immutable planner root is malformed")
    require_no_symlink_components_below(planner_root, pic_root)
    with PinnedDirectoryAncestry(planner_root, root=pic_root) as ancestry:
        files, directories = _scan_read_only_planner_tree(ancestry.descriptor)
        ancestry.require_same()
    try:
        inventory_payload = files["artifact_inventory.sha256"]
    except KeyError as error:
        raise ValueError("Immutable qualifying planner inventory is absent") from error
    if hashlib.sha256(inventory_payload).hexdigest() != inventory_sha256:
        raise ValueError("Immutable qualifying planner inventory checksum drifted")
    inventory = _planner_inventory(inventory_payload)
    measured = {
        path: hashlib.sha256(payload).hexdigest()
        for path, payload in files.items()
        if path != "artifact_inventory.sha256"
    }
    if measured != inventory or directories != _planner_required_directories(inventory):
        raise ValueError("Immutable qualifying planner tree closure drifted")
    supplied_receipt = _planner_binding(
        value["planner_materialization_receipt"],
        label="planner materialization receipt",
    )
    if supplied_receipt["path"] != "materialization_receipt.json":
        raise ValueError("Planner materialization receipt path drifted")
    receipt_binding, receipt_payload = _planner_member(
        files, supplied_receipt, label="planner materialization receipt"
    )
    receipt = _planner_json(receipt_payload, label="planner materialization receipt")
    if (
        set(receipt)
        != {
            "record_type",
            "schema_version",
            "plan_id",
            "campaign_plan",
            "helper_source_closure",
            "tree_inventory",
        }
        or receipt.get("record_type")
        != "q011_section54_qualifying_campaign_plan_materialization_receipt"
        or type(receipt.get("schema_version")) is not int
        or receipt["schema_version"] != 1
        or receipt.get("plan_id") != plan_id
    ):
        raise ValueError("Planner materialization receipt identity drifted")
    tree_inventory = receipt.get("tree_inventory")
    materialized_inventory_algorithm = (
        "sha256 of '<file_sha256>  <root-relative-path>\\n' entries ordered "
        "lexically by root-relative path"
    )
    if (
        not isinstance(tree_inventory, dict)
        or set(tree_inventory)
        != {
            "algorithm",
            "scope",
            "excludes",
            "sha256",
            "inventoried_file_count",
        }
        or tree_inventory.get("algorithm") != materialized_inventory_algorithm
        or tree_inventory.get("scope")
        != (
            "all materialized campaign-plan members before this receipt "
            "and recursive-freeze metadata"
        )
        or tree_inventory.get("excludes")
        != [
            "materialization_receipt.json",
            "freeze_receipt.json",
            "artifact_inventory.sha256",
        ]
        or not _is_lowercase_sha256(tree_inventory.get("sha256"))
        or type(tree_inventory.get("inventoried_file_count")) is not int
    ):
        raise ValueError("Planner materialization tree inventory drifted")
    pre_receipt_files = {
        path: payload
        for path, payload in files.items()
        if path
        not in {
            "materialization_receipt.json",
            "freeze_receipt.json",
            "artifact_inventory.sha256",
        }
    }
    pre_receipt_inventory = "".join(
        f"{hashlib.sha256(pre_receipt_files[path]).hexdigest()}  {path}\n"
        for path in sorted(pre_receipt_files)
    ).encode("utf-8")
    if (
        hashlib.sha256(pre_receipt_inventory).hexdigest()
        != tree_inventory["sha256"]
        or len(pre_receipt_files) != tree_inventory["inventoried_file_count"]
    ):
        raise ValueError("Planner materialization tree inventory checksum drifted")
    freeze_receipt = _planner_json(
        files.get("freeze_receipt.json", b""), label="planner freeze receipt"
    )
    if (
        freeze_receipt
        != {
            "schema_version": 1,
            "artifact_role": (
                "q011_section54_source_local_immutable_qualifying_campaign_plan"
            ),
            "qualification_effect": (
                "plan_only_no_execution_authorization_no_claim_closure"
            ),
            "inventory_excludes": "artifact_inventory.sha256",
            "freeze_policy": "remove all owner, group and other write bits recursively",
        }
        or type(freeze_receipt.get("schema_version")) is not int
    ):
        raise ValueError("Planner recursive-freeze receipt drifted")
    campaign_plan_binding, campaign_plan_payload = _planner_member(
        files, receipt["campaign_plan"], label="planner campaign plan"
    )
    if campaign_plan_binding["path"] != "campaign_plan.json":
        raise ValueError("Planner campaign-plan path drifted")
    plan = _planner_json(campaign_plan_payload, label="planner campaign plan")
    if (
        set(plan)
        != {
            "record_type",
            "schema_version",
            "plan_id",
            "artifact_role",
            "qualification_effect",
            "status",
            "authorized_orion_root",
            "authorized_orion_campaign_root",
            "selected_pressure",
            "candidate_binding",
            "source_bindings",
            "helper_source_closure",
            "campaign_matrix",
            "baseline_attempt_count",
            "baseline_attempt_descriptors",
            "restart_continuation_carrier",
            "independent_raw_artifact_recompute_plan",
            "nonauthorizing_policy_fragment",
            "execution_boundary",
            "preregistration_execution_boundary",
        }
        or plan.get("record_type") != "q011_section54_qualifying_campaign_execution_plan"
        or type(plan.get("schema_version")) is not int
        or plan["schema_version"] != 1
        or plan.get("plan_id") != plan_id
        or plan.get("artifact_role")
        != "q011_section54_source_local_immutable_qualifying_campaign_plan"
        or plan.get("qualification_effect")
        != "plan_only_no_execution_authorization_no_claim_closure"
        or plan.get("status") != "source_local_immutable_review_plan_only"
        or plan.get("authorized_orion_root") != str(pic_root)
        or plan.get("authorized_orion_campaign_root")
        != str(pic_root / "campaigns" / f"q011-section54-{plan_id}")
    ):
        raise ValueError("Planner campaign-plan identity drifted")
    source_bindings = plan.get("source_bindings")
    if not isinstance(source_bindings, dict):
        raise ValueError("Planner campaign-plan source bindings drifted")
    normalized_source_bindings = {
        name: _planner_binding(binding, label=f"planner source binding {name}")
        for name, binding in source_bindings.items()
    }
    if set(normalized_source_bindings) != {
        "pressure_selection_receipt",
        "clean_candidate_manifest",
        "environment_profile",
        "qualifying_preregistration",
        "restart_preregistration",
        "paper_deck",
    }:
        raise ValueError("Planner campaign-plan source-binding schema drifted")
    if {
        name: binding["path"] for name, binding in normalized_source_bindings.items()
    } != Q011_SECTION54_SOURCE_BINDING_PATHS:
        raise ValueError("Planner campaign-plan source-binding paths drifted")
    if expected_clean_candidate_manifest_sha256 is not None and (
        not _is_lowercase_sha256(expected_clean_candidate_manifest_sha256)
        or normalized_source_bindings["clean_candidate_manifest"]["sha256"]
        != expected_clean_candidate_manifest_sha256
    ):
        raise ValueError(
            "Planner clean-candidate manifest differs from submission binding"
        )
    for name, binding in normalized_source_bindings.items():
        _planner_member(files, binding, label=f"planner source binding {name}")
    source_archive, candidate = _planner_source_archive_and_candidate(
        files,
        normalized_source_bindings,
        plan.get("candidate_binding"),
        pic_root=pic_root,
    )
    for name, archive_relative in Q011_SECTION54_ARCHIVE_SOURCE_PATHS.items():
        retained = files[normalized_source_bindings[name]["path"]]
        if retained != _planner_archive_member(source_archive, archive_relative):
            raise ValueError(f"Planner reviewed source binding drifted: {name}")
    helper_binding, helper_payload = _planner_member(
        files, receipt["helper_source_closure"], label="planner helper-source closure"
    )
    if plan.get("helper_source_closure") != helper_binding:
        raise ValueError("Planner helper-source closure binding drifted")
    helper = _planner_json(helper_payload, label="planner helper-source closure")
    if (
        set(helper) != {"record_type", "schema_version", "plan_id", "sources"}
        or helper.get("record_type") != "q011_section54_helper_source_closure"
        or type(helper.get("schema_version")) is not int
        or helper["schema_version"] != 1
        or helper.get("plan_id") != plan_id
        or not isinstance(helper.get("sources"), list)
    ):
        raise ValueError("Planner helper-source closure identity drifted")
    expected_helper_sources = [
        {
            "path": relative,
            "sha256": hashlib.sha256(
                _planner_archive_member(source_archive, relative)
            ).hexdigest(),
        }
        for relative in Q011_SECTION54_HELPER_SOURCES
    ]
    if helper["sources"] != expected_helper_sources:
        raise ValueError("Planner helper-source closure drifted from reviewed archive bytes")
    matrix = _planner_expected_matrix()
    _planner_exact_primitive_types(
        plan.get("campaign_matrix"), matrix, label="planner campaign matrix"
    )
    if plan["campaign_matrix"] != matrix:
        raise ValueError("Planner campaign matrix drifted from reviewed Section 5.4 matrix")
    pressure_receipt_payload = files[
        normalized_source_bindings["pressure_selection_receipt"]["path"]
    ]
    pressure_receipt = _planner_normalized_pressure_receipt(
        _planner_json(
            pressure_receipt_payload,
            label="planner human pressure-selection receipt",
        ),
        pic_root=pic_root,
    )
    if pressure_receipt_payload != (
        json.dumps(pressure_receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8"):
        raise ValueError("Planner human pressure-selection receipt is not canonical JSON")
    try:
        pressure_reanalysis = consume_sealed_pressure_reanalysis_attestation(
            pressure_receipt["authoritative_reanalysis_attestation"],
            aggregate_receipt_binding=pressure_receipt[
                "published_pressure_pilot_receipt"
            ],
            packet_receipt_binding=pressure_receipt[
                "published_pressure_pilot_review_packet_receipt"
            ],
            pilot_bundle_manifest_sha256=pressure_receipt[
                "pilot_bundle_manifest_sha256"
            ],
            aggregate_pilot_analysis_sha256=pressure_receipt[
                "aggregate_pilot_analysis_sha256"
            ],
            authorized_pic_root=pic_root,
        )
        validate_pressure_reanalysis_source_snapshot(
            pressure_reanalysis,
            git_commit=candidate["git_commit"],
            source_archive_sha256=candidate["source_archive_sha256"],
            helper_source_closure=expected_helper_sources,
            authorized_pic_root=pic_root,
        )
    except ValueError as error:
        raise ValueError(
            f"Planner pressure reanalysis source-snapshot binding failed: {error}"
        ) from error
    selected_pressure = plan.get("selected_pressure")
    if (
        not isinstance(selected_pressure, dict)
        or set(selected_pressure) != {"selection_method", "selected_case", "receipt"}
        or selected_pressure.get("selection_method") != "human_review_only"
        or selected_pressure.get("receipt")
        != normalized_source_bindings["pressure_selection_receipt"]
        or pressure_receipt.get("selection_method") != "human_review_only"
        or pressure_receipt.get("selected_case") != selected_pressure.get("selected_case")
    ):
        raise ValueError("Planner selected-pressure binding drifted")
    selected_case = selected_pressure["selected_case"]
    if (
        not isinstance(selected_case, dict)
        or set(selected_case) != {"case_id", "problem_ps_p0"}
        or (selected_case.get("case_id"), selected_case.get("problem_ps_p0"))
        not in {
            ("ps_p0_1p00", 1.0),
            ("ps_p0_0p05", 0.05),
            ("ps_p0_0p10", 0.1),
            ("ps_p0_0p20", 0.2),
        }
        or type(selected_case["problem_ps_p0"]) is not float
    ):
        raise ValueError("Planner selected pressure case drifted")
    qualifying = _planner_json(
        files[normalized_source_bindings["qualifying_preregistration"]["path"]],
        label="planner qualifying preregistration",
    )
    if (
        not isinstance(qualifying.get("athenak_selected_release_criteria"), dict)
        or qualifying["athenak_selected_release_criteria"].get("campaign_matrix") != matrix
        or plan.get("preregistration_execution_boundary")
        != qualifying.get("qualifying_execution_bindings")
    ):
        raise ValueError("Planner qualifying preregistration binding drifted")
    expected_execution_boundary = {
        "mutates_live_policy": False,
        "scheduler_calls": False,
        "submits_jobs": False,
        "infers_pressure_selection": False,
        "launch_authorized": False,
        "frontier_execution_authorized": False,
        "claim_closure_authorized": False,
    }
    _planner_exact_primitive_types(
        plan.get("execution_boundary"),
        expected_execution_boundary,
        label="planner execution boundary",
    )
    if plan["execution_boundary"] != expected_execution_boundary:
        raise ValueError("Planner execution boundary drifted")
    basis = {
        "record_type": "q011_section54_qualifying_campaign_execution_plan",
        "schema_version": 1,
        "pressure_selection_receipt_sha256": normalized_source_bindings[
            "pressure_selection_receipt"
        ]["sha256"],
        "selected_case": selected_case,
        "candidate_binding": candidate,
        "source_binding_sha256": {
            name: binding["sha256"]
            for name, binding in normalized_source_bindings.items()
        },
        "helper_source_closure": expected_helper_sources,
        "campaign_matrix": matrix,
        "authorized_orion_root": str(pic_root),
    }
    if _planner_digest_value(basis) != plan_id:
        raise ValueError("Planner campaign-plan ID drifted from source-bound bytes")
    descriptors = plan.get("baseline_attempt_descriptors")
    if (
        not isinstance(descriptors, list)
        or plan.get("baseline_attempt_count") != 24
        or len(descriptors) != 24
    ):
        raise ValueError("Planner baseline-attempt descriptors drifted")
    campaign_root = pic_root / "campaigns" / f"q011-section54-{plan_id}"
    expected_attempts = []
    index = 0
    for variant, model_overrides in Q011_SECTION54_VARIANTS:
        for seed in Q011_SECTION54_SEEDS:
            index += 1
            expected_attempts.append(
                (index, _planner_attempt_id(index, variant, seed), variant, model_overrides, seed)
            )
    normalized_descriptors = [
        _planner_binding(binding, label="planner baseline-attempt descriptor")
        for binding in descriptors
    ]
    if [binding["path"] for binding in normalized_descriptors] != [
        f"attempts/baseline/{generated_attempt_id}.json"
        for _, generated_attempt_id, _, _, _ in expected_attempts
    ]:
        raise ValueError("Planner baseline descriptor path ordering drifted")
    selected_contract: dict[str, object] | None = None
    selected_attempt_root: Path | None = None
    restart_source_attempt: dict[str, object] | None = None
    contract_bindings: list[dict[str, str]] = []
    for descriptor_binding, (
        attempt_index,
        generated_attempt_id,
        variant,
        model_overrides,
        seed,
    ) in zip(normalized_descriptors, expected_attempts):
        _, descriptor_payload = _planner_member(
            files, descriptor_binding, label="planner baseline-attempt descriptor"
        )
        descriptor = _planner_json(
            descriptor_payload, label="planner baseline-attempt descriptor"
        )
        attempt_root = campaign_root / "baseline" / generated_attempt_id
        contract_path = f"launch_contracts/baseline/{generated_attempt_id}.json"
        contract_binding, contract_payload = _planner_member(
            files, descriptor.get("launch_contract"), label="planner baseline launch contract"
        )
        if contract_binding["path"] != contract_path:
            raise ValueError("Planner baseline launch-contract path drifted")
        contract_bindings.append(contract_binding)
        contract = _planner_json(contract_payload, label="planner baseline launch contract")
        expected_contract = _planner_expected_baseline_contract(
            attempt_id=generated_attempt_id,
            variant=variant,
            model_overrides=model_overrides,
            seed=seed,
            selected_ps_p0=selected_case["problem_ps_p0"],
            candidate=candidate,
            paper_deck_binding=normalized_source_bindings["paper_deck"],
            attempt_root=attempt_root,
        )
        _planner_exact_primitive_types(
            contract, expected_contract, label="planner baseline launch contract"
        )
        if contract != expected_contract:
            raise ValueError("Planner baseline launch contract drifted")
        expected_descriptor = _planner_expected_baseline_descriptor(
            index=attempt_index,
            attempt_id=generated_attempt_id,
            variant=variant,
            seed=seed,
            selected_ps_p0=selected_case["problem_ps_p0"],
            candidate=candidate,
            attempt_root=attempt_root,
            contract_path=contract_path,
            contract_payload=contract_payload,
        )
        _planner_exact_primitive_types(
            descriptor, expected_descriptor, label="planner baseline-attempt descriptor"
        )
        if descriptor != expected_descriptor:
            raise ValueError("Planner baseline-attempt descriptor drifted")
        if (
            variant == "three_level_amr_root_dx12_finest_dx3"
            and seed == Q011_SECTION54_SEEDS[0]
        ):
            restart_source_attempt = expected_descriptor
        if generated_attempt_id == attempt_id:
            if selected_contract is not None:
                raise ValueError("Planner retention selects duplicate immutable descriptors")
            selected_contract = contract
            selected_attempt_root = attempt_root
    if selected_contract is None or selected_attempt_root is None:
        raise ValueError("Planner retention does not select one immutable descriptor")
    if restart_source_attempt is None:
        raise ValueError("Planner restart-continuation source descriptor is absent")
    restart_preregistration = _planner_json(
        files[normalized_source_bindings["restart_preregistration"]["path"]],
        label="planner restart preregistration",
    )
    carrier_id = _planner_restart_carrier_id(
        restart_source_attempt["qualifying_seed"]
    )
    restart_attempt_root = campaign_root / "restart_continuation" / carrier_id
    restart_contract_path = f"launch_contracts/restart_continuation/{carrier_id}.json"
    restart_carrier_binding, restart_carrier_payload = _planner_member(
        files,
        plan.get("restart_continuation_carrier"),
        label="planner restart-continuation carrier",
    )
    if (
        restart_carrier_binding["path"]
        != "restart_continuation/amr_restart_continuation_carrier.json"
    ):
        raise ValueError("Planner restart-continuation carrier path drifted")
    restart_carrier = _planner_json(
        restart_carrier_payload, label="planner restart-continuation carrier"
    )
    restart_contract_binding, restart_contract_payload = _planner_member(
        files,
        restart_carrier.get("launch_contract"),
        label="planner restart-continuation launch contract",
    )
    if restart_contract_binding["path"] != restart_contract_path:
        raise ValueError("Planner restart-continuation launch-contract path drifted")
    restart_contract = _planner_json(
        restart_contract_payload, label="planner restart-continuation launch contract"
    )
    expected_restart_contract = _planner_expected_restart_contract(
        carrier_id=carrier_id,
        source_attempt=restart_source_attempt,
        restart_preregistration=restart_preregistration,
        candidate=candidate,
        paper_deck_binding=normalized_source_bindings["paper_deck"],
        attempt_root=restart_attempt_root,
    )
    _planner_exact_primitive_types(
        restart_contract,
        expected_restart_contract,
        label="planner restart-continuation launch contract",
    )
    if restart_contract != expected_restart_contract:
        raise ValueError("Planner restart-continuation launch contract drifted")
    expected_restart_carrier = _planner_expected_restart_carrier(
        carrier_id=carrier_id,
        source_attempt=restart_source_attempt,
        restart_preregistration=restart_preregistration,
        restart_preregistration_binding=normalized_source_bindings[
            "restart_preregistration"
        ],
        attempt_root=restart_attempt_root,
        contract_path=restart_contract_path,
        contract_payload=restart_contract_payload,
    )
    _planner_exact_primitive_types(
        restart_carrier,
        expected_restart_carrier,
        label="planner restart-continuation carrier",
    )
    if restart_carrier != expected_restart_carrier:
        raise ValueError("Planner restart-continuation carrier drifted")
    recompute_binding, recompute_payload = _planner_member(
        files,
        plan.get("independent_raw_artifact_recompute_plan"),
        label="planner independent raw-artifact recompute plan",
    )
    if recompute_binding["path"] != "independent_raw_artifact_recompute_plan.json":
        raise ValueError("Planner independent raw-artifact recompute-plan path drifted")
    recompute = _planner_json(
        recompute_payload, label="planner independent raw-artifact recompute plan"
    )
    expected_recompute = _planner_expected_independent_recompute_plan(
        plan_id=plan_id,
        campaign_root=campaign_root,
        qualifying_preregistration_binding=normalized_source_bindings[
            "qualifying_preregistration"
        ],
    )
    _planner_exact_primitive_types(
        recompute, expected_recompute, label="planner independent raw-artifact recompute plan"
    )
    if recompute != expected_recompute:
        raise ValueError("Planner independent raw-artifact recompute plan drifted")
    fragment_binding, fragment_payload = _planner_member(
        files,
        plan.get("nonauthorizing_policy_fragment"),
        label="planner nonauthorizing policy fragment",
    )
    if fragment_binding["path"] != "nonauthorizing_policy_fragment.json":
        raise ValueError("Planner nonauthorizing policy-fragment path drifted")
    fragment = _planner_json(fragment_payload, label="planner nonauthorizing policy fragment")
    expected_fragment = _planner_expected_policy_fragment(
        plan_id=plan_id,
        pic_root=pic_root,
        campaign_root=campaign_root,
        candidate=candidate,
        pressure_receipt_binding=normalized_source_bindings[
            "pressure_selection_receipt"
        ],
        contract_bindings=contract_bindings,
        restart_contract_binding=restart_contract_binding,
    )
    _planner_exact_primitive_types(
        fragment, expected_fragment, label="planner nonauthorizing policy fragment"
    )
    if fragment != expected_fragment:
        raise ValueError("Planner nonauthorizing policy fragment drifted")
    attempt_root = selected_attempt_root
    argv = selected_contract["argv"]
    raw_root = attempt_root / "raw"
    derived = {
        "schema_version": 1,
        "retention_role": Q011_PLANNER_RETENTION_ROLE,
        "planner_root": str(planner_root),
        "planner_inventory_sha256": inventory_sha256,
        "planner_plan_id": plan_id,
        "planner_materialization_receipt": receipt_binding,
        "attempt_id": attempt_id,
        "authorized_orion_attempt_root": str(attempt_root),
        "authorized_orion_raw_root": str(raw_root),
        "argv": list(argv),
    }
    if value != derived:
        raise ValueError("Planner-retention overlay differs from immutable planner bytes")
    return derived


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


_STRICT_STORAGE_PREFLIGHT_PROFILE = "strict_authenticated_mirror"
_HISTORICAL_RETIREMENT_STORAGE_PREFLIGHT_PROFILE = "historical_retirement_predecessor"
_EXACT_REVIEWED_STORAGE_PREFLIGHT_PREDECESSOR_PROFILE = (
    "exact_reviewed_storage_preflight_predecessor"
)


def _canonical_storage_preflight_bytes(value: dict[str, object]) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")


def _validate_storage_preflight_artifact(
    payload: bytes,
    *,
    probe_id: str,
    last_preflight_utc: str,
    orion_path: Path,
    project_home_path: Path,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    authorized_source_authentication: dict[str, object] | None = None,
) -> None:
    artifact = read_json_bytes(payload, label="storage-preflight evidence")
    if payload != _canonical_storage_preflight_bytes(artifact):
        raise ValueError("Storage-preflight evidence must use canonical JSON bytes")
    expected_schema_version = (
        2
        if authorized_source_authentication is None
        or "common_sha256" in authorized_source_authentication
        else 1
    )
    if (
        set(artifact)
        != {
            "completed_utc",
            "method",
            "probes",
            "probe_id",
            "publication",
            "record_type",
            "schema_version",
            "source_authentication",
            "started_utc",
            "status",
        }
        or type(artifact.get("schema_version")) is not int
        or artifact["schema_version"] != expected_schema_version
        or artifact.get("record_type") != "frontier_pic_storage_preflight_evidence"
        or artifact.get("probe_id") != probe_id
        or artifact.get("method") != AUTHORIZED_STORAGE_PREFLIGHT_METHOD
        or artifact.get("status") != "passed"
    ):
        raise ValueError("Storage-preflight evidence root schema is malformed")
    started = utc_datetime(
        artifact.get("started_utc"), field="storage_preflight_evidence.started_utc"
    )
    completed = utc_datetime(
        artifact.get("completed_utc"), field="storage_preflight_evidence.completed_utc"
    )
    if started > completed:
        raise ValueError("Storage-preflight evidence completion predates its start")
    if artifact["completed_utc"] != last_preflight_utc:
        raise ValueError(
            "Storage-preflight evidence completion differs from policy last_preflight_utc"
        )
    source = artifact.get("source_authentication")
    expected_source_keys = (
        set(authorized_source_authentication)
        if authorized_source_authentication is not None
        else {
            *AUTHORIZED_STORAGE_PREFLIGHT_CAPTURE_SOURCE_BLOBS,
            "common_sha256",
            "git_commit",
            "tracked_clean_head_blobs",
        }
    )
    if (
        not isinstance(source, dict)
        or set(source) != expected_source_keys
        or source.get("tracked_clean_head_blobs") is not True
        or not isinstance(source.get("git_commit"), str)
        or re.fullmatch(r"[0-9a-f]{40}", source["git_commit"]) is None
    ):
        raise ValueError("Storage-preflight evidence source authentication is not authorized")
    if authorized_source_authentication is not None:
        if source != authorized_source_authentication:
            raise ValueError(
                "Storage-preflight evidence source authentication is not authorized"
            )
    elif any(
        source.get(field) != digest
        for field, digest in AUTHORIZED_STORAGE_PREFLIGHT_CAPTURE_SOURCE_BLOBS.items()
    ) or source.get("common_sha256") != sha256(Path(__file__)):
        raise ValueError("Storage-preflight evidence source authentication is not authorized")
    # git_commit is retained provenance. The reviewed cycle-free tracked-file
    # digests, executing common-module digest and tracked-clean assertion are
    # the authorization boundary.
    if artifact.get("publication") != {
        "orion_path": str(orion_path),
        "project_home_path": str(project_home_path),
    }:
        raise ValueError("Storage-preflight evidence publication binding is malformed")
    probes = artifact.get("probes")
    expected_probes = [
        ("orion_simulation_root", Path(os.path.abspath(authorized_pic_root))),
        (
            "project_home_mirror_root",
            Path(os.path.abspath(authorized_project_home_root)),
        ),
    ]
    if not isinstance(probes, list) or len(probes) != len(expected_probes):
        raise ValueError("Storage-preflight evidence probes are malformed")
    for probe, (role, root) in zip(probes, expected_probes):
        if (
            not isinstance(probe, dict)
            or set(probe)
            != {
                "operations",
                "path",
                "payload_bytes",
                "payload_sha256",
                "role",
                "st_dev",
                "st_ino",
                "status",
            }
            or probe.get("role") != role
            or probe.get("path") != str(root)
            or probe.get("operations") != AUTHORIZED_STORAGE_PREFLIGHT_OPERATIONS
            or type(probe.get("payload_bytes")) is not int
            or probe["payload_bytes"] != 32
            or not _is_lowercase_sha256(probe.get("payload_sha256"))
            or type(probe.get("st_dev")) is not int
            or probe["st_dev"] < 0
            or type(probe.get("st_ino")) is not int
            or probe["st_ino"] < 1
            or probe.get("status") != "passed"
        ):
            raise ValueError("Storage-preflight evidence probe is malformed")


def _validate_storage_preflight_binding(
    storage: dict[str, object],
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    profile: str,
) -> None:
    if profile == _HISTORICAL_RETIREMENT_STORAGE_PREFLIGHT_PROFILE:
        if "storage_preflight_evidence" in storage:
            raise ValueError(
                "Historical storage-preflight compatibility accepts only absent evidence"
            )
        return
    if profile not in {
        _STRICT_STORAGE_PREFLIGHT_PROFILE,
        _EXACT_REVIEWED_STORAGE_PREFLIGHT_PREDECESSOR_PROFILE,
    }:
        raise ValueError("Unsupported storage-preflight validation profile")
    binding = storage.get("storage_preflight_evidence")
    if (
        not isinstance(binding, dict)
        or set(binding) != {"orion_path", "probe_id", "project_home_path", "sha256"}
    ):
        raise ValueError("Storage policy storage-preflight evidence binding is malformed")
    probe_id = binding.get("probe_id")
    if not isinstance(probe_id, str):
        raise ValueError("Storage policy storage-preflight evidence probe ID is malformed")
    try:
        parsed_probe_id = uuid.UUID(probe_id)
    except ValueError as error:
        raise ValueError(
            "Storage policy storage-preflight evidence probe ID is malformed"
        ) from error
    if str(parsed_probe_id) != probe_id:
        raise ValueError("Storage policy storage-preflight evidence probe ID is malformed")
    lexical_pic_root = Path(os.path.abspath(str(storage["orion_bulk_evidence_root"])))
    lexical_project_home_root = Path(
        os.path.abspath(str(storage["project_home_mirror_root"]))
    )
    evidence_relative = Path("policy") / "storage_preflight_evidence" / f"{probe_id}.json"
    orion_path = lexical_pic_root / evidence_relative
    project_home_path = lexical_project_home_root / evidence_relative
    if binding.get("orion_path") != str(orion_path) or binding.get(
        "project_home_path"
    ) != str(project_home_path):
        raise ValueError("Storage policy storage-preflight evidence paths are malformed")
    digest = binding.get("sha256")
    if not _is_lowercase_sha256(digest):
        raise ValueError("Storage policy storage-preflight evidence SHA-256 is malformed")
    orion_payload = read_stable_regular_file_below(
        orion_path, lexical_pic_root, require_read_only_mode=True
    )
    project_home_payload = read_stable_regular_file_below(
        project_home_path,
        lexical_project_home_root,
        require_read_only_mode=True,
    )
    if orion_payload != project_home_payload or sha256_bytes(orion_payload) != digest:
        raise ValueError("Storage-preflight evidence mirrored bytes differ")
    _validate_storage_preflight_artifact(
        orion_payload,
        probe_id=probe_id,
        last_preflight_utc=str(storage["last_preflight_utc"]),
        orion_path=orion_path,
        project_home_path=project_home_path,
        authorized_pic_root=lexical_pic_root,
        authorized_project_home_root=lexical_project_home_root,
        authorized_source_authentication=(
            AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_SOURCE_AUTHENTICATION
            if profile == _EXACT_REVIEWED_STORAGE_PREFLIGHT_PREDECESSOR_PROFILE
            else None
        ),
    )


def _validate_storage_policy(
    policy: dict[str, object],
    *,
    control_plane_version: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
    ledger_mirror_transport: str = AUTHORIZED_LEDGER_MIRROR_TRANSPORT,
    allow_pending_genesis: bool = False,
    storage_preflight_profile: str,
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
        "last_preflight_utc",
    }
    authenticated_storage_preflight_profiles = {
        _STRICT_STORAGE_PREFLIGHT_PROFILE,
        _EXACT_REVIEWED_STORAGE_PREFLIGHT_PREDECESSOR_PROFILE,
    }
    if storage_preflight_profile in authenticated_storage_preflight_profiles:
        if "storage_preflight_evidence" not in storage:
            raise ValueError(
                "Storage policy lacks authenticated storage-preflight evidence"
            )
        required_storage_keys.add("storage_preflight_evidence")
    elif storage_preflight_profile != _HISTORICAL_RETIREMENT_STORAGE_PREFLIGHT_PROFILE:
        raise ValueError("Unsupported storage-preflight validation profile")
    allowed_storage_keys = required_storage_keys | {
        "status",
        "storage_preflight_evidence",
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
    lexical_project_home_mirror_root = (
        Path(os.path.abspath(project_home_mirror_root))
        if isinstance(project_home_mirror_root, str)
        else None
    )
    if (
        not isinstance(project_home_mirror_root, str)
        or (
            storage_preflight_profile in authenticated_storage_preflight_profiles
            and (
                lexical_project_home_mirror_root
                != Path(os.path.abspath(authorized_project_home_root))
                or Path(os.path.abspath(authorized_project_home_root)).resolve()
                != Path(os.path.abspath(authorized_project_home_root))
            )
        )
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
        (
            "orion_simulation_root_preflight",
            Path(os.path.abspath(str(storage["orion_bulk_evidence_root"]))),
        ),
        (
            "project_home_preflight",
            Path(os.path.abspath(str(storage["project_home_mirror_root"]))),
        ),
    ]:
        record = storage.get(key)
        if (
            not isinstance(record, dict)
            or not set(record) <= {"status", "path", "method"}
            or record.get("status") != "passed"
            or record.get("method") != AUTHORIZED_STORAGE_PREFLIGHT_METHOD
        ):
            raise ValueError(f"Storage policy {key} has not passed")
        if (
            storage_preflight_profile
            == _HISTORICAL_RETIREMENT_STORAGE_PREFLIGHT_PROFILE
            and key == "project_home_preflight"
        ):
            if set(record) != {"status", "method"}:
                raise ValueError(
                    "Historical storage policy project_home_preflight shape is malformed"
                )
            continue
        if (
            set(record) != {"status", "path", "method"}
            or not isinstance(record["path"], str)
            or record["path"] != str(expected_path)
        ):
            raise ValueError(f"Storage policy {key}.path is not authorized")
    _validate_storage_preflight_binding(
        storage,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        profile=storage_preflight_profile,
    )
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
            Path(os.path.abspath(str(project_home_mirror_root)))
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
    """Validate one normal policy with authenticated mirrored storage evidence."""
    return _validate_storage_policy(
        policy,
        control_plane_version=control_plane_version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        ledger_mirror_transport=ledger_mirror_transport,
        allow_pending_genesis=allow_pending_genesis,
        storage_preflight_profile=_STRICT_STORAGE_PREFLIGHT_PROFILE,
    )


def require_policy_predecessor_snapshot_for_promotion(
    *,
    successor_policy: dict[str, object],
    successor_control_plane_version: str,
    permit_historical_retirement_predecessor: bool = False,
    permit_exact_reviewed_storage_preflight_predecessor: bool = False,
    permit_exact_authorized_clean_candidate_freeze_replacement: bool = False,
    permit_completed_q043_registered_slice_retirement: bool = False,
    permit_completed_q023_registered_slice_retirement: bool = False,
    q043_registered_matrix_path: Path | None = None,
    q043_registered_matrix_sha256: str | None = None,
    q023_registered_matrix_path: Path | None = None,
    q023_registered_matrix_sha256: str | None = None,
    expected_active_policy_sha256: str | None = None,
    expected_active_promotion_sha256: str | None = None,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
    ledger_mirror_transport: str = AUTHORIZED_LEDGER_MIRROR_TRANSPORT,
    predecessor_policy_bytes: bytes | None = None,
    predecessor_promotion_bytes: bytes | None = None,
    allow_active_promotion_transaction: bool = False,
) -> tuple[dict[str, object], dict[str, str]] | None:
    """Revalidate an active predecessor, with narrow one-use migration modes."""
    supplied_predecessor_bytes = (
        predecessor_policy_bytes is not None
        or predecessor_promotion_bytes is not None
    )
    if (
        supplied_predecessor_bytes
        and (
            type(predecessor_policy_bytes) is not bytes
            or type(predecessor_promotion_bytes) is not bytes
            or not allow_active_promotion_transaction
        )
    ):
        raise ValueError(
            "Prepared-transaction predecessor bytes require one complete locked "
            "recovery snapshot"
        )
    if not allow_active_promotion_transaction:
        require_no_active_promotion_transaction(
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
    transition_modes = [
        permit_historical_retirement_predecessor,
        permit_exact_reviewed_storage_preflight_predecessor,
        permit_exact_authorized_clean_candidate_freeze_replacement,
        permit_completed_q043_registered_slice_retirement,
        permit_completed_q023_registered_slice_retirement,
    ]
    if sum(transition_modes) > 1:
        raise ValueError("Policy predecessor transition modes are exclusive")
    exact_predecessor_mode = (
        permit_exact_authorized_clean_candidate_freeze_replacement
        or permit_completed_q043_registered_slice_retirement
        or permit_completed_q023_registered_slice_retirement
    )
    if exact_predecessor_mode:
        if not _is_lowercase_sha256(
            expected_active_policy_sha256
        ) or not _is_lowercase_sha256(expected_active_promotion_sha256):
            if permit_exact_authorized_clean_candidate_freeze_replacement:
                raise ValueError(
                    "Exact authorized clean-candidate freeze replacement requires "
                    "both exact active predecessor hashes"
                )
            campaign = (
                "Q043"
                if permit_completed_q043_registered_slice_retirement
                else "Q023"
            )
            raise ValueError(
                f"Completed {campaign} retirement requires both exact active "
                "predecessor hashes"
            )
    elif (
        expected_active_policy_sha256 is not None
        or expected_active_promotion_sha256 is not None
    ):
        raise ValueError(
            "Exact active predecessor hashes require an exact policy transition mode"
        )
    if permit_completed_q043_registered_slice_retirement:
        if not isinstance(q043_registered_matrix_path, Path) or not _is_lowercase_sha256(
            q043_registered_matrix_sha256
        ):
            raise ValueError(
                "Completed Q043 retirement requires the immutable matrix path and digest"
            )
    elif (
        q043_registered_matrix_path is not None
        or q043_registered_matrix_sha256 is not None
    ):
        raise ValueError(
            "Q043 registered matrix bindings require completed Q043 retirement mode"
        )
    if permit_completed_q023_registered_slice_retirement:
        if not isinstance(q023_registered_matrix_path, Path) or not _is_lowercase_sha256(
            q023_registered_matrix_sha256
        ):
            raise ValueError(
                "Completed Q023 retirement requires the immutable matrix path and digest"
            )
    elif (
        q023_registered_matrix_path is not None
        or q023_registered_matrix_sha256 is not None
    ):
        raise ValueError(
            "Q023 registered matrix bindings require completed Q023 retirement mode"
        )
    policy_path = canonical_policy_path(authorized_pic_root)
    mirror_policy_path = canonical_policy_path(authorized_project_home_root)
    promotion_path = active_promotion_path(authorized_pic_root)
    mirror_promotion_path = active_promotion_path(authorized_project_home_root)
    anchor_paths = [
        policy_path,
        mirror_policy_path,
        promotion_path,
        mirror_promotion_path,
    ]
    if supplied_predecessor_bytes:
        assert isinstance(predecessor_policy_bytes, bytes)
        assert isinstance(predecessor_promotion_bytes, bytes)
        policy_bytes = predecessor_policy_bytes
        promotion_bytes = predecessor_promotion_bytes
    else:
        anchor_exists = [path.exists() or path.is_symlink() for path in anchor_paths]
        if not any(anchor_exists):
            if (
                permit_historical_retirement_predecessor
                or permit_exact_reviewed_storage_preflight_predecessor
                or permit_exact_authorized_clean_candidate_freeze_replacement
                or permit_completed_q043_registered_slice_retirement
                or permit_completed_q023_registered_slice_retirement
            ):
                raise ValueError(
                    "Policy predecessor transition requires an active predecessor"
                )
            return None
        if not all(anchor_exists):
            raise ValueError("Active-policy predecessor anchors are incomplete")
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
        if promotion_bytes != artifacts[mirror_promotion_path]:
            raise ValueError("Orion and Project Home predecessor promotion bytes differ")
        policy_bytes = artifacts[policy_path]
        if policy_bytes != artifacts[mirror_policy_path]:
            raise ValueError("Orion and Project Home predecessor policy bytes differ")
    promotion = read_json_bytes(promotion_bytes, label=str(promotion_path))
    if exact_predecessor_mode and (
        sha256_bytes(policy_bytes) != expected_active_policy_sha256
        or sha256_bytes(promotion_bytes) != expected_active_promotion_sha256
    ):
        raise ValueError(
            "Exact policy transition active predecessor hashes changed"
        )
    predecessor_version = promotion.get("control_plane_version")
    if not _is_lowercase_sha256(predecessor_version):
        raise ValueError("Active-policy predecessor control-plane version is malformed")
    assert isinstance(predecessor_version, str)
    inventory: dict[str, object] | None = None
    for root in [authorized_pic_root, authorized_project_home_root]:
        lexical_root = Path(os.path.abspath(root))
        current = verify_historical_installed_control_plane(
            lexical_root / "control_plane" / predecessor_version,
            authorized_pic_root=root,
        )
        if inventory is not None and current != inventory:
            raise ValueError(
                "Historical Orion and Project Home control-plane inventories differ"
            )
        inventory = current
    policy = read_json_bytes(policy_bytes, label=str(policy_path))
    predecessor_storage = policy.get("olcf_side_storage")
    configured_project_home_root = (
        Path(os.path.abspath(str(predecessor_storage["project_home_mirror_root"])))
        if isinstance(predecessor_storage, dict)
        and isinstance(predecessor_storage.get("project_home_mirror_root"), str)
        else Path(os.path.abspath(authorized_project_home_root))
    )
    if (
        configured_project_home_root.resolve(strict=True)
        != Path(os.path.abspath(authorized_project_home_root)).resolve(strict=True)
    ):
        raise ValueError("Active-policy predecessor Project Home root is not authorized")
    expected = {
        "control_plane_version": predecessor_version,
        "policy_path": str(policy_path),
        "project_home_policy_path": str(
            canonical_policy_path(configured_project_home_root)
        ),
        "policy_sha256": sha256_bytes(policy_bytes),
    }
    slices = policy.get("registered_science_slices")
    if isinstance(slices, list) and slices:
        attestation_path = Path(
            str(promotion.get("pre_policy_promotion_attestation_path", ""))
        )
        authorization_id = promotion.get(
            "pre_policy_promotion_attestation_authorization_id"
        )
        digest = promotion.get("pre_policy_promotion_attestation_sha256")
        if (
            not isinstance(authorization_id, str)
            or not authorization_id
            or not _is_lowercase_sha256(digest)
            or attestation_path.name != "attestation.json"
            or attestation_path.parent.parent
            != Path(os.path.abspath(authorized_pic_root)) / "operator_attestations"
            or not attestation_path.parent.name.endswith("-pre_policy_promotion")
        ):
            raise ValueError("Active-policy predecessor attestation binding is malformed")
        attestation = validate_sealed_operator_attestation(
            attestation_path,
            authorization_id=authorization_id,
            phase="pre_policy_promotion",
            control_plane_version=predecessor_version,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=configured_project_home_root,
            enforce_freshness=False,
        )
        if attestation != {"path": str(attestation_path), "sha256": digest}:
            raise ValueError("Active-policy predecessor attestation binding differs")
        expected.update(
            {
                "pre_policy_promotion_attestation_authorization_id": authorization_id,
                "pre_policy_promotion_attestation_path": str(attestation_path),
                "pre_policy_promotion_attestation_sha256": digest,
            }
        )
    promotion_schema_version = promotion.get("schema_version")
    if (
        type(promotion_schema_version) is int
        and promotion_schema_version == 2
        and _is_canonical_uuid(promotion.get("promotion_id"))
    ):
        expected.update(
            {
                "schema_version": 2,
                "promotion_id": promotion["promotion_id"],
            }
        )
    elif (
        type(promotion_schema_version) is int
        and promotion_schema_version == 1
        and "promotion_id" not in promotion
        and (
            permit_historical_retirement_predecessor
            or permit_exact_reviewed_storage_preflight_predecessor
        )
    ):
        expected["schema_version"] = 1
    else:
        raise ValueError("Active-policy predecessor promotion record is malformed")
    if promotion != expected:
        raise ValueError("Active-policy predecessor promotion record is malformed")
    storage = policy.get("olcf_side_storage")
    has_evidence = isinstance(storage, dict) and "storage_preflight_evidence" in storage
    profile = _STRICT_STORAGE_PREFLIGHT_PROFILE
    if permit_historical_retirement_predecessor and has_evidence:
        raise ValueError(
            "Historical storage-preflight retirement flag requires a legacy predecessor"
        )
    if permit_exact_reviewed_storage_preflight_predecessor:
        if not has_evidence:
            raise ValueError(
                "Exact reviewed storage-preflight predecessor migration requires "
                "authenticated predecessor evidence"
            )
        if (
            sha256_bytes(policy_bytes)
            != AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_POLICY_SHA256
            or sha256_bytes(promotion_bytes)
            != AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_PROMOTION_SHA256
            or predecessor_version
            != AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_CONTROL_PLANE_VERSION
        ):
            raise ValueError(
                "Exact reviewed storage-preflight predecessor differs from the exact "
                "reviewed live anchors"
            )
        assert isinstance(storage, dict)
        predecessor_binding = storage["storage_preflight_evidence"]
        if (
            not isinstance(predecessor_binding, dict)
            or predecessor_binding.get("probe_id")
            != AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_PROBE_ID
            or predecessor_binding.get("sha256")
            != AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_EVIDENCE_SHA256
        ):
            raise ValueError(
                "Exact reviewed storage-preflight predecessor evidence differs from "
                "the exact reviewed live binding"
            )
        if predecessor_version == successor_control_plane_version:
            raise ValueError(
                "Exact reviewed storage-preflight predecessor migration requires "
                "a new control plane"
            )
        if successor_policy.get("registered_science_slices") != []:
            raise ValueError(
                "Exact reviewed storage-preflight predecessor migration requires one "
                "empty-allowlist replacement"
            )
        successor_storage = successor_policy.get("olcf_side_storage")
        if not isinstance(successor_storage, dict):
            raise ValueError(
                "Exact reviewed storage-preflight successor storage is malformed"
            )
        predecessor_comparable = {
            **policy,
            "olcf_side_storage": {
                **storage,
                "installed_control_plane_version": None,
                "staged_control_plane_candidate_version": None,
                "last_preflight_utc": None,
                "storage_preflight_evidence": None,
            },
        }
        successor_comparable = {
            **successor_policy,
            "olcf_side_storage": {
                **successor_storage,
                "installed_control_plane_version": None,
                "staged_control_plane_candidate_version": None,
                "last_preflight_utc": None,
                "storage_preflight_evidence": None,
            },
        }
        if not _exact_json_equal(successor_comparable, predecessor_comparable):
            raise ValueError(
                "Exact reviewed storage-preflight successor is not the exact "
                "authorized predecessor transformation"
            )
        successor_binding = successor_storage.get("storage_preflight_evidence")
        if (
            not isinstance(successor_binding, dict)
            or successor_binding == predecessor_binding
            or successor_binding.get("probe_id") == predecessor_binding.get("probe_id")
            or successor_binding.get("sha256") == predecessor_binding.get("sha256")
        ):
            raise ValueError(
                "Exact reviewed storage-preflight successor requires a different "
                "authenticated evidence binding"
            )
        if utc_datetime(
            successor_storage.get("last_preflight_utc"),
            field="successor.olcf_side_storage.last_preflight_utc",
        ) <= utc_datetime(
            storage.get("last_preflight_utc"),
            field="predecessor.olcf_side_storage.last_preflight_utc",
        ):
            raise ValueError(
                "Exact reviewed storage-preflight successor requires a newer "
                "authenticated storage preflight"
            )
        profile = _EXACT_REVIEWED_STORAGE_PREFLIGHT_PREDECESSOR_PROFILE
    if permit_completed_q043_registered_slice_retirement:
        if not isinstance(storage, dict):
            raise ValueError(
                "Completed Q043 retirement predecessor storage is malformed"
            )
        if predecessor_version == successor_control_plane_version:
            raise ValueError(
                "Completed Q043 retirement requires a new control plane"
            )
        predecessor_slices = policy.get("registered_science_slices")
        if (
            not isinstance(predecessor_slices, list)
            or len(predecessor_slices) != Q043_REGISTERED_CASE_COUNT
            or successor_policy.get("registered_science_slices") != []
        ):
            raise ValueError(
                "Completed Q043 retirement requires the exact nonempty Q043 "
                "predecessor and one empty-allowlist successor"
            )
        assert q043_registered_matrix_path is not None
        expected_matrix_path = (
            Path(os.path.abspath(authorized_pic_root))
            / Q043_REGISTERED_MATRIX_RELATIVE
        )
        if Path(os.path.abspath(q043_registered_matrix_path)) != expected_matrix_path:
            raise ValueError(
                "Completed Q043 retirement matrix path is not the canonical artifact"
            )
        matrix_bytes = read_stable_regular_file_below(
            q043_registered_matrix_path,
            authorized_pic_root,
            require_read_only_mode=True,
        )
        if sha256_bytes(matrix_bytes) != q043_registered_matrix_sha256:
            raise ValueError(
                "Completed Q043 retirement matrix digest changed"
            )
        matrix = read_json_bytes(
            matrix_bytes, label=str(q043_registered_matrix_path)
        )
        admissions = matrix.get("case_admissions")
        authorization = matrix.get("authorization")
        if (
            matrix.get("record_type") != Q043_REGISTERED_MATRIX_RECORD_TYPE
            or matrix.get("status") != Q043_REGISTERED_MATRIX_STATUS
            or matrix.get("campaign_id") != Q043_REGISTERED_CAMPAIGN
            or matrix.get("registered_execution_qualification_check_pass") is not True
            or matrix.get(
                "source_local_matrix_result_sufficient_for_downstream_qualification"
            )
            is not False
            or matrix.get("case_count") != Q043_REGISTERED_CASE_COUNT
            or not isinstance(admissions, list)
            or len(admissions) != Q043_REGISTERED_CASE_COUNT
            or not isinstance(authorization, dict)
            or not authorization
            or any(
                value is not False
                for key, value in authorization.items()
                if key.endswith("_authorized")
            )
        ):
            raise ValueError(
                "Completed Q043 retirement matrix qualification boundary drifted"
            )
        policy_pairs = []
        for index, item in enumerate(predecessor_slices, 1):
            expected_authorization = f"q043-re-{index:03d}-v1"
            if (
                not isinstance(item, dict)
                or item.get("authorization_id") != expected_authorization
                or item.get("status") != AUTHORIZED_REGISTERED_SCIENCE_SLICE_STATUS
                or item.get("campaign") != Q043_REGISTERED_CAMPAIGN
                or not isinstance(item.get("test_id"), str)
            ):
                raise ValueError(
                    "Completed Q043 retirement predecessor slice matrix drifted"
                )
            policy_pairs.append((expected_authorization, item["test_id"]))
        matrix_pairs = []
        for item in admissions:
            execution = item.get("execution_binding") if isinstance(item, dict) else None
            identity = (
                execution.get("registered_execution_identity")
                if isinstance(execution, dict)
                else None
            )
            if (
                not isinstance(item, dict)
                or not isinstance(item.get("case_id"), str)
                or item.get("campaign_id") != Q043_REGISTERED_CAMPAIGN
                or item.get("status") != Q043_REGISTERED_CASE_STATUS
                or not isinstance(identity, dict)
                or not isinstance(
                    identity.get("registered_science_authorization_id"), str
                )
            ):
                raise ValueError(
                    "Completed Q043 retirement matrix case binding is malformed"
                )
            matrix_pairs.append(
                (
                    identity["registered_science_authorization_id"],
                    item["case_id"],
                )
            )
        if (
            matrix_pairs != policy_pairs
            or len(set(policy_pairs)) != Q043_REGISTERED_CASE_COUNT
        ):
            raise ValueError(
                "Completed Q043 retirement matrix differs from active authorizations"
            )
        successor_storage = successor_policy.get("olcf_side_storage")
        if not isinstance(successor_storage, dict):
            raise ValueError("Completed Q043 retirement successor storage is malformed")
        predecessor_comparable = {
            **policy,
            "registered_science_slices": [],
            "olcf_side_storage": {
                **storage,
                "installed_control_plane_version": None,
                "staged_control_plane_candidate_version": None,
                "last_preflight_utc": None,
                "storage_preflight_evidence": None,
            },
        }
        successor_comparable = {
            **successor_policy,
            "olcf_side_storage": {
                **successor_storage,
                "installed_control_plane_version": None,
                "staged_control_plane_candidate_version": None,
                "last_preflight_utc": None,
                "storage_preflight_evidence": None,
            },
        }
        if not _exact_json_equal(successor_comparable, predecessor_comparable):
            raise ValueError(
                "Completed Q043 retirement changed unrelated policy fields"
            )
        predecessor_binding = storage.get("storage_preflight_evidence")
        successor_binding = successor_storage.get("storage_preflight_evidence")
        if (
            not isinstance(predecessor_binding, dict)
            or not isinstance(successor_binding, dict)
            or successor_binding == predecessor_binding
            or successor_binding.get("probe_id") == predecessor_binding.get("probe_id")
            or successor_binding.get("sha256") == predecessor_binding.get("sha256")
            or utc_datetime(
                successor_storage.get("last_preflight_utc"),
                field="successor.olcf_side_storage.last_preflight_utc",
            )
            <= utc_datetime(
                storage.get("last_preflight_utc"),
                field="predecessor.olcf_side_storage.last_preflight_utc",
            )
        ):
            raise ValueError(
                "Completed Q043 retirement requires newer authenticated storage evidence"
            )
    if permit_completed_q023_registered_slice_retirement:
        if not isinstance(storage, dict):
            raise ValueError(
                "Completed Q023 retirement predecessor storage is malformed"
            )
        if predecessor_version == successor_control_plane_version:
            raise ValueError(
                "Completed Q023 retirement requires a new control plane"
            )
        predecessor_slices = policy.get("registered_science_slices")
        if (
            not isinstance(predecessor_slices, list)
            or len(predecessor_slices) != Q023_REGISTERED_CASE_COUNT
            or successor_policy.get("registered_science_slices") != []
        ):
            raise ValueError(
                "Completed Q023 retirement requires the exact nonempty Q023 "
                "predecessor and one empty-allowlist successor"
            )
        assert q023_registered_matrix_path is not None
        expected_matrix_path = (
            Path(os.path.abspath(authorized_pic_root))
            / Q023_REGISTERED_MATRIX_RELATIVE
        )
        if Path(os.path.abspath(q023_registered_matrix_path)) != expected_matrix_path:
            raise ValueError(
                "Completed Q023 retirement matrix path is not the canonical artifact"
            )
        matrix_bytes = read_stable_regular_file_below(
            q023_registered_matrix_path,
            authorized_pic_root,
            require_read_only_mode=True,
        )
        if sha256_bytes(matrix_bytes) != q023_registered_matrix_sha256:
            raise ValueError("Completed Q023 retirement matrix digest changed")
        matrix = read_json_bytes(
            matrix_bytes, label=str(q023_registered_matrix_path)
        )
        admissions = matrix.get("case_admissions")
        authorization = matrix.get("authorization")
        if (
            matrix.get("record_type") != Q023_REGISTERED_MATRIX_RECORD_TYPE
            or matrix.get("status") != Q023_REGISTERED_MATRIX_STATUS
            or matrix.get("campaign_id") != Q023_REGISTERED_CAMPAIGN_ID
            or matrix.get("registered_execution_qualification_check_pass") is not True
            or matrix.get("registered_linear_qualification_pass") is not True
            or matrix.get("case_count") != Q023_REGISTERED_CASE_COUNT
            or not isinstance(admissions, list)
            or len(admissions) != Q023_REGISTERED_CASE_COUNT
            or not isinstance(authorization, dict)
            or not authorization
            or any(value is not False for value in authorization.values())
        ):
            raise ValueError(
                "Completed Q023 retirement matrix qualification boundary drifted"
            )
        policy_pairs = []
        for index, item in enumerate(predecessor_slices, 1):
            expected_authorization = f"q023-linear-{index:03d}-v1"
            if (
                not isinstance(item, dict)
                or item.get("authorization_id") != expected_authorization
                or item.get("status") != AUTHORIZED_REGISTERED_SCIENCE_SLICE_STATUS
                or item.get("campaign") != Q023_REGISTERED_POLICY_CAMPAIGN
                or not isinstance(item.get("test_id"), str)
            ):
                raise ValueError(
                    "Completed Q023 retirement predecessor slice matrix drifted"
                )
            policy_pairs.append((expected_authorization, item["test_id"]))
        matrix_pairs = []
        for item in admissions:
            identity = item.get("execution_identity") if isinstance(item, dict) else None
            if (
                not isinstance(item, dict)
                or not isinstance(item.get("member_id"), str)
                or item.get("campaign_id") != Q023_REGISTERED_CAMPAIGN_ID
                or item.get("status") != Q023_REGISTERED_CASE_STATUS
                or not isinstance(identity, dict)
                or not isinstance(
                    identity.get("registered_science_authorization_id"), str
                )
            ):
                raise ValueError(
                    "Completed Q023 retirement matrix case binding is malformed"
                )
            matrix_pairs.append(
                (
                    identity["registered_science_authorization_id"],
                    item["member_id"],
                )
            )
        if (
            matrix_pairs != policy_pairs
            or len(set(policy_pairs)) != Q023_REGISTERED_CASE_COUNT
        ):
            raise ValueError(
                "Completed Q023 retirement matrix differs from active authorizations"
            )
        successor_storage = successor_policy.get("olcf_side_storage")
        if not isinstance(successor_storage, dict):
            raise ValueError("Completed Q023 retirement successor storage is malformed")
        predecessor_comparable = {
            **policy,
            "registered_science_slices": [],
            "olcf_side_storage": {
                **storage,
                "installed_control_plane_version": None,
                "staged_control_plane_candidate_version": None,
                "last_preflight_utc": None,
                "storage_preflight_evidence": None,
            },
        }
        successor_comparable = {
            **successor_policy,
            "olcf_side_storage": {
                **successor_storage,
                "installed_control_plane_version": None,
                "staged_control_plane_candidate_version": None,
                "last_preflight_utc": None,
                "storage_preflight_evidence": None,
            },
        }
        if not _exact_json_equal(successor_comparable, predecessor_comparable):
            raise ValueError(
                "Completed Q023 retirement changed unrelated policy fields"
            )
        predecessor_binding = storage.get("storage_preflight_evidence")
        successor_binding = successor_storage.get("storage_preflight_evidence")
        if (
            not isinstance(predecessor_binding, dict)
            or not isinstance(successor_binding, dict)
            or successor_binding == predecessor_binding
            or successor_binding.get("probe_id") == predecessor_binding.get("probe_id")
            or successor_binding.get("sha256") == predecessor_binding.get("sha256")
            or utc_datetime(
                successor_storage.get("last_preflight_utc"),
                field="successor.olcf_side_storage.last_preflight_utc",
            )
            <= utc_datetime(
                storage.get("last_preflight_utc"),
                field="predecessor.olcf_side_storage.last_preflight_utc",
            )
        ):
            raise ValueError(
                "Completed Q023 retirement requires newer authenticated storage evidence"
            )
    predecessor_slices = policy.get("registered_science_slices")
    successor_slices = successor_policy.get("registered_science_slices")
    if (
        isinstance(predecessor_slices, list)
        and predecessor_slices
        and not _exact_json_equal(predecessor_slices, successor_slices)
        and not permit_completed_q043_registered_slice_retirement
        and not permit_completed_q023_registered_slice_retirement
    ):
        raise ValueError(
            "A nonempty registered-science allowlist can change only through an "
            "evidence-bound retirement transition"
        )
    predecessor_freeze = policy.get("science_submission_freeze")
    successor_freeze = successor_policy.get("science_submission_freeze")
    predecessor_has_authorized_freeze = (
        isinstance(predecessor_freeze, dict)
        and predecessor_freeze.get("status") == AUTHORIZED_CLEAN_CANDIDATE_FREEZE
    )
    successor_has_authorized_freeze = (
        isinstance(successor_freeze, dict)
        and successor_freeze.get("status") == AUTHORIZED_CLEAN_CANDIDATE_FREEZE
    )
    authorized_freeze_changed = (
        predecessor_has_authorized_freeze
        and not _exact_json_equal(predecessor_freeze, successor_freeze)
    )
    if authorized_freeze_changed and not (
        permit_exact_authorized_clean_candidate_freeze_replacement
        or permit_historical_retirement_predecessor
    ):
        raise ValueError(
            "Any transition away from an authorized clean-candidate freeze requires "
            "an exact active-predecessor transition mode"
        )
    if permit_exact_authorized_clean_candidate_freeze_replacement:
        expected_freeze_keys = {
            "status",
            "manifest_path",
            "manifest_sha256",
            "build_profile_control_plane_version",
        }
        if (
            predecessor_version != successor_control_plane_version
            or policy.get("registered_science_slices") != []
            or successor_policy.get("registered_science_slices") != []
        ):
            raise ValueError(
                "Exact authorized clean-candidate freeze replacement requires the "
                "same control plane and empty allowlists"
            )
        if (
            not isinstance(predecessor_freeze, dict)
            or not isinstance(successor_freeze, dict)
            or set(predecessor_freeze) != expected_freeze_keys
            or set(successor_freeze) != expected_freeze_keys
            or predecessor_freeze.get("status")
            != AUTHORIZED_CLEAN_CANDIDATE_FREEZE
            or successor_freeze.get("status") != AUTHORIZED_CLEAN_CANDIDATE_FREEZE
            or successor_freeze.get("build_profile_control_plane_version")
            != successor_control_plane_version
        ):
            raise ValueError(
                "Exact authorized clean-candidate freeze replacement requires exact "
                "authorized freeze bindings"
            )
        predecessor_comparable = {
            **policy,
            "science_submission_freeze": None,
        }
        successor_comparable = {
            **successor_policy,
            "science_submission_freeze": None,
        }
        if not _exact_json_equal(successor_comparable, predecessor_comparable):
            raise ValueError(
                "Exact authorized clean-candidate freeze replacement changed "
                "unrelated policy fields"
            )
        if (
            not authorized_freeze_changed
            or successor_freeze.get("manifest_path")
            == predecessor_freeze.get("manifest_path")
            or successor_freeze.get("manifest_sha256")
            == predecessor_freeze.get("manifest_sha256")
        ):
            raise ValueError(
                "Exact authorized clean-candidate freeze replacement requires a "
                "different manifest path and digest"
            )
    if not has_evidence:
        if not permit_historical_retirement_predecessor:
            raise ValueError(
                "Active-policy predecessor lacks authenticated storage-preflight evidence"
            )
        if predecessor_version == successor_control_plane_version:
            raise ValueError(
                "Historical storage-preflight retirement requires a new control plane"
            )
        if (
            successor_policy.get("registered_science_slices") != []
            or successor_policy.get("science_submission_freeze")
            != {"status": PENDING_CLEAN_CANDIDATE_FREEZE}
        ):
            raise ValueError(
                "Historical storage-preflight retirement requires one launch-prohibited "
                "pending-freeze replacement"
            )
        if (
            sha256_bytes(policy_bytes)
            != AUTHORIZED_HISTORICAL_STORAGE_PREFLIGHT_RETIREMENT_POLICY_SHA256
            or sha256_bytes(promotion_bytes)
            != AUTHORIZED_HISTORICAL_STORAGE_PREFLIGHT_RETIREMENT_PROMOTION_SHA256
        ):
            raise ValueError(
                "Historical storage-preflight retirement predecessor differs from "
                "the exact reviewed live anchors"
            )
        profile = _HISTORICAL_RETIREMENT_STORAGE_PREFLIGHT_PROFILE
    _validate_storage_policy(
        policy,
        control_plane_version=predecessor_version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        ledger_mirror_transport=ledger_mirror_transport,
        allow_pending_genesis=True,
        storage_preflight_profile=profile,
    )
    if not allow_active_promotion_transaction:
        require_no_active_promotion_transaction(
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
    return policy, {
        "active_policy_sha256": sha256_bytes(policy_bytes),
        "active_promotion_sha256": sha256_bytes(promotion_bytes),
    }


def require_storage_policy_unlock_snapshot(
    *,
    control_plane_version: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
    ledger_mirror_transport: str = AUTHORIZED_LEDGER_MIRROR_TRANSPORT,
    allow_pending_genesis: bool = False,
    allow_active_promotion_transaction: bool = False,
) -> tuple[dict[str, object], dict[str, str]]:
    if not allow_active_promotion_transaction:
        require_no_active_promotion_transaction(
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
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
        "schema_version": 2,
        "promotion_id": promotion.get("promotion_id"),
        "control_plane_version": control_plane_version,
        "policy_path": str(policy_path),
        "project_home_policy_path": str(mirror_policy_path),
        "policy_sha256": policy_sha256,
    }
    policy_value = read_json_bytes(policy_bytes, label=str(policy_path))
    registered_science_slices = policy_value.get("registered_science_slices")
    if isinstance(registered_science_slices, list) and registered_science_slices:
        attestation_path = Path(
            str(promotion.get("pre_policy_promotion_attestation_path", ""))
        )
        if (
            not isinstance(
                promotion.get("pre_policy_promotion_attestation_authorization_id"),
                str,
            )
            or not promotion["pre_policy_promotion_attestation_authorization_id"]
            or not _is_lowercase_sha256(
                promotion.get("pre_policy_promotion_attestation_sha256")
            )
            or attestation_path.name != "attestation.json"
            or attestation_path.parent.parent
            != Path(os.path.abspath(authorized_pic_root)) / "operator_attestations"
            or not attestation_path.parent.name.endswith("-pre_policy_promotion")
        ):
            raise ValueError("Active-policy promotion attestation binding is malformed")
        attestation = validate_sealed_operator_attestation(
            attestation_path,
            authorization_id=promotion[
                "pre_policy_promotion_attestation_authorization_id"
            ],
            phase="pre_policy_promotion",
            control_plane_version=control_plane_version,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=project_home_ledger_root(
                authorized_project_home_root
            ),
            enforce_freshness=False,
        )
        if attestation != {
            "path": str(attestation_path),
            "sha256": promotion["pre_policy_promotion_attestation_sha256"],
        }:
            raise ValueError("Active-policy promotion attestation binding differs")
        expected.update(
            {
                "pre_policy_promotion_attestation_authorization_id": promotion[
                    "pre_policy_promotion_attestation_authorization_id"
                ],
                "pre_policy_promotion_attestation_path": str(attestation_path),
                "pre_policy_promotion_attestation_sha256": promotion[
                    "pre_policy_promotion_attestation_sha256"
                ],
            }
        )
    if (
        type(promotion.get("schema_version")) is not int
        or not _is_canonical_uuid(promotion.get("promotion_id"))
        or promotion != expected
    ):
        raise ValueError("Active-policy promotion record is not anchored to this control plane")
    if sha256_bytes(mirror_policy_bytes) != promotion["policy_sha256"]:
        raise ValueError("Project Home active-policy mirror checksum differs")
    policy = validate_storage_policy(
        policy_value,
        control_plane_version=control_plane_version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        ledger_mirror_transport=ledger_mirror_transport,
        allow_pending_genesis=allow_pending_genesis,
    )
    if not allow_active_promotion_transaction:
        require_no_active_promotion_transaction(
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
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
            project_home_ledger_root(authorized_project_home_root)
            / "ledger"
            / "node_hours.jsonl",
            project_home_ledger_root(authorized_project_home_root),
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
