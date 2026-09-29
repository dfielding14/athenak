#!/usr/bin/env python3
"""Prepare isolated Stage I scaling pilots and select qualified node profiles.

``prepare-wave`` retains deterministic authenticated packets beneath one
qualification-only run root.  ``submit-all-waves`` submits those exact scripts
through bounded Slurm dependency waves.  ``audit-wave`` and ``select-profile``
require complete retained scheduler and case-aware scientific evidence before
recommending a production profile.
"""

from __future__ import annotations

import argparse
from array import array
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
import errno
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import stat
import struct
import subprocess
import sys
import tempfile


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
UTILITY_RELATIVE_PATH = Path("scripts/frontier/cgl_lf_stage_i_qualification.py")
DEFAULT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
DEFAULT_MATRIX = REPOSITORY_ROOT / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
ACCOUNT = "AST207"
PARTITION = "batch"
MAX_CONCURRENT_NODES = 10
RANKS_PER_NODE = 8
CPUS_PER_TASK = 7
RESTART_SMOKE_WALLTIME = "00:05:00"
MINIMUM_RUNTIME_SAVINGS_SECONDS = 600.0
MAXIMUM_NODE_HOUR_RATIO = 1.5
RUNTIME_TIE_FRACTION = 0.10
CROSS_PROFILE_RELATIVE_TOLERANCE = 5.0e-7
CROSS_PROFILE_ABSOLUTE_TOLERANCE = 1.0e-10
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
GIT_REVISION_PATTERN = re.compile(r"[0-9a-f]{40}")
QUALIFICATION_ROOT_PATTERN = re.compile(r"qualification-[A-Za-z0-9][A-Za-z0-9_.-]*")
JOB_ID_PATTERN = re.compile(r"[1-9][0-9]*")
LIVE_JOB_ID_PATTERN = re.compile(r"[1-9][0-9]*(?:_[0-9]+|\+[0-9]+)?")
SELF_DESCRIPTOR_ENV = "_CGL_LF_STAGE_I_QUALIFICATION_DESCRIPTOR"
SELF_SOURCE_ENV = "_CGL_LF_STAGE_I_QUALIFICATION_SOURCE"
REPOSITORY_ROOT_ENV = "_CGL_LF_STAGE_I_QUALIFICATION_REPOSITORY_ROOT"
SBATCH = Path("/usr/bin/sbatch")
SACCT = Path("/usr/bin/sacct")
SQUEUE = Path("/usr/bin/squeue")
SCONTROL = Path("/usr/bin/scontrol")
SCANCEL = Path("/usr/bin/scancel")
GIT = Path("/usr/lib/git/git")
GIT_EXEC_PATH = Path("/usr/lib/git")
SYSTEM_PYTHON = Path("/usr/bin/python3.11")
TRUSTED_SYSTEM_PATH = "/usr/bin:/bin"
EXPECTED_CASE_IDS = tuple(f"R{number:02d}" for number in range(2, 18))
FROZEN_SOURCE_REVISION = "9e07542281e4e6d125582f253df3ad2e3b8b154d"
FROZEN_MATRIX_SHA256 = "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
SLURM_TIME_FORMAT = "%Y-%m-%dT%H:%M:%S%z"
QUALIFICATION_JOB_PREFIX = "cglq_"
GLOBAL_LOCK_NAME = ".cgl_lf_stage_i_qualification.lock"
EXECUTION_EPOCH = "E03-forcing-policy"
EXECUTION_EPOCH_SLUG = "E03_forcing_policy"
R17_CASE_ID = "R17"
R17_OPERATIONAL_NODES = 8
R17_OPERATIONAL_RANKS = 64
R12_CASE_ID = "R12"
R12_FRESH_RERUN_PROFILE = {
    "case_id": R12_CASE_ID,
    "segment": "s01_rankio_t0_t0p12",
    "start_time": 0.0,
    "target_time": 0.12,
    "nodes": 4,
    "ranks_per_node": RANKS_PER_NODE,
    "total_ranks": 32,
    "walltime": "02:00:00",
    "athena_walltime": "01:50:00",
    "parent_job_id": None,
    "parent_result": None,
    "parent_segment": None,
    "restart_file": None,
    "restart_file_sha256": None,
    "restart_time": None,
}
STAGE_I_LOCK_NAME = f".mks24_stage_i_{EXECUTION_EPOCH_SLUG}.lock"
STAGE_I_RESERVATIONS_NAME = (
    f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_reservations.json"
)
STAGE_I_TRANSACTION_NAMES = (
    f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_transactions",
    f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_transactions",
)
EVIDENCE_MAX_AGE_SECONDS = 24 * 60 * 60
SELECTION_VALIDITY_SECONDS = 6 * 60 * 60
FUTURE_TIMESTAMP_TOLERANCE_SECONDS = 5 * 60
SUBMISSION_RECOVERY_WINDOW_SECONDS = 5 * 60
MAX_NORMALIZED_CT_DIVB_TEXT = "1e-12"
MAX_NORMALIZED_CT_DIVB = float(MAX_NORMALIZED_CT_DIVB_TEXT)
MINIMUM_FORCING_WORK_RATE = 1.0e-4
MINIMUM_PRESSURE_WORK_RATE = 1.0e-7
MINIMUM_LF_WORK_RATE = 1.0e-7
MINIMUM_FINAL_INTERVAL_WORK_FRACTION = 1.0e-4
STRICT_LF_FAILURE_COLUMNS = (
    "lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd",
)
REQUIRED_SNAPSHOT_VARIABLES = (
    "dens", "velx", "vely", "velz", "eint", "p_perp", "bcc1", "bcc2", "bcc3",
)
MAX_BINARY_HEADER_BYTES = 16 * 1024 * 1024
MAX_BINARY_PREHEADER_LINES = 64
MAX_RESTART_PARAMETER_DUMP_BYTES = 11 * 4096 + 1
QUALIFIED_RESTART_BINARY_ABIS = {
    (
        "9e07542281e4e6d125582f253df3ad2e3b8b154d",
        "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c",
    ): {
        "mesh_time_offset_after_parameter_dump": 232,
        "mesh_time_format": "<d",
        "allowed_marker_modes": frozenset({
            "full_precision", "legacy_default_precision",
        }),
    },
}
ENDPOINT_TOLERANCE = 1.0e-10
MASS_TOLERANCE = 1.0e-12
ACTIVE_FORCING_CLOSURE_TOLERANCE = 1.0e-8
NONTRIVIAL_WORK_TOLERANCE = 1.0e-12
HISTORY_CADENCE_TOLERANCE = 1.0e-10
CASE_POLICIES = {
    "R04": {
        "profile_class": "standard_192x192x384",
        "node_profiles": (1, 2, 4),
        "target_time": 0.25,
        "walltime": "02:00:00",
        "athena_walltime": "01:50:00",
        "case_name": "paper_standard_active_random_beta10",
        "input_basename": "cgl_lf_paper_standard_active_random_beta10",
        "input_relative_path": (
            "inputs/cgl_lf_paper/cgl_lf_paper_standard_active_random_beta10.athinput"
        ),
        "input_sha256": "571eea2ccec5d069ccb1b49d132ba4c5b8bdac7ce15ee8d8a04c92327c453137",
        "resolution": "192x192x384",
        "mesh_shape": (192, 192, 384),
        "meshblock_shape": (32, 32, 64),
        "mesh_bounds": ((0.0, 1.0), (0.0, 1.0), (0.0, 2.0)),
        "physics_contract": {
            "mhd/passive": "false",
            "mhd/cgl_heat_flux": "landau_fluid",
            "mhd/lf_k_parallel": "6.283185307179586",
            "mhd/limiter_nu_coll": "1.0e10",
            "mhd/limiter_hardwall": "true",
            "problem/passive_delta": "false",
            "problem/beta0": "10.0",
            "turb_driving/driving_type": "0",
            "turb_driving/projection_policy": "mks24_random_unprojected",
        },
    },
    "R12": {
        "profile_class": "standard_192x192x384",
        "node_profiles": (R12_FRESH_RERUN_PROFILE["nodes"],),
        "target_time": R12_FRESH_RERUN_PROFILE["target_time"],
        "walltime": R12_FRESH_RERUN_PROFILE["walltime"],
        "athena_walltime": R12_FRESH_RERUN_PROFILE["athena_walltime"],
        "fixed_fresh_profile": R12_FRESH_RERUN_PROFILE,
        "case_name": "paper_heat_flux_beta10_strong",
        "input_basename": "cgl_lf_paper_heat_flux_beta10_strong",
        "input_relative_path": (
            "inputs/cgl_lf_paper/cgl_lf_paper_heat_flux_beta10_strong.athinput"
        ),
        "input_sha256": "98ddea4b4f7fec18cc40abdbf5f7c8ba5b583a91f84f4dae00e1411f23d42e7c",
        "resolution": "192x192x384",
        "mesh_shape": (192, 192, 384),
        "meshblock_shape": (32, 32, 64),
        "mesh_bounds": ((0.0, 1.0), (0.0, 1.0), (0.0, 2.0)),
        "physics_contract": {
            "mhd/passive": "false",
            "mhd/cgl_heat_flux": "landau_fluid",
            "mhd/lf_k_parallel": "0.06283185307179586",
            "mhd/limiter_nu_coll": "1.0e10",
            "mhd/limiter_hardwall": "true",
            "problem/passive_delta": "false",
            "problem/beta0": "10.0",
            "turb_driving/driving_type": "1",
            "turb_driving/projection_policy": "mks24_alfvenic_perpendicular",
        },
    },
    "R16": {
        "profile_class": "scale_separation_96x96x192",
        "node_profiles": (1, 2),
        "target_time": 1.5,
        "walltime": "01:30:00",
        "athena_walltime": "01:15:00",
        "case_name": "paper_scale_separation_active_alfvenic_beta10_nperp96",
        "input_basename": "cgl_lf_paper_scale_separation_beta10_nperp96",
        "input_relative_path": (
            "inputs/cgl_lf_paper/cgl_lf_paper_scale_separation_beta10_nperp96.athinput"
        ),
        "input_sha256": "c0ac4b54248e8f8dfb0f5fd34c0cfb4414b5330529cbf2836961c5277af3f2d1",
        "resolution": "96x96x192",
        "mesh_shape": (96, 96, 192),
        "meshblock_shape": (32, 32, 64),
        "mesh_bounds": ((0.0, 1.0), (0.0, 1.0), (0.0, 2.0)),
        "physics_contract": {
            "mhd/passive": "false",
            "mhd/cgl_heat_flux": "landau_fluid",
            "mhd/lf_k_parallel": "6.283185307179586",
            "mhd/limiter_nu_coll": "1.0e10",
            "mhd/limiter_hardwall": "true",
            "problem/passive_delta": "false",
            "problem/beta0": "10.0",
            "turb_driving/driving_type": "1",
            "turb_driving/projection_policy": "mks24_alfvenic_perpendicular",
        },
    },
    "R17": {
        "profile_class": "scale_separation_384x384x768",
        "node_profiles": (R17_OPERATIONAL_NODES,),
        "operational_only": True,
        "required_ranks": R17_OPERATIONAL_RANKS,
        "required_meshblocks_per_rank": 27,
        "target_time": 0.25,
        "walltime": "02:00:00",
        "athena_walltime": "01:50:00",
        "case_name": "paper_scale_separation_active_alfvenic_beta10_nperp384",
        "input_basename": "cgl_lf_paper_scale_separation_beta10_nperp384",
        "input_relative_path": (
            "inputs/cgl_lf_paper/"
            "cgl_lf_paper_scale_separation_beta10_nperp384.athinput"
        ),
        "input_sha256": "cc1092404b82129f807308a64f7585a6da31f45f41f1d2263acad0c8d30a7e04",
        "resolution": "384x384x768",
        "mesh_shape": (384, 384, 768),
        "meshblock_shape": (32, 32, 64),
        "mesh_bounds": ((0.0, 1.0), (0.0, 1.0), (0.0, 2.0)),
        "physics_contract": {
            "mhd/passive": "false",
            "mhd/cgl_heat_flux": "landau_fluid",
            "mhd/lf_k_parallel": "6.283185307179586",
            "mhd/limiter_nu_coll": "1.0e10",
            "mhd/limiter_hardwall": "true",
            "problem/passive_delta": "false",
            "problem/beta0": "10.0",
            "turb_driving/driving_type": "1",
            "turb_driving/projection_policy": "mks24_alfvenic_perpendicular",
        },
    },
}
_AUTHENTICATED_SELF: tuple[Path, Path, str] | None = None
for active_case in CASE_POLICIES:
    CASE_POLICIES[active_case].update({
        "operational_only": CASE_POLICIES[active_case].get("operational_only", False),
        "fixed_fresh_profile": CASE_POLICIES[active_case].get("fixed_fresh_profile"),
        "scientific_policy": "active_hardwall",
        "history_dt": 0.02,
        "snapshot_dt": 0.25,
        "snapshot_location_size_bytes": 8,
        "snapshot_variable_size_bytes": 4,
        "minimum_restart_bytes": 1024 * 1024,
        "required_snapshot_variables": REQUIRED_SNAPSHOT_VARIABLES,
    })

COMMON_INPUT_CONTRACT = {
    "mesh/nghost": "2",
    "mesh/x1min": "0.0",
    "mesh/x1max": "1.0",
    "mesh/ix1_bc": "periodic",
    "mesh/ox1_bc": "periodic",
    "mesh/x2min": "0.0",
    "mesh/x2max": "1.0",
    "mesh/ix2_bc": "periodic",
    "mesh/ox2_bc": "periodic",
    "mesh/x3min": "0.0",
    "mesh/x3max": "2.0",
    "mesh/ix3_bc": "periodic",
    "mesh/ox3_bc": "periodic",
    "time/evolution": "dynamic",
    "time/integrator": "rk2",
    "time/sts_integrator": "rkl2",
    "time/tlim": "10.0",
    "mhd/eos": "cgl",
    "mhd/cgl_heat_flux_integrator": "sts",
    "mhd/cgl_lf_strict_admissibility": "true",
    "mhd/cgl_lf_record_pressure_work": "true",
    "mhd/lf_coefficient_mode": "local",
    "mhd/mirror_limiter": "true",
    "mhd/firehose_limiter": "true",
    "mhd/cgl_firehose_threshold": "parallel",
    "mhd/backup_limiters": "false",
    "mhd/reconstruct": "plm",
    "mhd/rsolver": "hlle",
    "problem/pgen_name": "cgl_lf_paper",
    "problem/user_hist": "true",
    "problem/paper_mode": "turbulence",
    "turb_driving/dedt": "0.32",
    "turb_driving/tcorr": "2.0",
    "turb_driving/record_injected_work": "true",
    "output1/file_type": "hst",
    "output1/dt": "0.02",
    "output2/file_type": "bin",
    "output2/variable": "mhd_w_bcc",
    "output2/dt": "0.25",
    "output2/single_file_per_rank": "true",
    "output3/file_type": "rst",
    "output3/dt": "1.0",
    "output3/single_file_per_rank": "true",
}


def unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject duplicate JSON keys instead of silently choosing one value."""

    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"duplicate JSON key: {key}")
        value[key] = item
    return value


def reject_json_constant(value: str) -> object:
    """Reject non-finite JSON numbers."""

    raise ValueError(f"invalid JSON numeric constant: {value}")


def load_json(path: Path, label: str) -> dict[str, object]:
    """Load one unambiguous JSON object."""

    try:
        value = json.loads(
            read_regular_bytes(path, label, single_link=True).decode("utf-8"),
            object_pairs_hook=unique_json_object,
            parse_constant=reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is invalid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object: {path}")
    return value


def stable_json(value: object) -> str:
    """Return canonical, finite JSON with a trailing newline."""

    return json.dumps(
        value, indent=2, sort_keys=True, allow_nan=False
    ) + "\n"


def stable_json_sha256(value: object) -> str:
    """Digest one canonical compact JSON value."""

    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def retained_inventory_sha256(value: object) -> str:
    """Return the exact inventory digest consumed by Stage I recost."""

    payload = (json.dumps(value, sort_keys=True, allow_nan=False) + "\n").encode(
        "utf-8"
    )
    return hashlib.sha256(payload).hexdigest()


def stable_stat(profile: os.stat_result) -> tuple[int, int, int, int, int]:
    """Return the file identity and mutation fields used by authenticated reads."""

    return (
        profile.st_dev, profile.st_ino, profile.st_mode, profile.st_size,
        profile.st_mtime_ns,
    )


@contextmanager
def open_regular_fd(path: Path, label: str, *, single_link: bool = False):
    """Open one unchanged regular file without following its leaf symlink."""

    absolute = path.expanduser().absolute()
    try:
        before = absolute.lstat()
    except FileNotFoundError as error:
        raise ValueError(f"{label} is missing: {absolute}") from error
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ValueError(f"{label} must be a regular non-symlink file: {absolute}")
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(absolute, flags)
    try:
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or before.st_dev != opened.st_dev
            or before.st_ino != opened.st_ino
            or (single_link and opened.st_nlink != 1)
        ):
            raise ValueError(f"{label} changed while it was opened: {absolute}")
        yield descriptor
        after = os.fstat(descriptor)
        try:
            rebound = absolute.lstat()
        except FileNotFoundError as error:
            raise ValueError(f"{label} disappeared during inspection: {absolute}") from error
        if (
            stable_stat(opened) != stable_stat(after)
            or after.st_dev != rebound.st_dev
            or after.st_ino != rebound.st_ino
            or stat.S_ISLNK(rebound.st_mode)
            or (single_link and rebound.st_nlink != 1)
        ):
            raise ValueError(f"{label} changed during inspection: {absolute}")
    finally:
        os.close(descriptor)


def read_regular_bytes(path: Path, label: str, *, single_link: bool = False) -> bytes:
    """Read exact bytes through one authenticated descriptor."""

    with open_regular_fd(path, label, single_link=single_link) as descriptor:
        with os.fdopen(os.dup(descriptor), "rb") as stream:
            return stream.read()


def read_regular_profile(path: Path, label: str, *,
                         single_link: bool = False
                         ) -> tuple[bytes, dict[str, object]]:
    """Read and profile one file through the same authenticated descriptor."""

    with open_regular_fd(path, label, single_link=single_link) as descriptor:
        profile = os.fstat(descriptor)
        with os.fdopen(os.dup(descriptor), "rb") as stream:
            payload = stream.read()
    return payload, {
        "path": str(path.expanduser().absolute().resolve(strict=True)),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": profile.st_size,
    }


def regular_file_profile(path: Path, label: str, *,
                         single_link: bool = False) -> dict[str, object]:
    """Hash one unchanged regular file through a single descriptor."""

    digest = hashlib.sha256()
    with open_regular_fd(path, label, single_link=single_link) as descriptor:
        profile = os.fstat(descriptor)
        with os.fdopen(os.dup(descriptor), "rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    return {
        "path": str(path.expanduser().absolute().resolve(strict=True)),
        "sha256": digest.hexdigest(),
        "size_bytes": profile.st_size,
    }


def sha256(path: Path) -> str:
    """Return one unchanged regular file's SHA-256 digest."""

    return str(regular_file_profile(path, "digest input", single_link=True)["sha256"])


def sha256_descriptor(descriptor: int) -> str:
    """Return the SHA-256 digest of one open descriptor without changing its offset."""

    offset = os.lseek(descriptor, 0, os.SEEK_CUR)
    try:
        os.lseek(descriptor, 0, os.SEEK_SET)
        digest = hashlib.sha256()
        while block := os.read(descriptor, 1024 * 1024):
            digest.update(block)
        return digest.hexdigest()
    finally:
        os.lseek(descriptor, offset, os.SEEK_SET)


def require_regular_file(path: Path, label: str) -> None:
    """Require a directly named regular file, never a symlink."""

    try:
        profile = path.lstat()
    except FileNotFoundError as error:
        raise ValueError(f"{label} is missing: {path}") from error
    if stat.S_ISLNK(profile.st_mode) or not stat.S_ISREG(profile.st_mode):
        raise ValueError(f"{label} must be a regular non-symlink file: {path}")


def require_directory(path: Path, label: str) -> None:
    """Require a directly named directory, never a symlink."""

    try:
        profile = path.lstat()
    except FileNotFoundError as error:
        raise ValueError(f"{label} is missing: {path}") from error
    if stat.S_ISLNK(profile.st_mode) or not stat.S_ISDIR(profile.st_mode):
        raise ValueError(f"{label} must be a non-symlink directory: {path}")


def require_owned_executable_profile(profile: os.stat_result, label: str) -> None:
    """Require one immutable single-link executable owned by the effective user."""

    mode = stat.S_IMODE(profile.st_mode)
    if (
        not stat.S_ISREG(profile.st_mode)
        or profile.st_uid != os.geteuid()
        or profile.st_nlink != 1
        or mode & 0o022
        or mode & 0o111 == 0
    ):
        raise ValueError(f"{label} does not have the trusted owned executable profile")


def require_system_executable_profile(profile: os.stat_result, label: str) -> None:
    """Require one root-owned, single-link, non-writable executable."""

    mode = stat.S_IMODE(profile.st_mode)
    if (
        not stat.S_ISREG(profile.st_mode)
        or profile.st_uid != 0
        or profile.st_nlink != 1
        or mode & 0o022
        or mode & 0o111 == 0
    ):
        raise ValueError(f"{label} does not have the trusted system executable profile")


def require_system_directory_profile(path: Path, label: str) -> None:
    """Require one directly named root-owned, non-writable system directory."""

    require_directory(path, label)
    profile = path.lstat()
    if profile.st_uid != 0 or stat.S_IMODE(profile.st_mode) & 0o022:
        raise ValueError(f"{label} does not have the trusted system directory profile")


def require_system_directory_descriptor_profile(profile: os.stat_result,
                                                label: str) -> None:
    """Require one root-owned directory without unsafe replacement permissions."""

    if (
        not stat.S_ISDIR(profile.st_mode)
        or profile.st_uid != 0
        or stat.S_IMODE(profile.st_mode) & 0o022
    ):
        raise ValueError(f"{label} does not have the trusted system directory profile")


@contextmanager
def open_system_executable_fd(path: Path, label: str):
    """Open one root-owned absolute executable through trusted directory descriptors."""

    if (
        not path.is_absolute()
        or path != Path(os.path.normpath(str(path)))
        or ".." in path.parts
        or len(path.parts) < 2
    ):
        raise ValueError(f"{label} path is not normalized and absolute")
    directory_flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        directory_flags |= os.O_DIRECTORY
    if hasattr(os, "O_CLOEXEC"):
        directory_flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        directory_flags |= os.O_NOFOLLOW
    executable_flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        executable_flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        executable_flags |= os.O_NOFOLLOW
    parent = os.open("/", directory_flags)
    opened = [parent]
    try:
        require_system_directory_descriptor_profile(
            os.fstat(parent), f"{label} root directory"
        )
        for component in path.parts[1:-1]:
            parent = os.open(component, directory_flags, dir_fd=parent)
            opened.append(parent)
            require_system_directory_descriptor_profile(
                os.fstat(parent), f"{label} directory"
            )
        descriptor = os.open(path.parts[-1], executable_flags, dir_fd=parent)
        opened.append(descriptor)
        profile = os.fstat(descriptor)
        require_system_executable_profile(profile, label)
        require_bound_entry_identity(parent, path.parts[-1], profile, label)
        yield descriptor
        retained = os.fstat(descriptor)
        require_system_executable_profile(retained, label)
        if stable_stat(retained) != stable_stat(profile):
            raise ValueError(f"{label} changed during execution")
        require_bound_entry_identity(parent, path.parts[-1], retained, label)
    finally:
        for descriptor in reversed(opened):
            os.close(descriptor)


def require_no_symlink_chain(path: Path, stop: Path, label: str,
                             *, require_exists: bool = True) -> Path:
    """Resolve one path while rejecting every symlink between it and ``stop``."""

    absolute = path.expanduser().absolute()
    boundary = stop.expanduser().absolute()
    if boundary.exists() and boundary.resolve(strict=True) != boundary:
        raise ValueError(f"{label} boundary traverses a symlink: {boundary}")
    try:
        absolute.relative_to(boundary)
    except ValueError as error:
        raise ValueError(f"{label} must be beneath {boundary}: {absolute}") from error
    current = absolute
    while True:
        try:
            profile = current.lstat()
        except FileNotFoundError:
            pass
        else:
            if stat.S_ISLNK(profile.st_mode):
                raise ValueError(f"{label} traverses a symlink: {current}")
        if current == boundary:
            break
        current = current.parent
    if require_exists and not absolute.exists():
        raise ValueError(f"{label} is missing: {absolute}")
    resolved = absolute.resolve(strict=require_exists)
    if resolved != absolute:
        raise ValueError(f"{label} traverses a symlink: {absolute}")
    return absolute


def require_owned_regular_file(path: Path, parent: Path, label: str) -> Path:
    """Require one single-link regular file below a symlink-free parent."""

    resolved = require_no_symlink_chain(path, parent, label)
    require_regular_file(resolved, label)
    if resolved.lstat().st_nlink != 1:
        raise ValueError(f"{label} must have exactly one filesystem link")
    return resolved


def utc_now() -> str:
    """Return a stable UTC timestamp."""

    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def parse_utc_timestamp(value: object, label: str) -> datetime:
    """Parse one timezone-aware timestamp and normalize it to UTC."""

    try:
        parsed = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} is not an ISO timestamp") from error
    if parsed.tzinfo is None:
        raise ValueError(f"{label} lacks an explicit UTC offset")
    return parsed.astimezone(timezone.utc)


def format_utc_timestamp(value: datetime) -> str:
    """Format one timezone-aware timestamp as stable UTC."""

    if value.tzinfo is None:
        raise ValueError("UTC timestamp formatter requires an explicit offset")
    return value.astimezone(timezone.utc).isoformat(timespec="seconds")


def require_recent_scheduler_end(value: object, *,
                                 now: datetime | None = None) -> datetime:
    """Require completed scheduler evidence to remain fresh for selection."""

    end = parse_utc_timestamp(value, "scheduler end time")
    current = now or datetime.now(timezone.utc)
    if current.tzinfo is None:
        raise ValueError("freshness reference time lacks an explicit UTC offset")
    age = (current.astimezone(timezone.utc) - end).total_seconds()
    if age < -FUTURE_TIMESTAMP_TOLERANCE_SECONDS:
        raise ValueError("scheduler evidence end time is implausibly in the future")
    if age > EVIDENCE_MAX_AGE_SECONDS:
        raise ValueError("qualification scheduler evidence is stale")
    return end


def require_entry_name(name: str, label: str) -> str:
    """Require one direct-child basename for descriptor-relative mutation."""

    if not name or name in {".", ".."} or Path(name).name != name or "/" in name:
        raise ValueError(f"{label} must be one direct-child basename")
    return name


def require_directory_descriptor_binding(path: Path, descriptor: int,
                                         expected: os.stat_result,
                                         label: str) -> os.stat_result:
    """Require a directory pathname to retain the opened descriptor identity."""

    current = os.fstat(descriptor)
    try:
        named = path.lstat()
    except FileNotFoundError as error:
        raise ValueError(f"{label} path disappeared before mutation") from error
    if (
        not stat.S_ISDIR(current.st_mode)
        or stat.S_ISLNK(named.st_mode)
        or not stat.S_ISDIR(named.st_mode)
        or (current.st_dev, current.st_ino) != (expected.st_dev, expected.st_ino)
        or (named.st_dev, named.st_ino) != (expected.st_dev, expected.st_ino)
    ):
        raise ValueError(f"{label} path changed before mutation")
    return current


def require_bound_entry_identity(parent: int, name: str,
                                 expected: os.stat_result,
                                 label: str) -> os.stat_result:
    """Require one descriptor-relative entry to retain its inode identity."""

    name = require_entry_name(name, label)
    try:
        named = os.stat(name, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError as error:
        raise ValueError(f"{label} disappeared during mutation") from error
    if (
        stat.S_ISLNK(named.st_mode)
        or (named.st_dev, named.st_ino) != (expected.st_dev, expected.st_ino)
    ):
        raise ValueError(f"{label} inode identity changed during mutation")
    return named


@contextmanager
def open_directory_fd(path: Path, label: str):
    """Yield one pathname-bound directory descriptor."""

    absolute = path.expanduser().absolute()
    try:
        before = absolute.lstat()
    except FileNotFoundError as error:
        raise ValueError(f"{label} is missing: {absolute}") from error
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
        raise ValueError(f"{label} must be a non-symlink directory: {absolute}")
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(absolute, flags)
    try:
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISDIR(opened.st_mode)
            or (opened.st_dev, opened.st_ino) != (before.st_dev, before.st_ino)
        ):
            raise ValueError(f"{label} changed while it was opened: {absolute}")
        yield descriptor, opened
        require_directory_descriptor_binding(absolute, descriptor, opened, label)
    finally:
        os.close(descriptor)


def write_exclusive(path: Path, payload: bytes, mode: int = 0o640,
                    mutation_guard=None) -> None:
    """Create one descriptor-bound retained artifact without replacing data."""

    path = path.expanduser().absolute()
    require_no_symlink_chain(path.parent, path.parent.parent, "artifact parent")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = None
    with open_directory_fd(path.parent, "artifact parent") as (
        parent_descriptor, parent_profile,
    ):
        require_directory_descriptor_binding(
            path.parent, parent_descriptor, parent_profile, "artifact parent"
        )
        if mutation_guard is not None:
            mutation_guard()
        require_directory_descriptor_binding(
            path.parent, parent_descriptor, parent_profile, "artifact parent"
        )
        try:
            descriptor = os.open(
                require_entry_name(path.name, "artifact"),
                flags,
                mode,
                dir_fd=parent_descriptor,
            )
            opened = os.fstat(descriptor)
            if not stat.S_ISREG(opened.st_mode) or opened.st_nlink != 1:
                raise ValueError(f"artifact does not have one regular link: {path}")
            require_bound_entry_identity(
                parent_descriptor, path.name, opened, "artifact"
            )
            if mutation_guard is not None:
                mutation_guard()
            require_directory_descriptor_binding(
                path.parent, parent_descriptor, parent_profile, "artifact parent"
            )
            os.fchmod(descriptor, mode)
            with os.fdopen(descriptor, "wb", closefd=False) as stream:
                stream.write(payload)
                stream.flush()
                os.fsync(stream.fileno())
            retained = os.fstat(descriptor)
            require_bound_entry_identity(
                parent_descriptor, path.name, retained, "artifact"
            )
            if (
                not stat.S_ISREG(retained.st_mode)
                or retained.st_nlink != 1
                or stat.S_IMODE(retained.st_mode) != mode
                or retained.st_size != len(payload)
            ):
                raise ValueError(f"artifact profile changed during mutation: {path}")
            if mutation_guard is not None:
                mutation_guard()
            require_directory_descriptor_binding(
                path.parent, parent_descriptor, parent_profile, "artifact parent"
            )
            os.fsync(parent_descriptor)
        finally:
            if descriptor is not None:
                os.close(descriptor)


def mkdir_exclusive(path: Path, mode: int = 0o750, mutation_guard=None) -> None:
    """Create one descriptor-bound directory without replacing an entry."""

    path = path.expanduser().absolute()
    require_no_symlink_chain(
        path.parent, path.parent.parent, "directory creation parent"
    )
    with open_directory_fd(path.parent, "directory creation parent") as (
        parent_descriptor, parent_profile,
    ):
        if mutation_guard is not None:
            mutation_guard()
        require_directory_descriptor_binding(
            path.parent, parent_descriptor, parent_profile,
            "directory creation parent",
        )
        os.mkdir(
            require_entry_name(path.name, "created directory"),
            mode,
            dir_fd=parent_descriptor,
        )
        flags = os.O_RDONLY
        if hasattr(os, "O_DIRECTORY"):
            flags |= os.O_DIRECTORY
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(path.name, flags, dir_fd=parent_descriptor)
        try:
            opened = os.fstat(descriptor)
            require_bound_entry_identity(
                parent_descriptor, path.name, opened, "created directory"
            )
            if not stat.S_ISDIR(opened.st_mode):
                raise ValueError(f"created directory profile differs: {path}")
            if mutation_guard is not None:
                mutation_guard()
            os.fchmod(descriptor, mode)
            os.fsync(descriptor)
            require_bound_entry_identity(
                parent_descriptor, path.name, os.fstat(descriptor),
                "created directory",
            )
            require_directory_descriptor_binding(
                path.parent, parent_descriptor, parent_profile,
                "directory creation parent",
            )
            os.fsync(parent_descriptor)
        finally:
            os.close(descriptor)


def chmod_directory(path: Path, mode: int, mutation_guard=None) -> None:
    """Change one descriptor-bound directory mode."""

    path = path.expanduser().absolute()
    require_no_symlink_chain(path.parent, path.parent.parent, "directory mode parent")
    with open_directory_fd(path.parent, "directory mode parent") as (
        parent_descriptor, parent_profile,
    ):
        flags = os.O_RDONLY
        if hasattr(os, "O_DIRECTORY"):
            flags |= os.O_DIRECTORY
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(
            require_entry_name(path.name, "directory mode target"),
            flags,
            dir_fd=parent_descriptor,
        )
        try:
            opened = os.fstat(descriptor)
            require_bound_entry_identity(
                parent_descriptor, path.name, opened, "directory mode target"
            )
            if mutation_guard is not None:
                mutation_guard()
            require_directory_descriptor_binding(
                path.parent, parent_descriptor, parent_profile,
                "directory mode parent",
            )
            os.fchmod(descriptor, mode)
            os.fsync(descriptor)
            retained = require_bound_entry_identity(
                parent_descriptor, path.name, os.fstat(descriptor),
                "directory mode target",
            )
            if not stat.S_ISDIR(retained.st_mode) or stat.S_IMODE(retained.st_mode) != mode:
                raise ValueError(f"directory mode mutation differs: {path}")
            if mutation_guard is not None:
                mutation_guard()
            os.fsync(parent_descriptor)
        finally:
            os.close(descriptor)


def write_json_exclusive(path: Path, value: object, mutation_guard=None) -> None:
    """Retain one immutable deterministic JSON artifact."""

    write_exclusive(
        path, stable_json(value).encode("utf-8"), mutation_guard=mutation_guard
    )


def copy_exclusive(source: Path, destination: Path, mode: int = 0o640,
                   mutation_guard=None) -> None:
    """Retain one exact file copy without replacing or following destinations."""

    write_exclusive(
        destination,
        read_regular_bytes(source, "copy source", single_link=True),
        mode,
        mutation_guard=mutation_guard,
    )


def require_exact_keys(value: object, keys: set[str], label: str) -> dict[str, object]:
    """Require one object with no missing or unreviewed columns."""

    if not isinstance(value, dict) or set(value) != keys:
        raise ValueError(f"{label} has invalid columns")
    return value


def require_sha256(value: object, label: str) -> str:
    """Require one lowercase SHA-256 digest."""

    if not isinstance(value, str) or SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return value


def require_nonempty_string(value: object, label: str) -> str:
    """Require one nonempty trimmed string."""

    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError(f"{label} must be a nonempty trimmed string")
    return value


def require_git_revision(value: object, label: str) -> str:
    """Require one full lowercase Git revision."""

    if not isinstance(value, str) or GIT_REVISION_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{label} must be a full lowercase Git revision")
    return value


def require_finite_positive(value: object, label: str) -> float:
    """Require one positive finite numeric value."""

    if isinstance(value, bool):
        raise ValueError(f"{label} must be positive and finite")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be positive and finite") from error
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{label} must be positive and finite")
    return number


def require_project_root(root_value: str | Path, allow_local_root: bool) -> Path:
    """Require the canonical project root except for explicit offline tests."""

    root_candidate = Path(root_value).expanduser().absolute()
    root = root_candidate.resolve(strict=True)
    require_directory(root, "project root")
    if root_candidate != root:
        raise ValueError("project root may not traverse a symlink")
    if root != DEFAULT_ROOT.resolve() and not allow_local_root:
        raise ValueError(
            f"Stage I qualification project root must be {DEFAULT_ROOT}; "
            "use --allow-local-root only for offline validation"
        )
    return root


def require_qualification_root(root: Path, value: str | Path) -> Path:
    """Confine every pilot product to one qualification-only direct child."""

    runs = root / "runs"
    require_directory(runs, "project runs directory")
    candidate = Path(value).expanduser().absolute()
    if not candidate.is_absolute():
        raise ValueError("qualification root must be an absolute path")
    if candidate.exists() or candidate.is_symlink():
        require_directory(candidate, "qualification root")
    qualification = require_no_symlink_chain(
        candidate, runs.resolve(strict=True), "qualification root",
        require_exists=candidate.exists(),
    )
    if qualification.parent != runs.resolve(strict=True):
        raise ValueError(
            "qualification root must be a direct child of the project runs directory"
        )
    if QUALIFICATION_ROOT_PATTERN.fullmatch(qualification.name) is None:
        raise ValueError("qualification root must use the qualification-* namespace")
    canonical = (runs / "mks24-stage-i").resolve(strict=False)
    if qualification == canonical or canonical in qualification.parents:
        raise ValueError("qualification output may not enter canonical Stage I lineages")
    return qualification


def require_beneath(path: Path, parent: Path, label: str) -> None:
    """Require one resolved path to remain below a reviewed parent."""

    try:
        path.relative_to(parent)
    except ValueError as error:
        raise ValueError(f"{label} must be beneath {parent}: {path}") from error


def require_qualification_path(path_value: str | Path, qualification_root: Path,
                               label: str, must_exist: bool = False) -> Path:
    """Require one path below the isolated qualification root."""

    candidate = Path(path_value).expanduser().absolute()
    if must_exist:
        path = require_owned_regular_file(candidate, qualification_root, label)
    else:
        path = require_no_symlink_chain(
            candidate, qualification_root, label, require_exists=False
        )
    if path == qualification_root:
        raise ValueError(f"{label} may not equal the qualification root")
    return path


def parse_walltime(value: str) -> int:
    """Parse one HH:MM:SS walltime."""

    match = re.fullmatch(r"([0-9]{2}):([0-5][0-9]):([0-5][0-9])", value)
    if match is None:
        raise ValueError(f"walltime must use HH:MM:SS: {value}")
    hours, minutes, seconds = (int(item) for item in match.groups())
    return hours * 3600 + minutes * 60 + seconds


def require_r12_fresh_profile_policy() -> dict[str, object]:
    """Require the exact measured parentless fresh R12 qualification profile."""

    policy = CASE_POLICIES[R12_CASE_ID]
    profile = require_exact_keys(
        policy.get("fixed_fresh_profile"),
        {
            "case_id", "segment", "start_time", "target_time", "nodes",
            "ranks_per_node", "total_ranks", "walltime", "athena_walltime",
            "parent_job_id", "parent_result", "parent_segment", "restart_file",
            "restart_file_sha256", "restart_time",
        },
        "R12 fixed fresh qualification profile",
    )
    null_lineage = (
        "parent_job_id", "parent_result", "parent_segment", "restart_file",
        "restart_file_sha256", "restart_time",
    )
    if (
        profile != R12_FRESH_RERUN_PROFILE
        or profile["case_id"] != R12_CASE_ID
        or profile["segment"] != "s01_rankio_t0_t0p12"
        or profile["start_time"] != 0.0
        or profile["target_time"] != 0.12
        or profile["nodes"] != 4
        or profile["ranks_per_node"] != RANKS_PER_NODE
        or profile["nodes"] * profile["ranks_per_node"] != profile["total_ranks"]
        or profile["walltime"] != "02:00:00"
        or profile["athena_walltime"] != "01:50:00"
        or any(profile[key] is not None for key in null_lineage)
        or policy["node_profiles"] != (profile["nodes"],)
        or policy["target_time"] != profile["target_time"]
        or policy["walltime"] != profile["walltime"]
        or policy["athena_walltime"] != profile["athena_walltime"]
        or policy["history_dt"] != 0.02
        or policy["snapshot_dt"] != 0.25
        or COMMON_INPUT_CONTRACT["output1/dt"] != "0.02"
        or COMMON_INPUT_CONTRACT["output2/dt"] != "0.25"
        or COMMON_INPUT_CONTRACT["output3/dt"] != "1.0"
    ):
        raise ValueError(
            "R12 qualification policy differs from the measured parentless "
            "s01_rankio_t0_t0p12 4-node/32-rank fresh profile"
        )
    return profile


def retained_file(path: Path, label: str) -> dict[str, object]:
    """Return deterministic provenance for one regular file."""

    return regular_file_profile(path, label, single_link=True)


def root_relative_binding(path: Path, root: Path, label: str) -> dict[str, str]:
    """Return one recost-compatible root-relative immutable input binding."""

    candidate = path.expanduser().absolute()
    path = require_owned_regular_file(candidate, root, label)
    require_beneath(path, root, label)
    return {
        "path": path.relative_to(root).as_posix(),
        "sha256": sha256(path),
    }


def root_relative_payload_binding(path: Path, payload: bytes, root: Path,
                                  label: str) -> dict[str, str]:
    """Bind not-yet-published canonical bytes to one root-relative path."""

    candidate = path.expanduser().absolute()
    require_beneath(candidate, root, label)
    if candidate == root:
        raise ValueError(f"{label} may not equal the project root")
    return {
        "path": candidate.relative_to(root).as_posix(),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def build_manifest_inventory(manifest: Path) -> list[dict[str, str]]:
    """Authenticate the complete flat retained build-manifest inventory."""

    require_directory(manifest, "build manifest")
    entries = []
    names = sorted(item.name for item in manifest.iterdir())
    if not names:
        raise ValueError("build manifest must not be empty")
    for name in names:
        path = require_owned_regular_file(
            manifest / name, manifest, f"build-manifest file {name}"
        )
        mode = stat.S_IMODE(path.stat().st_mode)
        if mode != 0o644:
            raise ValueError(f"build-manifest file mode differs: {path}")
        entries.append({
            "name": name,
            "mode": f"{mode:04o}",
            "sha256": sha256(path),
        })
    return entries


def validate_retained_file(record: object, qualification_root: Path,
                           label: str) -> Path:
    """Revalidate one qualification evidence file and its exact digest."""

    value = require_exact_keys(
        record, {"path", "sha256", "size_bytes"}, label
    )
    path = require_qualification_path(
        str(value["path"]), qualification_root, label, must_exist=True
    )
    require_sha256(value["sha256"], f"{label} sha256")
    if isinstance(value["size_bytes"], bool):
        raise ValueError(f"{label} size is invalid")
    try:
        size = int(value["size_bytes"])
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} size is invalid") from error
    observed = regular_file_profile(path, label, single_link=True)
    if size < 0 or observed != value:
        raise ValueError(f"{label} retained file binding differs: {path}")
    return path


def output_file_inventory(output: Path, qualification_root: Path
                          ) -> list[dict[str, object]]:
    """Hash the complete regular-file inventory below one pilot output root."""

    output = require_no_symlink_chain(output, qualification_root, "pilot output directory")
    require_directory(output, "pilot output directory")
    records = []
    for current, directories, names in os.walk(output, followlinks=False):
        directories.sort()
        names.sort()
        current_path = Path(current)
        for name in directories:
            directory = current_path / name
            require_no_symlink_chain(directory, qualification_root, "pilot output directory")
            require_directory(directory, "pilot output directory")
        for name in names:
            path = require_owned_regular_file(
                current_path / name, qualification_root, "pilot output binding file"
            )
            profile = retained_file(path, "pilot output binding file")
            records.append({
                "path": path.relative_to(output).as_posix(),
                "sha256": profile["sha256"],
                "size_bytes": profile["size_bytes"],
            })
    if not records:
        raise ValueError("pilot output binding inventory is empty")
    return records


def output_binding_record(packet: dict[str, object], qualification_root: Path,
                          provenance: dict[str, object], *, job_id: str,
                          completed_utc: str) -> dict[str, object]:
    """Build the exact end-of-job output-tree binding for one pilot."""

    if JOB_ID_PATTERN.fullmatch(job_id) is None:
        raise ValueError("output binding job ID is invalid")
    parse_utc_timestamp(completed_utc, "output binding completion time")
    intent = packet["execution_intent"]
    executable = require_exact_keys(
        provenance.get("executable"),
        {"path", "sha256", "size_bytes", "revision"},
        "qualification executable provenance",
    )
    output = require_no_symlink_chain(
        Path(str(intent["paths"]["output_dir"])),
        qualification_root,
        "pilot output directory",
    )
    files = output_file_inventory(output, qualification_root)
    return {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_qualification_output_binding",
        "job_id": job_id,
        "execution_contract_sha256": intent["execution_contract_sha256"],
        "executable_sha256": executable["sha256"],
        "completed_utc": completed_utc,
        "output_dir": str(output),
        "files": files,
        "output_tree_sha256": stable_json_sha256(files),
    }


def validate_output_binding(path: Path, *, packet: dict[str, object],
                            binding: dict[str, object],
                            scheduler: dict[str, object],
                            qualification_root: Path,
                            provenance: dict[str, object]) -> dict[str, object]:
    """Rehash every output byte and bind it to one authenticated Slurm interval."""

    expected_path = Path(str(packet["execution_intent"]["paths"]["output_binding"]))
    path = require_owned_regular_file(path, qualification_root, "pilot output binding")
    if path != expected_path.resolve(strict=True):
        raise ValueError("pilot output binding path differs from authenticated packet")
    if path.stat().st_mode & 0o222:
        raise ValueError("pilot output binding must be immutable")
    retained = require_exact_keys(
        load_json(path, "pilot output binding"),
        {
            "schema_version", "record_type", "job_id", "execution_contract_sha256",
            "executable_sha256", "completed_utc", "output_dir", "files",
            "output_tree_sha256",
        },
        "pilot output binding",
    )
    completed = parse_utc_timestamp(
        retained["completed_utc"], "output binding completion time"
    )
    started = parse_scheduler_time(str(scheduler["start_utc"]), "scheduler start")
    ended = parse_scheduler_time(str(scheduler["end_utc"]), "scheduler end")
    if not started <= completed <= ended:
        raise ValueError("pilot output binding completion falls outside Slurm job")
    reproduced = output_binding_record(
        packet,
        qualification_root,
        provenance,
        job_id=str(binding["job_id"]),
        completed_utc=str(retained["completed_utc"]),
    )
    if retained != reproduced:
        raise ValueError("pilot output binding differs from live output tree")
    return retained


def validate_live_file_record(record: object, label: str) -> Path:
    """Revalidate one exact retained file record without a path namespace."""

    value = require_exact_keys(
        record, {"path", "sha256", "size_bytes"}, label
    )
    candidate = Path(str(value["path"])).expanduser().absolute()
    require_regular_file(candidate, label)
    path = candidate.resolve(strict=True)
    if candidate != path:
        raise ValueError(f"{label} may not traverse a symlink: {candidate}")
    require_sha256(value["sha256"], f"{label} sha256")
    if isinstance(value["size_bytes"], bool):
        raise ValueError(f"{label} size is invalid")
    try:
        size = int(value["size_bytes"])
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} size is invalid") from error
    observed = regular_file_profile(path, label, single_link=True)
    if size < 0 or observed != value:
        raise ValueError(f"{label} retained file binding differs: {path}")
    return path


def validate_prepared_provenance(record: object, root: Path,
                                 selected_paths: list[Path]) -> None:
    """Require complete live provenance for one retained qualification wave."""

    provenance = require_exact_keys(
        record,
        {
            "source", "source_bundle", "matrix", "executable",
            "build_manifest", "qualification_helper",
        },
        "prepared qualification provenance",
    )
    source = require_exact_keys(
        provenance["source"], {"directory", "revision"}, "source provenance"
    )
    source_candidate = Path(str(source["directory"]))
    require_directory(source_candidate, "source provenance directory")
    source_directory = source_candidate.resolve(strict=True)
    if source_candidate.absolute() != source_directory:
        raise ValueError("prepared source directory may not traverse a symlink")
    source_revision = require_git_revision(
        source["revision"], "source provenance revision"
    )
    if source_revision != FROZEN_SOURCE_REVISION:
        raise ValueError("prepared source revision is not the reviewed frozen revision")
    matrix_path = validate_live_file_record(provenance["matrix"], "prepared matrix")
    if provenance["matrix"]["sha256"] != FROZEN_MATRIX_SHA256:
        raise ValueError("prepared matrix is not the reviewed frozen matrix")

    bundle = require_exact_keys(
        provenance["source_bundle"],
        {"path", "sha256", "size_bytes", "verified_revisions"},
        "prepared source bundle",
    )
    bundle_path = validate_live_file_record(
        {key: bundle[key] for key in ("path", "sha256", "size_bytes")},
        "prepared source bundle",
    )
    source_archives = root / "source-archives"
    require_directory(source_archives, "source-archive directory")
    require_beneath(
        bundle_path, source_archives.resolve(strict=True),
        "prepared source bundle",
    )

    executable = require_exact_keys(
        provenance["executable"],
        {"path", "sha256", "size_bytes", "revision"},
        "prepared executable",
    )
    executable_path = validate_live_file_record(
        {key: executable[key] for key in ("path", "sha256", "size_bytes")},
        "prepared executable",
    )
    require_beneath(executable_path, root, "prepared executable")
    executable_revision = require_git_revision(
        executable["revision"], "prepared executable revision"
    )
    if executable_revision != source_revision:
        raise ValueError("prepared executable revision differs from frozen source revision")

    helper = require_exact_keys(
        provenance["qualification_helper"],
        {"path", "sha256", "size_bytes", "revision", "committed"},
        "prepared qualification helper",
    )
    validate_live_file_record(
        {key: helper[key] for key in ("path", "sha256", "size_bytes")},
        "prepared qualification helper",
    )
    helper_revision = require_git_revision(
        helper["revision"], "prepared qualification helper revision"
    )
    if helper["committed"] is not True:
        raise ValueError("prepared qualification helper was not committed")

    revisions = bundle["verified_revisions"]
    expected_revisions = sorted({
        source_revision, executable_revision, helper_revision,
    })
    if revisions != expected_revisions:
        raise ValueError("prepared source bundle revisions are incomplete or ambiguous")

    build = require_exact_keys(
        provenance["build_manifest"],
        {
            "path", "athena_sha256", "environment", "inventory",
            "inventory_sha256",
        },
        "prepared build manifest",
    )
    build_candidate = Path(str(build["path"]))
    require_directory(build_candidate, "prepared build manifest")
    build_path = build_candidate.resolve(strict=True)
    if build_candidate.absolute() != build_path:
        raise ValueError("prepared build manifest may not traverse a symlink")
    require_beneath(build_path, root, "prepared build manifest")
    for key in ("athena_sha256", "environment"):
        path = validate_live_file_record(
            build[key], f"prepared build manifest {key}"
        )
        if path.parent != build_path:
            raise ValueError("prepared build-manifest files have inconsistent parents")
    inventory = build["inventory"]
    if (
        not isinstance(inventory, list)
        or inventory != build_manifest_inventory(build_path)
        or require_sha256(
            build["inventory_sha256"], "prepared build-manifest inventory digest"
        )
        != retained_inventory_sha256(inventory)
    ):
        raise ValueError("prepared build-manifest inventory differs from live build")
    if matrix_path != Path(str(provenance["matrix"]["path"])).resolve(strict=True):
        raise ValueError("prepared matrix path is ambiguous")
    observed_revision = committed_source_revision(
        source_directory, [matrix_path, *selected_paths]
    )
    if observed_revision != source_revision:
        raise ValueError("prepared source revision differs from live committed source")
    observed_helper = utility_provenance()
    if observed_helper != helper:
        raise ValueError("prepared qualification helper provenance differs from live helper")
    observed_build = build_provenance(executable_path, build_path)
    if observed_build["executable"] != executable:
        raise ValueError("prepared executable provenance differs from live build")
    if observed_build["build_manifest"] != build:
        raise ValueError("prepared build-manifest provenance differs from live build")
    observed_bundle = source_bundle_provenance(
        bundle_path, root, expected_revisions
    )
    if observed_bundle != bundle:
        raise ValueError("prepared source-bundle provenance differs from live bundle")


def hardened_child_environment(*forbidden_prefixes: str) -> dict[str, str]:
    """Strip private reexec, loader, interpreter, Git, and caller PATH controls."""

    private = {SELF_DESCRIPTOR_ENV, SELF_SOURCE_ENV, REPOSITORY_ROOT_ENV}
    forbidden_names = {
        "BASH_ENV", "BASHOPTS", "CDPATH", "ENV", "GLOBIGNORE", "IFS",
        "PROMPT_COMMAND", "SHELLOPTS",
    }
    environment = {
        key: value
        for key, value in os.environ.items()
        if key not in private
        and key not in forbidden_names
        and not key.startswith("_CGL_LF_")
        and not key.startswith("BASH_FUNC_")
        and not key.startswith("GIT_")
        and not key.startswith("LD_")
        and not key.startswith("PERL")
        and not key.startswith("PYTHON")
        and not key.startswith("RUBY")
        and not any(key.startswith(prefix) for prefix in forbidden_prefixes)
    }
    environment.update({"LC_ALL": "C", "PATH": TRUSTED_SYSTEM_PATH})
    return environment


def scheduler_environment(*, exact_timestamps: bool = False) -> dict[str, str]:
    """Return caller-independent routing and loader state for one Slurm child."""

    environment = hardened_child_environment("SBATCH_", "SLURM_")
    if exact_timestamps:
        environment["SLURM_TIME_FORMAT"] = SLURM_TIME_FORMAT
    return environment


def hardened_git_environment() -> dict[str, str]:
    """Return Git's complete caller-independent execution environment."""

    return {
        "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_CONFIG_SYSTEM": os.devnull,
        "GIT_EXEC_PATH": str(GIT_EXEC_PATH),
        "GIT_OPTIONAL_LOCKS": "0",
        "GIT_TERMINAL_PROMPT": "0",
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": TRUSTED_SYSTEM_PATH,
    }


def hardened_git_arguments(arguments: list[str]) -> list[str]:
    """Disable repository-local execution surfaces for one exact Git query."""

    return [
        "-c", "core.fsmonitor=false",
        "-c", f"core.hooksPath={os.devnull}",
        "-c", f"core.attributesFile={os.devnull}",
        "-c", f"core.excludesFile={os.devnull}",
        "-c", f"init.templateDir={os.devnull}",
        "-c", "init.defaultObjectFormat=sha1",
        *arguments,
    ]


def git_run(arguments: list[str], *, capture_output: bool = False,
            text: bool = False,
            pass_fds: tuple[int, ...] = ()) -> subprocess.CompletedProcess:
    """Run one authenticated root-owned Git descriptor with isolated configuration."""

    with open_system_executable_fd(
        GIT, "authenticated absolute Git binary"
    ) as git:
        try:
            return subprocess.run(
                [str(GIT), "--no-replace-objects", *hardened_git_arguments(arguments)],
                executable=f"/proc/self/fd/{git}",
                check=False,
                stdin=subprocess.DEVNULL,
                capture_output=capture_output,
                text=text,
                env=hardened_git_environment(),
                pass_fds=(git, *pass_fds),
                timeout=120,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            raise ValueError("authenticated absolute Git query failed to execute") from error


def git_output(repository: Path, arguments: list[str], label: str) -> str:
    """Run one read-only authenticated Git query and return stripped stdout."""

    result = git_run(
        ["-C", str(repository), *arguments], capture_output=True, text=True
    )
    if result.returncode != 0:
        raise ValueError(f"cannot determine {label}")
    return result.stdout.strip()


def committed_source_revision(source_dir: Path, paths: list[Path]) -> str:
    """Require selected source inputs and matrix to be tracked and committed."""

    require_directory(source_dir, "source directory")
    revision = require_git_revision(
        git_output(
            source_dir, ["rev-parse", "--verify", "HEAD"],
            "source revision",
        ),
        "source revision",
    )
    relative_paths = []
    before_profiles = []
    for path in paths:
        resolved = path.resolve(strict=True)
        if path.expanduser().absolute() != resolved:
            raise ValueError(f"source input may not traverse a symlink: {path}")
        require_beneath(resolved, source_dir, "source input")
        relative_paths.append(str(resolved.relative_to(source_dir)))
        before_profiles.append(
            regular_file_profile(resolved, "committed source input", single_link=True)
        )
    for relative in relative_paths:
        tracked = git_run(
            [
                "-C", str(source_dir), "ls-files", "--error-unmatch",
                "--", relative,
            ],
            capture_output=True,
            text=True,
        )
        if tracked.returncode != 0:
            raise ValueError(f"source provenance is untracked: {relative}")
    for diff_arguments in (
        ["diff", "--no-ext-diff", "--no-textconv", "--quiet", "--"],
        ["diff", "--cached", "--no-ext-diff", "--no-textconv", "--quiet", "--"],
    ):
        result = git_run(
            ["-C", str(source_dir), *diff_arguments, *relative_paths],
        )
        if result.returncode == 1:
            raise ValueError("selected source inputs and matrix must be committed")
        if result.returncode != 0:
            raise ValueError("cannot determine selected source worktree status")
    final_revision = require_git_revision(
        git_output(
            source_dir, ["rev-parse", "--verify", "HEAD"],
            "source revision",
        ),
        "source revision",
    )
    after_profiles = [
        regular_file_profile(path.resolve(strict=True), "committed source input",
                             single_link=True)
        for path in paths
    ]
    if final_revision != revision or after_profiles != before_profiles:
        raise ValueError("selected source inputs or revision changed during inspection")
    return revision


def initial_self_source_path() -> Path:
    """Return the named qualification source for initial descriptor reexecution."""

    source = Path(__file__).absolute()
    if source.is_symlink():
        raise ValueError("qualification utility source must not be a symlink")
    return source.resolve(strict=True)


def require_reexec_source_relationship(source: Path, repository: Path) -> None:
    """Bind private reexec metadata to the exact qualification utility path."""

    if source != repository / UTILITY_RELATIVE_PATH:
        raise ValueError("authenticated qualification source/repository relationship differs")


def inherited_reexec_path(name: str, label: str) -> Path:
    """Read one private reexec path only after authenticating the self descriptor."""

    retained = os.environ.get(name)
    if retained is None:
        raise ValueError(f"authenticated qualification reexecution lacks {label}")
    path = Path(retained)
    if (
        not path.is_absolute()
        or path != Path(os.path.normpath(retained))
        or ".." in path.parts
        or len(path.parts) < 2
    ):
        raise ValueError(
            f"authenticated qualification {label} is not normalized and absolute"
        )
    require_no_symlink_chain(path, Path(path.anchor), f"authenticated qualification {label}")
    return path


def reexec_environment(descriptor: int, source: Path,
                       repository: Path) -> dict[str, str]:
    """Return a qualification reexec environment without caller execution controls."""

    environment = hardened_child_environment("SBATCH_", "SLURM_")
    environment.update(
        {
            SELF_DESCRIPTOR_ENV: str(descriptor),
            SELF_SOURCE_ENV: str(source),
            REPOSITORY_ROOT_ENV: str(repository),
            "HOME": "/nonexistent",
        }
    )
    return environment


def reexec_authenticated_self(descriptor: int, source: Path, repository: Path,
                              argv: list[str]) -> None:
    """Reexecute source bytes through an authenticated isolated system Python."""

    with open_system_executable_fd(
        SYSTEM_PYTHON, "authenticated absolute system Python"
    ) as python_descriptor:
        os.set_inheritable(descriptor, True)
        os.set_inheritable(python_descriptor, True)
        os.execve(
            f"/proc/self/fd/{python_descriptor}",
            [
                str(SYSTEM_PYTHON), "-I", "-S", "-B",
                f"/proc/self/fd/{descriptor}", *argv,
            ],
            reexec_environment(descriptor, source, repository),
        )
    raise AssertionError("authenticated qualification reexecution unexpectedly returned")


def require_authenticated_reexec_runtime() -> None:
    """Require the inherited process to be the isolated authenticated Python."""

    if os.environ.get("HOME") != "/nonexistent":
        raise ValueError("authenticated qualification reexecution HOME is not sanitized")
    flags = sys.flags
    if (
        flags.isolated != 1
        or flags.ignore_environment != 1
        or flags.no_user_site != 1
        or flags.no_site != 1
        or flags.dont_write_bytecode != 1
    ):
        raise ValueError("authenticated qualification reexecution is not isolated")
    with open_system_executable_fd(
        SYSTEM_PYTHON, "authenticated absolute system Python"
    ) as descriptor:
        expected = os.fstat(descriptor)
        running = os.stat("/proc/self/exe")
        if (
            not stat.S_ISREG(running.st_mode)
            or (running.st_dev, running.st_ino) != (expected.st_dev, expected.st_ino)
        ):
            raise ValueError(
                "authenticated qualification reexecution interpreter differs"
            )


def authenticate_committed_self(source: Path, repository: Path,
                                expected_sha256: str) -> str:
    """Bind immutable executing bytes to the committed qualification utility."""

    revision = committed_source_revision(repository, [source])
    retained = git_run(
        [
            "-C", str(repository), "show",
            f"{revision}:{UTILITY_RELATIVE_PATH.as_posix()}",
        ],
        capture_output=True,
    )
    if (
        retained.returncode != 0
        or hashlib.sha256(retained.stdout).hexdigest() != expected_sha256
    ):
        raise ValueError("committed qualification utility bytes differ from execution")
    return revision


def authenticate_self(argv: list[str]) -> tuple[Path, Path, str]:
    """Reexecute immutable committed qualification bytes and retain their identity."""

    global _AUTHENTICATED_SELF
    inherited = os.environ.get(SELF_DESCRIPTOR_ENV)
    if inherited is not None:
        if re.fullmatch(r"[0-9]+", inherited) is None:
            raise ValueError("qualification descriptor marker is invalid")
        descriptor = int(inherited)
        if __file__ != f"/proc/self/fd/{descriptor}":
            raise ValueError(
                "qualification descriptor marker is not attached to this execution"
            )
        require_owned_executable_profile(
            os.fstat(descriptor), "authenticated qualification descriptor"
        )
        require_authenticated_reexec_runtime()
        source = inherited_reexec_path(SELF_SOURCE_ENV, "source path")
        repository = inherited_reexec_path(REPOSITORY_ROOT_ENV, "repository root")
        require_reexec_source_relationship(source, repository)
        digest = sha256_descriptor(descriptor)
        authenticate_committed_self(source, repository, digest)
        _AUTHENTICATED_SELF = (source, repository, digest)
        return _AUTHENTICATED_SELF

    orphaned = [
        name for name in (SELF_SOURCE_ENV, REPOSITORY_ROOT_ENV) if name in os.environ
    ]
    if orphaned:
        raise ValueError(
            "private qualification reexecution path is forbidden without an "
            "authenticated descriptor"
        )
    source = initial_self_source_path()
    repository = source.parents[2].resolve(strict=True)
    require_reexec_source_relationship(source, repository)
    with open_regular_fd(
        source, "retained qualification utility", single_link=True
    ) as descriptor:
        require_owned_executable_profile(
            os.fstat(descriptor), "retained qualification utility"
        )
        digest = sha256_descriptor(descriptor)
        authenticate_committed_self(source, repository, digest)
        reexec_authenticated_self(descriptor, source, repository, argv)
    raise AssertionError("qualification descriptor reexecution unexpectedly returned")


def utility_provenance() -> dict[str, object]:
    """Require the qualification helper itself to be tracked and committed."""

    if _AUTHENTICATED_SELF is None:
        path = Path(__file__).resolve(strict=True)
        repository = REPOSITORY_ROOT
    else:
        path, repository, _ = _AUTHENTICATED_SELF
    if path.relative_to(repository) != UTILITY_RELATIVE_PATH:
        raise ValueError("qualification helper path is inconsistent")
    revision = committed_source_revision(repository, [path])
    profile = retained_file(path, "qualification helper")
    if (
        committed_source_revision(repository, [path]) != revision
        or retained_file(path, "qualification helper") != profile
    ):
        raise ValueError("qualification helper changed during provenance inspection")
    return {
        **profile,
        "revision": revision,
        "committed": True,
    }


def build_provenance(executable: Path, manifest: Path) -> dict[str, object]:
    """Authenticate one executable against its archived build manifest."""

    executable_candidate = executable.expanduser().absolute()
    require_regular_file(executable_candidate, "Athena executable")
    executable = executable_candidate.resolve(strict=True)
    if executable_candidate != executable:
        raise ValueError("Athena executable may not traverse a symlink")
    if not os.access(executable, os.X_OK):
        raise ValueError(f"Athena executable is not executable: {executable}")
    manifest_candidate = manifest.expanduser().absolute()
    require_directory(manifest_candidate, "build manifest")
    manifest = manifest_candidate.resolve(strict=True)
    if manifest_candidate != manifest:
        raise ValueError("build manifest may not traverse a symlink")
    digest_file = manifest / "athena.sha256"
    environment_file = manifest / "environment.txt"
    require_regular_file(digest_file, "build-manifest digest")
    require_regular_file(environment_file, "build-manifest environment")
    try:
        digest_payload, digest_profile = read_regular_profile(
            digest_file, "build-manifest digest", single_link=True
        )
        fields = digest_payload.decode("utf-8").split()
    except UnicodeDecodeError as error:
        raise ValueError("build-manifest digest is not UTF-8") from error
    if not fields:
        raise ValueError("build manifest does not record an executable digest")
    recorded = require_sha256(fields[0], "build-manifest executable digest")
    executable_profile = regular_file_profile(
        executable, "Athena executable", single_link=True
    )
    if recorded != executable_profile["sha256"]:
        raise ValueError("executable digest does not match build manifest")
    try:
        environment_payload, environment_profile = read_regular_profile(
            environment_file, "build-manifest environment", single_link=True
        )
        environment_text = environment_payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("build-manifest environment is not UTF-8") from error
    match = re.search(r"(?m)^git_revision=([0-9a-f]{40})\s*$", environment_text)
    if match is None:
        raise ValueError("build manifest does not record a full git revision")
    inventory = build_manifest_inventory(manifest)
    return {
        "executable": {
            **executable_profile,
            "revision": match.group(1),
        },
        "build_manifest": {
            "path": str(manifest.resolve(strict=True)),
            "athena_sha256": digest_profile,
            "environment": environment_profile,
            "inventory": inventory,
            "inventory_sha256": retained_inventory_sha256(inventory),
        },
    }


def source_bundle_provenance(bundle: Path, root: Path,
                             revisions: list[str]) -> dict[str, object]:
    """Require a retained project bundle containing every launch revision."""

    candidate = bundle.expanduser().absolute()
    require_regular_file(candidate, "source bundle")
    bundle = candidate.resolve(strict=True)
    if candidate != bundle:
        raise ValueError("source bundle may not traverse a symlink")
    source_archives = root / "source-archives"
    require_directory(source_archives, "source-archive directory")
    require_beneath(bundle, source_archives.resolve(strict=True), "source bundle")
    unique_revisions = sorted(set(
        require_git_revision(revision, "bundle revision") for revision in revisions
    ))
    with open_regular_fd(bundle, "source bundle", single_link=True) as descriptor:
        profile = os.fstat(descriptor)
        digest = hashlib.sha256()
        with os.fdopen(os.dup(descriptor), "rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        os.lseek(descriptor, 0, os.SEEK_SET)
        with tempfile.TemporaryDirectory(
            prefix="cgl_lf_qualification_bundle_"
        ) as directory:
            repository = Path(directory) / "source.git"
            descriptor_path = f"/proc/self/fd/{descriptor}"
            cloned = git_run(
                ["clone", "--bare", "--quiet", descriptor_path, str(repository)],
                capture_output=True,
                text=True,
                pass_fds=(descriptor,),
            )
            if cloned.returncode != 0:
                raise ValueError(f"source bundle cannot be cloned: {bundle}")
            for revision in unique_revisions:
                present = git_run(
                    [
                        "-C", str(repository), "cat-file", "-e",
                        f"{revision}^{{commit}}",
                    ],
                    capture_output=True,
                    text=True,
                )
                if present.returncode != 0:
                    raise ValueError(
                        f"source bundle does not contain revision {revision}: {bundle}"
                    )
    bundle_profile = {
        "path": str(bundle),
        "sha256": digest.hexdigest(),
        "size_bytes": profile.st_size,
    }
    return {
        **bundle_profile,
        "verified_revisions": unique_revisions,
    }


def parse_athinput_text(text: str, label: str) -> dict[str, str]:
    """Parse one unambiguous Athena parameter dump into exact block/key values."""

    block = ""
    values: dict[str, str] = {}
    for original in text.splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            if not block:
                raise ValueError(f"{label} has an empty parameter block")
            continue
        key, separator, value = line.partition("=")
        if not block or not separator or not key.strip() or not value.strip():
            raise ValueError(f"{label} has an invalid parameter line: {original!r}")
        qualified = f"{block}/{key.strip()}"
        if qualified in values:
            raise ValueError(f"{label} duplicates parameter {qualified}")
        values[qualified] = value.strip()
    return values


def parse_athinput(path: Path, label: str) -> dict[str, str]:
    """Parse one descriptor-authenticated Athena input deck."""

    try:
        text = read_regular_bytes(path, label, single_link=True).decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} is not UTF-8") from error
    return parse_athinput_text(text, label)


def expected_case_configuration(case_id: str, *, basename: str,
                                tlim: str) -> dict[str, str]:
    """Return the exact reviewed physical and output configuration."""

    policy = CASE_POLICIES[case_id]
    return {
        **COMMON_INPUT_CONTRACT,
        **policy["physics_contract"],
        "job/basename": basename,
        "time/tlim": tlim,
        "mesh/nx1": str(policy["mesh_shape"][0]),
        "mesh/nx2": str(policy["mesh_shape"][1]),
        "mesh/nx3": str(policy["mesh_shape"][2]),
        "meshblock/nx1": str(policy["meshblock_shape"][0]),
        "meshblock/nx2": str(policy["meshblock_shape"][1]),
        "meshblock/nx3": str(policy["meshblock_shape"][2]),
    }


def validate_case_input_contract(case_id: str, path: Path) -> dict[str, str]:
    """Require the exact reviewed frozen deck and its physical configuration."""

    policy = CASE_POLICIES[case_id]
    observed_sha256 = sha256(path)
    if observed_sha256 != policy["input_sha256"]:
        raise ValueError(f"{case_id} input differs from the reviewed frozen deck")
    values = parse_athinput(path, f"{case_id} input")
    expected = expected_case_configuration(
        case_id, basename=str(policy["input_basename"]),
        tlim=COMMON_INPUT_CONTRACT["time/tlim"],
    )
    mismatches = {
        key: {"expected": expected_value, "observed": values.get(key)}
        for key, expected_value in expected.items()
        if values.get(key) != expected_value
    }
    if mismatches:
        raise ValueError(
            f"{case_id} input violates the reviewed physics contract: {mismatches}"
        )
    return values


def load_case_records(matrix_path: Path, source_dir: Path,
                      selected_case_ids: list[str]) -> tuple[dict[str, object],
                                                            list[dict[str, object]]]:
    """Load unique mapped cases and authenticate selected input decks."""

    if sha256(matrix_path) != FROZEN_MATRIX_SHA256:
        raise ValueError("Stage I matrix differs from the reviewed frozen matrix")
    matrix = load_json(matrix_path, "Stage I matrix")
    cases = matrix.get("cases")
    if not isinstance(cases, list):
        raise ValueError("Stage I matrix cases must be a list")
    ids = [str(case.get("id")) for case in cases if isinstance(case, dict)]
    if ids != list(EXPECTED_CASE_IDS):
        raise ValueError("Stage I matrix case identifiers are incomplete or reordered")
    if len(selected_case_ids) != len(set(selected_case_ids)):
        raise ValueError("qualification case IDs must be unique")
    if not selected_case_ids:
        raise ValueError("at least one qualification case ID is required")
    records = []
    for case_id in sorted(selected_case_ids):
        if case_id not in CASE_POLICIES:
            raise ValueError(
                "qualification supports only reviewed cases R04, R12, R16, and R17"
            )
        matches = [
            case for case in cases
            if isinstance(case, dict) and case.get("id") == case_id
        ]
        if len(matches) != 1:
            raise ValueError(f"matrix has missing or ambiguous case {case_id}")
        case = matches[0]
        name = case.get("name")
        input_value = case.get("input")
        resolution = case.get("resolution")
        policy = CASE_POLICIES[case_id]
        if (
            name != policy["case_name"]
            or input_value != policy["input_relative_path"]
            or resolution != policy["resolution"]
        ):
            raise ValueError(f"matrix case {case_id} differs from reviewed identity")
        input_path = (source_dir / input_value).resolve(strict=True)
        require_beneath(input_path, source_dir, f"{case_id} input")
        validate_case_input_contract(case_id, input_path)
        records.append({
            "case_id": case_id,
            "case_name": name,
            "input": retained_file(input_path, f"{case_id} input"),
        })
    return matrix, records


def quote(value: object) -> str:
    """Quote one shell literal."""

    return shlex.quote(str(value))


def format_target(value: float) -> str:
    """Return one stable target-time path token."""

    return f"{value:.12g}".replace(".", "p")


def pack_waves(pilots: list[dict[str, object]],
               max_concurrent_nodes: int) -> list[list[dict[str, object]]]:
    """First-fit-decreasing pack pilots into deterministic bounded waves."""

    if (
        isinstance(max_concurrent_nodes, bool)
        or max_concurrent_nodes < 1
        or max_concurrent_nodes > MAX_CONCURRENT_NODES
    ):
        raise ValueError(
            f"max concurrent nodes must be between 1 and {MAX_CONCURRENT_NODES}"
        )
    waves: list[list[dict[str, object]]] = []
    for pilot in sorted(
        pilots, key=lambda item: (-int(item["nodes"]), str(item["case_id"]))
    ):
        nodes = int(pilot["nodes"])
        if nodes > max_concurrent_nodes:
            raise ValueError(
                f"{pilot['case_id']} {nodes}-node profile exceeds wave ceiling "
                f"{max_concurrent_nodes}"
            )
        for wave in waves:
            if sum(int(item["nodes"]) for item in wave) + nodes <= max_concurrent_nodes:
                wave.append(pilot)
                break
        else:
            waves.append([pilot])
    return [
        sorted(wave, key=lambda item: (str(item["case_id"]), int(item["nodes"])))
        for wave in waves
    ]


def finalize_batch_script(script: str) -> tuple[str, str, str]:
    """Bind one retained script to its normalized and exact digests."""

    placeholder = "BATCH_SCRIPT_SHA256=\n"
    normalized_sha256 = hashlib.sha256(script.encode("utf-8")).hexdigest()
    if script.count(placeholder) != 1:
        raise ValueError("qualification batch script digest placeholder is ambiguous")
    finalized = script.replace(
        placeholder, f"BATCH_SCRIPT_SHA256={normalized_sha256}\n"
    )
    return finalized, normalized_sha256, hashlib.sha256(
        finalized.encode("utf-8")
    ).hexdigest()


def generated_job_script(intent: dict[str, object],
                         provenance: dict[str, object]) -> str:
    """Generate one fail-closed retained qualification batch script."""

    allocation = intent["allocation"]
    paths = intent["paths"]
    input_record = intent["input"]
    executable = provenance["executable"]
    bundle = provenance["source_bundle"]
    matrix = provenance["matrix"]
    helper = provenance["qualification_helper"]
    build_manifest = provenance["build_manifest"]
    build_manifest_requirements = "\n".join(
        "require_sha "
        f"{quote(record['sha256'])} "
        f"{quote(Path(str(build_manifest['path'])) / str(record['name']))} "
        f"{quote('build_manifest_' + str(record['name']))}"
        for record in build_manifest["inventory"]
    )
    build_manifest_names_sha256 = hashlib.sha256(
        "".join(
            f"{record['name']}\n" for record in build_manifest["inventory"]
        ).encode("utf-8")
    ).hexdigest()
    target = intent["target_time"]
    return f"""#!/bin/bash
#SBATCH -J {intent["job_name"]}
#SBATCH -A {ACCOUNT}
#SBATCH -o {paths["slurm_log"]}
#SBATCH -p {PARTITION}
#SBATCH -t {allocation["walltime"]}
#SBATCH -N {allocation["nodes"]}
#SBATCH --gpus-per-node=8
#SBATCH --threads-per-core=1

set -eCuo pipefail
umask 027
BATCH_SCRIPT_SHA256=
QUALIFICATION_ROOT={quote(intent["qualification_root"])}
EXECUTION_CONTRACT={quote(paths["execution_contract"])}
PREPARED_MANIFEST={quote(paths["prepared_manifest"])}
ARCHIVED_INPUT={quote(paths["archived_input"])}
OUT_DIR={quote(paths["output_dir"])}
ENV_LOG={quote(paths["environment_log"])}
SMOKE_DIR={quote(paths["restart_smoke_dir"])}
SMOKE_LOG={quote(paths["restart_smoke_log"])}
SMOKE_RESULT={quote(paths["restart_smoke_result"])}
OUTPUT_BINDING={quote(paths["output_binding"])}

require_sha() {{
  local expected="$1"
  local path="$2"
  local label="$3"
  test -f "${{path}}" || {{ echo "missing ${{label}}: ${{path}}" >&2; exit 1; }}
  test ! -L "${{path}}" || {{ echo "symlink ${{label}}: ${{path}}" >&2; exit 1; }}
  test "$(realpath -e "${{path}}")" = "${{path}}" || {{
    echo "non-canonical ${{label}} path: ${{path}}" >&2
    exit 1
  }}
  test "$(stat -c %h "${{path}}")" -eq 1 || {{
    echo "multiply linked ${{label}}: ${{path}}" >&2
    exit 1
  }}
  local actual
  actual="$(sha256sum "${{path}}" | awk '{{print $1}}')"
  test "${{actual}}" = "${{expected}}" || {{
    echo "checksum mismatch for ${{label}}: ${{path}}" >&2
    exit 1
  }}
}}
require_no_symlink_path() {{
  local root="$1"
  local path="$2"
  local label="$3"
  local current="${{path}}"
  test "$(realpath -e "${{root}}")" = "${{root}}" || {{
    echo "non-canonical ${{label}} root: ${{root}}" >&2
    exit 1
  }}
  while [[ "${{current}}" != "${{root}}" ]]; do
    test ! -L "${{current}}" || {{
      echo "symlink in ${{label}} path: ${{current}}" >&2
      exit 1
    }}
    current="$(dirname "${{current}}")"
  done
  test ! -L "${{root}}" || {{ echo "symlink root: ${{root}}" >&2; exit 1; }}
}}
normalized_script_sha256() {{
  sed -E 's/^BATCH_SCRIPT_SHA256=[0-9a-f]{{64}}$/BATCH_SCRIPT_SHA256=/' "$0" |
    sha256sum | awk '{{print $1}}'
}}
test "${{SLURM_NNODES:?Missing SLURM_NNODES}}" -eq {allocation["nodes"]}
test "${{SLURM_JOB_NAME:?Missing SLURM_JOB_NAME}}" = {quote(intent["job_name"])}
test "$(normalized_script_sha256)" = "${{BATCH_SCRIPT_SHA256}}"
require_no_symlink_path "${{QUALIFICATION_ROOT}}" {quote(paths["run_dir"])} run_dir
require_no_symlink_path "${{QUALIFICATION_ROOT}}" "${{OUT_DIR}}" output_dir
require_no_symlink_path "${{QUALIFICATION_ROOT}}" "${{PREPARED_MANIFEST}}" prepared_manifest
require_no_symlink_path "${{QUALIFICATION_ROOT}}" "${{EXECUTION_CONTRACT}}" execution_contract
require_no_symlink_path "${{QUALIFICATION_ROOT}}" "${{ARCHIVED_INPUT}}" archived_input
require_sha {quote(intent["execution_contract_sha256"])} "${{EXECUTION_CONTRACT}}" execution_contract
require_sha {quote(input_record["sha256"])} "${{ARCHIVED_INPUT}}" archived_input
require_sha {quote(executable["sha256"])} {quote(executable["path"])} executable
require_sha {quote(matrix["sha256"])} {quote(matrix["path"])} matrix
require_sha {quote(bundle["sha256"])} {quote(bundle["path"])} source_bundle
require_sha {quote(helper["sha256"])} {quote(helper["path"])} qualification_helper
test "$(find {quote(build_manifest["path"])} -mindepth 1 -maxdepth 1 -printf '%f\\n' | LC_ALL=C sort | sha256sum | awk '{{print $1}}')" = {quote(build_manifest_names_sha256)}
{build_manifest_requirements}
test -d {quote(paths["run_dir"])}
test ! -L {quote(paths["run_dir"])}
test -d "${{OUT_DIR}}"
test ! -L "${{OUT_DIR}}"
test -d {quote(Path(str(paths["environment_log"])).parent)}
test ! -L {quote(Path(str(paths["environment_log"])).parent)}
test ! -e {quote(paths["environment_log"])}
test ! -e "${{SMOKE_DIR}}"
test ! -e "${{SMOKE_LOG}}"
test ! -e "${{SMOKE_RESULT}}"
test -z "$(find "${{OUT_DIR}}" -mindepth 1 -maxdepth 1 -print -quit)"
module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime
export LD_LIBRARY_PATH=${{CRAY_LD_LIBRARY_PATH}}:${{LD_LIBRARY_PATH:-}}
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_OFI_NUM_CQ_ENTRIES=131072
export FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software
export HSA_XNACK=1
export OMP_NUM_THREADS=1
{{
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
  echo "slurm_job_id=${{SLURM_JOB_ID:?Missing SLURM_JOB_ID}}"
  echo "qualification_execution_contract_sha256={intent["execution_contract_sha256"]}"
  echo "prepared_manifest=${{PREPARED_MANIFEST}}"
  echo "nodes=${{SLURM_NNODES}}"
  echo "ranks=$((SLURM_NNODES * {allocation["ranks_per_node"]}))"
  module -t list 2>&1
}} > {quote(paths["environment_log"])}
srun -N "${{SLURM_NNODES}}" -n "$((SLURM_NNODES * {allocation["ranks_per_node"]}))" \
  --ntasks-per-node={allocation["ranks_per_node"]} -c {allocation["cpus_per_task"]} \
  --threads-per-core=1 --cpu-bind=threads --gpus-per-task=1 --gpu-bind=closest \
  {quote(executable["path"])} -i "${{ARCHIVED_INPUT}}" \
  -d "${{OUT_DIR}}" -t {quote(allocation["athena_walltime"])} \
  job/basename={quote(intent["run_basename"])} time/tlim={quote(target)}
TERMINAL_RESTART="$(
  find "${{OUT_DIR}}/rst/rank_00000000" -mindepth 1 -maxdepth 1 -type f -print |
    LC_ALL=C sort | tail -n 1
)"
test -n "${{TERMINAL_RESTART}}"
test -f "${{TERMINAL_RESTART}}"
test ! -L "${{TERMINAL_RESTART}}"
test "$(realpath -e "${{TERMINAL_RESTART}}")" = "${{TERMINAL_RESTART}}"
test "$(stat -c %h "${{TERMINAL_RESTART}}")" -eq 1
require_no_symlink_path "${{QUALIFICATION_ROOT}}" "${{TERMINAL_RESTART}}" terminal_restart
RESTART_SHA256="$(sha256sum "${{TERMINAL_RESTART}}" | awk '{{print $1}}')"
RESTART_SIZE_BYTES="$(stat -c %s "${{TERMINAL_RESTART}}")"
mkdir "${{SMOKE_DIR}}"
srun -N "${{SLURM_NNODES}}" -n "$((SLURM_NNODES * {allocation["ranks_per_node"]}))" \
  --ntasks-per-node={allocation["ranks_per_node"]} -c {allocation["cpus_per_task"]} \
  --threads-per-core=1 --cpu-bind=threads --gpus-per-task=1 --gpu-bind=closest \
  {quote(executable["path"])} -r "${{TERMINAL_RESTART}}" \
  -d "${{SMOKE_DIR}}" -t {RESTART_SMOKE_WALLTIME} \
  job/basename={quote(intent["run_basename"] + "_restart_smoke")} \
  time/tlim={quote(target)} time/nlim=0 > "${{SMOKE_LOG}}" 2>&1
COMPLETED_UTC="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
printf '{{"schema_version":1,"record_type":"cgl_lf_stage_i_restart_load_smoke",'\
'"execution_contract_sha256":"{intent["execution_contract_sha256"]}",'\
'"executable_sha256":"{executable["sha256"]}",'\
'"job_id":"%s","restart_path":"%s","restart_sha256":"%s",'\
'"restart_size_bytes":%s,"exit_code":0,"completed_utc":"%s"}}\n' \
  "${{SLURM_JOB_ID}}" "${{TERMINAL_RESTART}}" "${{RESTART_SHA256}}" \
  "${{RESTART_SIZE_BYTES}}" "${{COMPLETED_UTC}}" > "${{SMOKE_RESULT}}"
FINISHED_UTC="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
printf 'finished_utc=%s\n' "${{FINISHED_UTC}}" >> {quote(paths["environment_log"])}
python3 - "${{OUT_DIR}}" "${{OUTPUT_BINDING}}" "${{SLURM_JOB_ID}}" \
  {quote(intent["execution_contract_sha256"])} {quote(executable["sha256"])} \
  "${{FINISHED_UTC}}" <<'PY'
import hashlib
import json
import os
from pathlib import Path
import stat
import sys

output = Path(sys.argv[1]).resolve(strict=True)
destination = Path(sys.argv[2])
job_id, contract, executable, completed = sys.argv[3:]
files = []
for current, directories, names in os.walk(output, followlinks=False):
    directories.sort()
    names.sort()
    current_path = Path(current)
    for name in directories:
        profile = (current_path / name).lstat()
        if stat.S_ISLNK(profile.st_mode) or not stat.S_ISDIR(profile.st_mode):
            raise SystemExit(f"invalid output directory: {{current_path / name}}")
    for name in names:
        path = current_path / name
        profile = path.lstat()
        if (
            stat.S_ISLNK(profile.st_mode)
            or not stat.S_ISREG(profile.st_mode)
            or profile.st_nlink != 1
        ):
            raise SystemExit(f"invalid output file: {{path}}")
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        files.append({{
            "path": path.relative_to(output).as_posix(),
            "sha256": digest.hexdigest(),
            "size_bytes": profile.st_size,
        }})
tree_sha256 = hashlib.sha256(
    json.dumps(files, sort_keys=True, separators=(",", ":")).encode("utf-8")
).hexdigest()
record = {{
    "schema_version": 1,
    "record_type": "cgl_lf_stage_i_qualification_output_binding",
    "job_id": job_id,
    "execution_contract_sha256": contract,
    "executable_sha256": executable,
    "completed_utc": completed,
    "output_dir": str(output),
    "files": files,
    "output_tree_sha256": tree_sha256,
}}
payload = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode("utf-8")
flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
if hasattr(os, "O_NOFOLLOW"):
    flags |= os.O_NOFOLLOW
descriptor = os.open(destination, flags, 0o440)
try:
    with os.fdopen(descriptor, "wb", closefd=False) as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
finally:
    os.close(descriptor)
PY
"""


def make_pilot_packet(root: Path, qualification_root: Path,
                      case: dict[str, object], nodes: int,
                      provenance: dict[str, object],
                      provenance_sha256: str) -> dict[str, object]:
    """Build one deterministic, non-submitting qualification packet."""

    case_id = str(case["case_id"])
    policy = CASE_POLICIES[case_id]
    target = float(policy["target_time"])
    packet_id = f"{case_id}_n{nodes:02d}_t{format_target(target)}"
    run_dir = qualification_root / case_id / f"nodes-{nodes:02d}"
    for directory, label in (
        (qualification_root, "qualification root"),
        (run_dir.parent, f"{case_id} qualification directory"),
        (qualification_root / "logs", "qualification log directory"),
    ):
        if directory.exists() or directory.is_symlink():
            require_directory(directory, label)
    if run_dir.exists() or run_dir.is_symlink():
        raise ValueError(f"qualification pilot output already exists: {run_dir}")
    manifest_dir = run_dir / "manifest"
    paths = {
        "run_dir": str(run_dir),
        "output_dir": str(run_dir / "output"),
        "manifest_dir": str(manifest_dir),
        "environment_log": str(manifest_dir / "run_environment.txt"),
        "execution_contract": str(manifest_dir / "execution_contract.json"),
        "prepared_manifest": str(manifest_dir / "prepared_pilot.json"),
        "archived_input": str(manifest_dir / "submitted_input.athinput"),
        "batch_script": str(manifest_dir / "cgl_lf_stage_i_qualification.sbatch"),
        "scheduler_batch_script": str(
            manifest_dir / "scheduler_stored_batch_script.sbatch"
        ),
        "submission_receipt": str(manifest_dir / "submission_receipt.json"),
        "submission_recovery": str(manifest_dir / "submission_recovery.json"),
        "job_binding": str(manifest_dir / "job_binding.json"),
        "scheduler_raw": str(manifest_dir / "scheduler.sacct.txt"),
        "scheduler_evidence": str(manifest_dir / "scheduler_evidence.json"),
        "scientific_evidence": str(manifest_dir / "scientific_evidence.json"),
        "output_binding": str(manifest_dir / "output_binding.json"),
        "restart_smoke_dir": str(run_dir / "restart-smoke"),
        "restart_smoke_log": str(manifest_dir / "restart_load_smoke.log"),
        "restart_smoke_result": str(manifest_dir / "restart_load_smoke.json"),
        "slurm_log": str(qualification_root / "logs" / f"{packet_id}.%j.log"),
    }
    intent_core: dict[str, object] = {
        "packet_id": packet_id,
        "project_root": str(root),
        "qualification_root": str(qualification_root),
        "case_id": case_id,
        "case_name": case["case_name"],
        "profile_class": policy["profile_class"],
        "scientific_policy": policy["scientific_policy"],
        "target_time": target,
        "run_basename": f"qualification_{case_id}_n{nodes:02d}",
        "job_name": f"cglq_{case_id}_n{nodes:02d}",
        "input": case["input"],
        "provenance_sha256": provenance_sha256,
        "allocation": {
            "nodes": nodes,
            "walltime": policy["walltime"],
            "walltime_seconds": parse_walltime(str(policy["walltime"])),
            "athena_walltime": policy["athena_walltime"],
            "athena_walltime_seconds": parse_walltime(
                str(policy["athena_walltime"])
            ),
            "ranks_per_node": RANKS_PER_NODE,
            "cpus_per_task": CPUS_PER_TASK,
        },
        "paths": paths,
    }
    contract_sha256 = stable_json_sha256(intent_core)
    intent = {
        **intent_core,
        "execution_contract_sha256": contract_sha256,
    }
    script = generated_job_script(intent, provenance)
    script, normalized_script_sha256, exact_script_sha256 = finalize_batch_script(
        script
    )
    commands = {
        "batch_script_text": script,
        "batch_script_sha256": exact_script_sha256,
        "normalized_batch_script_sha256": normalized_script_sha256,
        "base_submission_argv": [
            str(SBATCH), "--parsable",
            f"--comment={contract_sha256}",
        ],
    }
    intent["authenticated_commands"] = commands
    intent_sha256 = stable_json_sha256(intent)
    intent["execution_intent_sha256"] = intent_sha256
    return {
        "execution_intent": intent,
        "execution_intent_sha256": intent_sha256,
    }


def build_prepared_wave(root: Path, qualification_root: Path,
                        cases: list[dict[str, object]],
                        provenance: dict[str, object],
                        max_concurrent_nodes: int) -> dict[str, object]:
    """Build a complete deterministic qualification-wave packet."""

    qualification_root = require_qualification_root(root, qualification_root)
    provenance_sha256 = stable_json_sha256(provenance)
    selected_case_ids = sorted(str(case.get("case_id", "")) for case in cases)
    operational_only = any(
        bool(CASE_POLICIES.get(case_id, {}).get("operational_only"))
        for case_id in selected_case_ids
    )
    if operational_only and (
        selected_case_ids != [R17_CASE_ID]
        or max_concurrent_nodes != R17_OPERATIONAL_NODES
    ):
        raise ValueError(
            "R17 operational qualification requires one exclusive eight-node wave"
        )
    pilots = []
    for case in sorted(cases, key=lambda item: str(item["case_id"])):
        case_id = str(case.get("case_id", ""))
        if case_id not in CASE_POLICIES:
            raise ValueError(
                "qualification supports only reviewed cases R04, R12, R16, and R17"
            )
        if case_id == R12_CASE_ID:
            require_r12_fresh_profile_policy()
        for nodes in CASE_POLICIES[case_id]["node_profiles"]:
            pilots.append(
                make_pilot_packet(
                    root, qualification_root, case, nodes, provenance,
                    provenance_sha256,
                )
            )
    packed = pack_waves(
        [
            {
                "case_id": packet["execution_intent"]["case_id"],
                "nodes": packet["execution_intent"]["allocation"]["nodes"],
                "packet": packet,
            }
            for packet in pilots
        ],
        max_concurrent_nodes,
    )
    waves = []
    for number, wave in enumerate(packed, start=1):
        packets = [item["packet"] for item in wave]
        waves.append({
            "wave": number,
            "total_nodes": sum(
                int(packet["execution_intent"]["allocation"]["nodes"])
                for packet in packets
            ),
            "packets": packets,
        })
    return {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_qualification_wave",
        "project_root": str(root),
        "qualification_root": str(qualification_root),
        "policy": {
            "canonical_acceptance_eligible": False,
            "submission_action_available": True,
            "submission_performed": False,
            "operational_only": operational_only,
            "max_concurrent_nodes": max_concurrent_nodes,
            "selection_minimum_runtime_savings_seconds": (
                MINIMUM_RUNTIME_SAVINGS_SECONDS
            ),
            "selection_maximum_node_hour_ratio": MAXIMUM_NODE_HOUR_RATIO,
            "selection_runtime_tie_fraction": RUNTIME_TIE_FRACTION,
        },
        "provenance": provenance,
        "provenance_sha256": provenance_sha256,
        "waves": waves,
    }


def prepare_wave(args: argparse.Namespace) -> dict[str, object]:
    """Collect live provenance and emit one isolated qualification wave."""

    root = require_project_root(args.root, args.allow_local_root)
    qualification_root = require_qualification_root(root, args.qualification_root)
    source_candidate = Path(args.source_dir).expanduser().absolute()
    matrix_candidate = Path(args.matrix).expanduser().absolute()
    executable_candidate = Path(args.executable).expanduser().absolute()
    build_manifest_candidate = Path(args.build_manifest).expanduser().absolute()
    source_dir = source_candidate.resolve(strict=True)
    matrix_path = matrix_candidate.resolve(strict=True)
    executable = executable_candidate.resolve(strict=True)
    build_manifest_path = build_manifest_candidate.resolve(strict=True)
    if any(
        candidate != resolved
        for candidate, resolved in (
            (source_candidate, source_dir),
            (matrix_candidate, matrix_path),
            (executable_candidate, executable),
            (build_manifest_candidate, build_manifest_path),
        )
    ):
        raise ValueError("qualification preparation paths may not traverse symlinks")
    require_beneath(executable, root, "Athena executable")
    require_beneath(build_manifest_path, root, "build manifest")
    _, cases = load_case_records(matrix_path, source_dir, args.case_id)
    selected_paths = [
        matrix_path,
        *(Path(str(case["input"]["path"])) for case in cases),
    ]
    source_revision = committed_source_revision(source_dir, selected_paths)
    if source_revision != FROZEN_SOURCE_REVISION:
        raise ValueError("qualification source must be the reviewed frozen revision")
    helper = utility_provenance()
    build = build_provenance(executable, build_manifest_path)
    if build["executable"]["revision"] != source_revision:
        raise ValueError("qualification executable must be built from frozen source")
    bundle = source_bundle_provenance(
        Path(args.source_bundle).expanduser(),
        root,
        [
            source_revision,
            str(helper["revision"]),
            str(build["executable"]["revision"]),
        ],
    )
    provenance = {
        "source": {
            "directory": str(source_dir),
            "revision": source_revision,
        },
        "source_bundle": bundle,
        "matrix": retained_file(matrix_path, "Stage I matrix"),
        **build,
        "qualification_helper": helper,
    }
    wave = build_prepared_wave(
        root, qualification_root, cases, provenance, args.max_concurrent_nodes
    )
    materialize_prepared_wave(wave)
    return wave


def intent_core(intent: dict[str, object]) -> dict[str, object]:
    """Return the immutable pre-command execution contract."""

    core = dict(intent)
    core.pop("execution_intent_sha256", None)
    core.pop("authenticated_commands", None)
    core.pop("execution_contract_sha256", None)
    return core


def execution_contract_bytes(intent: dict[str, object]) -> bytes:
    """Return exact bytes whose digest is the execution-contract digest."""

    return json.dumps(
        intent_core(intent), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def materialize_prepared_wave(wave: dict[str, object]) -> Path:
    """Retain exact inputs, contracts, manifests, and batch scripts."""

    root = Path(str(wave["project_root"])).resolve(strict=True)
    qualification_root = require_qualification_root(
        root, str(wave["qualification_root"])
    )
    if qualification_root.exists():
        raise ValueError(f"qualification root already exists: {qualification_root}")
    mkdir_exclusive(qualification_root, 0o750)
    mkdir_exclusive(qualification_root / "logs", 0o750)
    for wave_record in wave["waves"]:
        for packet in wave_record["packets"]:
            intent = packet["execution_intent"]
            paths = intent["paths"]
            case_dir = Path(str(paths["run_dir"])).parent
            if not case_dir.exists():
                mkdir_exclusive(case_dir, 0o750)
            run_dir = Path(str(paths["run_dir"]))
            output_dir = Path(str(paths["output_dir"]))
            manifest_dir = Path(str(paths["manifest_dir"]))
            for directory in (run_dir, output_dir, manifest_dir):
                require_no_symlink_chain(
                    directory, qualification_root, "qualification materialization",
                    require_exists=False,
                )
                mkdir_exclusive(directory, 0o750)
            contract = execution_contract_bytes(intent)
            if hashlib.sha256(contract).hexdigest() != intent["execution_contract_sha256"]:
                raise ValueError("execution contract bytes differ from prepared digest")
            write_exclusive(Path(str(paths["execution_contract"])), contract)
            input_source = validate_live_file_record(
                intent["input"], f"{intent['case_id']} source input"
            )
            copy_exclusive(input_source, Path(str(paths["archived_input"])))
            script = str(intent["authenticated_commands"]["batch_script_text"])
            write_exclusive(
                Path(str(paths["batch_script"])), script.encode("utf-8"), mode=0o750
            )
            write_json_exclusive(Path(str(paths["prepared_manifest"])), packet)
    wave_path = qualification_root / "prepared_wave.json"
    write_json_exclusive(wave_path, wave)
    return wave_path


def validate_materialized_packet(packet: dict[str, object],
                                 qualification_root: Path) -> dict[str, Path]:
    """Revalidate all immutable retained launch artifacts for one pilot."""

    intent = packet["execution_intent"]
    paths = intent["paths"]
    resolved = {
        key: require_owned_regular_file(
            Path(str(paths[key])), qualification_root, f"materialized {key}"
        )
        for key in (
            "execution_contract", "prepared_manifest", "archived_input", "batch_script",
        )
    }
    if read_regular_bytes(
        resolved["execution_contract"], "materialized execution contract",
        single_link=True,
    ) != execution_contract_bytes(intent):
        raise ValueError("materialized execution contract differs from prepared intent")
    retained_packet = load_json(resolved["prepared_manifest"], "prepared pilot manifest")
    if retained_packet != packet:
        raise ValueError("materialized prepared pilot manifest differs from wave")
    if sha256(resolved["archived_input"]) != intent["input"]["sha256"]:
        raise ValueError("materialized pilot input differs from authenticated source input")
    commands = intent["authenticated_commands"]
    if (
        read_regular_bytes(
            resolved["batch_script"], "materialized batch script", single_link=True
        ).decode("utf-8") != commands["batch_script_text"]
        or sha256(resolved["batch_script"]) != commands["batch_script_sha256"]
    ):
        raise ValueError("materialized batch script differs from authenticated command")
    if not os.access(resolved["batch_script"], os.X_OK):
        raise ValueError("materialized qualification batch script is not executable")
    return resolved


def validate_prepared_wave(wave: dict[str, object]) -> dict[tuple[str, int],
                                                            dict[str, object]]:
    """Validate a retained wave and return its unique pilot packets."""

    value = require_exact_keys(
        wave,
        {
            "schema_version", "record_type", "project_root", "qualification_root",
            "policy", "provenance", "provenance_sha256", "waves",
        },
        "prepared qualification wave",
    )
    if (
        value["schema_version"] != 1
        or value["record_type"] != "cgl_lf_stage_i_qualification_wave"
    ):
        raise ValueError("prepared qualification wave has wrong schema")
    root = require_project_root(str(value["project_root"]), allow_local_root=False)
    qualification_root = require_qualification_root(
        root, str(value["qualification_root"])
    )
    policy = require_exact_keys(
        value["policy"],
        {
            "canonical_acceptance_eligible", "submission_action_available",
            "submission_performed", "operational_only", "max_concurrent_nodes",
            "selection_minimum_runtime_savings_seconds",
            "selection_maximum_node_hour_ratio", "selection_runtime_tie_fraction",
        },
        "prepared qualification policy",
    )
    if (
        policy["canonical_acceptance_eligible"] is not False
        or policy["submission_action_available"] is not True
        or policy["submission_performed"] is not False
        or not isinstance(policy["operational_only"], bool)
        or policy["selection_minimum_runtime_savings_seconds"]
        != MINIMUM_RUNTIME_SAVINGS_SECONDS
        or policy["selection_maximum_node_hour_ratio"] != MAXIMUM_NODE_HOUR_RATIO
        or policy["selection_runtime_tie_fraction"] != RUNTIME_TIE_FRACTION
    ):
        raise ValueError("prepared qualification policy is not fail-closed")
    max_nodes = policy["max_concurrent_nodes"]
    if (
        isinstance(max_nodes, bool)
        or not isinstance(max_nodes, int)
        or not 1 <= max_nodes <= MAX_CONCURRENT_NODES
    ):
        raise ValueError("prepared qualification wave has invalid node ceiling")
    require_sha256(value["provenance_sha256"], "prepared provenance digest")
    if stable_json_sha256(value["provenance"]) != value["provenance_sha256"]:
        raise ValueError("prepared qualification provenance digest differs")
    waves = value["waves"]
    if not isinstance(waves, list) or not waves:
        raise ValueError("prepared qualification wave has no launch waves")
    packets: dict[tuple[str, int], dict[str, object]] = {}
    selected_input_paths: list[Path] = []
    for expected_number, wave_record in enumerate(waves, start=1):
        record = require_exact_keys(
            wave_record, {"wave", "total_nodes", "packets"},
            "prepared launch wave",
        )
        if record["wave"] != expected_number or not isinstance(record["packets"], list):
            raise ValueError("prepared launch waves are ambiguous or reordered")
        total_nodes = 0
        for packet_value in record["packets"]:
            packet = require_exact_keys(
                packet_value,
                {"execution_intent", "execution_intent_sha256"},
                "prepared pilot packet",
            )
            intent = require_exact_keys(
                packet["execution_intent"],
                {
                    "packet_id", "project_root", "qualification_root", "case_id",
                    "case_name", "profile_class", "target_time", "run_basename",
                    "job_name", "scientific_policy", "input", "provenance_sha256",
                    "allocation", "paths", "execution_contract_sha256",
                    "authenticated_commands", "execution_intent_sha256",
                },
                "prepared pilot execution intent",
            )
            retained_digest = require_sha256(
                packet["execution_intent_sha256"], "execution intent digest"
            )
            intent_without_digest = dict(intent)
            embedded_digest = intent_without_digest.pop(
                "execution_intent_sha256", None
            )
            if (
                embedded_digest != retained_digest
                or stable_json_sha256(intent_without_digest) != retained_digest
            ):
                raise ValueError("prepared pilot execution intent digest differs")
            case_id = str(intent.get("case_id", ""))
            if case_id not in CASE_POLICIES:
                raise ValueError("prepared pilot has invalid case or allocation")
            if case_id == R12_CASE_ID:
                require_r12_fresh_profile_policy()
            allocation = require_exact_keys(
                intent["allocation"],
                {
                    "nodes", "walltime", "walltime_seconds", "athena_walltime",
                    "athena_walltime_seconds", "ranks_per_node", "cpus_per_task",
                },
                "prepared pilot allocation",
            )
            nodes = allocation.get("nodes")
            if (
                isinstance(nodes, bool)
                or not isinstance(nodes, int)
                or nodes not in CASE_POLICIES[case_id]["node_profiles"]
            ):
                raise ValueError("prepared pilot has unauthorized node profile")
            if (
                intent.get("project_root") != str(root)
                or intent.get("qualification_root") != str(qualification_root)
                or intent.get("provenance_sha256") != value["provenance_sha256"]
                or float(intent.get("target_time", -1.0))
                != CASE_POLICIES[case_id]["target_time"]
            ):
                raise ValueError("prepared pilot intent differs from qualification policy")
            case_policy = CASE_POLICIES[case_id]
            if (
                allocation["walltime"] != case_policy["walltime"]
                or allocation["walltime_seconds"]
                != parse_walltime(str(case_policy["walltime"]))
                or allocation["athena_walltime"] != case_policy["athena_walltime"]
                or allocation["athena_walltime_seconds"]
                != parse_walltime(str(case_policy["athena_walltime"]))
                or allocation["ranks_per_node"] != RANKS_PER_NODE
                or allocation["cpus_per_task"] != CPUS_PER_TASK
                or intent["profile_class"] != case_policy["profile_class"]
                or intent["scientific_policy"] != case_policy["scientific_policy"]
                or intent["run_basename"] != f"qualification_{case_id}_n{nodes:02d}"
                or intent["job_name"] != f"cglq_{case_id}_n{nodes:02d}"
                or intent["packet_id"]
                != (
                    f"{case_id}_n{nodes:02d}_t"
                    f"{format_target(float(case_policy['target_time']))}"
                )
                or not isinstance(intent["case_name"], str)
                or not intent["case_name"]
            ):
                raise ValueError("prepared pilot allocation or identity differs")
            if case_policy["operational_only"] and (
                nodes != R17_OPERATIONAL_NODES
                or allocation["ranks_per_node"] * nodes
                != case_policy["required_ranks"]
            ):
                raise ValueError("R17 operational allocation is not exactly 8 nodes/64 ranks")
            selected_input_paths.append(
                validate_live_file_record(intent["input"], f"prepared {case_id} input")
            )
            paths = require_exact_keys(
                intent["paths"],
                {
                    "run_dir", "output_dir", "manifest_dir", "environment_log",
                    "execution_contract", "prepared_manifest", "archived_input",
                    "batch_script", "scheduler_batch_script", "submission_receipt",
                    "submission_recovery", "job_binding",
                    "scheduler_raw",
                    "scheduler_evidence", "scientific_evidence", "output_binding",
                    "restart_smoke_dir", "restart_smoke_log",
                    "restart_smoke_result", "slurm_log",
                },
                "prepared pilot paths",
            )
            expected_run_dir = qualification_root / case_id / f"nodes-{nodes:02d}"
            expected_manifest_dir = expected_run_dir / "manifest"
            expected_paths = {
                "run_dir": str(expected_run_dir),
                "output_dir": str(expected_run_dir / "output"),
                "manifest_dir": str(expected_manifest_dir),
                "environment_log": str(
                    expected_manifest_dir / "run_environment.txt"
                ),
                "execution_contract": str(expected_manifest_dir / "execution_contract.json"),
                "prepared_manifest": str(expected_manifest_dir / "prepared_pilot.json"),
                "archived_input": str(expected_manifest_dir / "submitted_input.athinput"),
                "batch_script": str(
                    expected_manifest_dir / "cgl_lf_stage_i_qualification.sbatch"
                ),
                "scheduler_batch_script": str(
                    expected_manifest_dir / "scheduler_stored_batch_script.sbatch"
                ),
                "submission_receipt": str(
                    expected_manifest_dir / "submission_receipt.json"
                ),
                "submission_recovery": str(
                    expected_manifest_dir / "submission_recovery.json"
                ),
                "job_binding": str(expected_manifest_dir / "job_binding.json"),
                "scheduler_raw": str(expected_manifest_dir / "scheduler.sacct.txt"),
                "scheduler_evidence": str(
                    expected_manifest_dir / "scheduler_evidence.json"
                ),
                "scientific_evidence": str(
                    expected_manifest_dir / "scientific_evidence.json"
                ),
                "output_binding": str(expected_manifest_dir / "output_binding.json"),
                "restart_smoke_dir": str(expected_run_dir / "restart-smoke"),
                "restart_smoke_log": str(
                    expected_manifest_dir / "restart_load_smoke.log"
                ),
                "restart_smoke_result": str(
                    expected_manifest_dir / "restart_load_smoke.json"
                ),
                "slurm_log": str(
                    qualification_root / "logs"
                    / f"{intent['packet_id']}.%j.log"
                ),
            }
            if paths != expected_paths:
                raise ValueError("prepared pilot paths differ from isolated layout")
            for label, path_value in paths.items():
                require_qualification_path(
                    str(path_value), qualification_root, f"prepared pilot {label}"
                )
            commands = require_exact_keys(
                intent["authenticated_commands"],
                {
                    "batch_script_text", "batch_script_sha256",
                    "normalized_batch_script_sha256", "base_submission_argv",
                },
                "prepared pilot commands",
            )
            intent_core = dict(intent)
            intent_core.pop("execution_intent_sha256")
            intent_core.pop("authenticated_commands")
            contract_sha256 = intent_core.pop("execution_contract_sha256")
            if (
                require_sha256(contract_sha256, "execution contract digest")
                != stable_json_sha256(intent_core)
            ):
                raise ValueError("prepared pilot execution contract digest differs")
            regeneration_intent = {
                **intent_core,
                "execution_contract_sha256": contract_sha256,
            }
            regenerated, normalized_sha256, exact_sha256 = finalize_batch_script(
                generated_job_script(regeneration_intent, value["provenance"])
            )
            expected_commands = {
                "batch_script_text": regenerated,
                "batch_script_sha256": exact_sha256,
                "normalized_batch_script_sha256": normalized_sha256,
                "base_submission_argv": [
                    str(SBATCH), "--parsable", f"--comment={contract_sha256}",
                ],
            }
            if commands != expected_commands:
                raise ValueError("prepared pilot commands differ from exact regeneration")
            key = (case_id, nodes)
            if key in packets:
                raise ValueError(f"prepared qualification profile is duplicated: {key}")
            packets[key] = packet
            total_nodes += nodes
        if record["total_nodes"] != total_nodes or total_nodes > max_nodes:
            raise ValueError("prepared launch wave exceeds its node ceiling")
    selected_case_ids = sorted({case_id for case_id, _ in packets})
    observed_operational_only = selected_case_ids == [R17_CASE_ID]
    if policy["operational_only"] != observed_operational_only:
        raise ValueError("prepared qualification operational policy differs")
    if observed_operational_only and (
        max_nodes != R17_OPERATIONAL_NODES
        or len(waves) != 1
        or len(waves[0]["packets"]) != 1
        or waves[0]["total_nodes"] != R17_OPERATIONAL_NODES
    ):
        raise ValueError("R17 operational qualification is not one exclusive wave")
    for case_id in selected_case_ids:
        expected_keys = {
            (case_id, nodes) for nodes in CASE_POLICIES[case_id]["node_profiles"]
        }
        if {key for key in packets if key[0] == case_id} != expected_keys:
            raise ValueError(f"prepared qualification profiles are incomplete for {case_id}")
    expected_packing = (
        [[{"case_id": R17_CASE_ID, "nodes": R17_OPERATIONAL_NODES}]]
        if observed_operational_only
        else pack_waves(
            [
                {"case_id": case_id, "nodes": nodes}
                for case_id, nodes in sorted(packets)
            ],
            max_nodes,
        )
    )
    observed_packing = [
        [
            (
                str(packet["execution_intent"]["case_id"]),
                int(packet["execution_intent"]["allocation"]["nodes"]),
            )
            for packet in wave_record["packets"]
        ]
        for wave_record in waves
    ]
    expected_packing_keys = [
        [(str(item["case_id"]), int(item["nodes"])) for item in wave_record]
        for wave_record in expected_packing
    ]
    if observed_packing != expected_packing_keys:
        raise ValueError("prepared qualification wave packing is not deterministic")
    provenance = value["provenance"]
    source_directory = Path(str(provenance["source"]["directory"])).resolve(strict=True)
    matrix_path = Path(str(provenance["matrix"]["path"])).resolve(strict=True)
    _, matrix_cases = load_case_records(
        matrix_path, source_directory, selected_case_ids
    )
    mapped_cases = {str(case["case_id"]): case for case in matrix_cases}
    for (case_id, _), packet in packets.items():
        intent = packet["execution_intent"]
        expected_case = mapped_cases[case_id]
        if (
            intent["case_name"] != expected_case["case_name"]
            or intent["input"] != expected_case["input"]
        ):
            raise ValueError("prepared pilot identity differs from committed matrix")
    validate_prepared_provenance(value["provenance"], root, selected_input_paths)
    return packets


def parse_job_id(output: str | bytes) -> str:
    """Parse one top-level numeric Slurm job identifier from exact sbatch output."""

    if isinstance(output, bytes):
        try:
            value = output.decode("utf-8").strip()
        except UnicodeDecodeError as error:
            raise ValueError("sbatch output is not UTF-8") from error
    else:
        value = output.strip()
    job_id, separator, cluster = value.partition(";")
    if (
        JOB_ID_PATTERN.fullmatch(job_id) is None
        or (separator and (not cluster or re.fullmatch(r"[A-Za-z0-9_.-]+", cluster) is None))
    ):
        raise ValueError(f"sbatch did not return one parsable job ID: {output!r}")
    return job_id


def collect_live_batch_script(job_id: str,
                              runner=subprocess.run) -> bytes:
    """Read the exact batch-script bytes retained by the Slurm controller."""

    if JOB_ID_PATTERN.fullmatch(job_id) is None:
        raise ValueError("scheduler batch-script query has invalid job ID")
    argv = [str(SCONTROL), "write", "batch_script", job_id, "-"]
    completed = runner(
        argv, check=True, capture_output=True, env=scheduler_environment()
    )
    stdout = completed.stdout
    if isinstance(stdout, str):
        stdout = stdout.encode("utf-8")
    if not isinstance(stdout, bytes) or not stdout:
        raise ValueError("Slurm returned an empty or invalid stored batch script")
    return stdout


def actual_submission_argv(packet: dict[str, object],
                           prior_wave_job_ids: list[str]) -> list[str]:
    """Return the exact authenticated submission argv for one dependency wave."""

    base = list(
        packet["execution_intent"]["authenticated_commands"]["base_submission_argv"]
    )
    if prior_wave_job_ids:
        dependency = "afterok:" + ":".join(prior_wave_job_ids)
        return [*base, f"--dependency={dependency}"]
    return base


def require_global_qualification_lock_profile(profile: os.stat_result,
                                              path: Path) -> None:
    """Require the owned immutable-profile qualification submission lock."""

    if (
        not stat.S_ISREG(profile.st_mode)
        or profile.st_uid != os.geteuid()
        or profile.st_nlink != 1
        or profile.st_size != 0
        or stat.S_IMODE(profile.st_mode) & 0o022
    ):
        raise ValueError(f"qualification lock profile differs: {path}")


class MutationLockBinding:
    """Continuously authenticate one held mutation-lock pathname and inode."""

    def __init__(self, path: Path, parent_path: Path, parent_descriptor: int,
                 parent_profile: os.stat_result, descriptor: int,
                 opened: os.stat_result, profile_validator):
        self.path = path
        self.parent_path = parent_path
        self.parent_descriptor = parent_descriptor
        self.parent_profile = parent_profile
        self.descriptor = descriptor
        self.opened = opened
        self.profile_validator = profile_validator

    def require_bound(self) -> None:
        """Fail unless both the lock entry and its parent retain their identities."""

        require_directory_descriptor_binding(
            self.parent_path,
            self.parent_descriptor,
            self.parent_profile,
            f"{self.path.name} parent",
        )
        current = os.fstat(self.descriptor)
        if (
            (current.st_dev, current.st_ino)
            != (self.opened.st_dev, self.opened.st_ino)
        ):
            raise ValueError(f"mutation lock descriptor changed while held: {self.path}")
        self.profile_validator(current, self.path)
        named = require_bound_entry_identity(
            self.parent_descriptor,
            self.path.name,
            current,
            "mutation lock",
        )
        self.profile_validator(named, self.path)


@contextmanager
def global_qualification_lock(root: Path):
    """Hold and continuously authenticate the cross-root submission lock."""

    runs = require_no_symlink_chain(root / "runs", root, "project runs directory")
    lock_path = runs / GLOBAL_LOCK_NAME
    flags = os.O_RDWR | os.O_CREAT
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    with open_directory_fd(runs, "project runs directory") as (
        runs_descriptor, runs_profile,
    ):
        descriptor = os.open(
            require_entry_name(lock_path.name, "qualification lock"),
            flags,
            0o640,
            dir_fd=runs_descriptor,
        )
        acquired = False
        binding = None
        try:
            profile = os.fstat(descriptor)
            require_global_qualification_lock_profile(profile, lock_path)
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                acquired = True
            except BlockingIOError as error:
                raise ValueError(
                    "another canonical qualification submission holds the lock"
                ) from error
            binding = MutationLockBinding(
                lock_path,
                runs,
                runs_descriptor,
                runs_profile,
                descriptor,
                profile,
                require_global_qualification_lock_profile,
            )
            binding.require_bound()
            yield binding
        finally:
            try:
                if acquired and binding is not None:
                    binding.require_bound()
            finally:
                try:
                    if acquired:
                        fcntl.flock(descriptor, fcntl.LOCK_UN)
                finally:
                    os.close(descriptor)


def require_stage_i_lock_profile(profile: os.stat_result, path: Path) -> None:
    """Require the existing canonical Stage I mutation-lock profile."""

    if (
        not stat.S_ISREG(profile.st_mode)
        or stat.S_IMODE(profile.st_mode) != 0o644
        or profile.st_uid != os.geteuid()
        or profile.st_nlink != 1
        or profile.st_size != 0
    ):
        raise ValueError(f"Stage I lock profile differs: {path}")


@contextmanager
def stage_i_exclusivity_lock(root: Path):
    """Acquire and continuously authenticate the existing Stage I lock."""

    path = require_no_symlink_chain(
        root / STAGE_I_LOCK_NAME, root, "Stage I lock"
    )
    flags = os.O_RDWR
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    with open_directory_fd(root, "Stage I project root") as (
        root_descriptor, root_profile,
    ):
        descriptor = os.open(
            require_entry_name(path.name, "Stage I lock"),
            flags,
            dir_fd=root_descriptor,
        )
        acquired = False
        binding = None
        try:
            opened = os.fstat(descriptor)
            require_stage_i_lock_profile(opened, path)
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                acquired = True
            except OSError as error:
                if error.errno not in (errno.EACCES, errno.EAGAIN):
                    raise
                raise ValueError(f"another Stage I mutation holds {path}") from error
            binding = MutationLockBinding(
                path,
                root,
                root_descriptor,
                root_profile,
                descriptor,
                opened,
                require_stage_i_lock_profile,
            )
            binding.require_bound()
            yield binding
        finally:
            try:
                if acquired and binding is not None:
                    binding.require_bound()
            finally:
                try:
                    if acquired:
                        fcntl.flock(descriptor, fcntl.LOCK_UN)
                finally:
                    os.close(descriptor)


def require_empty_canonical_stage_i_state(root: Path) -> dict[str, object]:
    """Require no active reservation or pending Stage I metadata transaction."""

    accounting = root / "accounting"
    require_directory(accounting, "Stage I accounting directory")
    reservations_path = require_owned_regular_file(
        accounting / STAGE_I_RESERVATIONS_NAME,
        accounting,
        "Stage I reservations",
    )
    try:
        reservations = json.loads(
            read_regular_bytes(
                reservations_path, "Stage I reservations", single_link=True
            ).decode("utf-8"),
            object_pairs_hook=unique_json_object,
            parse_constant=reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("Stage I reservations are invalid JSON") from error
    if not isinstance(reservations, list):
        raise ValueError("Stage I reservations must contain a list")
    active = []
    for index, reservation in enumerate(reservations):
        if not isinstance(reservation, dict):
            raise ValueError(f"Stage I reservation {index} is invalid")
        state = reservation.get("state")
        if state not in {"prepared", "submitted", "recorded", "cancelled"}:
            raise ValueError(f"Stage I reservation {index} has invalid state")
        if state in {"prepared", "submitted"}:
            active.append(index)
    if active:
        raise ValueError("canonical Stage I retains active reservations")
    transaction_stores = []
    for name in STAGE_I_TRANSACTION_NAMES:
        directory = require_no_symlink_chain(
            accounting / name, accounting, "Stage I transaction store"
        )
        require_directory(directory, "Stage I transaction store")
        entries = sorted(item.name for item in directory.iterdir())
        if entries:
            raise ValueError("canonical Stage I transaction store is not empty")
        transaction_stores.append(str(directory))
    return {
        "reservations": str(reservations_path),
        "reservation_count": len(reservations),
        "active_reservations": 0,
        "transaction_stores": transaction_stores,
        "transactions": 0,
    }


class R17SubmissionBoundary:
    """Bind retained canonical admission state to the live Stage I lock."""

    def __init__(self, lock: MutationLockBinding, stage_i_state: dict[str, object]):
        self.lock = lock
        self.stage_i_state = stage_i_state

    def require_bound(self) -> None:
        self.lock.require_bound()


@contextmanager
def r17_submission_exclusivity(root: Path, required: bool):
    """Hold and authenticate the canonical Stage I boundary for R17."""

    if not required:
        yield None
        return
    with stage_i_exclusivity_lock(root) as lock:
        state = require_empty_canonical_stage_i_state(root)
        lock.require_bound()
        yield R17SubmissionBoundary(lock, state)


SQUEUE_FORMAT = "%i|%F|%K|%j|%T|%D|%u"


def live_queue_jobs(runner=subprocess.run) -> list[dict[str, object]]:
    """Query exact live project-account queue state for global node admission."""

    argv = [
        str(SQUEUE), "-A", ACCOUNT, "--array", "-h",
        "-t", "PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED",
        "-o", SQUEUE_FORMAT,
    ]
    completed = runner(
        argv, check=True, capture_output=True, text=True, env=scheduler_environment()
    )
    records = []
    for line in completed.stdout.splitlines():
        if not line:
            continue
        fields = line.split("|")
        if len(fields) != 7:
            raise ValueError("live squeue output has invalid columns")
        (
            job_id, array_job_id, array_task_id, job_name, state, nodes_text,
            owner,
        ) = fields
        if LIVE_JOB_ID_PATTERN.fullmatch(job_id) is None or not owner:
            raise ValueError("live squeue output has invalid job identity")
        if JOB_ID_PATTERN.fullmatch(array_job_id) is None:
            raise ValueError("live squeue output has invalid array parent identity")
        if array_task_id == "N/A":
            if job_id != array_job_id:
                raise ValueError("live squeue output has ambiguous non-array identity")
        elif (
            re.fullmatch(r"[0-9]+", array_task_id) is None
            or job_id != f"{array_job_id}_{array_task_id}"
        ):
            raise ValueError("live squeue output contains a condensed or ambiguous array")
        try:
            nodes = int(nodes_text)
        except ValueError as error:
            raise ValueError("live squeue output has invalid node count") from error
        if nodes <= 0:
            raise ValueError("live squeue output has nonpositive node count")
        records.append({
            "job_id": job_id,
            "array_job_id": array_job_id,
            "array_task_id": array_task_id,
            "job_name": job_name,
            "state": state,
            "nodes": nodes,
            "owner": owner,
        })
    if len({str(record["job_id"]) for record in records}) != len(records):
        raise ValueError("live squeue output contains duplicate job identities")
    return records


def require_global_node_capacity(wave_record: dict[str, object],
                                 submitted_job_ids: set[str],
                                 queue_runner=subprocess.run, *,
                                 exclusive: bool = False,
                                 stage_i_state: dict[str, object] | None = None,
                                 ) -> dict[str, object]:
    """Require one wave plus every potentially overlapping user job to fit."""

    jobs = live_queue_jobs(queue_runner)
    if exclusive and jobs:
        raise ValueError("R17 operational qualification requires an empty account queue")
    foreign = [job for job in jobs if str(job["job_id"]) not in submitted_job_ids]
    active_qualification = [
        job for job in foreign if str(job["job_name"]).startswith(QUALIFICATION_JOB_PREFIX)
    ]
    if active_qualification:
        raise ValueError("another qualification root has live scheduler jobs")
    external_nodes = sum(int(job["nodes"]) for job in foreign)
    wave_nodes = int(wave_record["total_nodes"])
    if external_nodes + wave_nodes > MAX_CONCURRENT_NODES:
        raise ValueError(
            "live user queue plus qualification wave exceeds the ten-node ceiling"
        )
    result = {
        "queried_utc": utc_now(),
        "external_nodes": external_nodes,
        "wave_nodes": wave_nodes,
        "total_nodes_if_admitted": external_nodes + wave_nodes,
        "live_jobs": jobs,
    }
    if exclusive:
        result["exclusive_operational_admission"] = True
        result["canonical_stage_i_state"] = stage_i_state
    return result


def require_post_submission_node_capacity(
    wave_record: dict[str, object],
    submitted_job_ids: set[str],
    current_wave_job_ids: list[str],
    expected_current_jobs: dict[str, tuple[str, int]],
    queue_runner=subprocess.run,
    *,
    exclusive: bool = False,
) -> dict[str, object]:
    """Requery Slurm after submission and require exact account-wide capacity."""

    jobs = live_queue_jobs(queue_runner)
    by_id = {str(job["job_id"]): job for job in jobs}
    for job_id in current_wave_job_ids:
        expected = expected_current_jobs.get(job_id)
        observed = by_id.get(job_id)
        if expected is None or observed is None:
            raise ValueError("post-submission queue lacks a submitted qualification job")
        expected_name, expected_nodes = expected
        if (
            observed["array_task_id"] != "N/A"
            or observed["job_name"] != expected_name
            or observed["nodes"] != expected_nodes
        ):
            raise ValueError("post-submission queue differs from submitted qualification job")
    foreign = [job for job in jobs if str(job["job_id"]) not in submitted_job_ids]
    if exclusive and foreign:
        raise ValueError(
            "R17 operational qualification lost exclusive account ownership"
        )
    if any(
        str(job["job_name"]).startswith(QUALIFICATION_JOB_PREFIX)
        for job in foreign
    ):
        raise ValueError("another qualification root has live scheduler jobs")
    external_nodes = sum(int(job["nodes"]) for job in foreign)
    wave_nodes = int(wave_record["total_nodes"])
    running_nodes = sum(
        int(job["nodes"]) for job in jobs
        if str(job["state"]) != "PENDING"
    )
    if (
        external_nodes + wave_nodes > MAX_CONCURRENT_NODES
        or running_nodes > MAX_CONCURRENT_NODES
    ):
        raise ValueError(
            "post-submission account queue exceeds the ten-node ceiling"
        )
    result = {
        "queried_utc": utc_now(),
        "external_nodes": external_nodes,
        "wave_nodes": wave_nodes,
        "running_nodes": running_nodes,
        "total_nodes_if_wave_runs": external_nodes + wave_nodes,
        "submitted_job_ids": sorted(submitted_job_ids),
        "current_wave_job_ids": list(current_wave_job_ids),
        "live_jobs": jobs,
    }
    return result


def write_journal_record(directory: Path, sequence: int, kind: str,
                         value: dict[str, object],
                         previous_sha256: str | None,
                         mutation_guard=None) -> tuple[Path, dict[str, object]]:
    """Write one exclusive hash-chained submission journal record."""

    record = {
        "schema_version": 1,
        "record_type": f"cgl_lf_stage_i_qualification_submission_{kind}",
        "sequence": sequence,
        "previous_record_sha256": previous_sha256,
        **value,
    }
    path = directory / f"{sequence:06d}-{kind}.json"
    write_exclusive(
        path,
        stable_json(record).encode("utf-8"),
        mode=0o440,
        mutation_guard=mutation_guard,
    )
    return path, retained_file(path, f"submission {kind} record")


def completed_output(value: object) -> str:
    """Normalize one subprocess output stream without losing retained bytes."""

    if value is None:
        return ""
    if isinstance(value, bytes):
        try:
            return value.decode("utf-8")
        except UnicodeDecodeError as error:
            raise ValueError("scheduler command output is not UTF-8") from error
    return str(value)


def write_submission_receipt(packet: dict[str, object], wave_number: int,
                             boundary: dict[str, object], *,
                             submission_argv: list[str],
                             state: str, job_id: str | None,
                             stdout: str = "", stderr: str = "",
                             error: str | None = None,
                             mutation_guard=None) -> dict[str, object]:
    """Durably retain the scheduler outcome before any post-submit checks."""

    if state not in {
        "ambiguous_sbatch_outcome",
        "scheduler_job_id_observed",
        "scheduler_job_id_recovered_from_error",
    }:
        raise ValueError("submission receipt state is invalid")
    if state == "scheduler_job_id_observed":
        if job_id is None or JOB_ID_PATTERN.fullmatch(job_id) is None or error is not None:
            raise ValueError("known submission receipt has invalid job identity")
    elif state == "scheduler_job_id_recovered_from_error":
        if (
            job_id is None
            or JOB_ID_PATTERN.fullmatch(job_id) is None
            or not isinstance(error, str)
            or not error
        ):
            raise ValueError("recovered submission receipt has invalid job identity")
    elif job_id is not None:
        raise ValueError("ambiguous submission receipt may not assert a job ID")
    intent = packet["execution_intent"]
    record = {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_qualification_submission_receipt",
        "recorded_utc": utc_now(),
        "state": state,
        "wave": wave_number,
        "packet_id": intent["packet_id"],
        "execution_intent_sha256": packet["execution_intent_sha256"],
        "execution_contract_sha256": intent["execution_contract_sha256"],
        "submission_argv": submission_argv,
        "submission_stdin_sha256": intent["authenticated_commands"]["batch_script_sha256"],
        "submission_stdin_size_bytes": len(
            intent["authenticated_commands"]["batch_script_text"].encode("utf-8")
        ),
        "submission_boundary": boundary,
        "sbatch_stdout": stdout,
        "sbatch_stderr": stderr,
        "job_id": job_id,
        "error": error,
    }
    path = Path(str(intent["paths"]["submission_receipt"]))
    write_exclusive(
        path,
        stable_json(record).encode("utf-8"),
        mode=0o440,
        mutation_guard=mutation_guard,
    )
    return record


def validate_submission_receipt(record: object, packet: dict[str, object],
                                qualification_root: Path,
                                *, expected_job_id: str | None,
                                expected_submission_argv: list[str] | None = None,
                                ) -> tuple[dict[str, object], Path]:
    """Authenticate one immutable pre-verification scheduler receipt."""

    intent = packet["execution_intent"]
    value = require_exact_keys(
        record,
        {
            "schema_version", "record_type", "recorded_utc", "state", "wave",
            "packet_id", "execution_intent_sha256", "execution_contract_sha256",
            "submission_argv", "submission_stdin_sha256",
            "submission_stdin_size_bytes", "submission_boundary", "sbatch_stdout",
            "sbatch_stderr", "job_id", "error",
        },
        "submission receipt",
    )
    parse_utc_timestamp(value["recorded_utc"], "submission receipt time")
    path = require_owned_regular_file(
        Path(str(intent["paths"]["submission_receipt"])),
        qualification_root,
        "submission receipt",
    )
    if path.stat().st_mode & 0o222:
        raise ValueError("submission receipt must be immutable")
    boundary_path = validate_retained_file(
        value["submission_boundary"], qualification_root, "submission receipt boundary"
    )
    if (
        boundary_path.parent != qualification_root / "submission-journal"
        or boundary_path.stat().st_mode & 0o222
    ):
        raise ValueError("submission receipt boundary path differs")
    boundary = load_json(boundary_path, "submission receipt boundary")
    if (
        value["schema_version"] != 1
        or value["record_type"] != "cgl_lf_stage_i_qualification_submission_receipt"
        or value["packet_id"] != intent["packet_id"]
        or value["execution_intent_sha256"] != packet["execution_intent_sha256"]
        or value["execution_contract_sha256"] != intent["execution_contract_sha256"]
        or (
            expected_submission_argv is not None
            and value["submission_argv"] != expected_submission_argv
        )
        or value["submission_stdin_sha256"]
        != intent["authenticated_commands"]["batch_script_sha256"]
        or value["submission_stdin_size_bytes"]
        != len(intent["authenticated_commands"]["batch_script_text"].encode("utf-8"))
        or boundary.get("record_type")
        != "cgl_lf_stage_i_qualification_submission_boundary"
        or boundary.get("wave") != value["wave"]
        or boundary.get("execution_intent_sha256") != packet["execution_intent_sha256"]
        or boundary.get("execution_contract_sha256")
        != intent["execution_contract_sha256"]
        or boundary.get("submission_argv") != value["submission_argv"]
        or boundary.get("submission_stdin_sha256") != value["submission_stdin_sha256"]
        or boundary.get("submission_stdin_size_bytes")
        != value["submission_stdin_size_bytes"]
    ):
        raise ValueError("submission receipt differs from authenticated packet")
    job_id = value["job_id"]
    if expected_job_id is None:
        if (
            value["state"] != "ambiguous_sbatch_outcome"
            or job_id is not None
            or not isinstance(value["error"], str)
            or not value["error"]
        ):
            raise ValueError("submission receipt is not an ambiguous outcome")
    elif (
        value["state"] not in {
            "scheduler_job_id_observed",
            "scheduler_job_id_recovered_from_error",
        }
        or job_id != expected_job_id
        or JOB_ID_PATTERN.fullmatch(str(job_id)) is None
        or (
            value["state"] == "scheduler_job_id_observed"
            and value["error"] is not None
        )
        or (
            value["state"] == "scheduler_job_id_recovered_from_error"
            and (
                not isinstance(value["error"], str)
                or not value["error"]
            )
        )
        or parse_job_id(str(value["sbatch_stdout"])) != expected_job_id
    ):
        raise ValueError("submission receipt does not bind the observed scheduler job")
    return value, path


def cancel_known_submission(job_id: str, *,
                            cancel_runner=subprocess.run,
                            queue_runner=subprocess.run) -> dict[str, object]:
    """Request cancellation and prove the known job is absent from the live queue."""

    if JOB_ID_PATTERN.fullmatch(job_id) is None:
        raise ValueError("submission cancellation job ID is invalid")
    command = [str(SCANCEL), job_id]
    returncode = None
    stdout = ""
    stderr = ""
    command_error = None
    try:
        completed = cancel_runner(
            command, check=False, capture_output=True, text=True,
            env=scheduler_environment(),
        )
        returncode = int(completed.returncode)
        stdout = completed_output(completed.stdout)
        stderr = completed_output(completed.stderr)
    except (OSError, subprocess.SubprocessError, ValueError) as error:
        command_error = str(error)
    queue_error = None
    jobs: list[dict[str, object]] = []
    try:
        jobs = live_queue_jobs(queue_runner)
    except (OSError, subprocess.SubprocessError, ValueError) as error:
        queue_error = str(error)
    active = [job for job in jobs if str(job["job_id"]) == job_id]
    confirmed_absent = queue_error is None and not active
    return {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_qualification_cancellation_evidence",
        "checked_utc": utc_now(),
        "job_id": job_id,
        "command": command,
        "returncode": returncode,
        "stdout": stdout,
        "stderr": stderr,
        "command_error": command_error,
        "queue_error": queue_error,
        "active_matching_jobs": active,
        "confirmed_absent": confirmed_absent,
    }


def write_submission_recovery(packet: dict[str, object],
                              qualification_root: Path,
                              receipt: dict[str, object] | None, *,
                              boundary: dict[str, object] | None = None,
                              failure: str,
                              cancellation: dict[str, object] | None,
                              notes: str,
                              resolved_job_id: str | None = None,
                              mutation_guard=None) -> Path:
    """Retain immutable resolution evidence for a failed submission boundary."""

    if not notes.strip():
        raise ValueError("submission recovery notes must be nonempty")
    intent = packet["execution_intent"]
    receipt_value = None
    receipt_record = None
    if receipt is not None:
        receipt_value, receipt_path = validate_submission_receipt(
            receipt,
            packet,
            qualification_root,
            expected_job_id=(
                str(receipt["job_id"]) if receipt["job_id"] is not None else None
            ),
        )
        receipt_record = retained_file(receipt_path, "submission receipt")
        boundary = receipt_value["submission_boundary"]
    if boundary is None:
        raise ValueError("submission recovery lacks a retained submission boundary")
    boundary_path = validate_retained_file(
        boundary, qualification_root, "submission recovery boundary"
    )
    if (
        boundary_path.parent != qualification_root / "submission-journal"
        or boundary_path.stat().st_mode & 0o222
    ):
        raise ValueError("submission recovery boundary path differs")
    boundary_value = load_json(boundary_path, "submission recovery boundary")
    if (
        boundary_value.get("record_type")
        != "cgl_lf_stage_i_qualification_submission_boundary"
        or boundary_value.get("execution_intent_sha256")
        != packet["execution_intent_sha256"]
        or boundary_value.get("execution_contract_sha256")
        != intent["execution_contract_sha256"]
    ):
        raise ValueError("submission recovery boundary differs from packet")
    if resolved_job_id is not None:
        if JOB_ID_PATTERN.fullmatch(resolved_job_id) is None:
            raise ValueError("submission recovery resolved job ID is invalid")
        if (
            receipt_value is not None
            and receipt_value["job_id"] is not None
            and receipt_value["job_id"] != resolved_job_id
        ):
            raise ValueError("submission recovery resolved job ID differs from receipt")
    state = "ambiguous"
    if cancellation is not None:
        state = (
            "cancelled_confirmed"
            if cancellation.get("confirmed_absent") is True
            else "cancellation_unconfirmed"
        )
    record = {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_qualification_submission_recovery",
        "recorded_utc": utc_now(),
        "state": state,
        "packet_id": intent["packet_id"],
        "execution_intent_sha256": packet["execution_intent_sha256"],
        "execution_contract_sha256": intent["execution_contract_sha256"],
        "job_id": resolved_job_id or (
            receipt_value["job_id"] if receipt_value is not None else None
        ),
        "submission_boundary": boundary,
        "submission_receipt": receipt_record,
        "failure": failure,
        "cancellation": cancellation,
        "notes": notes,
    }
    path = Path(str(intent["paths"]["submission_recovery"]))
    write_exclusive(
        path,
        stable_json(record).encode("utf-8"),
        mode=0o440,
        mutation_guard=mutation_guard,
    )
    return path


def recover_ambiguous_submission(prepared_wave_path: Path, packet_id: str,
                                 job_id: str, notes: str, *,
                                 batch_script_runner=subprocess.run,
                                 cancel_runner=subprocess.run,
                                 queue_runner=subprocess.run) -> dict[str, object]:
    """Authenticate and cancel one operator-identified ambiguous sbatch outcome."""

    require_regular_file(prepared_wave_path, "prepared qualification wave")
    wave = load_json(prepared_wave_path, "prepared qualification wave")
    packets = validate_prepared_wave(wave)
    root = require_project_root(str(wave["project_root"]), allow_local_root=False)
    qualification_root = require_qualification_root(
        root, str(wave["qualification_root"])
    )
    expected_wave_path = qualification_root / "prepared_wave.json"
    if prepared_wave_path.resolve(strict=True) != expected_wave_path.resolve(strict=True):
        raise ValueError("submission recovery requires the canonical retained prepared wave")
    if JOB_ID_PATTERN.fullmatch(job_id) is None:
        raise ValueError("submission recovery job ID is invalid")
    matches = [
        packet for packet in packets.values()
        if packet["execution_intent"]["packet_id"] == packet_id
    ]
    if len(matches) != 1:
        raise ValueError("submission recovery packet ID is missing or ambiguous")
    packet = matches[0]
    intent = packet["execution_intent"]
    if Path(str(intent["paths"]["job_binding"])).exists():
        raise ValueError("submission recovery refuses an already bound job")
    recovery_path = Path(str(intent["paths"]["submission_recovery"]))
    if recovery_path.exists() or recovery_path.is_symlink():
        raise ValueError("submission recovery record already exists")
    receipt_path = require_owned_regular_file(
        Path(str(intent["paths"]["submission_receipt"])),
        qualification_root,
        "submission receipt",
    )
    receipt, _ = validate_submission_receipt(
        load_json(receipt_path, "submission receipt"),
        packet,
        qualification_root,
        expected_job_id=None,
    )
    recorded = parse_utc_timestamp(
        receipt["recorded_utc"], "ambiguous submission receipt time"
    )
    current = datetime.now(timezone.utc)
    if (
        recorded > current + timedelta(seconds=FUTURE_TIMESTAMP_TOLERANCE_SECONDS)
        or current - recorded > timedelta(seconds=SUBMISSION_RECOVERY_WINDOW_SECONDS)
    ):
        raise ValueError("ambiguous submission receipt is outside the recovery window")
    expected_script = intent["authenticated_commands"]["batch_script_text"].encode("utf-8")
    if collect_live_batch_script(job_id, batch_script_runner) != expected_script:
        raise ValueError(
            "operator-identified ambiguous job does not store the authenticated script"
        )
    cancellation = cancel_known_submission(
        job_id, cancel_runner=cancel_runner, queue_runner=queue_runner
    )
    path = write_submission_recovery(
        packet,
        qualification_root,
        receipt,
        failure="operator resolved ambiguous sbatch outcome",
        cancellation=cancellation,
        notes=notes,
        resolved_job_id=job_id,
    )
    if cancellation["confirmed_absent"] is not True:
        raise ValueError("ambiguous submission cancellation could not be confirmed")
    return {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_qualification_submission_recovery_result",
        "packet_id": packet_id,
        "job_id": job_id,
        "submission_recovery": retained_file(path, "submission recovery"),
        "cancellation": cancellation,
    }


def submit_all_waves(prepared_wave_path: Path,
                     runner=subprocess.run,
                     queue_runner=subprocess.run,
                     batch_script_runner=subprocess.run,
                     cancel_runner=subprocess.run) -> dict[str, object]:
    """Submit all bounded waves with Slurm afterok barriers and retain bindings."""

    require_regular_file(prepared_wave_path, "prepared qualification wave")
    wave = load_json(prepared_wave_path, "prepared qualification wave")
    packets = validate_prepared_wave(wave)
    root = require_project_root(str(wave["project_root"]), allow_local_root=False)
    qualification_root = require_qualification_root(
        root, str(wave["qualification_root"])
    )
    expected_wave_path = qualification_root / "prepared_wave.json"
    if prepared_wave_path.resolve(strict=True) != expected_wave_path.resolve(strict=True):
        raise ValueError("submission requires the canonical retained prepared wave")
    journal = qualification_root / "submission-journal"
    if journal.exists() or journal.is_symlink():
        raise ValueError("qualification submission journal already exists")
    for packet in packets.values():
        validate_materialized_packet(packet, qualification_root)
        for key in (
            "job_binding", "scheduler_batch_script", "submission_receipt",
            "submission_recovery", "output_binding",
        ):
            candidate = Path(str(packet["execution_intent"]["paths"][key]))
            if candidate.exists() or candidate.is_symlink():
                raise ValueError(
                    f"qualification submission artifact already exists: {candidate}"
                )
    previous_job_ids: list[str] = []
    retained_bindings = []
    retained_journal = []
    sequence = 0
    previous_record_sha256 = None
    submitted_job_ids: set[str] = set()
    r17_operational = set(packets) == {(R17_CASE_ID, R17_OPERATIONAL_NODES)}
    with global_qualification_lock(root) as qualification_lock, r17_submission_exclusivity(
        root, r17_operational
    ) as r17_boundary:
        def require_submission_mutation_boundary() -> None:
            qualification_lock.require_bound()
            if r17_boundary is not None:
                r17_boundary.require_bound()

        stage_i_state = (
            None if r17_boundary is None else r17_boundary.stage_i_state
        )
        mkdir_exclusive(
            journal, 0o750,
            mutation_guard=require_submission_mutation_boundary,
        )
        for wave_record in wave["waves"]:
            admission = require_global_node_capacity(
                wave_record,
                submitted_job_ids,
                queue_runner,
                exclusive=r17_operational,
                stage_i_state=stage_i_state,
            )
            current_job_ids = []
            expected_current_jobs: dict[str, tuple[str, int]] = {}
            for packet in wave_record["packets"]:
                intent = packet["execution_intent"]
                resolved = validate_materialized_packet(packet, qualification_root)
                script_bytes = read_regular_bytes(
                    resolved["batch_script"], "pilot batch script", single_link=True
                )
                if hashlib.sha256(script_bytes).hexdigest() != (
                    intent["authenticated_commands"]["batch_script_sha256"]
                ):
                    raise ValueError("submission stdin differs from authenticated script")
                argv = actual_submission_argv(packet, previous_job_ids)
                sequence += 1
                boundary_path, boundary_record = write_journal_record(
                    journal,
                    sequence,
                    "boundary",
                    {
                        "recorded_utc": utc_now(),
                        "wave": wave_record["wave"],
                        "execution_intent_sha256": packet["execution_intent_sha256"],
                        "execution_contract_sha256": intent["execution_contract_sha256"],
                        "submission_argv": argv,
                        "submission_stdin_sha256": hashlib.sha256(script_bytes).hexdigest(),
                        "submission_stdin_size_bytes": len(script_bytes),
                        "queue_admission": admission,
                    },
                    previous_record_sha256,
                    mutation_guard=require_submission_mutation_boundary,
                )
                retained_journal.append(boundary_record)
                previous_record_sha256 = str(boundary_record["sha256"])
                completed = None
                try:
                    require_submission_mutation_boundary()
                    completed = runner(
                        argv, input=script_bytes, check=True, capture_output=True,
                        env=scheduler_environment(),
                    )
                    stdout = completed_output(completed.stdout)
                    stderr = completed_output(completed.stderr)
                    job_id = parse_job_id(stdout)
                except Exception as error:
                    raw_stdout = getattr(completed, "stdout", None)
                    raw_stderr = getattr(completed, "stderr", None)
                    if raw_stdout is None:
                        raw_stdout = getattr(error, "stdout", None)
                    if raw_stdout is None:
                        raw_stdout = getattr(error, "output", "")
                    if raw_stderr is None:
                        raw_stderr = getattr(error, "stderr", "")
                    stdout = completed_output(raw_stdout)
                    stderr = completed_output(raw_stderr)
                    failure = f"{type(error).__name__}: {error}"
                    try:
                        recovered_job_id = parse_job_id(stdout)
                    except ValueError:
                        recovered_job_id = None
                    if recovered_job_id is not None:
                        try:
                            receipt = write_submission_receipt(
                                packet,
                                int(wave_record["wave"]),
                                boundary_record,
                                submission_argv=argv,
                                state="scheduler_job_id_recovered_from_error",
                                job_id=recovered_job_id,
                                stdout=stdout,
                                stderr=stderr,
                                error=failure,
                                mutation_guard=require_submission_mutation_boundary,
                            )
                        except Exception as receipt_error:
                            cancellation = cancel_known_submission(
                                recovered_job_id,
                                cancel_runner=cancel_runner,
                                queue_runner=queue_runner,
                            )
                            write_submission_recovery(
                                packet,
                                qualification_root,
                                None,
                                boundary=boundary_record,
                                failure=(
                                    f"{failure}; receipt failure: "
                                    f"{type(receipt_error).__name__}: {receipt_error}"
                                ),
                                cancellation=cancellation,
                                notes=(
                                    "automatic fail-closed recovery after sbatch "
                                    "failed but exposed an exact scheduler job ID"
                                ),
                                resolved_job_id=recovered_job_id,
                                mutation_guard=require_submission_mutation_boundary,
                            )
                            raise ValueError(
                                "qualification sbatch failed after exposing an exact "
                                "job ID; cancellation was requested but the receipt "
                                f"could not be retained: {receipt_error}"
                            ) from error
                        cancellation = cancel_known_submission(
                            recovered_job_id,
                            cancel_runner=cancel_runner,
                            queue_runner=queue_runner,
                        )
                        write_submission_recovery(
                            packet,
                            qualification_root,
                            receipt,
                            boundary=boundary_record,
                            failure=failure,
                            cancellation=cancellation,
                            notes=(
                                "automatic fail-closed recovery after sbatch failed "
                                "but exposed an exact scheduler job ID"
                            ),
                            mutation_guard=require_submission_mutation_boundary,
                        )
                        raise ValueError(
                            "qualification sbatch failed after exposing an exact job "
                            "ID; scheduler job cancellation was requested and retained: "
                            f"{error}"
                        ) from error
                    write_submission_receipt(
                        packet,
                        int(wave_record["wave"]),
                        boundary_record,
                        submission_argv=argv,
                        state="ambiguous_sbatch_outcome",
                        job_id=None,
                        stdout=stdout,
                        stderr=stderr,
                        error=failure,
                        mutation_guard=require_submission_mutation_boundary,
                    )
                    raise ValueError(
                        "qualification submission outcome is ambiguous; "
                        "recover the immutable submission receipt before retrying"
                    ) from error
                try:
                    receipt = write_submission_receipt(
                        packet,
                        int(wave_record["wave"]),
                        boundary_record,
                        submission_argv=argv,
                        state="scheduler_job_id_observed",
                        job_id=job_id,
                        stdout=stdout,
                        stderr=stderr,
                        mutation_guard=require_submission_mutation_boundary,
                    )
                except Exception as error:
                    cancellation = cancel_known_submission(
                        job_id,
                        cancel_runner=cancel_runner,
                        queue_runner=queue_runner,
                    )
                    write_submission_recovery(
                        packet,
                        qualification_root,
                        None,
                        boundary=boundary_record,
                        failure=f"{type(error).__name__}: {error}",
                        cancellation=cancellation,
                        notes=(
                            "automatic fail-closed recovery after the durable known-job "
                            "submission receipt could not be retained"
                        ),
                        resolved_job_id=job_id,
                        mutation_guard=require_submission_mutation_boundary,
                    )
                    raise ValueError(
                        "qualification known-job submission receipt failed; "
                        "scheduler job cancellation was requested and retained: "
                        f"{error}"
                    ) from error
                try:
                    scheduler_script = collect_live_batch_script(
                        job_id, batch_script_runner
                    )
                    if scheduler_script != script_bytes:
                        raise ValueError(
                            "Slurm-stored batch script differs from authenticated stdin"
                        )
                    scheduler_script_path = Path(
                        str(intent["paths"]["scheduler_batch_script"])
                    )
                    write_exclusive(
                        scheduler_script_path,
                        scheduler_script,
                        mode=0o440,
                        mutation_guard=require_submission_mutation_boundary,
                    )
                    scheduler_script_record = retained_file(
                        scheduler_script_path, "scheduler-stored batch script"
                    )
                    current_job_ids.append(job_id)
                    submitted_job_ids.add(job_id)
                    expected_current_jobs[job_id] = (
                        str(intent["job_name"]), int(intent["allocation"]["nodes"])
                    )
                    post_submission_queue = require_post_submission_node_capacity(
                        wave_record,
                        submitted_job_ids,
                        current_job_ids,
                        expected_current_jobs,
                        queue_runner,
                        exclusive=r17_operational,
                    )
                    sequence += 1
                    result_path, result_record = write_journal_record(
                        journal,
                        sequence,
                        "result",
                        {
                            "recorded_utc": utc_now(),
                            "wave": wave_record["wave"],
                            "execution_intent_sha256": packet["execution_intent_sha256"],
                            "execution_contract_sha256": intent["execution_contract_sha256"],
                            "submission_argv": argv,
                            "submission_stdin_sha256": hashlib.sha256(script_bytes).hexdigest(),
                            "submission_stdin_size_bytes": len(script_bytes),
                            "sbatch_stdout": stdout,
                            "sbatch_stderr": stderr,
                            "job_id": job_id,
                            "submission_receipt": retained_file(
                                Path(str(intent["paths"]["submission_receipt"])),
                                "submission receipt",
                            ),
                            "scheduler_batch_script": scheduler_script_record,
                            "post_submission_queue": post_submission_queue,
                            "boundary_sha256": boundary_record["sha256"],
                        },
                        previous_record_sha256,
                        mutation_guard=require_submission_mutation_boundary,
                    )
                    retained_journal.append(result_record)
                    previous_record_sha256 = str(result_record["sha256"])
                    binding = {
                        "schema_version": 4,
                        "record_type": "cgl_lf_stage_i_qualification_job_binding",
                        "submitted_utc": utc_now(),
                        "wave": wave_record["wave"],
                        "job_id": job_id,
                        "job_name": intent["job_name"],
                        "nodes": intent["allocation"]["nodes"],
                        "prior_wave_job_ids": list(previous_job_ids),
                        "execution_intent_sha256": packet["execution_intent_sha256"],
                        "execution_contract_sha256": intent["execution_contract_sha256"],
                        "submission_argv": argv,
                        "submission_stdin_sha256": hashlib.sha256(script_bytes).hexdigest(),
                        "submission_stdin_size_bytes": len(script_bytes),
                        "sbatch_stdout": stdout,
                        "sbatch_stderr": stderr,
                        "submission_receipt": retained_file(
                            Path(str(intent["paths"]["submission_receipt"])),
                            "submission receipt",
                        ),
                        "scheduler_batch_script": scheduler_script_record,
                        "post_submission_queue": post_submission_queue,
                        "submission_boundary": retained_file(
                            boundary_path, "submission boundary"
                        ),
                        "submission_result": retained_file(result_path, "submission result"),
                        "prepared_wave": retained_file(expected_wave_path, "prepared wave"),
                        "prepared_manifest": retained_file(
                            Path(str(intent["paths"]["prepared_manifest"])),
                            "prepared pilot manifest",
                        ),
                        "archived_input": retained_file(
                            Path(str(intent["paths"]["archived_input"])),
                            "archived pilot input",
                        ),
                        "batch_script": retained_file(
                            Path(str(intent["paths"]["batch_script"])), "pilot batch script"
                        ),
                    }
                    binding_path = Path(str(intent["paths"]["job_binding"]))
                    write_json_exclusive(
                        binding_path,
                        binding,
                        mutation_guard=require_submission_mutation_boundary,
                    )
                    retained_bindings.append(retained_file(binding_path, "job binding"))
                except Exception as error:
                    cancellation = cancel_known_submission(
                        job_id,
                        cancel_runner=cancel_runner,
                        queue_runner=queue_runner,
                    )
                    write_submission_recovery(
                        packet,
                        qualification_root,
                        receipt,
                        boundary=boundary_record,
                        failure=f"{type(error).__name__}: {error}",
                        cancellation=cancellation,
                        notes=(
                            "automatic fail-closed recovery after a known scheduler "
                            "job failed post-submit authentication"
                        ),
                        mutation_guard=require_submission_mutation_boundary,
                    )
                    raise ValueError(
                        "qualification post-submit authentication failed; "
                        "known scheduler job cancellation was requested and retained: "
                        f"{error}"
                    ) from error
            previous_job_ids = current_job_ids
        chmod_directory(
            journal, 0o550,
            mutation_guard=require_submission_mutation_boundary,
        )
    return {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_qualification_submission",
        "prepared_wave": retained_file(expected_wave_path, "prepared wave"),
        "submission_journal": retained_journal,
        "job_bindings": retained_bindings,
    }


def validate_submission_journal(wave: dict[str, object],
                                qualification_root: Path) -> None:
    """Require one complete exclusive hash-chained submission journal."""

    directory = require_no_symlink_chain(
        qualification_root / "submission-journal",
        qualification_root,
        "submission journal directory",
    )
    require_directory(directory, "submission journal directory")
    if directory.stat().st_mode & 0o222:
        raise ValueError("submission journal directory must be immutable")
    expected_count = 2 * sum(
        len(record["packets"]) for record in wave["waves"]
    )
    files = sorted(directory.iterdir())
    if len(files) != expected_count:
        raise ValueError("submission journal is incomplete or ambiguous")
    previous_sha256 = None
    for sequence, path in enumerate(files, start=1):
        expected_kind = "boundary" if sequence % 2 else "result"
        if path.name != f"{sequence:06d}-{expected_kind}.json":
            raise ValueError("submission journal ordering is ambiguous")
        path = require_owned_regular_file(path, qualification_root, "submission record")
        if path.stat().st_mode & 0o222:
            raise ValueError("submission journal records must be immutable")
        record = load_json(path, "submission record")
        if (
            record.get("schema_version") != 1
            or record.get("sequence") != sequence
            or record.get("previous_record_sha256") != previous_sha256
            or record.get("record_type")
            != f"cgl_lf_stage_i_qualification_submission_{expected_kind}"
        ):
            raise ValueError("submission journal hash chain differs")
        previous_sha256 = sha256(path)


def validate_post_submission_queue(record: object,
                                   current_wave_job_ids: list[str], *,
                                   expected_wave_nodes: int,
                                   expected_current_job: tuple[str, str, int],
                                   ) -> None:
    """Validate one retained post-submit account-wide queue snapshot."""

    value = require_exact_keys(
        record,
        {
            "queried_utc", "external_nodes", "wave_nodes", "running_nodes",
            "total_nodes_if_wave_runs", "submitted_job_ids",
            "current_wave_job_ids", "live_jobs",
        },
        "post-submission queue audit",
    )
    parse_utc_timestamp(value["queried_utc"], "post-submission queue query time")
    if (
        value["current_wave_job_ids"] != current_wave_job_ids
        or len(current_wave_job_ids) != len(set(current_wave_job_ids))
        or any(JOB_ID_PATTERN.fullmatch(job_id) is None for job_id in current_wave_job_ids)
    ):
        raise ValueError("post-submission queue job IDs differ")
    submitted_job_ids = value["submitted_job_ids"]
    if (
        not isinstance(submitted_job_ids, list)
        or submitted_job_ids != sorted(set(submitted_job_ids))
        or any(JOB_ID_PATTERN.fullmatch(str(job_id)) is None
               for job_id in submitted_job_ids)
        or any(job_id not in submitted_job_ids for job_id in current_wave_job_ids)
    ):
        raise ValueError("post-submission submitted job IDs are invalid")
    jobs = value["live_jobs"]
    if not isinstance(jobs, list):
        raise ValueError("post-submission queue jobs are invalid")
    numeric = {
        key: value[key]
        for key in (
            "external_nodes", "wave_nodes", "running_nodes",
            "total_nodes_if_wave_runs",
        )
    }
    if (
        any(isinstance(number, bool) or not isinstance(number, int)
            for number in numeric.values())
        or numeric["external_nodes"] < 0
        or numeric["running_nodes"] < 0
        or numeric["wave_nodes"] != expected_wave_nodes
        or numeric["wave_nodes"] <= 0
        or numeric["total_nodes_if_wave_runs"] <= 0
    ):
        raise ValueError("post-submission queue node totals are invalid")
    seen = set()
    external_nodes = 0
    running_nodes = 0
    current_jobs: dict[str, dict[str, object]] = {}
    for item in jobs:
        job = require_exact_keys(
            item,
            {
                "job_id", "array_job_id", "array_task_id", "job_name", "state",
                "nodes", "owner",
            },
            "post-submission queue job",
        )
        job_id = str(job["job_id"])
        array_job_id = str(job["array_job_id"])
        array_task_id = str(job["array_task_id"])
        if (
            LIVE_JOB_ID_PATTERN.fullmatch(job_id) is None
            or JOB_ID_PATTERN.fullmatch(array_job_id) is None
            or job_id in seen
        ):
            raise ValueError("post-submission queue job identity is invalid")
        if array_task_id == "N/A":
            if job_id != array_job_id:
                raise ValueError("post-submission queue has ambiguous non-array identity")
        elif (
            re.fullmatch(r"[0-9]+", array_task_id) is None
            or job_id != f"{array_job_id}_{array_task_id}"
        ):
            raise ValueError(
                "post-submission queue contains a condensed or ambiguous array"
            )
        if (
            not isinstance(job["job_name"], str)
            or not job["job_name"]
            or not isinstance(job["state"], str)
            or job["state"] not in {
                "PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED",
            }
            or not isinstance(job["owner"], str)
            or not job["owner"]
        ):
            raise ValueError("post-submission queue job fields are invalid")
        seen.add(job_id)
        nodes = job["nodes"]
        if isinstance(nodes, bool) or not isinstance(nodes, int) or nodes <= 0:
            raise ValueError("post-submission queue node count is invalid")
        if job_id not in submitted_job_ids:
            external_nodes += nodes
            if str(job["job_name"]).startswith(QUALIFICATION_JOB_PREFIX):
                raise ValueError(
                    "post-submission queue contains another qualification root"
                )
        if job_id in current_wave_job_ids:
            if array_task_id != "N/A":
                raise ValueError("submitted qualification job became an array")
            current_jobs[job_id] = job
        if job["state"] != "PENDING":
            running_nodes += nodes
    expected_job_id, expected_job_name, expected_job_nodes = expected_current_job
    observed_current = current_jobs.get(expected_job_id)
    if (
        any(job_id not in seen for job_id in current_wave_job_ids)
        or observed_current is None
        or observed_current["job_name"] != expected_job_name
        or observed_current["nodes"] != expected_job_nodes
        or value["external_nodes"] != external_nodes
        or value["running_nodes"] != running_nodes
        or value["total_nodes_if_wave_runs"]
        != external_nodes + int(value["wave_nodes"])
        or value["total_nodes_if_wave_runs"] > MAX_CONCURRENT_NODES
        or running_nodes > MAX_CONCURRENT_NODES
    ):
        raise ValueError("post-submission queue audit arithmetic differs")


def validate_r17_submission_admission(record: object, root: Path) -> None:
    """Validate the retained exclusive account and Stage I admission boundary."""

    value = require_exact_keys(
        record,
        {
            "queried_utc", "external_nodes", "wave_nodes",
            "total_nodes_if_admitted", "live_jobs",
            "exclusive_operational_admission", "canonical_stage_i_state",
        },
        "R17 submission admission",
    )
    parse_utc_timestamp(value["queried_utc"], "R17 admission query time")
    if (
        value["external_nodes"] != 0
        or value["wave_nodes"] != R17_OPERATIONAL_NODES
        or value["total_nodes_if_admitted"] != R17_OPERATIONAL_NODES
        or value["live_jobs"] != []
        or value["exclusive_operational_admission"] is not True
    ):
        raise ValueError("R17 retained account admission was not exclusive")
    state = require_exact_keys(
        value["canonical_stage_i_state"],
        {
            "reservations", "reservation_count", "active_reservations",
            "transaction_stores", "transactions",
        },
        "R17 canonical Stage I admission",
    )
    expected_accounting = root / "accounting"
    expected_transactions = [
        str(expected_accounting / name) for name in STAGE_I_TRANSACTION_NAMES
    ]
    if (
        state["reservations"] != str(expected_accounting / STAGE_I_RESERVATIONS_NAME)
        or isinstance(state["reservation_count"], bool)
        or not isinstance(state["reservation_count"], int)
        or state["reservation_count"] < 0
        or state["active_reservations"] != 0
        or state["transaction_stores"] != expected_transactions
        or state["transactions"] != 0
    ):
        raise ValueError("R17 retained canonical Stage I admission was not empty")


def validate_job_binding(packet: dict[str, object], wave_number: int,
                         previous_job_ids: list[str],
                         qualification_root: Path,
                         current_wave_job_ids: list[str] | None = None,
                         expected_wave_nodes: int | None = None,
                         batch_script_runner=subprocess.run) -> dict[str, object]:
    """Authenticate one exact submitted job against its retained launch packet."""

    intent = packet["execution_intent"]
    path = require_owned_regular_file(
        Path(str(intent["paths"]["job_binding"])), qualification_root, "job binding"
    )
    binding = require_exact_keys(
        load_json(path, "job binding"),
        {
            "schema_version", "record_type", "submitted_utc", "wave", "job_id",
            "job_name", "nodes", "prior_wave_job_ids", "execution_intent_sha256",
            "execution_contract_sha256", "submission_argv", "prepared_wave",
            "submission_stdin_sha256", "submission_stdin_size_bytes",
            "sbatch_stdout", "sbatch_stderr", "submission_boundary",
            "submission_result", "prepared_manifest", "archived_input",
            "batch_script", "scheduler_batch_script", "submission_receipt",
            "post_submission_queue",
        },
        "job binding",
    )
    try:
        submitted = datetime.fromisoformat(str(binding["submitted_utc"]))
    except ValueError as error:
        raise ValueError("job binding has invalid submission time") from error
    if submitted.tzinfo is None:
        raise ValueError("job binding submission time lacks an explicit UTC offset")
    expected_argv = actual_submission_argv(packet, previous_job_ids)
    if (
        binding["schema_version"] != 4
        or binding["record_type"] != "cgl_lf_stage_i_qualification_job_binding"
        or binding["wave"] != wave_number
        or JOB_ID_PATTERN.fullmatch(str(binding["job_id"])) is None
        or binding["job_name"] != intent["job_name"]
        or binding["nodes"] != intent["allocation"]["nodes"]
        or binding["prior_wave_job_ids"] != previous_job_ids
        or binding["execution_intent_sha256"] != packet["execution_intent_sha256"]
        or binding["execution_contract_sha256"] != intent["execution_contract_sha256"]
        or binding["submission_argv"] != expected_argv
        or binding["submission_stdin_sha256"]
        != intent["authenticated_commands"]["batch_script_sha256"]
        or binding["submission_stdin_size_bytes"]
        != len(intent["authenticated_commands"]["batch_script_text"].encode("utf-8"))
        or parse_job_id(str(binding["sbatch_stdout"])) != str(binding["job_id"])
    ):
        raise ValueError("job binding differs from authenticated submission")
    expected_records = {
        "prepared_wave": qualification_root / "prepared_wave.json",
        "prepared_manifest": Path(str(intent["paths"]["prepared_manifest"])),
        "archived_input": Path(str(intent["paths"]["archived_input"])),
        "batch_script": Path(str(intent["paths"]["batch_script"])),
        "scheduler_batch_script": Path(str(intent["paths"]["scheduler_batch_script"])),
        "submission_receipt": Path(str(intent["paths"]["submission_receipt"])),
    }
    for key, expected in expected_records.items():
        observed = validate_retained_file(binding[key], qualification_root, key)
        if observed != expected.resolve(strict=True):
            raise ValueError(f"job binding {key} path differs")
    boundary_path = validate_retained_file(
        binding["submission_boundary"], qualification_root, "submission boundary"
    )
    result_path = validate_retained_file(
        binding["submission_result"], qualification_root, "submission result"
    )
    if boundary_path.parent != qualification_root / "submission-journal":
        raise ValueError("job binding submission boundary path differs")
    if result_path.parent != qualification_root / "submission-journal":
        raise ValueError("job binding submission result path differs")
    boundary = load_json(boundary_path, "submission boundary")
    result = load_json(result_path, "submission result")
    receipt_path = validate_retained_file(
        binding["submission_receipt"], qualification_root, "submission receipt"
    )
    receipt, _ = validate_submission_receipt(
        load_json(receipt_path, "submission receipt"),
        packet,
        qualification_root,
        expected_job_id=str(binding["job_id"]),
        expected_submission_argv=expected_argv,
    )
    recovery_path = Path(str(intent["paths"]["submission_recovery"]))
    if recovery_path.exists() or recovery_path.is_symlink():
        raise ValueError("bound qualification job has a submission recovery record")
    current_ids = [*(current_wave_job_ids or []), str(binding["job_id"])]
    if expected_wave_nodes is None:
        raise ValueError("job binding validation lacks the reviewed wave node total")
    validate_post_submission_queue(
        binding["post_submission_queue"],
        current_ids,
        expected_wave_nodes=expected_wave_nodes,
        expected_current_job=(
            str(binding["job_id"]),
            str(intent["job_name"]),
            int(intent["allocation"]["nodes"]),
        ),
    )
    if intent["case_id"] == R17_CASE_ID:
        queue = binding["post_submission_queue"]
        if (
            queue["external_nodes"] != 0
            or queue["submitted_job_ids"] != [str(binding["job_id"])]
            or queue["current_wave_job_ids"] != [str(binding["job_id"])]
            or len(queue["live_jobs"]) != 1
        ):
            raise ValueError("R17 post-submission queue was not exclusive")
        validate_r17_submission_admission(
            boundary.get("queue_admission"),
            Path(str(intent["project_root"])),
        )
    common = {
        "wave": wave_number,
        "execution_intent_sha256": packet["execution_intent_sha256"],
        "execution_contract_sha256": intent["execution_contract_sha256"],
        "submission_argv": expected_argv,
        "submission_stdin_sha256": binding["submission_stdin_sha256"],
        "submission_stdin_size_bytes": binding["submission_stdin_size_bytes"],
    }
    if any(boundary.get(key) != value for key, value in common.items()):
        raise ValueError("submission boundary differs from authenticated job")
    if any(result.get(key) != value for key, value in common.items()):
        raise ValueError("submission result differs from authenticated job")
    if (
        result.get("job_id") != binding["job_id"]
        or result.get("sbatch_stdout") != binding["sbatch_stdout"]
        or result.get("sbatch_stderr") != binding["sbatch_stderr"]
        or result.get("submission_receipt") != binding["submission_receipt"]
        or result.get("scheduler_batch_script") != binding["scheduler_batch_script"]
        or result.get("post_submission_queue") != binding["post_submission_queue"]
        or result.get("boundary_sha256") != sha256(boundary_path)
        or result.get("previous_record_sha256") != sha256(boundary_path)
    ):
        raise ValueError("submission result does not bind the exact sbatch response")
    if (
        receipt["wave"] != wave_number
        or receipt["submission_boundary"] != binding["submission_boundary"]
    ):
        raise ValueError("submission receipt does not bind the exact submission boundary")
    validate_materialized_packet(packet, qualification_root)
    scheduler_path = validate_retained_file(
        binding["scheduler_batch_script"],
        qualification_root,
        "scheduler-stored batch script",
    )
    expected_script = intent["authenticated_commands"]["batch_script_text"].encode("utf-8")
    if (
        read_regular_bytes(
            scheduler_path, "scheduler-stored batch script", single_link=True
        ) != expected_script
        or collect_live_batch_script(
            str(binding["job_id"]), batch_script_runner
        ) != expected_script
    ):
        raise ValueError("Slurm-stored batch script differs from authenticated job")
    return binding


def load_bound_evidence_file(record: object, qualification_root: Path,
                             label: str) -> Path:
    """Validate one exact retained evidence record."""

    value = require_exact_keys(record, {"path", "sha256", "size_bytes"}, label)
    return validate_retained_file(value, qualification_root, label)


SACCT_FIELDS = (
    "JobIDRaw", "JobName", "State", "ExitCode", "NNodes", "ElapsedRaw",
    "Submit", "Start", "End", "Partition", "Account", "AllocTRES",
    "ReqTRES", "TimelimitRaw", "Comment", "SubmitLine",
)
SCHEDULER_HEADER = "|".join(SACCT_FIELDS)
ACCOUNT_SACCT_FIELDS = (
    "JobIDRaw", "JobName", "State", "ExitCode", "NNodes", "ElapsedRaw",
    "Submit", "Start", "End", "Partition", "Account", "User",
)
ACCOUNT_SCHEDULER_HEADER = "|".join(ACCOUNT_SACCT_FIELDS)
ACCOUNT_MISSING_TIMESTAMPS = frozenset({"", "Unknown", "N/A", "None"})


def parse_scheduler_time(value: str, label: str) -> datetime:
    """Parse one retained Slurm timestamp."""

    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise ValueError(f"{label} is not an ISO timestamp") from error
    if parsed.tzinfo is None:
        raise ValueError(f"{label} must include an explicit UTC offset")
    return parsed


def parse_tres(value: str, label: str) -> dict[str, str]:
    """Parse one unambiguous comma-delimited Slurm TRES value."""

    result = {}
    for item in value.split(","):
        if not item:
            continue
        key, separator, item_value = item.partition("=")
        if not separator or not key or not item_value or key in result:
            raise ValueError(f"{label} is ambiguous")
        result[key] = item_value
    return result


def parse_scheduler_text(text: str) -> dict[str, object]:
    """Parse exactly one top-level pipe-delimited sacct record."""

    lines = [line for line in text.splitlines() if line]
    if len(lines) != 2 or lines[0] != SCHEDULER_HEADER:
        raise ValueError("scheduler raw evidence must contain one exact sacct record")
    fields = lines[1].split("|")
    if len(fields) != len(SACCT_FIELDS):
        raise ValueError("scheduler raw evidence has invalid columns")
    (
        job_id, job_name, state_value, exit_code, nodes_text, elapsed_text,
        submit, start, end, partition, account, alloc_tres_text, req_tres_text,
        timelimit_text, comment, submit_line,
    ) = fields
    state = state_value.split()[0].split("+")[0]
    if JOB_ID_PATTERN.fullmatch(job_id) is None:
        raise ValueError("scheduler raw evidence has invalid job ID")
    try:
        nodes = int(nodes_text)
        elapsed = int(elapsed_text)
        timelimit_minutes = int(timelimit_text)
    except ValueError as error:
        raise ValueError("scheduler raw evidence has invalid numeric values") from error
    if nodes <= 0 or elapsed <= 0 or timelimit_minutes <= 0:
        raise ValueError("scheduler raw evidence has nonpositive numeric values")
    submit_time = parse_scheduler_time(submit, "scheduler submit time")
    start_time = parse_scheduler_time(start, "scheduler start time")
    end_time = parse_scheduler_time(end, "scheduler end time")
    if not submit_time <= start_time < end_time:
        raise ValueError("scheduler raw evidence has invalid time ordering")
    if abs((end_time - start_time).total_seconds() - elapsed) > 1.0:
        raise ValueError("scheduler elapsed seconds differ from retained interval")
    return {
        "job_id": job_id,
        "job_name": job_name,
        "state": state,
        "exit_code": exit_code,
        "nodes": nodes,
        "elapsed_seconds": elapsed,
        "submit_utc": submit,
        "start_utc": start,
        "end_utc": end,
        "partition": partition,
        "account": account,
        "alloc_tres": parse_tres(alloc_tres_text, "scheduler AllocTRES"),
        "req_tres": parse_tres(req_tres_text, "scheduler ReqTRES"),
        "timelimit_minutes": timelimit_minutes,
        "comment": comment,
        "submit_line": submit_line,
    }


def parse_scheduler_raw(path: Path) -> dict[str, object]:
    """Parse one descriptor-authenticated retained sacct record."""

    try:
        text = read_regular_bytes(
            path, "scheduler raw evidence", single_link=True
        ).decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("scheduler raw evidence is not UTF-8") from error
    return parse_scheduler_text(text)


def collect_live_scheduler_raw(job_id: str,
                               runner=subprocess.run) -> str:
    """Collect one live completed top-level Slurm record with exact timestamps."""

    if JOB_ID_PATTERN.fullmatch(job_id) is None:
        raise ValueError("live scheduler query has invalid job ID")
    argv = [
        str(SACCT), "-j", job_id, "-X", "-n", "-P",
        "-o", ",".join(SACCT_FIELDS),
    ]
    completed = runner(
        argv, check=True, capture_output=True, text=True,
        env=scheduler_environment(exact_timestamps=True),
    )
    stdout = completed.stdout
    if isinstance(stdout, bytes):
        try:
            stdout = stdout.decode("utf-8")
        except UnicodeDecodeError as error:
            raise ValueError("live sacct output is not UTF-8") from error
    lines = [line for line in str(stdout).splitlines() if line]
    if len(lines) != 1:
        raise ValueError("live sacct query must return exactly one top-level job record")
    text = SCHEDULER_HEADER + "\n" + lines[0] + "\n"
    parsed = parse_scheduler_text(text)
    if parsed["job_id"] != job_id:
        raise ValueError("live sacct query returned a different job ID")
    return text


def parse_optional_scheduler_time(value: str, label: str) -> datetime | None:
    """Parse one optional account-wide Slurm timestamp."""

    if value in ACCOUNT_MISSING_TIMESTAMPS:
        return None
    return parse_scheduler_time(value, label).astimezone(timezone.utc)


def parse_account_scheduler_text(text: str) -> list[dict[str, object]]:
    """Parse a complete all-users top-level Slurm account query."""

    lines = text.splitlines()
    if not lines or lines[0] != ACCOUNT_SCHEDULER_HEADER:
        raise ValueError("account scheduler raw evidence has an invalid header")
    records = []
    seen = set()
    for line in lines[1:]:
        if not line:
            continue
        fields = line.split("|")
        if len(fields) != len(ACCOUNT_SACCT_FIELDS):
            raise ValueError("account scheduler raw evidence has invalid columns")
        (
            job_id, job_name, state_value, exit_code, nodes_text, elapsed_text,
            submit, start, end, partition, account, owner,
        ) = fields
        state_parts = state_value.split()
        state = state_parts[0].split("+")[0] if state_parts else ""
        if (
            LIVE_JOB_ID_PATTERN.fullmatch(job_id) is None
            or job_id in seen
            or not job_name
            or not state
            or not exit_code
            or not partition
            or str(account).casefold() != ACCOUNT.casefold()
            or not owner
        ):
            raise ValueError("account scheduler raw evidence has invalid job identity")
        try:
            nodes = int(nodes_text)
            elapsed = int(elapsed_text)
        except ValueError as error:
            raise ValueError(
                "account scheduler raw evidence has invalid numeric values"
            ) from error
        if nodes < 0 or elapsed < 0:
            raise ValueError("account scheduler raw evidence has negative numeric values")
        submit_time = parse_optional_scheduler_time(
            submit, "account scheduler submit time"
        )
        start_time = parse_optional_scheduler_time(
            start, "account scheduler start time"
        )
        end_time = parse_optional_scheduler_time(end, "account scheduler end time")
        if (
            submit_time is None
            or (
                start_time is None
                and (
                    elapsed != 0
                    or state in {"RUNNING", "COMPLETING", "SUSPENDED", "COMPLETED"}
                    or (end_time is not None and submit_time > end_time)
                )
            )
            or (start_time is not None and nodes <= 0)
            or (
                start_time is not None
                and submit_time > start_time
            )
            or (
                start_time is not None
                and end_time is not None
                and (
                    end_time < start_time
                    or abs((end_time - start_time).total_seconds() - elapsed) > 1.0
                )
            )
        ):
            raise ValueError("account scheduler raw evidence has invalid time ordering")
        seen.add(job_id)
        records.append({
            "job_id": job_id,
            "job_name": job_name,
            "state": state,
            "exit_code": exit_code,
            "nodes": nodes,
            "elapsed_seconds": elapsed,
            "submit_utc": submit,
            "start_utc": start,
            "end_utc": end,
            "partition": partition,
            "account": account,
            "owner": owner,
        })
    if not records:
        raise ValueError("account scheduler raw evidence is empty")
    return sorted(records, key=lambda record: str(record["job_id"]))


def account_scheduler_query_contract(start: datetime,
                                     end: datetime) -> dict[str, object]:
    """Return the exact all-users account query spanning one execution interval."""

    if start.tzinfo is None or end.tzinfo is None:
        raise ValueError("account scheduler query interval lacks an explicit offset")
    if not start < end:
        raise ValueError("account scheduler query interval is invalid")
    query_start = start - timedelta(seconds=1)
    query_end = end + timedelta(seconds=1)
    return {
        "account": ACCOUNT,
        "all_users": True,
        "allocations_only": True,
        "expanded_arrays": True,
        "start_utc": format_utc_timestamp(query_start.astimezone(timezone.utc)),
        "end_utc": format_utc_timestamp(query_end.astimezone(timezone.utc)),
        "start_argument": query_start.strftime("%Y-%m-%dT%H:%M:%S"),
        "end_argument": query_end.strftime("%Y-%m-%dT%H:%M:%S"),
        "start_scheduler_offset": query_start.strftime("%z"),
        "end_scheduler_offset": query_end.strftime("%z"),
        "fields": list(ACCOUNT_SACCT_FIELDS),
    }


def collect_live_account_visibility(runner=subprocess.run) -> dict[str, object]:
    """Require Slurm accounting to expose every user's job allocation."""

    argv = [str(SCONTROL), "show", "config"]
    completed = runner(
        argv, check=True, capture_output=True, text=True, env=scheduler_environment()
    )
    stdout = completed.stdout
    if isinstance(stdout, bytes):
        try:
            stdout = stdout.decode("utf-8")
        except UnicodeDecodeError as error:
            raise ValueError("live Slurm configuration is not UTF-8") from error
    matches = re.findall(r"(?m)^\s*PrivateData\s*=\s*(\S.*?)\s*$", str(stdout))
    if len(matches) != 1:
        raise ValueError("live Slurm PrivateData configuration is ambiguous")
    private_data = matches[0]
    settings = {
        value.strip().casefold()
        for value in private_data.split(",")
        if value.strip()
    }
    if not settings:
        raise ValueError("live Slurm PrivateData configuration is empty")
    if "all" in settings or "jobs" in settings:
        raise ValueError("Slurm PrivateData prevents account-wide job evidence")
    return {
        "private_data": private_data,
        "all_users_job_visibility": True,
    }


def collect_live_account_scheduler_raw(start: datetime, end: datetime,
                                       runner=subprocess.run) -> tuple[str, dict[str, object]]:
    """Collect every top-level account job that can intersect one interval."""

    contract = account_scheduler_query_contract(start, end)
    argv = [
        str(SACCT), "-a", "-A", ACCOUNT, "-X", "--array", "-n", "-P",
        "-S", str(contract["start_argument"]), "-E", str(contract["end_argument"]),
        "-o", ",".join(ACCOUNT_SACCT_FIELDS),
    ]
    completed = runner(
        argv, check=True, capture_output=True, text=True,
        env=scheduler_environment(exact_timestamps=True),
    )
    stdout = completed.stdout
    if isinstance(stdout, bytes):
        try:
            stdout = stdout.decode("utf-8")
        except UnicodeDecodeError as error:
            raise ValueError("live account sacct output is not UTF-8") from error
    rows = [line for line in str(stdout).splitlines() if line]
    text = ACCOUNT_SCHEDULER_HEADER + "\n" + "\n".join(rows) + "\n"
    parse_account_scheduler_text(text)
    return text, contract


def build_r17_account_exclusivity_evidence(
    scheduler: dict[str, object],
    raw_text: str,
    query_contract: dict[str, object],
    visibility_contract: dict[str, object],
    measured_utc: str,
) -> dict[str, object]:
    """Prove that no other account allocation overlapped the R17 execution."""

    measured = parse_utc_timestamp(measured_utc, "R17 account-exclusivity measurement")
    scheduler_start = parse_scheduler_time(
        str(scheduler["start_utc"]), "R17 scheduler start"
    )
    scheduler_end = parse_scheduler_time(
        str(scheduler["end_utc"]), "R17 scheduler end"
    )
    start = scheduler_start.astimezone(timezone.utc)
    end = scheduler_end.astimezone(timezone.utc)
    expected_contract = account_scheduler_query_contract(
        scheduler_start, scheduler_end
    )
    if (
        query_contract != expected_contract
        or visibility_contract.get("all_users_job_visibility") is not True
        or set(visibility_contract) != {
            "private_data", "all_users_job_visibility",
        }
        or not isinstance(visibility_contract["private_data"], str)
        or not visibility_contract["private_data"]
        or any(
            value in {"all", "jobs"}
            for value in (
                item.strip().casefold()
                for item in visibility_contract["private_data"].split(",")
            )
        )
        or measured < end
    ):
        raise ValueError("R17 account scheduler query contract or chronology differs")
    records = parse_account_scheduler_text(raw_text)
    target_jobs = [
        record for record in records if record["job_id"] == scheduler["job_id"]
    ]
    if len(target_jobs) != 1:
        raise ValueError("R17 account scheduler evidence lacks one exact qualification job")
    target = target_jobs[0]
    if (
        target["job_name"] != scheduler["job_name"]
        or target["state"] != scheduler["state"]
        or target["exit_code"] != scheduler["exit_code"]
        or target["nodes"] != scheduler["nodes"]
        or target["elapsed_seconds"] != scheduler["elapsed_seconds"]
        or target["submit_utc"] != scheduler["submit_utc"]
        or target["start_utc"] != scheduler["start_utc"]
        or target["end_utc"] != scheduler["end_utc"]
        or target["partition"] != scheduler["partition"]
        or str(target["account"]).casefold()
        != str(scheduler["account"]).casefold()
    ):
        raise ValueError("R17 account scheduler evidence differs from qualification job")
    overlapping = []
    for record in records:
        record_start = parse_optional_scheduler_time(
            str(record["start_utc"]), "account scheduler start time"
        )
        record_end = parse_optional_scheduler_time(
            str(record["end_utc"]), "account scheduler end time"
        )
        if (
            record_start is not None
            and record_start < end
            and (record_end is None or record_end > start)
        ):
            overlapping.append(str(record["job_id"]))
    if overlapping != [str(scheduler["job_id"])]:
        raise ValueError(
            "R17 account scheduler evidence contains an overlapping account job"
        )
    target_scheduler = {
        key: scheduler[key]
        for key in (
            "job_id", "job_name", "state", "exit_code", "nodes",
            "elapsed_seconds", "submit_utc", "start_utc", "end_utc",
            "partition", "account",
        )
    }
    return {
        "schema_version": 1,
        "record_type": "stage-i-r17-account-exclusivity-evidence",
        "execution_epoch": EXECUTION_EPOCH,
        "measured_utc": measured_utc,
        "query_contract": query_contract,
        "visibility_contract": visibility_contract,
        "raw_account_scheduler_sha256": hashlib.sha256(
            raw_text.encode("utf-8")
        ).hexdigest(),
        "qualification_job": target_scheduler,
        "qualification_job_sha256": stable_json_sha256(target_scheduler),
        "account_jobs": records,
        "account_jobs_sha256": stable_json_sha256(records),
        "overlapping_job_ids": overlapping,
        "exclusive_entire_execution_interval": True,
    }


def validate_r17_account_exclusivity_evidence(
    value: object,
    scheduler: dict[str, object],
    raw_text: str,
) -> dict[str, object]:
    """Reproduce one retained full-interval R17 account-exclusivity proof."""

    evidence = require_exact_keys(
        value,
        {
            "schema_version", "record_type", "execution_epoch", "measured_utc",
            "query_contract", "visibility_contract",
            "raw_account_scheduler_sha256", "qualification_job",
            "qualification_job_sha256", "account_jobs", "account_jobs_sha256",
            "overlapping_job_ids", "exclusive_entire_execution_interval",
        },
        "R17 account exclusivity evidence",
    )
    reproduced = build_r17_account_exclusivity_evidence(
        scheduler,
        raw_text,
        evidence["query_contract"],
        evidence["visibility_contract"],
        str(evidence["measured_utc"]),
    )
    if (
        evidence["schema_version"] != 1
        or evidence["record_type"] != "stage-i-r17-account-exclusivity-evidence"
        or evidence["execution_epoch"] != EXECUTION_EPOCH
        or evidence != reproduced
    ):
        raise ValueError("R17 account exclusivity evidence differs")
    return evidence


def build_scheduler_evidence(packet: dict[str, object],
                             binding: dict[str, object],
                             raw_path: Path,
                             qualification_root: Path) -> dict[str, object]:
    """Build authenticated scheduler evidence from one retained sacct record."""

    raw = require_owned_regular_file(raw_path, qualification_root, "scheduler raw evidence")
    parsed = parse_scheduler_raw(raw)
    intent = packet["execution_intent"]
    if (
        parsed["job_id"] != binding["job_id"]
        or parsed["job_name"] != intent["job_name"]
        or parsed["nodes"] != intent["allocation"]["nodes"]
        or parsed["state"] != "COMPLETED"
        or parsed["exit_code"] != "0:0"
        or parsed["elapsed_seconds"] > intent["allocation"]["walltime_seconds"]
        or parsed["partition"] != PARTITION
        or str(parsed["account"]).casefold() != ACCOUNT.casefold()
        or parsed["timelimit_minutes"]
        != intent["allocation"]["walltime_seconds"] // 60
        or parsed["comment"] != intent["execution_contract_sha256"]
        or shlex.split(str(parsed["submit_line"])) != binding["submission_argv"]
        or parsed["alloc_tres"].get("node") != str(intent["allocation"]["nodes"])
        or parsed["req_tres"].get("node") != str(intent["allocation"]["nodes"])
    ):
        raise ValueError("scheduler raw evidence differs from authenticated job binding")
    return {
        "schema_version": 3,
        "record_type": "cgl_lf_stage_i_qualification_scheduler_evidence",
        **parsed,
        "execution_intent_sha256": packet["execution_intent_sha256"],
        "execution_contract_sha256": intent["execution_contract_sha256"],
        "job_binding_sha256": stable_json_sha256(binding),
        "raw_scheduler": retained_file(raw, "scheduler raw evidence"),
    }


def validate_scheduler_evidence(path: Path, *, packet: dict[str, object],
                                binding: dict[str, object],
                                qualification_root: Path,
                                runner=subprocess.run) -> dict[str, object]:
    """Reproduce retained evidence and bind it to a fresh live Slurm query."""

    retained = load_json(path, "qualification scheduler evidence")
    record = require_exact_keys(
        retained,
        {
            "schema_version", "record_type", "job_id", "job_name", "state",
            "exit_code", "nodes", "elapsed_seconds", "submit_utc", "start_utc",
            "end_utc", "partition", "account", "alloc_tres", "req_tres",
            "timelimit_minutes", "comment", "submit_line",
            "execution_intent_sha256", "execution_contract_sha256",
            "job_binding_sha256", "raw_scheduler",
        },
        "qualification scheduler evidence",
    )
    raw_path = validate_retained_file(
        record["raw_scheduler"], qualification_root, "raw scheduler evidence"
    )
    reproduced = build_scheduler_evidence(
        packet, binding, raw_path, qualification_root
    )
    if record != reproduced:
        raise ValueError("qualification scheduler evidence differs from retained sacct data")
    live = parse_scheduler_text(
        collect_live_scheduler_raw(str(binding["job_id"]), runner)
    )
    retained_scheduler = {
        key: value for key, value in record.items()
        if key not in {
            "schema_version", "record_type", "execution_intent_sha256",
            "execution_contract_sha256", "job_binding_sha256", "raw_scheduler",
        }
    }
    if live != retained_scheduler:
        raise ValueError("retained scheduler evidence differs from live Slurm")
    return record


def parse_history_bytes(payload: bytes, label: str) -> dict[str, list[float]]:
    """Parse one complete finite Athena history payload."""

    labels: list[str] = []
    rows: list[list[float]] = []
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} is not UTF-8 history data") from error
    for line in text.splitlines():
        if line.startswith("#"):
            found = re.findall(r"\[\d+\]=([^\s]+)", line)
            if found:
                labels = found
        elif line.strip():
            try:
                row = [float(value) for value in line.split()]
            except ValueError as error:
                raise ValueError(f"{label} has a nonnumeric row") from error
            rows.append(row)
    if (
        not labels
        or len(labels) != len(set(labels))
        or not rows
        or any(len(row) != len(labels) for row in rows)
    ):
        raise ValueError(f"{label} has incomplete labels or rows")
    if any(not math.isfinite(value) for row in rows for value in row):
        raise ValueError(f"{label} contains non-finite values")
    return {
        name: [row[index] for row in rows] for index, name in enumerate(labels)
    }


def authenticated_history(path: Path, qualification_root: Path,
                          label: str) -> tuple[dict[str, list[float]], dict[str, object]]:
    """Read and hash one unchanged retained history payload."""

    path = require_owned_regular_file(path, qualification_root, label)
    payload = read_regular_bytes(path, label, single_link=True)
    return parse_history_bytes(payload, label), {
        "path": str(path),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": len(payload),
    }


def expected_history_times(target: float, cadence: float) -> list[float]:
    """Return the exact reviewed start/cadence/final-time history schedule."""

    count = int(math.floor(target / cadence + HISTORY_CADENCE_TOLERANCE))
    result = [index * cadence for index in range(count + 1)]
    if abs(result[-1] - target) > HISTORY_CADENCE_TOLERANCE:
        result.append(target)
    else:
        result[-1] = target
    return result


def require_history_schedule(history: dict[str, list[float]], target: float,
                             cadence: float, label: str) -> None:
    """Require an exact zero start, strict monotonicity, cadence, and endpoint."""

    times = history["time"]
    expected = expected_history_times(target, cadence)
    if (
        len(times) != len(expected)
        or any(later <= earlier for earlier, later in zip(times, times[1:]))
        or any(
            abs(observed - required) > HISTORY_CADENCE_TOLERANCE
            for observed, required in zip(times, expected)
        )
    ):
        raise ValueError(f"{label} violates the exact start/cadence schedule")


def require_monotonic_columns(history: dict[str, list[float]],
                              names: tuple[str, ...], label: str) -> None:
    """Require reviewed cumulative diagnostic columns to never decrease."""

    for name in names:
        if any(
            later < earlier
            for earlier, later in zip(history[name], history[name][1:])
        ):
            raise ValueError(f"{label} cumulative column is nonmonotonic: {name}")


def read_binary_line(stream, path: Path, label: str) -> bytes:
    """Read one bounded newline-terminated Athena binary header line."""

    line = stream.readline(MAX_BINARY_HEADER_BYTES + 1)
    if not line or len(line) > MAX_BINARY_HEADER_BYTES or not line.endswith(b"\n"):
        raise ValueError(f"{label} is missing or oversized: {path}")
    return line


def binary_header_assignment(stream, path: Path, expected: str) -> str:
    """Read one exact UTF-8 ``key=value`` binary-header assignment."""

    try:
        line = read_binary_line(stream, path, f"Athena snapshot {expected}").decode(
            "utf-8"
        ).strip()
    except UnicodeDecodeError as error:
        raise ValueError(f"Athena snapshot header is not UTF-8: {path}") from error
    key, separator, value = line.partition("=")
    if not separator or key.strip() != expected or not value.strip():
        raise ValueError(f"invalid Athena snapshot {expected} header: {path}")
    return value.strip()


def snapshot_mesh_nghost(header: bytes, path: Path) -> int:
    """Read one unique nonnegative ``mesh/nghost`` parameter."""

    try:
        text = header.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"Athena snapshot parameter header is not UTF-8: {path}") from error
    block = ""
    values = []
    for original in text.splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            continue
        if block == "mesh" and "=" in line:
            key, value = line.split("=", 1)
            if key.strip() == "nghost":
                values.append(value.strip())
    if len(values) != 1:
        raise ValueError(f"Athena snapshot must contain one mesh/nghost value: {path}")
    try:
        nghost = int(values[0])
    except ValueError as error:
        raise ValueError(f"Athena snapshot mesh/nghost is invalid: {path}") from error
    if nghost < 0:
        raise ValueError(f"Athena snapshot mesh/nghost is negative: {path}")
    return nghost


def read_exact_binary(stream, size: int, path: Path, label: str) -> bytes:
    """Read one exact binary extent."""

    payload = stream.read(size)
    if len(payload) != size:
        raise ValueError(f"{label} is truncated: {path}")
    return payload


def structural_snapshot_profile(path: Path,
                                intent: dict[str, object] | None = None
                                ) -> dict[str, object]:
    """Stream-validate one complete physically admissible Athena snapshot."""

    with path.open("rb") as stream:
        file_size = os.fstat(stream.fileno()).st_size
        code_header = read_binary_line(
            stream, path, "Athena snapshot code header"
        ).split()
        if (
            not code_header
            or code_header[0] != b"Athena"
            or code_header[-1] != b"version=1.1"
        ):
            raise ValueError(f"invalid Athena v1.1 binary snapshot: {path}")
        try:
            preheader_count = int(
                binary_header_assignment(stream, path, "size of preheader")
            )
        except ValueError as error:
            raise ValueError(f"invalid Athena snapshot preheader: {path}") from error
        if not 2 <= preheader_count <= MAX_BINARY_PREHEADER_LINES:
            raise ValueError(f"invalid Athena snapshot preheader count: {path}")
        preheader: dict[str, str] = {}
        for _ in range(preheader_count - 1):
            try:
                line = read_binary_line(
                    stream, path, "Athena snapshot preheader"
                ).decode("utf-8").strip()
            except UnicodeDecodeError as error:
                raise ValueError(f"Athena snapshot preheader is not UTF-8: {path}") from error
            key, separator, value = line.partition("=")
            key = key.strip()
            if not separator or not key or not value.strip() or key in preheader:
                raise ValueError(f"invalid Athena snapshot preheader assignment: {path}")
            preheader[key] = value.strip()
        if set(preheader) != {
            "time", "cycle", "size of location", "size of variable",
        }:
            raise ValueError(f"Athena snapshot preheader metadata differs: {path}")
        try:
            snapshot_time = float(preheader["time"])
            cycle = int(preheader["cycle"])
            location_size = int(preheader["size of location"])
            variable_size = int(preheader["size of variable"])
        except (KeyError, ValueError) as error:
            raise ValueError(f"incomplete Athena snapshot preheader: {path}") from error
        if (
            not math.isfinite(snapshot_time)
            or cycle < 0
            or location_size not in (4, 8)
            or variable_size not in (4, 8)
        ):
            raise ValueError(f"invalid Athena snapshot preheader values: {path}")
        try:
            number_variables = int(
                binary_header_assignment(stream, path, "number of variables")
            )
        except ValueError as error:
            raise ValueError(f"invalid Athena snapshot variable count: {path}") from error
        try:
            variable_line = read_binary_line(
                stream, path, "Athena snapshot variable list"
            ).decode("utf-8").strip()
        except UnicodeDecodeError as error:
            raise ValueError(f"Athena snapshot variable list is not UTF-8: {path}") from error
        fields = variable_line.split()
        if (
            number_variables <= 0
            or number_variables > 4096
            or not fields
            or fields[0] != "variables:"
            or len(fields[1:]) != number_variables
            or len(set(fields[1:])) != number_variables
        ):
            raise ValueError(f"invalid Athena snapshot variable inventory: {path}")
        if tuple(fields[1:]) != REQUIRED_SNAPSHOT_VARIABLES:
            raise ValueError(f"Athena snapshot lacks the required variable inventory: {path}")
        try:
            header_size = int(binary_header_assignment(stream, path, "header offset"))
        except ValueError as error:
            raise ValueError(f"invalid Athena snapshot parameter-header size: {path}") from error
        if not 0 < header_size <= MAX_BINARY_HEADER_BYTES:
            raise ValueError(f"invalid Athena snapshot parameter-header size: {path}")
        header = read_exact_binary(
            stream, header_size, path, "Athena snapshot parameter header"
        )
        nghost = snapshot_mesh_nghost(header, path)
        try:
            parameter_values = parse_athinput_text(
                header.decode("utf-8"), "Athena snapshot parameter header"
            )
        except UnicodeDecodeError as error:
            raise ValueError(
                f"Athena snapshot parameter header is not UTF-8: {path}"
            ) from error
        parameter_contract_sha256 = None
        if intent is not None:
            expected_parameters = expected_case_configuration(
                str(intent["case_id"]),
                basename=str(intent["run_basename"]),
                tlim=str(intent["target_time"]),
            )
            mismatches = {
                key: {
                    "expected": expected_value,
                    "observed": parameter_values.get(key),
                }
                for key, expected_value in expected_parameters.items()
                if parameter_values.get(key) != expected_value
            }
            if mismatches:
                raise ValueError(
                    "Athena snapshot parameter header violates the reviewed "
                    f"contract: {mismatches}"
                )
            parameter_contract_sha256 = stable_json_sha256(expected_parameters)
        meshblocks = 0
        total_values = 0
        dimensions = set()
        logical_locations = set()
        meshblock_metadata = []
        positive_minima = {
            "dens": math.inf,
            "eint": math.inf,
            "p_perp": math.inf,
        }
        variable_indices = {
            name: fields[1:].index(name) for name in REQUIRED_SNAPSHOT_VARIABLES
        }
        variable_sums = {name: 0.0 for name in REQUIRED_SNAPSHOT_VARIABLES}
        cell_count = 0
        hard_bound_violation_cells = 0
        minimum_mirror_hard_margin = math.inf
        minimum_firehose_hard_margin = math.inf
        geometry_format = "<6d" if location_size == 8 else "<6f"
        while stream.tell() < file_size:
            indices = struct.unpack(
                "<6i",
                read_exact_binary(stream, 24, path, "Athena snapshot meshblock index"),
            )
            shape = (
                indices[1] - indices[0] + 1,
                indices[3] - indices[2] + 1,
                indices[5] - indices[4] + 1,
            )
            if any(value <= 0 for value in shape):
                raise ValueError(f"Athena snapshot meshblock dimensions are invalid: {path}")
            logical = struct.unpack(
                "<4i",
                read_exact_binary(
                    stream, 16, path, "Athena snapshot meshblock logical index"
                ),
            )
            if (
                logical[3] != 0
                or any(value < 0 for value in logical[:3])
                or logical in logical_locations
            ):
                raise ValueError(f"Athena snapshot logical meshblock is invalid: {path}")
            logical_locations.add(logical)
            geometry = struct.unpack(
                geometry_format,
                read_exact_binary(
                    stream,
                    6 * location_size,
                    path,
                    "Athena snapshot meshblock geometry",
                ),
            )
            if (
                any(not math.isfinite(value) for value in geometry)
                or any(
                    geometry[index + 1] <= geometry[index]
                    for index in (0, 2, 4)
                )
            ):
                raise ValueError(f"Athena snapshot meshblock geometry is non-finite: {path}")
            values = math.prod(shape) * number_variables
            data_bytes = values * variable_size
            data = read_exact_binary(
                stream, data_bytes, path, "Athena snapshot meshblock data"
            )
            decoded = array("f" if variable_size == 4 else "d")
            if decoded.itemsize != variable_size:
                raise ValueError("host floating-point ABI cannot decode Athena snapshot")
            decoded.frombytes(data)
            if len(decoded) != values or any(not math.isfinite(value) for value in decoded):
                raise ValueError(f"Athena snapshot payload contains non-finite values: {path}")
            cells = math.prod(shape)
            for name in positive_minima:
                variable = variable_indices[name]
                minimum = min(decoded[variable * cells:(variable + 1) * cells])
                if minimum <= 0.0:
                    raise ValueError(
                        f"Athena snapshot contains nonpositive {name}: {path}"
                    )
                positive_minima[name] = min(positive_minima[name], minimum)
            variables = {
                name: decoded[index * cells:(index + 1) * cells]
                for name, index in variable_indices.items()
            }
            for name, values_for_variable in variables.items():
                variable_sums[name] += math.fsum(values_for_variable)
            for ppar, pperp, bcc1, bcc2, bcc3 in zip(
                variables["eint"],
                variables["p_perp"],
                variables["bcc1"],
                variables["bcc2"],
                variables["bcc3"],
            ):
                bsqr = bcc1 * bcc1 + bcc2 * bcc2 + bcc3 * bcc3
                paniso = pperp - ppar
                mirror_margin = bsqr - paniso
                firehose_margin = paniso + 1.5 * bsqr
                minimum_mirror_hard_margin = min(
                    minimum_mirror_hard_margin, mirror_margin
                )
                minimum_firehose_hard_margin = min(
                    minimum_firehose_hard_margin, firehose_margin
                )
                if mirror_margin <= 0.0 or firehose_margin <= 0.0:
                    hard_bound_violation_cells += 1
            cell_count += cells
            meshblocks += 1
            total_values += values
            dimensions.add(shape)
            meshblock_metadata.append({
                "logical_location": list(logical),
                "active_indices": list(indices),
                "geometry": list(geometry),
            })
        if stream.tell() != file_size or meshblocks == 0:
            raise ValueError(f"Athena snapshot has no complete meshblock payload: {path}")
        if hard_bound_violation_cells:
            raise ValueError(
                f"Athena snapshot independently violates CGL hard bounds: {path}"
            )
    return {
        "format": "Athena binary output version=1.1",
        "time": snapshot_time,
        "cycle": cycle,
        "number_variables": number_variables,
        "variables": fields[1:],
        "location_size_bytes": location_size,
        "variable_size_bytes": variable_size,
        "parameter_header_size_bytes": header_size,
        "parameter_contract_sha256": parameter_contract_sha256,
        "mesh_nghost": nghost,
        "meshblock_count": meshblocks,
        "meshblock_dimensions": [list(value) for value in sorted(dimensions)],
        "logical_locations": [list(value) for value in sorted(logical_locations)],
        "meshblocks": sorted(
            meshblock_metadata, key=lambda value: value["logical_location"]
        ),
        "positive_variable_minima": positive_minima,
        "variable_sums": variable_sums,
        "cell_count": cell_count,
        "hard_bound_violation_cells": hard_bound_violation_cells,
        "minimum_mirror_hard_margin": minimum_mirror_hard_margin,
        "minimum_firehose_hard_margin": minimum_firehose_hard_margin,
        "total_variable_values": total_values,
    }


def restart_parameter_dump(path: Path) -> tuple[str, int]:
    """Read restart parameter text and return its exact binary payload offset."""

    marker = b"<par_end>\n"
    header = b""
    with path.open("rb") as stream:
        while marker not in header and len(header) < MAX_RESTART_PARAMETER_DUMP_BYTES:
            block = stream.read(
                min(4096, MAX_RESTART_PARAMETER_DUMP_BYTES - len(header))
            )
            if not block:
                break
            header += block
    end = header.find(marker)
    if end < 0:
        raise ValueError(f"restart lacks loadable <par_end> terminator: {path}")
    try:
        text = header[:end].decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"restart parameter dump is not UTF-8: {path}") from error
    return text, end + len(marker)


def restart_time_marker_text(path: Path) -> str:
    """Read the unique physical-time marker from a restart parameter dump."""

    text, _ = restart_parameter_dump(path)
    block = ""
    markers = []
    for original in text.splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            continue
        if block == "time" and "=" in line:
            key, value = line.split("=", 1)
            if key.strip() == "restart_time":
                markers.append(value.strip())
    if len(markers) != 1:
        raise ValueError(f"restart must contain one time/restart_time marker: {path}")
    return markers[0]


def qualified_restart_abi(executable: object) -> tuple[dict[str, object],
                                                       str, str]:
    """Return the reviewed binary-restart ABI for one exact executable."""

    record = require_exact_keys(
        executable, {"path", "sha256", "size_bytes", "revision"},
        "qualification executable provenance",
    )
    revision = require_git_revision(record["revision"], "qualification executable revision")
    digest = require_sha256(record["sha256"], "qualification executable sha256")
    abi = QUALIFIED_RESTART_BINARY_ABIS.get((revision, digest))
    if abi is None:
        raise ValueError("qualification executable has no qualified binary-restart ABI")
    return abi, revision, digest


def structural_restart_profile(path: Path, executable: object,
                               intent: dict[str, object]) -> dict[str, object]:
    """Authenticate restart completeness, configuration, ABI, and physical time."""

    abi, revision, digest = qualified_restart_abi(executable)
    case_id = str(intent["case_id"])
    policy = CASE_POLICIES[case_id]
    if path.stat().st_size < int(policy["minimum_restart_bytes"]):
        raise ValueError(f"restart is smaller than the reviewed completeness floor: {path}")
    parameter_text, payload_offset = restart_parameter_dump(path)
    parameters = parse_athinput_text(parameter_text, "restart parameter dump")
    expected = expected_case_configuration(
        case_id,
        basename=str(intent["run_basename"]),
        tlim=str(intent["target_time"]),
    )
    mismatches = {
        key: {"expected": expected_value, "observed": parameters.get(key)}
        for key, expected_value in expected.items()
        if parameters.get(key) != expected_value
    }
    if mismatches:
        raise ValueError(f"restart parameter dump violates the reviewed contract: {mismatches}")
    marker_text = parameters.get("time/restart_time", "")
    try:
        marker_time = float(marker_text)
    except ValueError as error:
        raise ValueError(f"restart time marker is not numeric: {path}") from error
    if not math.isfinite(marker_time):
        raise ValueError(f"restart time marker is not finite: {path}")
    value_format = str(abi["mesh_time_format"])
    value_size = struct.calcsize(value_format)
    binary_offset = payload_offset + int(abi["mesh_time_offset_after_parameter_dump"])
    with path.open("rb") as stream:
        stream.seek(binary_offset)
        encoded = stream.read(value_size)
    if len(encoded) != value_size:
        raise ValueError(f"restart binary header is truncated: {path}")
    binary_time = float(struct.unpack(value_format, encoded)[0])
    if not math.isfinite(binary_time):
        raise ValueError(f"restart binary physical time is not finite: {path}")
    if marker_text == format(binary_time, ".17g"):
        marker_mode = "full_precision"
    elif marker_text == format(binary_time, ".6g"):
        marker_mode = "legacy_default_precision"
    else:
        raise ValueError(f"restart marker does not authenticate binary physical time: {path}")
    if marker_mode not in abi["allowed_marker_modes"]:
        raise ValueError("restart marker precision is not qualified for this executable")
    return {
        "executable_revision": revision,
        "executable_sha256": digest,
        "parameter_dump_size_bytes": payload_offset,
        "binary_time_offset_bytes": binary_offset,
        "binary_time_format": value_format,
        "marker_mode": marker_mode,
        "time": binary_time,
        "parameter_contract_sha256": stable_json_sha256(expected),
    }


def authenticated_product(path: Path, qualification_root: Path, *, kind: str,
                          executable: object | None = None,
                          intent: dict[str, object] | None = None,
                          ) -> tuple[dict[str, object], float]:
    """Hash and structurally inspect one unchanged rank-local output product."""

    path = require_owned_regular_file(path, qualification_root, "pilot output product")
    with open_regular_fd(
        path, "pilot output product", single_link=True
    ) as descriptor:
        before = os.fstat(descriptor)
        descriptor_path = Path(f"/proc/self/fd/{descriptor}")
        if kind == "snapshot":
            first_profile = structural_snapshot_profile(descriptor_path, intent)
        elif kind == "restart":
            if intent is None:
                raise ValueError("restart inspection lacks execution intent")
            first_profile = structural_restart_profile(descriptor_path, executable, intent)
        else:
            raise ValueError(f"unsupported pilot output product kind: {kind}")
        digest = hashlib.sha256()
        with os.fdopen(os.dup(descriptor), "rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        if kind == "snapshot":
            second_profile = structural_snapshot_profile(descriptor_path, intent)
        else:
            second_profile = structural_restart_profile(
                descriptor_path, executable, intent
            )
        if first_profile != second_profile:
            raise ValueError(f"pilot output product changed during inspection: {path}")
    return {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "size_bytes": before.st_size,
        "inspection": first_profile,
    }, float(first_profile["time"])


def inspect_rank_products(directory: Path, expected_ranks: int,
                          qualification_root: Path, *, kind: str,
                          executable: object | None = None,
                          intent: dict[str, object] | None = None,
                          ) -> tuple[list[dict[str, object]], list[float]]:
    """Require a complete identical rank-local product inventory."""

    require_no_symlink_chain(directory, qualification_root, "rank product directory")
    require_directory(directory, "rank product directory")
    expected = [f"rank_{rank:08d}" for rank in range(expected_ranks)]
    rank_directories = sorted(item.name for item in directory.iterdir())
    if rank_directories != expected:
        raise ValueError("rank-local output directories are incomplete or ambiguous")
    inventories: list[list[str]] = []
    for name in expected:
        rank_dir = require_no_symlink_chain(
            directory / name, qualification_root, "rank output directory"
        )
        require_directory(rank_dir, "rank output directory")
        files = sorted(item.name for item in rank_dir.iterdir())
        if not files:
            raise ValueError("rank-local output directory is empty")
        for item in files:
            require_owned_regular_file(
                rank_dir / item, qualification_root, "rank-local output product"
            )
        inventories.append(files)
    if any(inventory != inventories[0] for inventory in inventories[1:]):
        raise ValueError("rank-local output inventories differ")
    records = []
    times = []
    for product_name in inventories[0]:
        product_records = []
        product_times = []
        for rank_name in expected:
            record, product_time = authenticated_product(
                directory / rank_name / product_name,
                qualification_root,
                kind=kind,
                executable=executable,
                intent=intent,
            )
            product_records.append(record)
            product_times.append(product_time)
        if any(
            abs(value - product_times[0]) > ENDPOINT_TOLERANCE
            for value in product_times[1:]
        ):
            raise ValueError(f"rank-local {kind} times disagree")
        if kind == "snapshot":
            policy = CASE_POLICIES[str(intent["case_id"])]
            expected_dimensions = [list(policy["meshblock_shape"])]
            expected_nghost = int(COMMON_INPUT_CONTRACT["mesh/nghost"])
            expected_indices = [
                value
                for size in policy["meshblock_shape"]
                for value in (expected_nghost, expected_nghost + int(size) - 1)
            ]
            expected_parameter_contract = stable_json_sha256(
                expected_case_configuration(
                    str(intent["case_id"]),
                    basename=str(intent["run_basename"]),
                    tlim=str(intent["target_time"]),
                )
            )
            expected_locations = {
                (lx1, lx2, lx3, 0)
                for lx1 in range(policy["mesh_shape"][0] // policy["meshblock_shape"][0])
                for lx2 in range(policy["mesh_shape"][1] // policy["meshblock_shape"][1])
                for lx3 in range(policy["mesh_shape"][2] // policy["meshblock_shape"][2])
            }
            expected_geometry = {}
            block_counts = [
                policy["mesh_shape"][axis] // policy["meshblock_shape"][axis]
                for axis in range(3)
            ]
            for logical in expected_locations:
                geometry = []
                for axis in range(3):
                    lower, upper = policy["mesh_bounds"][axis]
                    width = (upper - lower) / block_counts[axis]
                    geometry.extend([
                        lower + logical[axis] * width,
                        lower + (logical[axis] + 1) * width,
                    ])
                expected_geometry[logical] = geometry
            observed_locations = []
            observed_cycles = set()
            for record in product_records:
                inspection = record["inspection"]
                if (
                    inspection["variables"] != list(policy["required_snapshot_variables"])
                    or inspection["meshblock_dimensions"] != expected_dimensions
                    or inspection["mesh_nghost"] != expected_nghost
                    or inspection["location_size_bytes"]
                    != policy["snapshot_location_size_bytes"]
                    or inspection["variable_size_bytes"]
                    != policy["snapshot_variable_size_bytes"]
                    or inspection["parameter_contract_sha256"]
                    != expected_parameter_contract
                ):
                    raise ValueError("snapshot metadata differs from the reviewed case")
                required_per_rank = policy.get("required_meshblocks_per_rank")
                if (
                    required_per_rank is not None
                    and inspection["meshblock_count"] != required_per_rank
                ):
                    raise ValueError(
                        "R17 snapshot does not retain exactly 27 meshblocks per rank"
                    )
                observed_cycles.add(int(inspection["cycle"]))
                observed_locations.extend(
                    tuple(location) for location in inspection["logical_locations"]
                )
                for meshblock in inspection["meshblocks"]:
                    logical = tuple(meshblock["logical_location"])
                    geometry = expected_geometry.get(logical)
                    if (
                        meshblock["active_indices"] != expected_indices
                        or geometry is None
                        or any(
                            abs(float(observed) - expected) > ENDPOINT_TOLERANCE
                            for observed, expected in zip(
                                meshblock["geometry"], geometry
                            )
                        )
                    ):
                        raise ValueError(
                            "snapshot meshblock metadata or geometry differs "
                            "from the reviewed case"
                        )
            if (
                len(observed_cycles) != 1
                or len(observed_locations) != len(set(observed_locations))
                or set(observed_locations) != expected_locations
            ):
                raise ValueError("snapshot logical meshblock inventory is incomplete or ambiguous")
        records.append({"name": product_name, "rank_files": product_records})
        times.append(product_times[0])
    return records, times


def validate_restart_smoke(packet: dict[str, object], qualification_root: Path,
                           terminal_restart: dict[str, object]) -> dict[str, object]:
    """Require the exact terminal restart to have completed a bound load smoke."""

    intent = packet["execution_intent"]
    smoke_dir = require_no_symlink_chain(
        Path(str(intent["paths"]["restart_smoke_dir"])),
        qualification_root,
        "restart smoke directory",
    )
    require_directory(smoke_dir, "restart smoke directory")
    result_path = require_owned_regular_file(
        Path(str(intent["paths"]["restart_smoke_result"])),
        qualification_root,
        "restart smoke result",
    )
    log_path = require_owned_regular_file(
        Path(str(intent["paths"]["restart_smoke_log"])),
        qualification_root,
        "restart smoke log",
    )
    environment_path = require_owned_regular_file(
        Path(str(intent["paths"]["environment_log"])),
        qualification_root,
        "run environment log",
    )
    binding_path = require_owned_regular_file(
        Path(str(intent["paths"]["job_binding"])),
        qualification_root,
        "restart smoke job binding",
    )
    binding = load_json(binding_path, "restart smoke job binding")
    result = require_exact_keys(
        load_json(result_path, "restart smoke result"),
        {
            "schema_version", "record_type", "execution_contract_sha256",
            "executable_sha256", "job_id", "restart_path", "restart_sha256",
            "restart_size_bytes", "exit_code", "completed_utc",
        },
        "restart smoke result",
    )
    try:
        completed = datetime.fromisoformat(str(result["completed_utc"]).replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError("restart smoke result has invalid completion time") from error
    if completed.tzinfo is None:
        raise ValueError("restart smoke result lacks an explicit UTC offset")
    try:
        environment_text = read_regular_bytes(
            environment_path, "run environment log", single_link=True
        ).decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("run environment log is not UTF-8") from error
    required_environment_keys = {
        "started_utc", "finished_utc", "slurm_job_id",
        "qualification_execution_contract_sha256", "prepared_manifest",
        "nodes", "ranks",
    }
    environment_values = {}
    for line in environment_text.splitlines():
        key, separator, value = line.partition("=")
        if key in required_environment_keys:
            if not separator or not value or key in environment_values:
                raise ValueError("run environment log is incomplete or ambiguous")
            environment_values[key] = value
    if set(environment_values) != required_environment_keys:
        raise ValueError("run environment log is incomplete or ambiguous")
    try:
        started = datetime.fromisoformat(
            environment_values["started_utc"].replace("Z", "+00:00")
        )
        finished = datetime.fromisoformat(
            environment_values["finished_utc"].replace("Z", "+00:00")
        )
    except ValueError as error:
        raise ValueError("run environment log has invalid timestamps") from error
    if (
        started.tzinfo is None
        or finished.tzinfo is None
        or not started <= completed <= finished
    ):
        raise ValueError("run environment and restart smoke timestamps are inconsistent")
    if (
        result["schema_version"] != 1
        or result["record_type"] != "cgl_lf_stage_i_restart_load_smoke"
        or result["execution_contract_sha256"] != intent["execution_contract_sha256"]
        or result["executable_sha256"]
        != terminal_restart["inspection"]["executable_sha256"]
        or binding.get("schema_version") != 4
        or binding.get("execution_intent_sha256") != packet["execution_intent_sha256"]
        or binding.get("execution_contract_sha256") != intent["execution_contract_sha256"]
        or environment_values["slurm_job_id"] != binding.get("job_id")
        or environment_values["qualification_execution_contract_sha256"]
        != intent["execution_contract_sha256"]
        or environment_values["prepared_manifest"] != intent["paths"]["prepared_manifest"]
        or environment_values["nodes"] != str(intent["allocation"]["nodes"])
        or environment_values["ranks"]
        != str(
            int(intent["allocation"]["nodes"])
            * int(intent["allocation"]["ranks_per_node"])
        )
        or result["job_id"] != binding.get("job_id")
        or result["restart_path"] != terminal_restart["path"]
        or result["restart_sha256"] != terminal_restart["sha256"]
        or result["restart_size_bytes"] != terminal_restart["size_bytes"]
        or result["exit_code"] != 0
    ):
        raise ValueError("restart smoke result differs from the authenticated terminal restart")
    executable_sha256 = require_sha256(
        result["executable_sha256"], "restart smoke executable sha256"
    )
    return {
        "result": retained_file(result_path, "restart smoke result"),
        "log": retained_file(log_path, "restart smoke log"),
        "environment": retained_file(environment_path, "run environment log"),
        "restart_smoke_dir": str(smoke_dir),
        "job_id": result["job_id"],
        "restart_sha256": result["restart_sha256"],
        "executable_sha256": executable_sha256,
        "exit_code": 0,
        "started_utc": environment_values["started_utc"],
        "completed_utc": result["completed_utc"],
        "finished_utc": environment_values["finished_utc"],
    }


def require_columns(history: dict[str, list[float]], names: tuple[str, ...],
                    label: str) -> None:
    """Require all named history columns."""

    missing = [name for name in names if name not in history]
    if missing:
        raise ValueError(f"{label} lacks required columns: {missing}")


def relative_difference(first: float, second: float) -> float:
    """Return one stable normalized absolute difference."""

    return abs(first - second) / max(abs(first), abs(second), 1.0)


def compare_scientific_signatures(baseline: object, candidate: object, *,
                                  path: str = "signature",
                                  differences: list[tuple[float, float]] | None = None,
                                  ) -> list[tuple[float, float]]:
    """Require recursively identical structure and tolerance-bounded numeric values."""

    if differences is None:
        differences = []
    if isinstance(baseline, dict):
        if not isinstance(candidate, dict) or set(candidate) != set(baseline):
            raise ValueError(f"cross-profile scientific agreement structure differs at {path}")
        for key in sorted(baseline):
            compare_scientific_signatures(
                baseline[key],
                candidate[key],
                path=f"{path}.{key}",
                differences=differences,
            )
        return differences
    if isinstance(baseline, list):
        if not isinstance(candidate, list) or len(candidate) != len(baseline):
            raise ValueError(f"cross-profile scientific agreement structure differs at {path}")
        for index, (first, second) in enumerate(zip(baseline, candidate)):
            compare_scientific_signatures(
                first,
                second,
                path=f"{path}[{index}]",
                differences=differences,
            )
        return differences
    if type(baseline) is not type(candidate):
        raise ValueError(f"cross-profile scientific agreement type differs at {path}")
    if isinstance(baseline, bool):
        if baseline != candidate:
            raise ValueError(f"cross-profile scientific agreement differs at {path}")
        return differences
    if isinstance(baseline, int):
        if baseline != candidate:
            raise ValueError(f"cross-profile scientific agreement differs at {path}")
        return differences
    if isinstance(baseline, float):
        if not isinstance(candidate, float) or not math.isfinite(candidate):
            raise ValueError(f"cross-profile scientific agreement is invalid at {path}")
        absolute = abs(baseline - candidate)
        relative = relative_difference(baseline, candidate)
        differences.append((absolute, relative))
        if not math.isclose(
            baseline,
            candidate,
            rel_tol=CROSS_PROFILE_RELATIVE_TOLERANCE,
            abs_tol=CROSS_PROFILE_ABSOLUTE_TOLERANCE,
        ):
            raise ValueError(f"cross-profile scientific agreement failed at {path}")
        return differences
    if baseline != candidate:
        raise ValueError(f"cross-profile scientific agreement differs at {path}")
    return differences


def cross_profile_scientific_agreement(
    observed: dict[int, dict[str, object]], baseline_nodes: int = 1,
) -> dict[str, object]:
    """Require every candidate node profile to reproduce the reviewed baseline."""

    if baseline_nodes not in observed:
        raise ValueError("cross-profile scientific agreement lacks reviewed baseline")
    baseline = observed[baseline_nodes].get("scientific_agreement_signature")
    if not isinstance(baseline, dict):
        raise ValueError("baseline scientific agreement signature is invalid")
    comparisons = []
    for nodes in sorted(observed):
        if nodes == baseline_nodes:
            continue
        candidate = observed[nodes].get("scientific_agreement_signature")
        if not isinstance(candidate, dict):
            raise ValueError("candidate scientific agreement signature is invalid")
        differences = compare_scientific_signatures(baseline, candidate)
        comparisons.append({
            "nodes": nodes,
            "numeric_values_compared": len(differences),
            "maximum_absolute_difference": max(
                (value[0] for value in differences), default=0.0
            ),
            "maximum_relative_difference": max(
                (value[1] for value in differences), default=0.0
            ),
            "agrees_with_one_node": True,
        })
    return {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_cross_profile_scientific_agreement",
        "baseline_nodes": baseline_nodes,
        "relative_tolerance": CROSS_PROFILE_RELATIVE_TOLERANCE,
        "absolute_tolerance": CROSS_PROFILE_ABSOLUTE_TOLERANCE,
        "comparisons": comparisons,
        "all_profiles_agree": True,
    }


def rank_file_inventory(records: list[dict[str, object]], root: Path,
                        expected_ranks: int, label: str
                        ) -> tuple[list[dict[str, str]], str]:
    """Return one recost-compatible exact rank-local file inventory."""

    if len(records) != expected_ranks:
        raise ValueError(f"{label} does not contain exactly one file per rank")
    inventory = []
    ranks = set()
    for record in records:
        path = Path(str(record["path"])).resolve(strict=True)
        require_beneath(path, root, label)
        relative = path.relative_to(root)
        rank_parts = [
            part for part in relative.parts if re.fullmatch(r"rank_[0-9]{8}", part)
        ]
        if len(rank_parts) != 1:
            raise ValueError(f"{label} lacks one rank-local directory")
        ranks.add(rank_parts[0])
        inventory.append({
            "path": relative.as_posix(),
            "sha256": require_sha256(record["sha256"], f"{label} sha256"),
        })
    inventory.sort(key=lambda item: item["path"])
    expected = {f"rank_{rank:08d}" for rank in range(expected_ranks)}
    if ranks != expected or len({item["path"] for item in inventory}) != expected_ranks:
        raise ValueError(f"{label} rank inventory is incomplete or ambiguous")
    return inventory, retained_inventory_sha256(inventory)


def validate_r17_decomposition_evidence(
    value: object,
    terminal_output_inventory_sha256: str,
) -> dict[str, object]:
    """Authenticate the exact complete R17 logical-block/rank decomposition."""

    evidence = require_exact_keys(
        value,
        {
            "schema_version", "record_type", "resolution", "mesh_shape",
            "meshblock_shape", "logical_meshblock_grid", "logical_meshblocks",
            "ranks", "meshblocks_per_rank", "complete_block_rank_inventory",
            "complete_block_rank_inventory_sha256",
            "terminal_rank_local_output_inventory_sha256", "checks",
        },
        "R17 decomposition evidence",
    )
    inventory = evidence["complete_block_rank_inventory"]
    if not isinstance(inventory, list) or len(inventory) != R17_OPERATIONAL_RANKS:
        raise ValueError("R17 decomposition rank inventory is incomplete")
    expected_locations = {
        (lx1, lx2, lx3, 0)
        for lx1 in range(12)
        for lx2 in range(12)
        for lx3 in range(12)
    }
    observed_locations = []
    for rank, item in enumerate(inventory):
        record = require_exact_keys(
            item, {"rank", "rank_name", "logical_meshblocks"},
            "R17 decomposition rank record",
        )
        locations = record["logical_meshblocks"]
        if (
            record["rank"] != rank
            or record["rank_name"] != f"rank_{rank:08d}"
            or not isinstance(locations, list)
            or len(locations) != 27
            or locations != sorted(locations)
            or len({tuple(location) for location in locations}) != 27
            or any(
                not isinstance(location, list)
                or len(location) != 4
                or any(isinstance(index, bool) or not isinstance(index, int)
                       for index in location)
                or tuple(location) not in expected_locations
                for location in locations
            )
        ):
            raise ValueError("R17 decomposition does not retain 27 unique blocks per rank")
        observed_locations.extend(tuple(location) for location in locations)
    checks = require_exact_keys(
        evidence["checks"],
        {
            "exact_resolution", "exact_rank_count",
            "exact_meshblocks_per_rank", "complete_unique_logical_inventory",
        },
        "R17 decomposition checks",
    )
    if (
        evidence["schema_version"] != 1
        or evidence["record_type"] != "stage-i-r17-decomposition-evidence"
        or evidence["resolution"] != "384x384x768"
        or evidence["mesh_shape"] != [384, 384, 768]
        or evidence["meshblock_shape"] != [32, 32, 64]
        or evidence["logical_meshblock_grid"] != [12, 12, 12]
        or evidence["logical_meshblocks"] != 1728
        or evidence["ranks"] != R17_OPERATIONAL_RANKS
        or evidence["meshblocks_per_rank"] != 27
        or evidence["complete_block_rank_inventory_sha256"]
        != stable_json_sha256(inventory)
        or evidence["terminal_rank_local_output_inventory_sha256"]
        != require_sha256(
            terminal_output_inventory_sha256,
            "R17 terminal rank-local output inventory digest",
        )
        or len(observed_locations) != 1728
        or len(set(observed_locations)) != 1728
        or set(observed_locations) != expected_locations
        or checks != {
            "exact_resolution": True,
            "exact_rank_count": True,
            "exact_meshblocks_per_rank": True,
            "complete_unique_logical_inventory": True,
        }
    ):
        raise ValueError("R17 decomposition evidence differs from the reviewed contract")
    return evidence


def build_r17_decomposition_evidence(
    terminal_snapshot_files: list[dict[str, object]],
    terminal_output_inventory_sha256: str,
) -> dict[str, object]:
    """Build a deterministic complete R17 logical-block/rank inventory."""

    policy = CASE_POLICIES[R17_CASE_ID]
    logical_grid = [
        int(policy["mesh_shape"][axis]) // int(policy["meshblock_shape"][axis])
        for axis in range(3)
    ]
    if (
        policy["resolution"] != "384x384x768"
        or logical_grid != [12, 12, 12]
        or len(terminal_snapshot_files) != R17_OPERATIONAL_RANKS
    ):
        raise ValueError("R17 snapshot decomposition differs from the reviewed contract")
    inventory = []
    for rank, record in enumerate(terminal_snapshot_files):
        inspection = record["inspection"]
        locations = sorted(
            [list(location) for location in inspection["logical_locations"]]
        )
        inventory.append({
            "rank": rank,
            "rank_name": f"rank_{rank:08d}",
            "logical_meshblocks": locations,
        })
    evidence = {
        "schema_version": 1,
        "record_type": "stage-i-r17-decomposition-evidence",
        "resolution": "384x384x768",
        "mesh_shape": [384, 384, 768],
        "meshblock_shape": [32, 32, 64],
        "logical_meshblock_grid": [12, 12, 12],
        "logical_meshblocks": 1728,
        "ranks": R17_OPERATIONAL_RANKS,
        "meshblocks_per_rank": 27,
        "complete_block_rank_inventory": inventory,
        "complete_block_rank_inventory_sha256": stable_json_sha256(inventory),
        "terminal_rank_local_output_inventory_sha256": (
            terminal_output_inventory_sha256
        ),
        "checks": {
            "exact_resolution": True,
            "exact_rank_count": True,
            "exact_meshblocks_per_rank": True,
            "complete_unique_logical_inventory": True,
        },
    }
    return validate_r17_decomposition_evidence(
        evidence, terminal_output_inventory_sha256
    )


def inspect_pilot_output(packet: dict[str, object], qualification_root: Path,
                         provenance: dict[str, object]) -> dict[str, object]:
    """Reproduce a case-aware exact-endpoint pilot scientific inspection."""

    intent = packet["execution_intent"]
    if stable_json_sha256(provenance) != intent["provenance_sha256"]:
        raise ValueError("pilot scientific inspection provenance differs from intent")
    executable = require_exact_keys(
        provenance.get("executable"),
        {"path", "sha256", "size_bytes", "revision"},
        "qualification executable provenance",
    )
    output = require_no_symlink_chain(
        Path(str(intent["paths"]["output_dir"])), qualification_root,
        "pilot output directory",
    )
    require_directory(output, "pilot output directory")
    entries = {item.name for item in output.iterdir()}
    histories = sorted(output.glob("*.hst"))
    mhd_paths = [path for path in histories if path.name.endswith(".mhd.hst")]
    user_paths = [path for path in histories if path.name.endswith(".user.hst")]
    expected_entries = {"bin", "rst", *(path.name for path in histories)}
    if entries != expected_entries or len(mhd_paths) != 1 or len(user_paths) != 1:
        raise ValueError("pilot output top-level inventory is incomplete or ambiguous")
    mhd, mhd_record = authenticated_history(
        mhd_paths[0], qualification_root, "MHD history"
    )
    user, user_record = authenticated_history(
        user_paths[0], qualification_root, "user history"
    )
    require_columns(
        mhd,
        (
            "time", "mass", "tot-E", "lf_nstage", *STRICT_LF_FAILURE_COLUMNS,
            "lf_qface", "lf_qprcap", "lf_qpecap", "lf_qprwrk", "lf_qpewrk",
            "lf_hwproj", "lf_cpwrk", "lf_cawrk",
        ),
        "MHD history",
    )
    require_columns(
        user,
        ("time", "mass", "hard_vol", "force_pwr", "force_work", "max_ndiv"),
        "user history",
    )
    target = float(intent["target_time"])
    cadence = float(CASE_POLICIES[str(intent["case_id"])]["history_dt"])
    require_history_schedule(mhd, target, cadence, "MHD history")
    require_history_schedule(user, target, cadence, "user history")
    if (
        abs(mhd["time"][-1] - target) > ENDPOINT_TOLERANCE
        or abs(user["time"][-1] - target) > ENDPOINT_TOLERANCE
        or len(mhd["time"]) != len(user["time"])
        or any(
            abs(first - second) > ENDPOINT_TOLERANCE
            for first, second in zip(mhd["time"], user["time"])
        )
    ):
        raise ValueError("pilot histories do not share the exact target endpoint")
    if any(
        value != 0.0 for name in STRICT_LF_FAILURE_COLUMNS for value in mhd[name]
    ):
        raise ValueError("pilot strict LF failure counters are nonzero")
    if any(value != 0.0 for value in user["hard_vol"]):
        raise ValueError("pilot hard-bound volume is nonzero")
    if (
        any(value < 0.0 for value in user["max_ndiv"])
        or max(user["max_ndiv"]) >= MAX_NORMALIZED_CT_DIVB
    ):
        raise ValueError("pilot normalized CT-divB gate failed")
    if (
        any(
            relative_difference(mhd["mass"][0], value) > MASS_TOLERANCE
            for value in mhd["mass"]
        )
        or any(
            relative_difference(user["mass"][0], value) > MASS_TOLERANCE
            for value in user["mass"]
        )
        or any(
            relative_difference(first, second) > MASS_TOLERANCE
            for first, second in zip(mhd["mass"], user["mass"])
        )
    ):
        raise ValueError("pilot mass conservation or history agreement failed")
    if mhd["lf_nstage"][-1] <= mhd["lf_nstage"][0] or mhd["lf_qface"][-1] <= 0.0:
        raise ValueError("pilot did not exercise the LF update")
    require_monotonic_columns(
        mhd,
        ("lf_nstage", "lf_qface", "lf_qprcap", "lf_qpecap", "lf_hwproj"),
        "MHD history",
    )
    count_columns = ("lf_nstage", "lf_qface", "lf_qprcap", "lf_qpecap", "lf_hwproj")
    if any(
        value < 0.0 or not value.is_integer()
        for name in count_columns for value in mhd[name]
    ):
        raise ValueError("pilot LF diagnostic count history is not integral")
    if any(
        later <= earlier
        for name in ("lf_nstage", "lf_qface")
        for earlier, later in zip(mhd[name], mhd[name][1:])
    ):
        raise ValueError("pilot LF update counters have an inactive history interval")
    if any(
        value < 0.0 or value > qface
        for name in ("lf_qprcap", "lf_qpecap")
        for value, qface in zip(mhd[name], mhd["lf_qface"])
    ):
        raise ValueError("pilot LF cap counters are inconsistent")
    for name in ("lf_qprcap", "lf_qpecap"):
        for previous_cap, cap, previous_qface, qface in zip(
            mhd[name], mhd[name][1:], mhd["lf_qface"], mhd["lf_qface"][1:]
        ):
            if cap - previous_cap > qface - previous_qface:
                raise ValueError(
                    "pilot LF cap increment exceeds interval qface evaluations"
                )
    policy = str(intent["scientific_policy"])
    forcing_closure = None
    if policy == "passive_delta":
        if any(value != 0.0 for name in ("lf_cpwrk", "lf_cawrk") for value in mhd[name]):
            raise ValueError("passive-Delta pilot has active pressure work")
    else:
        energy_delta = mhd["tot-E"][-1] - mhd["tot-E"][0]
        forcing_delta = user["force_work"][-1] - user["force_work"][0]
        pressure_work_delta = sum(
            abs(mhd[name][-1] - mhd[name][0]) for name in ("lf_cpwrk", "lf_cawrk")
        )
        lf_work_delta = sum(
            abs(mhd[name][-1] - mhd[name][0]) for name in ("lf_qprwrk", "lf_qpewrk")
        )
        final_pressure_work = sum(
            abs(mhd[name][-1] - mhd[name][-2]) for name in ("lf_cpwrk", "lf_cawrk")
        )
        final_lf_work = sum(
            abs(mhd[name][-1] - mhd[name][-2]) for name in ("lf_qprwrk", "lf_qpewrk")
        )
        minimum_forcing_work = target * MINIMUM_FORCING_WORK_RATE
        minimum_pressure_work = target * MINIMUM_PRESSURE_WORK_RATE
        minimum_lf_work = target * MINIMUM_LF_WORK_RATE
        if (
            abs(forcing_delta) < minimum_forcing_work
            or max(abs(value) for value in user["force_pwr"])
            < MINIMUM_FORCING_WORK_RATE
            or pressure_work_delta < minimum_pressure_work
            or lf_work_delta < minimum_lf_work
            or abs(user["force_work"][-1] - user["force_work"][-2])
            < cadence * MINIMUM_FORCING_WORK_RATE
            or final_pressure_work
            < pressure_work_delta * MINIMUM_FINAL_INTERVAL_WORK_FRACTION
            or final_lf_work
            < lf_work_delta * MINIMUM_FINAL_INTERVAL_WORK_FRACTION
        ):
            raise ValueError("active-CGL pilot lacks nontrivial forcing or pressure work")
        forcing_closure = relative_difference(energy_delta, forcing_delta)
        if forcing_closure > ACTIVE_FORCING_CLOSURE_TOLERANCE:
            raise ValueError("active-CGL pilot forcing-energy closure failed")
    if policy == "active_finite_limiter":
        if any(value != 0.0 for value in mhd["lf_hwproj"]):
            raise ValueError("finite-limiter pilot has hardwall projections")
    else:
        if any(value < 0.0 for value in mhd["lf_hwproj"]) or any(
            later < earlier
            for earlier, later in zip(mhd["lf_hwproj"], mhd["lf_hwproj"][1:])
        ):
            raise ValueError("hardwall projection counter is invalid")
    expected_ranks = (
        int(intent["allocation"]["nodes"])
        * int(intent["allocation"]["ranks_per_node"])
    )
    snapshots, snapshot_times = inspect_rank_products(
        output / "bin", expected_ranks, qualification_root, kind="snapshot",
        intent=intent,
    )
    restarts, restart_times = inspect_rank_products(
        output / "rst", expected_ranks, qualification_root, kind="restart",
        executable=executable,
        intent=intent,
    )
    expected_snapshot_times = expected_history_times(
        target, float(CASE_POLICIES[str(intent["case_id"])]["snapshot_dt"])
    )
    if (
        len(snapshot_times) != len(expected_snapshot_times)
        or any(
            abs(observed - expected) > ENDPOINT_TOLERANCE
            for observed, expected in zip(snapshot_times, expected_snapshot_times)
        )
    ):
        raise ValueError("pilot snapshot products violate the exact output cadence")
    terminal = [
        value for value in snapshot_times if abs(value - target) <= ENDPOINT_TOLERANCE
    ]
    if len(terminal) != 1 or any(value > target + ENDPOINT_TOLERANCE for value in snapshot_times):
        raise ValueError("pilot lacks one unique exact terminal snapshot")
    terminal_snapshot_groups = [
        record for record, value in zip(snapshots, snapshot_times)
        if abs(value - target) <= ENDPOINT_TOLERANCE
    ]
    terminal_snapshot_files = terminal_snapshot_groups[0]["rank_files"]
    project_root = Path(str(intent["project_root"])).resolve(strict=True)
    terminal_outputs, terminal_output_sha256 = rank_file_inventory(
        terminal_snapshot_files,
        project_root,
        expected_ranks,
        "terminal rank-local outputs",
    )
    r17_decomposition = None
    if intent["case_id"] == R17_CASE_ID:
        r17_decomposition = build_r17_decomposition_evidence(
            terminal_snapshot_files, terminal_output_sha256
        )
    terminal_snapshot_cell_count = sum(
        int(record["inspection"]["cell_count"]) for record in terminal_snapshot_files
    )
    terminal_snapshot_variable_sums = {
        name: math.fsum(
            float(record["inspection"]["variable_sums"][name])
            for record in terminal_snapshot_files
        )
        for name in REQUIRED_SNAPSHOT_VARIABLES
    }
    if (
        terminal_snapshot_cell_count <= 0
        or any(
            record["inspection"]["hard_bound_violation_cells"] != 0
            for record in terminal_snapshot_files
        )
    ):
        raise ValueError("terminal snapshot independent CGL hard-bound gate failed")
    mhd_agreement_columns = (
        "mass", "tot-E", "lf_nstage", "lf_qface", "lf_qprcap",
        "lf_qpecap", "lf_qprwrk", "lf_qpewrk", "lf_hwproj",
        "lf_cpwrk", "lf_cawrk",
    )
    agreement_signature = {
        "schema_version": 1,
        "mhd_history": {
            name: (
                [int(value) for value in mhd[name]]
                if name in count_columns
                else list(mhd[name])
            )
            for name in mhd_agreement_columns
        },
        "user_history": {
            name: list(user[name])
            for name in ("mass", "hard_vol", "force_work", "max_ndiv")
        },
        "terminal_snapshot": {
            "cell_count": terminal_snapshot_cell_count,
            "variable_sums": terminal_snapshot_variable_sums,
        },
    }
    terminal_restarts = [
        value for value in restart_times if abs(value - target) <= ENDPOINT_TOLERANCE
    ]
    if len(terminal_restarts) != 1 or any(
        value > target + ENDPOINT_TOLERANCE for value in restart_times
    ):
        raise ValueError("pilot lacks one unique exact terminal restart")
    terminal_restart_groups = [
        record for record, value in zip(restarts, restart_times)
        if abs(value - target) <= ENDPOINT_TOLERANCE
    ]
    terminal_restart_files = terminal_restart_groups[0]["rank_files"]
    terminal_restarts_inventory, terminal_restart_sha256 = rank_file_inventory(
        terminal_restart_files,
        project_root,
        expected_ranks,
        "terminal rank-local restarts",
    )
    terminal_restart = terminal_restart_files[0]
    restart_smoke = validate_restart_smoke(
        packet, qualification_root, terminal_restart
    )
    mass_relative_drift_max = max(
        [
            relative_difference(mhd["mass"][0], value) for value in mhd["mass"]
        ]
        + [
            relative_difference(user["mass"][0], value) for value in user["mass"]
        ]
    )
    mhd_user_mass_mismatch_max = max(
        relative_difference(first, second)
        for first, second in zip(mhd["mass"], user["mass"])
    )
    normalized_ct_divb_max = max(user["max_ndiv"])
    operational_only = bool(
        CASE_POLICIES[str(intent["case_id"])]["operational_only"]
    )
    return {
        "schema_version": 6,
        "record_type": "cgl_lf_stage_i_qualification_scientific_evidence",
        "case_id": intent["case_id"],
        "nodes": intent["allocation"]["nodes"],
        "target_time": target,
        "scientific_policy": policy,
        "executable": {
            "revision": executable["revision"],
            "sha256": executable["sha256"],
        },
        "execution_intent_sha256": packet["execution_intent_sha256"],
        "execution_contract_sha256": intent["execution_contract_sha256"],
        "forcing_closure_normalized_residual": forcing_closure,
        "mhd_history": mhd_record,
        "user_history": user_record,
        "snapshots": snapshots,
        "snapshot_times": snapshot_times,
        "restarts": restarts,
        "restart_times": restart_times,
        "restart_load_smoke": restart_smoke,
        "terminal_rank_local_outputs": terminal_outputs,
        "terminal_rank_local_output_inventory_sha256": terminal_output_sha256,
        "r17_decomposition": r17_decomposition,
        "terminal_rank_local_restarts": terminal_restarts_inventory,
        "terminal_rank_local_restart_inventory_sha256": terminal_restart_sha256,
        "physics_measurements": {
            "finite_rank_outputs": expected_ranks,
            "mass_relative_drift_max": format(mass_relative_drift_max, ".17g"),
            "mhd_user_mass_mismatch_max": format(
                mhd_user_mass_mismatch_max, ".17g"
            ),
            "lf_bad_counts_total": 0,
            "normalized_ct_divb_max": format(normalized_ct_divb_max, ".17g"),
            "normalized_ct_divb_threshold": MAX_NORMALIZED_CT_DIVB_TEXT,
            "normalized_ct_divb_below_threshold": True,
        },
        "scientific_agreement_signature": agreement_signature,
        "checks": {
            "exact_endpoint": True,
            "complete_rank_inventory": True,
            "finite_synchronized_histories": True,
            "mass_conserved": True,
            "strict_lf_failure_counters_zero": True,
            "hard_volume_zero": True,
            "snapshot_hard_bounds_independently_verified": True,
            "normalized_ct_divb_below_threshold": True,
            "interval_cap_counts_valid": True,
            "nontrivial_forcing_and_pressure_work": True,
            "case_aware_policy_passed": True,
            "snapshot_cadence_complete": True,
            "terminal_snapshot_unique": True,
            "snapshots_structurally_complete": True,
            "restart_headers_authenticated": True,
            "terminal_restart_unique": True,
            "restart_load_smoke_passed": True,
        },
        "accepted_for_profile_selection": not operational_only,
        "accepted_for_operational_qualification": operational_only,
    }


def validate_scientific_evidence(path: Path, *, packet: dict[str, object],
                                 qualification_root: Path,
                                 provenance: dict[str, object]) -> dict[str, object]:
    """Reproduce and require exact case-aware scientific evidence."""

    retained = load_json(path, "qualification scientific evidence")
    reproduced = inspect_pilot_output(packet, qualification_root, provenance)
    if retained != reproduced:
        raise ValueError("qualification scientific evidence differs from live products")
    return retained


def audit_retained_wave(wave: dict[str, object],
                        qualification_root: Path,
                        scheduler_runner=subprocess.run,
                        batch_script_runner=None,
                        now: datetime | None = None,
                        ) -> dict[tuple[str, int], dict[str, object]]:
    """Reproduce the complete scheduler, launch, and scientific wave audit."""

    validate_prepared_wave(wave)
    validate_submission_journal(wave, qualification_root)
    if batch_script_runner is None:
        batch_script_runner = scheduler_runner
    observed: dict[tuple[str, int], dict[str, object]] = {}
    seen_job_ids = set()
    previous_job_ids: list[str] = []
    previous_wave_end: datetime | None = None
    for wave_record in wave["waves"]:
        current_job_ids = []
        current_starts = []
        current_ends = []
        current_nodes = 0
        for packet in wave_record["packets"]:
            intent = packet["execution_intent"]
            binding = validate_job_binding(
                packet,
                int(wave_record["wave"]),
                previous_job_ids,
                qualification_root,
                current_wave_job_ids=current_job_ids,
                expected_wave_nodes=int(wave_record["total_nodes"]),
                batch_script_runner=batch_script_runner,
            )
            job_id = str(binding["job_id"])
            if job_id in seen_job_ids:
                raise ValueError("qualification scheduler job IDs are ambiguous")
            seen_job_ids.add(job_id)
            scheduler_path = require_owned_regular_file(
                Path(str(intent["paths"]["scheduler_evidence"])),
                qualification_root,
                "scheduler evidence",
            )
            scientific_path = require_owned_regular_file(
                Path(str(intent["paths"]["scientific_evidence"])),
                qualification_root,
                "scientific evidence",
            )
            scheduler = validate_scheduler_evidence(
                scheduler_path,
                packet=packet,
                binding=binding,
                qualification_root=qualification_root,
                runner=scheduler_runner,
            )
            require_recent_scheduler_end(scheduler["end_utc"], now=now)
            output_binding_path = Path(str(intent["paths"]["output_binding"]))
            validate_output_binding(
                output_binding_path,
                packet=packet,
                binding=binding,
                scheduler=scheduler,
                qualification_root=qualification_root,
                provenance=wave["provenance"],
            )
            scientific = validate_scientific_evidence(
                scientific_path,
                packet=packet,
                qualification_root=qualification_root,
                provenance=wave["provenance"],
            )
            start = parse_scheduler_time(str(scheduler["start_utc"]), "scheduler start")
            end = parse_scheduler_time(str(scheduler["end_utc"]), "scheduler end")
            smoke = scientific["restart_load_smoke"]
            smoke_started = parse_scheduler_time(
                str(smoke["started_utc"]), "run environment start"
            )
            smoke_completed = parse_scheduler_time(
                str(smoke["completed_utc"]), "restart smoke completion"
            )
            smoke_finished = parse_scheduler_time(
                str(smoke["finished_utc"]), "run environment finish"
            )
            if not start <= smoke_started <= smoke_completed <= smoke_finished <= end:
                raise ValueError("run environment and restart smoke fall outside Slurm job")
            current_starts.append(start)
            current_ends.append(end)
            current_nodes += int(intent["allocation"]["nodes"])
            current_job_ids.append(job_id)
            observed[(str(intent["case_id"]), int(intent["allocation"]["nodes"]))] = {
                "nodes": intent["allocation"]["nodes"],
                "elapsed_seconds": scheduler["elapsed_seconds"],
                "job_id": job_id,
                "scheduler_end_utc": scheduler["end_utc"],
                "execution_intent_sha256": packet["execution_intent_sha256"],
                "job_binding": retained_file(
                    Path(str(intent["paths"]["job_binding"])), "job binding"
                ),
                "scheduler_evidence": retained_file(
                    scheduler_path, "scheduler evidence"
                ),
                "scientific_evidence": retained_file(
                    scientific_path, "scientific evidence"
                ),
                "output_binding": retained_file(
                    output_binding_path, "pilot output binding"
                ),
                "scientific_agreement_signature": scientific[
                    "scientific_agreement_signature"
                ],
            }
        if current_nodes > int(wave["policy"]["max_concurrent_nodes"]):
            raise ValueError("retained qualification wave exceeds node ceiling")
        if previous_wave_end is not None and min(current_starts) < previous_wave_end:
            raise ValueError("retained qualification waves overlap")
        previous_wave_end = max(current_ends)
        previous_job_ids = current_job_ids
    return observed


def retain_wave_audit(prepared_wave_path: Path,
                      scheduler_runner=subprocess.run,
                      batch_script_runner=None,
                      now: datetime | None = None) -> list[Path]:
    """Collect live Slurm records and reproduce scientific evidence for every pilot."""

    wave = load_json(prepared_wave_path, "prepared qualification wave")
    packets = validate_prepared_wave(wave)
    root = require_project_root(str(wave["project_root"]), allow_local_root=False)
    qualification_root = require_qualification_root(
        root, str(wave["qualification_root"])
    )
    if prepared_wave_path.resolve(strict=True) != (
        qualification_root / "prepared_wave.json"
    ).resolve(strict=True):
        raise ValueError("audit requires the canonical retained prepared wave")
    previous_job_ids: list[str] = []
    if batch_script_runner is None:
        batch_script_runner = scheduler_runner
    validate_submission_journal(wave, qualification_root)
    for wave_record in wave["waves"]:
        current_job_ids = []
        for packet in wave_record["packets"]:
            intent = packet["execution_intent"]
            binding = validate_job_binding(
                packet,
                int(wave_record["wave"]),
                previous_job_ids,
                qualification_root,
                current_wave_job_ids=current_job_ids,
                expected_wave_nodes=int(wave_record["total_nodes"]),
                batch_script_runner=batch_script_runner,
            )
            job_id = str(binding["job_id"])
            raw_path = Path(str(intent["paths"]["scheduler_raw"]))
            write_exclusive(
                raw_path,
                collect_live_scheduler_raw(job_id, scheduler_runner).encode("utf-8"),
                mode=0o440,
            )
            scheduler = build_scheduler_evidence(
                packet, binding, raw_path, qualification_root
            )
            validate_output_binding(
                Path(str(intent["paths"]["output_binding"])),
                packet=packet,
                binding=binding,
                scheduler=scheduler,
                qualification_root=qualification_root,
                provenance=wave["provenance"],
            )
            write_json_exclusive(
                Path(str(intent["paths"]["scheduler_evidence"])), scheduler
            )
            write_json_exclusive(
                Path(str(intent["paths"]["scientific_evidence"])),
                inspect_pilot_output(packet, qualification_root, wave["provenance"]),
            )
            current_job_ids.append(job_id)
        previous_job_ids = current_job_ids
    audited = audit_retained_wave(
        wave,
        qualification_root,
        scheduler_runner,
        batch_script_runner=batch_script_runner,
        now=now,
    )
    evidence_paths = []
    for case_id in sorted({case_id for case_id, _ in audited}):
        results = [
            audited[(case_id, nodes)]
            for nodes in CASE_POLICIES[case_id]["node_profiles"]
        ]
        evidence = {
            "schema_version": 2,
            "record_type": "cgl_lf_stage_i_qualification_evidence",
            "project_root": str(root),
            "qualification_root": str(qualification_root),
            "prepared_wave": retained_file(prepared_wave_path, "prepared wave"),
            "case_id": case_id,
            "target_time": CASE_POLICIES[case_id]["target_time"],
            "results": results,
        }
        path = qualification_root / f"{case_id}.qualification_evidence.json"
        write_json_exclusive(path, evidence)
        evidence_paths.append(path)
    return evidence_paths


def validate_r17_build_manifest_inventory(
    value: object,
    inventory_sha256: object,
) -> list[dict[str, object]]:
    """Authenticate the exact complete build-manifest inventory contract."""

    if not isinstance(value, list) or not value:
        raise ValueError("R17 build-manifest inventory is empty")
    inventory = []
    names = set()
    for item in value:
        record = require_exact_keys(
            item, {"name", "mode", "sha256"}, "R17 build-manifest inventory record"
        )
        name = require_nonempty_string(record["name"], "R17 build-manifest file name")
        if Path(name).name != name or name in names or record["mode"] != "0644":
            raise ValueError("R17 build-manifest inventory is incomplete or ambiguous")
        names.add(name)
        inventory.append({
            "name": name,
            "mode": "0644",
            "sha256": require_sha256(
                record["sha256"], "R17 build-manifest file sha256"
            ),
        })
    if (
        inventory != sorted(inventory, key=lambda item: str(item["name"]))
        or retained_inventory_sha256(inventory)
        != require_sha256(
            inventory_sha256, "R17 build-manifest inventory digest"
        )
    ):
        raise ValueError("R17 build-manifest inventory digest or ordering differs")
    return inventory


def validate_r17_root_relative_binding(
    value: object,
    root: Path,
    label: str,
    *,
    expected_path: Path | None = None,
    expected_mode: int | None = None,
) -> Path:
    """Authenticate one exact root-relative R17 evidence binding."""

    binding = require_exact_keys(value, {"path", "sha256"}, label)
    relative_text = require_nonempty_string(binding["path"], f"{label} path")
    relative = Path(relative_text)
    if (
        relative.is_absolute()
        or ".." in relative.parts
        or relative.as_posix() != relative_text
    ):
        raise ValueError(f"{label} path is not canonical root-relative")
    path = root / relative
    observed = root_relative_binding(path, root, label)
    if observed != binding:
        raise ValueError(f"{label} binding differs")
    if expected_path is not None and path != expected_path:
        raise ValueError(f"{label} path differs")
    if expected_mode is not None and stat.S_IMODE(path.stat().st_mode) != expected_mode:
        raise ValueError(f"{label} mode differs")
    return path


def validate_r17_rank_inventory_contract(
    value: object,
    inventory_sha256: object,
    root: Path,
    label: str,
) -> list[dict[str, str]]:
    """Authenticate exactly one terminal file for each of the 64 R17 ranks."""

    if not isinstance(value, list) or len(value) != R17_OPERATIONAL_RANKS:
        raise ValueError(f"{label} must retain exactly 64 files")
    inventory = []
    ranks = set()
    paths = set()
    for item in value:
        record = require_exact_keys(item, {"path", "sha256"}, f"{label} record")
        path = validate_r17_root_relative_binding(record, root, f"{label} file")
        rank_parts = [
            part for part in path.relative_to(root).parts
            if re.fullmatch(r"rank_[0-9]{8}", part)
        ]
        if len(rank_parts) != 1 or path in paths:
            raise ValueError(f"{label} rank inventory is incomplete or ambiguous")
        ranks.add(rank_parts[0])
        paths.add(path)
        inventory.append({
            "path": path.relative_to(root).as_posix(),
            "sha256": require_sha256(record["sha256"], f"{label} sha256"),
        })
    if (
        ranks != {f"rank_{rank:08d}" for rank in range(R17_OPERATIONAL_RANKS)}
        or inventory != sorted(inventory, key=lambda item: item["path"])
        or retained_inventory_sha256(inventory)
        != require_sha256(inventory_sha256, f"{label} inventory digest")
    ):
        raise ValueError(f"{label} rank inventory digest or ordering differs")
    return inventory


def validate_r17_frozen_science_build_contract(value: object) -> dict[str, object]:
    """Validate the exact frozen science, source, execution, and build contract."""

    contract = require_exact_keys(
        value,
        {
            "case_id", "case_name", "profile_class", "resolution", "mesh_shape",
            "meshblock_shape", "target_time", "scientific_policy", "run_basename",
            "source_revision", "source_bundle_sha256", "matrix_sha256",
            "input_sha256", "provenance_sha256", "execution_intent_sha256",
            "execution_contract_sha256", "parameter_contract",
            "parameter_contract_sha256", "executable_revision",
            "executable_sha256", "build_manifest_inventory_sha256",
        },
        "R17 frozen science/build contract",
    )
    policy = CASE_POLICIES[R17_CASE_ID]
    run_basename = f"qualification_{R17_CASE_ID}_n{R17_OPERATIONAL_NODES:02d}"
    parameter_contract = expected_case_configuration(
        R17_CASE_ID,
        basename=run_basename,
        tlim=str(policy["target_time"]),
    )
    expected_identity = {
        "case_id": R17_CASE_ID,
        "case_name": policy["case_name"],
        "profile_class": policy["profile_class"],
        "resolution": "384x384x768",
        "mesh_shape": [384, 384, 768],
        "meshblock_shape": [32, 32, 64],
        "target_time": policy["target_time"],
        "scientific_policy": policy["scientific_policy"],
        "run_basename": run_basename,
        "source_revision": FROZEN_SOURCE_REVISION,
        "matrix_sha256": FROZEN_MATRIX_SHA256,
        "input_sha256": policy["input_sha256"],
        "parameter_contract": parameter_contract,
        "parameter_contract_sha256": stable_json_sha256(parameter_contract),
    }
    executable_revision = require_git_revision(
        contract["executable_revision"], "R17 executable revision"
    )
    executable_sha256 = require_sha256(
        contract["executable_sha256"], "R17 executable sha256"
    )
    if (
        any(contract[key] != expected for key, expected in expected_identity.items())
        or executable_revision != contract["source_revision"]
        or (executable_revision, executable_sha256) not in QUALIFIED_RESTART_BINARY_ABIS
    ):
        raise ValueError("R17 frozen science/build contract differs")
    for key in (
        "source_bundle_sha256", "provenance_sha256", "execution_intent_sha256",
        "execution_contract_sha256", "build_manifest_inventory_sha256",
    ):
        require_sha256(contract[key], f"R17 frozen contract {key}")
    return contract


def build_r17_frozen_science_build_contract(
    wave: dict[str, object],
    packet: dict[str, object],
    build_manifest_inventory_sha256: str,
) -> dict[str, object]:
    """Build the exact frozen contract that one R17 qualification authenticates."""

    provenance = require_exact_keys(
        wave["provenance"],
        {
            "source", "source_bundle", "matrix", "executable",
            "build_manifest", "qualification_helper",
        },
        "prepared R17 provenance",
    )
    source = require_exact_keys(
        provenance["source"], {"directory", "revision"}, "prepared R17 source"
    )
    bundle = require_exact_keys(
        provenance["source_bundle"],
        {"path", "sha256", "size_bytes", "verified_revisions"},
        "prepared R17 source bundle",
    )
    matrix = require_exact_keys(
        provenance["matrix"], {"path", "sha256", "size_bytes"}, "prepared R17 matrix"
    )
    executable = require_exact_keys(
        provenance["executable"],
        {"path", "sha256", "size_bytes", "revision"},
        "prepared R17 executable",
    )
    intent = packet["execution_intent"]
    input_record = require_exact_keys(
        intent["input"], {"path", "sha256", "size_bytes"}, "prepared R17 input"
    )
    parameter_contract = expected_case_configuration(
        R17_CASE_ID,
        basename=str(intent["run_basename"]),
        tlim=str(intent["target_time"]),
    )
    if (
        wave["provenance_sha256"] != stable_json_sha256(provenance)
        or packet["execution_intent_sha256"] != intent["execution_intent_sha256"]
        or intent["case_id"] != R17_CASE_ID
    ):
        raise ValueError("prepared R17 execution binding differs")
    contract = {
        "case_id": R17_CASE_ID,
        "case_name": intent["case_name"],
        "profile_class": intent["profile_class"],
        "resolution": "384x384x768",
        "mesh_shape": [384, 384, 768],
        "meshblock_shape": [32, 32, 64],
        "target_time": intent["target_time"],
        "scientific_policy": intent["scientific_policy"],
        "run_basename": intent["run_basename"],
        "source_revision": source["revision"],
        "source_bundle_sha256": bundle["sha256"],
        "matrix_sha256": matrix["sha256"],
        "input_sha256": input_record["sha256"],
        "provenance_sha256": wave["provenance_sha256"],
        "execution_intent_sha256": packet["execution_intent_sha256"],
        "execution_contract_sha256": intent["execution_contract_sha256"],
        "parameter_contract": parameter_contract,
        "parameter_contract_sha256": stable_json_sha256(parameter_contract),
        "executable_revision": executable["revision"],
        "executable_sha256": executable["sha256"],
        "build_manifest_inventory_sha256": build_manifest_inventory_sha256,
    }
    return validate_r17_frozen_science_build_contract(contract)


def build_r17_independent_review_contract(
    qualification_path: Path,
    root: Path,
    measured_by: str,
    measured_utc: str,
) -> dict[str, object]:
    """Declare the exact external review required without circular self-binding."""

    reviewer_exclusions = [
        require_nonempty_string(measured_by, "R17 measurement author")
    ]
    parse_utc_timestamp(measured_utc, "R17 measurement timestamp")
    review_path = qualification_path.with_name(
        f"{qualification_path.name}.independent_review.json"
    )
    require_beneath(qualification_path, root, "R17 operational qualification")
    require_beneath(review_path, root, "R17 independent review")
    return {
        "required": True,
        "path": review_path.relative_to(root).as_posix(),
        "mode": "0444",
        "schema_version": 1,
        "record_type": "stage-i-r17-operational-qualification-independent-review",
        "execution_epoch": EXECUTION_EPOCH,
        "decision": "approved",
        "candidate_path": str(qualification_path),
        "candidate_sha256_required": True,
        "reviewer_must_differ_from": reviewer_exclusions,
        "reviewed_after_utc": measured_utc,
    }


def validate_r17_operational_qualification_contract(
    value: object,
    root: Path,
    qualification_path: Path,
) -> dict[str, object]:
    """Authenticate the one canonical schema-2 R17 producer contract."""

    qualification = require_exact_keys(
        value,
        {
            "schema_version", "record_type", "execution_epoch", "completed_utc",
            "measured_utc", "measured_by", "job_id", "state", "exit_code",
            "nodes", "ranks", "prepared_wave", "qualification_evidence",
            "scientific_evidence", "scheduler_evidence",
            "account_scheduler_evidence", "account_exclusivity_evidence",
            "executable_sha256", "build_manifest_inventory",
            "build_manifest_inventory_sha256", "rank_local_outputs",
            "rank_local_output_inventory_sha256", "rank_local_restarts",
            "rank_local_restart_inventory_sha256", "decomposition_evidence",
            "restart_load_evidence", "physics_validation_evidence",
            "frozen_science_build_contract", "independent_review_contract",
            "authority",
        },
        "R17 operational qualification",
    )
    measured_by = require_nonempty_string(
        qualification["measured_by"], "R17 measurement author"
    )
    completed = parse_scheduler_time(
        str(qualification["completed_utc"]), "R17 qualification completion"
    )
    measured = parse_utc_timestamp(
        qualification["measured_utc"], "R17 qualification measurement"
    )
    job_id = str(qualification["job_id"])
    expected_qualification_path = root / "accounting" / (
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_operational_qualification.json"
    )
    if (
        qualification_path != expected_qualification_path
        or qualification["schema_version"] != 2
        or qualification["record_type"] != "stage-i-r17-operational-qualification"
        or qualification["execution_epoch"] != EXECUTION_EPOCH
        or JOB_ID_PATTERN.fullmatch(job_id) is None
        or qualification["state"] != "COMPLETED"
        or qualification["exit_code"] != "0:0"
        or qualification["nodes"] != R17_OPERATIONAL_NODES
        or qualification["ranks"] != R17_OPERATIONAL_RANKS
        or measured < completed.astimezone(timezone.utc)
    ):
        raise ValueError("R17 operational qualification identity differs")

    prepared_wave_path = validate_r17_root_relative_binding(
        qualification["prepared_wave"], root, "R17 prepared wave"
    )
    qualification_evidence_path = validate_r17_root_relative_binding(
        qualification["qualification_evidence"], root, "R17 qualification evidence"
    )
    scientific_path = validate_r17_root_relative_binding(
        qualification["scientific_evidence"], root, "R17 scientific evidence"
    )
    scheduler_path = validate_r17_root_relative_binding(
        qualification["scheduler_evidence"],
        root,
        "R17 scheduler evidence",
        expected_path=root / "accounting" / f"{job_id}.r17_qualification.sacct.txt",
        expected_mode=0o444,
    )
    account_scheduler_path = validate_r17_root_relative_binding(
        qualification["account_scheduler_evidence"],
        root,
        "R17 account scheduler evidence",
        expected_path=(
            root / "accounting" / f"{job_id}.r17_qualification.account.sacct.txt"
        ),
        expected_mode=0o444,
    )
    account_exclusivity_path = validate_r17_root_relative_binding(
        qualification["account_exclusivity_evidence"],
        root,
        "R17 account exclusivity evidence",
        expected_path=root / "accounting" / (
            f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_account_exclusivity_evidence.json"
        ),
        expected_mode=0o444,
    )
    restart_path = validate_r17_root_relative_binding(
        qualification["restart_load_evidence"],
        root,
        "R17 restart-load evidence",
        expected_path=root / "accounting" / (
            f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_restart_load_evidence.json"
        ),
        expected_mode=0o444,
    )
    physics_path = validate_r17_root_relative_binding(
        qualification["physics_validation_evidence"],
        root,
        "R17 physics-validation evidence",
        expected_path=root / "accounting" / (
            f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_physics_validation_evidence.json"
        ),
        expected_mode=0o444,
    )

    wave = load_json(prepared_wave_path, "R17 prepared wave")
    wave_packets = []
    for wave_record in wave.get("waves", []):
        if not isinstance(wave_record, dict):
            raise ValueError("R17 prepared wave has an invalid wave record")
        for packet in wave_record.get("packets", []):
            if (
                isinstance(packet, dict)
                and isinstance(packet.get("execution_intent"), dict)
                and packet["execution_intent"].get("case_id") == R17_CASE_ID
            ):
                wave_packets.append(packet)
    if len(wave_packets) != 1:
        raise ValueError("R17 prepared wave does not retain one exact R17 packet")
    packet = wave_packets[0]
    frozen = build_r17_frozen_science_build_contract(
        wave, packet, str(qualification["build_manifest_inventory_sha256"])
    )
    if qualification["frozen_science_build_contract"] != frozen:
        raise ValueError("R17 frozen science/build binding differs from prepared wave")

    evidence = require_exact_keys(
        load_json(qualification_evidence_path, "R17 qualification evidence"),
        {
            "schema_version", "record_type", "project_root", "qualification_root",
            "prepared_wave", "case_id", "target_time", "results",
        },
        "R17 qualification evidence",
    )
    prepared_evidence_binding = require_exact_keys(
        evidence["prepared_wave"], {"path", "sha256", "size_bytes"},
        "R17 qualification evidence prepared wave",
    )
    if (
        evidence["schema_version"] != 2
        or evidence["record_type"] != "cgl_lf_stage_i_qualification_evidence"
        or evidence["project_root"] != str(root)
        or evidence["case_id"] != R17_CASE_ID
        or evidence["target_time"] != CASE_POLICIES[R17_CASE_ID]["target_time"]
        or prepared_evidence_binding["path"] != str(prepared_wave_path)
        or prepared_evidence_binding["sha256"] != qualification["prepared_wave"]["sha256"]
    ):
        raise ValueError("R17 qualification evidence binding differs")

    scientific = load_json(scientific_path, "R17 scientific evidence")
    output_inventory = validate_r17_rank_inventory_contract(
        qualification["rank_local_outputs"],
        qualification["rank_local_output_inventory_sha256"],
        root,
        "R17 terminal rank-local outputs",
    )
    restart_inventory = validate_r17_rank_inventory_contract(
        qualification["rank_local_restarts"],
        qualification["rank_local_restart_inventory_sha256"],
        root,
        "R17 terminal rank-local restarts",
    )
    build_inventory = validate_r17_build_manifest_inventory(
        qualification["build_manifest_inventory"],
        qualification["build_manifest_inventory_sha256"],
    )
    decomposition = validate_r17_decomposition_evidence(
        qualification["decomposition_evidence"],
        str(qualification["rank_local_output_inventory_sha256"]),
    )
    if (
        scientific.get("schema_version") != 6
        or scientific.get("record_type")
        != "cgl_lf_stage_i_qualification_scientific_evidence"
        or scientific.get("case_id") != R17_CASE_ID
        or scientific.get("nodes") != R17_OPERATIONAL_NODES
        or scientific.get("execution_intent_sha256")
        != frozen["execution_intent_sha256"]
        or scientific.get("execution_contract_sha256")
        != frozen["execution_contract_sha256"]
        or scientific.get("terminal_rank_local_outputs") != output_inventory
        or scientific.get("terminal_rank_local_output_inventory_sha256")
        != qualification["rank_local_output_inventory_sha256"]
        or scientific.get("terminal_rank_local_restarts") != restart_inventory
        or scientific.get("terminal_rank_local_restart_inventory_sha256")
        != qualification["rank_local_restart_inventory_sha256"]
        or scientific.get("r17_decomposition") != decomposition
        or scientific.get("accepted_for_operational_qualification") is not True
        or scientific.get("accepted_for_profile_selection") is not False
        or qualification["executable_sha256"] != frozen["executable_sha256"]
        or qualification["build_manifest_inventory_sha256"]
        != frozen["build_manifest_inventory_sha256"]
        or build_inventory != qualification["build_manifest_inventory"]
    ):
        raise ValueError("R17 scientific/build evidence differs from canonical contract")

    try:
        scheduler_text = read_regular_bytes(
            scheduler_path, "R17 scheduler evidence", single_link=True
        ).decode("utf-8")
        account_scheduler_text = read_regular_bytes(
            account_scheduler_path, "R17 account scheduler evidence", single_link=True
        ).decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("R17 retained scheduler evidence is not UTF-8") from error
    account_exclusivity = load_json(
        account_exclusivity_path, "R17 account exclusivity evidence"
    )
    account_job = account_exclusivity.get("qualification_job")
    if not isinstance(account_job, dict):
        raise ValueError("R17 account exclusivity evidence lacks the qualification job")
    validate_r17_account_exclusivity_evidence(
        account_exclusivity, account_job, account_scheduler_text
    )
    scheduler_lines = scheduler_text.splitlines()
    scheduler_expected = "|".join([
        job_id,
        str(account_job.get("job_name")),
        "COMPLETED",
        "0:0",
        str(R17_OPERATIONAL_NODES),
        str(account_job.get("elapsed_seconds")),
        str(account_job.get("submit_utc")),
        str(account_job.get("end_utc")),
    ])
    if (
        scheduler_lines != [scheduler_expected]
        or account_job.get("job_id") != job_id
        or account_job.get("state") != "COMPLETED"
        or account_job.get("exit_code") != "0:0"
        or account_job.get("nodes") != R17_OPERATIONAL_NODES
        or account_job.get("end_utc") != qualification["completed_utc"]
        or account_exclusivity.get("measured_utc") != qualification["measured_utc"]
        or account_exclusivity.get("exclusive_entire_execution_interval") is not True
        or account_exclusivity.get("overlapping_job_ids") != [job_id]
        or account_exclusivity.get("query_contract", {}).get("all_users") is not True
        or account_exclusivity.get("visibility_contract", {}).get(
            "all_users_job_visibility"
        )
        is not True
    ):
        raise ValueError("R17 full-interval account exclusivity contract differs")

    restart = require_exact_keys(
        load_json(restart_path, "R17 restart-load evidence"),
        {
            "schema_version", "record_type", "execution_epoch", "measured_utc",
            "measured_by", "job_id", "executable_sha256",
            "build_manifest_inventory_sha256", "passed", "measurements",
        },
        "R17 restart-load evidence",
    )
    physics = require_exact_keys(
        load_json(physics_path, "R17 physics-validation evidence"),
        {
            "schema_version", "record_type", "execution_epoch", "measured_utc",
            "measured_by", "job_id", "executable_sha256",
            "build_manifest_inventory_sha256", "passed", "measurements",
        },
        "R17 physics-validation evidence",
    )
    expected_restart_measurements = {
        "rank_local_restart_inventory_sha256": qualification[
            "rank_local_restart_inventory_sha256"
        ],
        "loaded_rank_count": R17_OPERATIONAL_RANKS,
        "load_state": "COMPLETED",
        "load_exit_code": "0:0",
    }
    expected_physics_measurements = {
        "rank_local_output_inventory_sha256": qualification[
            "rank_local_output_inventory_sha256"
        ],
        **scientific["physics_measurements"],
    }
    for record, record_type, measurements in (
        (restart, "stage-i-r17-restart-load-evidence", expected_restart_measurements),
        (physics, "stage-i-r17-physics-validation-evidence", expected_physics_measurements),
    ):
        if (
            record["schema_version"] != 1
            or record["record_type"] != record_type
            or record["execution_epoch"] != EXECUTION_EPOCH
            or record["measured_utc"] != qualification["measured_utc"]
            or record["measured_by"] != measured_by
            or record["job_id"] != job_id
            or record["executable_sha256"] != qualification["executable_sha256"]
            or record["build_manifest_inventory_sha256"]
            != qualification["build_manifest_inventory_sha256"]
            or record["passed"] is not True
            or record["measurements"] != measurements
        ):
            raise ValueError("R17 retained validation evidence differs")

    expected_review = build_r17_independent_review_contract(
        qualification_path, root, measured_by, str(qualification["measured_utc"])
    )
    if (
        qualification["independent_review_contract"] != expected_review
        or qualification["authority"]
        != {
            "r17_launch_authorized": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
    ):
        raise ValueError("R17 independent-review or authority contract differs")
    return qualification


def publish_r17_operational_qualification(
    root: Path,
    artifacts: list[tuple[Path, bytes]],
    qualification_path: Path,
    qualification: dict[str, object],
    review_path: Path,
) -> None:
    """Publish one descriptor-bound canonical R17 evidence set under the Stage I lock."""

    accounting = root / "accounting"
    require_no_symlink_chain(accounting, root, "Stage I accounting directory")
    with stage_i_exclusivity_lock(root) as lock, open_directory_fd(
        accounting, "Stage I accounting directory"
    ) as (accounting_descriptor, accounting_profile):
        def require_canonical_mutation_boundary() -> None:
            lock.require_bound()
            require_directory_descriptor_binding(
                accounting,
                accounting_descriptor,
                accounting_profile,
                "Stage I accounting directory",
            )

        require_empty_canonical_stage_i_state(root)
        require_canonical_mutation_boundary()
        destinations = [path for path, _ in artifacts]
        if (
            qualification_path.parent != accounting
            or review_path.parent != accounting
            or any(path.parent != accounting for path in destinations)
        ):
            raise ValueError("R17 operational evidence escaped canonical accounting")
        for path in [*destinations, qualification_path, review_path]:
            if path.exists() or path.is_symlink():
                raise ValueError(f"R17 operational evidence already exists: {path}")
        require_canonical_mutation_boundary()
        for path, payload in artifacts:
            write_exclusive(
                path,
                payload,
                mode=0o444,
                mutation_guard=require_canonical_mutation_boundary,
            )
        validate_r17_operational_qualification_contract(
            qualification, root, qualification_path
        )
        write_exclusive(
            qualification_path,
            stable_json(qualification).encode("utf-8"),
            mode=0o444,
            mutation_guard=require_canonical_mutation_boundary,
        )


def retain_r17_operational_qualification(
    evidence_path: Path,
    measured_by: str,
    scheduler_runner=subprocess.run,
    batch_script_runner=None,
    now: datetime | None = None,
) -> dict[str, object]:
    """Retain measured, non-self-approved R17 operational qualification evidence."""

    measured_by = require_nonempty_string(measured_by, "R17 measurement author")
    current = now or datetime.now(timezone.utc)
    if current.tzinfo is None:
        raise ValueError("R17 operational qualification time lacks an explicit offset")
    current = current.astimezone(timezone.utc)
    if batch_script_runner is None:
        batch_script_runner = scheduler_runner
    evidence = require_exact_keys(
        load_json(evidence_path, "R17 qualification evidence"),
        {
            "schema_version", "record_type", "project_root", "qualification_root",
            "prepared_wave", "case_id", "target_time", "results",
        },
        "R17 qualification evidence",
    )
    if (
        evidence["schema_version"] != 2
        or evidence["record_type"] != "cgl_lf_stage_i_qualification_evidence"
        or evidence["case_id"] != R17_CASE_ID
    ):
        raise ValueError("operational qualification requires exact R17 evidence")
    root = require_project_root(str(evidence["project_root"]), allow_local_root=False)
    qualification_root = require_qualification_root(
        root, str(evidence["qualification_root"])
    )
    expected_evidence = qualification_root / f"{R17_CASE_ID}.qualification_evidence.json"
    if require_owned_regular_file(
        evidence_path, qualification_root, "R17 qualification evidence"
    ) != expected_evidence:
        raise ValueError("R17 qualification evidence path differs")
    wave_path = load_bound_evidence_file(
        evidence["prepared_wave"], qualification_root, "prepared R17 wave"
    )
    wave = load_json(wave_path, "prepared R17 wave")
    packets = validate_prepared_wave(wave)
    if set(packets) != {(R17_CASE_ID, R17_OPERATIONAL_NODES)}:
        raise ValueError("operational qualification wave is not exact R17")
    audited = audit_retained_wave(
        wave,
        qualification_root,
        scheduler_runner,
        batch_script_runner=batch_script_runner,
        now=current,
    )
    result = audited.get((R17_CASE_ID, R17_OPERATIONAL_NODES))
    if evidence["results"] != [result]:
        raise ValueError("R17 qualification evidence differs from reproduced audit")
    packet = packets[(R17_CASE_ID, R17_OPERATIONAL_NODES)]
    intent = packet["execution_intent"]
    scheduler = load_json(
        Path(str(intent["paths"]["scheduler_evidence"])), "R17 scheduler evidence"
    )
    scientific = load_json(
        Path(str(intent["paths"]["scientific_evidence"])), "R17 scientific evidence"
    )
    if (
        scientific.get("accepted_for_operational_qualification") is not True
        or scientific.get("accepted_for_profile_selection") is not False
        or scheduler.get("state") != "COMPLETED"
        or scheduler.get("exit_code") != "0:0"
        or scheduler.get("nodes") != R17_OPERATIONAL_NODES
        or int(intent["allocation"]["nodes"]) * int(intent["allocation"]["ranks_per_node"])
        != R17_OPERATIONAL_RANKS
    ):
        raise ValueError("R17 evidence does not prove exact operational readiness")
    completed = parse_scheduler_time(
        str(scheduler["end_utc"]), "R17 scheduler completion"
    ).astimezone(timezone.utc)
    if current < completed:
        raise ValueError("R17 operational evidence predates scheduler completion")
    measured_utc = format_utc_timestamp(current)
    build = require_exact_keys(
        wave["provenance"]["build_manifest"],
        {
            "path", "athena_sha256", "environment", "inventory",
            "inventory_sha256",
        },
        "R17 build manifest",
    )
    inventory = build["inventory"]
    if (
        not isinstance(inventory, list)
        or build["inventory_sha256"] != retained_inventory_sha256(inventory)
    ):
        raise ValueError("R17 build-manifest inventory digest differs")
    executable = wave["provenance"]["executable"]
    output_inventory = scientific["terminal_rank_local_outputs"]
    restart_inventory = scientific["terminal_rank_local_restarts"]
    output_inventory_sha256 = str(
        scientific["terminal_rank_local_output_inventory_sha256"]
    )
    restart_inventory_sha256 = str(
        scientific["terminal_rank_local_restart_inventory_sha256"]
    )
    decomposition_evidence = validate_r17_decomposition_evidence(
        scientific.get("r17_decomposition"),
        output_inventory_sha256,
    )
    if (
        len(output_inventory) != R17_OPERATIONAL_RANKS
        or len(restart_inventory) != R17_OPERATIONAL_RANKS
        or output_inventory_sha256 != retained_inventory_sha256(output_inventory)
        or restart_inventory_sha256 != retained_inventory_sha256(restart_inventory)
    ):
        raise ValueError("R17 terminal rank-local inventory differs")
    physics_measurements = require_exact_keys(
        scientific["physics_measurements"],
        {
            "finite_rank_outputs", "mass_relative_drift_max",
            "mhd_user_mass_mismatch_max", "lf_bad_counts_total",
            "normalized_ct_divb_max", "normalized_ct_divb_threshold",
            "normalized_ct_divb_below_threshold",
        },
        "R17 physics measurements",
    )
    account_visibility = collect_live_account_visibility(scheduler_runner)
    account_scheduler_raw, account_query_contract = collect_live_account_scheduler_raw(
        parse_scheduler_time(str(scheduler["start_utc"]), "R17 scheduler start"),
        parse_scheduler_time(str(scheduler["end_utc"]), "R17 scheduler end"),
        scheduler_runner,
    )
    account_exclusivity_evidence = build_r17_account_exclusivity_evidence(
        scheduler,
        account_scheduler_raw,
        account_query_contract,
        account_visibility,
        measured_utc,
    )
    restart_evidence = {
        "schema_version": 1,
        "record_type": "stage-i-r17-restart-load-evidence",
        "execution_epoch": EXECUTION_EPOCH,
        "measured_utc": measured_utc,
        "measured_by": measured_by,
        "job_id": scheduler["job_id"],
        "executable_sha256": executable["sha256"],
        "build_manifest_inventory_sha256": build["inventory_sha256"],
        "passed": True,
        "measurements": {
            "rank_local_restart_inventory_sha256": scientific[
                "terminal_rank_local_restart_inventory_sha256"
            ],
            "loaded_rank_count": R17_OPERATIONAL_RANKS,
            "load_state": "COMPLETED",
            "load_exit_code": "0:0",
        },
    }
    physics_evidence = {
        "schema_version": 1,
        "record_type": "stage-i-r17-physics-validation-evidence",
        "execution_epoch": EXECUTION_EPOCH,
        "measured_utc": measured_utc,
        "measured_by": measured_by,
        "job_id": scheduler["job_id"],
        "executable_sha256": executable["sha256"],
        "build_manifest_inventory_sha256": build["inventory_sha256"],
        "passed": True,
        "measurements": {
            "rank_local_output_inventory_sha256": scientific[
                "terminal_rank_local_output_inventory_sha256"
            ],
            **physics_measurements,
        },
    }
    accounting = root / "accounting"
    require_directory(accounting, "Stage I accounting directory")
    scheduler_path = accounting / f"{scheduler['job_id']}.r17_qualification.sacct.txt"
    account_scheduler_path = accounting / (
        f"{scheduler['job_id']}.r17_qualification.account.sacct.txt"
    )
    account_exclusivity_path = accounting / (
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_account_exclusivity_evidence.json"
    )
    restart_path = accounting / (
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_restart_load_evidence.json"
    )
    physics_path = accounting / (
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_physics_validation_evidence.json"
    )
    qualification_path = accounting / (
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_operational_qualification.json"
    )
    review_path = qualification_path.with_name(
        f"{qualification_path.name}.independent_review.json"
    )
    scheduler_row = "|".join([
        str(scheduler["job_id"]),
        str(scheduler["job_name"]),
        str(scheduler["state"]),
        str(scheduler["exit_code"]),
        str(scheduler["nodes"]),
        str(scheduler["elapsed_seconds"]),
        str(scheduler["submit_utc"]),
        str(scheduler["end_utc"]),
    ]) + "\n"
    scheduler_payload = scheduler_row.encode("utf-8")
    account_scheduler_payload = account_scheduler_raw.encode("utf-8")
    account_exclusivity_payload = stable_json(
        account_exclusivity_evidence
    ).encode("utf-8")
    restart_payload = stable_json(restart_evidence).encode("utf-8")
    physics_payload = stable_json(physics_evidence).encode("utf-8")
    qualification = {
        "schema_version": 2,
        "record_type": "stage-i-r17-operational-qualification",
        "execution_epoch": EXECUTION_EPOCH,
        "completed_utc": scheduler["end_utc"],
        "measured_utc": measured_utc,
        "measured_by": measured_by,
        "job_id": scheduler["job_id"],
        "state": scheduler["state"],
        "exit_code": scheduler["exit_code"],
        "nodes": R17_OPERATIONAL_NODES,
        "ranks": R17_OPERATIONAL_RANKS,
        "prepared_wave": root_relative_binding(
            wave_path, root, "prepared R17 wave"
        ),
        "qualification_evidence": root_relative_binding(
            evidence_path, root, "R17 qualification evidence"
        ),
        "scientific_evidence": root_relative_binding(
            Path(str(intent["paths"]["scientific_evidence"])),
            root,
            "R17 scientific evidence",
        ),
        "scheduler_evidence": root_relative_payload_binding(
            scheduler_path, scheduler_payload, root, "R17 scheduler evidence"
        ),
        "account_scheduler_evidence": root_relative_payload_binding(
            account_scheduler_path,
            account_scheduler_payload,
            root,
            "R17 account scheduler evidence",
        ),
        "account_exclusivity_evidence": root_relative_payload_binding(
            account_exclusivity_path,
            account_exclusivity_payload,
            root,
            "R17 account exclusivity evidence",
        ),
        "executable_sha256": executable["sha256"],
        "build_manifest_inventory": inventory,
        "build_manifest_inventory_sha256": build["inventory_sha256"],
        "rank_local_outputs": output_inventory,
        "rank_local_output_inventory_sha256": output_inventory_sha256,
        "rank_local_restarts": restart_inventory,
        "rank_local_restart_inventory_sha256": restart_inventory_sha256,
        "decomposition_evidence": decomposition_evidence,
        "restart_load_evidence": root_relative_payload_binding(
            restart_path, restart_payload, root, "R17 restart-load evidence"
        ),
        "physics_validation_evidence": root_relative_payload_binding(
            physics_path, physics_payload, root, "R17 physics-validation evidence"
        ),
        "frozen_science_build_contract": build_r17_frozen_science_build_contract(
            wave, packet, str(build["inventory_sha256"])
        ),
        "independent_review_contract": build_r17_independent_review_contract(
            qualification_path, root, measured_by, measured_utc
        ),
        "authority": {
            "r17_launch_authorized": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        },
    }
    publish_r17_operational_qualification(
        root,
        [
            (scheduler_path, scheduler_payload),
            (account_scheduler_path, account_scheduler_payload),
            (account_exclusivity_path, account_exclusivity_payload),
            (restart_path, restart_payload),
            (physics_path, physics_payload),
        ],
        qualification_path,
        qualification,
        review_path,
    )
    return {
        "schema_version": 2,
        "record_type": "cgl_lf_stage_i_r17_operational_qualification_retention",
        "operational_qualification": root_relative_binding(
            qualification_path, root, "R17 operational qualification"
        ),
        "restart_load_evidence": qualification["restart_load_evidence"],
        "physics_validation_evidence": qualification["physics_validation_evidence"],
        "scheduler_evidence": qualification["scheduler_evidence"],
        "account_scheduler_evidence": qualification["account_scheduler_evidence"],
        "account_exclusivity_evidence": qualification["account_exclusivity_evidence"],
        "required_independent_review": qualification["independent_review_contract"],
        "authority": qualification["authority"],
        "independent_review_created": False,
        "self_approved": False,
        "launch_authorized": False,
    }


def select_profile(evidence_path: Path,
                   scheduler_runner=subprocess.run,
                   batch_script_runner=None,
                   now: datetime | None = None) -> dict[str, object]:
    """Select one node profile only from a reproduced complete wave audit."""

    current = now or datetime.now(timezone.utc)
    if current.tzinfo is None:
        raise ValueError("selection reference time lacks an explicit UTC offset")
    current = current.astimezone(timezone.utc)
    if batch_script_runner is None:
        batch_script_runner = scheduler_runner
    evidence = load_json(evidence_path, "qualification evidence")
    value = require_exact_keys(
        evidence,
        {
            "schema_version", "record_type", "project_root", "qualification_root",
            "prepared_wave", "case_id", "target_time", "results",
        },
        "qualification evidence",
    )
    if (
        value["schema_version"] != 2
        or value["record_type"] != "cgl_lf_stage_i_qualification_evidence"
    ):
        raise ValueError("qualification evidence has wrong schema")
    root = require_project_root(str(value["project_root"]), allow_local_root=False)
    qualification_root = require_qualification_root(
        root, str(value["qualification_root"])
    )
    require_qualification_path(
        evidence_path, qualification_root, "qualification evidence", must_exist=True
    )
    wave_path = load_bound_evidence_file(
        value["prepared_wave"], qualification_root, "prepared wave"
    )
    wave = load_json(wave_path, "prepared wave")
    audited = audit_retained_wave(
        wave,
        qualification_root,
        scheduler_runner,
        batch_script_runner=batch_script_runner,
        now=current,
    )
    if (
        wave["project_root"] != str(root)
        or wave["qualification_root"] != str(qualification_root)
    ):
        raise ValueError("qualification evidence and prepared wave roots differ")
    case_id = str(value["case_id"])
    if case_id not in CASE_POLICIES:
        raise ValueError("qualification evidence has an unsupported case")
    if CASE_POLICIES[case_id]["operational_only"]:
        raise ValueError("R17 operational qualification is not profile-selection evidence")
    target = require_finite_positive(value["target_time"], "evidence target time")
    if target != CASE_POLICIES[case_id]["target_time"]:
        raise ValueError("qualification evidence target time differs from policy")
    expected_evidence_path = qualification_root / f"{case_id}.qualification_evidence.json"
    if evidence_path.resolve(strict=True) != expected_evidence_path.resolve(strict=True):
        raise ValueError("qualification evidence path differs from isolated layout")
    expected_nodes = set(CASE_POLICIES[case_id]["node_profiles"])
    if not isinstance(value["results"], list) or len(value["results"]) != len(expected_nodes):
        raise ValueError("qualification evidence lacks one result per node profile")
    observed: dict[int, dict[str, object]] = {}
    for result in value["results"]:
        if not isinstance(result, dict):
            raise ValueError("qualification profile result is invalid")
        nodes = result.get("nodes")
        if isinstance(nodes, bool) or not isinstance(nodes, int) or nodes not in expected_nodes:
            raise ValueError("qualification result has unauthorized node profile")
        if nodes in observed or result != audited.get((case_id, nodes)):
            raise ValueError("qualification result differs from reproduced wave audit")
        observed[nodes] = {
            "nodes": nodes,
            "elapsed_seconds": require_finite_positive(
                result["elapsed_seconds"], f"{nodes}-node elapsed seconds"
            ),
            "job_id": result["job_id"],
            "execution_intent_sha256": result["execution_intent_sha256"],
            "scheduler_end_utc": result["scheduler_end_utc"],
            "scientific_agreement_signature": result[
                "scientific_agreement_signature"
            ],
        }
    if set(observed) != expected_nodes:
        raise ValueError("qualification evidence node profiles are incomplete")
    fixed_fresh_profile = (
        require_r12_fresh_profile_policy() if case_id == R12_CASE_ID else None
    )
    baseline_nodes = (
        int(fixed_fresh_profile["nodes"]) if fixed_fresh_profile is not None else 1
    )
    scientific_agreement = cross_profile_scientific_agreement(
        observed, baseline_nodes=baseline_nodes
    )
    baseline = observed[baseline_nodes]["elapsed_seconds"]
    evaluated = []
    eligible = []
    for nodes in sorted(observed):
        runtime = observed[nodes]["elapsed_seconds"]
        savings = baseline - runtime
        ratio = nodes * runtime / (baseline_nodes * baseline)
        larger_eligible = (
            nodes > baseline_nodes
            and savings >= MINIMUM_RUNTIME_SAVINGS_SECONDS
            and ratio <= MAXIMUM_NODE_HOUR_RATIO
        )
        record = {
            **observed[nodes],
            "runtime_savings_seconds_vs_n1": savings,
            "node_hour_ratio_vs_n1": ratio,
            "larger_profile_eligible": larger_eligible,
        }
        evaluated.append(record)
        if larger_eligible:
            eligible.append(record)
    if fixed_fresh_profile is not None:
        selected = next(
            item for item in evaluated if item["nodes"] == fixed_fresh_profile["nodes"]
        )
        reason = "selected exact measured parentless fresh R12 profile"
    elif not eligible:
        selected = next(item for item in evaluated if item["nodes"] == baseline_nodes)
        reason = "no larger profile satisfies both promotion thresholds"
    else:
        fastest = min(item["elapsed_seconds"] for item in eligible)
        tied = [
            item for item in eligible
            if item["elapsed_seconds"] <= fastest * (1.0 + RUNTIME_TIE_FRACTION)
        ]
        selected = min(
            tied, key=lambda item: (item["nodes"], item["elapsed_seconds"])
        )
        reason = (
            "selected fastest eligible larger profile; runtimes within 10% "
            "prefer fewer nodes"
        )
    evidence_ends = [
        require_recent_scheduler_end(item["scheduler_end_utc"], now=current)
        for item in observed.values()
    ]
    expires = min(
        current + timedelta(seconds=SELECTION_VALIDITY_SECONDS),
        *(end + timedelta(seconds=EVIDENCE_MAX_AGE_SECONDS) for end in evidence_ends),
    )
    case_policy = CASE_POLICIES[case_id]
    selection_core = {
        "schema_version": 3,
        "record_type": "cgl_lf_stage_i_qualification_selection",
        "generated_utc": format_utc_timestamp(current),
        "expires_utc": format_utc_timestamp(expires),
        "project_root": str(root),
        "qualification_root": str(qualification_root),
        "case_id": case_id,
        "target_time": target,
        "prepared_wave": value["prepared_wave"],
        "qualification_evidence": retained_file(
            evidence_path, "qualification evidence"
        ),
        "policy": {
            "minimum_runtime_savings_seconds": MINIMUM_RUNTIME_SAVINGS_SECONDS,
            "maximum_node_hour_ratio": MAXIMUM_NODE_HOUR_RATIO,
            "runtime_tie_fraction": RUNTIME_TIE_FRACTION,
            "evidence_max_age_seconds": EVIDENCE_MAX_AGE_SECONDS,
            "selection_validity_seconds": SELECTION_VALIDITY_SECONDS,
            "cross_profile_relative_tolerance": CROSS_PROFILE_RELATIVE_TOLERANCE,
            "cross_profile_absolute_tolerance": CROSS_PROFILE_ABSOLUTE_TOLERANCE,
        },
        "profiles": evaluated,
        "cross_profile_scientific_agreement": scientific_agreement,
        "selected_nodes": selected["nodes"],
        "selection_reason": reason,
        "consumption_contract": {
            "schema_version": 1,
            "record_type": (
                "cgl_lf_stage_i_non_authorizing_profile_consumption_contract"
            ),
            "advisory_only": True,
            "canonical_acceptance_eligible": False,
            "production_authorization": False,
            "production_controller_consumption_implemented": False,
            "reviewed_promotion_required": True,
            "cross_profile_scientific_agreement_required": True,
            "fixed_fresh_profile": fixed_fresh_profile,
            "selected_profile": {
                "case_id": case_id,
                "profile_class": case_policy["profile_class"],
                "nodes": selected["nodes"],
                "ranks_per_node": RANKS_PER_NODE,
                "cpus_per_task": CPUS_PER_TASK,
                "target_time": target,
                "walltime": case_policy["walltime"],
                "athena_walltime": case_policy["athena_walltime"],
            },
        },
    }
    return {
        **selection_core,
        "selection_contract_sha256": stable_json_sha256(selection_core),
    }


def retain_profile_selection(evidence_path: Path,
                             scheduler_runner=subprocess.run,
                             batch_script_runner=None,
                             now: datetime | None = None) -> dict[str, object]:
    """Retain one immutable expiring advisory selection artifact."""

    selection = select_profile(
        evidence_path,
        scheduler_runner,
        batch_script_runner=batch_script_runner,
        now=now,
    )
    root = require_project_root(str(selection["project_root"]), allow_local_root=False)
    qualification_root = require_qualification_root(
        root, str(selection["qualification_root"])
    )
    path = qualification_root / f"{selection['case_id']}.qualification_selection.json"
    write_exclusive(path, stable_json(selection).encode("utf-8"), mode=0o440)
    return {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_retained_qualification_selection",
        "selection": selection,
        "retained_selection": retained_file(path, "retained qualification selection"),
    }


def validate_profile_selection(selection_path: Path,
                               scheduler_runner=subprocess.run,
                               batch_script_runner=None,
                               now: datetime | None = None) -> dict[str, object]:
    """Reproduce one durable selection while preserving non-authorizing status."""

    current = now or datetime.now(timezone.utc)
    if current.tzinfo is None:
        raise ValueError("selection verification time lacks an explicit UTC offset")
    current = current.astimezone(timezone.utc)
    if batch_script_runner is None:
        batch_script_runner = scheduler_runner
    selection = require_exact_keys(
        load_json(selection_path, "retained qualification selection"),
        {
            "schema_version", "record_type", "generated_utc", "expires_utc",
            "project_root", "qualification_root", "case_id", "target_time",
            "prepared_wave", "qualification_evidence", "policy", "profiles",
            "cross_profile_scientific_agreement", "selected_nodes",
            "selection_reason", "consumption_contract",
            "selection_contract_sha256",
        },
        "retained qualification selection",
    )
    retained_digest = require_sha256(
        selection["selection_contract_sha256"], "selection contract digest"
    )
    digest_core = dict(selection)
    digest_core.pop("selection_contract_sha256")
    if stable_json_sha256(digest_core) != retained_digest:
        raise ValueError("retained qualification selection digest differs")
    generated = parse_utc_timestamp(selection["generated_utc"], "selection generation time")
    expires = parse_utc_timestamp(selection["expires_utc"], "selection expiry time")
    if (
        expires <= generated
        or current < generated - timedelta(seconds=FUTURE_TIMESTAMP_TOLERANCE_SECONDS)
        or current > expires
    ):
        raise ValueError("retained qualification selection is expired or future-dated")
    root = require_project_root(str(selection["project_root"]), allow_local_root=False)
    qualification_root = require_qualification_root(
        root, str(selection["qualification_root"])
    )
    expected_path = qualification_root / (
        f"{selection['case_id']}.qualification_selection.json"
    )
    selection_path = require_owned_regular_file(
        selection_path, qualification_root, "retained qualification selection"
    )
    if (
        selection_path != expected_path
        or selection_path.stat().st_mode & 0o222
    ):
        raise ValueError("retained qualification selection path differs")
    evidence_path = validate_retained_file(
        selection["qualification_evidence"],
        qualification_root,
        "selection qualification evidence",
    )
    reproduced = select_profile(
        evidence_path,
        scheduler_runner,
        batch_script_runner=batch_script_runner,
        now=generated,
    )
    if selection != reproduced:
        raise ValueError("retained qualification selection differs from reproduced audit")
    return selection


def parser() -> argparse.ArgumentParser:
    """Build the isolated qualification lifecycle command-line parser."""

    command = argparse.ArgumentParser(description=__doc__)
    actions = command.add_subparsers(dest="action", required=True)
    prepare = actions.add_parser(
        "prepare-wave",
        help="Retain deterministic authenticated qualification packets.",
    )
    prepare.add_argument("--root", default=str(DEFAULT_ROOT))
    prepare.add_argument("--allow-local-root", action="store_true")
    prepare.add_argument("--qualification-root", required=True)
    prepare.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    prepare.add_argument("--source-dir", default=str(REPOSITORY_ROOT))
    prepare.add_argument("--source-bundle", required=True)
    prepare.add_argument("--executable", required=True)
    prepare.add_argument("--build-manifest", required=True)
    prepare.add_argument("--case-id", action="append", required=True)
    prepare.add_argument("--max-concurrent-nodes", type=int, default=10)
    submit = actions.add_parser(
        "submit-all-waves",
        help="Submit retained waves with exact afterok dependency barriers.",
    )
    submit.add_argument("--prepared-wave", required=True)
    recover = actions.add_parser(
        "recover-submit",
        help="Authenticate and cancel an operator-identified ambiguous submission.",
    )
    recover.add_argument("--prepared-wave", required=True)
    recover.add_argument("--packet-id", required=True)
    recover.add_argument("--job-id", required=True)
    recover.add_argument("--notes", required=True)
    audit = actions.add_parser(
        "audit-wave",
        help="Collect live sacct records and reproduce every scientific inspection.",
    )
    audit.add_argument("--prepared-wave", required=True)
    select = actions.add_parser(
        "select-profile",
        help="Retain an expiring non-authorizing profile recommendation.",
    )
    select.add_argument("--evidence", required=True)
    operational = actions.add_parser(
        "retain-r17-operational-qualification",
        help="Retain recost-compatible measured R17 operational evidence.",
    )
    operational.add_argument("--evidence", required=True)
    operational.add_argument("--measured-by", required=True)
    verify = actions.add_parser(
        "verify-selection",
        help="Reproduce and verify a retained advisory profile selection.",
    )
    verify.add_argument("--selection", required=True)
    return command


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point."""

    retained_argv = sys.argv[1:] if argv is None else argv
    try:
        authenticate_self(retained_argv)
        args = parser().parse_args(retained_argv)
        if args.action == "prepare-wave":
            result = prepare_wave(args)
        elif args.action == "submit-all-waves":
            result = submit_all_waves(Path(args.prepared_wave).expanduser())
        elif args.action == "recover-submit":
            result = recover_ambiguous_submission(
                Path(args.prepared_wave).expanduser(),
                args.packet_id,
                args.job_id,
                args.notes,
            )
        elif args.action == "audit-wave":
            paths = retain_wave_audit(Path(args.prepared_wave).expanduser())
            result = {
                "schema_version": 1,
                "record_type": "cgl_lf_stage_i_qualification_audit",
                "evidence": [retained_file(path, "qualification evidence") for path in paths],
            }
        elif args.action == "select-profile":
            result = retain_profile_selection(Path(args.evidence).expanduser())
        elif args.action == "retain-r17-operational-qualification":
            result = retain_r17_operational_qualification(
                Path(args.evidence).expanduser(), args.measured_by
            )
        elif args.action == "verify-selection":
            result = validate_profile_selection(Path(args.selection).expanduser())
        else:
            raise ValueError(f"unsupported action: {args.action}")
        sys.stdout.write(stable_json(result))
        return 0
    except (
        KeyError,
        OSError,
        TypeError,
        ValueError,
        subprocess.CalledProcessError,
    ) as error:
        print(f"Stage I qualification utility failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
