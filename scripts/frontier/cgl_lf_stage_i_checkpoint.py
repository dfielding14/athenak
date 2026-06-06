#!/opt/cray/pe/python/3.11.7/bin/python3.11 -I
"""Verify and publish retained Stage I recost checkpoints.

The companion is intentionally generic: callers supply the retained artifact
name, authenticated file digests, one legacy launch packet or schema-2 strict
non-authorizing recommendation packet, expected reconcile counts, and JSON
pointers locating those bindings in the artifact.  Publishing recommendation
evidence never grants launch authority.  Canonical publication uses a journaled
same-directory link/fsync/copy/exchange/forensic-retirement/fsync sequence
under the Stage I lock.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import ctypes
from datetime import datetime, timedelta, timezone
import errno
import fcntl
import hashlib
from importlib.machinery import SourceFileLoader
import importlib.util
import json
import os
from pathlib import Path
import pwd
import re
import stat
import subprocess
import sys
import tempfile
import uuid


DEFAULT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
CANONICAL_REPOSITORY_ROOT = Path("/autofs/nccs-svm1_home2/dfielding/athenak-df")
SQUEUE = Path("/usr/bin/squeue")
GIT = Path("/usr/lib/git/git")
GIT_EXEC_PATH = Path("/usr/lib/git")
STAGE_I_PYTHON = Path("/usr/bin/python3.11")
TRUSTED_SYSTEM_PATH = "/usr/bin:/bin"
OFFLINE_ALLOWED_PREFIXES = (Path("/tmp"),)
EXECUTION_EPOCH = "E03-forcing-policy"
EXECUTION_EPOCH_SLUG = "E03_forcing_policy"
STAGE_I_RELATIVE = Path("scripts/frontier/cgl_lf_stage_i.py")
UTILITY_RELATIVE = Path("scripts/frontier/cgl_lf_stage_i_checkpoint.py")
RECOST_RELATIVE = Path("scripts/frontier/cgl_lf_stage_i_recost.py")
F113_RELATIVE = Path(
    "accounting/"
    "mks24_stage_i_E03_forcing_policy_F113_controller_transition_evidence.json"
)
F113_PUBLICATION_AUDIT_RELATIVE = Path(f"{F113_RELATIVE}.publication_audit.json")
F116_RELATIVE = Path(
    "accounting/"
    "mks24_stage_i_E03_forcing_policy_F116_current_source_authority_supersession_evidence.json"
)
F116_PUBLICATION_AUDIT_RELATIVE = Path(f"{F116_RELATIVE}.publication_audit.json")
F116_PROVENANCE_REVIEW_RELATIVE = Path(f"{F116_RELATIVE}.provenance_security_review.json")
F116_PLASMA_REVIEW_RELATIVE = Path(f"{F116_RELATIVE}.plasma_scientific_review.json")
F116_CANONICAL_SHA256 = {
    "evidence": "6cdbf9e4d10f1282744c6274aa3ef08afec4c510420296837fdbbdfcefe30a2a",
    "provenance_review": "e9731aab8305505e058c68ae8bb61c9ec5ff4885bde1bee3162719c41ab9bafd",
    "plasma_review": "bc4da6897263843f5d233f996539ff047d6ef465b775a9b5cdc8f864de96a6b8",
    "publication_audit": "3a6168e3039c02656b38ebdfcadffc07a2f151b1a474084307a80a341ba83096",
}
F118_RELATIVE = Path(
    "accounting/"
    "mks24_stage_i_E03_forcing_policy_F118_current_source_authority_supersession_evidence.json"
)
F118_PUBLICATION_AUDIT_RELATIVE = Path(f"{F118_RELATIVE}.publication_audit.json")
F118_PROVENANCE_REVIEW_RELATIVE = Path(f"{F118_RELATIVE}.provenance_security_review.json")
F118_PLASMA_REVIEW_RELATIVE = Path(f"{F118_RELATIVE}.plasma_scientific_review.json")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
REVISION_PATTERN = re.compile(r"[0-9a-f]{40}")
EXIT_CODE_PATTERN = re.compile(r"[0-9]+:[0-9]+")
ARTIFACT_NAME_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,199}\.json")
V2_RECOST_ARTIFACT_PATTERN = re.compile(
    r"mks24_stage_i_E03_forcing_policy_F([0-9]+)_recost_evidence\.json"
)
SELF_DESCRIPTOR_ENV = "_CGL_LF_RECOST_UTILITY_DESCRIPTOR"
PYTHON_DESCRIPTOR_ENV = "_CGL_LF_RECOST_UTILITY_PYTHON_DESCRIPTOR"
SELF_SOURCE_ENV = "_CGL_LF_RECOST_UTILITY_SOURCE"
ROOT_DIR_ENV = "_CGL_LF_RECOST_REPOSITORY_ROOT"
STAGE_I_SELF_DESCRIPTOR_ENV = "_CGL_LF_STAGE_I_CONTROLLER_DESCRIPTOR"
STAGE_I_PYTHON_DESCRIPTOR_ENV = "_CGL_LF_STAGE_I_CONTROLLER_PYTHON_DESCRIPTOR"
STAGE_I_SOURCE_ENV = "_CGL_LF_STAGE_I_CONTROLLER_SOURCE"
STAGE_I_REPOSITORY_ROOT_ENV = "_CGL_LF_STAGE_I_CONTROLLER_REPOSITORY_ROOT"
_AUTHENTICATED_SOURCE_PATH: Path | None = None
F116_REQUIRED_TOOLS = {
    "scripts/frontier/cgl_lf_stage_i.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_checkpoint.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_qualification.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_recost.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_source_authority.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_validate_segment.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_wave_plan.py": "0644",
}
F116_PUBLISHER_RELATIVE = "scripts/frontier/cgl_lf_stage_i_source_authority.py"
F118_REQUIRED_TOOLS = F116_REQUIRED_TOOLS
F118_PUBLISHER_RELATIVE = F116_PUBLISHER_RELATIVE
F118_AUTHORIZATION = {
    "current_source_selection_authorized": True,
    "source_authority_publication_authorized": True,
    "prepare_authorized": False,
    "submit_authorized": False,
    "direct_sbatch_authorized": False,
    "scheduler_mutation_authorized": False,
    "stage_i_execution_state_mutation_authorized": False,
    "scientific_configuration_change_authorized": False,
    "historical_manifest_rebinding_authorized": False,
}
F118_PRESERVES = [
    "The immutable F-116 evidence, reviews, publication audit, selected source bundle, and nested F-115 historical authority.",
    "The qualified executable, frozen source revision, inputs, matrix, restart lineages, targets, resources, qualification, and Stage I budget policy.",
    "Every prior active source-archive checksum-ledger entry and the corrupt-C7 incident-evidence exclusion.",
]
F118_DOES_NOT_AUTHORIZE = [
    "prepare",
    "submit",
    "direct sbatch",
    "scheduler mutation",
    "Stage I execution-state mutation",
    "scientific configuration change",
    "historical manifest rebinding",
]
F118_PUBLICATION_REQUIREMENTS = {
    "published_evidence_mode": "0444",
    "published_review_mode": "0444",
    "published_audit_mode": "0444",
    "published_links": 1,
    "publication_audit_is_authority_commit_marker": True,
    "recovery_required_after_interruption": True,
}
F118_VALIDATION_CLAIMS = {
    "historical_f116_chain": "passed",
    "historical_f115_chain": "passed",
    "bridge_bundle_complete_history": "passed",
    "predecessor_current_source_bundle_complete_history": "passed",
    "final_bundle_complete_history": "passed",
    "final_bundle_single_head_tip": "passed",
    "final_bundle_required_revisions": "passed",
    "committed_tool_bytes": "passed",
    "corrupt_c7_exclusion_preserved": True,
}
LEGACY_F114_RECOST_NAME = (
    "mks24_stage_i_E03_forcing_policy_R03_s00_clean_partial_recost_evidence.json"
)
LEGACY_F114_RECOST_SHA256 = (
    "27dd154b8594693afe0308f4de63a7faaed6daf81a396ab31166bb70a8c43b86"
)
LEGACY_F114_PUBLICATION_AUDIT_SHA256 = (
    "3bf50ef1359f7f1798da945185d42b61488cb98779d5a6942839c03e3eb6bc73"
)
F117_FAILED_ATTEMPT_RELATIVES = {
    "packet": Path(
        "accounting/mks24_stage_i_E03_forcing_policy_F117_recost_draft_packet.json"
    ),
    "request": Path(
        "accounting/mks24_stage_i_E03_forcing_policy_F117_recost_request.json"
    ),
    "reconciliation": Path(
        "accounting/mks24_stage_i_E03_forcing_policy_F117_reconciliation_evidence.json"
    ),
    "storage": Path(
        "accounting/mks24_stage_i_E03_forcing_policy_F117_storage_evidence.json"
    ),
}
F117_FAILED_ATTEMPT_SHA256 = {
    "packet": "8abceb6b2e3031a21b87f95b19053f2bf2b35f30d642c7a86bc4d9c6614071ad",
    "request": "5dea85b34b4f09cb3a926927ad73b4443e2f8ce4fdd2b2a1b1f84b70294bc632",
    "reconciliation": "677547889992d700e05559e322d95e7df6078a736bb2fbdff27ada38410c83e6",
    "storage": "637e43a3cd12f3051e1eae1d3102a2106d0445e26ac500b8946f8fea2064f218",
}
F119_ARTIFACT_NAME = "mks24_stage_i_E03_forcing_policy_F119_recost_evidence.json"
F119_LEGACY_BOOTSTRAP = (
    "exact-retained-legacy-F114-after-authenticated-F117-failed-attempt"
)
STORAGE_PROJECTION_METHOD = "observed-stage-i-output-byte-rate-v1"
INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION = (
    "Reviewer roles, agent identifiers, and process separation are retained "
    "declarations; exact artifact digests authenticate reviewed bytes but do not "
    "cryptographically authenticate a human or agent identity."
)
RENAME_NOREPLACE = 1
RENAME_EXCHANGE = 2
AT_FDCWD = -100
AT_SYMLINK_FOLLOW = 0x400
RENAMEAT2_UNSUPPORTED_ERRNOS = frozenset(
    {
        errno.EINVAL,
        errno.ENOSYS,
        getattr(errno, "ENOTSUP", errno.EINVAL),
        getattr(errno, "EOPNOTSUPP", errno.EINVAL),
    }
)
_ACTIVE_MUTATION_LOCK = None
_ACTIVE_PUBLIC_ROOT = None
_BOUND_MUTATION_PARENTS: dict[int, "BoundParentGuard"] = {}

# Same-user interposition inside one raw syscall cannot be prevented in user
# space.  Mutations below therefore issue that syscall once, durably classify
# the observable result, and never compensate after lock or namespace loss.

class NoRecostJournal(ValueError):
    """Report an authenticated empty recost transaction store."""

FRESH_R12_RERUN = {
    "case_id": "R12",
    "segment": "s01_rankio_t0_t0p12",
    "nodes": 4,
    "time_tlim_target": 0.12,
}
HISTORICAL_R12_CLEAN_PARTIAL = {
    "case_id": "R12",
    "segment": "s00_rankio_t0_t0p25",
    "job_id": "4766856",
    "result": "clean_partial",
}
COUNT_KEYS = (
    "transactions",
    "reservations",
    "active_reservations",
    "ledger_rows",
    "manifests",
)
CGL_JOB_NAME_PREFIX = "cgl_"
TERMINAL_SCHEDULER_STATES = frozenset(
    {
        "BOOT_FAIL",
        "CANCELLED",
        "COMPLETED",
        "DEADLINE",
        "FAILED",
        "NODE_FAIL",
        "OUT_OF_MEMORY",
        "PREEMPTED",
        "REVOKED",
        "SPECIAL_EXIT",
        "TIMEOUT",
    }
)
SAFE_ID_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,199}")
LOWER_RESOLUTION_CASE_PATTERN = re.compile(r"R(?:0[2-9]|1[0-6])")
MAX_BOUNDED_WAVE_LANES = 4
MAX_BOUNDED_WAVE_NODES = 10
MAX_STAGE_I_TIME = 10.0
MAX_WALLTIME_SECONDS = 2 * 60 * 60
SCHEDULER_TIME_TOLERANCE = timedelta(minutes=5)
V2_AUTHORIZATION_MAX_LIFETIME = timedelta(hours=24)
BOUNDED_PROFILE_KEYS = frozenset(
    {
        "acceptance_criterion",
        "acceptance_policy",
        "athena_walltime",
        "build_manifest",
        "build_manifest_sha256",
        "case_id",
        "controller_walltime_max_seconds",
        "cpus_per_task",
        "estimated_storage_bytes",
        "executable",
        "executable_revision",
        "executable_sha256",
        "input_file",
        "input_revision",
        "input_sha256",
        "nodes",
        "output_layout",
        "segment",
        "parent_job_id",
        "parent_result",
        "parent_segment",
        "restart_file",
        "restart_file_sha256",
        "restart_time",
        "ranks_per_node",
        "time_tlim_target",
        "walltime",
        "source_bundle",
        "source_bundle_sha256",
    }
)
SOLE_COMPATIBLE_PROFILE_KEYS = (
    "athena_walltime",
    "case_id",
    "nodes",
    "parent_job_id",
    "parent_result",
    "parent_segment",
    "restart_file",
    "restart_time",
    "segment",
    "source_bundle",
    "source_bundle_sha256",
    "time_tlim_target",
    "walltime",
)
SCHEDULER_EVIDENCE_KEYS = frozenset(
    {
        "path",
        "sha256",
        "job_id",
        "job_name",
        "state",
        "exit_code",
        "nodes",
        "elapsed_seconds",
        "submitted_utc",
        "completed_utc",
    }
)
V2_ARTIFACT_KEYS = frozenset(
    {
        "schema_version",
        "record_type",
        "checkpoint",
        "artifact_name",
        "execution_epoch",
        "generated_utc",
        "expires_utc",
        "scope",
        "predecessor_recost",
        "authorization",
        "barrier",
        "budget",
        "storage",
        "ledger",
        "reservations",
        "manifests",
        "r17_readiness",
        "promoted_f113",
        "reconcile",
        "provenance",
    }
)
SCHEMA2_RECOST_ARTIFACT_KEYS = frozenset(
    {
        "schema_version",
        "record_type",
        "checkpoint",
        "artifact_name",
        "execution_epoch",
        "generated_utc",
        "expires_utc",
        "requested_by",
        "scope",
        "predecessor_recost",
        "authority",
        "publication_requirements",
        "recommendations",
        "barrier",
        "budget",
        "storage",
        "ledger",
        "reservations",
        "manifests",
        "r17_readiness",
        "promoted_f113",
        "reconcile",
        "provenance",
    }
)
V2_PROVENANCE_KEYS = frozenset(
    {
        "request_sha256",
        "generator_sha256",
        "generator_revision",
        "stage_i_helper_sha256",
        "stage_i_helper_revision",
        "matrix_sha256",
        "matrix_revision",
        "source_bundle_sha256",
        "source_bundle_verified_revisions",
        "ceiling_evidence_sha256",
        "ceiling_publication_audit_sha256",
        "storage_evidence_sha256",
        "reconciliation_sha256",
        "ledger_sha256",
        "reservations_sha256",
        "scheduler_evidence",
        "predecessor_recost_sha256",
        "predecessor_recost_publication_audit_sha256",
        "authenticated_lineages_sha256",
        "computed_projection_sha256",
        "r17_readiness_evidence_sha256",
    }
)
SCHEMA2_RECOST_PROVENANCE_KEYS = frozenset(
    {
        *V2_PROVENANCE_KEYS,
        "request_independent_review_sha256",
        "source_authority",
        "qualification_approval_sha256",
        "f113_historical_helper_revision",
        "f113_historical_helper_sha256",
        "predecessor_recost_independent_review_sha256",
    }
)


class LinkAttemptError(RuntimeError):
    """Report an os.link attempt whose post-state requires recovery review."""


def utc_now() -> str:
    """Return one stable UTC timestamp."""

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def require_utc_timestamp(value: object, label: str) -> None:
    """Require one UTC ISO-8601 timestamp."""

    if not isinstance(value, str):
        raise ValueError(f"{label} must be an ISO-8601 timestamp")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} must be an ISO-8601 timestamp") from error
    if parsed.tzinfo is None:
        raise ValueError(f"{label} must include a timezone")
    if parsed.utcoffset() != timedelta(0):
        raise ValueError(f"{label} must use UTC")


def parse_utc_timestamp(value: object, label: str) -> datetime:
    """Parse and return one required UTC ISO-8601 timestamp."""

    require_utc_timestamp(value, label)
    assert isinstance(value, str)
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


def sha256_bytes(value: bytes) -> str:
    """Return the SHA-256 digest for retained bytes."""

    return hashlib.sha256(value).hexdigest()


def read_descriptor_bytes(descriptor: int) -> bytes:
    """Read one descriptor without changing its caller-visible offset."""

    offset = os.lseek(descriptor, 0, os.SEEK_CUR)
    try:
        os.lseek(descriptor, 0, os.SEEK_SET)
        blocks = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                return b"".join(blocks)
            blocks.append(block)
    finally:
        os.lseek(descriptor, offset, os.SEEK_SET)


def write_descriptor_bytes(descriptor: int, value: bytes) -> None:
    """Write all selected bytes to one descriptor."""

    view = memoryview(value)
    while view:
        written = os.write(descriptor, view)
        if written <= 0:
            raise OSError("descriptor write made no progress")
        view = view[written:]


def sha256_descriptor(descriptor: int) -> str:
    """Return the digest of one open descriptor."""

    return sha256_bytes(read_descriptor_bytes(descriptor))


def sha256(path: Path) -> str:
    """Return the digest of one regular retained file."""

    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        return sha256_descriptor(descriptor)
    finally:
        os.close(descriptor)


def require_sha256(value: object, label: str) -> str:
    """Require a lowercase SHA-256 digest."""

    if not isinstance(value, str) or SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return value


def sha256_arg(value: str) -> str:
    """Parse a SHA-256 command argument."""

    try:
        return require_sha256(value, "value")
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def revision_arg(value: str) -> str:
    """Parse one exact lowercase Git revision."""

    if REVISION_PATTERN.fullmatch(value) is None:
        raise argparse.ArgumentTypeError("revision must be a lowercase 40-hex commit ID")
    return value


def mode_arg(value: str) -> int:
    """Parse a four-digit retained-file mode."""

    if re.fullmatch(r"0?[0-7]{3}", value) is None:
        raise argparse.ArgumentTypeError("mode must be an octal value such as 0644")
    return int(value, 8)


def capped_mode_arg(value: str, maximum: int, label: str) -> int:
    """Parse a retained-file mode with no permissions beyond one ceiling."""

    mode = mode_arg(value)
    if mode & ~maximum:
        raise argparse.ArgumentTypeError(
            f"{label} mode must not exceed {maximum:04o}"
        )
    return mode


def artifact_mode_arg(value: str) -> int:
    """Parse a retained recost-artifact mode."""

    return capped_mode_arg(value, 0o644, "artifact")


def generator_mode_arg(value: str) -> int:
    """Parse a retained recost-generator mode."""

    return capped_mode_arg(value, 0o755, "generator")


def scheduler_mode_arg(value: str) -> int:
    """Parse a retained scheduler-evidence mode."""

    return capped_mode_arg(value, 0o644, "scheduler evidence")


def source_bundle_mode_arg(value: str) -> int:
    """Parse a retained source-bundle mode."""

    return capped_mode_arg(value, 0o644, "source bundle")


def json_object_arg(value: str) -> dict[str, object]:
    """Parse one explicit JSON object argument."""

    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as error:
        raise argparse.ArgumentTypeError("value must be a JSON object") from error
    if not isinstance(parsed, dict) or not parsed:
        raise argparse.ArgumentTypeError("value must be a nonempty JSON object")
    return parsed


def json_array_arg(value: str) -> list[object]:
    """Parse one explicit nonempty JSON array argument."""

    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as error:
        raise argparse.ArgumentTypeError("value must be a JSON array") from error
    if not isinstance(parsed, list) or not parsed:
        raise argparse.ArgumentTypeError("value must be a nonempty JSON array")
    return parsed


def revision_array_arg(value: str) -> list[str]:
    """Parse one nonempty JSON array of unique lowercase Git revisions."""

    parsed = json_array_arg(value)
    if (
        any(not isinstance(item, str) or REVISION_PATTERN.fullmatch(item) is None
            for item in parsed)
        or len(set(parsed)) != len(parsed)
    ):
        raise argparse.ArgumentTypeError(
            "value must be a nonempty JSON array of unique lowercase 40-hex revisions"
        )
    return parsed


def nonnegative_int_arg(value: str) -> int:
    """Parse one non-negative retained count."""

    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("count must be a non-negative integer") from error
    if parsed < 0:
        raise argparse.ArgumentTypeError("count must be a non-negative integer")
    return parsed


def expected_self_sha256(argv: list[str]) -> str:
    """Read the externally supplied self digest before general CLI parsing."""

    option = "--expected-utility-sha256"
    matches = [
        argv[index + 1]
        for index, value in enumerate(argv[:-1])
        if value == option
    ]
    if len(matches) != 1:
        raise ValueError(f"{option} must be supplied exactly once")
    return require_sha256(matches[0], option)


def initial_source_path() -> Path:
    """Return the named source path for initial descriptor re-execution."""

    if _AUTHENTICATED_SOURCE_PATH is not None:
        return _AUTHENTICATED_SOURCE_PATH
    source = Path(__file__).expanduser().absolute()
    if source.is_symlink():
        raise ValueError("utility source must not be a symbolic link")
    return source.resolve()


def repository_root(source: Path) -> Path:
    """Return the named repository root for initial descriptor re-execution."""

    try:
        return source.parents[2]
    except IndexError as error:
        raise ValueError(f"utility source path is invalid: {source}") from error


def require_reexec_source_relationship(source: Path, root: Path) -> None:
    """Bind private reexec metadata to the utility's repository-relative path."""

    if source != root / UTILITY_RELATIVE:
        raise ValueError("authenticated utility source/repository relationship differs")


def inherited_reexec_path(name: str, label: str, *, directory: bool = False) -> Path:
    """Read one private reexec path only after authenticating the self descriptor."""

    retained = os.environ.get(name)
    if retained is None:
        raise ValueError(f"authenticated utility reexecution lacks {label}")
    path = Path(retained)
    if (
        not path.is_absolute()
        or path != Path(os.path.normpath(retained))
        or ".." in path.parts
        or len(path.parts) < 2
    ):
        raise ValueError(f"authenticated utility {label} is not normalized and absolute")
    flags = os.O_RDONLY | (os.O_DIRECTORY if directory else 0)
    with absolute_descriptor(path, f"authenticated utility {label}", flags=flags):
        pass
    return path


def reexec_environment(descriptor: int, python_descriptor: int,
                       source: Path, root: Path) -> dict[str, str]:
    """Return the complete caller-independent authenticated reexec environment."""

    return {
        SELF_DESCRIPTOR_ENV: str(descriptor),
        PYTHON_DESCRIPTOR_ENV: str(python_descriptor),
        SELF_SOURCE_ENV: str(source),
        ROOT_DIR_ENV: str(root),
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": TRUSTED_SYSTEM_PATH,
        "PYTHONDONTWRITEBYTECODE": "1",
        "XDG_CONFIG_HOME": "/nonexistent",
    }


def profile_identity(profile: os.stat_result) -> tuple[int, int]:
    """Return one stable filesystem object identity."""

    return profile.st_dev, profile.st_ino


def profile_security_binding(profile: os.stat_result) -> tuple[int, ...]:
    """Return the mutation-stable security profile for one filesystem object."""

    return (
        profile.st_dev,
        profile.st_ino,
        stat.S_IFMT(profile.st_mode),
        stat.S_IMODE(profile.st_mode),
        profile.st_nlink,
        profile.st_uid,
        profile.st_gid,
        profile.st_size,
    )


def require_authenticated_python_descriptor() -> int:
    """Authenticate the inherited interpreter descriptor against this process."""

    retained = os.environ.get(PYTHON_DESCRIPTOR_ENV)
    if retained is None or re.fullmatch(r"[0-9]+", retained) is None:
        raise ValueError("authenticated utility Python descriptor marker is invalid")
    descriptor = int(retained)
    profile = os.fstat(descriptor)
    require_system_executable_profile(profile, "authenticated utility Python interpreter")
    if profile_identity(profile) != profile_identity(os.stat("/proc/self/exe")):
        raise ValueError("authenticated utility Python descriptor is not this interpreter")
    if not sys.flags.isolated:
        raise ValueError("authenticated utility Python interpreter is not isolated")
    return descriptor


def authenticate_retained_self_source(source: Path, expected: str) -> None:
    """Require the retained source profile even after descriptor re-execution."""

    with absolute_descriptor(source, "retained utility", flags=os.O_RDONLY) as descriptor:
        require_regular_profile(
            os.fstat(descriptor),
            source,
            "retained utility",
            expected_mode=0o755,
            expected_links=1,
        )
        if sha256_descriptor(descriptor) != expected:
            raise ValueError("retained utility checksum has changed")


def authenticate_self(expected: str) -> tuple[Path, Path]:
    """Authenticate this utility and re-execute immutable descriptor bytes."""

    global _AUTHENTICATED_SOURCE_PATH
    expected = require_sha256(expected, "utility SHA-256")
    inherited = os.environ.get(SELF_DESCRIPTOR_ENV)
    if inherited is not None:
        require_authenticated_python_descriptor()
        if re.fullmatch(r"[0-9]+", inherited) is None:
            raise ValueError("utility descriptor marker is invalid")
        try:
            descriptor = int(inherited)
        except ValueError as error:
            raise ValueError("utility descriptor marker is invalid") from error
        if __file__ != f"/proc/self/fd/{descriptor}":
            raise ValueError("utility descriptor marker is not attached to this execution")
        require_regular_profile(
            os.fstat(descriptor), Path(__file__), "authenticated utility descriptor"
        )
        if sha256_descriptor(descriptor) != expected:
            raise ValueError("authenticated utility descriptor checksum has changed")
        source = inherited_reexec_path(SELF_SOURCE_ENV, "source path")
        root = inherited_reexec_path(ROOT_DIR_ENV, "repository root", directory=True)
        require_reexec_source_relationship(source, root)
        authenticate_retained_self_source(source, expected)
        _AUTHENTICATED_SOURCE_PATH = source
        return source, root

    orphaned = [
        name for name in (PYTHON_DESCRIPTOR_ENV, SELF_SOURCE_ENV, ROOT_DIR_ENV)
        if name in os.environ
    ]
    if orphaned:
        raise ValueError(
            "private utility reexecution path is forbidden without an authenticated "
            "descriptor"
        )
    source = initial_source_path()
    root = repository_root(source)
    require_reexec_source_relationship(source, root)
    authenticate_retained_self_source(source, expected)
    value = require_file_sha256(
        source,
        expected,
        "retained utility",
        expected_mode=0o755,
        expected_links=1,
    )
    if not hasattr(os, "memfd_create"):
        raise ValueError("descriptor re-execution requires os.memfd_create")
    retained = os.memfd_create("cgl_lf_stage_i_checkpoint")
    os.write(retained, value)
    os.lseek(retained, 0, os.SEEK_SET)
    interpreter = Path(sys.executable).resolve(strict=True)
    with absolute_descriptor(
        interpreter, "authenticated utility Python interpreter", flags=os.O_RDONLY
    ) as python_descriptor:
        require_system_executable_profile(
            os.fstat(python_descriptor), "authenticated utility Python interpreter"
        )
        if profile_identity(os.fstat(python_descriptor)) != profile_identity(
            os.stat("/proc/self/exe")
        ):
            raise ValueError("selected utility Python interpreter is not this process")
        os.set_inheritable(retained, True)
        os.set_inheritable(python_descriptor, True)
        os.execve(
            f"/proc/self/fd/{python_descriptor}",
            [
                str(interpreter),
                "-I",
                "-B",
                f"/proc/self/fd/{retained}",
                *sys.argv[1:],
            ],
            reexec_environment(retained, python_descriptor, source, root),
        )
    raise AssertionError("descriptor re-execution returned unexpectedly")


def fsync_directory(path: Path) -> None:
    """Persist directory-entry changes."""

    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def require_entry_name(name: str, label: str) -> str:
    """Require one direct-child basename for descriptor-relative mutation."""

    if not name or name in {".", ".."} or Path(name).name != name or "/" in name:
        raise ValueError(f"{label} is not a direct-child entry name")
    return name


def require_bound_entry_identity(parent: int, name: str, expected: os.stat_result,
                                 label: str) -> os.stat_result:
    """Require one directory entry to retain an authenticated inode identity."""

    name = require_entry_name(name, label)
    try:
        named = os.stat(name, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError as error:
        raise ValueError(f"{label} disappeared before mutation") from error
    if profile_identity(named) != profile_identity(expected):
        raise ValueError(f"{label} inode identity changed before mutation")
    return named


def require_bound_entry_profile(parent: int, name: str, expected: os.stat_result,
                                label: str) -> os.stat_result:
    """Require one entry to retain its complete mutation-stable profile."""

    named = require_bound_entry_identity(parent, name, expected, label)
    if profile_security_binding(named) != profile_security_binding(expected):
        raise ValueError(f"{label} security profile changed during mutation")
    return named


def require_bound_entry_security(
    parent: int,
    name: str,
    expected: os.stat_result,
    label: str,
    *,
    expected_sha256: str | None = None,
    allow_mode_zero: bool = False,
) -> str | None:
    """Authenticate one bound regular entry's profile and readable content."""

    named = require_bound_entry_profile(parent, name, expected, label)
    require_regular_profile(named, Path(name), label)
    if stat.S_IMODE(named.st_mode) == 0:
        if not allow_mode_zero:
            raise ValueError(f"{label} mode-0000 content cannot be authenticated")
        if expected_sha256 is not None:
            raise ValueError(f"{label} mode-0000 content cannot satisfy a digest binding")
        return None
    try:
        descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent)
    except OSError as error:
        raise ValueError(f"{label} cannot be opened for content authentication") from error
    try:
        opened = os.fstat(descriptor)
        if profile_security_binding(opened) != profile_security_binding(expected):
            raise ValueError(f"{label} descriptor security profile changed during mutation")
        digest = sha256_descriptor(descriptor)
        if profile_security_binding(os.fstat(descriptor)) != profile_security_binding(
            expected
        ):
            raise ValueError(f"{label} descriptor security profile changed during mutation")
    finally:
        os.close(descriptor)
    require_bound_entry_profile(parent, name, expected, label)
    if expected_sha256 is not None and digest != require_sha256(
        expected_sha256, f"{label} SHA-256"
    ):
        raise ValueError(f"{label} content digest changed during mutation")
    return digest


def require_bound_descriptor_security(
    parent: int,
    name: str,
    descriptor: int,
    expected: os.stat_result,
    label: str,
    *,
    expected_sha256: str,
) -> None:
    """Authenticate one open descriptor and its selected public name."""

    expected_sha256 = require_sha256(expected_sha256, f"{label} SHA-256")
    require_bound_entry_profile(parent, name, expected, label)
    opened = os.fstat(descriptor)
    require_regular_profile(opened, Path(name), label)
    if profile_security_binding(opened) != profile_security_binding(expected):
        raise ValueError(f"{label} descriptor security profile changed during mutation")
    if sha256_descriptor(descriptor) != expected_sha256:
        raise ValueError(f"{label} content digest changed during mutation")
    if profile_security_binding(os.fstat(descriptor)) != profile_security_binding(
        expected
    ):
        raise ValueError(f"{label} descriptor security profile changed during mutation")
    require_bound_entry_profile(parent, name, expected, label)


class MutationLockGuard:
    """Continuously rebind mutation boundaries to the selected lock inode."""

    def __init__(self, parent: int, name: str, descriptor: int,
                 profile: os.stat_result, label: str):
        self.parent = parent
        self.name = require_entry_name(name, label)
        self.descriptor = descriptor
        self.identity = profile_identity(profile)
        self.label = label

    def assert_bound(self) -> None:
        opened = os.fstat(self.descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or profile_identity(opened) != self.identity
            or opened.st_nlink != 1
        ):
            raise ValueError(f"{self.label} descriptor changed while mutation is active")
        try:
            named = os.stat(self.name, dir_fd=self.parent, follow_symlinks=False)
        except FileNotFoundError as error:
            raise ValueError(f"{self.label} path changed while mutation is active") from error
        if (
            not stat.S_ISREG(named.st_mode)
            or profile_identity(named) != self.identity
            or named.st_nlink != 1
        ):
            raise ValueError(f"{self.label} path changed while mutation is active")


class BoundParentGuard:
    """Continuously bind one mutation descriptor to its trusted public path."""

    def __init__(self, descriptor: int, path: Path, profile: os.stat_result, label: str):
        self.descriptor = descriptor
        self.path = absolute_path(path)
        self.label = label
        self.boundary = (
            *profile_identity(profile),
            stat.S_IMODE(profile.st_mode),
            profile.st_uid,
            profile.st_gid,
        )
        self.ancestors = mutation_ancestor_boundaries(self.path, self.label)

    def assert_bound(self) -> None:
        require_mutation_ancestor_boundaries(self.ancestors, self.label)
        opened = os.fstat(self.descriptor)
        boundary = (
            *profile_identity(opened),
            stat.S_IMODE(opened.st_mode),
            opened.st_uid,
            opened.st_gid,
        )
        if boundary != self.boundary:
            raise ValueError(f"{self.label} parent path changed during mutation")
        require_trusted_directory_profile(opened, self.path, self.label)
        with absolute_descriptor(
            self.path, f"{self.label} parent binding", flags=os.O_RDONLY | os.O_DIRECTORY
        ) as named:
            named_profile = os.fstat(named)
            named_boundary = (
                *profile_identity(named_profile),
                stat.S_IMODE(named_profile.st_mode),
                named_profile.st_uid,
                named_profile.st_gid,
            )
        if named_boundary != self.boundary:
            raise ValueError(f"{self.label} parent path changed during mutation")
        require_trusted_directory_profile(named_profile, self.path, self.label)


def directory_boundary(profile: os.stat_result) -> tuple[int, ...]:
    """Return the identity and replacement-relevant profile of one directory."""

    return (
        *profile_identity(profile),
        stat.S_IMODE(profile.st_mode),
        profile.st_uid,
        profile.st_gid,
    )


def require_mutation_ancestor_profile(value: os.stat_result, path: Path,
                                      label: str) -> None:
    """Apply the leaf trust profile below an authenticated public boundary."""

    active_root = _ACTIVE_PUBLIC_ROOT.path if _ACTIVE_PUBLIC_ROOT is not None else None
    offline_boundary = next(
        (
            absolute_path(prefix)
            for prefix in OFFLINE_ALLOWED_PREFIXES
            if path == absolute_path(prefix) or absolute_path(prefix) in path.parents
        ),
        None,
    )
    if active_root is not None and (path == active_root or active_root in path.parents):
        require_trusted_directory_profile(value, path, label)
    elif offline_boundary is not None and path != offline_boundary:
        require_trusted_directory_profile(value, path, label)
    else:
        require_trusted_ancestor_profile(value, path, label)


def mutation_ancestor_boundaries(path: Path, label: str
                                 ) -> list[tuple[Path, tuple[int, ...]]]:
    """Capture every ancestor selecting one mutation directory."""

    retained = []
    for ancestor in reversed(path.parents):
        if ancestor == Path("/"):
            profile = os.stat("/", follow_symlinks=False)
        else:
            with absolute_descriptor(
                ancestor,
                f"{label} ancestor",
                flags=os.O_RDONLY | os.O_DIRECTORY,
            ) as descriptor:
                profile = os.fstat(descriptor)
        require_mutation_ancestor_profile(profile, ancestor, f"{label} ancestor")
        retained.append((ancestor, directory_boundary(profile)))
    return retained


def require_mutation_ancestor_boundaries(
    retained: list[tuple[Path, tuple[int, ...]]], label: str
) -> None:
    """Revalidate every exact ancestor selecting one mutation directory."""

    for ancestor, expected in retained:
        if ancestor == Path("/"):
            profile = os.stat("/", follow_symlinks=False)
        else:
            with absolute_descriptor(
                ancestor,
                f"{label} ancestor binding",
                flags=os.O_RDONLY | os.O_DIRECTORY,
            ) as descriptor:
                profile = os.fstat(descriptor)
        if directory_boundary(profile) != expected:
            raise ValueError(f"{label} ancestor path changed during mutation: {ancestor}")
        require_mutation_ancestor_profile(profile, ancestor, f"{label} ancestor")


def require_trusted_ancestor_profile(value: os.stat_result, path: Path, label: str) -> None:
    """Accept only an exact safe ancestor or the deployed trusted-project boundary."""

    if not stat.S_ISDIR(value.st_mode):
        raise ValueError(f"{label} is not a directory: {path}")
    mode = stat.S_IMODE(value.st_mode)
    offline_sticky_boundary = (
        any(path == absolute_path(prefix) for prefix in OFFLINE_ALLOWED_PREFIXES)
        and mode & stat.S_ISVTX
    )
    if mode & stat.S_IWOTH:
        if not offline_sticky_boundary:
            raise ValueError(f"{label} is world-writable and untrusted: {path}")
    if mode & stat.S_IWGRP and not offline_sticky_boundary:
        # Orion's canonical project boundary is intentionally group-writable;
        # bind its exact deployed path, owner, group, and mode instead of
        # treating all group writability as equivalent.
        canonical_project_boundary = (
            path == DEFAULT_ROOT.parents[1]
            and value.st_uid == 0
            and value.st_gid == 31114
            and mode == 0o2770
        )
        if not canonical_project_boundary:
            raise ValueError(f"{label} is group-writable and untrusted: {path}")


class BoundPublicRootGuard:
    """Hold and revalidate every component selecting the public root."""

    def __init__(self, path: Path, components: list[tuple[Path, int, os.stat_result]]):
        self.path = absolute_path(path)
        self.components = [
            (component, descriptor, directory_boundary(profile))
            for component, descriptor, profile in components
        ]

    @property
    def descriptor(self) -> int:
        return self.components[-1][1]

    def assert_bound(self) -> None:
        for component, descriptor, expected in self.components:
            opened = os.fstat(descriptor)
            if directory_boundary(opened) != expected:
                raise ValueError(f"Stage I public namespace changed: {component}")
            if component == self.path:
                require_trusted_directory_profile(opened, component, "Stage I root")
            else:
                require_trusted_ancestor_profile(
                    opened, component, "Stage I root ancestor"
                )
            if component == Path("/"):
                named_boundary = directory_boundary(os.stat("/", follow_symlinks=False))
            else:
                with absolute_descriptor(
                    component,
                    "Stage I public namespace binding",
                    flags=os.O_RDONLY | os.O_DIRECTORY,
                ) as named:
                    named_boundary = directory_boundary(os.fstat(named))
            if named_boundary != expected:
                raise ValueError(f"Stage I public namespace changed: {component}")


@contextmanager
def bound_public_root(path: Path):
    """Bind every no-symlink component selecting one public root."""

    path = absolute_path(path)
    descriptor = os.open("/", os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    opened = [descriptor]
    components = [(Path("/"), descriptor, os.fstat(descriptor))]
    current = Path("/")
    try:
        for name in path.parts[1:]:
            descriptor = os.open(
                name,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=descriptor,
            )
            opened.append(descriptor)
            current = current / name
            profile = os.fstat(descriptor)
            if current == path:
                require_trusted_directory_profile(profile, current, "Stage I root")
            else:
                require_trusted_ancestor_profile(
                    profile, current, "Stage I root ancestor"
                )
            components.append((current, descriptor, profile))
        guard = BoundPublicRootGuard(path, components)
        guard.assert_bound()
        yield guard
        guard.assert_bound()
    finally:
        for retained in reversed(opened):
            os.close(retained)


def require_active_mutation_lock_bound() -> None:
    """Fail before and after irreversible writes if the active lock moved."""

    if _ACTIVE_PUBLIC_ROOT is not None:
        _ACTIVE_PUBLIC_ROOT.assert_bound()
    if _ACTIVE_MUTATION_LOCK is not None:
        _ACTIVE_MUTATION_LOCK.assert_bound()


def require_mutation_authority_bound(parent: int) -> None:
    """Require the active lock and selected mutation parent to remain authorized."""

    require_active_mutation_lock_bound()
    guards = list(_BOUND_MUTATION_PARENTS.values())
    for guard in guards:
        guard.assert_bound()
    if parent not in _BOUND_MUTATION_PARENTS:
        require_trusted_directory_profile(
            os.fstat(parent), Path(f"/proc/self/fd/{parent}"), "mutation parent"
        )


def renameat2_between(source_parent: int, source: str, target_parent: int,
                      target: str, flags: int, label: str) -> None:
    """Perform one descriptor-relative Linux renameat2 operation."""

    source = require_entry_name(source, label)
    target = require_entry_name(target, f"{label} target")
    try:
        operation = ctypes.CDLL(None, use_errno=True).renameat2
    except AttributeError as error:
        raise ValueError("descriptor-bound publication requires renameat2") from error
    operation.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    operation.restype = ctypes.c_int
    if operation(
        source_parent,
        os.fsencode(source),
        target_parent,
        os.fsencode(target),
        flags,
    ) != 0:
        retained_errno = ctypes.get_errno()
        raise OSError(
            retained_errno, os.strerror(retained_errno), f"{source} -> {target}"
        )


def renameat2(parent: int, source: str, target: str, flags: int, label: str) -> None:
    """Perform and durably record one raw same-directory renameat2 attempt."""

    try:
        renameat2_between(parent, source, parent, target, flags, label)
    finally:
        os.fsync(parent)


def entry_security_matches(parent: int, name: str, expected: os.stat_result,
                           label: str, expected_sha256: str | None, *,
                           allow_mode_zero: bool = False) -> bool:
    """Return whether one name exactly retains an authenticated regular entry."""

    try:
        require_bound_entry_security(
            parent,
            name,
            expected,
            label,
            expected_sha256=expected_sha256,
            allow_mode_zero=allow_mode_zero,
        )
    except (FileNotFoundError, OSError, ValueError):
        return False
    return True


def entry_absent(parent: int, name: str) -> bool:
    """Return whether one direct-child name is absent without following links."""

    try:
        os.stat(name, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError:
        return True
    return False


def entry_identity_matches(parent: int, name: str, expected: os.stat_result) -> bool:
    """Return whether one direct child selects the expected inode identity."""

    try:
        profile = os.stat(name, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError:
        return False
    return profile_identity(profile) == profile_identity(expected)


def fsync_descriptors(*descriptors: int) -> BaseException | None:
    """Attempt every selected durability barrier and retain the first failure."""

    failure: BaseException | None = None
    for descriptor in descriptors:
        try:
            os.fsync(descriptor)
        except BaseException as error:
            if failure is None:
                failure = error
    return failure


def renameat2_is_unsupported(error: BaseException | None) -> bool:
    """Return whether one renameat2 failure permits the reviewed Lustre fallback."""

    return isinstance(error, OSError) and error.errno in RENAMEAT2_UNSUPPORTED_ERRNOS


def profile_security_binding_without_links(profile: os.stat_result) -> tuple[int, ...]:
    """Return the stable regular-file profile while a hard-link commit is active."""

    return (
        profile.st_dev,
        profile.st_ino,
        stat.S_IFMT(profile.st_mode),
        stat.S_IMODE(profile.st_mode),
        profile.st_uid,
        profile.st_gid,
        profile.st_size,
    )


def authenticate_descriptor_with_links(
    descriptor: int,
    expected: os.stat_result,
    expected_links: int,
    label: str,
    *,
    expected_sha256: str | None = None,
) -> tuple[os.stat_result, str]:
    """Authenticate exact readable bytes while permitting one selected link count."""

    before = os.fstat(descriptor)
    require_regular_profile(before, Path(f"/proc/self/fd/{descriptor}"), label)
    if (
        profile_security_binding_without_links(before)
        != profile_security_binding_without_links(expected)
        or before.st_nlink != expected_links
    ):
        raise ValueError(f"{label} security profile changed during hard-link commit")
    if stat.S_IMODE(before.st_mode) == 0:
        raise ValueError(f"{label} mode-0000 content cannot be exactly authenticated")
    digest = sha256_descriptor(descriptor)
    after = os.fstat(descriptor)
    if (
        profile_security_binding_without_links(after)
        != profile_security_binding_without_links(expected)
        or after.st_nlink != expected_links
    ):
        raise ValueError(f"{label} descriptor changed during hard-link authentication")
    if expected_sha256 is not None and digest != require_sha256(
        expected_sha256, f"{label} SHA-256"
    ):
        raise ValueError(f"{label} content digest changed during hard-link commit")
    return after, digest


def open_bound_entry_with_links(
    parent: int,
    name: str,
    expected: os.stat_result,
    expected_links: int,
    label: str,
    *,
    expected_sha256: str | None = None,
) -> tuple[int, os.stat_result, str]:
    """Open and authenticate one exact direct child during a hard-link commit."""

    name = require_entry_name(name, label)
    named = os.stat(name, dir_fd=parent, follow_symlinks=False)
    if (
        profile_security_binding_without_links(named)
        != profile_security_binding_without_links(expected)
        or named.st_nlink != expected_links
    ):
        raise ValueError(f"{label} name changed during hard-link commit")
    descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent)
    try:
        opened, digest = authenticate_descriptor_with_links(
            descriptor,
            expected,
            expected_links,
            label,
            expected_sha256=expected_sha256,
        )
        current = os.stat(name, dir_fd=parent, follow_symlinks=False)
        if (
            profile_security_binding_without_links(current)
            != profile_security_binding_without_links(opened)
            or current.st_nlink != expected_links
        ):
            raise ValueError(f"{label} name changed during hard-link authentication")
        return descriptor, opened, digest
    except BaseException:
        os.close(descriptor)
        raise


def classify_bound_hard_link_move(
    source_parent: int,
    source: str,
    target_parent: int,
    target: str,
    expected: os.stat_result,
    label: str,
    *,
    expected_sha256: str | None = None,
) -> tuple[str, str, int]:
    """Authenticate and classify every durable state of one hard-link move."""

    source = require_entry_name(source, label)
    target = require_entry_name(target, f"{label} target")
    try:
        source_profile = os.stat(source, dir_fd=source_parent, follow_symlinks=False)
    except FileNotFoundError:
        source_profile = None
    try:
        target_profile = os.stat(target, dir_fd=target_parent, follow_symlinks=False)
    except FileNotFoundError:
        target_profile = None
    for profile, selected_label in (
        (source_profile, label),
        (target_profile, f"{label} target"),
    ):
        if profile is not None and profile_security_binding_without_links(
            profile
        ) != profile_security_binding_without_links(expected):
            raise ValueError(f"{selected_label} differs during hard-link commit")
    if source_profile is not None and target_profile is None:
        state = "source-only"
        links = source_profile.st_nlink
        if expected.st_nlink != links:
            raise ValueError(f"{label} link count changed during hard-link commit")
    elif source_profile is not None and target_profile is not None:
        if profile_identity(source_profile) != profile_identity(target_profile):
            raise ValueError(f"{label} names select different inodes during hard-link commit")
        state = "linked"
        links = source_profile.st_nlink
        if (
            target_profile.st_nlink != links
            or links < 2
            or expected.st_nlink not in {links - 1, links}
        ):
            raise ValueError(f"{label} link count changed during hard-link commit")
    elif source_profile is None and target_profile is not None:
        state = "target-only"
        links = target_profile.st_nlink
        if expected.st_nlink not in {links, links + 1}:
            raise ValueError(f"{label} link count changed during hard-link commit")
    else:
        raise ValueError(f"{label} hard-link commit lost both selected names")
    digest = expected_sha256
    for parent, name, profile, selected_label in (
        (source_parent, source, source_profile, label),
        (target_parent, target, target_profile, f"{label} target"),
    ):
        if profile is None:
            continue
        descriptor, _, observed_digest = open_bound_entry_with_links(
            parent,
            name,
            expected,
            links,
            selected_label,
            expected_sha256=digest,
        )
        os.close(descriptor)
        if digest is None:
            digest = observed_digest
    assert digest is not None
    return state, digest, links


def linkat_descriptor_noreplace(
    source: int, target_parent: int, target: str, label: str
) -> None:
    """Create one descriptor-backed hard link without replacing a direct child."""

    target = require_entry_name(target, label)
    try:
        operation = ctypes.CDLL(None, use_errno=True).linkat
    except AttributeError as error:
        raise ValueError("descriptor-bound publication requires linkat") from error
    operation.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
    )
    operation.restype = ctypes.c_int
    if operation(
        AT_FDCWD,
        os.fsencode(f"/proc/self/fd/{source}"),
        target_parent,
        os.fsencode(target),
        AT_SYMLINK_FOLLOW,
    ) != 0:
        retained_errno = ctypes.get_errno()
        raise OSError(retained_errno, os.strerror(retained_errno), target)


def complete_bound_hard_link_move(
    source_parent: int,
    source: str,
    target_parent: int,
    target: str,
    expected: os.stat_result,
    label: str,
    *,
    expected_sha256: str | None,
) -> None:
    """Resume or complete one exact no-clobber hard-link move."""

    require_mutation_authority_bound(source_parent)
    require_mutation_authority_bound(target_parent)
    source = require_entry_name(source, label)
    target = require_entry_name(target, f"{label} target")
    state, digest, links = classify_bound_hard_link_move(
        source_parent,
        source,
        target_parent,
        target,
        expected,
        label,
        expected_sha256=expected_sha256,
    )
    if state == "target-only":
        require_mutation_authority_bound(source_parent)
        require_mutation_authority_bound(target_parent)
        return
    if state == "source-only":
        descriptor, _, retained_digest = open_bound_entry_with_links(
            source_parent,
            source,
            expected,
            links,
            label,
            expected_sha256=digest,
        )
        operation_error: BaseException | None = None
        try:
            require_mutation_authority_bound(source_parent)
            require_mutation_authority_bound(target_parent)
            linkat_descriptor_noreplace(descriptor, target_parent, target, label)
        except BaseException as error:
            operation_error = error
        finally:
            os.close(descriptor)
        durability_error = fsync_descriptors(source_parent, target_parent)
        try:
            state, digest, _ = classify_bound_hard_link_move(
                source_parent,
                source,
                target_parent,
                target,
                expected,
                label,
                expected_sha256=retained_digest,
            )
        except BaseException:
            require_mutation_authority_bound(source_parent)
            require_mutation_authority_bound(target_parent)
            if durability_error is not None:
                raise durability_error
            if isinstance(operation_error, FileExistsError):
                raise operation_error
            raise
        require_mutation_authority_bound(source_parent)
        require_mutation_authority_bound(target_parent)
        if durability_error is not None:
            raise durability_error
        if state == "source-only" and operation_error is not None:
            raise operation_error
        if state == "source-only":
            raise ValueError(f"{label} hard-link no-replace did not publish the target")
        if state == "target-only":
            return
    if state != "linked":
        raise ValueError(f"{label} hard-link commit entered an invalid state")
    require_mutation_authority_bound(source_parent)
    require_mutation_authority_bound(target_parent)
    operation_error: BaseException | None = None
    try:
        os.unlink(source, dir_fd=source_parent)
    except BaseException as error:
        operation_error = error
    durability_error = fsync_descriptors(source_parent, target_parent)
    state, _, _ = classify_bound_hard_link_move(
        source_parent,
        source,
        target_parent,
        target,
        expected,
        label,
        expected_sha256=digest,
    )
    require_mutation_authority_bound(source_parent)
    require_mutation_authority_bound(target_parent)
    if durability_error is not None:
        raise durability_error
    if state == "target-only":
        return
    if state == "linked" and operation_error is not None:
        raise operation_error
    if state != "linked":
        raise ValueError(f"{label} hard-link source unlink changed the selected names")
    raise ValueError(f"{label} hard-link source unlink did not complete")


def deterministic_retirement_name(
    parent: int, name: str, expected: os.stat_result
) -> str:
    """Return the discoverable forensic name for one exact retirement."""

    name = require_entry_name(name, "retirement source")
    parent_identity = profile_identity(os.fstat(parent))
    token = sha256_bytes(
        (
            f"{parent_identity[0]}:{parent_identity[1]}\0{name}\0"
            f"{expected.st_dev}:{expected.st_ino}"
        ).encode()
    )
    return f".cgl-checkpoint-retired-{token}.forensic"


def link_descriptor_noreplace(source: int, parent: int, target: str,
                              expected_identity: tuple[int, int], label: str) -> None:
    """Create one hard link and durably classify the raw linkat outcome."""

    require_mutation_authority_bound(parent)
    target = require_entry_name(target, label)
    source_profile = os.fstat(source)
    if profile_identity(source_profile) != expected_identity:
        raise ValueError(f"{label} source descriptor inode changed")
    require_regular_profile(source_profile, Path(target), f"{label} source")
    source_digest = sha256_descriptor(source)
    stable_source_profile = (
        profile_identity(source_profile),
        stat.S_IMODE(source_profile.st_mode),
        source_profile.st_uid,
        source_profile.st_gid,
        source_profile.st_size,
    )
    try:
        operation = ctypes.CDLL(None, use_errno=True).linkat
    except AttributeError as error:
        raise ValueError("descriptor-bound publication requires linkat") from error
    operation.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
    )
    operation.restype = ctypes.c_int
    operation_error: BaseException | None = None
    try:
        if operation(
            AT_FDCWD,
            os.fsencode(f"/proc/self/fd/{source}"),
            parent,
            os.fsencode(target),
            AT_SYMLINK_FOLLOW,
        ) != 0:
            retained_errno = ctypes.get_errno()
            raise OSError(retained_errno, os.strerror(retained_errno), target)
    except BaseException as error:
        operation_error = error
    durability_error = fsync_descriptors(parent)
    source_after = os.fstat(source)
    if (
        (
            profile_identity(source_after),
            stat.S_IMODE(source_after.st_mode),
            source_after.st_uid,
            source_after.st_gid,
            source_after.st_size,
        )
        != stable_source_profile
        or sha256_descriptor(source) != source_digest
    ):
        raise ValueError(f"{label} source descriptor changed during link attempt")
    try:
        target_profile = os.stat(target, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError:
        linked = False
    else:
        linked = profile_identity(target_profile) == expected_identity
        if linked:
            require_regular_profile(target_profile, Path(target), label)
            target_descriptor = os.open(
                target, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent
            )
            try:
                if (
                    profile_identity(os.fstat(target_descriptor)) != expected_identity
                    or sha256_descriptor(target_descriptor) != source_digest
                ):
                    linked = False
            finally:
                os.close(target_descriptor)
    if not linked and not entry_absent(parent, target):
        raise ValueError(f"{label} target differs after descriptor link attempt")
    require_mutation_authority_bound(parent)
    if durability_error is not None:
        raise durability_error
    if operation_error is not None:
        raise operation_error
    if not linked:
        raise ValueError(f"{label} link attempt did not publish the target")


def unlink_bound_entry(parent: int, name: str, expected: os.stat_result,
                       label: str, *, expected_sha256: str | None = None) -> None:
    """Retire one exact direct child without leaving an undiscoverable state."""

    require_mutation_authority_bound(parent)
    name = require_entry_name(name, label)
    retired = deterministic_retirement_name(parent, name, expected)
    forensic_parent = os.open("..", os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=parent)
    require_trusted_directory_profile(
        os.fstat(forensic_parent), Path(f"/proc/self/fd/{forensic_parent}"),
        f"{label} forensic parent",
    )

    try:
        readable = stat.S_IMODE(expected.st_mode) != 0
        if readable:
            try:
                state, expected_digest, _ = classify_bound_hard_link_move(
                    parent,
                    name,
                    forensic_parent,
                    retired,
                    expected,
                    f"{label} retirement",
                    expected_sha256=expected_sha256,
                )
            except ValueError as error:
                raise ValueError(
                    f"{label} changed during retirement; deterministic forensic "
                    f"name is {retired}"
                ) from error
            if state != "source-only":
                complete_bound_hard_link_move(
                    parent,
                    name,
                    forensic_parent,
                    retired,
                    expected,
                    f"{label} retirement",
                    expected_sha256=expected_digest,
                )
                return
        else:
            expected_digest = require_bound_entry_security(
                parent,
                name,
                expected,
                label,
                expected_sha256=expected_sha256,
                allow_mode_zero=True,
            )
            if not entry_absent(forensic_parent, retired):
                raise ValueError(
                    f"{label} deterministic forensic target already exists: {retired}"
                )
        operation_error: BaseException | None = None
        try:
            require_mutation_authority_bound(parent)
            require_trusted_directory_profile(
                os.fstat(forensic_parent),
                Path(f"/proc/self/fd/{forensic_parent}"),
                f"{label} forensic parent",
            )
            renameat2_between(
                parent,
                name,
                forensic_parent,
                retired,
                RENAME_NOREPLACE,
                f"{label} retirement",
            )
        except BaseException as error:
            operation_error = error
        durability_error = fsync_descriptors(parent, forensic_parent)
        if readable:
            source_matches = entry_security_matches(
                parent, name, expected, label, expected_digest
            )
            retired_matches = entry_security_matches(
                forensic_parent,
                retired,
                expected,
                f"{label} retirement",
                expected_digest,
            )
            if retired_matches and entry_absent(parent, name):
                state = "retired"
            elif retired_matches:
                raise ValueError(f"{label} public name reappeared during retirement")
            elif source_matches and entry_absent(forensic_parent, retired):
                state = "unchanged"
            else:
                raise ValueError(
                    f"{label} changed during atomic retirement; retained as {retired}"
                )
            if state == "unchanged" and renameat2_is_unsupported(operation_error):
                if durability_error is not None:
                    raise durability_error
                complete_bound_hard_link_move(
                    parent,
                    name,
                    forensic_parent,
                    retired,
                    expected,
                    f"{label} retirement",
                    expected_sha256=expected_digest,
                )
                return
            require_mutation_authority_bound(parent)
            if durability_error is not None:
                raise durability_error
            if operation_error is not None:
                raise operation_error
            if state != "retired":
                raise ValueError(f"{label} retirement did not remove the public name")
            return
        source_matches = entry_security_matches(
            parent,
            name,
            expected,
            label,
            expected_digest,
            allow_mode_zero=expected_digest is None,
        )
        retired_matches = entry_security_matches(
            forensic_parent,
            retired,
            expected,
            f"{label} retirement",
            expected_digest,
            allow_mode_zero=expected_digest is None,
        )
        if retired_matches and entry_absent(parent, name):
            state = "retired"
        elif retired_matches:
            raise ValueError(f"{label} public name reappeared during retirement")
        elif source_matches and entry_absent(forensic_parent, retired):
            state = "unchanged"
        else:
            raise ValueError(
                f"{label} changed during atomic retirement; retained as {retired}"
            )
        if state == "unchanged" and renameat2_is_unsupported(operation_error):
            if durability_error is not None:
                raise durability_error
            raise ValueError(
                f"{label} Lustre retirement requires exact readable-byte "
                "authentication; mode-0000 source left unchanged"
            ) from operation_error
        require_mutation_authority_bound(parent)
        if durability_error is not None:
            raise durability_error
        if operation_error is not None:
            raise operation_error
        if state != "retired":
            raise ValueError(f"{label} retirement did not remove the public name")
    finally:
        os.close(forensic_parent)


def rename_bound_noreplace(parent: int, source: str, target: str,
                           expected: os.stat_result, label: str, *,
                           expected_sha256: str | None = None) -> None:
    """Move one authenticated entry through a retryable no-clobber protocol."""

    require_mutation_authority_bound(parent)
    try:
        state, digest, _ = classify_bound_hard_link_move(
            parent,
            source,
            parent,
            target,
            expected,
            label,
            expected_sha256=expected_sha256,
        )
    except ValueError as error:
        if not entry_absent(parent, target):
            raise ValueError(f"{label} target already exists") from error
        raise
    if state != "source-only":
        complete_bound_hard_link_move(
            parent,
            source,
            parent,
            target,
            expected,
            label,
            expected_sha256=digest,
        )
        return
    operation_error: BaseException | None = None
    try:
        renameat2(parent, source, target, RENAME_NOREPLACE, label)
    except BaseException as error:
        operation_error = error
    try:
        state, _, _ = classify_bound_hard_link_move(
            parent,
            source,
            parent,
            target,
            expected,
            label,
            expected_sha256=digest,
        )
    except ValueError as error:
        if isinstance(operation_error, FileExistsError):
            raise ValueError(f"{label} target already exists") from operation_error
        if "content digest changed" in str(error):
            raise ValueError(f"{label} content digest changed during mutation") from error
        raise ValueError(f"{label} names changed during atomic no-replace") from error
    if state in {"source-only", "linked"} and renameat2_is_unsupported(operation_error):
        try:
            complete_bound_hard_link_move(
                parent,
                source,
                parent,
                target,
                expected,
                label,
                expected_sha256=digest,
            )
        except FileExistsError as error:
            raise ValueError(f"{label} target already exists") from error
        return
    require_mutation_authority_bound(parent)
    if operation_error is not None:
        raise operation_error
    if state != "target-only":
        raise ValueError(f"{label} no-replace did not publish the target")


def exchange_bound_entries(parent: int, source: str, target: str,
                           source_expected: os.stat_result,
                           target_expected: os.stat_result, label: str, *,
                           source_sha256: str | None = None,
                           target_sha256: str | None = None) -> None:
    """Atomically exchange two names and durably classify the resulting state."""

    require_mutation_authority_bound(parent)
    source_digest = require_bound_entry_security(
        parent,
        source,
        source_expected,
        label,
        expected_sha256=source_sha256,
    )
    target_digest = require_bound_entry_security(
        parent,
        target,
        target_expected,
        f"{label} target",
        expected_sha256=target_sha256,
    )
    assert source_digest is not None and target_digest is not None
    operation_error: BaseException | None = None
    try:
        renameat2(parent, source, target, RENAME_EXCHANGE, label)
    except BaseException as error:
        operation_error = error
    exchanged = (
        entry_security_matches(
            parent, target, source_expected, label, source_digest
        )
        and entry_security_matches(
            parent, source, target_expected, f"{label} predecessor", target_digest
        )
    )
    unchanged = (
        entry_security_matches(parent, source, source_expected, label, source_digest)
        and entry_security_matches(
            parent, target, target_expected, f"{label} target", target_digest
        )
    )
    if unchanged and renameat2_is_unsupported(operation_error):
        require_mutation_authority_bound(parent)
        raise ValueError(
            f"{label} occupied-target exchange requires RENAME_EXCHANGE support; "
            "both authenticated names were left unchanged"
        ) from operation_error
    if not exchanged and not unchanged:
        raise ValueError(f"{label} names changed after atomic exchange")
    require_mutation_authority_bound(parent)
    if operation_error is not None:
        raise operation_error
    if not exchanged:
        raise ValueError(f"{label} exchange did not change the selected names")


def replace_bound_entry_forward(
    parent: int,
    source: str,
    target: str,
    source_expected: os.stat_result,
    target_expected: os.stat_result,
    label: str,
    *,
    source_sha256: str,
    target_sha256: str,
) -> None:
    """Forward-replace one target through a durable authenticated commit marker."""

    source_sha256 = require_sha256(source_sha256, f"{label} source SHA-256")
    target_sha256 = require_sha256(target_sha256, f"{label} predecessor SHA-256")
    require_bound_entry_security(
        parent,
        source,
        source_expected,
        label,
        expected_sha256=source_sha256,
    )
    require_bound_entry_security(
        parent,
        target,
        target_expected,
        f"{label} predecessor",
        expected_sha256=target_sha256,
    )
    unlink_bound_entry(
        parent,
        target,
        target_expected,
        f"{label} predecessor",
        expected_sha256=target_sha256,
    )
    require_mutation_authority_bound(parent)
    rename_bound_noreplace(
        parent,
        source,
        target,
        source_expected,
        label,
        expected_sha256=source_sha256,
    )
    require_bound_entry_security(
        parent,
        target,
        source_expected,
        label,
        expected_sha256=source_sha256,
    )
    os.fsync(parent)
    require_mutation_authority_bound(parent)


@contextmanager
def bound_parent_descriptor(path: Path, label: str):
    """Yield a stable trusted no-symlink descriptor for one entry's parent."""

    with absolute_descriptor(
        path.parent, f"{label} parent", flags=os.O_RDONLY | os.O_DIRECTORY
    ) as parent:
        profile = os.fstat(parent)
        require_trusted_directory_profile(profile, path.parent, f"{label} parent")
        guard = BoundParentGuard(parent, path.parent, profile, label)
        guard.assert_bound()
        if parent in _BOUND_MUTATION_PARENTS:
            raise ValueError(f"{label} parent descriptor is already a mutation boundary")
        _BOUND_MUTATION_PARENTS[parent] = guard
        try:
            yield parent
        finally:
            try:
                guard.assert_bound()
            finally:
                _BOUND_MUTATION_PARENTS.pop(parent, None)


@contextmanager
def bound_directory_descriptor(path: Path, label: str):
    """Yield one trusted directory continuously bound to its public pathname."""

    path = absolute_path(path)
    with absolute_descriptor(
        path, label, flags=os.O_RDONLY | os.O_DIRECTORY
    ) as descriptor:
        profile = os.fstat(descriptor)
        require_trusted_directory_profile(profile, path, label)
        guard = BoundParentGuard(descriptor, path, profile, label)
        guard.assert_bound()
        if descriptor in _BOUND_MUTATION_PARENTS:
            raise ValueError(f"{label} descriptor is already a mutation boundary")
        _BOUND_MUTATION_PARENTS[descriptor] = guard
        try:
            yield descriptor
        finally:
            try:
                guard.assert_bound()
            finally:
                _BOUND_MUTATION_PARENTS.pop(descriptor, None)


def directory_entry_names(descriptor: int) -> list[str]:
    """Return sorted direct-child names through one bound directory descriptor."""

    return sorted(os.listdir(descriptor))


def durably_classify_created_regular(
    parent: int,
    name: str,
    descriptor: int | None,
    label: str,
    *,
    expected_mode: int,
    operation_error: BaseException | None = None,
) -> os.stat_result:
    """Persist and authenticate one exact O_CREAT outcome before continuing."""

    name = require_entry_name(name, label)
    durability_error = fsync_descriptors(
        *(() if descriptor is None else (descriptor,)),
        parent,
    )
    retained = descriptor
    close_retained = False
    try:
        if retained is None:
            try:
                retained = os.open(
                    name,
                    getattr(os, "O_PATH", os.O_RDONLY) | os.O_NOFOLLOW,
                    dir_fd=parent,
                )
                close_retained = True
            except FileNotFoundError:
                require_mutation_authority_bound(parent)
                if durability_error is not None:
                    raise durability_error
                if operation_error is not None:
                    raise operation_error
                raise ValueError(f"{label} create outcome disappeared")
        profile = os.fstat(retained)
        require_regular_profile(
            profile,
            Path(name),
            label,
            expected_mode=expected_mode,
            expected_links=1,
            expected_uid=os.geteuid(),
        )
        if profile.st_size != 0:
            raise ValueError(f"{label} create outcome is not empty: {name}")
        require_bound_entry_profile(parent, name, profile, label)
    finally:
        if close_retained and retained is not None:
            os.close(retained)
    require_mutation_authority_bound(parent)
    if durability_error is not None:
        raise durability_error
    if operation_error is not None:
        raise operation_error
    return profile


def write_bound_exclusive(parent: int, name: str, payload: bytes, mode: int,
                          label: str) -> os.stat_result:
    """Create exact bytes as one authenticated direct child."""

    require_mutation_authority_bound(parent)
    name = require_entry_name(name, label)
    expected_digest = sha256_bytes(payload)
    descriptor = None
    operation_error = None
    previous = os.umask(0)
    try:
        try:
            descriptor = os.open(
                name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o600,
                dir_fd=parent,
            )
        except BaseException as error:
            operation_error = error
    finally:
        try:
            os.umask(previous)
        except BaseException as error:
            if operation_error is None:
                operation_error = error
    try:
        durably_classify_created_regular(
            parent,
            name,
            descriptor,
            label,
            expected_mode=0o600,
            operation_error=operation_error,
        )
    except BaseException:
        if descriptor is not None:
            os.close(descriptor)
        raise
    assert descriptor is not None
    try:
        write_descriptor_bytes(descriptor, payload)
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
        profile = os.fstat(descriptor)
        require_regular_profile(
            profile,
            Path(name),
            label,
            expected_mode=mode,
            expected_links=1,
            expected_uid=os.geteuid(),
        )
    except BaseException:
        profile = os.fstat(descriptor)
        os.close(descriptor)
        try:
            unlink_bound_entry(parent, name, profile, label)
        except (FileNotFoundError, ValueError):
            pass
        raise
    os.close(descriptor)
    require_bound_entry_security(
        parent, name, profile, label, expected_sha256=expected_digest
    )
    os.fsync(parent)
    require_bound_entry_security(
        parent, name, profile, label, expected_sha256=expected_digest
    )
    require_mutation_authority_bound(parent)
    return profile


def create_exclusive(path: Path, flags: int, mode: int) -> int:
    """Create one exact-mode entry without allowing the caller umask to narrow it."""

    previous = os.umask(0)
    try:
        return os.open(path, flags | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, mode)
    finally:
        os.umask(previous)


def json_temporary_target_name(name: str) -> str | None:
    """Return the exact public target encoded by one JSON atomic temporary."""

    name = require_entry_name(name, "JSON atomic temporary")
    match = re.fullmatch(
        r"\.(?P<target>.+)\.[0-9]+\.[0-9a-f]{32}\.tmp",
        name,
    )
    if match is None:
        return None
    return require_entry_name(match.group("target"), "JSON atomic temporary target")


def json_temporary_names_in_directory(parent: int, target: str) -> list[str]:
    """Return every narrowly named atomic temporary for one public target."""

    target = require_entry_name(target, "JSON atomic target")
    return [
        name
        for name in directory_entry_names(parent)
        if json_temporary_target_name(name) == target
    ]


def recover_linked_json_publication_in_directory(
    parent: int,
    directory: Path,
    target: str,
    label: str,
    *,
    expected_mode: int | None = None,
    expected_sha256: str | None = None,
) -> bool:
    """Complete only one exact post-link/pre-unlink JSON publication state."""

    require_mutation_authority_bound(parent)
    target = require_entry_name(target, label)
    temporaries = json_temporary_names_in_directory(parent, target)
    temporary_profiles: dict[str, os.stat_result] = {}
    for temporary in temporaries:
        try:
            temporary_profiles[temporary] = os.stat(
                temporary, dir_fd=parent, follow_symlinks=False
            )
        except FileNotFoundError as error:
            raise ValueError(f"{label} temporary namespace changed during recovery") from error
    try:
        public_profile = os.stat(target, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError:
        public_profile = None

    linked_temporaries = [
        name for name, profile in temporary_profiles.items() if profile.st_nlink != 1
    ]
    if public_profile is None:
        if linked_temporaries:
            raise ValueError(
                f"{label} has a linked temporary without its exact public name"
            )
        return False
    if public_profile.st_nlink == 1:
        if linked_temporaries:
            raise ValueError(
                f"{label} temporary link state differs from its public name"
            )
        return False
    if public_profile.st_nlink != 2:
        raise ValueError(
            f"{label} public name has {public_profile.st_nlink} links, expected 1 or "
            "one exact recoverable pair"
        )
    if len(temporaries) != 1:
        raise ValueError(
            f"{label} two-link public name requires exactly one correctly named temporary"
        )

    temporary = temporaries[0]
    temporary_profile = temporary_profiles[temporary]
    if (
        temporary_profile.st_nlink != 2
        or profile_identity(temporary_profile) != profile_identity(public_profile)
    ):
        raise ValueError(
            f"{label} correctly named temporary and public name do not select the "
            "same exact two-link inode"
        )
    require_regular_profile(
        temporary_profile,
        directory / temporary,
        f"{label} temporary",
        expected_mode=expected_mode,
        expected_links=2,
        expected_uid=os.geteuid(),
    )
    require_regular_profile(
        public_profile,
        directory / target,
        f"{label} public name",
        expected_mode=expected_mode,
        expected_links=2,
        expected_uid=os.geteuid(),
    )
    state, digest, links = classify_bound_hard_link_move(
        parent,
        temporary,
        parent,
        target,
        temporary_profile,
        label,
        expected_sha256=expected_sha256,
    )
    if state != "linked" or links != 2:
        raise ValueError(f"{label} is not the exact recoverable two-link state")
    complete_bound_hard_link_move(
        parent,
        temporary,
        parent,
        target,
        temporary_profile,
        label,
        expected_sha256=digest,
    )
    if not entry_absent(parent, temporary):
        raise ValueError(f"{label} temporary remains after recovery")
    current = os.stat(target, dir_fd=parent, follow_symlinks=False)
    if profile_identity(current) != profile_identity(public_profile):
        raise ValueError(f"{label} public inode changed during recovery")
    require_regular_profile(
        current,
        directory / target,
        f"{label} public name",
        expected_mode=expected_mode,
        expected_links=1,
        expected_uid=os.geteuid(),
    )
    require_bound_entry_security(
        parent,
        target,
        current,
        f"{label} public name",
        expected_sha256=digest,
    )
    require_mutation_authority_bound(parent)
    return True


def recover_linked_json_publication(
    path: Path,
    label: str,
    *,
    expected_mode: int | None = None,
    expected_sha256: str | None = None,
) -> bool:
    """Recover one exact JSON post-link/pre-unlink state through its bound parent."""

    with bound_parent_descriptor(path, label) as parent:
        return recover_linked_json_publication_in_directory(
            parent,
            path.parent,
            path.name,
            label,
            expected_mode=expected_mode,
            expected_sha256=expected_sha256,
        )


def write_json(path: Path, value: object, mode: int = 0o644, *,
               simulate_interruption_before_directory_fsync: bool = False) -> None:
    """Atomically replace one exact JSON entry through its stable parent."""

    payload = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    payload_digest = sha256_bytes(payload)
    temporary = f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp"
    with bound_parent_descriptor(path, str(path)) as parent:
        if recover_linked_json_publication_in_directory(
            parent,
            path.parent,
            path.name,
            "JSON atomic write retry",
            expected_mode=mode,
            expected_sha256=payload_digest,
        ):
            return
        try:
            target_profile = os.stat(path.name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            target_profile = None
            target_digest = None
        if target_profile is not None:
            require_regular_profile(
                target_profile,
                path,
                "JSON predecessor",
                expected_links=1,
                expected_uid=os.geteuid(),
            )
            target_digest = require_bound_entry_security(
                parent,
                path.name,
                target_profile,
                "JSON predecessor",
            )
        temporary_profile = write_bound_exclusive(
            parent, temporary, payload, mode, "JSON atomic temporary"
        )
        try:
            if target_profile is None:
                rename_bound_noreplace(
                    parent,
                    temporary,
                    path.name,
                    temporary_profile,
                    "JSON atomic write",
                    expected_sha256=payload_digest,
                )
            else:
                assert target_digest is not None
                replace_bound_entry_forward(
                    parent,
                    temporary,
                    path.name,
                    temporary_profile,
                    target_profile,
                    "JSON atomic write",
                    source_sha256=payload_digest,
                    target_sha256=target_digest,
                )
                if simulate_interruption_before_directory_fsync:
                    raise ValueError("simulated interruption before JSON directory fsync")
            if target_profile is None and simulate_interruption_before_directory_fsync:
                raise ValueError("simulated interruption before JSON directory fsync")
            require_bound_entry_security(
                parent,
                path.name,
                temporary_profile,
                "JSON atomic write",
                expected_sha256=payload_digest,
            )
            os.fsync(parent)
            require_bound_entry_security(
                parent,
                path.name,
                temporary_profile,
                "JSON atomic write",
                expected_sha256=payload_digest,
            )
            require_mutation_authority_bound(parent)
        except BaseException:
            raise


def json_temporary_entries(path: Path) -> list[Path]:
    """Return narrowly matched atomic-write temporaries for one JSON path."""

    return [
        entry for entry in directory_entries(path.parent)
        if json_temporary_target_name(entry.name) == path.name
    ]


def read_bound_json_object(parent: int, name: str, path: Path, label: str, *,
                           expected_mode: int | None = None
                           ) -> tuple[os.stat_result, str, bytes, dict[str, object]]:
    """Read one exact descriptor-relative JSON object with stable profile/content."""

    name = require_entry_name(name, label)
    profile = os.stat(name, dir_fd=parent, follow_symlinks=False)
    require_regular_mode_subset(
        profile,
        path,
        label,
        maximum_mode=0o644,
        expected_links=1,
        expected_uid=os.geteuid(),
    )
    if expected_mode is not None and stat.S_IMODE(profile.st_mode) != expected_mode:
        raise ValueError(
            f"{label} mode is {stat.S_IMODE(profile.st_mode):04o}, "
            f"expected {expected_mode:04o}: {path}"
        )
    digest = require_bound_entry_security(parent, name, profile, label)
    assert digest is not None
    descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent)
    try:
        retained = read_descriptor_bytes(descriptor)
        if sha256_bytes(retained) != digest:
            raise ValueError(f"{label} content changed while reading: {path}")
    finally:
        os.close(descriptor)
    require_bound_entry_security(
        parent, name, profile, label, expected_sha256=digest
    )
    try:
        record = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON: {path}") from error
    if not isinstance(record, dict):
        raise ValueError(f"{label} must be a JSON object: {path}")
    return profile, digest, retained, record


def remove_json_temporaries_in_directory(parent: int, directory: Path,
                                         names: list[str], label: str) -> None:
    """Retire invalid partial JSON temporaries through one bound directory."""

    for name in names:
        path = directory / require_entry_name(name, label)
        try:
            profile = os.stat(name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            continue
        require_regular_mode_subset(
            profile,
            path,
            label,
            maximum_mode=0o644,
            expected_links=1,
            expected_uid=os.geteuid(),
        )
        valid = False
        if stat.S_IMODE(profile.st_mode) != 0:
            try:
                read_bound_json_object(parent, name, path, label)
            except ValueError as error:
                if "not valid JSON" not in str(error) and "must be a JSON object" not in str(
                    error
                ):
                    raise
            else:
                valid = True
        if valid:
            raise ValueError(
                f"{label} is a valid JSON recovery candidate and was preserved: {path}"
            )
        unlink_bound_entry(parent, name, profile, label)
    require_mutation_authority_bound(parent)


def remove_json_temporaries(entries: list[Path], label: str) -> None:
    """Retire only partial JSON temporaries; preserve valid recovery bytes."""

    grouped: dict[Path, list[str]] = {}
    for entry in entries:
        grouped.setdefault(entry.parent, []).append(entry.name)
    for directory, names in grouped.items():
        with bound_directory_descriptor(directory, label) as parent:
            remove_json_temporaries_in_directory(parent, directory, names, label)


def mkdir_durable(path: Path) -> None:
    """Create one directory tree through stable parent descriptors."""

    missing = []
    current = path
    while not current.exists():
        missing.append(current)
        current = current.parent
    for directory in reversed(missing):
        with bound_parent_descriptor(directory, str(directory)) as parent:
            require_mutation_authority_bound(parent)
            parent_profile = os.fstat(parent)
            expected_mode = 0o755 | (
                stat.S_ISGID if stat.S_IMODE(parent_profile.st_mode) & stat.S_ISGID else 0
            )
            operation_error = None
            previous = os.umask(0)
            try:
                try:
                    os.mkdir(directory.name, mode=0o755, dir_fd=parent)
                except BaseException as error:
                    operation_error = error
            finally:
                try:
                    os.umask(previous)
                except BaseException as error:
                    if operation_error is None:
                        operation_error = error
            durability_error = fsync_descriptors(parent)
            descriptor = None
            try:
                try:
                    descriptor = os.open(
                        directory.name,
                        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                        dir_fd=parent,
                    )
                except BaseException:
                    descriptor = os.open(
                        directory.name,
                        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                        dir_fd=parent,
                    )
                profile = os.fstat(descriptor)
                require_trusted_directory_profile(profile, directory, str(directory))
                if stat.S_IMODE(profile.st_mode) != expected_mode:
                    raise ValueError(
                        f"{directory} mode is {stat.S_IMODE(profile.st_mode):04o}, "
                        f"expected {expected_mode:04o}: {directory}"
                    )
                require_bound_entry_profile(
                    parent, directory.name, profile, str(directory)
                )
                followup_error = fsync_descriptors(descriptor, parent)
                require_bound_entry_profile(
                    parent, directory.name, profile, str(directory)
                )
            finally:
                if descriptor is not None:
                    os.close(descriptor)
            require_mutation_authority_bound(parent)
            if durability_error is not None:
                raise durability_error
            if followup_error is not None:
                raise followup_error
            if operation_error is not None:
                raise operation_error


def unlink_durable(path: Path) -> None:
    """Remove only the exact inode authenticated through a stable parent."""

    with bound_parent_descriptor(path, str(path)) as parent:
        descriptor = os.open(path.name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent)
        try:
            profile = os.fstat(descriptor)
            require_regular_profile(
                profile, path, str(path), expected_links=1, expected_uid=os.geteuid()
            )
            require_bound_entry_identity(parent, path.name, profile, str(path))
        finally:
            os.close(descriptor)
        unlink_bound_entry(parent, path.name, profile, str(path))


def entry_exists(path: Path) -> bool:
    """Return whether one directory entry exists without following symlinks."""

    try:
        path.lstat()
    except FileNotFoundError:
        return False
    return True


def require_regular_profile(value: os.stat_result, path: Path, label: str, *,
                            expected_mode: int | None = None,
                            expected_links: int | None = None,
                            expected_uid: int | None = None) -> None:
    """Require one descriptor-backed regular-file profile."""

    if not stat.S_ISREG(value.st_mode):
        raise ValueError(f"{label} is not a regular file: {path}")
    mode = stat.S_IMODE(value.st_mode)
    if expected_mode is not None and mode != expected_mode:
        raise ValueError(
            f"{label} mode is {mode:04o}, expected {expected_mode:04o}: {path}"
        )
    if expected_links is not None and value.st_nlink != expected_links:
        raise ValueError(
            f"{label} has {value.st_nlink} links, expected {expected_links}: {path}"
        )
    if expected_uid is not None and value.st_uid != expected_uid:
        raise ValueError(
            f"{label} owner is UID {value.st_uid}, expected {expected_uid}: {path}"
        )


def require_regular_mode_subset(value: os.stat_result, path: Path, label: str, *,
                                maximum_mode: int,
                                expected_links: int,
                                expected_uid: int | None = None) -> None:
    """Require a retained file with no permissions beyond a selected cap."""

    require_regular_profile(
        value,
        path,
        label,
        expected_links=expected_links,
        expected_uid=expected_uid,
    )
    mode = stat.S_IMODE(value.st_mode)
    if mode & ~maximum_mode:
        raise ValueError(
            f"{label} mode is {mode:04o}, exceeds {maximum_mode:04o}: {path}"
        )


@contextmanager
def regular_descriptor(path: Path, label: str, *,
                       flags: int = os.O_RDONLY,
                       expected_mode: int | None = None,
                       expected_links: int | None = None):
    """Yield one O_NOFOLLOW descriptor after validating its regular profile."""

    try:
        descriptor = os.open(path, flags | os.O_NOFOLLOW)
    except FileNotFoundError as error:
        raise ValueError(f"{label} is missing: {path}") from error
    try:
        require_regular_profile(
            os.fstat(descriptor),
            path,
            label,
            expected_mode=expected_mode,
            expected_links=expected_links,
        )
        yield descriptor
    finally:
        os.close(descriptor)


def require_file_sha256(path: Path, expected: str, label: str, *,
                        expected_mode: int | None = None,
                        expected_links: int | None = None) -> bytes:
    """Return descriptor-read bytes after authenticating one retained file."""

    with regular_descriptor(
        path,
        label,
        expected_mode=expected_mode,
        expected_links=expected_links,
    ) as descriptor:
        value = read_descriptor_bytes(descriptor)
    if sha256_bytes(value) != require_sha256(expected, f"{label} SHA-256"):
        raise ValueError(f"{label} checksum has changed: {path}")
    return value


def safe_relative_path(value: str, label: str) -> Path:
    """Require a relative path without traversal."""

    path = Path(value)
    if path.is_absolute() or not path.parts or ".." in path.parts:
        raise ValueError(f"{label} must be a relative path without '..': {value}")
    return path


def absolute_path(value: Path) -> Path:
    """Return one lexical absolute path without resolving symlinks."""

    return Path(os.path.abspath(os.path.expanduser(str(value))))


def require_offline_path(path: Path, label: str) -> Path:
    """Require offline fixture paths beneath an explicit local prefix."""

    path = absolute_path(path)
    for allowed in OFFLINE_ALLOWED_PREFIXES:
        allowed = absolute_path(allowed)
        if path == allowed or allowed in path.parents:
            return path
    allowed = ", ".join(str(absolute_path(path)) for path in OFFLINE_ALLOWED_PREFIXES)
    raise ValueError(f"{label} must be under local fixture storage {allowed}: {path}")


@contextmanager
def absolute_descriptor(path: Path, label: str, *, flags: int):
    """Open an absolute path without accepting symlinks in any component."""

    path = absolute_path(path)
    if not path.is_absolute() or len(path.parts) < 2:
        raise ValueError(f"{label} path is invalid: {path}")
    parent = os.open("/", os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    opened = [parent]
    try:
        for component in path.parts[1:-1]:
            parent = os.open(
                component,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=parent,
            )
            opened.append(parent)
        descriptor = os.open(path.parts[-1], flags | os.O_NOFOLLOW, dir_fd=parent)
        opened.append(descriptor)
        yield descriptor
    finally:
        for descriptor in reversed(opened):
            os.close(descriptor)


def read_confined_file_sha256(root: Path, relative_value: str, expected: str,
                              label: str, *, expected_mode: int
                              ) -> tuple[Path, bytes]:
    """Read one authenticated in-root file without accepting symlink traversal."""

    relative = safe_relative_path(relative_value, f"{label} relative path")
    parent = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    opened = [parent]
    try:
        for component in relative.parts[:-1]:
            parent = os.open(
                component,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=parent,
            )
            opened.append(parent)
        descriptor = os.open(
            relative.parts[-1],
            os.O_RDONLY | os.O_NOFOLLOW,
            dir_fd=parent,
        )
        opened.append(descriptor)
        require_regular_profile(
            os.fstat(descriptor),
            root / relative,
            label,
            expected_mode=expected_mode,
            expected_links=1,
        )
        retained = read_descriptor_bytes(descriptor)
        if sha256_bytes(retained) != require_sha256(expected, f"{label} SHA-256"):
            raise ValueError(f"{label} checksum has changed: {root / relative}")
    finally:
        for descriptor in reversed(opened):
            os.close(descriptor)
    return root / relative, retained


def require_confined_file_sha256(root: Path, relative_value: str, expected: str,
                                 label: str, *, expected_mode: int) -> Path:
    """Authenticate one in-root file without accepting symlink traversal."""

    path, _ = read_confined_file_sha256(
        root, relative_value, expected, label, expected_mode=expected_mode
    )
    return path


def require_artifact_name(value: str) -> str:
    """Require a basename suitable for exact namespace checks."""

    if ARTIFACT_NAME_PATTERN.fullmatch(value) is None:
        raise ValueError(
            "artifact name must be a JSON basename using ASCII letters, digits, "
            "dots, underscores, or hyphens"
        )
    return value


def require_root(root: Path, allow_local_root: bool) -> tuple[Path, bool]:
    """Select the canonical root or an explicitly enabled offline fixture."""

    requested = absolute_path(root)
    canonical = absolute_path(DEFAULT_ROOT)
    if requested == canonical:
        if allow_local_root:
            raise ValueError("--allow-local-root may not relax the canonical root")
        if root.expanduser().resolve() != canonical:
            raise ValueError(f"canonical Stage I root resolves unexpectedly: {root}")
        return canonical, False
    if not allow_local_root:
        raise ValueError(
            f"Stage I root must be {canonical}; use --allow-local-root only "
            "for offline fixtures"
        )
    requested = require_offline_path(requested, "offline fixture root")
    with absolute_descriptor(
        requested, "offline fixture root", flags=os.O_RDONLY | os.O_DIRECTORY
    ):
        pass
    return requested, True


def require_canonical_repository(root_dir: Path, offline: bool) -> None:
    """Pin canonical publication and verification to the retained source repository."""

    if not offline and root_dir != CANONICAL_REPOSITORY_ROOT:
        raise ValueError(
            "canonical use requires the retained repository root "
            f"{CANONICAL_REPOSITORY_ROOT}: {root_dir}"
        )


def layout(root: Path, artifact_name: str) -> dict[str, Path]:
    """Return companion and Stage I retained paths."""

    accounting = root / "accounting"
    canonical = accounting / require_artifact_name(artifact_name)
    return {
        "root": root,
        "accounting": accounting,
        "stage_i_transactions": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_transactions"
        ),
        "recost_transactions": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_transactions"
        ),
        "recost_forensics_root": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_forensics"
        ),
        "recost_forensics": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_forensics"
            / canonical.name
        ),
        "lock": root / f".mks24_stage_i_{EXECUTION_EPOCH_SLUG}.lock",
        "canonical": canonical,
        "staged": accounting / f"{artifact_name}.staged",
        "audit": accounting / f"{artifact_name}.publication_audit.json",
    }


def require_directory(path: Path, label: str) -> None:
    """Require one retained directory without accepting a symlink."""

    try:
        value = path.lstat()
    except FileNotFoundError as error:
        raise ValueError(f"{label} is missing: {path}") from error
    if not stat.S_ISDIR(value.st_mode):
        raise ValueError(f"{label} is not a directory: {path}")


def require_trusted_directory(path: Path, label: str) -> None:
    """Require one owned directory without unsafe replacement permissions."""

    require_directory(path, label)
    require_trusted_directory_profile(path.lstat(), path, label)


def require_trusted_directory_profile(value: os.stat_result, path: Path,
                                      label: str) -> None:
    """Require one descriptor- or pathname-backed trusted directory profile."""

    if not stat.S_ISDIR(value.st_mode):
        raise ValueError(f"{label} is not a directory: {path}")
    if value.st_uid != os.geteuid():
        raise ValueError(
            f"{label} owner is UID {value.st_uid}, expected {os.geteuid()}: {path}"
        )
    mode = stat.S_IMODE(value.st_mode)
    allowed = 0o755 | stat.S_ISGID
    if mode & ~allowed:
        raise ValueError(
            f"{label} mode is {mode:04o}, exceeds trusted profile 0755: {path}"
        )


def require_managed_directory(path: Path, label: str) -> None:
    """Require one companion-managed trusted directory."""

    require_trusted_directory(path, label)


def directory_entries(path: Path) -> list[Path]:
    """Return sorted direct children of one directory."""

    return sorted(path.iterdir(), key=lambda item: item.name)


def require_empty_directory(path: Path, label: str, *,
                            allow_absent: bool = False,
                            managed: bool = False,
                            trusted: bool = False) -> None:
    """Require an empty retained directory."""

    if not entry_exists(path) and allow_absent:
        return
    if managed or trusted:
        require_managed_directory(path, label)
    else:
        require_directory(path, label)
    entries = directory_entries(path)
    if entries:
        raise ValueError(
            f"{label} is not empty: " + ", ".join(str(item) for item in entries)
        )


def require_exact_directory_entries(path: Path, label: str,
                                    expected: set[Path], *,
                                    managed: bool = False) -> None:
    """Require one retained directory to contain exactly the selected entries."""

    if managed:
        require_managed_directory(path, label)
    else:
        require_directory(path, label)
    retained = set(directory_entries(path))
    if retained != expected:
        raise ValueError(
            f"{label} entries differ; found: "
            + ", ".join(str(item) for item in sorted(retained))
        )


def is_artifact_namespace_name(name: str, candidate: str) -> bool:
    """Return whether one name is companion-managed artifact state."""

    if candidate == f"{name}.independent_review.json":
        return False
    return (
        candidate == name
        or candidate.startswith(f"{name}.")
        or candidate == f".{name}"
        or candidate.startswith(f".{name}.")
    )


def artifact_namespace_entries(paths: dict[str, Path]) -> list[Path]:
    """Return direct accounting entries related to one artifact basename."""

    accounting = paths["accounting"]
    require_trusted_directory(accounting, "accounting directory")
    name = paths["canonical"].name
    return [
        path for path in directory_entries(accounting)
        if is_artifact_namespace_name(name, path.name)
    ]


def require_artifact_namespace(paths: dict[str, Path], state: str) -> None:
    """Require an exact staged-only or canonical-only artifact namespace."""

    if state == "staged":
        expected = {paths["staged"]}
    elif state == "promoted":
        expected = {paths["canonical"], paths["audit"]}
    elif state == "linked-pair":
        expected = {paths["staged"], paths["canonical"]}
    elif state == "canonical-pending-audit":
        expected = {paths["canonical"]}
    else:
        raise ValueError(f"unsupported namespace state: {state}")
    retained = set(artifact_namespace_entries(paths))
    if retained != expected:
        raise ValueError(
            f"recost namespace is not {state}-only; found: "
            + ", ".join(str(path) for path in sorted(retained))
        )


def json_pointer(value: object, pointer: str, label: str) -> object:
    """Resolve one RFC 6901 JSON pointer."""

    if pointer == "":
        return value
    if not pointer.startswith("/"):
        raise ValueError(f"{label} JSON pointer must start with '/': {pointer}")
    current = value
    for raw in pointer[1:].split("/"):
        token = raw.replace("~1", "/").replace("~0", "~")
        if isinstance(current, dict) and token in current:
            current = current[token]
            continue
        if isinstance(current, list) and re.fullmatch(r"0|[1-9][0-9]*", token):
            index = int(token)
            if index < len(current):
                current = current[index]
                continue
        raise ValueError(f"{label} JSON pointer is absent: {pointer}")
    return current


def require_exact_keys(value: object, keys: set[str] | frozenset[str],
                       label: str) -> dict[str, object]:
    """Require one object with an exact key set."""

    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    if frozenset(value) != frozenset(keys):
        raise ValueError(
            f"{label} schema differs; expected {sorted(keys)}, "
            f"found {sorted(value)}"
        )
    return value


def require_nonempty_string(value: object, label: str) -> str:
    """Require one nonempty trimmed string."""

    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError(f"{label} must be a nonempty trimmed string")
    return value


def require_safe_id(value: object, label: str) -> str:
    """Require one filesystem-safe identifier."""

    retained = require_nonempty_string(value, label)
    if SAFE_ID_PATTERN.fullmatch(retained) is None:
        raise ValueError(f"{label} is not a safe identifier")
    return retained


def require_integer(value: object, label: str, *, minimum: int = 0) -> int:
    """Require one non-boolean integer."""

    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


def require_finite_float(value: object, label: str) -> float:
    """Require one finite non-boolean number."""

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a finite number")
    retained = float(value)
    if not (float("-inf") < retained < float("inf")):
        raise ValueError(f"{label} must be a finite number")
    return retained


def walltime_seconds(value: object, label: str) -> int:
    """Parse one strict Stage I HH:MM:SS walltime."""

    retained = require_nonempty_string(value, label)
    match = re.fullmatch(r"([0-9]{2}):([0-5][0-9]):([0-5][0-9])", retained)
    if match is None:
        raise ValueError(f"{label} must use HH:MM:SS")
    hours, minutes, seconds = (int(item) for item in match.groups())
    total = hours * 3600 + minutes * 60 + seconds
    if total <= 0 or total > MAX_WALLTIME_SECONDS:
        raise ValueError(f"{label} must be within the two-hour Stage I limit")
    return total


def schema2_evidence_mode(args: argparse.Namespace) -> bool:
    """Return whether this invocation publishes schema-2 recommendation evidence."""

    return getattr(args, "recost_recommendations_json", None) is not None


def v2_packet(args: argparse.Namespace) -> dict[str, object] | None:
    """Return the selected generalized packet without relabeling evidence authority."""

    if schema2_evidence_mode(args):
        return args.recost_recommendations_json
    return args.authorized_bounded_wave_json


def publication_mode(args: argparse.Namespace) -> str:
    """Require one complete, non-mixed publication packet mode."""

    sole = args.authorized_next_segment_profile_json
    bounded = v2_packet(args)
    if sole is not None:
        if args.scheduler_relative_path is None or args.expected_scheduler_sha256 is None:
            raise ValueError(
                "sole-profile mode requires one scheduler relative path and SHA-256"
            )
        if (
            args.expected_source_bundle_verified_revisions_json is None
            or not args.expected_source_bundle_verified_revisions_json
        ):
            raise ValueError(
                "sole-profile mode requires source-bundle verified revisions"
            )
        if (
            not isinstance(sole, dict)
            or sole.get("source_bundle") is None
            or sole.get("source_bundle_sha256") is None
        ):
            raise ValueError(
                "sole-profile mode requires a complete legacy source-bundle binding"
            )
        if (
            args.bounded_wave_scheduler_evidence_json is not None
            or args.source_bundle_relative_path is not None
            or args.expected_source_bundle_sha256 is not None
            or args.expected_stage_i_revision is not None
            or args.expected_generator_revision is not None
            or args.recost_request_relative_path is not None
            or args.expected_request_sha256 is not None
        ):
            raise ValueError("sole-profile mode rejects V2 bindings")
        return "sole-profile"
    if bounded is not None:
        if args.scheduler_relative_path is not None or args.expected_scheduler_sha256 is not None:
            raise ValueError("V2 mode rejects sole-profile scheduler bindings")
        if args.bounded_wave_scheduler_evidence_json is None:
            raise ValueError("V2 mode requires scheduler-evidence JSON")
        if (
            args.source_bundle_relative_path is None
            or args.expected_source_bundle_sha256 is None
            or args.expected_stage_i_revision is None
            or args.expected_generator_revision is None
            or args.expected_source_bundle_verified_revisions_json is None
            or args.recost_request_relative_path is None
            or args.expected_request_sha256 is None
        ):
            raise ValueError("V2 mode requires exact request and source provenance")
        if (
            args.expected_stage_i_revision not in
            args.expected_source_bundle_verified_revisions_json
            or args.expected_generator_revision not in
            args.expected_source_bundle_verified_revisions_json
        ):
            raise ValueError(
                "bounded-wave source bundle revisions must cover helper and generator"
            )
        expected_pointers = {
            "artifact_epoch_pointer": "/execution_epoch",
            "artifact_authorization_pointer": (
                "/recommendations" if schema2_evidence_mode(args) else "/authorization"
            ),
            "artifact_counts_pointer": "/reconcile/counts",
            "artifact_stage_i_sha256_pointer": "/provenance/stage_i_helper_sha256",
            "artifact_generator_sha256_pointer": "/provenance/generator_sha256",
            "artifact_stage_i_revision_pointer": "/provenance/stage_i_helper_revision",
            "artifact_generator_revision_pointer": "/provenance/generator_revision",
            "artifact_scheduler_evidence_pointer": "/provenance/scheduler_evidence",
            "artifact_barrier_scheduler_evidence_pointer": (
                "/barrier/scheduler_evidence"
            ),
            "artifact_source_bundle_sha256_pointer": (
                "/provenance/source_bundle_sha256"
            ),
            "artifact_source_bundle_revisions_pointer": (
                "/provenance/source_bundle_verified_revisions"
            ),
        }
        for attribute, expected in expected_pointers.items():
            if getattr(args, attribute) != expected:
                raise ValueError(
                    f"V2 mode requires canonical JSON pointer {expected}"
                )
        generator = safe_relative_path(
            args.generator_relative_path, "recost generator relative path"
        )
        if generator.parent != Path("accounting/utilities"):
            raise ValueError(
                "V2 recost generator must be retained under accounting/utilities"
            )
        source = safe_relative_path(
            args.source_bundle_relative_path, "source bundle relative path"
        )
        if source.parent != Path("source-archives"):
            raise ValueError(
                "V2 source bundle must be a direct child of source-archives"
            )
        request = safe_relative_path(
            args.recost_request_relative_path, "recost request relative path"
        )
        if request.parent != Path("accounting"):
            raise ValueError("V2 recost request must be a direct child of accounting")
        if schema2_evidence_mode(args):
            if args.expected_artifact_mode != 0o444:
                raise ValueError("schema-2 recost evidence must be staged and published mode 0444")
            if (
                args.independent_review_relative_path is None
                or args.expected_independent_review_sha256 is None
            ):
                raise ValueError(
                    "schema-2 recost evidence requires one exact accounting review binding"
                )
            review = safe_relative_path(
                args.independent_review_relative_path,
                "recost evidence independent-review relative path",
            )
            if review.parent != Path("accounting"):
                raise ValueError(
                    "schema-2 recost evidence requires one exact accounting review binding"
                )
        return "v2"
    raise ValueError("one sole-profile or bounded-wave authorization is required")


def require_root_member(root: Path, value: object, label: str) -> Path:
    """Require one normalized absolute path lexically beneath the root."""

    retained = require_nonempty_string(value, label)
    path = Path(retained)
    if not path.is_absolute() or path != path.absolute():
        raise ValueError(f"{label} must be an absolute normalized path")
    try:
        path.relative_to(root)
    except ValueError as error:
        raise ValueError(f"{label} must be beneath the Stage I root") from error
    return path


def is_fresh_r12_rerun_profile(profile: dict[str, object]) -> bool:
    """Return whether one profile is the sole authorized fresh R12 rerun shape."""

    parent_keys = (
        "parent_job_id",
        "parent_result",
        "parent_segment",
        "restart_file",
        "restart_file_sha256",
        "restart_time",
    )
    return (
        all(profile.get(key) is None for key in parent_keys)
        and all(profile.get(key) == value for key, value in FRESH_R12_RERUN.items())
    )


def validate_fresh_r12_rerun_evidence(profiles: object, barrier: object) -> None:
    """Require exact historical non-authorizing evidence for a fresh R12 rerun."""

    if not isinstance(profiles, list) or not isinstance(barrier, list):
        raise ValueError("fresh R12 rerun bindings are invalid")
    fresh = [
        profile
        for profile in profiles
        if isinstance(profile, dict)
        and profile.get("case_id") == FRESH_R12_RERUN["case_id"]
        and profile.get("segment") == FRESH_R12_RERUN["segment"]
    ]
    if not fresh:
        return
    if len(fresh) != 1 or not is_fresh_r12_rerun_profile(fresh[0]):
        raise ValueError("fresh R12 rerun profile differs from the exact t=0 decision")
    if barrier.count(HISTORICAL_R12_CLEAN_PARTIAL) != 1:
        raise ValueError(
            "fresh R12 rerun requires exact historical s00/4766856 clean_partial evidence"
        )


def validate_bounded_wave_authorization(value: object, root: Path,
                                        args: argparse.Namespace
                                        ) -> dict[str, object]:
    """Validate one exact generalized V2 authorization."""

    if not isinstance(value, dict):
        raise ValueError("V2 authorization must be an object")
    mode = value.get("mode")
    common_keys = {
        "mode",
        "authorizing",
        "authorized_next_profiles",
        "bounded_concurrency",
        "controller_consumption_state",
    }
    if mode == "bounded-wave":
        authorization = require_exact_keys(
            value,
            {*common_keys, "non_authorizing_reason"},
            "bounded-wave V2 authorization",
        )
        if authorization["authorizing"] is not False:
            raise ValueError("bounded-wave V2 artifact must remain non-authorizing")
        state = require_nonempty_string(
            authorization["controller_consumption_state"],
            "bounded-wave controller consumption state",
        ).lower()
        if "non-authorizing" not in state or "pending" not in state:
            raise ValueError(
                "bounded-wave V2 artifact does not disclose pending controller enforcement"
            )
        require_nonempty_string(
            authorization["non_authorizing_reason"],
            "bounded-wave non-authorizing reason",
        )
    elif mode == "sole-next-profile":
        authorization = require_exact_keys(
            value,
            {*common_keys, "sole_next_segment_profile"},
            "sole-profile V2 authorization",
        )
        if authorization["authorizing"] is not True:
            raise ValueError("sole-profile V2 authorization must be authorizing")
        require_nonempty_string(
            authorization["controller_consumption_state"],
            "sole-profile controller consumption state",
        )
    else:
        raise ValueError("V2 authorization mode differs")

    profiles = authorization["authorized_next_profiles"]
    if not isinstance(profiles, list) or not profiles:
        raise ValueError("V2 authorization profiles must be a nonempty list")
    if mode == "bounded-wave" and not 2 <= len(profiles) <= MAX_BOUNDED_WAVE_LANES:
        raise ValueError(
            f"bounded-wave authorization requires 2-{MAX_BOUNDED_WAVE_LANES} profiles"
        )
    if mode == "sole-next-profile" and len(profiles) != 1:
        raise ValueError("sole-profile V2 authorization requires exactly one profile")
    if (
        mode == "sole-next-profile"
        and authorization["sole_next_segment_profile"] != {
            key: profiles[0][key] for key in SOLE_COMPATIBLE_PROFILE_KEYS
        }
    ):
        raise ValueError("sole-profile V2 authorization profile bindings differ")
    concurrency = require_exact_keys(
        authorization["bounded_concurrency"],
        {"max_active_segments", "max_wave_nodes", "r17_exclusive_and_last"},
        "V2 bounded concurrency",
    )
    max_active = require_integer(
        concurrency["max_active_segments"],
        "V2 max active segments",
        minimum=1,
    )
    if max_active > MAX_BOUNDED_WAVE_LANES or len(profiles) > max_active:
        raise ValueError(
            f"V2 authorization exceeds the {MAX_BOUNDED_WAVE_LANES}-lane ceiling"
        )
    max_nodes = require_integer(
        concurrency["max_wave_nodes"], "V2 max wave nodes", minimum=1
    )
    if max_nodes > MAX_BOUNDED_WAVE_NODES:
        raise ValueError(
            f"V2 authorization exceeds the {MAX_BOUNDED_WAVE_NODES}-node ceiling"
        )
    if concurrency["r17_exclusive_and_last"] is not True:
        raise ValueError("V2 authorization does not preserve R17 exclusive/last")

    source_relative = safe_relative_path(
        args.source_bundle_relative_path, "source bundle relative path"
    )
    source_path = root / source_relative
    source_sha256 = require_sha256(
        args.expected_source_bundle_sha256, "source bundle SHA-256"
    )
    retained_profiles: list[dict[str, object]] = []
    cases: set[str] = set()
    segments: set[tuple[str, str]] = set()
    total_nodes = 0
    for index, item in enumerate(profiles):
        profile = require_exact_keys(item, BOUNDED_PROFILE_KEYS, f"V2 profile {index}")
        case_id = require_safe_id(profile["case_id"], f"V2 profile {index} case")
        if mode == "bounded-wave" and case_id == "R17":
            raise ValueError("R17 is forbidden in bounded-wave mode; it must remain exclusive")
        if case_id != "R17" and LOWER_RESOLUTION_CASE_PATTERN.fullmatch(case_id) is None:
            raise ValueError(f"V2 profile {index} is not a Stage I case")
        if case_id in cases:
            raise ValueError(f"V2 profiles duplicate lane {case_id}")
        cases.add(case_id)
        segment = require_safe_id(profile["segment"], f"V2 profile {index} segment")
        if (
            case_id == HISTORICAL_R12_CLEAN_PARTIAL["case_id"]
            and segment == HISTORICAL_R12_CLEAN_PARTIAL["segment"]
        ):
            raise ValueError(
                "historical R12 s00/4766856 clean_partial evidence is inventory-only"
            )
        identity = (case_id, segment)
        if identity in segments:
            raise ValueError("V2 profiles duplicate a segment")
        segments.add(identity)
        nodes = require_integer(profile["nodes"], f"V2 profile {index} nodes", minimum=1)
        total_nodes += nodes
        slurm_seconds = walltime_seconds(profile["walltime"], f"V2 profile {index} walltime")
        if require_integer(
            profile["controller_walltime_max_seconds"],
            f"V2 profile {index} controller walltime maximum",
            minimum=1,
        ) != MAX_WALLTIME_SECONDS:
            raise ValueError(f"V2 profile {index} controller walltime maximum differs")
        athena_seconds = walltime_seconds(
            profile["athena_walltime"], f"V2 profile {index} Athena walltime"
        )
        if slurm_seconds - athena_seconds < 600:
            raise ValueError(f"V2 profile {index} lacks the ten-minute shutdown margin")
        if require_integer(
            profile["ranks_per_node"], f"V2 profile {index} ranks per node", minimum=1
        ) != 8:
            raise ValueError(f"V2 profile {index} ranks per node differ")
        if require_integer(
            profile["cpus_per_task"], f"V2 profile {index} CPUs per task", minimum=1
        ) != 7:
            raise ValueError(f"V2 profile {index} CPUs per task differ")
        target = require_finite_float(
            profile["time_tlim_target"], f"V2 profile {index} target"
        )
        if target <= 0 or target > MAX_STAGE_I_TIME:
            raise ValueError(f"V2 profile {index} target is outside Stage I")
        require_integer(
            profile["estimated_storage_bytes"],
            f"V2 profile {index} estimated storage bytes",
            minimum=1,
        )
        if profile["source_bundle"] != str(source_path):
            raise ValueError(f"V2 profile {index} source bundle path differs")
        if profile["source_bundle_sha256"] != source_sha256:
            raise ValueError(f"V2 profile {index} source bundle SHA-256 differs")
        if profile["output_layout"] != "rank-local":
            raise ValueError(f"V2 profile {index} output layout differs")
        require_nonempty_string(
            profile["acceptance_policy"], f"V2 profile {index} acceptance policy"
        )
        require_nonempty_string(
            profile["acceptance_criterion"], f"V2 profile {index} acceptance criterion"
        )
        for key in ("executable_sha256", "build_manifest_sha256", "input_sha256"):
            require_sha256(profile[key], f"V2 profile {index} {key}")
        for key in ("executable_revision", "input_revision"):
            if (
                not isinstance(profile[key], str)
                or REVISION_PATTERN.fullmatch(profile[key]) is None
            ):
                raise ValueError(f"V2 profile {index} {key} is invalid")

        parent_values = [
            profile["parent_job_id"],
            profile["parent_result"],
            profile["parent_segment"],
            profile["restart_file"],
            profile["restart_file_sha256"],
            profile["restart_time"],
        ]
        if all(item is None for item in parent_values):
            if not segment.startswith("s00_") and not is_fresh_r12_rerun_profile(profile):
                raise ValueError(f"fresh V2 profile {index} must use an s00 segment")
        elif any(item is None for item in parent_values):
            raise ValueError(f"V2 profile {index} has incomplete parent bindings")
        else:
            require_safe_id(profile["parent_job_id"], f"V2 profile {index} parent job")
            require_safe_id(profile["parent_result"], f"V2 profile {index} parent result")
            require_safe_id(profile["parent_segment"], f"V2 profile {index} parent segment")
            if (
                case_id == HISTORICAL_R12_CLEAN_PARTIAL["case_id"]
                and (
                    profile["parent_job_id"] == HISTORICAL_R12_CLEAN_PARTIAL["job_id"]
                    or profile["parent_segment"]
                    == HISTORICAL_R12_CLEAN_PARTIAL["segment"]
                )
            ):
                raise ValueError(
                    "historical R12 s00/4766856 clean_partial evidence is non-authorizing"
                )
            require_root_member(root, profile["restart_file"], f"V2 profile {index} restart")
            require_sha256(profile["restart_file_sha256"], f"V2 profile {index} restart SHA-256")
            restart_time = require_finite_float(
                profile["restart_time"], f"V2 profile {index} restart time"
            )
            if restart_time < 0 or target <= restart_time:
                raise ValueError(f"V2 profile {index} restart/target times differ")
        retained_profiles.append(dict(profile))

    if profiles != sorted(
        retained_profiles, key=lambda item: (str(item["case_id"]), str(item["segment"]))
    ):
        raise ValueError("V2 profiles are not in canonical case/segment order")
    if total_nodes != max_nodes:
        raise ValueError("V2 profile nodes differ from exact max_wave_nodes")
    if "R17" in cases and (mode != "sole-next-profile" or len(profiles) != 1 or max_nodes != 8):
        raise ValueError("R17 V2 authorization must be exclusive on eight nodes")
    return authorization


def validate_recost_recommendations(value: object, root: Path,
                                    args: argparse.Namespace
                                    ) -> dict[str, object]:
    """Validate one exact schema-2 non-authorizing recommendation packet."""

    if not isinstance(value, dict):
        raise ValueError("schema-2 recost recommendations must be an object")
    mode = value.get("mode")
    keys = {
        "mode",
        "authorizing",
        "recommended_next_profiles",
        "bounded_concurrency",
        "controller_consumption_state",
        "non_authorizing_reason",
    }
    if mode == "sole-next-profile":
        keys.add("sole_next_segment_recommendation")
    elif mode != "bounded-wave":
        raise ValueError("schema-2 recost recommendation mode differs")
    recommendations = require_exact_keys(value, keys, "schema-2 recost recommendations")
    if recommendations["authorizing"] is not False:
        raise ValueError("schema-2 recost recommendations must remain non-authorizing")
    state = require_nonempty_string(
        recommendations["controller_consumption_state"],
        "schema-2 controller consumption state",
    ).lower()
    if "non-authorizing" not in state:
        raise ValueError("schema-2 controller consumption state must remain non-authorizing")
    require_nonempty_string(
        recommendations["non_authorizing_reason"],
        "schema-2 non-authorizing reason",
    )
    profiles = recommendations["recommended_next_profiles"]
    if not isinstance(profiles, list) or not profiles:
        raise ValueError("schema-2 recommended profiles must be a nonempty list")
    projected_profiles = []
    for index, item in enumerate(profiles):
        profile = require_exact_keys(
            item,
            {*BOUNDED_PROFILE_KEYS, "recommendation_basis"},
            f"schema-2 recommended profile {index}",
        )
        basis = profile["recommendation_basis"]
        if not isinstance(basis, dict) or not basis:
            raise ValueError(f"schema-2 recommended profile {index} basis is invalid")
        projected_profiles.append(
            {key: profile[key] for key in BOUNDED_PROFILE_KEYS}
        )
    projected = {
        "mode": mode,
        "authorizing": False if mode == "bounded-wave" else True,
        "authorized_next_profiles": projected_profiles,
        "bounded_concurrency": recommendations["bounded_concurrency"],
        "controller_consumption_state": recommendations["controller_consumption_state"],
    }
    if mode == "bounded-wave":
        projected["non_authorizing_reason"] = recommendations["non_authorizing_reason"]
    else:
        compatible = {
            key: profiles[0][key] for key in SOLE_COMPATIBLE_PROFILE_KEYS
        }
        if recommendations["sole_next_segment_recommendation"] != compatible:
            raise ValueError("schema-2 sole recommendation compatibility binding differs")
        projected["sole_next_segment_profile"] = compatible
    validate_bounded_wave_authorization(projected, root, args)
    return recommendations


def validate_bounded_scheduler_evidence(value: object) -> list[dict[str, object]]:
    """Validate the exact sorted bounded-wave scheduler-evidence list."""

    if not isinstance(value, list) or not value:
        raise ValueError("bounded-wave scheduler evidence must be a nonempty list")
    retained: list[dict[str, object]] = []
    paths: set[str] = set()
    jobs: set[str] = set()
    for index, item in enumerate(value):
        evidence = require_exact_keys(
            item, SCHEDULER_EVIDENCE_KEYS, f"bounded scheduler evidence {index}"
        )
        relative = safe_relative_path(
            require_nonempty_string(
                evidence["path"], f"bounded scheduler evidence {index} path"
            ),
            f"bounded scheduler evidence {index} path",
        )
        if relative.parent != Path("accounting"):
            raise ValueError("bounded scheduler evidence must be a direct child of accounting")
        if relative.as_posix() in paths:
            raise ValueError("bounded scheduler evidence duplicates a path")
        paths.add(relative.as_posix())
        job_id = require_safe_id(
            evidence["job_id"], f"bounded scheduler evidence {index} job ID"
        )
        if job_id in jobs:
            raise ValueError(f"bounded scheduler evidence duplicates job {job_id}")
        jobs.add(job_id)
        job_name = require_safe_id(
            evidence["job_name"], f"bounded scheduler evidence {index} job name"
        )
        if not job_name.startswith(CGL_JOB_NAME_PREFIX):
            raise ValueError("bounded scheduler evidence job name is not a CGL job")
        state = require_safe_id(
            evidence["state"], f"bounded scheduler evidence {index} state"
        )
        if state not in TERMINAL_SCHEDULER_STATES:
            raise ValueError(f"bounded scheduler evidence {index} state is not terminal")
        exit_code = require_nonempty_string(
            evidence["exit_code"], f"bounded scheduler evidence {index} exit code"
        )
        if EXIT_CODE_PATTERN.fullmatch(exit_code) is None:
            raise ValueError(f"bounded scheduler evidence {index} exit code is invalid")
        require_sha256(
            evidence["sha256"], f"bounded scheduler evidence {index} SHA-256"
        )
        require_integer(
            evidence["nodes"], f"bounded scheduler evidence {index} nodes", minimum=1
        )
        submitted = parse_utc_timestamp(
            evidence["submitted_utc"],
            f"bounded scheduler evidence {index} submit time",
        )
        completed = parse_utc_timestamp(
            evidence["completed_utc"],
            f"bounded scheduler evidence {index} completion time",
        )
        if completed <= submitted:
            raise ValueError(
                f"bounded scheduler evidence {index} chronology is invalid"
            )
        elapsed = require_integer(
            evidence["elapsed_seconds"],
            f"bounded scheduler evidence {index} elapsed seconds",
        )
        if elapsed > int((completed - submitted).total_seconds()) + int(
            SCHEDULER_TIME_TOLERANCE.total_seconds()
        ):
            raise ValueError(
                f"bounded scheduler evidence {index} elapsed time exceeds its chronology"
            )
        retained.append(dict(evidence))
    return retained


def validate_bounded_barrier(value: object,
                             scheduler: list[dict[str, object]]) -> None:
    """Bind every scheduler evidence item to one exact barrier segment."""

    if not isinstance(value, list) or not value:
        raise ValueError("bounded-wave barrier must contain recorded segments")
    barriers: dict[str, dict[str, object]] = {}
    for index, item in enumerate(value):
        barrier = require_exact_keys(
            item,
            {"case_id", "segment", "job_id", "result"},
            f"bounded barrier segment {index}",
        )
        job_id = require_safe_id(barrier["job_id"], f"bounded barrier segment {index} job")
        if job_id in barriers:
            raise ValueError(f"bounded barrier duplicates job {job_id}")
        require_safe_id(barrier["case_id"], f"bounded barrier segment {index} case")
        require_safe_id(barrier["segment"], f"bounded barrier segment {index} segment")
        require_safe_id(barrier["result"], f"bounded barrier segment {index} result")
        barriers[job_id] = barrier
    barrier_jobs = [str(item["job_id"]) for item in value]
    scheduler_jobs = [str(item["job_id"]) for item in scheduler]
    if barrier_jobs != scheduler_jobs:
        raise ValueError("bounded scheduler evidence job order differs from the barrier")
    for item in scheduler:
        barrier = barriers[str(item["job_id"])]
        expected_name = (
            f"cgl_mks24_{EXECUTION_EPOCH_SLUG}_"
            f"{barrier['case_id']}_{barrier['segment']}"
        )
        if item["job_name"] != expected_name:
            raise ValueError(
                f"bounded scheduler evidence job {item['job_id']} name differs from barrier"
            )


def require_v2_input_binding(value: object, label: str) -> dict[str, str]:
    """Require one exact root-relative path and digest binding."""

    binding = require_exact_keys(value, {"path", "sha256"}, label)
    path = safe_relative_path(
        require_nonempty_string(binding["path"], f"{label} path"),
        f"{label} path",
    ).as_posix()
    return {
        "path": path,
        "sha256": require_sha256(binding["sha256"], f"{label} SHA-256"),
    }


def validate_f118_source_authority_binding(
    value: object, root: Path, args: argparse.Namespace
) -> dict[str, object]:
    """Require the finalized F118 binding and exact committed seven-tool vector."""

    retained = require_exact_keys(
        value,
        {
            "checkpoint",
            "evidence",
            "provenance_review",
            "plasma_review",
            "publication_audit",
            "final_source_bundle",
        },
        "schema-2 F118 source authority binding",
    )
    if retained["checkpoint"] != "F-118":
        raise ValueError("schema-2 current source authority binding is not F-118")
    paths = {
        "evidence": F118_RELATIVE,
        "provenance_review": F118_PROVENANCE_REVIEW_RELATIVE,
        "plasma_review": F118_PLASMA_REVIEW_RELATIVE,
        "publication_audit": F118_PUBLICATION_AUDIT_RELATIVE,
    }
    normalized: dict[str, object] = {"checkpoint": "F-118"}
    payloads: dict[str, bytes] = {}
    digests: dict[str, str] = {}
    for key, expected_path in paths.items():
        binding = require_v2_input_binding(
            retained[key], f"schema-2 F118 {key} binding"
        )
        if binding["path"] != expected_path.as_posix():
            raise ValueError(f"schema-2 F118 {key} path differs")
        _, payload = read_confined_file_sha256(
            root,
            binding["path"],
            binding["sha256"],
            f"schema-2 F118 {key}",
            expected_mode=0o444,
        )
        payloads[key] = payload
        digests[key] = binding["sha256"]
        normalized[key] = binding
    final = require_exact_keys(
        retained["final_source_bundle"],
        {"path", "sha256", "verified_revisions"},
        "schema-2 F118 final source bundle binding",
    )
    final_binding = {
        "path": safe_relative_path(
            require_nonempty_string(
                final["path"], "schema-2 F118 final source bundle path"
            ),
            "schema-2 F118 final source bundle path",
        ).as_posix(),
        "sha256": require_sha256(
            final["sha256"], "schema-2 F118 final source bundle SHA-256"
        ),
        "verified_revisions": final["verified_revisions"],
    }
    if final_binding != {
        "path": args.source_bundle_relative_path,
        "sha256": args.expected_source_bundle_sha256,
        "verified_revisions": args.expected_source_bundle_verified_revisions_json,
    }:
        raise ValueError("schema-2 F118 final source bundle binding differs")
    historical = validate_f118_committed_tools(payloads["evidence"], final_binding, root)
    validate_complete_f118_source_authority(
        payloads, digests, final_binding, root, historical
    )
    normalized["final_source_bundle"] = final_binding
    return normalized


def git_revision_subject(repository: Path, revision: str, label: str) -> str:
    """Return one exact retained Git commit subject."""

    completed = git_run(
        repository, ["show", "-s", "--format=%s", revision], capture_output=True
    )
    if completed.returncode:
        raise ValueError(f"cannot inspect {label} subject")
    try:
        subject = completed.stdout.decode("utf-8").rstrip("\n")
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} subject is not UTF-8") from error
    return require_nonempty_string(subject, f"{label} subject")


def validate_source_authority_tool_vector(
    tools: object,
    publisher_value: object,
    head: str,
    repository: Path,
    label: str,
) -> list[dict[str, object]]:
    """Authenticate the exact seven-tool source-authority contract at one revision."""

    expected_paths = sorted(F118_REQUIRED_TOOLS)
    if not isinstance(tools, list) or len(tools) != len(expected_paths):
        raise ValueError(f"{label} committed_tools must contain exactly seven tools")
    normalized = []
    for index, expected_path in enumerate(expected_paths):
        tool = require_exact_keys(
            tools[index],
            {"path", "revision", "sha256", "mode"},
            f"{label} committed tool {index}",
        )
        expected_mode = F118_REQUIRED_TOOLS[expected_path]
        if (
            tool["path"] != expected_path
            or tool["revision"] != head
            or tool["mode"] != expected_mode
        ):
            raise ValueError(
                f"{label} committed tool {index} path, revision, or mode differs"
            )
        digest = require_sha256(tool["sha256"], f"{label} committed tool {index} SHA-256")
        require_committed_revision_bytes(
            repository,
            Path(expected_path),
            head,
            digest,
            f"{label} committed tool {index}",
        )
        tree = git_run(
            repository,
            ["ls-tree", head, "--", expected_path],
            capture_output=True,
        )
        expected_git_mode = "100755" if expected_mode == "0755" else "100644"
        try:
            tree_fields = tree.stdout.decode("ascii").strip().split()
        except UnicodeDecodeError as error:
            raise ValueError(
                f"{label} committed tool {index} tree mode is not ASCII"
            ) from error
        if (
            tree.returncode
            or len(tree_fields) != 4
            or tree_fields[0] != expected_git_mode
            or tree_fields[1] != "blob"
            or tree_fields[3] != expected_path
        ):
            raise ValueError(f"{label} committed tool {index} Git mode differs")
        normalized.append(
            {
                "path": expected_path,
                "revision": head,
                "sha256": digest,
                "mode": expected_mode,
            }
        )
    publisher = require_exact_keys(
        publisher_value,
        {"path", "revision", "sha256", "mode"},
        f"{label} publisher",
    )
    publisher_tool = normalized[expected_paths.index(F118_PUBLISHER_RELATIVE)]
    if publisher != publisher_tool:
        raise ValueError(f"{label} publisher differs from committed_tools")
    return normalized


def validate_f118_committed_tools(
    payload: bytes, final_binding: dict[str, object], root: Path
) -> dict[str, object]:
    """Validate F118 current tools and its immutable exact F116 predecessor."""

    try:
        evidence = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("schema-2 F118 evidence is not valid JSON") from error
    if (json.dumps(evidence, indent=2, sort_keys=True) + "\n").encode() != payload:
        raise ValueError("schema-2 F118 evidence must be stable canonical JSON")
    retained = require_exact_keys(
        evidence,
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "execution_epoch",
            "generated_utc",
            "scope",
            "predecessor_authorities",
            "implementation",
            "source_archive_catalog",
            "authorization",
            "validation",
            "publication_requirements",
        },
        "schema-2 F118 evidence",
    )
    if (
        retained["schema_version"] != 1
        or retained["record_type"]
        != "stage-i-current-source-authority-supersession-evidence"
        or retained["checkpoint"] != "F-118"
        or retained["execution_epoch"] != EXECUTION_EPOCH
    ):
        raise ValueError("schema-2 F118 evidence identity differs")
    predecessors = require_exact_keys(
        retained["predecessor_authorities"],
        {"historical_f116"},
        "schema-2 F118 predecessor authorities",
    )
    implementation = require_exact_keys(
        retained["implementation"],
        {
            "publisher",
            "committed_tools",
            "intermediate_36140_bundle",
            "predecessor_current_source_bundle",
            "current_source_bundle",
        },
        "schema-2 F118 implementation",
    )
    current = require_exact_keys(
        implementation["current_source_bundle"],
        {
            "path",
            "sha256",
            "complete_history",
            "head",
            "advertised_tip",
            "verified_revisions",
            "selected_as_current",
            "candidate_path",
            "subject",
        },
        "schema-2 F118 current source bundle",
    )
    head = require_nonempty_string(current["head"], "schema-2 F118 final HEAD")
    if REVISION_PATTERN.fullmatch(head) is None:
        raise ValueError("schema-2 F118 final HEAD is invalid")
    if {
        "path": current["path"],
        "sha256": current["sha256"],
        "verified_revisions": current["verified_revisions"],
    } != final_binding:
        raise ValueError("schema-2 F118 evidence final source bundle differs")

    historical_bindings = require_exact_keys(
        predecessors["historical_f116"],
        {"evidence", "publication_audit", "provenance_review", "plasma_review"},
        "schema-2 F118 historical F116 bindings",
    )
    historical_payload = b""
    for key, expected_path in {
        "evidence": F116_RELATIVE,
        "publication_audit": F116_PUBLICATION_AUDIT_RELATIVE,
        "provenance_review": F116_PROVENANCE_REVIEW_RELATIVE,
        "plasma_review": F116_PLASMA_REVIEW_RELATIVE,
    }.items():
        binding = require_v2_input_binding(
            historical_bindings[key], f"schema-2 historical F116 {key}"
        )
        if binding["path"] != expected_path.as_posix():
            raise ValueError(f"schema-2 historical F116 {key} path differs")
        if root == DEFAULT_ROOT and binding["sha256"] != F116_CANONICAL_SHA256[key]:
            raise ValueError(f"schema-2 historical F116 {key} canonical digest differs")
        _, retained_payload = read_confined_file_sha256(
            root,
            binding["path"],
            binding["sha256"],
            f"schema-2 historical F116 {key}",
            expected_mode=0o444,
        )
        if key == "evidence":
            historical_payload = retained_payload
    try:
        historical = json.loads(historical_payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("schema-2 historical F116 evidence is not valid JSON") from error
    historical = require_exact_keys(
        historical,
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "execution_epoch",
            "generated_utc",
            "scope",
            "predecessor_authorities",
            "implementation",
            "source_archive_catalog",
            "authorization",
            "validation",
            "publication_requirements",
        },
        "schema-2 historical F116 evidence",
    )
    if (
        historical["schema_version"] != 1
        or historical["record_type"]
        != "stage-i-current-source-authority-supersession-evidence"
        or historical["checkpoint"] != "F-116"
        or historical["execution_epoch"] != EXECUTION_EPOCH
    ):
        raise ValueError("schema-2 historical F116 evidence identity differs")
    historical_implementation = require_exact_keys(
        historical["implementation"],
        {
            "publisher",
            "committed_tools",
            "intermediate_36140_bundle",
            "current_source_bundle",
        },
        "schema-2 historical F116 implementation",
    )
    historical_current = require_exact_keys(
        historical_implementation["current_source_bundle"],
        {
            "path",
            "sha256",
            "complete_history",
            "head",
            "advertised_tip",
            "verified_revisions",
            "selected_as_current",
            "candidate_path",
            "subject",
        },
        "schema-2 historical F116 current source bundle",
    )
    historical_bridge = require_exact_keys(
        historical_implementation["intermediate_36140_bundle"],
        {
            "path",
            "sha256",
            "complete_history",
            "head",
            "advertised_tip",
            "verified_revisions",
            "selected_as_current",
            "role",
        },
        "schema-2 historical F116 bridge bundle",
    )
    if implementation["intermediate_36140_bundle"] != historical_bridge:
        raise ValueError("schema-2 F118 bridge does not exactly preserve F116")
    predecessor = require_exact_keys(
        implementation["predecessor_current_source_bundle"],
        {
            "path",
            "sha256",
            "complete_history",
            "head",
            "advertised_tip",
            "verified_revisions",
            "selected_as_current",
            "role",
            "subject",
        },
        "schema-2 F118 predecessor current source bundle",
    )
    expected_predecessor = dict(historical_current)
    expected_predecessor.pop("candidate_path")
    expected_predecessor["selected_as_current"] = False
    expected_predecessor["role"] = "retained-non-current-predecessor"
    if predecessor != expected_predecessor:
        raise ValueError("schema-2 F118 predecessor does not exactly preserve F116")

    repository = repository_root(initial_source_path())
    validate_source_authority_tool_vector(
        historical_implementation["committed_tools"],
        historical_implementation["publisher"],
        require_nonempty_string(
            historical_current["head"], "schema-2 historical F116 current HEAD"
        ),
        repository,
        "schema-2 historical F116",
    )
    validate_source_authority_tool_vector(
        implementation["committed_tools"],
        implementation["publisher"],
        head,
        repository,
        "schema-2 F118",
    )
    if current["subject"] != git_revision_subject(repository, head, "schema-2 F118 final HEAD"):
        raise ValueError("schema-2 F118 final HEAD subject differs")
    revisions = current["verified_revisions"]
    if (
        not isinstance(revisions, list)
        or not revisions
        or len(set(revisions)) != len(revisions)
        or any(
            not isinstance(revision, str) or REVISION_PATTERN.fullmatch(revision) is None
            for revision in revisions
        )
    ):
        raise ValueError("schema-2 F118 final verified revisions differ")
    required_revisions = {
        head,
        require_nonempty_string(historical_current["head"], "historical F116 current HEAD"),
        require_nonempty_string(historical_bridge["head"], "historical F116 bridge HEAD"),
        *historical_current["verified_revisions"],
        *historical_bridge["verified_revisions"],
    }
    if not required_revisions.issubset(set(revisions)):
        raise ValueError("schema-2 F118 final verified revisions omit F116 history")
    historical_catalog = require_exact_keys(
        historical["source_archive_catalog"],
        {"before", "after"},
        "schema-2 historical F116 source-archive catalog",
    )
    historical_after = require_exact_keys(
        historical_catalog["after"],
        {
            "readme_sha256",
            "sha256sums_sha256",
            "bridge_listed_exactly_once",
            "final_bundle_listed_exactly_once",
            "corrupt_c7_listed",
            "historical_f115_preserved",
            "sole_current_source_bundle",
        },
        "schema-2 historical F116 source-archive catalog after",
    )
    for key in ("readme_sha256", "sha256sums_sha256"):
        require_sha256(historical_after[key], f"schema-2 historical F116 catalog {key}")
    if (
        historical_after["bridge_listed_exactly_once"] is not True
        or historical_after["final_bundle_listed_exactly_once"] is not True
        or historical_after["corrupt_c7_listed"] is not False
        or historical_after["historical_f115_preserved"] is not True
        or historical_after["sole_current_source_bundle"] != historical_current["path"]
    ):
        raise ValueError("schema-2 historical F116 source-archive catalog differs")
    return {
        "bridge": historical_bridge,
        "current": historical_current,
        "catalog_after": historical_after,
    }


def canonical_json_object(payload: bytes, label: str) -> dict[str, object]:
    """Parse one stable canonical JSON object."""

    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON") from error
    if (
        not isinstance(value, dict)
        or (json.dumps(value, indent=2, sort_keys=True) + "\n").encode() != payload
    ):
        raise ValueError(f"{label} must be stable canonical JSON")
    return value


def exact_publication_binding(
    value: object,
    path: Path,
    digest: str,
    label: str,
    *,
    mode: str,
) -> dict[str, object]:
    """Require one exact single-link publication binding."""

    retained = require_exact_keys(
        value, {"path", "sha256", "mode", "links"}, label
    )
    if retained != {
        "path": str(path),
        "sha256": digest,
        "mode": mode,
        "links": 1,
    }:
        raise ValueError(f"{label} differs")
    return retained


def validate_f118_catalog_snapshot(value: object, label: str, *, after: bool
                                   ) -> dict[str, object]:
    """Require the exact F118 before/after source-catalog claim schema."""

    keys = {
        "readme_sha256",
        "sha256sums_sha256",
        "bridge_listed_exactly_once",
        "predecessor_current_source_bundle_listed_exactly_once",
        "corrupt_c7_listed",
        "historical_f115_preserved",
    }
    if after:
        keys |= {
            "final_bundle_listed_exactly_once",
            "historical_f116_preserved",
            "all_prior_checksum_entries_preserved",
            "sole_current_source_bundle",
        }
    else:
        keys.add("final_bundle_listed")
    retained = require_exact_keys(value, keys, label)
    require_sha256(retained["readme_sha256"], f"{label} README SHA-256")
    require_sha256(retained["sha256sums_sha256"], f"{label} SHA256SUMS SHA-256")
    expected = {
        "bridge_listed_exactly_once": True,
        "predecessor_current_source_bundle_listed_exactly_once": True,
        "corrupt_c7_listed": False,
        "historical_f115_preserved": True,
    }
    if after:
        expected |= {
            "final_bundle_listed_exactly_once": True,
            "historical_f116_preserved": True,
            "all_prior_checksum_entries_preserved": True,
        }
        require_nonempty_string(
            retained["sole_current_source_bundle"],
            f"{label} sole current source bundle",
        )
    else:
        expected["final_bundle_listed"] = False
    if any(retained[key] != expected_value for key, expected_value in expected.items()):
        raise ValueError(f"{label} authority claims differ")
    return retained


def parse_source_archive_sha256sums(payload: bytes, label: str) -> list[tuple[str, str]]:
    """Parse the exact structured source-archive checksum ledger."""

    if not payload.endswith(b"\n"):
        raise ValueError(f"{label} must end with newline")
    try:
        lines = payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} is not UTF-8") from error
    entries: list[tuple[str, str]] = []
    names: set[str] = set()
    for line in lines:
        match = re.fullmatch(r"([0-9a-f]{64})  ([A-Za-z0-9_.-]+)", line)
        if match is None:
            raise ValueError(f"{label} row is malformed: {line}")
        digest, name = match.groups()
        if name in names:
            raise ValueError(f"{label} duplicates {name}")
        names.add(name)
        entries.append((digest, name))
    if not entries:
        raise ValueError(f"{label} is empty")
    return entries


def f118_catalog_readme_block(current: dict[str, object]) -> bytes:
    """Return the exact source-authority F118 README append block."""

    name = Path(str(current["path"])).name
    return (
        f"`{name}` records complete history through commit `{current['head']}` "
        f"(`{current['subject']}`). It is the sole current Stage I source-selection "
        "bundle after independently reviewed F-118 publication. It preserves the "
        "immutable four-part F-116 authority and every prior active source-archive "
        "catalog entry; it does not itself authorize prepare or submission. "
        f"Its SHA-256 is `{current['sha256']}`.\n\n"
    ).encode()


def validate_f118_catalog_transition(
    root: Path,
    current: dict[str, object],
    bridge: dict[str, object],
    predecessor: dict[str, object],
    before: dict[str, object],
    after: dict[str, object],
) -> None:
    """Prove the live catalog is exactly F116 plus the deterministic F118 append."""

    _, readme = read_confined_file_sha256(
        root,
        "source-archives/README.md",
        str(after["readme_sha256"]),
        "schema-2 F118 live source-archive README",
        expected_mode=0o644,
    )
    _, sums = read_confined_file_sha256(
        root,
        "source-archives/SHA256SUMS",
        str(after["sha256sums_sha256"]),
        "schema-2 F118 live source-archive SHA256SUMS",
        expected_mode=0o644,
    )
    marker = b"## AthenaK\n\n"
    block = f118_catalog_readme_block(current)
    if readme.count(marker) != 1 or readme.count(block) != 1:
        raise ValueError("schema-2 F118 live README lacks the exact single append block")
    old_readme = readme.replace(marker + block, marker, 1)
    if (
        old_readme == readme
        or sha256_bytes(old_readme) != before["readme_sha256"]
    ):
        raise ValueError("schema-2 F118 live README does not derive from F116 catalog_after")
    bridge_name = Path(str(bridge["path"])).name
    predecessor_name = Path(str(predecessor["path"])).name
    current_name = Path(str(current["path"])).name
    if (
        bridge_name.encode() not in old_readme
        or predecessor_name.encode() not in old_readme
        or current_name.encode() in old_readme
    ):
        raise ValueError("schema-2 F118 predecessor README entries differ")

    entries = parse_source_archive_sha256sums(sums, "schema-2 F118 live SHA256SUMS")
    expected_final = (str(current["sha256"]), current_name)
    if entries[-1] != expected_final:
        raise ValueError("schema-2 F118 current bundle is not the exact final checksum entry")
    final_line = f"{expected_final[0]}  {expected_final[1]}\n".encode()
    if not sums.endswith(final_line):
        raise ValueError("schema-2 F118 checksum-ledger append bytes differ")
    old_sums = sums[: -len(final_line)]
    if sha256_bytes(old_sums) != before["sha256sums_sha256"]:
        raise ValueError("schema-2 F118 live SHA256SUMS does not derive from F116 catalog_after")
    old_entries = parse_source_archive_sha256sums(
        old_sums, "schema-2 F118 predecessor SHA256SUMS"
    )
    old_names = [name for _, name in old_entries]
    if (
        old_names.count(bridge_name) != 1
        or old_names.count(predecessor_name) != 1
        or current_name in old_names
    ):
        raise ValueError("schema-2 F118 predecessor checksum entries differ")
    for digest, name in old_entries:
        read_confined_file_sha256(
            root,
            f"source-archives/{name}",
            digest,
            f"schema-2 F118 preserved source archive {name}",
            expected_mode=0o644,
        )


def validate_complete_f118_source_authority(
    payloads: dict[str, bytes],
    digests: dict[str, str],
    final_binding: dict[str, object],
    root: Path,
    historical: dict[str, object],
) -> None:
    """Validate the complete F118 evidence, independent reviews, and audit."""

    if set(payloads) != {"evidence", "provenance_review", "plasma_review", "publication_audit"}:
        raise ValueError("schema-2 F118 retained record set is incomplete")
    if set(digests) != set(payloads):
        raise ValueError("schema-2 F118 retained digest set is incomplete")
    loaded = {
        key: canonical_json_object(payload, f"schema-2 F118 {key}")
        for key, payload in payloads.items()
    }
    evidence = require_exact_keys(
        loaded["evidence"],
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "execution_epoch",
            "generated_utc",
            "scope",
            "predecessor_authorities",
            "implementation",
            "source_archive_catalog",
            "authorization",
            "validation",
            "publication_requirements",
        },
        "schema-2 F118 evidence",
    )
    generated = parse_utc_timestamp(evidence["generated_utc"], "schema-2 F118 generation")
    if (
        evidence["schema_version"] != 1
        or evidence["record_type"]
        != "stage-i-current-source-authority-supersession-evidence"
        or evidence["checkpoint"] != "F-118"
        or evidence["execution_epoch"] != EXECUTION_EPOCH
        or evidence["authorization"] != F118_AUTHORIZATION
        or evidence["validation"] != F118_VALIDATION_CLAIMS
        or evidence["publication_requirements"] != F118_PUBLICATION_REQUIREMENTS
    ):
        raise ValueError("schema-2 F118 evidence identity or authority differs")
    scope = require_exact_keys(
        evidence["scope"],
        {"relationship", "summary", "preserves", "does_not_authorize"},
        "schema-2 F118 scope",
    )
    if (
        scope["relationship"] != "current-source-selection-only-supersession"
        or not require_nonempty_string(scope["summary"], "schema-2 F118 scope summary")
        or scope["preserves"] != F118_PRESERVES
        or scope["does_not_authorize"] != F118_DOES_NOT_AUTHORIZE
    ):
        raise ValueError("schema-2 F118 scope differs or broadens authority")
    implementation = require_exact_keys(
        evidence["implementation"],
        {
            "publisher",
            "committed_tools",
            "intermediate_36140_bundle",
            "predecessor_current_source_bundle",
            "current_source_bundle",
        },
        "schema-2 F118 implementation",
    )
    current = require_exact_keys(
        implementation["current_source_bundle"],
        {
            "path",
            "sha256",
            "complete_history",
            "head",
            "advertised_tip",
            "verified_revisions",
            "selected_as_current",
            "candidate_path",
            "subject",
        },
        "schema-2 F118 current source bundle",
    )
    current_head = require_nonempty_string(current["head"], "schema-2 F118 current head")
    if (
        current["complete_history"] is not True
        or current["selected_as_current"] is not True
        or current_head not in current["verified_revisions"]
        or not require_nonempty_string(current["candidate_path"], "schema-2 F118 candidate")
        or not require_nonempty_string(current["subject"], "schema-2 F118 subject")
        or require_exact_keys(
            current["advertised_tip"],
            {"revision", "name"},
            "schema-2 F118 advertised tip",
        )
        != {"revision": current_head, "name": "HEAD"}
        or {
            "path": current["path"],
            "sha256": current["sha256"],
            "verified_revisions": current["verified_revisions"],
        }
        != final_binding
    ):
        raise ValueError("schema-2 F118 current source bundle authority differs")
    bridge = require_exact_keys(
        implementation["intermediate_36140_bundle"],
        {
            "path",
            "sha256",
            "complete_history",
            "head",
            "advertised_tip",
            "verified_revisions",
            "selected_as_current",
            "role",
        },
        "schema-2 F118 bridge bundle",
    )
    predecessor = require_exact_keys(
        implementation["predecessor_current_source_bundle"],
        {
            "path",
            "sha256",
            "complete_history",
            "head",
            "advertised_tip",
            "verified_revisions",
            "selected_as_current",
            "role",
            "subject",
        },
        "schema-2 F118 predecessor source bundle",
    )
    for bundle, role, tip_name, label in (
        (
            bridge,
            "retained-non-current-bridge",
            "refs/heads/feature/cgl-landau-fluid",
            "bridge",
        ),
        (predecessor, "retained-non-current-predecessor", "HEAD", "predecessor"),
    ):
        if (
            bundle["complete_history"] is not True
            or bundle["selected_as_current"] is not False
            or bundle["role"] != role
            or require_exact_keys(
                bundle["advertised_tip"],
                {"revision", "name"},
                f"schema-2 F118 {label} advertised tip",
            )
            != {"revision": bundle["head"], "name": tip_name}
        ):
            raise ValueError(f"schema-2 F118 {label} source bundle authority differs")
        require_sha256(bundle["sha256"], f"schema-2 F118 {label} bundle SHA-256")
    catalog = require_exact_keys(
        evidence["source_archive_catalog"],
        {"before", "after"},
        "schema-2 F118 evidence source-archive catalog",
    )
    before = validate_f118_catalog_snapshot(
        catalog["before"], "schema-2 F118 catalog before", after=False
    )
    after = validate_f118_catalog_snapshot(
        catalog["after"], "schema-2 F118 catalog after", after=True
    )
    historical_after = require_exact_keys(
        historical["catalog_after"],
        {
            "readme_sha256",
            "sha256sums_sha256",
            "bridge_listed_exactly_once",
            "final_bundle_listed_exactly_once",
            "corrupt_c7_listed",
            "historical_f115_preserved",
            "sole_current_source_bundle",
        },
        "schema-2 authenticated historical F116 catalog_after",
    )
    expected_before = {
        "readme_sha256": historical_after["readme_sha256"],
        "sha256sums_sha256": historical_after["sha256sums_sha256"],
        "bridge_listed_exactly_once": historical_after["bridge_listed_exactly_once"],
        "predecessor_current_source_bundle_listed_exactly_once": (
            historical_after["final_bundle_listed_exactly_once"]
        ),
        "final_bundle_listed": False,
        "corrupt_c7_listed": historical_after["corrupt_c7_listed"],
        "historical_f115_preserved": historical_after["historical_f115_preserved"],
    }
    expected_predecessor = dict(historical["current"])
    expected_predecessor.pop("candidate_path")
    expected_predecessor["selected_as_current"] = False
    expected_predecessor["role"] = "retained-non-current-predecessor"
    if (
        before != expected_before
        or bridge != historical["bridge"]
        or predecessor != expected_predecessor
        or after["sole_current_source_bundle"] != current["path"]
    ):
        raise ValueError("schema-2 F118 catalog predecessor authority differs")
    validate_f118_catalog_transition(root, current, bridge, predecessor, before, after)

    verified = {
        "authorization_broadening": False,
        "bridge_selected_as_current": False,
        "predecessor_current_source_bundle_selected_as_current": False,
        "corrupt_c7_excluded": True,
        "current_source_selection_only": True,
        "final_bundle_sha256": current["sha256"],
        "final_head": current_head,
        "historical_f115_preserved": True,
        "historical_f116_preserved": True,
    }
    reviewers: set[str] = set()
    reviewed_times: list[datetime] = []
    reviewed_candidate_path = None
    for key, review_kind, decision in (
        ("provenance_review", "provenance-security", "approved-for-publication"),
        ("plasma_review", "plasma-scientific-continuation", "approved"),
    ):
        review = require_exact_keys(
            loaded[key],
            {
                "schema_version",
                "record_type",
                "checkpoint",
                "execution_epoch",
                "review_kind",
                "decision",
                "reviewed_candidate",
                "published_f118",
                "reviewer",
                "reviewed_utc",
                "findings",
                "limitations",
                "verified",
            },
            f"schema-2 F118 {key}",
        )
        candidate = require_exact_keys(
            review["reviewed_candidate"],
            {"path", "sha256"},
            f"schema-2 F118 {key} candidate",
        )
        reviewer = require_exact_keys(
            review["reviewer"],
            {"agent_id", "identity"},
            f"schema-2 F118 {key} reviewer",
        )
        candidate_path = require_nonempty_string(
            candidate["path"], f"schema-2 F118 {key} candidate path"
        )
        selected_candidate = Path(candidate_path)
        if (
            not selected_candidate.is_absolute()
            or selected_candidate != absolute_path(selected_candidate)
        ):
            raise ValueError(f"schema-2 F118 {key} candidate path is not absolute normalized")
        if reviewed_candidate_path is None:
            reviewed_candidate_path = candidate_path
        if (
            review["schema_version"] != 1
            or review["record_type"]
            != "stage-i-current-source-authority-supersession-independent-review"
            or review["checkpoint"] != "F-118"
            or review["execution_epoch"] != EXECUTION_EPOCH
            or review["review_kind"] != review_kind
            or review["decision"] != decision
            or candidate_path == str(root / F118_RELATIVE)
            or candidate_path != reviewed_candidate_path
            or candidate["sha256"] != digests["evidence"]
            or review["published_f118"]
            != {"path": str(root / F118_RELATIVE), "sha256": digests["evidence"]}
            or review["verified"] != verified
        ):
            raise ValueError(f"schema-2 F118 {key} identity or decision differs")
        agent = require_nonempty_string(
            reviewer["agent_id"], f"schema-2 F118 {key} reviewer agent"
        )
        require_nonempty_string(
            reviewer["identity"], f"schema-2 F118 {key} reviewer identity"
        )
        if agent in reviewers:
            raise ValueError("schema-2 F118 independent reviews reuse one reviewer")
        reviewers.add(agent)
        for field in ("findings", "limitations"):
            if (
                not isinstance(review[field], list)
                or not review[field]
                or any(not isinstance(item, str) or not item for item in review[field])
            ):
                raise ValueError(f"schema-2 F118 {key} {field} differ")
        if INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION not in review["limitations"]:
            raise ValueError(
                f"schema-2 F118 {key} omits the reviewer identity limitation"
            )
        reviewed = parse_utc_timestamp(review["reviewed_utc"], f"schema-2 F118 {key}")
        if reviewed < generated:
            raise ValueError(f"schema-2 F118 {key} predates evidence")
        reviewed_times.append(reviewed)

    audit = require_exact_keys(
        loaded["publication_audit"],
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "execution_epoch",
            "published_utc",
            "artifact",
            "independent_reviews",
            "historical_f116_authority",
            "source_archive_catalog",
            "authority_and_enforcement",
            "publication",
        },
        "schema-2 F118 publication audit",
    )
    historical = require_exact_keys(
        evidence["predecessor_authorities"],
        {"historical_f116"},
        "schema-2 F118 predecessor authorities",
    )["historical_f116"]
    historical_bindings = require_exact_keys(
        historical,
        {"evidence", "publication_audit", "provenance_review", "plasma_review"},
        "schema-2 F118 historical F116 bindings",
    )
    historical_digests = {
        f"{key}_sha256": require_v2_input_binding(
            historical_bindings[key], f"schema-2 F118 historical F116 {key}"
        )["sha256"]
        for key in historical_bindings
    }
    published = parse_utc_timestamp(audit["published_utc"], "schema-2 F118 publication")
    if (
        audit["schema_version"] != 1
        or audit["record_type"]
        != "stage-i-current-source-authority-supersession-publication-audit"
        or audit["checkpoint"] != "F-118"
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or audit["historical_f116_authority"] != historical_digests
        or audit["authority_and_enforcement"] != F118_AUTHORIZATION
        or audit["publication"]
        != "recoverable-forward-transaction-with-publication-audit-commit-marker-under-stage-i-lock"
        or published < generated
        or any(published < reviewed for reviewed in reviewed_times)
    ):
        raise ValueError("schema-2 F118 publication audit identity or authority differs")
    exact_publication_binding(
        audit["artifact"],
        root / F118_RELATIVE,
        digests["evidence"],
        "schema-2 F118 publication artifact",
        mode="0444",
    )
    reviews = require_exact_keys(
        audit["independent_reviews"],
        {
            "reviews_bind_exact_published_f118_sha256",
            "provenance_security",
            "plasma_scientific_continuation",
        },
        "schema-2 F118 publication reviews",
    )
    if reviews["reviews_bind_exact_published_f118_sha256"] != digests["evidence"]:
        raise ValueError("schema-2 F118 publication review digest binding differs")
    exact_publication_binding(
        reviews["provenance_security"],
        root / F118_PROVENANCE_REVIEW_RELATIVE,
        digests["provenance_review"],
        "schema-2 F118 provenance review publication",
        mode="0444",
    )
    exact_publication_binding(
        reviews["plasma_scientific_continuation"],
        root / F118_PLASMA_REVIEW_RELATIVE,
        digests["plasma_review"],
        "schema-2 F118 plasma review publication",
        mode="0444",
    )
    audit_catalog = require_exact_keys(
        audit["source_archive_catalog"],
        {
            "readme",
            "sha256sums",
            "bridge_bundle",
            "predecessor_current_source_bundle",
            "current_source_bundle",
            "corrupt_c7_absent_from_active_checksum_ledger",
            "sole_current_source_bundle",
        },
        "schema-2 F118 publication source-archive catalog",
    )
    for key, relative, digest_key in (
        ("readme", Path("source-archives/README.md"), "readme_sha256"),
        ("sha256sums", Path("source-archives/SHA256SUMS"), "sha256sums_sha256"),
    ):
        exact_publication_binding(
            audit_catalog[key],
            root / relative,
            str(after[digest_key]),
            f"schema-2 F118 publication catalog {key}",
            mode="0644",
        )
        read_confined_file_sha256(
            root,
            relative.as_posix(),
            str(after[digest_key]),
            f"schema-2 F118 publication catalog {key}",
            expected_mode=0o644,
        )
    expected_bundles = (
        (
            "bridge_bundle",
            bridge,
            {"role": "retained-non-current-bridge", "selected_as_current": False},
        ),
        (
            "predecessor_current_source_bundle",
            predecessor,
            {"role": "retained-non-current-predecessor", "selected_as_current": False},
        ),
        ("current_source_bundle", current, {"selected_as_current": True}),
    )
    for key, bundle, extra in expected_bundles:
        expected = {
            "path": str(root / safe_relative_path(str(bundle["path"]), f"F118 {key} path")),
            "sha256": bundle["sha256"],
            "mode": "0644",
            "links": 1,
            "head": bundle["head"],
            **extra,
        }
        if audit_catalog[key] != expected:
            raise ValueError(f"schema-2 F118 publication catalog {key} differs")
        read_confined_file_sha256(
            root,
            str(bundle["path"]),
            str(bundle["sha256"]),
            f"schema-2 F118 publication catalog {key}",
            expected_mode=0o644,
        )
    if (
        audit_catalog["corrupt_c7_absent_from_active_checksum_ledger"] is not True
        or audit_catalog["sole_current_source_bundle"] != str(root / Path(str(current["path"])))
    ):
        raise ValueError("schema-2 F118 publication catalog authority differs")


def require_v2_repository_binding(value: object, label: str) -> dict[str, str]:
    """Require one exact repository path, revision, and digest binding."""

    binding = require_exact_keys(value, {"path", "revision", "sha256"}, label)
    path = safe_relative_path(
        require_nonempty_string(binding["path"], f"{label} path"),
        f"{label} path",
    ).as_posix()
    revision = require_nonempty_string(binding["revision"], f"{label} revision")
    if REVISION_PATTERN.fullmatch(revision) is None:
        raise ValueError(f"{label} revision is invalid")
    return {
        "path": path,
        "revision": revision,
        "sha256": require_sha256(binding["sha256"], f"{label} SHA-256"),
    }


def v2_request(root: Path, args: argparse.Namespace) -> tuple[Path, dict[str, object]]:
    """Authenticate and parse the exact reviewed V2 recost request."""

    request_path, payload = read_confined_file_sha256(
        root,
        args.recost_request_relative_path,
        args.expected_request_sha256,
        "V2 recost request",
        expected_mode=0o644,
    )
    try:
        request = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("V2 recost request is not valid JSON") from error
    return request_path, require_exact_keys(
        request,
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "artifact_name",
            "execution_epoch",
            "generated_utc",
            "expires_utc",
            "scope",
            "barrier",
            "inputs",
            "authorization",
        },
        "V2 recost request",
    )


def declared_process_independence_assurance(
    role_agents: dict[str, str], label: str
) -> dict[str, object]:
    """Represent strict declared process independence without identity overclaim."""

    if len(role_agents) < 2:
        raise ValueError(f"{label} must declare at least two distinct process roles")
    retained = {
        require_nonempty_string(role, f"{label} role"): require_nonempty_string(
            agent, f"{label} agent"
        )
        for role, agent in role_agents.items()
    }
    if len(retained) != len(role_agents) or len(set(retained.values())) != len(retained):
        raise ValueError(f"{label} roles or agents are not strictly distinct")
    return {
        "basis": "declared-process-independence",
        "strict_distinct_role_and_agent_declarations": True,
        "declared_role_agents": retained,
        "cryptographic_identity_verified": False,
        "non_cryptographic_limitation": INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION,
    }


def validate_schema2_independent_review(payload: bytes, root: Path,
                                        artifact_name: str,
                                        args: argparse.Namespace, *,
                                        artifact_generated_utc: object,
                                        publication_published_utc: object | None = None
                                        ) -> dict[str, object]:
    """Validate review bytes and return explicit declared-process assurance."""

    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("schema-2 recost evidence independent review is not valid JSON") from error
    review = require_exact_keys(
        value,
        {
            "schema_version",
            "record_type",
            "execution_epoch",
            "reviewed_utc",
            "decision",
            "reviewer",
            "candidate",
            "scope",
        },
        "schema-2 recost evidence independent review",
    )
    reviewer = require_exact_keys(
        review["reviewer"],
        {"agent_id", "independent_from_generator"},
        "schema-2 recost evidence reviewer",
    )
    reviewer_agent = require_nonempty_string(
        reviewer["agent_id"], "schema-2 recost evidence reviewer"
    )
    if (
        review["schema_version"] != 1
        or review["record_type"] != "stage-i-recost-recommendation-independent-review"
        or review["execution_epoch"] != EXECUTION_EPOCH
        or review["decision"] != "approved-for-publication"
        or reviewer["independent_from_generator"] is not True
        or review["candidate"]
        != {
            "path": str(root / "accounting" / artifact_name),
            "sha256": args.expected_artifact_sha256,
        }
        or review["scope"] != {"non_authorizing": True}
    ):
        raise ValueError("schema-2 recost evidence independent review differs")
    generated = parse_utc_timestamp(
        artifact_generated_utc, "schema-2 recost evidence artifact generation timestamp"
    )
    reviewed = parse_utc_timestamp(
        review["reviewed_utc"], "schema-2 recost evidence review timestamp"
    )
    if reviewed < generated:
        raise ValueError(
            "schema-2 recost evidence independent review predates artifact generation"
        )
    if publication_published_utc is None:
        if reviewed > datetime.now(timezone.utc):
            raise ValueError(
                "schema-2 recost evidence independent review is in the future"
            )
    else:
        published = parse_utc_timestamp(
            publication_published_utc,
            "schema-2 recost evidence publication timestamp",
        )
        if reviewed > published:
            raise ValueError(
                "schema-2 recost evidence publication predates independent review"
            )
    generator_declaration = getattr(args, "expected_generator_revision", None)
    if not isinstance(generator_declaration, str) or not generator_declaration:
        generator_declaration = getattr(args, "expected_generator_sha256", None)
    if not isinstance(generator_declaration, str) or not generator_declaration:
        generator_declaration = "declared-recost-generator-process"
    return declared_process_independence_assurance(
        {
            "recost-generator": f"generator:{generator_declaration}",
            "independent-reviewer": reviewer_agent,
        },
        "schema-2 recost independent-review process",
    )


def schema2_independent_review(root: Path, artifact_name: str,
                               args: argparse.Namespace, *,
                               artifact_generated_utc: object,
                               publication_published_utc: object | None = None,
                               retained_review: tuple[Path, bytes] | None = None,
                               ) -> dict[str, object]:
    """Authenticate the independent review required before evidence publication."""

    if retained_review is None:
        review_path, payload = read_confined_file_sha256(
            root,
            args.independent_review_relative_path,
            args.expected_independent_review_sha256,
            "schema-2 recost evidence independent review",
            expected_mode=0o444,
        )
    else:
        review_path, payload = retained_review
        expected = (
            root
            / safe_relative_path(
                args.independent_review_relative_path,
                "schema-2 recost evidence independent-review relative path",
            )
        )
        if review_path != expected:
            raise ValueError("schema-2 recost evidence independent review path differs")
        if sha256_bytes(payload) != args.expected_independent_review_sha256:
            raise ValueError("schema-2 recost evidence independent review checksum has changed")
    validate_schema2_independent_review(
        payload,
        root,
        artifact_name,
        args,
        artifact_generated_utc=artifact_generated_utc,
        publication_published_utc=publication_published_utc,
    )
    return {
        "path": str(review_path),
        "sha256": args.expected_independent_review_sha256,
        "mode": "0444",
        "links": 1,
    }


def authenticate_f117_failed_attempt_for_f119(root: Path) -> dict[str, object]:
    """Authenticate the exact unpromoted F117 quartet superseded only by F119."""

    allowed_names = {relative.name for relative in F117_FAILED_ATTEMPT_RELATIVES.values()}
    with absolute_descriptor(
        root / "accounting",
        "schema-2 F119 accounting directory",
        flags=os.O_RDONLY | os.O_DIRECTORY,
    ) as descriptor:
        names = sorted(entry.name for entry in os.scandir(descriptor))
        unexpected = sorted(
            name
            for name in names
            if name.startswith(f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_F117_")
            and name not in allowed_names
        )
        prior_schema2 = sorted(
            name
            for name in names
            if name.endswith("_recost_evidence.json.publication_audit.json")
            and name
            not in {
                f"{LEGACY_F114_RECOST_NAME}.publication_audit.json",
                f"{F119_ARTIFACT_NAME}.publication_audit.json",
            }
        )
    if unexpected:
        raise ValueError(
            "schema-2 F119 supersession requires no F117 artifact, review, audit, "
            f"or extra prerequisite: {unexpected}"
        )
    if prior_schema2:
        raise ValueError(
            "schema-2 F119 retained F114 predecessor was already superseded by "
            f"another publication: {prior_schema2}"
        )

    records: dict[str, dict[str, object]] = {}
    for key, relative in F117_FAILED_ATTEMPT_RELATIVES.items():
        _, payload = read_confined_file_sha256(
            root,
            relative.as_posix(),
            F117_FAILED_ATTEMPT_SHA256[key],
            f"schema-2 F119 retained failed F117 {key}",
            expected_mode=0o644,
        )
        records[key] = canonical_json_object(
            payload, f"schema-2 F119 retained failed F117 {key}"
        )

    packet = require_exact_keys(
        records["packet"],
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "artifact_name",
            "execution_epoch",
            "generated_utc",
            "expires_utc",
            "requested_by",
            "scope",
            "barrier",
            "inputs",
            "recommendations",
            "draft_policy",
        },
        "schema-2 F119 retained failed F117 packet",
    )
    request = require_exact_keys(
        records["request"],
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "artifact_name",
            "execution_epoch",
            "generated_utc",
            "expires_utc",
            "requested_by",
            "scope",
            "barrier",
            "inputs",
            "recommendations",
        },
        "schema-2 F119 retained failed F117 request",
    )
    failed_artifact = "mks24_stage_i_E03_forcing_policy_F117_recost_evidence.json"
    if (
        packet["schema_version"] != 1
        or packet["record_type"] != "stage-i-recost-request-draft-packet"
        or request["schema_version"] != 2
        or request["record_type"] != "stage-i-recost-recommendation-request"
        or packet["checkpoint"] != "F-117"
        or request["checkpoint"] != "F-117"
        or packet["artifact_name"] != failed_artifact
        or request["artifact_name"] != failed_artifact
        or packet["execution_epoch"] != EXECUTION_EPOCH
        or request["execution_epoch"] != EXECUTION_EPOCH
    ):
        raise ValueError("schema-2 F119 retained failed F117 identity differs")
    policy = require_exact_keys(
        packet["draft_policy"],
        {
            "independent_review_created",
            "self_approved",
            "scheduler_mutation_authorized",
            "canonical_mutation_authorized",
            "required_storage_safety_bytes",
        },
        "schema-2 F119 retained failed F117 draft policy",
    )
    if (
        policy["independent_review_created"] is not False
        or policy["self_approved"] is not False
        or policy["scheduler_mutation_authorized"] is not False
        or policy["canonical_mutation_authorized"] is not False
    ):
        raise ValueError("schema-2 F119 retained failed F117 packet over-authorizes")
    require_integer(
        policy["required_storage_safety_bytes"],
        "schema-2 F119 retained failed F117 safety bytes",
        minimum=1,
    )
    for key in (
        "checkpoint",
        "artifact_name",
        "execution_epoch",
        "generated_utc",
        "expires_utc",
        "requested_by",
        "scope",
        "barrier",
        "recommendations",
    ):
        if packet[key] != request[key]:
            raise ValueError(
                f"schema-2 F119 retained failed F117 packet/request {key} differs"
            )
    parse_utc_timestamp(packet["generated_utc"], "schema-2 F119 failed F117 generation")
    parse_utc_timestamp(packet["expires_utc"], "schema-2 F119 failed F117 expiry")
    require_nonempty_string(packet["requested_by"], "schema-2 F119 failed F117 author")
    require_nonempty_string(packet["scope"], "schema-2 F119 failed F117 scope")
    require_exact_keys(
        packet["barrier"], {"recorded_segments"}, "schema-2 F119 failed F117 barrier"
    )
    require_exact_keys(
        packet["recommendations"],
        {"mode", "max_wave_nodes", "profiles"},
        "schema-2 F119 failed F117 recommendations",
    )
    packet_inputs = packet["inputs"]
    request_inputs = request["inputs"]
    expected_input_keys = {
        "reconciliation",
        "ledger",
        "reservations",
        "manifests",
        "scheduler_evidence",
        "storage_evidence",
        "source_bundle",
        "matrix",
        "stage_i_helper",
        "ceiling_evidence",
        "ceiling_publication_audit",
        "source_authority",
        "qualification_approval",
        "predecessor_recost",
        "predecessor_recost_independent_review",
        "predecessor_recost_publication_audit",
        "r17_readiness_evidence",
        "r17_readiness_independent_review",
        "r17_readiness_publication_audit",
    }
    packet_inputs = require_exact_keys(
        packet_inputs, expected_input_keys, "schema-2 F119 failed F117 packet inputs"
    )
    request_inputs = require_exact_keys(
        request_inputs, expected_input_keys, "schema-2 F119 failed F117 request inputs"
    )
    live_keys = {
        "reconciliation",
        "ledger",
        "reservations",
        "manifests",
        "scheduler_evidence",
        "storage_evidence",
    }
    if any(packet_inputs[key] is not None for key in live_keys):
        raise ValueError("schema-2 F119 failed F117 packet retains live bindings")
    if any(
        packet_inputs[key] != request_inputs[key]
        for key in expected_input_keys - live_keys
    ):
        raise ValueError("schema-2 F119 failed F117 authority bindings differ")
    expected_legacy = {
        "path": f"accounting/{LEGACY_F114_RECOST_NAME}",
        "sha256": LEGACY_F114_RECOST_SHA256,
    }
    expected_legacy_audit = {
        "path": f"accounting/{LEGACY_F114_RECOST_NAME}.publication_audit.json",
        "sha256": LEGACY_F114_PUBLICATION_AUDIT_SHA256,
    }
    if (
        request_inputs["predecessor_recost"] != expected_legacy
        or request_inputs["predecessor_recost_independent_review"] is not None
        or request_inputs["predecessor_recost_publication_audit"]
        != expected_legacy_audit
        or not isinstance(request_inputs["source_authority"], dict)
        or request_inputs["source_authority"].get("checkpoint") != "F-116"
    ):
        raise ValueError("schema-2 F119 failed F117 predecessor or authority differs")
    for input_key, record_key in (
        ("reconciliation", "reconciliation"),
        ("storage_evidence", "storage"),
    ):
        if request_inputs[input_key] != {
            "path": F117_FAILED_ATTEMPT_RELATIVES[record_key].as_posix(),
            "sha256": F117_FAILED_ATTEMPT_SHA256[record_key],
        }:
            raise ValueError(f"schema-2 F119 failed F117 {input_key} binding differs")

    reconciliation = records["reconciliation"]
    counts = reconciliation.get("counts")
    if (
        reconciliation.get("execution_epoch") != EXECUTION_EPOCH
        or reconciliation.get("root") != str(root)
        or reconciliation.get("consistent") is not True
        or reconciliation.get("issues") != []
        or not isinstance(counts, dict)
        or any(
            isinstance(counts.get(key), bool)
            or not isinstance(counts.get(key), int)
            or counts.get(key) < 0
            for key in COUNT_KEYS
        )
        or counts.get("transactions") != 0
        or counts.get("active_reservations") != 0
    ):
        raise ValueError("schema-2 F119 retained failed F117 reconciliation is not clean")
    storage = require_exact_keys(
        records["storage"],
        {
            "schema_version",
            "record_type",
            "execution_epoch",
            "root",
            "measured_utc",
            "available_bytes",
            "retained_stage_i_bytes",
            "required_safety_bytes",
            "projected_authorized_wave_growth_bytes",
            "projection_method",
            "profile_projections_sha256",
        },
        "schema-2 F119 retained failed F117 storage",
    )
    if (
        storage["schema_version"] != 1
        or storage["record_type"] != "stage-i-storage-evidence"
        or storage["execution_epoch"] != EXECUTION_EPOCH
        or storage["root"] != str(root)
        or storage["measured_utc"] != request["generated_utc"]
        or storage["projection_method"] != STORAGE_PROJECTION_METHOD
    ):
        raise ValueError("schema-2 F119 retained failed F117 storage identity differs")
    require_integer(storage["available_bytes"], "schema-2 F119 failed F117 available bytes", minimum=1)
    require_integer(storage["retained_stage_i_bytes"], "schema-2 F119 failed F117 retained bytes")
    require_integer(storage["required_safety_bytes"], "schema-2 F119 failed F117 safety bytes", minimum=1)
    require_integer(
        storage["projected_authorized_wave_growth_bytes"],
        "schema-2 F119 failed F117 projected growth",
        minimum=1,
    )
    require_sha256(
        storage["profile_projections_sha256"],
        "schema-2 F119 failed F117 profile projection SHA-256",
    )
    return {
        "checkpoint": "F-117",
        **{
            key: {
                "path": relative.as_posix(),
                "sha256": F117_FAILED_ATTEMPT_SHA256[key],
            }
            for key, relative in F117_FAILED_ATTEMPT_RELATIVES.items()
        },
        "status": "authenticated-unpromoted-failed-attempt",
    }


def validate_f119_failed_f117_predecessor(
    value: object,
    request_inputs: dict[str, object],
    provenance: dict[str, object],
    root: Path,
    artifact_name: str,
    checkpoint: str,
    request_generated_utc: object,
) -> None:
    """Require the exact retained-F114/failed-F117 predecessor only for F119."""

    predecessor_path = root / "accounting" / LEGACY_F114_RECOST_NAME
    audit_path = predecessor_path.with_name(f"{predecessor_path.name}.publication_audit.json")
    expected_request_predecessor = {
        "path": predecessor_path.relative_to(root).as_posix(),
        "sha256": LEGACY_F114_RECOST_SHA256,
    }
    expected_request_audit = {
        "path": audit_path.relative_to(root).as_posix(),
        "sha256": LEGACY_F114_PUBLICATION_AUDIT_SHA256,
    }
    has_failed_f117_marker = (
        isinstance(value, dict)
        and (
            value.get("bootstrap") == F119_LEGACY_BOOTSTRAP
            or "superseded_failed_attempt" in value
        )
    )
    has_legacy_f114_identity = (
        isinstance(value, dict)
        and (
            value.get("path") == str(predecessor_path)
            or value.get("artifact_name") == LEGACY_F114_RECOST_NAME
            or value.get("sha256") == LEGACY_F114_RECOST_SHA256
            or value.get("publication_audit_path") == str(audit_path)
            or value.get("publication_audit_sha256")
            == LEGACY_F114_PUBLICATION_AUDIT_SHA256
            or value.get("checkpoint") == "F-114"
        )
    ) or (
        isinstance(request_inputs, dict)
        and (
            request_inputs.get("predecessor_recost") == expected_request_predecessor
            or request_inputs.get("predecessor_recost_publication_audit")
            == expected_request_audit
        )
    ) or (
        isinstance(provenance, dict)
        and (
            provenance.get("predecessor_recost_sha256") == LEGACY_F114_RECOST_SHA256
            or provenance.get("predecessor_recost_publication_audit_sha256")
            == LEGACY_F114_PUBLICATION_AUDIT_SHA256
        )
    )
    if checkpoint != "F-119" or artifact_name != F119_ARTIFACT_NAME:
        if has_failed_f117_marker or has_legacy_f114_identity:
            raise ValueError("legacy F114 predecessor consumption is restricted to exact F119")
        return
    request_inputs = require_exact_keys(
        request_inputs,
        {
            "reconciliation",
            "ledger",
            "reservations",
            "manifests",
            "scheduler_evidence",
            "storage_evidence",
            "source_bundle",
            "matrix",
            "stage_i_helper",
            "ceiling_evidence",
            "ceiling_publication_audit",
            "source_authority",
            "qualification_approval",
            "predecessor_recost",
            "predecessor_recost_independent_review",
            "predecessor_recost_publication_audit",
            "r17_readiness_evidence",
            "r17_readiness_independent_review",
            "r17_readiness_publication_audit",
        },
        "schema-2 F119 request inputs",
    )
    predecessor = require_exact_keys(
        value,
        {
            "path",
            "artifact_name",
            "sha256",
            "publication_audit_path",
            "publication_audit_sha256",
            "independent_review_path",
            "independent_review_sha256",
            "checkpoint",
            "generated_utc",
            "published_utc",
            "bootstrap",
            "superseded_failed_attempt",
        },
        "schema-2 F119 predecessor recost",
    )
    if (
        predecessor["path"] != str(predecessor_path)
        or predecessor["artifact_name"] != LEGACY_F114_RECOST_NAME
        or predecessor["sha256"] != LEGACY_F114_RECOST_SHA256
        or predecessor["publication_audit_path"] != str(audit_path)
        or predecessor["publication_audit_sha256"]
        != LEGACY_F114_PUBLICATION_AUDIT_SHA256
        or predecessor["independent_review_path"] is not None
        or predecessor["independent_review_sha256"] is not None
        or predecessor["checkpoint"] != "F-114"
        or predecessor["bootstrap"] != F119_LEGACY_BOOTSTRAP
        or request_inputs.get("predecessor_recost") != expected_request_predecessor
        or request_inputs.get("predecessor_recost_independent_review") is not None
        or request_inputs.get("predecessor_recost_publication_audit")
        != expected_request_audit
        or provenance.get("predecessor_recost_sha256") != LEGACY_F114_RECOST_SHA256
        or provenance.get("predecessor_recost_independent_review_sha256") is not None
        or provenance.get("predecessor_recost_publication_audit_sha256")
        != LEGACY_F114_PUBLICATION_AUDIT_SHA256
    ):
        raise ValueError("schema-2 F119 predecessor identity or request binding differs")
    _, predecessor_payload = read_confined_file_sha256(
        root,
        predecessor_path.relative_to(root).as_posix(),
        LEGACY_F114_RECOST_SHA256,
        "schema-2 F119 retained legacy F114 predecessor",
        expected_mode=0o644,
    )
    _, audit_payload = read_confined_file_sha256(
        root,
        audit_path.relative_to(root).as_posix(),
        LEGACY_F114_PUBLICATION_AUDIT_SHA256,
        "schema-2 F119 retained legacy F114 predecessor audit",
        expected_mode=0o644,
    )
    try:
        retained_predecessor = json.loads(predecessor_payload)
        retained_audit = json.loads(audit_payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("schema-2 F119 retained legacy F114 predecessor is not JSON") from error
    if (
        not isinstance(retained_predecessor, dict)
        or retained_predecessor.get("schema_version") != 1
        or retained_predecessor.get("record_type")
        != "stage-i-clean-partial-recost-checkpoint"
        or retained_predecessor.get("checkpoint") != "F-114"
        or retained_predecessor.get("execution_epoch") != EXECUTION_EPOCH
    ):
        raise ValueError("schema-2 F119 retained legacy F114 predecessor identity differs")
    generated = parse_utc_timestamp(
        retained_predecessor.get("generated_utc"),
        "schema-2 F119 retained legacy F114 generation",
    )
    request_generated = parse_utc_timestamp(
        request_generated_utc, "schema-2 F119 request generation"
    )
    if (
        not isinstance(retained_audit, dict)
        or retained_audit.get("schema_version") != 1
        or retained_audit.get("record_type") != "observed-publication"
        or retained_audit.get("execution_epoch") != EXECUTION_EPOCH
        or retained_audit.get("artifact")
        != {
            "path": str(predecessor_path),
            "sha256": LEGACY_F114_RECOST_SHA256,
            "mode": "0644",
            "links": 1,
        }
    ):
        raise ValueError("schema-2 F119 retained legacy F114 publication audit differs")
    published = parse_utc_timestamp(
        retained_audit.get("published_utc"),
        "schema-2 F119 retained legacy F114 publication",
    )
    if (
        predecessor["generated_utc"] != generated.isoformat()
        or predecessor["published_utc"] != published.isoformat()
        or published < generated
        or request_generated < published
    ):
        raise ValueError("schema-2 F119 retained legacy F114 chronology differs")
    failed_f117 = authenticate_f117_failed_attempt_for_f119(root)
    if predecessor["superseded_failed_attempt"] != failed_f117:
        raise ValueError("schema-2 F119 failed-F117 supersession binding differs")


def schema2_publication_context(value: dict[str, object], path: Path,
                                args: argparse.Namespace,
                                scheduler: list[dict[str, object]], *,
                                allow_expired: bool = False,
                                retained_review: tuple[Path, bytes] | None = None,
                                ) -> dict[str, object]:
    """Validate and return the strict non-authorizing schema-2 publication context."""

    root = path.parent.parent
    request_path, payload = read_confined_file_sha256(
        root,
        args.recost_request_relative_path,
        args.expected_request_sha256,
        "schema-2 recost request",
        expected_mode=0o644,
    )
    try:
        request_value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("schema-2 recost request is not valid JSON") from error
    request = require_exact_keys(
        request_value,
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "artifact_name",
            "execution_epoch",
            "generated_utc",
            "expires_utc",
            "requested_by",
            "scope",
            "barrier",
            "inputs",
            "recommendations",
        },
        "schema-2 recost request",
    )
    if (
        request["schema_version"] != 2
        or request["record_type"] != "stage-i-recost-recommendation-request"
        or request["execution_epoch"] != EXECUTION_EPOCH
    ):
        raise ValueError("schema-2 recost request identity differs")
    artifact_name = require_artifact_name(
        require_nonempty_string(request["artifact_name"], "schema-2 request artifact name")
    )
    checkpoint = require_safe_id(request["checkpoint"], "schema-2 request checkpoint")
    match = V2_RECOST_ARTIFACT_PATTERN.fullmatch(artifact_name)
    if (
        match is None
        or checkpoint != f"F-{match.group(1)}"
        or artifact_name != args.artifact_name
        or value.get("artifact_name") != artifact_name
        or value.get("checkpoint") != checkpoint
    ):
        raise ValueError("schema-2 artifact name and checkpoint binding differs")
    generated = parse_utc_timestamp(request["generated_utc"], "schema-2 request generation")
    expires = parse_utc_timestamp(request["expires_utc"], "schema-2 request expiry")
    now = datetime.now(timezone.utc)
    if generated > now + timedelta(minutes=5):
        raise ValueError("schema-2 recost request generation is in the future")
    if expires <= generated or expires - generated > V2_AUTHORIZATION_MAX_LIFETIME:
        raise ValueError("schema-2 recost request lifetime is invalid")
    if now > expires and not allow_expired:
        raise ValueError("schema-2 recost request has expired")
    for key in ("generated_utc", "expires_utc", "requested_by", "scope"):
        if value.get(key) != request[key]:
            raise ValueError(f"schema-2 artifact {key} binding differs")
    recommendations = validate_recost_recommendations(
        args.recost_recommendations_json, root, args
    )
    if value.get("recommendations") != recommendations:
        raise ValueError("schema-2 artifact recommendation packet differs")
    requested = require_exact_keys(
        request["recommendations"],
        {"mode", "max_wave_nodes", "profiles"},
        "schema-2 request recommendations",
    )
    if (
        requested["mode"] != recommendations["mode"]
        or requested["max_wave_nodes"]
        != recommendations["bounded_concurrency"]["max_wave_nodes"]
        or requested["profiles"] != recommendations["recommended_next_profiles"]
    ):
        raise ValueError("schema-2 request and artifact recommendations differ")
    authority = require_exact_keys(
        value.get("authority"),
        {
            "authorizing",
            "action_authority",
            "scheduler_mutation_authorized",
            "canonical_mutation_authorized",
        },
        "schema-2 recost authority",
    )
    if authority != {
        "authorizing": False,
        "action_authority": "none-until-independent-review-and-publication",
        "scheduler_mutation_authorized": False,
        "canonical_mutation_authorized": False,
    }:
        raise ValueError("schema-2 recost evidence improperly grants authority")
    publication_requirements = require_exact_keys(
        value.get("publication_requirements"),
        {
            "independent_review_required",
            "publication_audit_required",
            "published_mode",
            "published_links",
            "controller_consumption_requires_exact_published_sha256",
        },
        "schema-2 publication requirements",
    )
    if publication_requirements != {
        "independent_review_required": True,
        "publication_audit_required": True,
        "published_mode": "0444",
        "published_links": 1,
        "controller_consumption_requires_exact_published_sha256": True,
    }:
        raise ValueError("schema-2 publication requirements differ")
    barrier = require_exact_keys(
        value.get("barrier"),
        {"job_ids", "recorded_segments", "scheduler_evidence"},
        "schema-2 artifact barrier",
    )
    request_barrier = require_exact_keys(
        request["barrier"], {"recorded_segments"}, "schema-2 request barrier"
    )
    validate_bounded_barrier(barrier["recorded_segments"], scheduler)
    validate_fresh_r12_rerun_evidence(
        recommendations["recommended_next_profiles"],
        barrier["recorded_segments"],
    )
    if (
        barrier["recorded_segments"] != request_barrier["recorded_segments"]
        or barrier["scheduler_evidence"] != scheduler
        or barrier["job_ids"] != sorted(str(item["job_id"]) for item in scheduler)
    ):
        raise ValueError("schema-2 artifact barrier bindings differ")
    provenance = require_exact_keys(
        value.get("provenance"),
        SCHEMA2_RECOST_PROVENANCE_KEYS
        | ({"scheduler_sha256"} if recommendations["mode"] == "sole-next-profile" else set()),
        "schema-2 artifact provenance",
    )
    request_inputs = request["inputs"]
    if not isinstance(request_inputs, dict):
        raise ValueError("schema-2 request inputs must be an object")
    source_authority = validate_f118_source_authority_binding(
        request_inputs.get("source_authority"), root, args
    )
    if provenance["source_authority"] != source_authority:
        raise ValueError("schema-2 artifact F118 source authority binding differs")
    validate_f119_failed_f117_predecessor(
        value.get("predecessor_recost"),
        request_inputs,
        provenance,
        root,
        artifact_name,
        checkpoint,
        request["generated_utc"],
    )
    for key, expected in (
        ("request_sha256", args.expected_request_sha256),
        ("generator_sha256", args.expected_generator_sha256),
        ("generator_revision", args.expected_generator_revision),
        ("stage_i_helper_sha256", args.expected_stage_i_sha256),
        ("stage_i_helper_revision", args.expected_stage_i_revision),
        ("source_bundle_sha256", args.expected_source_bundle_sha256),
        ("source_bundle_verified_revisions", args.expected_source_bundle_verified_revisions_json),
        ("scheduler_evidence", scheduler),
    ):
        if provenance.get(key) != expected:
            raise ValueError(f"schema-2 artifact provenance {key} binding differs")
    require_selected_counts(
        value.get("reconcile", {}).get("counts") if isinstance(value.get("reconcile"), dict) else None,
        expected_counts(args),
        "schema-2 artifact reconciliation counts",
    )
    review = schema2_independent_review(
        root,
        artifact_name,
        args,
        artifact_generated_utc=value.get("generated_utc"),
        retained_review=retained_review,
    )
    return {
        "schema_version": 2,
        "artifact_name": artifact_name,
        "checkpoint": checkpoint,
        "generated_utc": request["generated_utc"],
        "expires_utc": request["expires_utc"],
        "request": {"path": str(request_path), "sha256": args.expected_request_sha256},
        "artifact": {"basename": artifact_name, "sha256": args.expected_artifact_sha256},
        "independent_review": review,
        "authority": authority,
        "publication_requirements": publication_requirements,
        "recommendations": recommendations,
        "barrier": barrier,
        "provenance": provenance,
        "predecessor_recost": value.get("predecessor_recost"),
        "reconciliation": value.get("reconcile"),
        "budget": value.get("budget"),
        "storage": value.get("storage"),
        "r17_readiness": value.get("r17_readiness"),
        "controller_enforcement": {
            "state": recommendations["controller_consumption_state"],
            "launch_authority": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        },
    }


def v2_publication_context(value: dict[str, object], path: Path,
                           args: argparse.Namespace,
                           scheduler: list[dict[str, object]], *,
                           allow_expired: bool = False,
                           retained_review: tuple[Path, bytes] | None = None,
                           ) -> dict[str, object]:
    """Validate and return every generalized V2 publication binding."""

    if schema2_evidence_mode(args):
        return schema2_publication_context(
            value,
            path,
            args,
            scheduler,
            allow_expired=allow_expired,
            retained_review=retained_review,
        )
    root = path.parent.parent
    request_path, request = v2_request(root, args)
    if request["schema_version"] != 1:
        raise ValueError("V2 recost request schema version differs")
    if request["record_type"] != "stage-i-barrier-recost-request":
        raise ValueError("V2 recost request record type differs")
    if request["execution_epoch"] != EXECUTION_EPOCH:
        raise ValueError("V2 recost request execution epoch differs")
    artifact_name = require_artifact_name(
        require_nonempty_string(request["artifact_name"], "V2 request artifact name")
    )
    if artifact_name != args.artifact_name or value.get("artifact_name") != artifact_name:
        raise ValueError("V2 artifact basename binding differs")
    checkpoint = require_safe_id(request["checkpoint"], "V2 request checkpoint")
    match = V2_RECOST_ARTIFACT_PATTERN.fullmatch(artifact_name)
    if match is None or checkpoint != f"F-{match.group(1)}":
        raise ValueError("V2 artifact basename and checkpoint ID differ")
    if value.get("checkpoint") != checkpoint:
        raise ValueError("V2 artifact checkpoint binding differs")
    generated = parse_utc_timestamp(request["generated_utc"], "V2 request generation")
    expires = parse_utc_timestamp(request["expires_utc"], "V2 request expiry")
    now = datetime.now(timezone.utc)
    if generated > now + timedelta(minutes=5):
        raise ValueError("V2 recost request generation is in the future")
    if expires <= generated or expires - generated > V2_AUTHORIZATION_MAX_LIFETIME:
        raise ValueError("V2 recost request lifetime is invalid")
    if now > expires and not allow_expired:
        raise ValueError("V2 recost request has expired")
    if value.get("generated_utc") != request["generated_utc"]:
        raise ValueError("V2 artifact generation timestamp differs")
    if value.get("expires_utc") != request["expires_utc"]:
        raise ValueError("V2 artifact expiry timestamp differs")
    require_nonempty_string(request["scope"], "V2 request scope")
    if value.get("scope") != request["scope"]:
        raise ValueError("V2 artifact scope differs")

    request_authorization = require_exact_keys(
        request["authorization"],
        {"mode", "max_wave_nodes", "profiles"},
        "V2 request authorization",
    )
    authorization = validate_bounded_wave_authorization(
        args.authorized_bounded_wave_json, root, args
    )
    if request_authorization["mode"] != authorization["mode"]:
        raise ValueError("V2 request and artifact authorization modes differ")
    if request_authorization["profiles"] != authorization["authorized_next_profiles"]:
        raise ValueError("V2 request and artifact profile vectors differ")
    if (
        request_authorization["max_wave_nodes"]
        != authorization["bounded_concurrency"]["max_wave_nodes"]
    ):
        raise ValueError("V2 request and artifact wave-node bindings differ")

    inputs = require_exact_keys(
        request["inputs"],
        {
            "reconciliation",
            "ledger",
            "reservations",
            "manifests",
            "scheduler_evidence",
            "storage_evidence",
            "source_bundle",
            "matrix",
            "stage_i_helper",
            "ceiling_evidence",
            "ceiling_publication_audit",
            "predecessor_recost",
            "predecessor_recost_publication_audit",
            "r17_readiness_evidence",
        },
        "V2 request inputs",
    )
    simple_bindings = {
        key: require_v2_input_binding(inputs[key], f"V2 {key}")
        for key in (
            "reconciliation",
            "ledger",
            "reservations",
            "storage_evidence",
            "ceiling_evidence",
            "ceiling_publication_audit",
            "predecessor_recost",
            "predecessor_recost_publication_audit",
        )
    }
    readiness_binding = None
    if inputs["r17_readiness_evidence"] is not None:
        readiness_binding = require_v2_input_binding(
            inputs["r17_readiness_evidence"], "V2 R17 readiness"
        )
    manifests = inputs["manifests"]
    if not isinstance(manifests, list):
        raise ValueError("V2 manifest binding vector must be a list")
    manifest_bindings = [
        require_v2_input_binding(item, f"V2 manifest binding {index}")
        for index, item in enumerate(manifests)
    ]
    retained_input_payloads: dict[str, bytes] = {}
    for key, binding in simple_bindings.items():
        _, retained_input_payloads[key] = read_confined_file_sha256(
            root,
            binding["path"],
            binding["sha256"],
            f"V2 {key}",
            expected_mode=0o644,
        )
    for index, binding in enumerate(manifest_bindings):
        read_confined_file_sha256(
            root,
            binding["path"],
            binding["sha256"],
            f"V2 manifest binding {index}",
            expected_mode=0o644,
        )
    scheduler_inputs = inputs["scheduler_evidence"]
    if not isinstance(scheduler_inputs, list) or not scheduler_inputs:
        raise ValueError("V2 scheduler input vector must be nonempty")
    scheduler_bindings = [
        require_v2_input_binding(item, f"V2 scheduler binding {index}")
        for index, item in enumerate(scheduler_inputs)
    ]
    scheduler_projection = [
        {"path": str(item["path"]), "sha256": str(item["sha256"])}
        for item in scheduler
    ]
    if scheduler_bindings != scheduler_projection:
        raise ValueError("V2 request and artifact scheduler vectors differ")

    source_bundle = require_exact_keys(
        inputs["source_bundle"],
        {"path", "sha256", "verified_revisions"},
        "V2 source bundle binding",
    )
    source_binding = {
        "path": safe_relative_path(
            require_nonempty_string(source_bundle["path"], "V2 source bundle path"),
            "V2 source bundle path",
        ).as_posix(),
        "sha256": require_sha256(source_bundle["sha256"], "V2 source bundle SHA-256"),
        "verified_revisions": source_bundle["verified_revisions"],
    }
    if source_binding != {
        "path": args.source_bundle_relative_path,
        "sha256": args.expected_source_bundle_sha256,
        "verified_revisions": args.expected_source_bundle_verified_revisions_json,
    }:
        raise ValueError("V2 request source-bundle binding differs")
    matrix_binding = require_v2_repository_binding(inputs["matrix"], "V2 matrix")
    repository = repository_root(initial_source_path())
    matrix_relative = safe_relative_path(matrix_binding["path"], "V2 matrix path")
    matrix_path = repository / matrix_relative
    require_file_sha256(
        matrix_path,
        matrix_binding["sha256"],
        "V2 matrix",
        expected_links=1,
    )
    require_committed_revision_bytes(
        repository,
        matrix_relative,
        matrix_binding["revision"],
        matrix_binding["sha256"],
        "V2 matrix",
    )
    helper_binding = require_exact_keys(
        inputs["stage_i_helper"], {"revision", "sha256"}, "V2 Stage I helper"
    )
    helper = {
        "revision": require_nonempty_string(
            helper_binding["revision"], "V2 Stage I helper revision"
        ),
        "sha256": require_sha256(
            helper_binding["sha256"], "V2 Stage I helper SHA-256"
        ),
    }
    if helper != {
        "revision": args.expected_stage_i_revision,
        "sha256": args.expected_stage_i_sha256,
    }:
        raise ValueError("V2 request Stage I helper binding differs")

    provenance_keys = set(V2_PROVENANCE_KEYS)
    if authorization["mode"] == "sole-next-profile":
        provenance_keys.add("scheduler_sha256")
    provenance = require_exact_keys(
        value.get("provenance"), provenance_keys, "V2 artifact provenance"
    )
    provenance_bindings = {
        "request_sha256": args.expected_request_sha256,
        "generator_sha256": args.expected_generator_sha256,
        "generator_revision": args.expected_generator_revision,
        "stage_i_helper_sha256": args.expected_stage_i_sha256,
        "stage_i_helper_revision": args.expected_stage_i_revision,
        "matrix_sha256": matrix_binding["sha256"],
        "matrix_revision": matrix_binding["revision"],
        "source_bundle_sha256": args.expected_source_bundle_sha256,
        "source_bundle_verified_revisions": (
            args.expected_source_bundle_verified_revisions_json
        ),
        "ceiling_evidence_sha256": simple_bindings["ceiling_evidence"]["sha256"],
        "ceiling_publication_audit_sha256": (
            simple_bindings["ceiling_publication_audit"]["sha256"]
        ),
        "storage_evidence_sha256": simple_bindings["storage_evidence"]["sha256"],
        "reconciliation_sha256": simple_bindings["reconciliation"]["sha256"],
        "ledger_sha256": simple_bindings["ledger"]["sha256"],
        "reservations_sha256": simple_bindings["reservations"]["sha256"],
        "scheduler_evidence": scheduler,
        "predecessor_recost_sha256": simple_bindings["predecessor_recost"]["sha256"],
        "predecessor_recost_publication_audit_sha256": (
            simple_bindings["predecessor_recost_publication_audit"]["sha256"]
        ),
    }
    for key, expected in provenance_bindings.items():
        if provenance.get(key) != expected:
            raise ValueError(f"V2 artifact provenance {key} binding differs")
    if authorization["mode"] == "sole-next-profile" and (
        len(scheduler) != 1
        or provenance["scheduler_sha256"] != scheduler[0]["sha256"]
    ):
        raise ValueError("sole-profile V2 scheduler compatibility binding differs")

    request_barrier = require_exact_keys(
        request["barrier"], {"recorded_segments"}, "V2 request barrier"
    )
    barrier = require_exact_keys(
        value.get("barrier"),
        {"job_ids", "recorded_segments", "scheduler_evidence"},
        "V2 artifact barrier",
    )
    if barrier["recorded_segments"] != request_barrier["recorded_segments"]:
        raise ValueError("V2 request and artifact barrier vectors differ")
    validate_bounded_barrier(barrier["recorded_segments"], scheduler)
    validate_fresh_r12_rerun_evidence(
        authorization["authorized_next_profiles"],
        barrier["recorded_segments"],
    )
    for index, item in enumerate(scheduler):
        completed = parse_utc_timestamp(
            item["completed_utc"],
            f"bounded scheduler evidence {index} completion time",
        )
        if completed > generated + SCHEDULER_TIME_TOLERANCE:
            raise ValueError(
                f"bounded scheduler evidence {index} completes after the request"
            )
    barrier_jobs = sorted(str(item["job_id"]) for item in scheduler)
    if barrier["job_ids"] != barrier_jobs or barrier["scheduler_evidence"] != scheduler:
        raise ValueError("V2 artifact barrier ledger-tail bindings differ")

    promoted_f113 = require_exact_keys(
        value.get("promoted_f113"),
        {
            "path",
            "sha256",
            "publication_audit_path",
            "publication_audit_sha256",
            "publication_audit",
        },
        "V2 promoted F113",
    )
    if (
        promoted_f113["path"] != str(root / F113_RELATIVE)
        or simple_bindings["ceiling_evidence"]["path"] != F113_RELATIVE.as_posix()
        or promoted_f113["sha256"] != simple_bindings["ceiling_evidence"]["sha256"]
        or promoted_f113["publication_audit_path"]
        != str(root / F113_PUBLICATION_AUDIT_RELATIVE)
        or simple_bindings["ceiling_publication_audit"]["path"]
        != F113_PUBLICATION_AUDIT_RELATIVE.as_posix()
        or promoted_f113["publication_audit_sha256"]
        != simple_bindings["ceiling_publication_audit"]["sha256"]
    ):
        raise ValueError("V2 promoted F113 binding differs")
    try:
        retained_f113_audit = json.loads(
            retained_input_payloads["ceiling_publication_audit"]
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("V2 promoted F113 publication audit is not valid JSON") from error
    if promoted_f113["publication_audit"] != retained_f113_audit:
        raise ValueError("V2 promoted F113 publication audit object differs")

    artifact_ledger = require_exact_keys(
        value.get("ledger"),
        {"rows", "sha256", "cumulative_stage_i_node_hours"},
        "V2 artifact ledger",
    )
    if (
        artifact_ledger["rows"] != args.expected_ledger_rows
        or artifact_ledger["sha256"] != simple_bindings["ledger"]["sha256"]
    ):
        raise ValueError("V2 artifact ledger binding differs")
    artifact_reservations = require_exact_keys(
        value.get("reservations"),
        {"rows", "sha256", "active"},
        "V2 artifact reservations",
    )
    if (
        artifact_reservations["rows"] != args.expected_reservations
        or artifact_reservations["active"] != args.expected_active_reservations
        or artifact_reservations["sha256"] != simple_bindings["reservations"]["sha256"]
    ):
        raise ValueError("V2 artifact reservation binding differs")
    artifact_manifests = require_exact_keys(
        value.get("manifests"),
        {"rows", "bindings", "authenticated_lineages_sha256"},
        "V2 artifact manifests",
    )
    if artifact_manifests["rows"] != len(manifest_bindings):
        raise ValueError("V2 artifact manifest count differs")
    if artifact_manifests["bindings"] != manifest_bindings:
        raise ValueError("V2 artifact manifest vector differs")
    if provenance.get("authenticated_lineages_sha256") != artifact_manifests[
        "authenticated_lineages_sha256"
    ]:
        raise ValueError("V2 authenticated-lineage digest differs")
    require_sha256(
        artifact_manifests["authenticated_lineages_sha256"],
        "V2 authenticated-lineage digest",
    )

    predecessor = require_exact_keys(
        value.get("predecessor_recost"),
        {
            "artifact_name",
            "path",
            "sha256",
            "publication_audit_path",
            "publication_audit_sha256",
            "checkpoint",
            "generated_utc",
            "published_utc",
        },
        "V2 artifact predecessor recost",
    )
    if (
        predecessor["path"] != str(root / simple_bindings["predecessor_recost"]["path"])
        or predecessor["sha256"] != simple_bindings["predecessor_recost"]["sha256"]
        or predecessor["publication_audit_path"]
        != str(root / simple_bindings["predecessor_recost_publication_audit"]["path"])
        or predecessor["publication_audit_sha256"]
        != simple_bindings["predecessor_recost_publication_audit"]["sha256"]
    ):
        raise ValueError("V2 predecessor recost binding differs")
    predecessor_artifact_name = require_artifact_name(
        require_nonempty_string(
            predecessor["artifact_name"], "V2 predecessor artifact name"
        )
    )
    if (
        Path(str(predecessor["path"])).name != predecessor_artifact_name
        or Path(str(predecessor["publication_audit_path"])).name
        != f"{predecessor_artifact_name}.publication_audit.json"
    ):
        raise ValueError("V2 predecessor artifact-name binding differs")
    predecessor_checkpoint = require_safe_id(
        predecessor["checkpoint"], "V2 predecessor checkpoint"
    )
    predecessor_artifact_match = V2_RECOST_ARTIFACT_PATTERN.fullmatch(
        predecessor_artifact_name
    )
    predecessor_match = re.fullmatch(r"F-([0-9]+)", predecessor_checkpoint)
    checkpoint_match = re.fullmatch(r"F-([0-9]+)", checkpoint)
    if (
        predecessor_match is None
        or checkpoint_match is None
        or int(predecessor_match.group(1)) >= int(checkpoint_match.group(1))
    ):
        raise ValueError("V2 predecessor checkpoint does not precede publication")
    if (
        predecessor_artifact_match is not None
        and predecessor_checkpoint != f"F-{predecessor_artifact_match.group(1)}"
    ):
        raise ValueError("V2 predecessor artifact name and checkpoint ID differ")
    parse_utc_timestamp(predecessor["generated_utc"], "V2 predecessor generation")
    parse_utc_timestamp(predecessor["published_utc"], "V2 predecessor publication")

    reconcile = value.get("reconcile")
    if not isinstance(reconcile, dict):
        raise ValueError("V2 artifact reconciliation must be an object")
    try:
        retained_reconcile = json.loads(retained_input_payloads["reconciliation"])
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("V2 retained reconciliation is not valid JSON") from error
    if reconcile != retained_reconcile:
        raise ValueError("V2 artifact reconciliation differs from retained evidence")
    require_selected_counts(
        reconcile.get("counts"), expected_counts(args), "V2 artifact reconciliation counts"
    )
    budget = value.get("budget")
    if not isinstance(budget, dict):
        raise ValueError("V2 artifact budget must be an object")
    projection_sha256 = sha256_bytes(
        (json.dumps(budget, sort_keys=True) + "\n").encode()
    )
    if provenance.get("computed_projection_sha256") != projection_sha256:
        raise ValueError("V2 artifact budget projection digest differs")
    storage = require_exact_keys(
        value.get("storage"),
        {
            "available_bytes",
            "retained_stage_i_bytes",
            "required_safety_bytes",
            "projected_authorized_wave_growth_bytes",
            "headroom_after_authorized_wave_and_safety_bytes",
        },
        "V2 artifact storage",
    )
    available = require_integer(storage["available_bytes"], "V2 available storage")
    require_integer(storage["retained_stage_i_bytes"], "V2 retained Stage I storage")
    safety = require_integer(
        storage["required_safety_bytes"], "V2 required storage safety", minimum=1
    )
    projected_growth = require_integer(
        storage["projected_authorized_wave_growth_bytes"],
        "V2 projected authorized-wave growth",
    )
    headroom = require_integer(
        storage["headroom_after_authorized_wave_and_safety_bytes"],
        "V2 storage headroom",
    )
    if headroom != available - safety - projected_growth:
        raise ValueError("V2 artifact storage arithmetic differs")
    profiles = authorization["authorized_next_profiles"]
    assert isinstance(profiles, list)
    if projected_growth != sum(
        require_integer(
            profile["estimated_storage_bytes"],
            f"V2 profile {index} estimated storage bytes",
            minimum=1,
        )
        for index, profile in enumerate(profiles)
        if isinstance(profile, dict)
    ):
        raise ValueError("V2 artifact storage projection differs from authorization")
    contains_r17 = any(
        isinstance(profile, dict) and profile.get("case_id") == "R17"
        for profile in profiles
    )
    readiness = value.get("r17_readiness")
    readiness_sha = provenance.get("r17_readiness_evidence_sha256")
    if contains_r17:
        if (
            not isinstance(readiness, dict)
            or readiness_binding is None
            or readiness_sha != readiness_binding["sha256"]
        ):
            raise ValueError("V2 R17 authorization lacks reviewed readiness")
        read_confined_file_sha256(
            root,
            readiness_binding["path"],
            readiness_binding["sha256"],
            "V2 R17 readiness",
            expected_mode=0o644,
        )
    elif readiness is not None or readiness_sha is not None or readiness_binding is not None:
        raise ValueError("non-R17 V2 artifact retains R17 readiness")

    context = {
        "schema_version": 2,
        "artifact_name": artifact_name,
        "checkpoint": checkpoint,
        "generated_utc": request["generated_utc"],
        "expires_utc": request["expires_utc"],
        "request": {
            "path": str(request_path),
            "sha256": args.expected_request_sha256,
        },
        "artifact": {
            "basename": artifact_name,
            "sha256": args.expected_artifact_sha256,
        },
        "authorization": authorization,
        "barrier": barrier,
        "ledger_tail_job_ids": barrier_jobs,
        "inputs": inputs,
        "provenance": provenance,
        "predecessor_recost": predecessor,
        "promoted_f113": promoted_f113,
        "ledger": artifact_ledger,
        "reservations": artifact_reservations,
        "manifests": artifact_manifests,
        "reconciliation": reconcile,
        "budget": budget,
        "storage": storage,
        "r17_readiness": readiness,
        "controller_enforcement": {
            "state": authorization["controller_consumption_state"],
            "bounded_wave_authorizing": (
                authorization["authorizing"]
                if authorization["mode"] == "bounded-wave"
                else None
            ),
        },
    }
    return context


def selected_v2_publication_context(paths: dict[str, Path],
                                    args: argparse.Namespace, *,
                                    allow_expired: bool = False
                                    ) -> dict[str, object]:
    """Return the full V2 context from the selected staged or canonical name."""

    if publication_mode(args) != "v2":
        raise ValueError("generalized publication context is available only in V2 mode")
    artifact = paths["canonical"] if entry_exists(paths["canonical"]) else paths["staged"]
    retained = require_file_sha256(
        artifact,
        args.expected_artifact_sha256,
        "V2 recost artifact",
        expected_mode=args.expected_artifact_mode,
    )
    value = validate_artifact_payload(
        retained, artifact, args, allow_expired_v2=allow_expired
    )
    scheduler = validate_bounded_scheduler_evidence(
        args.bounded_wave_scheduler_evidence_json
    )
    return v2_publication_context(
        value, artifact, args, scheduler, allow_expired=allow_expired
    )


def expected_counts(args: argparse.Namespace) -> dict[str, int]:
    """Return the explicitly authorized reconcile counts."""

    counts = {
        key: int(getattr(args, f"expected_{key}"))
        for key in COUNT_KEYS
    }
    if counts["transactions"] != 0:
        raise ValueError("expected Stage I transaction count must be zero")
    if counts["active_reservations"] != 0:
        raise ValueError("expected active reservation count must be zero")
    return counts


def require_selected_counts(value: object, expected: dict[str, int],
                            label: str) -> None:
    """Require the selected artifact or reconcile report counts."""

    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    selected = {key: value.get(key) for key in COUNT_KEYS}
    if selected != expected:
        raise ValueError(f"{label} differ: expected {expected}, found {selected}")


def validate_artifact_payload(retained: bytes, path: Path,
                              args: argparse.Namespace, *,
                              allow_expired_v2: bool = False,
                              retained_review: tuple[Path, bytes] | None = None,
                              ) -> dict[str, object]:
    """Validate authenticated recost-artifact bytes."""

    try:
        value = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"recost artifact is not valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"recost artifact must be a JSON object: {path}")
    mode = publication_mode(args)
    if mode == "v2":
        if schema2_evidence_mode(args):
            require_exact_keys(
                value, SCHEMA2_RECOST_ARTIFACT_KEYS, "schema-2 recost evidence"
            )
            if value["schema_version"] != 2:
                raise ValueError("schema-2 recost evidence schema version differs")
            if value["record_type"] != "stage-i-recost-recommendation-evidence":
                raise ValueError("schema-2 recost evidence record type differs")
        else:
            require_exact_keys(value, V2_ARTIFACT_KEYS, "V2 recost artifact")
            if value["schema_version"] != 1:
                raise ValueError("V2 recost artifact schema version differs")
            if value["record_type"] != "stage-i-barrier-recost-checkpoint":
                raise ValueError("V2 recost artifact record type differs")
    bindings = [
        ("execution epoch", args.artifact_epoch_pointer, EXECUTION_EPOCH),
        (
            "Stage I helper SHA-256",
            args.artifact_stage_i_sha256_pointer,
            args.expected_stage_i_sha256,
        ),
        (
            "generator SHA-256",
            args.artifact_generator_sha256_pointer,
            args.expected_generator_sha256,
        ),
    ]
    if mode == "sole-profile":
        bindings[1:1] = [
            (
                "sole next-segment profile",
                args.artifact_profile_pointer,
                args.authorized_next_segment_profile_json,
            )
        ]
        bindings.append(
            (
                "scheduler SHA-256",
                args.artifact_scheduler_sha256_pointer,
                args.expected_scheduler_sha256,
            )
        )
    else:
        root = path.parent.parent
        packet = (
            validate_recost_recommendations(args.recost_recommendations_json, root, args)
            if schema2_evidence_mode(args)
            else validate_bounded_wave_authorization(
                args.authorized_bounded_wave_json, root, args
            )
        )
        scheduler = validate_bounded_scheduler_evidence(
            args.bounded_wave_scheduler_evidence_json
        )
        validate_bounded_barrier(
            json_pointer(value, "/barrier/recorded_segments", "bounded barrier"),
            scheduler,
        )
        bindings.extend(
            [
                (
                    "V2 authorization",
                    args.artifact_authorization_pointer,
                    packet,
                ),
                (
                    "V2 scheduler evidence",
                    args.artifact_scheduler_evidence_pointer,
                    scheduler,
                ),
                (
                    "V2 barrier scheduler evidence",
                    args.artifact_barrier_scheduler_evidence_pointer,
                    scheduler,
                ),
                (
                    "source bundle SHA-256",
                    args.artifact_source_bundle_sha256_pointer,
                    args.expected_source_bundle_sha256,
                ),
                (
                    "Stage I helper revision",
                    args.artifact_stage_i_revision_pointer,
                    args.expected_stage_i_revision,
                ),
                (
                    "generator revision",
                    args.artifact_generator_revision_pointer,
                    args.expected_generator_revision,
                ),
                (
                    "source bundle verified revisions",
                    args.artifact_source_bundle_revisions_pointer,
                    args.expected_source_bundle_verified_revisions_json,
                ),
            ]
        )
    for label, pointer, expected in bindings:
        selected = json_pointer(value, pointer, label)
        if selected != expected:
            raise ValueError(
                f"recost artifact {label} differs: expected {expected!r}, "
                f"found {selected!r}"
            )
    require_selected_counts(
        json_pointer(value, args.artifact_counts_pointer, "counts"),
        expected_counts(args),
        "recost artifact counts",
    )
    if mode == "v2":
        count_bindings = (
            ("reservation rows", "/reservations/rows", args.expected_reservations),
            (
                "active reservations",
                "/reservations/active",
                args.expected_active_reservations,
            ),
            ("ledger rows", "/ledger/rows", args.expected_ledger_rows),
            ("manifest rows", "/manifests/rows", args.expected_manifests),
        )
        for label, pointer, expected in count_bindings:
            selected = json_pointer(value, pointer, label)
            if selected != expected:
                raise ValueError(
                    f"recost artifact {label} differs: expected {expected!r}, "
                    f"found {selected!r}"
                )
        v2_publication_context(
            value,
            path,
            args,
            scheduler,
            allow_expired=allow_expired_v2,
            retained_review=retained_review,
        )
    return value


def read_artifact(path: Path, args: argparse.Namespace, *,
                  expected_links: int = 1,
                  allow_expired_v2: bool = False) -> dict[str, object]:
    """Authenticate and validate one staged or promoted recost artifact."""

    retained = require_file_sha256(
        path,
        args.expected_artifact_sha256,
        "recost artifact",
        expected_mode=args.expected_artifact_mode,
        expected_links=expected_links,
    )
    return validate_artifact_payload(
        retained, path, args, allow_expired_v2=allow_expired_v2
    )


def require_system_executable_profile(value: os.stat_result, label: str) -> None:
    """Require one root-owned, single-link, non-writable executable."""

    mode = stat.S_IMODE(value.st_mode)
    if (
        not stat.S_ISREG(value.st_mode)
        or value.st_uid != 0
        or value.st_nlink != 1
        or mode & 0o022
        or mode & 0o111 == 0
    ):
        raise ValueError(f"{label} does not have the trusted system executable profile")


def require_system_directory_profile(value: os.stat_result, label: str) -> None:
    """Require one root-owned directory without unsafe replacement permissions."""

    if (
        not stat.S_ISDIR(value.st_mode)
        or value.st_uid != 0
        or stat.S_IMODE(value.st_mode) & 0o022
    ):
        raise ValueError(f"{label} does not have the trusted system directory profile")


@contextmanager
def authenticated_python_interpreter():
    """Yield the exact root-owned interpreter inode running this process."""

    interpreter = Path(sys.executable).resolve(strict=True)
    with absolute_descriptor(
        interpreter, "authenticated Python interpreter", flags=os.O_RDONLY
    ) as descriptor:
        profile = os.fstat(descriptor)
        require_system_executable_profile(profile, "authenticated Python interpreter")
        if profile_identity(profile) != profile_identity(os.stat("/proc/self/exe")):
            raise ValueError("selected Python interpreter is not this process")
        os.set_inheritable(descriptor, True)
        yield interpreter, descriptor
        if profile_identity(os.fstat(descriptor)) != profile_identity(profile):
            raise ValueError("authenticated Python interpreter descriptor changed")


@contextmanager
def authenticated_stage_i_python_interpreter():
    """Yield the controller's fixed root-owned Python interpreter descriptor."""

    with absolute_descriptor(
        STAGE_I_PYTHON, "Stage I Python interpreter", flags=os.O_RDONLY
    ) as descriptor:
        profile = os.fstat(descriptor)
        require_system_executable_profile(profile, "Stage I Python interpreter")
        if profile_identity(profile) != profile_identity(STAGE_I_PYTHON.lstat()):
            raise ValueError("Stage I Python interpreter pathname changed")
        os.set_inheritable(descriptor, True)
        yield STAGE_I_PYTHON, descriptor
        if (
            profile_identity(os.fstat(descriptor)) != profile_identity(profile)
            or profile_identity(STAGE_I_PYTHON.lstat()) != profile_identity(profile)
        ):
            raise ValueError("Stage I Python interpreter changed")


def hardened_python_child_environment() -> dict[str, str]:
    """Return a complete caller-independent environment for Python children."""

    return {
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": TRUSTED_SYSTEM_PATH,
        "PYTHONDONTWRITEBYTECODE": "1",
        "XDG_CONFIG_HOME": "/nonexistent",
    }


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


def hardened_git_arguments(root: Path, arguments: list[str]) -> list[str]:
    """Pin repository mode/worktree and disable local-config execution surfaces."""

    require_trusted_directory(root, "Git repository")
    metadata = root / ".git"
    retained = list(arguments)
    initializing = bool(retained) and retained[0] == "init"
    if entry_exists(metadata):
        require_trusted_directory(metadata, "Git metadata directory")
        repository_location = [
            f"--git-dir={metadata}",
            f"--work-tree={root}",
        ]
    else:
        repository_location = [f"--git-dir={root}"]
    if initializing:
        repository_location = []
    if retained and retained[0] == "diff":
        retained[1:1] = ["--no-ext-diff", "--no-textconv"]
    return [
        *repository_location,
        "-c",
        "core.fsmonitor=false",
        "-c",
        f"core.hooksPath={os.devnull}",
        "-c",
        f"core.attributesFile={os.devnull}",
        "-c",
        f"core.excludesFile={os.devnull}",
        "-c",
        f"init.templateDir={os.devnull}",
        "-c",
        "init.defaultObjectFormat=sha1",
        *retained,
    ]


def git_run(root: Path, arguments: list[str], *,
            capture_output: bool = False,
            pass_fds: tuple[int, ...] = ()) -> subprocess.CompletedProcess:
    """Run one exact authenticated Git descriptor with isolated configuration."""

    with absolute_descriptor(
        GIT_EXEC_PATH, "authenticated Git execution path",
        flags=os.O_RDONLY | os.O_DIRECTORY,
    ) as execution_path:
        require_system_directory_profile(
            os.fstat(execution_path), "authenticated Git execution path"
        )
        with absolute_descriptor(
            GIT, "authenticated absolute Git binary", flags=os.O_RDONLY
        ) as git:
            require_system_executable_profile(
                os.fstat(git), "authenticated absolute Git binary"
            )
            try:
                return subprocess.run(
                    [
                        str(GIT),
                        "--no-replace-objects",
                        "-C",
                        str(root),
                        *hardened_git_arguments(root, arguments),
                    ],
                    executable=f"/proc/self/fd/{git}",
                    check=False,
                    stdin=subprocess.DEVNULL,
                    capture_output=capture_output,
                    env=hardened_git_environment(),
                    pass_fds=(git, *pass_fds),
                    timeout=120,
                )
            except (OSError, subprocess.TimeoutExpired) as error:
                raise ValueError("authenticated absolute Git query failed to execute") from error


def git_descriptor_run(root: Path, descriptor: int, arguments: list[str], *,
                       capture_output: bool = False) -> subprocess.CompletedProcess:
    """Run one Git query against an inherited immutable descriptor path."""

    return git_run(
        root,
        arguments,
        capture_output=capture_output,
        pass_fds=(descriptor,),
    )


def authenticate_committed_stage_i(root_dir: Path, expected: str) -> Path:
    """Authenticate the committed live Stage I helper bytes."""

    expected = require_sha256(expected, "Stage I helper SHA-256")
    path = root_dir / STAGE_I_RELATIVE
    require_file_sha256(
        path, expected, "Stage I helper", expected_mode=0o644, expected_links=1
    )
    relative = STAGE_I_RELATIVE.as_posix()
    if git_run(root_dir, ["ls-files", "--error-unmatch", "--", relative]).returncode:
        raise ValueError("Stage I helper is not tracked by Git")
    for arguments in (
        ["diff", "--quiet", "--", relative],
        ["diff", "--cached", "--quiet", "--", relative],
    ):
        if git_run(root_dir, arguments).returncode:
            raise ValueError("Stage I helper must be committed before recost publication")
    retained = git_run(
        root_dir, ["show", f"HEAD:{relative}"], capture_output=True
    )
    if retained.returncode or sha256_bytes(retained.stdout) != expected:
        raise ValueError("committed Stage I helper checksum differs from live bytes")
    return path


def authenticate_committed_utility(root_dir: Path, expected: str) -> Path:
    """Authenticate the tracked, clean companion before canonical use."""

    expected = require_sha256(expected, "utility SHA-256")
    path = root_dir / UTILITY_RELATIVE
    require_file_sha256(
        path, expected, "retained utility", expected_mode=0o755, expected_links=1
    )
    relative = UTILITY_RELATIVE.as_posix()
    if git_run(root_dir, ["ls-files", "--error-unmatch", "--", relative]).returncode:
        raise ValueError("retained utility is not tracked by Git")
    for arguments in (
        ["diff", "--quiet", "--", relative],
        ["diff", "--cached", "--quiet", "--", relative],
    ):
        if git_run(root_dir, arguments).returncode:
            raise ValueError("retained utility must be committed before canonical use")
    retained = git_run(
        root_dir, ["show", f"HEAD:{relative}"], capture_output=True
    )
    if retained.returncode or sha256_bytes(retained.stdout) != expected:
        raise ValueError("committed utility checksum differs from live bytes")
    return path


def require_committed_revision_bytes(root_dir: Path, relative: Path,
                                     revision: str, expected: str,
                                     label: str) -> None:
    """Require one exact revision to retain the selected file bytes."""

    retained = git_run(
        root_dir,
        ["show", f"{revision}:{relative.as_posix()}"],
        capture_output=True,
    )
    if retained.returncode or sha256_bytes(retained.stdout) != expected:
        raise ValueError(f"{label} requested revision bytes differ")


def authenticate_bounded_source_provenance(root_dir: Path,
                                           args: argparse.Namespace, *,
                                           require_generator_at_head: bool = True
                                           ) -> None:
    """Authenticate bounded-wave generator and helper revision provenance."""

    if publication_mode(args) != "v2":
        return
    generator = root_dir / RECOST_RELATIVE
    require_file_sha256(
        generator,
        args.expected_generator_sha256,
        "Stage I recost generator",
        expected_mode=0o755,
        expected_links=1,
    )
    relative = RECOST_RELATIVE.as_posix()
    if git_run(root_dir, ["ls-files", "--error-unmatch", "--", relative]).returncode:
        raise ValueError("Stage I recost generator is not tracked by Git")
    for arguments in (
        ["diff", "--quiet", "--", relative],
        ["diff", "--cached", "--quiet", "--", relative],
    ):
        if git_run(root_dir, arguments).returncode:
            raise ValueError(
                "Stage I recost generator must be committed before publication"
            )
    require_committed_revision_bytes(
        root_dir,
        RECOST_RELATIVE,
        args.expected_generator_revision,
        args.expected_generator_sha256,
        "Stage I recost generator",
    )
    require_committed_revision_bytes(
        root_dir,
        STAGE_I_RELATIVE,
        args.expected_stage_i_revision,
        args.expected_stage_i_sha256,
        "Stage I helper",
    )
    if require_generator_at_head:
        head = git_run(root_dir, ["rev-parse", "HEAD"], capture_output=True)
        if (
            head.returncode
            or head.stdout.decode("utf-8").strip() != args.expected_generator_revision
        ):
            raise ValueError(
                "Stage I recost generator revision differs from repository HEAD"
            )


def authenticate_publication_sources(root_dir: Path,
                                     args: argparse.Namespace, *,
                                     require_generator_at_head: bool = True) -> Path:
    """Authenticate every mode-selected committed source file and revision."""

    helper = authenticate_committed_stage_i(root_dir, args.expected_stage_i_sha256)
    authenticate_bounded_source_provenance(
        root_dir, args, require_generator_at_head=require_generator_at_head
    )
    return helper


@contextmanager
def authenticated_recost_module(root_dir: Path, args: argparse.Namespace):
    """Load the exact committed V2 recost generator from authenticated bytes."""

    path = root_dir / RECOST_RELATIVE
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    module_name = f"_cgl_lf_stage_i_recost_v2_{uuid.uuid4().hex}"
    try:
        require_regular_profile(
            os.fstat(descriptor),
            path,
            "Stage I recost generator",
            expected_mode=0o755,
            expected_links=1,
        )
        if sha256_descriptor(descriptor) != args.expected_generator_sha256:
            raise ValueError("Stage I recost generator descriptor checksum differs")
        loader = SourceFileLoader(module_name, f"/proc/self/fd/{descriptor}")
        spec = importlib.util.spec_from_loader(module_name, loader)
        if spec is None or spec.loader is None:
            raise ValueError("Stage I recost generator module loader is unavailable")
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        for name in (
            "build_payload",
            "run_authenticated_reconcile",
            "require_empty_transaction_stores",
            "require_live_storage_boundary",
            "require_directory_measurement_boundaries",
        ):
            if not callable(getattr(module, name, None)):
                raise ValueError(f"Stage I recost generator lacks required V2 API: {name}")
        yield module
        if sha256_descriptor(descriptor) != args.expected_generator_sha256:
            raise ValueError("Stage I recost generator changed during V2 revalidation")
    finally:
        sys.modules.pop(module_name, None)
        os.close(descriptor)


def reexecute_v2_recost(paths: dict[str, Path], root_dir: Path, root: Path,
                        args: argparse.Namespace, artifact: Path, *,
                        allow_recost_journal: bool = False,
                        stage_i_lock_held: bool = False) -> None:
    """Rebuild one V2 artifact from every live evidence binding and compare bytes."""

    if publication_mode(args) != "v2":
        return
    retained = require_file_sha256(
        artifact,
        args.expected_artifact_sha256,
        "V2 recost artifact",
        expected_mode=args.expected_artifact_mode,
    )
    request = root / safe_relative_path(
        args.recost_request_relative_path, "recost request relative path"
    )
    namespace = argparse.Namespace(
        request=request,
        expected_request_sha256=args.expected_request_sha256,
        output=paths["staged"],
    )
    with authenticated_recost_module(root_dir, args) as module:
        original_empty_check = module.require_empty_transaction_stores
        if allow_recost_journal:
            module.require_empty_transaction_stores = lambda selected_root: None
        try:
            build = module.build_payload(
                namespace,
                root,
                root_dir / RECOST_RELATIVE,
                root_dir,
                args.expected_generator_sha256,
                stage_i_lock_held=stage_i_lock_held,
            )
        finally:
            module.require_empty_transaction_stores = original_empty_check
        for name in (
            "payload",
            "tracker",
            "reconcile",
            "helper_path",
            "helper_sha256",
            "storage_available_bytes",
            "storage_retained_stage_i_bytes",
            "storage_required_safety_bytes",
            "projected_storage_bytes",
            "storage_measurements",
        ):
            if not hasattr(build, name):
                raise ValueError(f"Stage I recost generator build lacks V2 binding: {name}")
        if not callable(getattr(build.tracker, "reauthenticate_all", None)):
            raise ValueError("Stage I recost generator tracker lacks V2 reauthentication")
        if not isinstance(build.payload, bytes) or build.payload != retained:
            raise ValueError("V2 recost artifact differs from authenticated re-execution")
        build.tracker.reauthenticate_all()
        live_reconcile = module.run_authenticated_reconcile(
            build.helper_path,
            build.helper_sha256,
            root,
            stage_i_lock_held=stage_i_lock_held,
        )
        if live_reconcile != build.reconcile:
            raise ValueError("V2 live reconciliation changed during re-execution")
        if not allow_recost_journal:
            module.require_empty_transaction_stores(root)
        module.require_live_storage_boundary(
            root,
            build.storage_available_bytes,
            build.storage_retained_stage_i_bytes,
            build.storage_required_safety_bytes,
            build.projected_storage_bytes,
        )
        module.require_directory_measurement_boundaries(build.storage_measurements)
    require_file_sha256(
        artifact,
        args.expected_artifact_sha256,
        "V2 recost artifact",
        expected_mode=args.expected_artifact_mode,
    )


@contextmanager
def stage_i_descriptor(root_dir: Path, expected: str):
    """Yield an authenticated descriptor for the committed Stage I helper."""

    path = authenticate_committed_stage_i(root_dir, expected)
    with regular_descriptor(
        path, "Stage I helper", expected_mode=0o644, expected_links=1
    ) as descriptor:
        if sha256_descriptor(descriptor) != expected:
            raise ValueError("Stage I helper descriptor checksum has changed")
        os.set_inheritable(descriptor, True)
        yield descriptor
        if sha256_descriptor(descriptor) != expected:
            raise ValueError("Stage I helper descriptor checksum changed during reconcile")


def stage_i_reconcile_environment(root_dir: Path, descriptor: int,
                                  python_descriptor: int) -> dict[str, str]:
    """Return the controller's exact inherited descriptor-reexec environment."""

    environment = hardened_python_child_environment()
    environment.update(
        {
            STAGE_I_SELF_DESCRIPTOR_ENV: str(descriptor),
            STAGE_I_PYTHON_DESCRIPTOR_ENV: str(python_descriptor),
            STAGE_I_SOURCE_ENV: str(root_dir / STAGE_I_RELATIVE),
            STAGE_I_REPOSITORY_ROOT_ENV: str(root_dir),
        }
    )
    return environment


def run_reconcile(root_dir: Path, root: Path, offline: bool,
                  args: argparse.Namespace, *,
                  stage_i_lock_held: bool = False) -> dict[str, object]:
    """Execute reconcile from the authenticated Stage I helper descriptor."""

    with stage_i_descriptor(root_dir, args.expected_stage_i_sha256) as descriptor:
        if stage_i_lock_held:
            module_name = (
                f"_cgl_lf_checkpoint_locked_reconcile_{os.getpid()}_{descriptor}_{id(root)}"
            )
            loader = SourceFileLoader(module_name, f"/proc/self/fd/{descriptor}")
            spec = importlib.util.spec_from_loader(module_name, loader)
            if spec is None or spec.loader is None:
                raise ValueError("Stage I helper module loader is unavailable")
            module = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = module
            try:
                spec.loader.exec_module(module)
                reconcile_report = getattr(module, "reconcile_report", None)
                if not callable(reconcile_report):
                    raise ValueError("Stage I helper lacks read-only reconcile_report")
                report = reconcile_report(root)
            finally:
                sys.modules.pop(module_name, None)
            if not isinstance(report, dict) or report.get("consistent") is not True:
                raise ValueError(
                    f"authenticated Stage I reconcile is inconsistent: {report!r}"
                )
            if report.get("execution_epoch") != EXECUTION_EPOCH:
                raise ValueError("authenticated Stage I reconcile returned the wrong epoch")
            require_selected_counts(
                report.get("counts"), expected_counts(args), "Stage I reconcile counts"
            )
            return report
        with authenticated_stage_i_python_interpreter() as (
            interpreter,
            python_descriptor,
        ):
            command = [
                str(interpreter),
                "-I",
                "-S",
                "-B",
                f"/proc/self/fd/{descriptor}",
                "--root",
                str(root),
            ]
            if offline:
                command.append("--allow-local-root")
            command.append("reconcile")
            completed = subprocess.run(
                command,
                executable=f"/proc/self/fd/{python_descriptor}",
                check=False,
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
                env=stage_i_reconcile_environment(
                    root_dir, descriptor, python_descriptor
                ),
                pass_fds=(python_descriptor, descriptor),
                timeout=120,
            )
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise ValueError(f"authenticated Stage I reconcile failed: {detail}")
    try:
        report = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise ValueError("authenticated Stage I reconcile did not return JSON") from error
    if not isinstance(report, dict) or report.get("consistent") is not True:
        raise ValueError(f"authenticated Stage I reconcile is inconsistent: {report!r}")
    if report.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError("authenticated Stage I reconcile returned the wrong epoch")
    require_selected_counts(
        report.get("counts"), expected_counts(args), "Stage I reconcile counts"
    )
    return report


def scheduler_environment() -> dict[str, str]:
    """Return the complete caller-independent scheduler child environment."""

    return {
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": TRUSTED_SYSTEM_PATH,
        "XDG_CONFIG_HOME": "/nonexistent",
    }


def require_empty_queue(root: Path, offline: bool, queue_file: str | None) -> None:
    """Require no queued CGL workflow job, with fixture injection only offline."""

    if offline and queue_file is None:
        raise ValueError("offline fixture root requires --squeue-file")
    if queue_file is not None:
        if not offline:
            raise ValueError("--squeue-file is permitted only with a local fixture root")
        path = require_offline_path(Path(queue_file), "queue fixture")
        with absolute_descriptor(path, "queue fixture", flags=os.O_RDONLY) as descriptor:
            require_regular_profile(os.fstat(descriptor), path, "queue fixture")
            output = read_descriptor_bytes(descriptor).decode("utf-8")
    else:
        try:
            user = pwd.getpwuid(os.geteuid()).pw_name
        except KeyError as error:
            raise ValueError("effective UID has no account name; cannot check the queue") from error
        try:
            with absolute_descriptor(
                SQUEUE, "authenticated squeue binary", flags=os.O_RDONLY
            ) as descriptor:
                require_system_executable_profile(
                    os.fstat(descriptor), "authenticated squeue binary"
                )
                output = subprocess.run(
                    [str(SQUEUE), "-h", "-u", user, "-o", "%i|%j|%T"],
                    executable=f"/proc/self/fd/{descriptor}",
                    check=True,
                    stdin=subprocess.DEVNULL,
                    capture_output=True,
                    text=True,
                    env=scheduler_environment(),
                    pass_fds=(descriptor,),
                    timeout=120,
                ).stdout
        except (
            OSError,
            subprocess.CalledProcessError,
            subprocess.TimeoutExpired,
        ) as error:
            raise ValueError("squeue is unavailable; refusing recost publication") from error
    queued = []
    for line in output.splitlines():
        if not line.strip():
            continue
        if line.strip() != line:
            raise ValueError("squeue output row has surrounding whitespace: " + line)
        fields = line.split("|")
        if len(fields) != 3:
            raise ValueError("squeue output row is not exactly three fields: " + line)
        if any(not field or field.strip() != field for field in fields):
            raise ValueError("squeue output row has malformed fields: " + line)
        if fields[1].startswith(CGL_JOB_NAME_PREFIX):
            queued.append(line)
    if queued:
        raise ValueError("another CGL job is queued: " + "; ".join(queued))


def validate_fixture_options(args: argparse.Namespace, offline: bool) -> None:
    """Keep mutation simulation and queue injection out of canonical use."""

    fixture_values = (
        getattr(args, "squeue_file", None),
        getattr(args, "pre_link_squeue_file", None),
        getattr(args, "pre_unlink_squeue_file", None),
        getattr(args, "simulate_post_link_failure", False),
        getattr(args, "simulate_link_return_failure", False),
        getattr(args, "simulate_single_link_post_exchange_failure", False),
        getattr(args, "simulate_adoption_interruption_after_journal", False),
        getattr(
            args,
            "simulate_adoption_interruption_before_forensic_directory_fsync",
            False,
        ),
        getattr(args, "simulate_adoption_interruption_after_forensic", False),
        getattr(
            args,
            "simulate_adoption_interruption_before_audit_directory_fsync",
            False,
        ),
        getattr(args, "simulate_adoption_interruption_after_audit", False),
    )
    if any(value not in (None, False) for value in fixture_values) and not offline:
        raise ValueError("fixture injection options are permitted only for a local root")


def authenticate_root_member_sha256(root: Path, path: Path, expected: str,
                                    label: str) -> None:
    """Authenticate one absolute in-root product without symlink traversal."""

    require_root_member(root, str(path), label)
    with absolute_descriptor(path, label, flags=os.O_RDONLY) as descriptor:
        require_regular_profile(
            os.fstat(descriptor), path, label, expected_links=1
        )
        if sha256_descriptor(descriptor) != require_sha256(expected, f"{label} SHA-256"):
            raise ValueError(f"{label} checksum has changed: {path}")


def authenticate_source_bundle_coverage(repository: Path, bundle: Path,
                                        expected_sha256: str,
                                        revisions: list[str], *,
                                        expected_mode: int) -> None:
    """Require a descriptor-bound self-contained bundle covering exact revisions."""

    expected = require_sha256(expected_sha256, "source bundle SHA-256")
    with absolute_descriptor(bundle, "source bundle", flags=os.O_RDONLY) as descriptor:
        profile = os.fstat(descriptor)
        require_regular_profile(
            profile,
            bundle,
            "source bundle",
            expected_mode=expected_mode,
            expected_links=1,
        )
        if sha256_descriptor(descriptor) != expected:
            raise ValueError("source bundle checksum differs before verification")
        descriptor_path = f"/proc/self/fd/{descriptor}"
        header = bytearray()
        os.lseek(descriptor, 0, os.SEEK_SET)
        while b"\n\n" not in header:
            block = os.read(descriptor, 4096)
            if not block or len(header) + len(block) > 1024 * 1024:
                raise ValueError("source bundle header is invalid or unexpectedly large")
            header.extend(block)
        try:
            lines = bytes(header).split(b"\n\n", 1)[0].decode().splitlines()
        except UnicodeDecodeError as error:
            raise ValueError("source bundle header is not UTF-8") from error
        if not lines or lines[0] not in {"# v2 git bundle", "# v3 git bundle"}:
            raise ValueError("source bundle header version is invalid")
        if any(line.startswith("-") for line in lines[1:]):
            raise ValueError("source bundle must be self-contained without prerequisites")
        advertised = []
        for line in lines[1:]:
            if line.startswith("@"):
                continue
            match = re.fullmatch(r"([0-9a-f]{40}) (.+)", line)
            if match is None:
                raise ValueError("source bundle advertised reference is malformed")
            advertised.append(match.group(1))
        if not advertised:
            raise ValueError("source bundle advertises no retained revision")
        if git_descriptor_run(
            repository,
            descriptor,
            ["bundle", "verify", descriptor_path],
            capture_output=True,
        ).returncode:
            raise ValueError("source bundle verification failed")
        with tempfile.TemporaryDirectory(prefix="cgl-lf-checkpoint-bundle-") as directory:
            os.chmod(directory, 0o700)
            previous_umask = os.umask(0o077)
            try:
                isolated = Path(directory) / "repository.git"
                if git_run(
                    repository,
                    ["init", "--bare", str(isolated)],
                    capture_output=True,
                ).returncode:
                    raise ValueError(
                        "cannot initialize isolated source-bundle verification"
                    )
                if git_descriptor_run(
                    isolated,
                    descriptor,
                    ["bundle", "unbundle", descriptor_path],
                    capture_output=True,
                ).returncode:
                    raise ValueError("source bundle cannot be reconstructed in isolation")
                if git_run(
                    isolated,
                    ["fsck", "--full", "--strict", "--no-reflogs", *advertised],
                    capture_output=True,
                ).returncode:
                    raise ValueError("source bundle advertised history is incomplete")
                for revision in revisions:
                    if git_run(
                        isolated, ["cat-file", "-e", f"{revision}^{{commit}}"]
                    ).returncode:
                        raise ValueError(
                            f"source bundle does not contain requested revision: {revision}"
                        )
                    if not any(
                        git_run(
                            isolated, ["merge-base", "--is-ancestor", revision, head]
                        ).returncode == 0
                        for head in advertised
                    ):
                        raise ValueError(
                            f"source bundle does not cover requested revision: {revision}"
                        )
            finally:
                os.umask(previous_umask)
        with absolute_descriptor(bundle, "source bundle", flags=os.O_RDONLY) as named:
            named_profile = os.fstat(named)
            require_regular_profile(
                named_profile,
                bundle,
                "source bundle",
                expected_mode=expected_mode,
                expected_links=1,
            )
            if (
                (named_profile.st_dev, named_profile.st_ino)
                != (profile.st_dev, profile.st_ino)
                or sha256_descriptor(named) != expected
            ):
                raise ValueError("source bundle pathname changed during verification")


def authenticate_bounded_scheduler_files(root: Path, args: argparse.Namespace
                                         ) -> tuple[Path, ...]:
    """Authenticate and parse every exact bounded-wave scheduler item."""

    retained = []
    evidence = validate_bounded_scheduler_evidence(
        args.bounded_wave_scheduler_evidence_json
    )
    for index, item in enumerate(evidence):
        path, payload = read_confined_file_sha256(
            root,
            str(item["path"]),
            str(item["sha256"]),
            f"bounded scheduler evidence {index}",
            expected_mode=args.expected_scheduler_mode,
        )
        try:
            lines = payload.decode("utf-8").splitlines()
        except UnicodeDecodeError as error:
            raise ValueError("bounded scheduler evidence is not UTF-8") from error
        if len(lines) != 1:
            raise ValueError("bounded scheduler evidence must contain exactly one row")
        fields = lines[0].split("|")
        if len(fields) != 8 or any(
            not field or field.strip() != field for field in fields
        ):
            raise ValueError("bounded scheduler evidence row is malformed")
        job_id, job_name, state, exit_code, nodes, elapsed, submitted, completed = fields
        if path.name != f"{job_id}.stage_i.sacct.txt":
            raise ValueError("bounded scheduler evidence filename differs from its job ID")
        expected_fields = (
            str(item["job_id"]),
            str(item["job_name"]),
            str(item["state"]),
            str(item["exit_code"]),
            str(item["nodes"]),
            str(item["elapsed_seconds"]),
        )
        if (job_id, job_name, state, exit_code, nodes, elapsed) != expected_fields:
            raise ValueError(
                f"bounded scheduler evidence differs for job {item['job_id']}"
            )
        try:
            submitted_time = datetime.fromisoformat(submitted.replace("Z", "+00:00"))
            completed_time = datetime.fromisoformat(completed.replace("Z", "+00:00"))
        except ValueError as error:
            raise ValueError("bounded scheduler evidence timestamp is invalid") from error
        if submitted_time.tzinfo is None:
            submitted_time = submitted_time.replace(tzinfo=timezone.utc)
        if completed_time.tzinfo is None:
            completed_time = completed_time.replace(tzinfo=timezone.utc)
        if (
            submitted_time.utcoffset() != timedelta(0)
            or completed_time.utcoffset() != timedelta(0)
            or submitted_time.isoformat() != item["submitted_utc"]
            or completed_time.isoformat() != item["completed_utc"]
        ):
            raise ValueError(
                f"bounded scheduler evidence timestamp differs for job {item['job_id']}"
            )
        retained.append(path)
    return tuple(retained)


def authenticate_legacy_profile_source_bundle(root: Path,
                                              args: argparse.Namespace) -> None:
    """Authenticate an embedded legacy F114-compatible source-bundle binding."""

    profile = args.authorized_next_segment_profile_json
    if not isinstance(profile, dict):
        raise ValueError("sole-profile authorization must be an object")
    path_value = profile.get("source_bundle")
    sha_value = profile.get("source_bundle_sha256")
    if path_value is None or sha_value is None:
        raise ValueError("legacy sole-profile source-bundle binding is mandatory")
    revisions = args.expected_source_bundle_verified_revisions_json
    if not isinstance(revisions, list) or not revisions:
        raise ValueError(
            "legacy sole-profile source-bundle revisions are mandatory"
        )
    repository = repository_root(initial_source_path())
    helper_revision_bound = False
    for revision in revisions:
        retained = git_run(
            repository,
            ["show", f"{revision}:{STAGE_I_RELATIVE.as_posix()}"],
            capture_output=True,
        )
        if (
            retained.returncode == 0
            and sha256_bytes(retained.stdout) == args.expected_stage_i_sha256
        ):
            helper_revision_bound = True
            break
    if not helper_revision_bound:
        raise ValueError(
            "legacy source-bundle revisions do not bind the Stage I helper"
        )
    bundle = require_root_member(root, path_value, "legacy sole-profile source bundle")
    if bundle.parent != root / "source-archives":
        raise ValueError(
            "legacy sole-profile source bundle must be a direct child of source-archives"
        )
    expected = require_sha256(sha_value, "legacy sole-profile source bundle SHA-256")
    require_confined_file_sha256(
        root,
        bundle.relative_to(root).as_posix(),
        expected,
        "legacy sole-profile source bundle",
        expected_mode=args.expected_source_bundle_mode,
    )
    authenticate_source_bundle_coverage(
        repository,
        bundle,
        expected,
        revisions,
        expected_mode=args.expected_source_bundle_mode,
    )


def authenticate_external_files(root: Path, args: argparse.Namespace
                                ) -> tuple[Path, Path | tuple[Path, ...]]:
    """Authenticate retained generator and all mode-selected external evidence."""

    generator = require_confined_file_sha256(
        root,
        args.generator_relative_path,
        args.expected_generator_sha256,
        "recost generator",
        expected_mode=args.expected_generator_mode,
    )
    if publication_mode(args) == "sole-profile":
        scheduler: Path | tuple[Path, ...] = require_confined_file_sha256(
            root,
            args.scheduler_relative_path,
            args.expected_scheduler_sha256,
            "scheduler evidence",
            expected_mode=args.expected_scheduler_mode,
        )
        authenticate_legacy_profile_source_bundle(root, args)
    else:
        packet = (
            validate_recost_recommendations(args.recost_recommendations_json, root, args)
            if schema2_evidence_mode(args)
            else validate_bounded_wave_authorization(
                args.authorized_bounded_wave_json, root, args
            )
        )
        source_bundle = require_confined_file_sha256(
            root,
            args.source_bundle_relative_path,
            args.expected_source_bundle_sha256,
            "source bundle",
            expected_mode=args.expected_source_bundle_mode,
        )
        authenticate_source_bundle_coverage(
            repository_root(initial_source_path()),
            source_bundle,
            args.expected_source_bundle_sha256,
            args.expected_source_bundle_verified_revisions_json,
            expected_mode=args.expected_source_bundle_mode,
        )
        profiles = packet[
            "recommended_next_profiles"
            if schema2_evidence_mode(args)
            else "authorized_next_profiles"
        ]
        assert isinstance(profiles, list)
        for index, profile in enumerate(profiles):
            assert isinstance(profile, dict)
            restart = profile["restart_file"]
            if restart is not None:
                authenticate_root_member_sha256(
                    root,
                    Path(str(restart)),
                    str(profile["restart_file_sha256"]),
                    f"bounded profile {index} restart",
                )
        scheduler = authenticate_bounded_scheduler_files(root, args)
    return generator, scheduler


@contextmanager
def promotion_lock(paths: dict[str, Path]):
    """Take the nonblocking cooperative Stage I root lock."""

    global _ACTIVE_MUTATION_LOCK, _ACTIVE_PUBLIC_ROOT
    root = absolute_path(paths["root"])
    lock = absolute_path(paths["lock"])
    expected_lock = root / f".mks24_stage_i_{EXECUTION_EPOCH_SLUG}.lock"
    if lock != expected_lock:
        raise ValueError(f"Stage I lock path differs: {lock}")
    if _ACTIVE_MUTATION_LOCK is not None or _ACTIVE_PUBLIC_ROOT is not None:
        raise ValueError("nested checkpoint mutation locks are forbidden")
    with bound_public_root(root) as root_guard:
        root_descriptor = root_guard.descriptor
        root_guard.assert_bound()
        created = False
        descriptor = None
        operation_error = None
        previous = os.umask(0)
        try:
            try:
                descriptor = os.open(
                    lock.name,
                    os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                    0o644,
                    dir_fd=root_descriptor,
                )
                created = True
            except FileExistsError:
                descriptor = os.open(
                    lock.name,
                    os.O_RDWR | os.O_NOFOLLOW,
                    dir_fd=root_descriptor,
                )
            except BaseException as error:
                operation_error = error
        finally:
            try:
                os.umask(previous)
            except BaseException as error:
                if operation_error is None:
                    operation_error = error
        if operation_error is not None:
            try:
                durably_classify_created_regular(
                    root_descriptor,
                    lock.name,
                    descriptor,
                    "Stage I lock",
                    expected_mode=0o644,
                    operation_error=operation_error,
                )
            finally:
                if descriptor is not None:
                    os.close(descriptor)
        assert descriptor is not None
        if created:
            try:
                durably_classify_created_regular(
                    root_descriptor,
                    lock.name,
                    descriptor,
                    "Stage I lock",
                    expected_mode=0o644,
                )
            except BaseException:
                os.close(descriptor)
                raise
        try:
            profile = os.fstat(descriptor)
            require_regular_mode_subset(
                profile,
                lock,
                "Stage I lock",
                maximum_mode=0o644,
                expected_links=1,
                expected_uid=os.geteuid(),
            )
            if profile.st_size != 0:
                raise ValueError(f"Stage I lock must be empty: {lock}")
        except BaseException:
            os.close(descriptor)
            raise
        root_guard.assert_bound()
        stream = os.fdopen(descriptor, "r+", encoding="utf-8")
        try:
            try:
                fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError as error:
                if error.errno not in (errno.EACCES, errno.EAGAIN):
                    raise
                raise ValueError(f"another Stage I mutation holds {lock}") from error
            opened_profile = os.fstat(stream.fileno())
            try:
                named_descriptor = os.open(
                    lock.name,
                    os.O_RDWR | os.O_NOFOLLOW,
                    dir_fd=root_descriptor,
                )
            except OSError as error:
                raise ValueError(
                    f"Stage I lock path changed while locking: {lock}"
                ) from error
            try:
                named_profile = os.fstat(named_descriptor)
            finally:
                os.close(named_descriptor)
            if (
                not stat.S_ISREG(named_profile.st_mode)
                or (named_profile.st_dev, named_profile.st_ino)
                != (opened_profile.st_dev, opened_profile.st_ino)
            ):
                raise ValueError(f"Stage I lock path changed while locking: {lock}")
            require_regular_mode_subset(
                named_profile,
                lock,
                "Stage I lock",
                maximum_mode=0o644,
                expected_links=1,
                expected_uid=os.geteuid(),
            )
            if named_profile.st_size != 0:
                raise ValueError(f"Stage I lock must be empty: {lock}")
            guard = MutationLockGuard(
                root_descriptor,
                lock.name,
                stream.fileno(),
                opened_profile,
                "Stage I lock",
            )
            guard.assert_bound()
            root_guard.assert_bound()
            _ACTIVE_PUBLIC_ROOT = root_guard
            _ACTIVE_MUTATION_LOCK = guard
            if stat.S_IMODE(opened_profile.st_mode) != 0o644:
                guard.assert_bound()
                os.fchmod(stream.fileno(), 0o644)
                os.fsync(stream.fileno())
                os.fsync(root_descriptor)
                guard.assert_bound()
            try:
                yield
            finally:
                try:
                    guard.assert_bound()
                    root_guard.assert_bound()
                finally:
                    _ACTIVE_MUTATION_LOCK = None
                    _ACTIVE_PUBLIC_ROOT = None
        finally:
            try:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
            finally:
                stream.close()


def copy_forensic(source: Path, destination: Path, expected: str, *,
                  expected_mode: int,
                  source_label: str = "staged recost artifact",
                  simulate_interruption_before_directory_fsync: bool = False) -> None:
    """Retain one independent authenticated forensic copy."""

    expected = require_sha256(expected, "recost forensic copy SHA-256")
    remove_forensic_temporaries(
        destination.parent,
        "recost forensic temporary",
    )
    with absolute_descriptor(source, source_label, flags=os.O_RDONLY) as source_descriptor:
        source_profile = os.fstat(source_descriptor)
        require_regular_profile(
            source_profile,
            source,
            source_label,
            expected_mode=expected_mode,
            expected_links=1,
        )
        with bound_parent_descriptor(destination, "recost forensic copy") as parent:
            require_mutation_authority_bound(parent)
            temporary = forensic_temporary_name(destination)
            destination_descriptor = None
            operation_error = None
            previous = os.umask(0)
            try:
                try:
                    destination_descriptor = os.open(
                        temporary,
                        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                        0o600,
                        dir_fd=parent,
                    )
                except BaseException as error:
                    operation_error = error
            finally:
                try:
                    os.umask(previous)
                except BaseException as error:
                    if operation_error is None:
                        operation_error = error
            try:
                durably_classify_created_regular(
                    parent,
                    temporary,
                    destination_descriptor,
                    "recost forensic copy temporary",
                    expected_mode=0o600,
                    operation_error=operation_error,
                )
            except BaseException:
                if destination_descriptor is not None:
                    os.close(destination_descriptor)
                raise
            assert destination_descriptor is not None
            digest = hashlib.sha256()
            destination_profile = os.fstat(destination_descriptor)
            try:
                while True:
                    block = os.read(source_descriptor, 1024 * 1024)
                    if not block:
                        break
                    digest.update(block)
                    write_descriptor_bytes(destination_descriptor, block)
                os.fchmod(destination_descriptor, 0o444)
                os.fsync(destination_descriptor)
                destination_profile = os.fstat(destination_descriptor)
                require_regular_profile(
                    destination_profile,
                    destination,
                    "recost forensic copy",
                    expected_mode=0o444,
                    expected_links=1,
                    expected_uid=os.geteuid(),
                )
            finally:
                os.close(destination_descriptor)
            if profile_security_binding(os.fstat(source_descriptor)) != (
                profile_security_binding(source_profile)
            ):
                raise ValueError("staged recost artifact security profile changed during copy")
            if digest.hexdigest() != expected:
                unlink_bound_entry(
                    parent,
                    temporary,
                    destination_profile,
                    "recost forensic copy",
                )
                raise ValueError("staged recost artifact checksum changed during forensic copy")
            require_bound_entry_security(
                parent,
                temporary,
                destination_profile,
                "recost forensic copy",
                expected_sha256=expected,
            )
            rename_bound_noreplace(
                parent,
                temporary,
                destination.name,
                destination_profile,
                "recost forensic publication",
                expected_sha256=expected,
            )
            if simulate_interruption_before_directory_fsync:
                raise ValueError("simulated interruption before forensic directory fsync")
            require_bound_entry_security(
                parent,
                destination.name,
                destination_profile,
                "recost forensic copy",
                expected_sha256=expected,
            )
            os.fsync(parent)
            require_mutation_authority_bound(parent)
            require_bound_entry_security(
                parent,
                destination.name,
                destination_profile,
                "recost forensic copy",
                expected_sha256=expected,
            )


def forensic_temporary_name(destination: Path) -> str:
    """Return the deterministic crash-recovery name for one forensic copy."""

    name = require_entry_name(destination.name, "recost forensic destination")
    return f".cgl-checkpoint-forensic-{sha256_bytes(name.encode())}.tmp"


def forensic_temporary_entries(directory: Path) -> list[Path]:
    """Return narrowly matched current and legacy forensic-copy temporaries."""

    require_managed_directory(directory, "recost forensic directory")
    pattern = re.compile(
        r"\.cgl-checkpoint-forensic-(?:[0-9a-f]{32}|[0-9a-f]{64}\.tmp)"
    )
    return [
        entry for entry in directory_entries(directory)
        if pattern.fullmatch(entry.name) is not None
    ]


def remove_forensic_temporaries(directory: Path, label: str) -> None:
    """Forensically retire authenticated partial forensic-copy temporaries."""

    for entry in forensic_temporary_entries(directory):
        with bound_parent_descriptor(entry, label) as parent:
            try:
                profile = os.stat(entry.name, dir_fd=parent, follow_symlinks=False)
            except FileNotFoundError:
                continue
            require_regular_profile(
                profile,
                entry,
                label,
                expected_links=1,
                expected_uid=os.geteuid(),
            )
            mode = stat.S_IMODE(profile.st_mode)
            if mode not in {0o000, 0o444, 0o600}:
                raise ValueError(
                    f"{label} mode is {mode:04o}, expected 0000, 0444, or 0600: {entry}"
                )
            unlink_bound_entry(parent, entry.name, profile, label)


def recover_or_copy_forensic(source: Path, destination: Path, expected: str, *,
                             expected_mode: int,
                             simulate_interruption_before_directory_fsync: bool
                             ) -> None:
    """Repair and durably bind one interrupted provisional forensic copy."""

    expected = require_sha256(expected, "recost forensic copy SHA-256")

    def authenticate_existing() -> bool:
        with bound_parent_descriptor(destination, "recost forensic recovery") as parent:
            try:
                descriptor = os.open(
                    destination.name,
                    os.O_RDONLY | os.O_NOFOLLOW,
                    dir_fd=parent,
                )
            except FileNotFoundError:
                return False
            try:
                profile = os.fstat(descriptor)
                require_regular_profile(
                    profile,
                    destination,
                    "recost forensic copy",
                    expected_mode=0o444,
                    expected_links=1,
                    expected_uid=os.geteuid(),
                )
                digest = sha256_descriptor(descriptor)
                if profile_security_binding(os.fstat(descriptor)) != (
                    profile_security_binding(profile)
                ):
                    raise ValueError(
                        "recost forensic copy descriptor security profile changed "
                        "during recovery"
                    )
                require_bound_entry_profile(
                    parent, destination.name, profile, "recost forensic copy"
                )
                if digest != expected:
                    unlink_bound_entry(
                        parent,
                        destination.name,
                        profile,
                        "recost forensic copy",
                        expected_sha256=digest,
                    )
                    return False
                require_bound_descriptor_security(
                    parent,
                    destination.name,
                    descriptor,
                    profile,
                    "recost forensic copy",
                    expected_sha256=expected,
                )
                os.fsync(descriptor)
                require_bound_descriptor_security(
                    parent,
                    destination.name,
                    descriptor,
                    profile,
                    "recost forensic copy",
                    expected_sha256=expected,
                )
                os.fsync(parent)
                require_bound_descriptor_security(
                    parent,
                    destination.name,
                    descriptor,
                    profile,
                    "recost forensic copy",
                    expected_sha256=expected,
                )
                require_mutation_authority_bound(parent)
                return True
            finally:
                os.close(descriptor)

    if not authenticate_existing():
        copy_forensic(
            source,
            destination,
            expected,
            expected_mode=expected_mode,
            source_label="canonical recost artifact",
            simulate_interruption_before_directory_fsync=(
                simulate_interruption_before_directory_fsync
            ),
        )
    if not authenticate_existing():
        raise ValueError("recost forensic copy disappeared after durable publication")


def descriptor_identity(descriptor: int) -> tuple[int, int]:
    """Return the stable device/inode identity for one open descriptor."""

    value = os.fstat(descriptor)
    return value.st_dev, value.st_ino


def require_named_artifact_identity(path: Path, args: argparse.Namespace, *,
                                    expected_links: int,
                                    expected_identity: tuple[int, int] | None,
                                    allow_expired_v2: bool = False) -> None:
    """Require one artifact name still references the selected inode."""

    read_artifact(
        path,
        args,
        expected_links=expected_links,
        allow_expired_v2=allow_expired_v2,
    )
    with regular_descriptor(
        path,
        "recost artifact",
        expected_mode=args.expected_artifact_mode,
        expected_links=expected_links,
    ) as descriptor:
        if (
            expected_identity is not None
            and descriptor_identity(descriptor) != expected_identity
        ):
            raise ValueError(f"recost artifact inode identity has changed: {path}")


def link_after_empty_queue(root: Path, offline: bool, queue_file: str | None,
                           directory_descriptor: int,
                           staged_descriptor: int,
                           paths: dict[str, Path],
                           args: argparse.Namespace,
                           expected_identity: tuple[int, int]) -> None:
    """Check the live queue and immediately attempt same-directory publication."""

    require_artifact_namespace_in_directory(directory_descriptor, paths, "staged")
    require_named_artifact_identity_in_directory(
        directory_descriptor,
        paths["staged"].name,
        paths["staged"],
        args,
        expected_links=1,
        expected_identity=expected_identity,
    )
    require_empty_queue(root, offline, queue_file)
    try:
        link_descriptor_noreplace(
            staged_descriptor,
            directory_descriptor,
            paths["canonical"].name,
            expected_identity,
            "canonical recost publication link",
        )
    except OSError as error:
        raise LinkAttemptError(f"recost publication link failed: {error}") from error
    if args.simulate_link_return_failure:
        raise LinkAttemptError("simulated failure immediately after recost publication link")


def single_link_replacement_name(transaction_id: str) -> str:
    """Return the transaction-bound temporary name for one canonical copy."""

    if not transaction_id or "/" in transaction_id or Path(transaction_id).name != transaction_id:
        raise ValueError("recost transaction ID is invalid for single-link publication")
    return f".cgl-checkpoint-canonical-{transaction_id}.copy"


def identity_binding(identity: tuple[int, int]) -> dict[str, int]:
    """Return one stable JSON inode-identity binding."""

    return {"device": identity[0], "inode": identity[1]}


def recovery_identity(record: dict[str, object], key: str, label: str) -> tuple[int, int]:
    """Return one exact inode identity from a recovery journal."""

    value = record.get(key)
    if not isinstance(value, dict) or frozenset(value) != {"device", "inode"}:
        raise ValueError(f"{label} identity is invalid")
    return (
        require_integer(value["device"], f"{label} device"),
        require_integer(value["inode"], f"{label} inode", minimum=1),
    )


def prepare_single_link_copy(
    directory_descriptor: int,
    paths: dict[str, Path],
    args: argparse.Namespace,
    replacement: str,
    *,
    expected_staged_identity: tuple[int, int],
    allow_expired_v2: bool = False,
) -> tuple[int, int]:
    """Persist one exact single-link copy before its journaled exchange."""

    require_linked_pair_in_directory(
        directory_descriptor,
        paths,
        args,
        expected_identity=expected_staged_identity,
        allow_expired_v2=allow_expired_v2,
    )
    replacement = require_entry_name(replacement, "canonical recost single-link copy")
    try:
        os.stat(replacement, dir_fd=directory_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        pass
    else:
        raise ValueError("canonical recost single-link copy already exists")
    staged_descriptor = os.open(
        paths["staged"].name,
        os.O_RDONLY | os.O_NOFOLLOW,
        dir_fd=directory_descriptor,
    )
    try:
        staged_profile = os.fstat(staged_descriptor)
        require_regular_profile(
            staged_profile,
            paths["staged"],
            "staged recost publication link",
            expected_mode=args.expected_artifact_mode,
            expected_links=2,
            expected_uid=os.geteuid(),
        )
        retained = read_descriptor_bytes(staged_descriptor)
    finally:
        os.close(staged_descriptor)
    if sha256_bytes(retained) != args.expected_artifact_sha256:
        raise ValueError("staged recost artifact checksum changed before canonical copy")
    replacement_profile = write_bound_exclusive(
        directory_descriptor,
        replacement,
        retained,
        args.expected_artifact_mode,
        "canonical recost single-link copy",
    )
    return profile_identity(replacement_profile)


def validate_single_link_transition(
    directory_descriptor: int,
    paths: dict[str, Path],
    args: argparse.Namespace,
    replacement: str,
    staged_identity: tuple[int, int],
    canonical_identity: tuple[int, int],
    *,
    allow_expired_v2: bool = False,
) -> str:
    """Authenticate every durable state of the journaled no-exchange transition."""

    replacement = require_entry_name(replacement, "canonical recost single-link copy")

    def selected(name: str) -> os.stat_result | None:
        try:
            return os.stat(name, dir_fd=directory_descriptor, follow_symlinks=False)
        except FileNotFoundError:
            return None

    staged = selected(paths["staged"].name)
    canonical = selected(paths["canonical"].name)
    copy = selected(replacement)

    def require_selected(
        profile: os.stat_result | None,
        name: str,
        path: Path,
        identity: tuple[int, int],
        links: set[int],
        label: str,
    ) -> None:
        if profile is None or profile_identity(profile) != identity:
            raise ValueError(f"{label} does not select its journaled inode")
        if profile.st_nlink not in links:
            raise ValueError(
                f"{label} has {profile.st_nlink} links, expected one of "
                + ", ".join(str(value) for value in sorted(links))
            )
        require_named_artifact_identity_in_directory(
            directory_descriptor,
            name,
            path,
            args,
            expected_links=profile.st_nlink,
            expected_identity=identity,
            allow_expired_v2=allow_expired_v2,
        )

    visible = (
        None if staged is None else profile_identity(staged),
        None if canonical is None else profile_identity(canonical),
        None if copy is None else profile_identity(copy),
    )
    if staged is not None and visible[0] != staged_identity:
        raise ValueError(
            "staged recost publication link does not select its journaled inode"
        )
    if canonical is not None and visible[1] not in {
        staged_identity,
        canonical_identity,
    }:
        raise ValueError(
            "canonical recost publication does not select a journaled inode"
        )
    if copy is not None and visible[2] != canonical_identity:
        raise ValueError(
            "canonical recost single-link copy does not select its journaled inode"
        )
    if visible == (staged_identity, staged_identity, canonical_identity):
        require_selected(
            staged,
            paths["staged"].name,
            paths["staged"],
            staged_identity,
            {2, 3},
            "staged recost publication link",
        )
        require_selected(
            canonical,
            paths["canonical"].name,
            paths["canonical"],
            staged_identity,
            {2, 3},
            "canonical recost predecessor",
        )
        if staged is None or canonical is None or staged.st_nlink != canonical.st_nlink:
            raise ValueError("canonical predecessor links changed during retirement")
        require_selected(
            copy,
            replacement,
            paths["accounting"] / replacement,
            canonical_identity,
            {1},
            "canonical recost single-link copy",
        )
        return "copy-prepared"
    if visible == (staged_identity, None, canonical_identity):
        require_selected(
            staged,
            paths["staged"].name,
            paths["staged"],
            staged_identity,
            {2},
            "staged recost publication link",
        )
        require_selected(
            copy,
            replacement,
            paths["accounting"] / replacement,
            canonical_identity,
            {1},
            "canonical recost single-link copy",
        )
        return "canonical-retired"
    if visible == (staged_identity, canonical_identity, canonical_identity):
        require_selected(
            staged,
            paths["staged"].name,
            paths["staged"],
            staged_identity,
            {2},
            "staged recost publication link",
        )
        require_selected(
            canonical,
            paths["canonical"].name,
            paths["canonical"],
            canonical_identity,
            {2},
            "canonical recost publication",
        )
        require_selected(
            copy,
            replacement,
            paths["accounting"] / replacement,
            canonical_identity,
            {2},
            "canonical recost single-link copy",
        )
        return "canonical-linked"
    if visible == (staged_identity, canonical_identity, None):
        require_selected(
            staged,
            paths["staged"].name,
            paths["staged"],
            staged_identity,
            {2, 3},
            "staged recost publication link",
        )
        require_selected(
            canonical,
            paths["canonical"].name,
            paths["canonical"],
            canonical_identity,
            {1},
            "canonical recost publication",
        )
        return "canonical-published"
    if visible == (None, canonical_identity, None):
        require_selected(
            canonical,
            paths["canonical"].name,
            paths["canonical"],
            canonical_identity,
            {1},
            "canonical recost publication",
        )
        return "complete"
    raise ValueError("canonical recost no-exchange transition state is invalid")


def complete_single_link_transition(
    directory_descriptor: int,
    paths: dict[str, Path],
    args: argparse.Namespace,
    replacement: str,
    staged_identity: tuple[int, int],
    canonical_identity: tuple[int, int],
    *,
    allow_expired_v2: bool = False,
    simulate_post_exchange_failure: bool = False,
) -> None:
    """Complete one journaled forward-only canonical publication transition."""

    state = validate_single_link_transition(
        directory_descriptor,
        paths,
        args,
        replacement,
        staged_identity,
        canonical_identity,
        allow_expired_v2=allow_expired_v2,
    )
    if state == "copy-prepared":
        predecessor = os.stat(
            paths["canonical"].name,
            dir_fd=directory_descriptor,
            follow_symlinks=False,
        )
        unlink_bound_entry(
            directory_descriptor,
            paths["canonical"].name,
            predecessor,
            "canonical recost predecessor",
            expected_sha256=args.expected_artifact_sha256,
        )
        state = validate_single_link_transition(
            directory_descriptor,
            paths,
            args,
            replacement,
            staged_identity,
            canonical_identity,
            allow_expired_v2=allow_expired_v2,
        )
    if state in {"canonical-retired", "canonical-linked"}:
        replacement_profile = os.stat(
            replacement, dir_fd=directory_descriptor, follow_symlinks=False
        )
        rename_bound_noreplace(
            directory_descriptor,
            replacement,
            paths["canonical"].name,
            replacement_profile,
            "canonical recost single-link publication",
            expected_sha256=args.expected_artifact_sha256,
        )
        state = validate_single_link_transition(
            directory_descriptor,
            paths,
            args,
            replacement,
            staged_identity,
            canonical_identity,
            allow_expired_v2=allow_expired_v2,
        )
    if state == "canonical-published" and simulate_post_exchange_failure:
        raise ValueError("simulated interruption after canonical copy publication")
    if state == "canonical-published":
        staged_profile = os.stat(
            paths["staged"].name,
            dir_fd=directory_descriptor,
            follow_symlinks=False,
        )
        unlink_bound_entry(
            directory_descriptor,
            paths["staged"].name,
            staged_profile,
            "staged recost publication link",
            expected_sha256=args.expected_artifact_sha256,
        )
    state = validate_single_link_transition(
        directory_descriptor,
        paths,
        args,
        replacement,
        staged_identity,
        canonical_identity,
        allow_expired_v2=allow_expired_v2,
    )
    if state != "complete":
        raise ValueError("canonical recost no-exchange transition did not complete")
    os.fsync(directory_descriptor)
    validate_single_link_transition(
        directory_descriptor,
        paths,
        args,
        replacement,
        staged_identity,
        canonical_identity,
        allow_expired_v2=allow_expired_v2,
    )
    require_mutation_authority_bound(directory_descriptor)


def scheduler_audit_binding(root: Path, args: argparse.Namespace,
                            scheduler: Path | tuple[Path, ...]) -> object:
    """Return the exact mode-selected scheduler audit binding."""

    if publication_mode(args) == "sole-profile":
        if not isinstance(scheduler, Path):
            raise ValueError("sole-profile scheduler binding is invalid")
        return {
            "path": str(scheduler),
            "sha256": args.expected_scheduler_sha256,
            "mode": f"{args.expected_scheduler_mode:04o}",
        }
    if not isinstance(scheduler, tuple):
        raise ValueError("bounded-wave scheduler bindings are invalid")
    evidence = validate_bounded_scheduler_evidence(
        args.bounded_wave_scheduler_evidence_json
    )
    if len(scheduler) != len(evidence):
        raise ValueError("bounded-wave scheduler path count differs")
    return [
        {
            **item,
            "path": str(path),
            "mode": f"{args.expected_scheduler_mode:04o}",
        }
        for item, path in zip(evidence, scheduler)
    ]


def source_bundle_audit_binding(root: Path, args: argparse.Namespace
                                ) -> dict[str, object]:
    """Return the exact mode-selected source-bundle audit binding."""

    if publication_mode(args) == "sole-profile":
        profile = args.authorized_next_segment_profile_json
        assert isinstance(profile, dict)
        return {
            "path": str(profile["source_bundle"]),
            "sha256": profile["source_bundle_sha256"],
            "mode": f"{args.expected_source_bundle_mode:04o}",
            "verified_revisions": args.expected_source_bundle_verified_revisions_json,
        }
    return {
        "path": str(root / safe_relative_path(
            args.source_bundle_relative_path, "source bundle relative path"
        )),
        "sha256": args.expected_source_bundle_sha256,
        "mode": f"{args.expected_source_bundle_mode:04o}",
        "verified_revisions": args.expected_source_bundle_verified_revisions_json,
    }


def publication_audit_mode(args: argparse.Namespace) -> int:
    """Return the immutable audit mode required by the selected lifecycle."""

    return 0o444 if schema2_evidence_mode(args) else 0o644


def audit_record(paths: dict[str, Path], args: argparse.Namespace,
                 source: Path, generator: Path,
                 scheduler: Path | tuple[Path, ...],
                 forensic: Path, transaction_id: str, *,
                 allow_expired_v2: bool = False,
                 published_utc: object | None = None) -> dict[str, object]:
    """Build the finalized durable publication audit."""

    if published_utc is None:
        publication_timestamp = datetime.now(timezone.utc).isoformat()
    else:
        require_utc_timestamp(published_utc, "publication audit timestamp")
        assert isinstance(published_utc, str)
        publication_timestamp = published_utc
    generalized_context = (
        selected_v2_publication_context(
            paths, args, allow_expired=allow_expired_v2
        )
        if publication_mode(args) == "v2"
        else None
    )
    record = {
        "schema_version": 1,
        "record_type": (
            "stage-i-recost-recommendation-publication-audit"
            if schema2_evidence_mode(args)
            else "observed-publication"
        ),
        "execution_epoch": EXECUTION_EPOCH,
        "transaction_id": transaction_id,
        "published_utc": publication_timestamp,
        "artifact": {
            "path": str(paths["canonical"]),
            "sha256": args.expected_artifact_sha256,
            "mode": f"{args.expected_artifact_mode:04o}",
            "links": 1,
        },
    }
    if publication_mode(args) == "sole-profile":
        record["authorized_sole_next_segment_profile"] = (
            args.authorized_next_segment_profile_json
        )
    elif schema2_evidence_mode(args):
        record["recost_recommendations"] = args.recost_recommendations_json
        assert generalized_context is not None
        record["independent_review"] = schema2_independent_review(
            paths["root"],
            paths["canonical"].name,
            args,
            artifact_generated_utc=generalized_context["generated_utc"],
            publication_published_utc=publication_timestamp,
        )
        record["authority"] = {
            "action_authority": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
    else:
        record["authorized_bounded_wave"] = args.authorized_bounded_wave_json
    record.update({
        "counts": expected_counts(args),
        "generator": {
            "path": str(generator),
            "sha256": args.expected_generator_sha256,
            "mode": f"{args.expected_generator_mode:04o}",
        },
        "scheduler_evidence": scheduler_audit_binding(
            paths["root"], args, scheduler
        ),
    })
    record["source_bundle"] = source_bundle_audit_binding(paths["root"], args)
    record.update({
        "stage_i_helper": {
            "path": str(source.parent / "cgl_lf_stage_i.py"),
            "sha256": args.expected_stage_i_sha256,
            "committed": True,
            "reconcile_execution": "descriptor",
        },
        "utility": {
            "path": str(source),
            "sha256": args.expected_utility_sha256,
            "execution": "authenticated-descriptor",
            "committed": not args.allow_local_root,
        },
        "forensic_copy": {
            "path": str(forensic),
            "sha256": args.expected_artifact_sha256,
            "mode": "0444",
            "links": 1,
        },
        "publication": (
            "same-directory-link-fsync-copy-exchange-forensic-retirement-fsync"
        ),
    })
    if publication_mode(args) == "v2":
        record["generator"]["revision"] = args.expected_generator_revision
        record["stage_i_helper"]["revision"] = args.expected_stage_i_revision
        assert generalized_context is not None
        record["generalized_publication_context"] = generalized_context
        record["forensic_copy"]["generalized_vectors"] = "exact-artifact-payload"
    return record


def adoption_audit_record(paths: dict[str, Path], args: argparse.Namespace,
                          source: Path, generator: Path,
                          scheduler: Path | tuple[Path, ...],
                          forensic: Path, transaction_id: str, *,
                          allow_expired_v2: bool = False) -> dict[str, object]:
    """Build a present-time audit for a legacy canonical-only artifact."""

    record = {
        "schema_version": 1,
        "record_type": "legacy-canonical-adoption",
        "execution_epoch": EXECUTION_EPOCH,
        "transaction_id": transaction_id,
        "adopted_utc": utc_now(),
        "artifact": {
            "path": str(paths["canonical"]),
            "sha256": args.expected_artifact_sha256,
            "mode": f"{args.expected_artifact_mode:04o}",
            "links": 1,
        },
    }
    if publication_mode(args) == "sole-profile":
        record["authorized_sole_next_segment_profile"] = (
            args.authorized_next_segment_profile_json
        )
    else:
        record["authorized_bounded_wave"] = args.authorized_bounded_wave_json
    record.update({
        "counts": expected_counts(args),
        "generator": {
            "path": str(generator),
            "sha256": args.expected_generator_sha256,
            "mode": f"{args.expected_generator_mode:04o}",
        },
        "scheduler_evidence": scheduler_audit_binding(
            paths["root"], args, scheduler
        ),
    })
    record["source_bundle"] = source_bundle_audit_binding(paths["root"], args)
    record.update({
        "stage_i_helper": {
            "path": str(source.parent / "cgl_lf_stage_i.py"),
            "sha256": args.expected_stage_i_sha256,
            "committed": True,
            "reconcile_execution": "descriptor",
        },
        "utility": {
            "path": str(source),
            "sha256": args.expected_utility_sha256,
            "execution": "authenticated-descriptor",
            "committed": not args.allow_local_root,
        },
        "forensic_copy": {
            "path": str(forensic),
            "sha256": args.expected_artifact_sha256,
            "mode": "0444",
            "links": 1,
        },
        "original_publication_transition_observed": False,
        "original_publication_method": "unknown",
        "original_publisher": "unknown",
        "present_authentication": {
            "canonical_artifact": "authenticated-under-stage-i-lock",
            "no_queued_cgl_jobs": "required-under-stage-i-lock",
            "reconciliation": "clean-required-under-stage-i-lock",
        },
    })
    if publication_mode(args) == "v2":
        record["generator"]["revision"] = args.expected_generator_revision
        record["stage_i_helper"]["revision"] = args.expected_stage_i_revision
        record["generalized_publication_context"] = selected_v2_publication_context(
            paths, args, allow_expired=allow_expired_v2
        )
        record["forensic_copy"]["generalized_vectors"] = "exact-artifact-payload"
    return record


def validate_audit(paths: dict[str, Path], args: argparse.Namespace,
                   source: Path, generator: Path,
                   scheduler: Path | tuple[Path, ...], *,
                   allow_expired_v2: bool = False) -> None:
    """Validate the finalized audit and its retained forensic copy."""

    with regular_descriptor(
        paths["audit"],
        "recost publication audit",
        expected_mode=publication_audit_mode(args),
        expected_links=1,
    ) as descriptor:
        retained = read_descriptor_bytes(descriptor)
    try:
        audit = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("recost publication audit is not valid JSON") from error
    if not isinstance(audit, dict):
        raise ValueError("recost publication audit must be an object")
    forensic_record = audit.get("forensic_copy")
    if not isinstance(forensic_record, dict):
        raise ValueError("recost publication audit forensic binding is invalid")
    transaction_id = audit.get("transaction_id")
    if not isinstance(transaction_id, str) or not transaction_id:
        raise ValueError("recost publication audit transaction ID is invalid")
    forensic = Path(str(forensic_record.get("path", "")))
    expected_forensic = (
        paths["recost_forensics"]
        / f"{transaction_id}.{paths['canonical'].name}.forensic"
    )
    if forensic != expected_forensic:
        raise ValueError("recost forensic copy path does not match its transaction")
    record_type = audit.get("record_type")
    if record_type in {
        "observed-publication",
        "stage-i-recost-recommendation-publication-audit",
    }:
        expected = audit_record(
            paths,
            args,
            source,
            generator,
            scheduler,
            forensic,
            transaction_id,
            allow_expired_v2=allow_expired_v2,
            published_utc=audit.get("published_utc"),
        )
        timestamp = "published_utc"
        timestamp_label = "publication audit timestamp"
    elif record_type == "legacy-canonical-adoption":
        expected = adoption_audit_record(
            paths,
            args,
            source,
            generator,
            scheduler,
            forensic,
            transaction_id,
            allow_expired_v2=allow_expired_v2,
        )
        timestamp = "adopted_utc"
        timestamp_label = "adoption audit timestamp"
    else:
        raise ValueError("recost publication audit record type is invalid")
    if frozenset(audit) != frozenset(expected):
        raise ValueError("recost publication audit schema differs")
    require_utc_timestamp(audit.get(timestamp), timestamp_label)
    for key in expected:
        if key == timestamp:
            continue
        if audit.get(key) != expected[key]:
            raise ValueError(f"recost publication audit {key} binding differs")
    require_exact_directory_entries(
        paths["recost_forensics"],
        "recost forensic directory",
        {forensic},
        managed=True,
    )
    require_file_sha256(
        forensic,
        args.expected_artifact_sha256,
        "recost forensic copy",
        expected_mode=0o444,
        expected_links=1,
    )


def verify_staged_state(paths: dict[str, Path], root_dir: Path, root: Path,
                        offline: bool, args: argparse.Namespace, *,
                        stage_i_lock_held: bool = False) -> None:
    """Validate a staged-only recost boundary without mutation."""

    validate_fixture_options(args, offline)
    require_empty_queue(root, offline, args.squeue_file)
    require_empty_directory(
        paths["stage_i_transactions"], "Stage I transaction directory", trusted=True
    )
    require_empty_directory(
        paths["recost_transactions"],
        "recost transaction directory",
        allow_absent=True,
        managed=True,
    )
    require_artifact_namespace(paths, "staged")
    authenticate_publication_sources(root_dir, args)
    authenticate_external_files(root, args)
    read_artifact(paths["staged"], args)
    reexecute_v2_recost(
        paths,
        root_dir,
        root,
        args,
        paths["staged"],
        stage_i_lock_held=stage_i_lock_held,
    )
    run_reconcile(
        root_dir, root, offline, args, stage_i_lock_held=stage_i_lock_held
    )


def verify_retirable_staged_state(paths: dict[str, Path], root_dir: Path, root: Path,
                                  offline: bool, args: argparse.Namespace) -> None:
    """Validate staged-only cleanup without consuming freshness or HEAD equality."""

    validate_fixture_options(args, offline)
    require_empty_queue(root, offline, args.squeue_file)
    require_empty_directory(
        paths["stage_i_transactions"], "Stage I transaction directory", trusted=True
    )
    require_empty_directory(
        paths["recost_transactions"],
        "recost transaction directory",
        allow_absent=True,
        managed=True,
    )
    require_artifact_namespace(paths, "staged")
    authenticate_publication_sources(
        root_dir, args, require_generator_at_head=False
    )
    authenticate_external_files(root, args)
    read_artifact(paths["staged"], args, allow_expired_v2=True)
    run_reconcile(root_dir, root, offline, args)


def verify_promoted_state(paths: dict[str, Path], root_dir: Path, root: Path,
                          offline: bool, args: argparse.Namespace) -> None:
    """Validate canonical history after pre-link freshness has been consumed."""

    validate_fixture_options(args, offline)
    require_empty_queue(root, offline, args.squeue_file)
    require_empty_directory(
        paths["stage_i_transactions"], "Stage I transaction directory", trusted=True
    )
    require_empty_directory(
        paths["recost_transactions"], "recost transaction directory", managed=True
    )
    require_managed_directory(paths["recost_forensics_root"], "recost forensic root directory")
    require_managed_directory(paths["recost_forensics"], "recost forensic directory")
    require_artifact_namespace(paths, "promoted")
    authenticate_publication_sources(
        root_dir, args, require_generator_at_head=False
    )
    generator, scheduler = authenticate_external_files(root, args)
    read_artifact(paths["canonical"], args, allow_expired_v2=True)
    validate_audit(
        paths,
        args,
        initial_source_path(),
        generator,
        scheduler,
        allow_expired_v2=True,
    )
    run_reconcile(root_dir, root, offline, args)


def remove_bound_recost_journal(parent: int, journal: Path,
                                expected_record: dict[str, object],
                                label: str) -> None:
    """Retire one exact journal through its continuously bound transaction store."""

    for temporary in json_temporary_names_in_directory(parent, journal.name):
        temporary_profile = os.stat(temporary, dir_fd=parent, follow_symlinks=False)
        if temporary_profile.st_nlink != 1:
            recover_linked_json_publication_in_directory(
                parent,
                journal.parent,
                journal.name,
                f"{label} publication recovery",
                expected_mode=0o644,
            )
            break
    names = directory_entry_names(parent)
    if journal.name not in names:
        if names:
            raise ValueError(f"{label} transaction directory changed before cleanup")
        require_mutation_authority_bound(parent)
        return
    if names != [journal.name]:
        raise ValueError(f"{label} transaction directory changed before cleanup")
    profile, digest, _, record = read_bound_json_object(
        parent, journal.name, journal, label, expected_mode=0o644
    )
    if record != expected_record:
        raise ValueError(f"{label} changed before cleanup")
    unlink_bound_entry(
        parent, journal.name, profile, label, expected_sha256=digest
    )
    if directory_entry_names(parent):
        raise ValueError(f"{label} transaction directory changed after cleanup")
    require_mutation_authority_bound(parent)


def remove_prelink_records_bound(parent: int, journal: Path, forensic: Path,
                                 record: dict[str, object]) -> None:
    """Remove retry-safe prelink records while retaining transaction authority."""

    if entry_exists(forensic):
        unlink_durable(forensic)
    remove_bound_recost_journal(
        parent, journal, record, "recost prelink transaction journal"
    )


def recost_journal(paths: dict[str, Path], transaction_descriptor: int | None = None
                   ) -> tuple[Path, dict[str, object]]:
    """Authenticate or exactly recover the sole journal through one directory fd."""

    directory = paths["recost_transactions"]
    if transaction_descriptor is None:
        with bound_directory_descriptor(
            directory, "recost transaction directory"
        ) as selected:
            return recost_journal(paths, selected)
    parent = transaction_descriptor
    require_mutation_authority_bound(parent)
    temporary_pattern = re.compile(
        r"\.(?P<target>.+\.json)\.[0-9]+\.[0-9a-f]{32}\.tmp"
    )
    names = directory_entry_names(parent)
    public_names = [name for name in names if name.endswith(".json")]
    temporary_names = [
        name for name in names if temporary_pattern.fullmatch(name) is not None
    ]
    for name in temporary_names:
        target_name = json_temporary_target_name(name)
        assert target_name is not None
        temporary_profile = os.stat(name, dir_fd=parent, follow_symlinks=False)
        if temporary_profile.st_nlink != 1:
            recover_linked_json_publication_in_directory(
                parent,
                directory,
                target_name,
                "recost transaction journal publication recovery",
                expected_mode=0o644,
            )
    names = directory_entry_names(parent)
    public_names = [name for name in names if name.endswith(".json")]
    temporary_names = [
        name for name in names if temporary_pattern.fullmatch(name) is not None
    ]
    unknown = sorted(set(names) - set(public_names) - set(temporary_names))
    if len(public_names) > 1 or unknown:
        raise ValueError(
            "linked-pair recovery requires exactly one recost journal: "
            + ", ".join(str(directory / name) for name in names)
        )

    def parse_candidate(name: str, target_name: str, label: str, *,
                        expected_mode: int | None = None
                        ) -> tuple[os.stat_result, str, bytes, dict[str, object]]:
        candidate = read_bound_json_object(
            parent, name, directory / name, label, expected_mode=expected_mode
        )
        if candidate[3].get("transaction_id") != Path(target_name).stem:
            raise ValueError("recost transaction journal ID differs from its filename")
        return candidate

    valid_temporaries: list[
        tuple[str, str, os.stat_result, str, bytes, dict[str, object]]
    ] = []
    invalid_temporaries: list[str] = []
    for name in temporary_names:
        match = temporary_pattern.fullmatch(name)
        assert match is not None
        target_name = match.group("target")
        try:
            profile, digest, retained, record = parse_candidate(
                name, target_name, "recost transaction recovery temporary"
            )
        except ValueError as error:
            if (
                "not valid JSON" in str(error)
                or "mode-0000 content cannot be authenticated" in str(error)
            ):
                invalid_temporaries.append(name)
                continue
            raise
        valid_temporaries.append(
            (name, target_name, profile, digest, retained, record)
        )
    if len(valid_temporaries) > 1:
        raise ValueError("multiple valid recost transaction recovery temporaries were preserved")

    public_name = public_names[0] if public_names else None
    public_candidate = None
    public_error: BaseException | None = None
    if public_name is not None:
        try:
            public_candidate = parse_candidate(
                public_name,
                public_name,
                "recost transaction journal",
                expected_mode=0o644,
            )
        except BaseException as error:
            public_error = error

    if public_candidate is not None:
        public_profile, public_digest, public_retained, record = public_candidate
        if valid_temporaries:
            temporary = valid_temporaries[0]
            if temporary[1] != public_name:
                raise ValueError(
                    "valid recost transaction recovery temporary targets a "
                    "different authenticated public journal"
                )
            if temporary[4] == public_retained:
                unlink_bound_entry(
                    parent,
                    temporary[0],
                    temporary[2],
                    "duplicate recost transaction recovery temporary",
                    expected_sha256=temporary[3],
                )
            else:
                replace_bound_entry_forward(
                    parent,
                    temporary[0],
                    public_name,
                    temporary[2],
                    public_profile,
                    "recost transaction journal recovery",
                    source_sha256=temporary[3],
                    target_sha256=public_digest,
                )
                public_candidate = parse_candidate(
                    public_name,
                    public_name,
                    "recost transaction journal",
                    expected_mode=0o644,
                )
                public_profile, public_digest, public_retained, record = public_candidate
        remove_json_temporaries_in_directory(
            parent,
            directory,
            invalid_temporaries,
            "recost transaction temporary",
        )
        require_bound_entry_security(
            parent,
            public_name,
            public_profile,
            "recost transaction journal",
            expected_sha256=public_digest,
        )
        if directory_entry_names(parent) != [public_name]:
            raise ValueError("recost transaction directory changed during journal recovery")
        require_mutation_authority_bound(parent)
        return directory / public_name, record

    if len(valid_temporaries) != 1:
        if public_error is not None:
            raise public_error
        remove_json_temporaries_in_directory(
            parent,
            directory,
            invalid_temporaries,
            "recost transaction temporary",
        )
        if directory_entry_names(parent):
            raise ValueError("linked-pair recovery requires exactly one recost journal")
        require_mutation_authority_bound(parent)
        raise NoRecostJournal("recost transaction directory has no durable journal")

    temporary_name, target_name, temporary_profile, temporary_digest, _, _ = (
        valid_temporaries[0]
    )
    if public_name is not None and public_name != target_name:
        raise ValueError("recost recovery temporary targets a different public journal")
    if public_name is None:
        rename_bound_noreplace(
            parent,
            temporary_name,
            target_name,
            temporary_profile,
            "recost transaction journal recovery",
            expected_sha256=temporary_digest,
        )
    else:
        public_profile = os.stat(public_name, dir_fd=parent, follow_symlinks=False)
        require_regular_mode_subset(
            public_profile,
            directory / public_name,
            "forged recost transaction journal",
            maximum_mode=0o644,
            expected_links=public_profile.st_nlink,
            expected_uid=os.geteuid(),
        )
        public_digest = require_bound_entry_security(
            parent,
            public_name,
            public_profile,
            "forged recost transaction journal",
        )
        assert public_digest is not None
        replace_bound_entry_forward(
            parent,
            temporary_name,
            public_name,
            temporary_profile,
            public_profile,
            "recost transaction journal recovery",
            source_sha256=temporary_digest,
            target_sha256=public_digest,
        )
        recovered = parse_candidate(
            public_name,
            public_name,
            "recost transaction journal",
            expected_mode=0o644,
        )
        os.fsync(parent)
        require_bound_entry_security(
            parent,
            public_name,
            recovered[0],
            "recost transaction journal",
            expected_sha256=recovered[1],
        )
        require_mutation_authority_bound(parent)
    remove_json_temporaries_in_directory(
        parent,
        directory,
        invalid_temporaries,
        "recost transaction temporary",
    )
    recovered = parse_candidate(
        target_name,
        target_name,
        "recost transaction journal",
        expected_mode=0o644,
    )
    if directory_entry_names(parent) != [target_name]:
        raise ValueError("recost transaction directory changed during journal recovery")
    require_mutation_authority_bound(parent)
    return directory / target_name, recovered[3]


def recovery_bindings(paths: dict[str, Path], args: argparse.Namespace,
                      transaction_id: str, forensic: Path, *,
                      allow_expired_v2: bool = False,
                      staged_identity: tuple[int, int] | None = None
                      ) -> dict[str, object]:
    """Return every authorization-relevant immutable journal binding."""

    record = {
        "schema_version": 1,
        "execution_epoch": EXECUTION_EPOCH,
        "transaction_id": transaction_id,
        "staged_path": str(paths["staged"]),
        "canonical_path": str(paths["canonical"]),
        "artifact_sha256": args.expected_artifact_sha256,
        "artifact_mode": f"{args.expected_artifact_mode:04o}",
        "artifact_epoch_pointer": args.artifact_epoch_pointer,
        "single_link_replacement_name": single_link_replacement_name(transaction_id),
    }
    if staged_identity is not None:
        record["staged_inode_identity"] = {
            "device": staged_identity[0],
            "inode": staged_identity[1],
        }
    if publication_mode(args) == "sole-profile":
        record.update({
            "artifact_profile_pointer": args.artifact_profile_pointer,
            "artifact_counts_pointer": args.artifact_counts_pointer,
            "artifact_stage_i_sha256_pointer": args.artifact_stage_i_sha256_pointer,
            "artifact_generator_sha256_pointer": args.artifact_generator_sha256_pointer,
            "artifact_scheduler_sha256_pointer": args.artifact_scheduler_sha256_pointer,
            "authorized_sole_next_segment_profile": (
                args.authorized_next_segment_profile_json
            ),
        })
    else:
        record.update({
            "artifact_authorization_pointer": args.artifact_authorization_pointer,
            "artifact_counts_pointer": args.artifact_counts_pointer,
            "artifact_stage_i_sha256_pointer": args.artifact_stage_i_sha256_pointer,
            "artifact_generator_sha256_pointer": args.artifact_generator_sha256_pointer,
            "artifact_stage_i_revision_pointer": args.artifact_stage_i_revision_pointer,
            "artifact_generator_revision_pointer": args.artifact_generator_revision_pointer,
            "artifact_scheduler_evidence_pointer": (
                args.artifact_scheduler_evidence_pointer
            ),
            "artifact_barrier_scheduler_evidence_pointer": (
                args.artifact_barrier_scheduler_evidence_pointer
            ),
            "artifact_source_bundle_sha256_pointer": (
                args.artifact_source_bundle_sha256_pointer
            ),
            "artifact_source_bundle_revisions_pointer": (
                args.artifact_source_bundle_revisions_pointer
            ),
            "bounded_wave_scheduler_evidence": (
                args.bounded_wave_scheduler_evidence_json
            ),
            "stage_i_helper_revision": args.expected_stage_i_revision,
            "generator_revision": args.expected_generator_revision,
        })
        if schema2_evidence_mode(args):
            record.update({
                "recost_recommendations": args.recost_recommendations_json,
                "independent_review_relative_path": args.independent_review_relative_path,
                "independent_review_sha256": args.expected_independent_review_sha256,
            })
        else:
            record["authorized_bounded_wave"] = args.authorized_bounded_wave_json
    record.update({
        "counts": expected_counts(args),
        "forensic_path": str(forensic),
        "utility_sha256": args.expected_utility_sha256,
        "stage_i_helper_sha256": args.expected_stage_i_sha256,
        "generator_relative_path": args.generator_relative_path,
        "generator_sha256": args.expected_generator_sha256,
        "generator_mode": f"{args.expected_generator_mode:04o}",
    })
    if publication_mode(args) == "sole-profile":
        record.update({
            "scheduler_relative_path": args.scheduler_relative_path,
            "scheduler_sha256": args.expected_scheduler_sha256,
            "scheduler_mode": f"{args.expected_scheduler_mode:04o}",
            "source_bundle_mode": f"{args.expected_source_bundle_mode:04o}",
            "source_bundle_verified_revisions": (
                args.expected_source_bundle_verified_revisions_json
            ),
        })
    else:
        record.update({
            "source_bundle_relative_path": args.source_bundle_relative_path,
            "source_bundle_sha256": args.expected_source_bundle_sha256,
            "source_bundle_mode": f"{args.expected_source_bundle_mode:04o}",
            "source_bundle_verified_revisions": (
                args.expected_source_bundle_verified_revisions_json
            ),
            "scheduler_mode": f"{args.expected_scheduler_mode:04o}",
            "recost_request_relative_path": args.recost_request_relative_path,
            "recost_request_sha256": args.expected_request_sha256,
            "generalized_publication_context": selected_v2_publication_context(
                paths, args, allow_expired=allow_expired_v2
            ),
            "forensic_generalized_vectors": "exact-artifact-payload",
        })
    return record


def adoption_bindings(paths: dict[str, Path], args: argparse.Namespace,
                      transaction_id: str, forensic: Path) -> dict[str, object]:
    """Return immutable bindings for one legacy canonical adoption."""

    return {
        **recovery_bindings(paths, args, transaction_id, forensic),
        "operation": "legacy-canonical-adoption",
    }


def validate_adoption_journal(paths: dict[str, Path], args: argparse.Namespace,
                              record: dict[str, object]) -> tuple[str, Path, str]:
    """Validate one resumable legacy canonical adoption journal."""

    transaction_id = record.get("transaction_id")
    if not isinstance(transaction_id, str) or not transaction_id:
        raise ValueError("legacy adoption journal transaction ID is invalid")
    forensic = (
        paths["recost_forensics"]
        / f"{transaction_id}.{paths['canonical'].name}.forensic"
    )
    expected = adoption_bindings(paths, args, transaction_id, forensic)
    state = record.get("state")
    if state not in {"preparing", "forensic-copied", "audit-written"}:
        raise ValueError("legacy adoption journal is not in an allowed recovery state")
    if frozenset(record) != frozenset({*expected, "state", "created_utc"}):
        raise ValueError("legacy adoption journal schema differs")
    for key, value in expected.items():
        if record.get(key) != value:
            raise ValueError(f"legacy adoption journal {key} binding differs")
    require_utc_timestamp(record.get("created_utc"), "legacy adoption journal timestamp")
    return transaction_id, forensic, state


def recovery_staged_identity(record: dict[str, object]) -> tuple[int, int]:
    """Return the exact original staged inode identity bound by a journal."""

    return recovery_identity(
        record,
        "staged_inode_identity",
        "recost recovery journal staged inode",
    )


def validate_recovery_journal(paths: dict[str, Path], args: argparse.Namespace,
                              record: dict[str, object], *,
                              allowed_states: set[str] | None = None,
                              require_forensic_sha256: bool = True,
                              allow_expired_v2: bool = False
                              ) -> tuple[str, Path]:
    """Validate immutable recovery bindings and its forensic copy."""

    transaction_id = record.get("transaction_id")
    if not isinstance(transaction_id, str) or not transaction_id:
        raise ValueError("recost recovery journal transaction ID is invalid")
    forensic = (
        paths["recost_forensics"]
        / f"{transaction_id}.{paths['canonical'].name}.forensic"
    )
    staged_identity = recovery_staged_identity(record)
    expected = recovery_bindings(
        paths,
        args,
        transaction_id,
        forensic,
        allow_expired_v2=allow_expired_v2,
        staged_identity=staged_identity,
    )
    state = record.get("state")
    if allowed_states is None:
        allowed_states = {
            "link-pending",
            "ambiguous-after-link-attempt",
            "single-link-copy-prepared",
        }
    if state not in allowed_states:
        raise ValueError("recost recovery journal is not in an allowed recovery state")
    keys = {*expected, "state", "created_utc"}
    if state == "ambiguous-after-link-attempt":
        keys.add("ambiguity_recorded_utc")
    if state == "single-link-copy-prepared":
        keys.add("canonical_inode_identity")
    if frozenset(record) != frozenset(keys):
        raise ValueError("recost recovery journal schema differs")
    for key, value in expected.items():
        if record.get(key) != value:
            raise ValueError(f"recost recovery journal {key} binding differs")
    require_utc_timestamp(record.get("created_utc"), "recost journal creation timestamp")
    if state == "ambiguous-after-link-attempt":
        require_utc_timestamp(
            record.get("ambiguity_recorded_utc"), "recost journal ambiguity timestamp"
        )
    replacement = require_entry_name(
        str(record["single_link_replacement_name"]),
        "recost recovery journal single-link replacement",
    )
    if state == "single-link-copy-prepared":
        canonical_identity = recovery_identity(
            record,
            "canonical_inode_identity",
            "recost recovery journal canonical inode",
        )
        with bound_parent_descriptor(
            paths["canonical"], "recost single-link recovery"
        ) as directory_descriptor:
            validate_single_link_transition(
                directory_descriptor,
                paths,
                args,
                replacement,
                staged_identity,
                canonical_identity,
                allow_expired_v2=allow_expired_v2,
            )
    else:
        retained_artifacts = [
            path for path in (paths["staged"], paths["canonical"]) if entry_exists(path)
        ]
        if not retained_artifacts:
            raise ValueError("recost recovery journal has no retained artifact inode")
        expected_links = 2 if len(retained_artifacts) == 2 else 1
        for retained_artifact in retained_artifacts:
            require_named_artifact_identity(
                retained_artifact,
                args,
                expected_links=expected_links,
                expected_identity=staged_identity,
                allow_expired_v2=allow_expired_v2,
            )
        # Before the journal binds a canonical-copy inode, this deterministic
        # name is an untrusted crash orphan. Recovery retires it forensically
        # before recreating the copy; strict artifact validation is premature.
    require_managed_directory(paths["recost_forensics_root"], "recost forensic root directory")
    require_exact_directory_entries(
        paths["recost_forensics"],
        "recost forensic directory",
        {forensic},
        managed=True,
    )
    if require_forensic_sha256:
        require_file_sha256(
            forensic,
            args.expected_artifact_sha256,
            "recost forensic copy",
            expected_mode=0o444,
            expected_links=1,
        )
    else:
        with regular_descriptor(
            forensic,
            "recost forensic copy",
            expected_mode=0o444,
            expected_links=1,
        ):
            pass
    return transaction_id, forensic


def require_linked_pair(paths: dict[str, Path], args: argparse.Namespace, *,
                        expected_identity: tuple[int, int] | None = None,
                        allow_expired_v2: bool = False) -> None:
    """Require staged and canonical names to be one authenticated two-link inode."""

    require_artifact_namespace(paths, "linked-pair")
    read_artifact(
        paths["staged"],
        args,
        expected_links=2,
        allow_expired_v2=allow_expired_v2,
    )
    read_artifact(
        paths["canonical"],
        args,
        expected_links=2,
        allow_expired_v2=allow_expired_v2,
    )
    with regular_descriptor(
        paths["staged"], "staged recost artifact", expected_links=2
    ) as staged_descriptor:
        with regular_descriptor(
            paths["canonical"], "canonical recost artifact", expected_links=2
        ) as canonical_descriptor:
            staged = os.fstat(staged_descriptor)
            canonical = os.fstat(canonical_descriptor)
            if (staged.st_dev, staged.st_ino) != (canonical.st_dev, canonical.st_ino):
                raise ValueError(
                    "staged and canonical recost entries are not one linked pair"
                )
            if (
                expected_identity is not None
                and (staged.st_dev, staged.st_ino) != expected_identity
            ):
                raise ValueError("linked recost pair differs from authenticated staged inode")


def require_artifact_namespace_in_directory(directory_descriptor: int,
                                            paths: dict[str, Path],
                                            state: str) -> None:
    """Require one artifact namespace through a stable accounting descriptor."""

    if state == "staged":
        expected = {paths["staged"].name}
    elif state == "linked-pair":
        expected = {paths["staged"].name, paths["canonical"].name}
    else:
        raise ValueError(f"unsupported descriptor-relative namespace state: {state}")
    name = paths["canonical"].name
    retained = {
        entry for entry in os.listdir(directory_descriptor)
        if is_artifact_namespace_name(name, entry)
    }
    if retained != expected:
        raise ValueError(
            f"descriptor-relative recost namespace is not {state}-only; found: "
            + ", ".join(sorted(retained))
        )


def require_named_artifact_identity_in_directory(
    directory_descriptor: int,
    name: str,
    path: Path,
    args: argparse.Namespace,
    *,
    expected_links: int,
    expected_identity: tuple[int, int] | None = None,
    allow_expired_v2: bool = False,
) -> tuple[int, int]:
    """Authenticate one artifact name through a stable accounting descriptor."""

    descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=directory_descriptor)
    try:
        value = os.fstat(descriptor)
        require_regular_profile(
            value,
            path,
            "recost artifact",
            expected_mode=args.expected_artifact_mode,
            expected_links=expected_links,
        )
        retained = read_descriptor_bytes(descriptor)
        if sha256_bytes(retained) != args.expected_artifact_sha256:
            raise ValueError(f"recost artifact checksum has changed: {path}")
        validate_artifact_payload(
            retained, path, args, allow_expired_v2=allow_expired_v2
        )
        identity = descriptor_identity(descriptor)
        if expected_identity is not None and identity != expected_identity:
            raise ValueError(f"recost artifact inode identity has changed: {path}")
        return identity
    finally:
        os.close(descriptor)


def require_linked_pair_in_directory(directory_descriptor: int,
                                     paths: dict[str, Path],
                                     args: argparse.Namespace, *,
                                     expected_identity: tuple[int, int] | None = None,
                                     allow_expired_v2: bool = False) -> None:
    """Require one authenticated linked pair through a stable directory."""

    require_artifact_namespace_in_directory(directory_descriptor, paths, "linked-pair")
    staged_identity = require_named_artifact_identity_in_directory(
        directory_descriptor,
        paths["staged"].name,
        paths["staged"],
        args,
        expected_links=2,
        expected_identity=expected_identity,
        allow_expired_v2=allow_expired_v2,
    )
    canonical_identity = require_named_artifact_identity_in_directory(
        directory_descriptor,
        paths["canonical"].name,
        paths["canonical"],
        args,
        expected_links=2,
        expected_identity=expected_identity,
        allow_expired_v2=allow_expired_v2,
    )
    if staged_identity != canonical_identity:
        raise ValueError("staged and canonical recost entries are not one linked pair")


def artifact_review_target(paths: dict[str, Path],
                           args: argparse.Namespace) -> Path:
    """Return the sole managed accounting target for one schema-2 review."""

    if not schema2_evidence_mode(args) or publication_mode(args) != "v2":
        raise ValueError("artifact-review installation requires schema-2 recost evidence")
    relative = safe_relative_path(
        args.independent_review_relative_path,
        "schema-2 recost evidence independent-review relative path",
    )
    target = paths["root"] / relative
    expected = paths["accounting"] / (
        f"{paths['canonical'].name}.independent_review.json"
    )
    if target != expected:
        raise ValueError(
            "artifact-review installation target must be the exact managed accounting name"
        )
    return target


@contextmanager
def authenticated_artifact_review_candidate(candidate: Path, root: Path,
                                            expected: str):
    """Open one immutable external review candidate without trusting its pathname."""

    if not candidate.is_absolute() or candidate != absolute_path(candidate):
        raise ValueError(
            "artifact-review candidate must be an absolute normalized path"
        )
    root = absolute_path(root)
    if candidate == root or root in candidate.parents:
        raise ValueError("artifact-review candidate must be external to the Stage I root")
    expected = require_sha256(expected, "artifact-review candidate SHA-256")
    with absolute_descriptor(
        candidate, "artifact-review candidate", flags=os.O_RDONLY
    ) as descriptor:
        require_regular_profile(
            os.fstat(descriptor),
            candidate,
            "artifact-review candidate",
            expected_mode=0o444,
            expected_links=1,
            expected_uid=os.geteuid(),
        )
        retained = read_descriptor_bytes(descriptor)
        if sha256_bytes(retained) != expected:
            raise ValueError(f"artifact-review candidate checksum has changed: {candidate}")
        yield candidate, descriptor, retained


def require_artifact_review_candidate_stable(candidate: Path, descriptor: int,
                                             expected: str) -> bytes:
    """Reauthenticate one external candidate descriptor and its selected pathname."""

    expected = require_sha256(expected, "artifact-review candidate SHA-256")
    profile = os.fstat(descriptor)
    require_regular_profile(
        profile,
        candidate,
        "artifact-review candidate",
        expected_mode=0o444,
        expected_links=1,
        expected_uid=os.geteuid(),
    )
    retained = read_descriptor_bytes(descriptor)
    if sha256_bytes(retained) != expected:
        raise ValueError(f"artifact-review candidate checksum has changed: {candidate}")
    with absolute_descriptor(
        candidate, "artifact-review candidate pathname", flags=os.O_RDONLY
    ) as named:
        named_profile = os.fstat(named)
        require_regular_profile(
            named_profile,
            candidate,
            "artifact-review candidate pathname",
            expected_mode=0o444,
            expected_links=1,
            expected_uid=os.geteuid(),
        )
        if (
            descriptor_identity(named) != descriptor_identity(descriptor)
            or read_descriptor_bytes(named) != retained
        ):
            raise ValueError(
                f"artifact-review candidate pathname changed during installation: {candidate}"
            )
    return retained


def read_bound_artifact_review(directory_descriptor: int, target: Path,
                               expected: str, *,
                               expected_identity: tuple[int, int] | None = None
                               ) -> tuple[bytes, tuple[int, int]]:
    """Read one installed review through a stable accounting descriptor."""

    expected = require_sha256(expected, "installed artifact-review SHA-256")
    descriptor = os.open(
        target.name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=directory_descriptor
    )
    try:
        profile = os.fstat(descriptor)
        require_regular_profile(
            profile,
            target,
            "installed artifact review",
            expected_mode=0o444,
            expected_links=1,
            expected_uid=os.geteuid(),
        )
        identity = descriptor_identity(descriptor)
        if expected_identity is not None and identity != expected_identity:
            raise ValueError(f"installed artifact-review inode changed: {target}")
        retained = read_descriptor_bytes(descriptor)
        if sha256_bytes(retained) != expected:
            raise ValueError(f"installed artifact-review checksum has changed: {target}")
    finally:
        os.close(descriptor)
    require_bound_entry_security(
        directory_descriptor,
        target.name,
        profile,
        "installed artifact review",
        expected_sha256=expected,
    )
    return retained, identity


def remove_created_artifact_review(directory_descriptor: int, target: Path,
                                   expected_identity: tuple[int, int]) -> None:
    """Remove only the exact review inode created by a failed transaction."""

    try:
        profile = os.stat(
            target.name, dir_fd=directory_descriptor, follow_symlinks=False
        )
    except FileNotFoundError:
        return
    if (
        not stat.S_ISREG(profile.st_mode)
        or (profile.st_dev, profile.st_ino) != expected_identity
    ):
        raise ValueError(
            f"artifact-review target changed during failed installation: {target}"
        )
    unlink_bound_entry(
        directory_descriptor, target.name, profile, "installed artifact review"
    )


def remove_mode_zero_artifact_review_temporaries(directory_descriptor: int) -> None:
    """Durably discard authenticated legacy and owner-readable review remnants."""

    pattern = re.compile(r"\.cgl-checkpoint-review-[0-9a-f]{32}")
    names = [
        name for name in directory_entry_names(directory_descriptor)
        if pattern.fullmatch(name) is not None
    ]
    for name in names:
        try:
            profile = os.stat(
                name, dir_fd=directory_descriptor, follow_symlinks=False
            )
        except FileNotFoundError:
            continue
        require_regular_mode_subset(
            profile,
            Path(name),
            "artifact-review temporary",
            maximum_mode=0o644,
            expected_links=1,
            expected_uid=os.geteuid(),
        )
        mode = stat.S_IMODE(profile.st_mode)
        if mode not in {0o000, 0o444, 0o600}:
            raise ValueError(
                f"artifact-review temporary mode is {mode:04o}, "
                f"expected 0000, 0444, or 0600: {name}"
            )
        if mode == 0o000:
            if profile.st_size != 0:
                raise ValueError(
                    f"artifact-review mode-0000 temporary is not empty: {name}"
                )
            digest = None
            require_bound_entry_security(
                directory_descriptor,
                name,
                profile,
                "artifact-review mode-0000 temporary",
                allow_mode_zero=True,
            )
        else:
            digest = require_bound_entry_security(
                directory_descriptor,
                name,
                profile,
                "artifact-review temporary",
            )
            assert digest is not None
        require_mutation_authority_bound(directory_descriptor)
        operation_error = None
        try:
            os.unlink(name, dir_fd=directory_descriptor)
        except BaseException as error:
            operation_error = error
        durability_error = fsync_descriptors(directory_descriptor)
        removed = entry_absent(directory_descriptor, name)
        unchanged = entry_security_matches(
            directory_descriptor,
            name,
            profile,
            "artifact-review temporary",
            digest,
            allow_mode_zero=digest is None,
        )
        if not removed and not unchanged:
            raise ValueError(
                f"artifact-review temporary changed during removal: {name}"
            )
        require_mutation_authority_bound(directory_descriptor)
        if durability_error is not None:
            raise durability_error
        if operation_error is not None:
            raise operation_error
        if not removed:
            raise ValueError(
                f"artifact-review temporary removal did not complete: {name}"
            )
    remaining = [
        name for name in directory_entry_names(directory_descriptor)
        if pattern.fullmatch(name) is not None
    ]
    if remaining:
        raise ValueError(
            "artifact-review temporary entries appeared during recovery: "
            + ", ".join(remaining)
        )
    require_mutation_authority_bound(directory_descriptor)


def create_artifact_review(directory_descriptor: int, target: Path,
                           retained: bytes, expected: str
                           ) -> tuple[tuple[int, int], bool]:
    """Publish or authenticate one exact immutable review with atomic no-clobber."""

    expected = require_sha256(expected, "artifact-review candidate SHA-256")
    if sha256_bytes(retained) != expected:
        raise ValueError("artifact-review candidate checksum differs before installation")
    remove_mode_zero_artifact_review_temporaries(directory_descriptor)
    try:
        installed, installed_identity = read_bound_artifact_review(
            directory_descriptor, target, expected
        )
    except FileNotFoundError:
        pass
    except (OSError, ValueError) as error:
        raise ValueError(
            f"artifact-review target already exists or changed: {target}"
        ) from error
    else:
        if installed != retained:
            raise ValueError(
                f"artifact-review target already exists or changed: {target}"
            )
        remove_mode_zero_artifact_review_temporaries(directory_descriptor)
        require_mutation_authority_bound(directory_descriptor)
        return installed_identity, False
    temporary = f".cgl-checkpoint-review-{uuid.uuid4().hex}"
    identity: tuple[int, int] | None = None
    profile = write_bound_exclusive(
        directory_descriptor,
        temporary,
        retained,
        0o444,
        "installed artifact review",
    )
    identity = profile_identity(profile)
    try:
        rename_bound_noreplace(
            directory_descriptor,
            temporary,
            target.name,
            profile,
            "artifact-review publication",
            expected_sha256=expected,
        )
    except BaseException as publication_error:
        require_mutation_authority_bound(directory_descriptor)
        try:
            installed, installed_identity = read_bound_artifact_review(
                directory_descriptor, target, expected
            )
        except (FileNotFoundError, OSError, ValueError):
            installed = None
            installed_identity = None
        if installed == retained and installed_identity is not None:
            created = installed_identity == identity
            try:
                temporary_profile = os.stat(
                    temporary, dir_fd=directory_descriptor, follow_symlinks=False
                )
            except FileNotFoundError:
                pass
            else:
                if profile_identity(temporary_profile) == identity:
                    unlink_bound_entry(
                        directory_descriptor,
                        temporary,
                        temporary_profile,
                        "artifact-review temporary",
                        expected_sha256=expected,
                    )
            remove_mode_zero_artifact_review_temporaries(directory_descriptor)
            os.fsync(directory_descriptor)
            require_mutation_authority_bound(directory_descriptor)
            return installed_identity, created
        try:
            temporary_profile = os.stat(
                temporary, dir_fd=directory_descriptor, follow_symlinks=False
            )
        except FileNotFoundError:
            pass
        else:
            if identity is not None and profile_identity(temporary_profile) == identity:
                unlink_bound_entry(
                    directory_descriptor,
                    temporary,
                    temporary_profile,
                    "artifact-review temporary",
                    expected_sha256=expected,
                )
        raise ValueError(
            f"artifact-review target already exists or changed: {target}"
        ) from publication_error
    profile = os.stat(target.name, dir_fd=directory_descriptor, follow_symlinks=False)
    if profile_identity(profile) != identity:
        raise ValueError(f"artifact-review target changed during installation: {target}")
    remove_mode_zero_artifact_review_temporaries(directory_descriptor)
    os.fsync(directory_descriptor)
    require_mutation_authority_bound(directory_descriptor)
    return identity, True


def read_bound_staged_artifact_for_review_install(
    directory_descriptor: int,
    paths: dict[str, Path],
    args: argparse.Namespace,
    retained_review: tuple[Path, bytes],
    *,
    expected_identity: tuple[int, int] | None = None,
) -> tuple[int, int]:
    """Validate the exact staged artifact against descriptor-backed review bytes."""

    require_artifact_namespace_in_directory(directory_descriptor, paths, "staged")
    descriptor = os.open(
        paths["staged"].name,
        os.O_RDONLY | os.O_NOFOLLOW,
        dir_fd=directory_descriptor,
    )
    try:
        profile = os.fstat(descriptor)
        require_regular_profile(
            profile,
            paths["staged"],
            "staged recost artifact",
            expected_mode=args.expected_artifact_mode,
            expected_links=1,
            expected_uid=os.geteuid(),
        )
        identity = descriptor_identity(descriptor)
        if expected_identity is not None and identity != expected_identity:
            raise ValueError("staged recost artifact changed during review installation")
        retained = read_descriptor_bytes(descriptor)
        if sha256_bytes(retained) != args.expected_artifact_sha256:
            raise ValueError("staged recost artifact checksum changed during review installation")
        validate_artifact_payload(
            retained,
            paths["staged"],
            args,
            retained_review=retained_review,
        )
        return identity
    finally:
        os.close(descriptor)


def install_artifact_review(paths: dict[str, Path], root: Path, offline: bool,
                            args: argparse.Namespace) -> Path:
    """Install one exact external F117 review under canonical accounting."""

    validate_fixture_options(args, offline)
    target = artifact_review_target(paths, args)
    candidate_value = getattr(args, "artifact_review_candidate", None)
    if candidate_value is None:
        raise ValueError("artifact-review installation requires one external candidate")
    with authenticated_artifact_review_candidate(
        Path(candidate_value), root, args.expected_independent_review_sha256
    ) as (candidate, candidate_descriptor, initially_retained):
        with promotion_lock(paths):
            with bound_directory_descriptor(
                paths["accounting"],
                "accounting directory",
            ) as accounting_descriptor:
                retained = require_artifact_review_candidate_stable(
                    candidate,
                    candidate_descriptor,
                    args.expected_independent_review_sha256,
                )
                if retained != initially_retained:
                    raise ValueError("artifact-review candidate changed before installation")
                staged_identity = read_bound_staged_artifact_for_review_install(
                    accounting_descriptor,
                    paths,
                    args,
                    (target, retained),
                )
                created_identity, created = create_artifact_review(
                    accounting_descriptor,
                    target,
                    retained,
                    args.expected_independent_review_sha256,
                )
                try:
                    installed, installed_identity = read_bound_artifact_review(
                        accounting_descriptor,
                        target,
                        args.expected_independent_review_sha256,
                        expected_identity=created_identity,
                    )
                    if installed != retained:
                        raise ValueError(
                            "installed artifact-review bytes differ from the external candidate"
                        )
                    require_artifact_review_candidate_stable(
                        candidate,
                        candidate_descriptor,
                        args.expected_independent_review_sha256,
                    )
                    read_bound_staged_artifact_for_review_install(
                        accounting_descriptor,
                        paths,
                        args,
                        (target, installed),
                        expected_identity=staged_identity,
                    )
                    if installed_identity != created_identity:
                        raise ValueError("installed artifact-review inode identity differs")
                    with absolute_descriptor(
                        paths["accounting"],
                        "accounting directory pathname",
                        flags=os.O_RDONLY | os.O_DIRECTORY,
                    ) as named_accounting:
                        require_trusted_directory_profile(
                            os.fstat(named_accounting),
                            paths["accounting"],
                            "accounting directory pathname",
                        )
                        if descriptor_identity(named_accounting) != descriptor_identity(
                            accounting_descriptor
                        ):
                            raise ValueError(
                                "accounting directory pathname changed during "
                                "artifact-review installation"
                            )
                    os.fsync(accounting_descriptor)
                except BaseException:
                    if created:
                        remove_created_artifact_review(
                            accounting_descriptor, target, created_identity
                        )
                    raise
    return target


def link_failure_is_retry_safe(paths: dict[str, Path], args: argparse.Namespace,
                               expected_identity: tuple[int, int]) -> bool:
    """Return whether a failed link attempt provably left staged-only state."""

    try:
        require_artifact_namespace(paths, "staged")
        require_named_artifact_identity(
            paths["staged"],
            args,
            expected_links=1,
            expected_identity=expected_identity,
        )
    except (OSError, TypeError, ValueError):
        return False
    return True


def retire_unjournaled_single_link_copy(
    paths: dict[str, Path],
    args: argparse.Namespace,
    record: dict[str, object],
    *,
    allow_expired_v2: bool = False,
) -> None:
    """Forensically retire a copy created before its journal update completed."""

    replacement = require_entry_name(
        str(record["single_link_replacement_name"]),
        "recost recovery journal single-link replacement",
    )
    replacement_path = paths["accounting"] / replacement
    if not entry_exists(replacement_path):
        return
    with bound_parent_descriptor(
        paths["canonical"], "recost unjournaled single-link copy"
    ) as directory_descriptor:
        profile = os.stat(
            replacement, dir_fd=directory_descriptor, follow_symlinks=False
        )
        require_regular_mode_subset(
            profile,
            replacement_path,
            "unjournaled canonical recost single-link copy",
            maximum_mode=args.expected_artifact_mode,
            expected_links=1,
            expected_uid=os.geteuid(),
        )
        unlink_bound_entry(
            directory_descriptor,
            replacement,
            profile,
            "unjournaled canonical recost single-link copy",
            expected_sha256=(
                None
                if stat.S_IMODE(profile.st_mode) == 0
                else args.expected_artifact_sha256
            ),
        )


def prepare_journaled_single_link_transition(
    paths: dict[str, Path],
    root: Path,
    offline: bool,
    queue_file: str | None,
    args: argparse.Namespace,
    journal: Path,
    record: dict[str, object],
    *,
    allow_expired_v2: bool = False,
) -> None:
    """Create and durably bind the fresh canonical inode before exchange."""

    retire_unjournaled_single_link_copy(
        paths, args, record, allow_expired_v2=allow_expired_v2
    )
    staged_identity = recovery_staged_identity(record)
    replacement = require_entry_name(
        str(record["single_link_replacement_name"]),
        "recost recovery journal single-link replacement",
    )
    require_empty_queue(root, offline, queue_file)
    with bound_parent_descriptor(
        paths["canonical"], "recost single-link publication"
    ) as directory_descriptor:
        require_trusted_directory_profile(
            os.fstat(directory_descriptor),
            paths["accounting"],
            "accounting directory",
        )
        canonical_identity = prepare_single_link_copy(
            directory_descriptor,
            paths,
            args,
            replacement,
            expected_staged_identity=staged_identity,
            allow_expired_v2=allow_expired_v2,
        )
    record["state"] = "single-link-copy-prepared"
    record.pop("ambiguity_recorded_utc", None)
    record["canonical_inode_identity"] = identity_binding(canonical_identity)
    write_json(journal, record)
    validate_recovery_journal(
        paths,
        args,
        record,
        allow_expired_v2=allow_expired_v2,
    )


def complete_journaled_single_link_transition(
    paths: dict[str, Path],
    root: Path,
    offline: bool,
    queue_file: str | None,
    args: argparse.Namespace,
    record: dict[str, object],
    *,
    allow_expired_v2: bool = False,
) -> None:
    """Complete a durably journaled canonical copy transition."""

    staged_identity = recovery_staged_identity(record)
    canonical_identity = recovery_identity(
        record,
        "canonical_inode_identity",
        "recost recovery journal canonical inode",
    )
    replacement = require_entry_name(
        str(record["single_link_replacement_name"]),
        "recost recovery journal single-link replacement",
    )
    require_empty_queue(root, offline, queue_file)
    with bound_parent_descriptor(
        paths["canonical"], "recost single-link publication"
    ) as directory_descriptor:
        require_trusted_directory_profile(
            os.fstat(directory_descriptor),
            paths["accounting"],
            "accounting directory",
        )
        complete_single_link_transition(
            directory_descriptor,
            paths,
            args,
            replacement,
            staged_identity,
            canonical_identity,
            allow_expired_v2=allow_expired_v2,
            simulate_post_exchange_failure=(
                getattr(args, "simulate_single_link_post_exchange_failure", False)
            ),
        )


def finalize_linked_pair_bound_transactions(
    transaction_descriptor: int,
    paths: dict[str, Path],
    root: Path,
    offline: bool,
    args: argparse.Namespace,
    generator: Path,
    scheduler: Path | tuple[Path, ...],
) -> None:
    """Finalize linked-pair recovery while continuously binding its journal store."""

    remove_forensic_temporaries(
        paths["recost_forensics"],
        "recost forensic temporary",
    )
    journal, record = recost_journal(paths, transaction_descriptor)
    transaction_id, forensic = validate_recovery_journal(
        paths,
        args,
        record,
        allow_expired_v2=True,
    )
    recover_linked_json_publication(
        paths["audit"],
        "recost publication audit recovery",
        expected_mode=publication_audit_mode(args),
    )
    remove_json_temporaries(
        json_temporary_entries(paths["audit"]),
        "recost publication audit temporary",
    )
    retained = set(artifact_namespace_entries(paths))
    if record["state"] == "single-link-copy-prepared":
        complete_journaled_single_link_transition(
            paths,
            root,
            offline,
            args.pre_unlink_squeue_file or args.squeue_file,
            args,
            record,
            allow_expired_v2=True,
        )
    elif retained == {paths["staged"], paths["canonical"]}:
        require_linked_pair(paths, args, allow_expired_v2=True)
        prepare_journaled_single_link_transition(
            paths,
            root,
            offline,
            args.pre_unlink_squeue_file or args.squeue_file,
            args,
            journal,
            record,
            allow_expired_v2=True,
        )
        complete_journaled_single_link_transition(
            paths,
            root,
            offline,
            args.pre_unlink_squeue_file or args.squeue_file,
            args,
            record,
            allow_expired_v2=True,
        )
    elif retained == {paths["canonical"]}:
        read_artifact(paths["canonical"], args, allow_expired_v2=True)
        fsync_directory(paths["accounting"])
    elif retained == {paths["canonical"], paths["audit"]}:
        read_artifact(paths["canonical"], args, allow_expired_v2=True)
        validate_audit(
            paths,
            args,
            initial_source_path(),
            generator,
            scheduler,
            allow_expired_v2=True,
        )
    else:
        raise ValueError(
            "recost recovery namespace is not a linked pair or canonical-only: "
            + ", ".join(str(path) for path in sorted(retained))
        )
    read_artifact(paths["canonical"], args, allow_expired_v2=True)
    if not paths["audit"].exists():
        require_artifact_namespace_after_unlink(paths)
        write_json(
            paths["audit"],
            audit_record(
                paths,
                args,
                initial_source_path(),
                generator,
                scheduler,
                forensic,
                transaction_id,
                allow_expired_v2=True,
            ),
            mode=publication_audit_mode(args),
        )
    validate_audit(
        paths,
        args,
        initial_source_path(),
        generator,
        scheduler,
        allow_expired_v2=True,
    )
    if directory_entry_names(transaction_descriptor) != [journal.name]:
        raise ValueError("recost transaction directory changed before final retirement")
    journal_profile, journal_digest, _, current_record = read_bound_json_object(
        transaction_descriptor,
        journal.name,
        journal,
        "recost transaction journal",
        expected_mode=0o644,
    )
    if current_record != record:
        raise ValueError("recost transaction journal changed before final retirement")
    unlink_bound_entry(
        transaction_descriptor,
        journal.name,
        journal_profile,
        "recost transaction journal",
        expected_sha256=journal_digest,
    )
    if directory_entry_names(transaction_descriptor):
        raise ValueError("recost transaction directory changed after final retirement")
    require_mutation_authority_bound(transaction_descriptor)


def finalize_linked_pair(paths: dict[str, Path], root_dir: Path, root: Path,
                         offline: bool, args: argparse.Namespace) -> None:
    """Finish an authenticated publication left ambiguous after os.link."""

    validate_fixture_options(args, offline)
    with promotion_lock(paths):
        with bound_directory_descriptor(
            paths["recost_transactions"], "recost transaction directory"
        ) as transaction_descriptor:
            require_empty_queue(root, offline, args.squeue_file)
            require_empty_directory(
                paths["stage_i_transactions"], "Stage I transaction directory", trusted=True
            )
            # The durable journal now authorizes recovery; exact historical source
            # bytes remain mandatory, but request freshness and HEAD equality do not.
            authenticate_publication_sources(
                root_dir, args, require_generator_at_head=False
            )
            generator, scheduler = authenticate_external_files(root, args)
            run_reconcile(root_dir, root, offline, args, stage_i_lock_held=True)
            finalize_linked_pair_bound_transactions(
                transaction_descriptor,
                paths,
                root,
                offline,
                args,
                generator,
                scheduler,
            )
    verify_promoted_state(paths, root_dir, root, offline, args)


def retire_preparing_bound_transactions(transaction_descriptor: int,
                                        paths: dict[str, Path], root_dir: Path,
                                        root: Path, offline: bool,
                                        args: argparse.Namespace) -> bool:
    """Retire one prepublication journal while its store remains descriptor-bound."""

    try:
        journal, record = recost_journal(paths, transaction_descriptor)
    except NoRecostJournal:
        if directory_entry_names(transaction_descriptor):
            raise ValueError("recost transaction directory changed while proving empty")
        require_empty_directory(
            paths["recost_forensics"], "recost forensic directory", managed=True
        )
        require_empty_queue(root, offline, args.squeue_file)
        verify_retirable_staged_state(paths, root_dir, root, offline, args)
        require_mutation_authority_bound(transaction_descriptor)
        if directory_entry_names(transaction_descriptor):
            raise ValueError("recost transaction directory changed while proving empty")
        return True
    transaction_id = record.get("transaction_id")
    if not isinstance(transaction_id, str) or not transaction_id:
        raise ValueError("recost recovery journal transaction ID is invalid")
    forensic = (
        paths["recost_forensics"]
        / f"{transaction_id}.{paths['canonical'].name}.forensic"
    )
    staged_identity = recovery_staged_identity(record)
    expected = recovery_bindings(
        paths,
        args,
        transaction_id,
        forensic,
        allow_expired_v2=True,
        staged_identity=staged_identity,
    )
    if frozenset(record) != frozenset({*expected, "state", "created_utc"}):
        raise ValueError("recost prepublication journal schema differs")
    for key, value in expected.items():
        if record.get(key) != value:
            raise ValueError(f"recost prepublication journal {key} binding differs")
    state = record.get("state")
    if state not in {"preparing", "link-pending"}:
        raise ValueError("recost transaction is not an interrupted prepublication state")
    require_utc_timestamp(record.get("created_utc"), "recost journal creation timestamp")
    require_named_artifact_identity(
        paths["staged"],
        args,
        expected_links=1,
        expected_identity=staged_identity,
        allow_expired_v2=True,
    )
    retire_unjournaled_single_link_copy(
        paths,
        args,
        record,
        allow_expired_v2=True,
    )
    retained = set(directory_entries(paths["recost_forensics"]))
    if retained not in (set(), {forensic}):
        raise ValueError(
            "recost forensic directory entries differ; found: "
            + ", ".join(str(item) for item in sorted(retained))
        )
    if state == "link-pending":
        if entry_exists(forensic):
            require_file_sha256(
                forensic,
                args.expected_artifact_sha256,
                "recost forensic copy",
                expected_mode=0o444,
                expected_links=1,
            )
    elif entry_exists(forensic):
        with regular_descriptor(
            forensic,
            "recost forensic copy",
            expected_links=1,
        ) as descriptor:
            require_regular_mode_subset(
                os.fstat(descriptor),
                forensic,
                "recost forensic copy",
                maximum_mode=0o444,
                expected_links=1,
                expected_uid=os.geteuid(),
            )
    require_empty_queue(root, offline, args.squeue_file)
    if entry_exists(forensic):
        unlink_durable(forensic)
    if directory_entry_names(transaction_descriptor) != [journal.name]:
        raise ValueError("recost transaction directory changed before journal retirement")
    journal_profile, journal_digest, _, current_record = read_bound_json_object(
        transaction_descriptor,
        journal.name,
        journal,
        "recost transaction journal",
        expected_mode=0o644,
    )
    if current_record != record:
        raise ValueError("recost transaction journal changed before retirement")
    unlink_bound_entry(
        transaction_descriptor,
        journal.name,
        journal_profile,
        "recost transaction journal",
        expected_sha256=journal_digest,
    )
    if directory_entry_names(transaction_descriptor):
        raise ValueError("recost transaction directory changed after journal retirement")
    require_mutation_authority_bound(transaction_descriptor)
    return False


def retire_preparing(paths: dict[str, Path], root_dir: Path, root: Path,
                     offline: bool, args: argparse.Namespace) -> None:
    """Retire a staged-only transaction interrupted before publication."""

    validate_fixture_options(args, offline)
    with promotion_lock(paths):
        with bound_directory_descriptor(
            paths["recost_transactions"], "recost transaction directory"
        ) as transaction_descriptor:
            require_empty_queue(root, offline, args.squeue_file)
            require_empty_directory(
                paths["stage_i_transactions"], "Stage I transaction directory", trusted=True
            )
            authenticate_publication_sources(
                root_dir, args, require_generator_at_head=False
            )
            authenticate_external_files(root, args)
            run_reconcile(root_dir, root, offline, args, stage_i_lock_held=True)
            require_artifact_namespace(paths, "staged")
            read_artifact(paths["staged"], args, allow_expired_v2=True)
            require_managed_directory(
                paths["recost_forensics_root"], "recost forensic root directory"
            )
            require_managed_directory(
                paths["recost_forensics"], "recost forensic directory"
            )
            remove_forensic_temporaries(
                paths["recost_forensics"],
                "recost forensic temporary",
            )
            empty = retire_preparing_bound_transactions(
                transaction_descriptor, paths, root_dir, root, offline, args
            )
            if empty:
                return
    verify_retirable_staged_state(paths, root_dir, root, offline, args)


def promote_bound_transactions(transaction_descriptor: int, paths: dict[str, Path],
                               root: Path, offline: bool,
                               args: argparse.Namespace) -> None:
    """Publish one staged artifact while continuously binding its journal store."""

    if directory_entry_names(transaction_descriptor):
        raise ValueError("recost transaction directory is not empty")
    if not paths["recost_forensics_root"].exists():
        mkdir_durable(paths["recost_forensics_root"])
    require_managed_directory(paths["recost_forensics_root"], "recost forensic root directory")
    if not paths["recost_forensics"].exists():
        mkdir_durable(paths["recost_forensics"])
    remove_forensic_temporaries(
        paths["recost_forensics"],
        "recost forensic temporary",
    )
    require_empty_directory(
        paths["recost_forensics"], "recost forensic directory", managed=True
    )
    transaction_id = f"{utc_now().replace(':', '')}-{uuid.uuid4().hex}"
    journal = paths["recost_transactions"] / f"{transaction_id}.json"
    forensic = (
        paths["recost_forensics"]
        / f"{transaction_id}.{paths['canonical'].name}.forensic"
    )
    with regular_descriptor(
        paths["staged"],
        "staged recost artifact",
        expected_mode=args.expected_artifact_mode,
        expected_links=1,
    ) as staged_descriptor:
        if sha256_descriptor(staged_descriptor) != args.expected_artifact_sha256:
            raise ValueError("staged recost artifact checksum changed before journal")
        staged_identity = descriptor_identity(staged_descriptor)
    record = {
        **recovery_bindings(
            paths,
            args,
            transaction_id,
            forensic,
            staged_identity=staged_identity,
        ),
        "state": "preparing",
        "created_utc": utc_now(),
    }
    write_json(journal, record)
    try:
        copy_forensic(
            paths["staged"],
            forensic,
            args.expected_artifact_sha256,
            expected_mode=args.expected_artifact_mode,
        )
        record["state"] = "link-pending"
        write_json(journal, record)
        queue_file = args.pre_link_squeue_file or args.squeue_file
        with regular_descriptor(
            paths["staged"],
            "staged recost artifact",
            expected_mode=args.expected_artifact_mode,
            expected_links=1,
        ) as staged_descriptor:
            if sha256_descriptor(staged_descriptor) != args.expected_artifact_sha256:
                raise ValueError("staged recost artifact checksum changed before link")
            if descriptor_identity(staged_descriptor) != staged_identity:
                raise ValueError("staged recost artifact inode changed before link")
            with bound_parent_descriptor(
                paths["canonical"], "recost publication accounting"
            ) as directory_descriptor:
                link_after_empty_queue(
                    root,
                    offline,
                    queue_file,
                    directory_descriptor,
                    staged_descriptor,
                    paths,
                    args,
                    staged_identity,
                )
                require_linked_pair_in_directory(
                    directory_descriptor,
                    paths,
                    args,
                    expected_identity=staged_identity,
                )
                os.fsync(directory_descriptor)
                if args.simulate_post_link_failure:
                    raise ValueError("simulated post-link publication failure")
            prepare_journaled_single_link_transition(
                paths,
                root,
                offline,
                queue_file,
                args,
                journal,
                record,
            )
            complete_journaled_single_link_transition(
                paths,
                root,
                offline,
                queue_file,
                args,
                record,
            )
            canonical_identity = recovery_identity(
                record,
                "canonical_inode_identity",
                "recost recovery journal canonical inode",
            )
            require_named_artifact_identity(
                paths["canonical"],
                args,
                expected_links=1,
                expected_identity=canonical_identity,
            )
        require_artifact_namespace_after_unlink(paths)
        generator, scheduler = authenticate_external_files(root, args)
        write_json(
            paths["audit"],
            audit_record(
                paths,
                args,
                initial_source_path(),
                generator,
                scheduler,
                forensic,
                transaction_id,
            ),
            mode=publication_audit_mode(args),
        )
        validate_audit(
            paths,
            args,
            initial_source_path(),
            generator,
            scheduler,
        )
        remove_bound_recost_journal(
            transaction_descriptor, journal, record, "recost transaction journal"
        )
    except BaseException:
        if link_failure_is_retry_safe(paths, args, staged_identity):
            remove_prelink_records_bound(
                transaction_descriptor, journal, forensic, record
            )
        elif record.get("state") == "link-pending":
            record["state"] = "ambiguous-after-link-attempt"
            record["ambiguity_recorded_utc"] = utc_now()
            write_json(journal, record)
        elif record.get("state") != "single-link-copy-prepared":
            remove_prelink_records_bound(
                transaction_descriptor, journal, forensic, record
            )
        raise


def promote(paths: dict[str, Path], root_dir: Path, root: Path, offline: bool,
            args: argparse.Namespace) -> None:
    """Publish one staged recost artifact exactly once."""

    validate_fixture_options(args, offline)
    with promotion_lock(paths):
        transaction_boundary = (
            directory_boundary(paths["recost_transactions"].lstat())
            if entry_exists(paths["recost_transactions"])
            else None
        )
        mkdir_durable(paths["recost_transactions"])
        with bound_directory_descriptor(
            paths["recost_transactions"], "recost transaction directory"
        ) as transaction_descriptor:
            if (
                transaction_boundary is not None
                and directory_boundary(os.fstat(transaction_descriptor))
                != transaction_boundary
            ):
                raise ValueError("recost transaction directory changed before binding")
            verify_staged_state(
                paths, root_dir, root, offline, args, stage_i_lock_held=True
            )
            reexecute_v2_recost(
                paths,
                root_dir,
                root,
                args,
                paths["staged"],
                stage_i_lock_held=True,
            )
            promote_bound_transactions(
                transaction_descriptor, paths, root, offline, args
            )
    verify_promoted_state(paths, root_dir, root, offline, args)


def adopt_legacy_bound_transactions(
    transaction_descriptor: int,
    paths: dict[str, Path],
    root: Path,
    offline: bool,
    args: argparse.Namespace,
    generator: Path,
    scheduler: Path | tuple[Path, ...],
    *,
    transactions_existed: bool,
) -> None:
    """Adopt one canonical artifact while continuously binding its journal store."""

    initial_names = directory_entry_names(transaction_descriptor)
    if not transactions_existed and initial_names:
        raise ValueError(
            "legacy adoption transaction entries appeared before directory binding"
        )
    if transactions_existed and not initial_names:
        raise ValueError(
            "legacy adoption requires an absent transaction directory "
            "or one resumable journal"
        )
    temporary_pattern = re.compile(r"\..+\.json\.[0-9]+\.[0-9a-f]{32}\.tmp")
    recovering_prejournal_temporary = bool(initial_names) and all(
        temporary_pattern.fullmatch(name) is not None for name in initial_names
    )
    recover_linked_json_publication(
        paths["audit"],
        "legacy adoption publication audit recovery",
        expected_mode=publication_audit_mode(args),
    )
    remove_forensic_temporaries(
        paths["recost_forensics"],
        "recost forensic temporary",
    )
    remove_json_temporaries(
        json_temporary_entries(paths["audit"]),
        "recost publication audit temporary",
    )
    try:
        journal, record = recost_journal(paths, transaction_descriptor)
    except NoRecostJournal:
        if transactions_existed and not recovering_prejournal_temporary:
            raise ValueError("legacy adoption transaction journal disappeared")
        require_artifact_namespace(paths, "canonical-pending-audit")
        require_empty_directory(
            paths["recost_forensics"], "recost forensic directory", managed=True
        )
        transaction_id = f"{utc_now().replace(':', '')}-{uuid.uuid4().hex}"
        journal = paths["recost_transactions"] / f"{transaction_id}.json"
        forensic = (
            paths["recost_forensics"]
            / f"{transaction_id}.{paths['canonical'].name}.forensic"
        )
        record = {
            **adoption_bindings(paths, args, transaction_id, forensic),
            "state": "preparing",
            "created_utc": utc_now(),
        }
        write_json(journal, record)
        state = "preparing"
        if args.simulate_adoption_interruption_after_journal:
            raise ValueError("simulated interruption after legacy adoption journal")
    else:
        transaction_id, forensic, state = validate_adoption_journal(
            paths, args, record
        )

    if state == "preparing":
        require_artifact_namespace(paths, "canonical-pending-audit")
        retained = set(directory_entries(paths["recost_forensics"]))
        if retained not in (set(), {forensic}):
            raise ValueError(
                "recost forensic directory entries differ; found: "
                + ", ".join(str(item) for item in sorted(retained))
            )
        recover_or_copy_forensic(
            paths["canonical"],
            forensic,
            args.expected_artifact_sha256,
            expected_mode=args.expected_artifact_mode,
            simulate_interruption_before_directory_fsync=(
                args.simulate_adoption_interruption_before_forensic_directory_fsync
            ),
        )
        record["state"] = "forensic-copied"
        write_json(journal, record)
        state = "forensic-copied"
        if args.simulate_adoption_interruption_after_forensic:
            raise ValueError("simulated interruption after legacy adoption forensic copy")

    if state == "forensic-copied":
        require_exact_directory_entries(
            paths["recost_forensics"],
            "recost forensic directory",
            {forensic},
            managed=True,
        )
        require_file_sha256(
            forensic,
            args.expected_artifact_sha256,
            "recost forensic copy",
            expected_mode=0o444,
            expected_links=1,
        )
        if not entry_exists(paths["audit"]):
            require_artifact_namespace(paths, "canonical-pending-audit")
            read_artifact(paths["canonical"], args)
            require_empty_queue(root, offline, args.squeue_file)
            write_json(
                paths["audit"],
                adoption_audit_record(
                    paths,
                    args,
                    initial_source_path(),
                    generator,
                    scheduler,
                    forensic,
                    transaction_id,
                ),
                mode=publication_audit_mode(args),
                simulate_interruption_before_directory_fsync=(
                    args.simulate_adoption_interruption_before_audit_directory_fsync
                ),
            )
            if args.simulate_adoption_interruption_after_audit:
                raise ValueError("simulated interruption after legacy adoption audit")
        fsync_directory(paths["accounting"])
        require_artifact_namespace(paths, "promoted")
        validate_audit(paths, args, initial_source_path(), generator, scheduler)
        record["state"] = "audit-written"
        write_json(journal, record)
        state = "audit-written"

    if state == "audit-written":
        require_artifact_namespace(paths, "promoted")
        validate_audit(paths, args, initial_source_path(), generator, scheduler)
    require_empty_queue(root, offline, args.squeue_file)
    if directory_entry_names(transaction_descriptor) != [journal.name]:
        raise ValueError("legacy adoption transaction directory changed before retirement")
    journal_profile, journal_digest, _, current_record = read_bound_json_object(
        transaction_descriptor,
        journal.name,
        journal,
        "recost transaction journal",
        expected_mode=0o644,
    )
    if current_record != record:
        raise ValueError("legacy adoption transaction journal changed before retirement")
    unlink_bound_entry(
        transaction_descriptor,
        journal.name,
        journal_profile,
        "recost transaction journal",
        expected_sha256=journal_digest,
    )
    if directory_entry_names(transaction_descriptor):
        raise ValueError("legacy adoption transaction directory changed after retirement")
    require_mutation_authority_bound(transaction_descriptor)


def adopt_legacy_canonical(paths: dict[str, Path], root_dir: Path, root: Path,
                           offline: bool, args: argparse.Namespace) -> None:
    """Attest an authenticated canonical-only artifact without publication claims."""

    validate_fixture_options(args, offline)
    with promotion_lock(paths):
        transactions_existed = entry_exists(paths["recost_transactions"])
        transactions_boundary = None
        if transactions_existed:
            transactions_boundary = directory_boundary(
                paths["recost_transactions"].lstat()
            )
        mkdir_durable(paths["recost_transactions"])
        with bound_directory_descriptor(
            paths["recost_transactions"], "recost transaction directory"
        ) as transaction_descriptor:
            if (
                transactions_boundary is not None
                and directory_boundary(os.fstat(transaction_descriptor))
                != transactions_boundary
            ):
                raise ValueError(
                    "legacy adoption transaction directory changed before binding"
                )
            require_empty_queue(root, offline, args.squeue_file)
            require_empty_directory(
                paths["stage_i_transactions"], "Stage I transaction directory", trusted=True
            )
            authenticate_publication_sources(root_dir, args)
            generator, scheduler = authenticate_external_files(root, args)
            read_artifact(paths["canonical"], args)
            reexecute_v2_recost(
                paths,
                root_dir,
                root,
                args,
                paths["canonical"],
                stage_i_lock_held=True,
            )
            run_reconcile(root_dir, root, offline, args, stage_i_lock_held=True)
            if not transactions_existed:
                require_artifact_namespace(paths, "canonical-pending-audit")
                if entry_exists(paths["recost_forensics_root"]):
                    raise ValueError(
                        "legacy adoption requires an absent recost forensic root "
                        "before initial attestation"
                    )
            if not paths["recost_forensics_root"].exists():
                mkdir_durable(paths["recost_forensics_root"])
            require_managed_directory(
                paths["recost_forensics_root"], "recost forensic root directory"
            )
            if not paths["recost_forensics"].exists():
                mkdir_durable(paths["recost_forensics"])
            require_managed_directory(
                paths["recost_forensics"], "recost forensic directory"
            )
            adopt_legacy_bound_transactions(
                transaction_descriptor,
                paths,
                root,
                offline,
                args,
                generator,
                scheduler,
                transactions_existed=transactions_existed,
            )
    verify_promoted_state(paths, root_dir, root, offline, args)


def require_artifact_namespace_after_unlink(paths: dict[str, Path]) -> None:
    """Require only the canonical artifact before its audit is written."""

    retained = set(artifact_namespace_entries(paths))
    expected = {paths["canonical"]}
    if retained != expected:
        raise ValueError(
            "recost namespace is not canonical-only after publication; found: "
            + ", ".join(str(path) for path in sorted(retained))
        )


def add_artifact_arguments(command: argparse.ArgumentParser) -> None:
    """Add explicit generic recost binding arguments."""

    command.add_argument("--artifact-name", required=True)
    command.add_argument("--expected-artifact-sha256", type=sha256_arg, required=True)
    command.add_argument("--expected-artifact-mode", type=artifact_mode_arg, default=0o644)
    command.add_argument("--generator-relative-path", required=True)
    command.add_argument("--expected-generator-sha256", type=sha256_arg, required=True)
    command.add_argument("--expected-generator-mode", type=generator_mode_arg, default=0o755)
    command.add_argument("--scheduler-relative-path")
    command.add_argument("--expected-scheduler-sha256", type=sha256_arg)
    command.add_argument("--expected-scheduler-mode", type=scheduler_mode_arg, default=0o644)
    authorization = command.add_mutually_exclusive_group(required=True)
    authorization.add_argument(
        "--authorized-next-segment-profile-json",
        type=json_object_arg,
    )
    authorization.add_argument(
        "--authorized-v2-json",
        "--authorized-bounded-wave-json",
        dest="authorized_bounded_wave_json",
        type=json_object_arg,
    )
    authorization.add_argument(
        "--recost-recommendations-json",
        dest="recost_recommendations_json",
        type=json_object_arg,
        help="Exact schema-2 non-authorizing recommendation packet.",
    )
    command.add_argument(
        "--v2-scheduler-evidence-json",
        "--bounded-wave-scheduler-evidence-json",
        dest="bounded_wave_scheduler_evidence_json",
        type=json_array_arg,
    )
    command.add_argument("--recost-request-relative-path")
    command.add_argument("--expected-request-sha256", type=sha256_arg)
    command.add_argument("--independent-review-relative-path")
    command.add_argument("--expected-independent-review-sha256", type=sha256_arg)
    command.add_argument("--source-bundle-relative-path")
    command.add_argument("--expected-source-bundle-sha256", type=sha256_arg)
    command.add_argument("--expected-stage-i-revision", type=revision_arg)
    command.add_argument("--expected-generator-revision", type=revision_arg)
    command.add_argument(
        "--expected-source-bundle-verified-revisions-json",
        type=revision_array_arg,
    )
    command.add_argument(
        "--expected-source-bundle-mode",
        type=source_bundle_mode_arg,
        default=0o644,
    )
    for key in COUNT_KEYS:
        command.add_argument(
            f"--expected-{key.replace('_', '-')}",
            dest=f"expected_{key}",
            type=nonnegative_int_arg,
            required=True,
        )
    command.add_argument("--artifact-epoch-pointer", default="/execution_epoch")
    command.add_argument(
        "--artifact-profile-pointer",
        default="/authorization/sole_next_segment_profile",
    )
    command.add_argument(
        "--artifact-authorization-pointer",
        default="/authorization",
    )
    command.add_argument("--artifact-counts-pointer", default="/reconcile/counts")
    command.add_argument(
        "--artifact-stage-i-sha256-pointer",
        default="/provenance/stage_i_helper_sha256",
    )
    command.add_argument(
        "--artifact-generator-sha256-pointer",
        default="/provenance/generator_sha256",
    )
    command.add_argument(
        "--artifact-stage-i-revision-pointer",
        default="/provenance/stage_i_helper_revision",
    )
    command.add_argument(
        "--artifact-generator-revision-pointer",
        default="/provenance/generator_revision",
    )
    command.add_argument(
        "--artifact-scheduler-sha256-pointer",
        default="/provenance/scheduler_sha256",
    )
    command.add_argument(
        "--artifact-scheduler-evidence-pointer",
        default="/provenance/scheduler_evidence",
    )
    command.add_argument(
        "--artifact-barrier-scheduler-evidence-pointer",
        default="/barrier/scheduler_evidence",
    )
    command.add_argument(
        "--artifact-source-bundle-sha256-pointer",
        default="/provenance/source_bundle_sha256",
    )
    command.add_argument(
        "--artifact-source-bundle-revisions-pointer",
        default="/provenance/source_bundle_verified_revisions",
    )


def parser() -> argparse.ArgumentParser:
    """Build the retained companion CLI."""

    command = argparse.ArgumentParser(prog=UTILITY_RELATIVE.name, description=__doc__)
    command.add_argument("--root", default=str(DEFAULT_ROOT))
    command.add_argument("--allow-local-root", action="store_true")
    command.add_argument("--expected-utility-sha256", type=sha256_arg, required=True)
    command.add_argument("--expected-stage-i-sha256", type=sha256_arg, required=True)
    command.add_argument(
        "--squeue-file",
        help="Offline local-root fixture input in place of querying squeue.",
    )
    actions = command.add_subparsers(dest="action", required=True)
    for name in (
        "audit-recost",
        "verify-staged-recost",
        "install-artifact-review",
        "promote-recost",
        "verify-promoted-recost",
        "finalize-linked-pair",
        "retire-preparing",
        "adopt-legacy-canonical",
    ):
        action = actions.add_parser(name)
        add_artifact_arguments(action)
        if name == "install-artifact-review":
            action.add_argument(
                "--artifact-review-candidate",
                type=Path,
                required=True,
                help="External immutable review candidate copied into canonical accounting.",
            )
        if name == "audit-recost":
            action.add_argument(
                "--publication-state",
                choices=("staged", "promoted"),
                required=True,
            )
        if name == "promote-recost":
            action.add_argument(
                "--pre-link-squeue-file",
                help="Offline queue fixture read immediately before publication.",
            )
            action.add_argument(
                "--simulate-post-link-failure",
                action="store_true",
                help="Offline fixture hook retaining an ambiguous link attempt.",
            )
            action.add_argument(
                "--simulate-link-return-failure",
                action="store_true",
                help="Offline fixture hook interrupting immediately after os.link().",
            )
            action.add_argument(
                "--simulate-single-link-post-exchange-failure",
                action="store_true",
                help="Offline fixture hook interrupting after canonical-copy exchange.",
            )
        if name == "finalize-linked-pair":
            action.add_argument(
                "--pre-unlink-squeue-file",
                help="Offline queue fixture read immediately before staged unlink.",
            )
        if name == "adopt-legacy-canonical":
            action.add_argument(
                "--simulate-adoption-interruption-after-journal",
                action="store_true",
                help="Offline fixture hook interrupting after the adoption journal.",
            )
            action.add_argument(
                "--simulate-adoption-interruption-before-forensic-directory-fsync",
                action="store_true",
                help="Offline fixture hook interrupting before forensic-directory fsync.",
            )
            action.add_argument(
                "--simulate-adoption-interruption-after-forensic",
                action="store_true",
                help="Offline fixture hook interrupting after the adoption forensic copy.",
            )
            action.add_argument(
                "--simulate-adoption-interruption-before-audit-directory-fsync",
                action="store_true",
                help="Offline fixture hook interrupting before audit-directory fsync.",
            )
            action.add_argument(
                "--simulate-adoption-interruption-after-audit",
                action="store_true",
                help="Offline fixture hook interrupting after the adoption audit.",
            )
    return command


def main() -> int:
    """Run one authenticated retained-companion operation."""

    try:
        utility_sha = expected_self_sha256(sys.argv[1:])
        _, root_dir = authenticate_self(utility_sha)
        args = parser().parse_args()
        if args.action == "adopt-legacy-canonical" and schema2_evidence_mode(args):
            raise ValueError(
                "schema-2 recost evidence requires an observed normal publication transition"
            )
        root, offline = require_root(Path(args.root), args.allow_local_root)
        require_trusted_directory(root, "Stage I root")
        if not offline:
            require_canonical_repository(root_dir, offline)
            authenticate_committed_utility(root_dir, args.expected_utility_sha256)
        paths = layout(root, args.artifact_name)
        if args.action == "audit-recost":
            if args.publication_state == "staged":
                verify_staged_state(paths, root_dir, root, offline, args)
                print(f"Audited staged-only recost artifact: {paths['staged']}")
            else:
                verify_promoted_state(paths, root_dir, root, offline, args)
                print(f"Audited canonical-only recost artifact: {paths['canonical']}")
            return 0
        if args.action == "verify-staged-recost":
            verify_staged_state(paths, root_dir, root, offline, args)
            print(f"Verified staged-only recost artifact: {paths['staged']}")
            return 0
        if args.action == "install-artifact-review":
            target = install_artifact_review(paths, root, offline, args)
            print(f"Installed authenticated artifact review: {target}")
            return 0
        if args.action == "promote-recost":
            promote(paths, root_dir, root, offline, args)
            print(f"Promoted and verified recost artifact: {paths['canonical']}")
            return 0
        if args.action == "verify-promoted-recost":
            verify_promoted_state(paths, root_dir, root, offline, args)
            print(f"Verified canonical-only recost artifact: {paths['canonical']}")
            return 0
        if args.action == "finalize-linked-pair":
            finalize_linked_pair(paths, root_dir, root, offline, args)
            print(f"Finalized and verified recost artifact: {paths['canonical']}")
            return 0
        if args.action == "retire-preparing":
            retire_preparing(paths, root_dir, root, offline, args)
            print(f"Retired interrupted preparing transaction: {paths['staged']}")
            return 0
        if args.action == "adopt-legacy-canonical":
            adopt_legacy_canonical(paths, root_dir, root, offline, args)
            print(f"Adopted and verified legacy canonical artifact: {paths['canonical']}")
            return 0
        raise ValueError(f"unsupported action: {args.action}")
    except (KeyError, OSError, TypeError, ValueError,
            subprocess.CalledProcessError) as error:
        print(f"Stage I recost checkpoint utility failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
