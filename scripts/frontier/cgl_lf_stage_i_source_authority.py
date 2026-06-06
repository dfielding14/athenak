#!/opt/cray/pe/python/3.11.7/bin/python3.11 -I
"""Draft and promote the exact F-118 Stage I current-source authority.

Authenticated drafting emits deterministic evidence and audit candidates but
never creates reviews, self-approves, or publishes.  Promotion changes source
selection only.  It never authorizes prepare, submit, scheduler mutation,
scientific changes, or historical manifest rebinding.  Publication is a
recoverable forward transaction under the canonical Stage-I lock; the
single-link exact-digest 0444 F-118 publication audit is the sole authority
commit marker.
Before that marker, the exact complete transaction matching the requested audit
is a recovery input; only complete inert staged transactions with no publication
state may coexist as non-authoritative prior-attempt debris.  Incomplete,
non-inert, or unsafe transactions block.  All retained transactions and bounded
private-publication remnants are non-authoritative recovery debris after the
marker.  Normal forward completion removes authenticated incomplete private
remnants before committing; private names observed after commit are never
mutated.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack, contextmanager
import ctypes
from datetime import datetime, timezone
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
import sys
import tempfile
import uuid
from typing import Iterator


DEFAULT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
CANONICAL_REPOSITORY_ROOT = Path("/autofs/nccs-svm1_home2/dfielding/athenak-df")
CANONICAL_TRUSTED_PROJECT_BOUNDARY = Path("/lustre/orion/ast207/proj-shared")
CANONICAL_TRUSTED_PROJECT_GID = 31114
CANONICAL_OWNER_UID = 18664
CANONICAL_PUBLIC_NAMESPACE_PROFILES = (
    (CANONICAL_TRUSTED_PROJECT_BOUNDARY, 0o2770, 0, CANONICAL_TRUSTED_PROJECT_GID),
    (CANONICAL_TRUSTED_PROJECT_BOUNDARY / "dfielding", 0o2755, CANONICAL_OWNER_UID,
     CANONICAL_TRUSTED_PROJECT_GID),
    (DEFAULT_ROOT, 0o2755, CANONICAL_OWNER_UID, CANONICAL_TRUSTED_PROJECT_GID),
)
EXECUTION_EPOCH = "E03-forcing-policy"
EXECUTION_EPOCH_SLUG = "E03_forcing_policy"
CHECKPOINT = "F-118"
PUBLISHER_RELATIVE = PurePosixPath(
    "scripts/frontier/cgl_lf_stage_i_source_authority.py"
)
F118_NAME = (
    "mks24_stage_i_E03_forcing_policy_"
    "F118_current_source_authority_supersession_evidence.json"
)
F116_NAME = (
    "mks24_stage_i_E03_forcing_policy_"
    "F116_current_source_authority_supersession_evidence.json"
)
F115_NAME = (
    "mks24_stage_i_E03_forcing_policy_"
    "F115_source_bundle_recovery_supersession_evidence.json"
)
F115_PATHS = {
    "evidence": PurePosixPath("accounting") / F115_NAME,
    "publication_audit": PurePosixPath("accounting") / f"{F115_NAME}.publication_audit.json",
    "provenance_review": (
        PurePosixPath("accounting") / f"{F115_NAME}.provenance_security_review.json"
    ),
    "plasma_review": (
        PurePosixPath("accounting") / f"{F115_NAME}.plasma_scientific_review.json"
    ),
}
F116_PATHS = {
    "evidence": PurePosixPath("accounting") / F116_NAME,
    "provenance_review": (
        PurePosixPath("accounting") / f"{F116_NAME}.provenance_security_review.json"
    ),
    "plasma_review": (
        PurePosixPath("accounting") / f"{F116_NAME}.plasma_scientific_review.json"
    ),
    "publication_audit": (
        PurePosixPath("accounting") / f"{F116_NAME}.publication_audit.json"
    ),
}
F118_PATHS = {
    "evidence": PurePosixPath("accounting") / F118_NAME,
    "provenance_review": (
        PurePosixPath("accounting") / f"{F118_NAME}.provenance_security_review.json"
    ),
    "plasma_review": (
        PurePosixPath("accounting") / f"{F118_NAME}.plasma_scientific_review.json"
    ),
    "publication_audit": (
        PurePosixPath("accounting") / f"{F118_NAME}.publication_audit.json"
    ),
}
F115_CANONICAL_SHA256 = {
    "evidence": "cb50beb064678a9446ac33801a0023547d06c432d8c59b6bd3fbe34b11cf0391",
    "publication_audit": "5923e3872b1d4a84a147bcd1d81fddc20ee79b72bfc410683d083039781f5b1f",
    "provenance_review": "6fcd19f9267f36332742f1f103968098216bd8e6b42fa9821964ab8df704bacd",
    "plasma_review": "a78357ed90e593809b1a82a641d2b40d651b16569ec943fc1781940e569b8440",
}
F116_CANONICAL_SHA256 = {
    "evidence": "6cdbf9e4d10f1282744c6274aa3ef08afec4c510420296837fdbbdfcefe30a2a",
    "provenance_review": "e9731aab8305505e058c68ae8bb61c9ec5ff4885bde1bee3162719c41ab9bafd",
    "plasma_review": "bc4da6897263843f5d233f996539ff047d6ef465b775a9b5cdc8f864de96a6b8",
    "publication_audit": "3a6168e3039c02656b38ebdfcadffc07a2f151b1a474084307a80a341ba83096",
}
BRIDGE_REVISION = "36140ea825cb853b298714c27720440fdab60b9e"
BRIDGE_SHA256 = "2c2f57a166877387244dd5bb6bdf87beb12492ea075a7431939b78e5df7307a0"
BRIDGE_NAME = "athenak-feature-cgl-through-36140ea82.bundle"
CORRUPT_C7_NAME = "athenak-feature-cgl-through-c7e4fa30e.bundle"
PRODUCTION_REQUIRED_REVISIONS = frozenset(
    {
        "9e07542281e4e6d125582f253df3ad2e3b8b154d",
        "c7e4fa30ea7162e4d5ce070a45dfec5b56ca2052",
        "b0d3e8d526d8f3b4e5977333000db24c09d4eab1",
        "38aedd2a65c3c11855721858f5dbcce20bae11e4",
        "5834a91e448a69ec0df5d011b7be3fe786666806",
        "469d38a841ef25d5c044713071adccd55ff90fef",
        "e1f4f4a0b62c3649b4d80a25a188b991f856ebbe",
        BRIDGE_REVISION,
    }
)
REQUIRED_TOOLS = {
    "scripts/frontier/cgl_lf_stage_i.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_checkpoint.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_qualification.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_recost.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_source_authority.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_validate_segment.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_wave_plan.py": "0644",
}
AUTHORIZATION = {
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
F116_PRESERVES = [
    "The immutable F-115 evidence, reviews, publication audit, and historical R03 s02 authority.",
    "The qualified executable, frozen source revision, inputs, matrix, restart lineages, targets, resources, qualification, and Stage I budget policy.",
    "Every prior active source-archive checksum-ledger entry and the corrupt-C7 incident-evidence exclusion.",
]
PRESERVES = [
    "The immutable F-116 evidence, reviews, publication audit, selected source bundle, and nested F-115 historical authority.",
    "The qualified executable, frozen source revision, inputs, matrix, restart lineages, targets, resources, qualification, and Stage I budget policy.",
    "Every prior active source-archive checksum-ledger entry and the corrupt-C7 incident-evidence exclusion.",
]
DOES_NOT_AUTHORIZE = [
    "prepare",
    "submit",
    "direct sbatch",
    "scheduler mutation",
    "Stage I execution-state mutation",
    "scientific configuration change",
    "historical manifest rebinding",
]
PUBLICATION_REQUIREMENTS = {
    "published_evidence_mode": "0444",
    "published_review_mode": "0444",
    "published_audit_mode": "0444",
    "published_links": 1,
    "publication_audit_is_authority_commit_marker": True,
    "recovery_required_after_interruption": True,
}
VALIDATION_CLAIMS = {
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
F116_VALIDATION_CLAIMS = {
    "historical_f115_chain": "passed",
    "bridge_bundle_complete_history": "passed",
    "final_bundle_complete_history": "passed",
    "final_bundle_single_head_tip": "passed",
    "final_bundle_required_revisions": "passed",
    "committed_tool_bytes": "passed",
    "corrupt_c7_exclusion_preserved": True,
}
REVIEW_VERIFIED_KEYS = frozenset(
    {
        "authorization_broadening",
        "bridge_selected_as_current",
        "predecessor_current_source_bundle_selected_as_current",
        "corrupt_c7_excluded",
        "current_source_selection_only",
        "final_bundle_sha256",
        "final_head",
        "historical_f115_preserved",
        "historical_f116_preserved",
    }
)
F116_TRANSACTION_ROOT_NAME = (
    f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_source_authority_transactions"
)
# F118 has a disjoint recovery namespace; retained F116 transactions are immutable history.
TRANSACTION_ROOT_NAME = (
    f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_F118_source_authority_transactions"
)
RETIRED_TRANSACTION_SUFFIX = ".retired"
STAGING_TRANSACTION_SUFFIX = ".staging"
PRIVATE_PUBLICATION_SLOT_COUNT = 2
PRIVATE_PUBLICATION_SUFFIX = ".cgl-source-authority.private"
FORENSIC_ENTRY_RE = re.compile(
    r"\.cgl-source-authority-retired-(?:directory-)?[0-9a-f]{32}\.forensic"
)
OTHER_TRANSACTION_NAMES = (
    f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_transactions",
    f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_transactions",
)
SELF_DESCRIPTOR_ENV = "_CGL_LF_SOURCE_AUTHORITY_DESCRIPTOR"
PYTHON_DESCRIPTOR_ENV = "_CGL_LF_SOURCE_AUTHORITY_PYTHON_DESCRIPTOR"
SELF_SOURCE_ENV = "_CGL_LF_SOURCE_AUTHORITY_SOURCE"
REPOSITORY_ROOT_ENV = "_CGL_LF_SOURCE_AUTHORITY_REPOSITORY_ROOT"
SHA256_RE = re.compile(r"[0-9a-f]{64}")
REVISION_RE = re.compile(r"[0-9a-f]{40}")
BUNDLE_NAME_RE = re.compile(r"athenak-feature-cgl-through-([0-9a-f]{9})\.bundle")
JOURNAL_TEMP_RE = re.compile(r"\.journal\.json\.[0-9]+\.[0-9a-f]{32}\.tmp")
JOURNAL_RECOVERY_NAME = ".journal.json.recovery.tmp"
JOURNAL_RECOVERY_ALTERNATE_NAME = ".journal.json.recovery.alternate.tmp"
JOURNAL_RECOVERY_NAMES = (JOURNAL_RECOVERY_NAME, JOURNAL_RECOVERY_ALTERNATE_NAME)
TRANSACTION_ID_RE = re.compile(
    r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{6}\+0000-[0-9a-f]{32}"
)
TRANSACTION_ID_TIMESTAMP_FORMAT = "%Y-%m-%dT%H%M%S+0000"
GIT = Path("/usr/lib/git/git")
GIT_EXEC_PATH = Path("/usr/lib/git")
TRUSTED_SYSTEM_PATH = "/usr/bin:/bin"
INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION = (
    "Reviewer roles, agent identifiers, and process separation are retained "
    "declarations; exact artifact digests authenticate reviewed bytes but do not "
    "cryptographically authenticate a human or agent identity."
)
RENAME_NOREPLACE = 1
RENAME_EXCHANGE = 2
_ACTIVE_MUTATION_LOCK = None
_ACTIVE_CANONICAL_PUBLIC_NAMESPACE = None
_ACTIVE_F118_PUBLIC_AUTHORITY_LEASE = None
_BOUND_DIRECTORY_GUARDS = {}

# Threat model: bound namespace identities/profiles may change immediately
# before, during, or after any syscall. Mutations revalidate immediately before
# the syscall; possible post-syscall completion permits only durable fsync and
# authenticated classification, never a blind cleanup or second mutation.


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def is_publisher_transaction_id(value: object) -> bool:
    """Return whether a value has the exact transaction ID grammar we generate."""

    if not isinstance(value, str) or TRANSACTION_ID_RE.fullmatch(value) is None:
        return False
    timestamp, _, _ = value.rpartition("-")
    try:
        datetime.strptime(timestamp, TRANSACTION_ID_TIMESTAMP_FORMAT)
    except ValueError:
        return False
    return True


def require_utc(value: object, label: str) -> datetime:
    if not isinstance(value, str) or not value.endswith(("+00:00", "Z")):
        raise ValueError(f"{label} must use UTC")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} is invalid") from error
    if parsed.utcoffset() is None or parsed.utcoffset().total_seconds() != 0:
        raise ValueError(f"{label} must use UTC")
    retained = parsed.astimezone(timezone.utc)
    if (retained - datetime.now(timezone.utc)).total_seconds() > 300:
        raise ValueError(f"{label} is more than five minutes in the future")
    return retained


def require_nonempty(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError(f"{label} must be a nonempty normalized string")
    return value


def require_sha256(value: object, label: str) -> str:
    retained = require_nonempty(value, label)
    if SHA256_RE.fullmatch(retained) is None:
        raise ValueError(f"{label} must be a lowercase SHA-256")
    return retained


def require_revision(value: object, label: str) -> str:
    retained = require_nonempty(value, label)
    if REVISION_RE.fullmatch(retained) is None:
        raise ValueError(f"{label} must be a lowercase full Git revision")
    return retained


def require_exact_keys(value: object, keys: set[str] | frozenset[str], label: str
                       ) -> dict[str, object]:
    if not isinstance(value, dict) or set(value) != set(keys):
        found = sorted(value) if isinstance(value, dict) else type(value).__name__
        raise ValueError(f"{label} schema differs; expected {sorted(keys)}, found {found}")
    return value


def require_relative(value: object, label: str) -> PurePosixPath:
    retained = require_nonempty(value, label)
    path = PurePosixPath(retained)
    if path.is_absolute() or path.as_posix() != retained or ".." in path.parts:
        raise ValueError(f"{label} must be a normalized relative path")
    return path


def normalized_absolute(path: Path, label: str) -> Path:
    path = path.absolute()
    if not path.is_absolute() or path != Path(os.path.normpath(str(path))) or ".." in path.parts:
        raise ValueError(f"{label} must be an absolute normalized path")
    return path


def require_no_symlinks(path: Path, label: str, *, include_leaf: bool = True) -> None:
    path = normalized_absolute(path, label)
    current = Path(path.parts[0])
    parts = path.parts[1:] if include_leaf else path.parts[1:-1]
    for part in parts:
        current /= part
        try:
            profile = os.lstat(current)
        except FileNotFoundError:
            continue
        if stat.S_ISLNK(profile.st_mode):
            raise ValueError(f"{label} path contains a symbolic link: {current}")


@contextmanager
def absolute_descriptor(path: Path, label: str, *, flags: int) -> Iterator[int]:
    """Open one absolute path without accepting symlinks in any component."""

    path = normalized_absolute(path, label)
    if len(path.parts) < 2:
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


def require_file_profile(profile: os.stat_result, label: str, *,
                         mode: int | None = None, links: int = 1) -> None:
    if not stat.S_ISREG(profile.st_mode):
        raise ValueError(f"{label} must be a regular file")
    if profile.st_uid != os.geteuid():
        raise ValueError(f"{label} must be owned by the effective user")
    if profile.st_nlink != links:
        raise ValueError(f"{label} has {profile.st_nlink} links, expected {links}")
    retained_mode = stat.S_IMODE(profile.st_mode)
    if mode is not None and retained_mode != mode:
        raise ValueError(f"{label} mode is {retained_mode:04o}, expected {mode:04o}")
    if retained_mode & 0o022:
        raise ValueError(f"{label} must not be group- or world-writable")


def require_system_executable_profile(profile: os.stat_result, label: str) -> None:
    """Require one root-owned, single-link, non-writable executable."""

    retained_mode = stat.S_IMODE(profile.st_mode)
    if (
        not stat.S_ISREG(profile.st_mode)
        or profile.st_uid != 0
        or profile.st_nlink != 1
        or retained_mode & 0o022
        or retained_mode & 0o111 == 0
    ):
        raise ValueError(f"{label} does not have the trusted system executable profile")


def require_system_directory_profile(profile: os.stat_result, label: str) -> None:
    """Require one root-owned directory without unsafe replacement permissions."""

    retained_mode = stat.S_IMODE(profile.st_mode)
    if (
        not stat.S_ISDIR(profile.st_mode)
        or profile.st_uid != 0
        or retained_mode & 0o022
    ):
        raise ValueError(f"{label} does not have the trusted system directory profile")


def require_directory_profile(profile: os.stat_result, label: str, *,
                              mode: int | None = None) -> None:
    """Require one owned directory with a trusted replacement boundary."""

    if not stat.S_ISDIR(profile.st_mode) or profile.st_uid != os.geteuid():
        raise ValueError(f"{label} must be an owned directory")
    retained_mode = stat.S_IMODE(profile.st_mode)
    if mode is not None and retained_mode not in {mode, mode | stat.S_ISGID}:
        raise ValueError(
            f"{label} mode is {retained_mode:04o}, expected {mode:04o} "
            "with optional inherited setgid"
        )
    if retained_mode & 0o022:
        raise ValueError(f"{label} must not be group- or world-writable")


def require_directory(path: Path, label: str, *, mode: int | None = None) -> None:
    require_no_symlinks(path, label)
    with absolute_descriptor(
        path, label, flags=os.O_RDONLY | os.O_DIRECTORY
    ) as descriptor:
        profile = os.fstat(descriptor)
        require_directory_profile(profile, label, mode=mode)
        named = os.stat(path, follow_symlinks=False)
        if (named.st_dev, named.st_ino) != (profile.st_dev, profile.st_ino):
            raise ValueError(f"{label} path changed while authenticating")


def read_descriptor(descriptor: int, label: str, *, expected: str | None = None,
                    mode: int | None = None, links: int = 1
                    ) -> tuple[bytes, str, os.stat_result]:
    profile = os.fstat(descriptor)
    require_file_profile(profile, label, mode=mode, links=links)
    os.lseek(descriptor, 0, os.SEEK_SET)
    chunks = []
    digest = hashlib.sha256()
    while block := os.read(descriptor, 1024 * 1024):
        chunks.append(block)
        digest.update(block)
    retained = digest.hexdigest()
    if expected is not None and retained != require_sha256(expected, f"{label} SHA-256"):
        raise ValueError(f"{label} checksum differs")
    return b"".join(chunks), retained, profile


def read_file(path: Path, label: str, *, expected: str | None = None,
              mode: int | None = None, links: int = 1) -> tuple[bytes, str]:
    path = normalized_absolute(path, label)
    require_no_symlinks(path, label)
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        payload, retained, _ = read_descriptor(
            descriptor, label, expected=expected, mode=mode, links=links
        )
    finally:
        os.close(descriptor)
    return payload, retained


def read_json(path: Path, label: str, *, expected: str | None = None,
              mode: int | None = None) -> tuple[dict[str, object], bytes, str]:
    payload, digest = read_file(path, label, expected=expected, mode=mode)
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object")
    if canonical_json(value) != payload:
        raise ValueError(f"{label} must be stable sorted/indented JSON ending newline")
    return value, payload, digest


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fsync_descriptor(descriptor: int) -> None:
    """Persist directory-entry changes through an authenticated descriptor."""

    os.fsync(descriptor)
    require_bound_directory_descriptor(descriptor)


def fsync_bound_namespace(*descriptors: int) -> None:
    """Persist every affected directory and then report any authority loss."""

    failures: list[BaseException] = []
    for descriptor in dict.fromkeys(descriptors):
        try:
            fsync_descriptor(descriptor)
        except BaseException as error:
            failures.append(error)
    try:
        require_active_mutation_lock_bound()
    except BaseException as error:
        failures.append(error)
    if failures:
        raise failures[0]


def preserve_namespace_ambiguity(label: str, error: BaseException,
                                 *descriptors: int) -> None:
    """Durably retain a possibly completed namespace syscall before raising."""

    try:
        fsync_bound_namespace(*descriptors)
    except BaseException as fsync_error:
        raise ValueError(
            f"{label} namespace mutation is durably ambiguous after authority loss"
        ) from fsync_error
    raise ValueError(
        f"{label} namespace mutation completed ambiguously after durable directory fsync"
    ) from error


def require_entry_name(name: str, label: str) -> str:
    """Require one direct-child basename for descriptor-relative operations."""

    if not name or name in {".", ".."} or Path(name).name != name or "/" in name:
        raise ValueError(f"{label} must be one direct-child basename")
    return name


def require_bound_entry_identity(parent: int, name: str, expected: os.stat_result,
                                 label: str) -> os.stat_result:
    """Require one directory entry to retain an authenticated inode identity."""

    name = require_entry_name(name, label)
    try:
        named = os.stat(name, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError as error:
        raise ValueError(f"{label} disappeared before mutation") from error
    if (named.st_dev, named.st_ino) != (expected.st_dev, expected.st_ino):
        raise ValueError(f"{label} inode identity changed before mutation")
    return named


def require_bound_entry_absent(parent: int, name: str, label: str) -> None:
    """Require one direct-child name to remain absent."""

    name = require_entry_name(name, label)
    try:
        os.stat(name, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError:
        return
    raise ValueError(f"{label} public name reappeared during mutation")


def profile_identity(profile: os.stat_result) -> tuple[int, int]:
    """Return one stable filesystem object identity."""

    return profile.st_dev, profile.st_ino


def file_profile_binding(profile: os.stat_result) -> tuple[int, ...]:
    """Return immutable security/content metadata for one authenticated file."""

    return (
        stat.S_IFMT(profile.st_mode),
        stat.S_IMODE(profile.st_mode),
        profile.st_uid,
        profile.st_gid,
        profile.st_nlink,
        profile.st_size,
        profile.st_mtime_ns,
        profile.st_ctime_ns,
    )


def renamed_file_profile_binding(profile: os.stat_result) -> tuple[int, ...]:
    """Return file metadata that one same-filesystem rename must preserve."""

    return file_profile_binding(profile)[:-1]


def directory_security_identity(profile: os.stat_result) -> tuple[int, int, int]:
    """Return stable directory security metadata, excluding mutable timestamps."""

    return stat.S_IMODE(profile.st_mode), profile.st_uid, profile.st_gid


def require_bound_entry_profile(parent: int, name: str, expected: os.stat_result,
                                label: str) -> os.stat_result:
    """Require one bound entry to retain its complete mutation-relevant profile."""

    named = require_bound_entry_identity(parent, name, expected, label)
    if stat.S_ISREG(expected.st_mode):
        if not stat.S_ISREG(named.st_mode) or file_profile_binding(named) != file_profile_binding(
            expected
        ):
            raise ValueError(f"{label} file security/content profile changed before mutation")
    elif stat.S_ISDIR(expected.st_mode):
        if (
            not stat.S_ISDIR(named.st_mode)
            or directory_security_identity(named) != directory_security_identity(expected)
        ):
            raise ValueError(f"{label} directory security profile changed before mutation")
    else:
        raise ValueError(f"{label} has an unsupported mutation profile")
    return named


def require_renamed_bound_entry_profile(parent: int, name: str,
                                        expected: os.stat_result,
                                        label: str) -> os.stat_result:
    """Bind one just-renamed entry while allowing only rename-induced ctime drift."""

    named = require_bound_entry_identity(parent, name, expected, label)
    if stat.S_ISREG(expected.st_mode):
        if (
            not stat.S_ISREG(named.st_mode)
            or renamed_file_profile_binding(named)
            != renamed_file_profile_binding(expected)
        ):
            raise ValueError(f"{label} file security/content profile changed during mutation")
    elif stat.S_ISDIR(expected.st_mode):
        if (
            not stat.S_ISDIR(named.st_mode)
            or directory_security_identity(named) != directory_security_identity(expected)
        ):
            raise ValueError(f"{label} directory security profile changed during mutation")
    else:
        raise ValueError(f"{label} has an unsupported mutation profile")
    return named


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
        require_bound_directory_descriptor(self.parent)
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


class DirectoryBindingGuard:
    """Rebind one open directory descriptor to its authenticated pathname."""

    def __init__(self, descriptor: int, profile: os.stat_result, label: str, *,
                 path: Path | None = None, parent: int | None = None,
                 name: str | None = None):
        if (path is None) == (parent is None or name is None):
            raise ValueError(f"{label} directory binding is invalid")
        self.descriptor = descriptor
        self.identity = profile_identity(profile)
        self.security_identity = directory_security_identity(profile)
        self.label = label
        self.path = normalized_absolute(path, label) if path is not None else None
        self.parent = parent
        self.name = require_entry_name(name, label) if name is not None else None

    def assert_bound(self) -> None:
        opened = os.fstat(self.descriptor)
        if (
            not stat.S_ISDIR(opened.st_mode)
            or profile_identity(opened) != self.identity
        ):
            raise ValueError(f"{self.label} descriptor changed while mutation is active")
        if directory_security_identity(opened) != self.security_identity:
            raise ValueError(
                f"{self.label} descriptor security metadata changed while mutation is active"
            )
        if self.path is not None:
            require_no_symlinks(self.path, self.label)
            try:
                named = os.stat(self.path, follow_symlinks=False)
            except FileNotFoundError as error:
                raise ValueError(
                    f"{self.label} path changed while mutation is active"
                ) from error
        else:
            if self.parent is None or self.name is None:
                raise ValueError(f"{self.label} directory binding is invalid")
            require_bound_directory_descriptor(self.parent)
            try:
                named = os.stat(
                    self.name, dir_fd=self.parent, follow_symlinks=False
                )
            except FileNotFoundError as error:
                raise ValueError(
                    f"{self.label} path changed while mutation is active"
                ) from error
        if (
            not stat.S_ISDIR(named.st_mode)
            or profile_identity(named) != self.identity
        ):
            raise ValueError(f"{self.label} path changed while mutation is active")
        if directory_security_identity(named) != self.security_identity:
            raise ValueError(
                f"{self.label} path security metadata changed while mutation is active"
            )


class CanonicalPublicNamespaceGuard:
    """Continuously bind the exact trusted canonical public namespace."""

    def __init__(self, bindings: list[tuple[Path, int, os.stat_result]]):
        self.bindings = bindings

    def assert_bound(self) -> None:
        for path, descriptor, expected in self.bindings:
            require_no_symlinks(path, f"canonical public namespace {path}")
            opened = os.fstat(descriptor)
            require_canonical_public_namespace_profile(path, opened)
            if (
                profile_identity(opened) != profile_identity(expected)
                or directory_security_identity(opened)
                != directory_security_identity(expected)
            ):
                raise ValueError(f"canonical public namespace descriptor changed: {path}")
            named = os.stat(path, follow_symlinks=False)
            require_canonical_public_namespace_profile(path, named)
            if (
                profile_identity(named) != profile_identity(expected)
                or directory_security_identity(named)
                != directory_security_identity(expected)
            ):
                raise ValueError(f"canonical public namespace path changed: {path}")


class F118PublicAuthorityLease:
    """Retain exact public F118 files through one authority-bearing return."""

    def __init__(self, declarations: list[tuple[Path, str, int, str]]):
        if not declarations:
            raise ValueError("F118 public-authority lease declarations are empty")
        self.members: list[dict[str, object]] = []
        self.closed = False
        seen: set[Path] = set()
        try:
            for declared_path, declared_digest, mode, label in declarations:
                path = normalized_absolute(declared_path, label)
                if path in seen:
                    raise ValueError(
                        f"F118 public-authority lease duplicates a path: {path}"
                    )
                seen.add(path)
                digest = require_sha256(
                    declared_digest, f"{label} public-authority lease SHA-256"
                )
                require_no_symlinks(path, label)
                parent = os.open(
                    path.parent,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                )
                descriptor = None
                try:
                    parent_profile = os.fstat(parent)
                    require_directory_profile(parent_profile, f"{label} parent")
                    parent_guard = DirectoryBindingGuard(
                        parent, parent_profile, f"{label} parent", path=path.parent
                    )
                    parent_guard.assert_bound()
                    descriptor = os.open(
                        path.name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent
                    )
                    _, _, profile = read_descriptor(
                        descriptor, label, expected=digest, mode=mode
                    )
                    require_bound_entry_profile(parent, path.name, profile, label)
                    parent_guard.assert_bound()
                except BaseException:
                    if descriptor is not None:
                        os.close(descriptor)
                    os.close(parent)
                    raise
                self.members.append(
                    {
                        "path": path,
                        "digest": digest,
                        "mode": mode,
                        "label": label,
                        "parent": parent,
                        "parent_guard": parent_guard,
                        "descriptor": descriptor,
                        "profile": profile,
                    }
                )
            self.assert_bound()
        except BaseException:
            self.close()
            raise

    def assert_bound(self) -> None:
        """Require every retained public name, inode, profile, and byte stream."""

        if self.closed:
            raise ValueError("F118 public-authority lease is closed")
        for member in self.members:
            label = str(member["label"])
            try:
                parent_guard = member["parent_guard"]
                if not isinstance(parent_guard, DirectoryBindingGuard):
                    raise ValueError("F118 public-authority lease parent guard is invalid")
                parent_guard.assert_bound()
                descriptor = int(member["descriptor"])
                _, _, opened = read_descriptor(
                    descriptor,
                    label,
                    expected=str(member["digest"]),
                    mode=int(member["mode"]),
                )
                expected = member["profile"]
                if not isinstance(expected, os.stat_result):
                    raise ValueError("F118 public-authority lease profile is invalid")
                if (
                    profile_identity(opened) != profile_identity(expected)
                    or file_profile_binding(opened) != file_profile_binding(expected)
                ):
                    raise ValueError("retained descriptor profile changed")
                require_bound_entry_profile(
                    int(member["parent"]),
                    Path(member["path"]).name,
                    expected,
                    label,
                )
                parent_guard.assert_bound()
            except BaseException as error:
                raise ValueError(
                    f"{label} changed while F118 public-authority lease is active"
                ) from error

    def close(self) -> None:
        """Release every retained F118 public-authority descriptor."""

        if self.closed:
            return
        self.closed = True
        for member in reversed(self.members):
            os.close(int(member["descriptor"]))
            os.close(int(member["parent"]))


@contextmanager
def active_f118_public_authority_lease(
    lease: F118PublicAuthorityLease,
) -> Iterator[F118PublicAuthorityLease]:
    """Bind one aggregate F118 lease to mutation boundaries and final return."""

    global _ACTIVE_F118_PUBLIC_AUTHORITY_LEASE

    lease.assert_bound()
    if _ACTIVE_F118_PUBLIC_AUTHORITY_LEASE is not None:
        raise ValueError("nested F118 public-authority leases are forbidden")
    _ACTIVE_F118_PUBLIC_AUTHORITY_LEASE = lease
    try:
        yield lease
        lease.assert_bound()
    finally:
        if _ACTIVE_F118_PUBLIC_AUTHORITY_LEASE is not lease:
            raise ValueError("F118 public-authority lease changed while active")
        _ACTIVE_F118_PUBLIC_AUTHORITY_LEASE = None


@contextmanager
def f118_public_authority_closure(
    declarations: list[tuple[Path, str, int, str]],
) -> Iterator[F118PublicAuthorityLease]:
    """Acquire, activate, and release one aggregate F118 public closure."""

    lease = F118PublicAuthorityLease(declarations)
    try:
        with active_f118_public_authority_lease(lease):
            yield lease
    finally:
        lease.close()


def require_bound_directory_descriptor(descriptor: int) -> None:
    """Require a registered directory descriptor to remain pathname-bound."""

    guard = _BOUND_DIRECTORY_GUARDS.get(descriptor)
    if guard is not None:
        guard.assert_bound()


def require_active_mutation_lock_bound() -> None:
    """Fail if canonical authority or the active mutation lock moved."""

    if _ACTIVE_CANONICAL_PUBLIC_NAMESPACE is not None:
        _ACTIVE_CANONICAL_PUBLIC_NAMESPACE.assert_bound()
    if _ACTIVE_F118_PUBLIC_AUTHORITY_LEASE is not None:
        _ACTIVE_F118_PUBLIC_AUTHORITY_LEASE.assert_bound()
    if _ACTIVE_MUTATION_LOCK is not None:
        _ACTIVE_MUTATION_LOCK.assert_bound()
    if _ACTIVE_F118_PUBLIC_AUTHORITY_LEASE is not None:
        _ACTIVE_F118_PUBLIC_AUTHORITY_LEASE.assert_bound()
    if _ACTIVE_CANONICAL_PUBLIC_NAMESPACE is not None:
        _ACTIVE_CANONICAL_PUBLIC_NAMESPACE.assert_bound()


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
    require_active_mutation_lock_bound()
    require_bound_directory_descriptor(source_parent)
    require_bound_directory_descriptor(target_parent)
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
    """Perform one same-directory descriptor-relative renameat2 operation."""

    renameat2_between(parent, source, parent, target, flags, label)


def rename_bound_noreplace(parent: int, source: str, target: str,
                           expected: os.stat_result, label: str) -> os.stat_result:
    """Move one authenticated entry without replacing a raced target."""

    require_active_mutation_lock_bound()
    require_bound_directory_descriptor(parent)
    source = require_entry_name(source, label)
    target = require_entry_name(target, f"{label} target")
    require_bound_entry_profile(parent, source, expected, label)
    try:
        renameat2(parent, source, target, RENAME_NOREPLACE, label)
    except FileExistsError as error:
        fsync_bound_namespace(parent)
        raise ValueError(f"{label} target already exists") from error
    except BaseException as error:
        preserve_namespace_ambiguity(label, error, parent)
    fsync_bound_namespace(parent)
    renamed = require_renamed_bound_entry_profile(parent, target, expected, label)
    require_bound_entry_absent(parent, source, label)
    fsync_bound_namespace(parent)
    require_bound_entry_profile(parent, target, renamed, label)
    require_bound_entry_absent(parent, source, label)
    require_bound_directory_descriptor(parent)
    require_active_mutation_lock_bound()
    return renamed


def exchange_bound_entries(parent: int, source: str, target: str,
                           source_expected: os.stat_result,
                           target_expected: os.stat_result, label: str, *,
                           source_content: tuple[str, int, bytes | None, int] | None = None,
                           target_content: tuple[str, int, bytes | None, int] | None = None,
                           ) -> tuple[os.stat_result, os.stat_result]:
    """Exchange exact descriptor-bound inodes and retain all names on ambiguity."""

    require_active_mutation_lock_bound()
    require_bound_directory_descriptor(parent)
    source = require_entry_name(source, label)
    target = require_entry_name(target, f"{label} target")
    if source_content is not None:
        source_expected = authenticate_bound_file_content(
            parent, source, source_expected, label, source_content
        )
    if target_content is not None:
        target_expected = authenticate_bound_file_content(
            parent, target, target_expected, f"{label} target", target_content
        )
    require_bound_entry_profile(parent, source, source_expected, label)
    require_bound_entry_profile(parent, target, target_expected, f"{label} target")
    try:
        renameat2(parent, source, target, RENAME_EXCHANGE, label)
    except BaseException as error:
        preserve_namespace_ambiguity(label, error, parent)
    fsync_bound_namespace(parent)
    published = require_renamed_bound_entry_profile(parent, target, source_expected, label)
    predecessor = require_renamed_bound_entry_profile(
        parent, source, target_expected, f"{label} predecessor"
    )
    if source_content is not None:
        published = authenticate_bound_file_content(
            parent, target, published, label, source_content
        )
    if target_content is not None:
        predecessor = authenticate_bound_file_content(
            parent, source, predecessor, f"{label} predecessor", target_content
        )
    fsync_bound_namespace(parent)
    require_bound_entry_profile(parent, target, published, label)
    require_bound_entry_profile(parent, source, predecessor, f"{label} predecessor")
    require_bound_directory_descriptor(parent)
    require_active_mutation_lock_bound()
    return published, predecessor


def preserve_bound_created_file_ambiguity(
    parent: int,
    name: str,
    label: str,
    error: BaseException,
    expected: os.stat_result | None,
) -> None:
    """Durably classify one ambiguous create without unlinking its pathname."""

    try:
        fsync_bound_namespace(parent)
    except BaseException as authority_error:
        raise ValueError(f"{label} create is durably ambiguous after authority loss") from (
            authority_error
        )
    try:
        observed = os.stat(name, dir_fd=parent, follow_symlinks=False)
        require_file_profile(observed, f"{label} retained create")
        require_bound_entry_profile(
            parent,
            name,
            expected if expected is not None else observed,
            f"{label} retained create",
        )
    except FileNotFoundError as missing:
        raise ValueError(f"{label} create failed without a durable retained entry") from missing
    except BaseException as profile_error:
        raise ValueError(f"{label} create left an unauthenticated durable entry") from (
            profile_error
        )
    try:
        fsync_bound_namespace(parent)
    except BaseException as authority_error:
        raise ValueError(f"{label} create is durably ambiguous after authority loss") from (
            authority_error
        )
    try:
        require_bound_entry_profile(parent, name, observed, f"{label} retained create")
    except BaseException as profile_error:
        raise ValueError(f"{label} create retained entry changed after durable fsync") from (
            profile_error
        )
    raise ValueError(
        f"{label} create completed ambiguously; retained authenticated recovery entry"
    ) from error


def preserve_bound_created_directory_ambiguity(
    parent: int, name: str, label: str, error: BaseException
) -> None:
    """Durably classify one ambiguous mkdir without removing its namespace."""

    try:
        fsync_bound_namespace(parent)
    except BaseException as authority_error:
        raise ValueError(f"{label} mkdir is durably ambiguous after authority loss") from (
            authority_error
        )
    try:
        require_active_mutation_lock_bound()
        with bound_child_directory(
            parent, name, f"{label} retained mkdir"
        ) as (_, observed):
            fsync_bound_namespace(parent)
            require_bound_entry_profile(
                parent, name, observed, f"{label} retained mkdir"
            )
    except BaseException as profile_error:
        raise ValueError(f"{label} mkdir left an unauthenticated durable entry") from (
            profile_error
        )
    raise ValueError(
        f"{label} mkdir completed ambiguously; retained authenticated recovery directory"
    ) from error


def mkdir_bound_exclusive(parent: int, name: str, mode: int,
                          label: str) -> os.stat_result:
    """Create and durably authenticate one direct-child directory."""

    require_active_mutation_lock_bound()
    require_bound_directory_descriptor(parent)
    name = require_entry_name(name, label)
    previous = os.umask(0)
    try:
        require_active_mutation_lock_bound()
        require_bound_directory_descriptor(parent)
        try:
            os.mkdir(name, mode=mode, dir_fd=parent)
        except BaseException as error:
            preserve_bound_created_directory_ambiguity(parent, name, label, error)
    finally:
        os.umask(previous)
    try:
        require_active_mutation_lock_bound()
        require_bound_directory_descriptor(parent)
        with bound_child_directory(parent, name, label, mode=mode) as (_, profile):
            fsync_bound_namespace(parent)
            require_bound_entry_profile(parent, name, profile, label)
        require_bound_directory_descriptor(parent)
        require_active_mutation_lock_bound()
        return profile
    except BaseException as error:
        preserve_bound_created_directory_ambiguity(parent, name, label, error)


def write_bound_exclusive(parent: int, name: str, payload: bytes, mode: int,
                          label: str) -> os.stat_result:
    """Create exact bytes while retaining any authenticated ambiguous remnant."""

    require_active_mutation_lock_bound()
    require_bound_directory_descriptor(parent)
    name = require_entry_name(name, label)
    descriptor: int | None = None
    retained_profile: os.stat_result | None = None
    previous = os.umask(0)
    try:
        require_active_mutation_lock_bound()
        require_bound_directory_descriptor(parent)
        try:
            descriptor = os.open(
                name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o000,
                dir_fd=parent,
            )
        except BaseException as error:
            preserve_bound_created_file_ambiguity(parent, name, label, error, None)
    finally:
        os.umask(previous)
    try:
        if descriptor is None:
            raise ValueError(f"{label} create did not return a descriptor")
        retained_profile = os.fstat(descriptor)
        require_file_profile(retained_profile, label, mode=0o000)
        require_bound_entry_profile(parent, name, retained_profile, label)
        require_active_mutation_lock_bound()
        require_bound_directory_descriptor(parent)
        offset = 0
        while offset < len(payload):
            require_active_mutation_lock_bound()
            require_bound_directory_descriptor(parent)
            written = os.write(descriptor, payload[offset:])
            require_active_mutation_lock_bound()
            require_bound_directory_descriptor(parent)
            if written <= 0:
                raise ValueError(f"{label} write made no progress")
            offset += written
        require_active_mutation_lock_bound()
        os.fchmod(descriptor, mode)
        require_active_mutation_lock_bound()
        os.fsync(descriptor)
        retained_profile = os.fstat(descriptor)
        require_file_profile(retained_profile, label, mode=mode)
        os.close(descriptor)
        descriptor = None
        require_bound_entry_profile(parent, name, retained_profile, label)
        fsync_bound_namespace(parent)
        require_bound_entry_profile(parent, name, retained_profile, label)
        require_bound_directory_descriptor(parent)
        require_active_mutation_lock_bound()
        return retained_profile
    except BaseException as error:
        if descriptor is not None:
            try:
                os.fsync(descriptor)
            except BaseException:
                pass
            try:
                retained_profile = os.fstat(descriptor)
            except BaseException:
                pass
            try:
                os.close(descriptor)
            except BaseException:
                pass
        preserve_bound_created_file_ambiguity(
            parent, name, label, error, retained_profile
        )


@contextmanager
def bound_directory(path: Path, label: str, *,
                    mode: int | None = None) -> Iterator[tuple[int, os.stat_result]]:
    """Yield one pathname-bound authenticated directory descriptor."""

    require_no_symlinks(path, label)
    with absolute_descriptor(
        path, label, flags=os.O_RDONLY | os.O_DIRECTORY
    ) as descriptor:
        profile = os.fstat(descriptor)
        require_directory_profile(profile, label, mode=mode)
        named = os.stat(path, follow_symlinks=False)
        if (named.st_dev, named.st_ino) != (profile.st_dev, profile.st_ino):
            raise ValueError(f"{label} path changed while authenticating")
        guard = DirectoryBindingGuard(descriptor, profile, label, path=path)
        if descriptor in _BOUND_DIRECTORY_GUARDS:
            raise ValueError(f"{label} directory descriptor is already bound")
        guard.assert_bound()
        _BOUND_DIRECTORY_GUARDS[descriptor] = guard
        try:
            yield descriptor, profile
            guard.assert_bound()
        finally:
            if _BOUND_DIRECTORY_GUARDS.get(descriptor) is guard:
                del _BOUND_DIRECTORY_GUARDS[descriptor]


@contextmanager
def bound_child_directory(parent: int, name: str, label: str, *,
                          mode: int | None = None
                          ) -> Iterator[tuple[int, os.stat_result]]:
    """Yield one authenticated direct-child directory descriptor."""

    name = require_entry_name(name, label)
    try:
        descriptor = os.open(
            name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=parent
        )
    except OSError as error:
        raise ValueError(f"{label} is unavailable") from error
    try:
        profile = os.fstat(descriptor)
        require_directory_profile(profile, label, mode=mode)
        require_bound_entry_identity(parent, name, profile, label)
        guard = DirectoryBindingGuard(
            descriptor, profile, label, parent=parent, name=name
        )
        if descriptor in _BOUND_DIRECTORY_GUARDS:
            raise ValueError(f"{label} directory descriptor is already bound")
        guard.assert_bound()
        _BOUND_DIRECTORY_GUARDS[descriptor] = guard
        try:
            yield descriptor, profile
            guard.assert_bound()
        finally:
            if _BOUND_DIRECTORY_GUARDS.get(descriptor) is guard:
                del _BOUND_DIRECTORY_GUARDS[descriptor]
    finally:
        os.close(descriptor)


def read_bound_file(parent: int, name: str, label: str, *,
                    expected: str | None = None, mode: int | None = None,
                    links: int = 1) -> tuple[bytes, str, os.stat_result]:
    """Read one authenticated direct-child file and retain its inode identity."""

    name = require_entry_name(name, label)
    try:
        descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent)
    except OSError as error:
        raise ValueError(f"{label} is unavailable") from error
    try:
        payload, digest, profile = read_descriptor(
            descriptor, label, expected=expected, mode=mode, links=links
        )
        require_bound_entry_profile(parent, name, profile, label)
        return payload, digest, profile
    finally:
        os.close(descriptor)


def read_bound_json(parent: int, name: str, label: str, *,
                    expected: str | None = None, mode: int | None = None
                    ) -> tuple[dict[str, object], bytes, str, os.stat_result]:
    """Read one authenticated direct-child stable JSON file."""

    payload, digest, profile = read_bound_file(
        parent, name, label, expected=expected, mode=mode
    )
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON") from error
    if not isinstance(value, dict) or canonical_json(value) != payload:
        raise ValueError(f"{label} must be stable sorted/indented JSON ending newline")
    return value, payload, digest, profile


def authenticate_bound_file_content(
    parent: int,
    name: str,
    expected_profile: os.stat_result,
    label: str,
    binding: tuple[str, int, bytes | None, int],
) -> os.stat_result:
    """Bind one file's complete profile and reviewed content at a mutation boundary."""

    expected_digest, mode, expected_payload, links = binding
    require_bound_entry_profile(parent, name, expected_profile, label)
    payload, _, authenticated = read_bound_file(
        parent,
        name,
        label,
        expected=expected_digest,
        mode=mode,
        links=links,
    )
    if expected_payload is not None and payload != expected_payload:
        raise ValueError(f"{label} payload differs")
    require_bound_entry_profile(parent, name, authenticated, label)
    return authenticated


def durably_authenticate_bound_file(
    parent: int,
    name: str,
    label: str,
    *,
    expected: str,
    mode: int,
    payload: bytes | None,
    links: int = 1,
) -> os.stat_result:
    """Authenticate exact bytes on both sides of the publication fsync."""

    first_payload, _, first = read_bound_file(
        parent, name, label, expected=expected, mode=mode, links=links
    )
    if payload is not None and first_payload != payload:
        raise ValueError(f"{label} payload differs")
    fsync_descriptor(parent)
    second_payload, _, second = read_bound_file(
        parent, name, label, expected=expected, mode=mode, links=links
    )
    if payload is not None and second_payload != payload:
        raise ValueError(f"{label} payload differs")
    if file_profile_binding(second) != file_profile_binding(first):
        raise ValueError(f"{label} file security/content profile changed")
    require_bound_directory_descriptor(parent)
    require_active_mutation_lock_bound()
    return second


@contextmanager
def bound_forensic_parent(parent: int, label: str) -> Iterator[int]:
    """Yield the exact bound namespace that owns one mutation parent."""

    guard = _BOUND_DIRECTORY_GUARDS.get(parent)
    if guard is None:
        raise ValueError(f"{label} requires an authenticated mutation parent")
    require_bound_directory_descriptor(parent)
    if guard.parent is not None:
        forensic_parent = guard.parent
        if forensic_parent not in _BOUND_DIRECTORY_GUARDS:
            raise ValueError(f"{label} forensic parent is not independently bound")
        require_bound_directory_descriptor(forensic_parent)
        try:
            yield forensic_parent
        finally:
            require_bound_directory_descriptor(forensic_parent)
        return
    if guard.path is None:
        raise ValueError(f"{label} mutation parent binding is incomplete")
    with bound_directory(guard.path.parent, f"{label} forensic parent") as (
        forensic_parent,
        _,
    ):
        require_bound_directory_descriptor(parent)
        try:
            yield forensic_parent
        finally:
            require_bound_directory_descriptor(forensic_parent)


def unlink_bound_entry(parent: int, name: str, expected: os.stat_result,
                       label: str) -> None:
    """Atomically retire one exact child to a no-clobber forensic sibling."""

    require_active_mutation_lock_bound()
    require_bound_directory_descriptor(parent)
    name = require_entry_name(name, label)
    require_bound_entry_profile(parent, name, expected, label)
    retired = f".cgl-source-authority-retired-{uuid.uuid4().hex}.forensic"
    with bound_forensic_parent(parent, label) as forensic_parent:
        require_bound_entry_profile(parent, name, expected, label)
        try:
            renameat2_between(
                parent,
                name,
                forensic_parent,
                retired,
                RENAME_NOREPLACE,
                f"{label} retirement",
            )
        except FileExistsError as error:
            fsync_bound_namespace(parent, forensic_parent)
            raise ValueError(f"{label} forensic retirement target already exists") from error
        except BaseException as error:
            preserve_namespace_ambiguity(f"{label} retirement", error, parent, forensic_parent)
        fsync_bound_namespace(parent, forensic_parent)
        try:
            retired_profile = require_renamed_bound_entry_profile(
                forensic_parent, retired, expected, f"{label} retired inode"
            )
            require_bound_entry_absent(parent, name, label)
            fsync_bound_namespace(parent, forensic_parent)
            require_bound_entry_profile(
                forensic_parent, retired, retired_profile, f"{label} retired inode"
            )
            require_bound_entry_absent(parent, name, label)
            require_bound_directory_descriptor(parent)
            require_active_mutation_lock_bound()
        except BaseException as error:
            os.fsync(parent)
            os.fsync(forensic_parent)
            require_bound_directory_descriptor(forensic_parent)
            require_active_mutation_lock_bound()
            raise ValueError(
                f"{label} changed during atomic retirement; retained as {retired}"
            ) from error


def retire_bound_recovery_file(parent: int, name: str, observed: os.stat_result,
                               label: str, *, links: int = 1) -> None:
    """Forensically retire one trusted deterministic recovery-file remnant."""

    require_file_profile(observed, label, links=links)
    require_bound_entry_profile(parent, name, observed, label)
    unlink_bound_entry(parent, name, observed, label)


def authenticate_or_retire_recovery_file(
    parent: int,
    name: str,
    observed: os.stat_result | None,
    label: str,
    *,
    expected: str,
    mode: int,
    payload: bytes | None,
    links: int = 1,
) -> os.stat_result | None:
    """Authenticate one deterministic recovery file or preserve it forensics."""

    if observed is None:
        return None
    require_file_profile(observed, label, links=links)
    require_bound_entry_profile(parent, name, observed, label)
    try:
        recovered_payload, _, recovered = read_bound_file(
            parent, name, label, expected=expected, mode=mode, links=links
        )
        require_bound_entry_profile(parent, name, recovered, label)
        if file_profile_binding(recovered) != file_profile_binding(observed):
            raise ValueError(f"{label} file security/content profile changed")
        if payload is not None and recovered_payload != payload:
            raise ValueError(f"{label} payload differs")
    except ValueError:
        # The deterministic name can retain mode-000, short, or otherwise
        # incomplete bytes after a crash. Preserve only the exact observed
        # owned regular inode; pathname substitutions fail retirement closed.
        retire_bound_recovery_file(
            parent, name, observed, f"{label} incomplete", links=links
        )
        return None
    return recovered


def publish_bound_file_noreplace(parent: int, source: str, target: str,
                                 expected: os.stat_result, label: str, *,
                                 content: tuple[str, int, bytes | None, int] | None = None,
                                 ) -> os.stat_result:
    """Publish one exact file inode by no-clobber atomic rename."""

    source = require_entry_name(source, label)
    target = require_entry_name(target, f"{label} target")
    if content is not None:
        expected = authenticate_bound_file_content(parent, source, expected, label, content)
    published = rename_bound_noreplace(parent, source, target, expected, label)
    if content is not None:
        published = authenticate_bound_file_content(
            parent, target, published, label, content
        )
    require_bound_entry_profile(parent, target, published, label)
    require_bound_entry_absent(parent, source, label)
    fsync_descriptor(parent)
    require_bound_entry_profile(parent, target, published, label)
    require_bound_entry_absent(parent, source, label)
    require_bound_directory_descriptor(parent)
    require_active_mutation_lock_bound()
    return published


def replace_bound_entry(parent: int, source: str, target: str,
                        expected: os.stat_result, label: str) -> None:
    """Publish one exact direct-child file without replacing another entry."""

    publish_bound_file_noreplace(parent, source, target, expected, label)


def rename_bound_entry(parent: int, source: str, target: str,
                       expected: os.stat_result, label: str) -> None:
    """Rename one authenticated direct child without replacing another entry."""

    rename_bound_noreplace(parent, source, target, expected, label)


def directory_content_bindings(parent: int, label: str) -> dict[str, tuple[object, ...]]:
    """Bind exact direct-child names, identities, and security/content profiles."""

    bindings: dict[str, tuple[object, ...]] = {}
    for name in sorted(os.listdir(parent)):
        profile = os.stat(name, dir_fd=parent, follow_symlinks=False)
        entry_label = f"{label} entry {name}"
        require_bound_entry_profile(parent, name, profile, entry_label)
        if stat.S_ISREG(profile.st_mode):
            binding: tuple[object, ...] = ("file", *file_profile_binding(profile))
        elif stat.S_ISDIR(profile.st_mode):
            binding = ("directory", *directory_security_identity(profile))
        else:
            raise ValueError(f"{entry_label} has an unsupported content profile")
        bindings[name] = binding
    return bindings


def rollback_bound_directory_retirement(
    parent: int,
    forensic_parent: int,
    name: str,
    retired: str,
    retired_profile: os.stat_result,
    label: str,
) -> None:
    """Restore a changed retired directory only while both namespaces remain bound."""

    require_active_mutation_lock_bound()
    require_bound_directory_descriptor(forensic_parent)
    require_bound_directory_descriptor(parent)
    require_bound_entry_profile(
        forensic_parent,
        retired,
        retired_profile,
        f"{label} changed retired directory",
    )
    require_bound_entry_absent(parent, name, label)
    try:
        renameat2_between(
            forensic_parent,
            retired,
            parent,
            name,
            RENAME_NOREPLACE,
            f"{label} contents rollback",
        )
    except FileExistsError as error:
        fsync_bound_namespace(forensic_parent, parent)
        raise ValueError(f"{label} contents rollback target already exists") from error
    except BaseException as error:
        preserve_namespace_ambiguity(
            f"{label} contents rollback", error, forensic_parent, parent
        )
    fsync_bound_namespace(forensic_parent, parent)
    require_renamed_bound_entry_profile(parent, name, retired_profile, label)
    raise ValueError(f"{label} contents changed during atomic retirement; restored")


def rmdir_bound_path(path: Path, expected: os.stat_result, label: str, *,
                     expected_contents: dict[str, tuple[object, ...]] | None = None
                     ) -> None:
    """Retire one authenticated directory entry and persist its parent."""

    with bound_directory(path.parent, f"{label} parent") as (parent, _):
        rmdir_bound_entry(
            parent, path.name, expected, label, expected_contents=expected_contents
        )


def rmdir_bound_entry(parent: int, name: str, expected: os.stat_result,
                      label: str, *,
                      expected_contents: dict[str, tuple[object, ...]] | None = None
                      ) -> None:
    """Atomically retire one exact directory to a no-clobber forensic sibling."""

    require_active_mutation_lock_bound()
    require_bound_directory_descriptor(parent)
    name = require_entry_name(name, label)
    require_bound_entry_profile(parent, name, expected, label)
    contents_descriptor: int | None = None
    if expected_contents is not None:
        try:
            contents_descriptor = os.open(
                name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=parent
            )
        except OSError as error:
            raise ValueError(f"{label} contents are unavailable") from error
        opened = os.fstat(contents_descriptor)
        require_directory_profile(opened, label)
        if (
            profile_identity(opened) != profile_identity(expected)
            or directory_security_identity(opened) != directory_security_identity(expected)
            or directory_content_bindings(contents_descriptor, label) != expected_contents
        ):
            os.close(contents_descriptor)
            raise ValueError(f"{label} contents changed before atomic retirement")
    retired = f".cgl-source-authority-retired-directory-{uuid.uuid4().hex}.forensic"
    try:
        with bound_forensic_parent(parent, label) as forensic_parent:
            require_bound_entry_profile(parent, name, expected, label)
            try:
                renameat2_between(
                    parent,
                    name,
                    forensic_parent,
                    retired,
                    RENAME_NOREPLACE,
                    f"{label} retirement",
                )
            except FileExistsError as error:
                fsync_bound_namespace(parent, forensic_parent)
                raise ValueError(f"{label} forensic retirement target already exists") from error
            except BaseException as error:
                preserve_namespace_ambiguity(
                    f"{label} retirement", error, parent, forensic_parent
                )
            fsync_bound_namespace(parent, forensic_parent)
            try:
                retired_profile = require_renamed_bound_entry_profile(
                    forensic_parent, retired, expected, f"{label} retired directory"
                )
                require_bound_entry_absent(parent, name, label)
                if (
                    contents_descriptor is not None
                    and directory_content_bindings(contents_descriptor, label)
                    != expected_contents
                ):
                    rollback_bound_directory_retirement(
                        parent,
                        forensic_parent,
                        name,
                        retired,
                        retired_profile,
                        label,
                    )
                fsync_bound_namespace(parent, forensic_parent)
                require_bound_entry_profile(
                    forensic_parent,
                    retired,
                    retired_profile,
                    f"{label} retired directory",
                )
                require_bound_entry_absent(parent, name, label)
                if (
                    contents_descriptor is not None
                    and directory_content_bindings(contents_descriptor, label)
                    != expected_contents
                ):
                    rollback_bound_directory_retirement(
                        parent,
                        forensic_parent,
                        name,
                        retired,
                        retired_profile,
                        label,
                    )
                require_bound_directory_descriptor(parent)
                require_active_mutation_lock_bound()
            except BaseException as error:
                os.fsync(parent)
                os.fsync(forensic_parent)
                require_bound_directory_descriptor(forensic_parent)
                require_active_mutation_lock_bound()
                if isinstance(error, ValueError) and str(error).endswith("; restored"):
                    raise
                raise ValueError(
                    f"{label} changed during atomic retirement; retained as {retired}"
                ) from error
    finally:
        if contents_descriptor is not None:
            os.close(contents_descriptor)


def atomic_recovery_names(path: Path) -> tuple[str, str]:
    """Return the two deterministic recovery slots for one atomic target."""

    return f".{path.name}.recovery.tmp", f".{path.name}.recovery.alternate.tmp"


def prepare_atomic_recovery_copy(
    parent: int,
    names: tuple[str, str],
    payload: bytes,
    expected: str,
    mode: int,
    label: str,
) -> tuple[str, os.stat_result]:
    """Install one exact recovery before retiring any older recovery slot."""

    valid: dict[str, tuple[bytes, str, os.stat_result]] = {}
    invalid: dict[str, os.stat_result] = {}
    absent: list[str] = []
    for name in names:
        try:
            observed = os.stat(name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            absent.append(name)
            continue
        require_file_profile(observed, f"{label} recovery slot")
        require_bound_entry_profile(parent, name, observed, f"{label} recovery slot")
        try:
            recovered_payload, recovered_digest, recovered = read_bound_file(
                parent, name, f"{label} recovery slot", mode=mode
            )
            if file_profile_binding(recovered) != file_profile_binding(observed):
                raise ValueError(f"{label} recovery slot profile changed")
        except ValueError:
            invalid[name] = observed
            continue
        valid[name] = (recovered_payload, recovered_digest, recovered)

    selected = next(
        (
            (name, recovered)
            for name, (recovered_payload, recovered_digest, recovered) in valid.items()
            if recovered_payload == payload and recovered_digest == expected
        ),
        None,
    )
    if selected is None:
        if absent:
            selected_name = absent[0]
        elif invalid and valid:
            selected_name = next(iter(invalid))
            retained_name, (retained_payload, retained_digest, _) = next(iter(valid.items()))
            durably_authenticate_bound_file(
                parent,
                retained_name,
                f"{label} retained recovery",
                expected=retained_digest,
                mode=mode,
                payload=retained_payload,
            )
            retire_bound_recovery_file(
                parent,
                selected_name,
                invalid[selected_name],
                f"{label} invalid recovery slot",
            )
        else:
            raise ValueError(f"{label} has no safe deterministic recovery slot")
        selected_profile = write_bound_exclusive(
            parent, selected_name, payload, mode, f"{label} recovery"
        )
    else:
        selected_name, selected_profile = selected

    selected_profile = durably_authenticate_bound_file(
        parent,
        selected_name,
        f"{label} selected recovery",
        expected=expected,
        mode=mode,
        payload=payload,
    )
    for name in names:
        if name == selected_name:
            continue
        try:
            observed = os.stat(name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            continue
        retire_bound_recovery_file(parent, name, observed, f"{label} superseded recovery")
        selected_profile = durably_authenticate_bound_file(
            parent,
            selected_name,
            f"{label} selected recovery",
            expected=expected,
            mode=mode,
            payload=payload,
        )
    return selected_name, selected_profile


def atomic_write(path: Path, payload: bytes, mode: int) -> None:
    """Atomically write through a stable parent without clobbering a raced inode."""

    with bound_directory(path.parent, f"{path.parent} directory") as (parent, _):
        temporary = f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp"
        recovery_names = atomic_recovery_names(path)
        try:
            target_profile = os.stat(path.name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            target_profile = None
        if target_profile is not None:
            predecessor_payload, predecessor_digest, target_profile = read_bound_file(
                parent, path.name, f"{path.name} atomic predecessor", mode=mode
            )
        else:
            predecessor_payload = payload
            predecessor_digest = sha256_bytes(payload)
        recovery, recovery_profile = prepare_atomic_recovery_copy(
            parent,
            recovery_names,
            predecessor_payload,
            predecessor_digest,
            mode,
            f"{path.name} atomic predecessor",
        )
        if target_profile is not None and predecessor_payload == payload:
            durably_authenticate_bound_file(
                parent,
                path.name,
                str(path),
                expected=predecessor_digest,
                mode=mode,
                payload=predecessor_payload,
            )
            return
        temporary_profile = write_bound_exclusive(
            parent, temporary, payload, mode, f"{path.name} atomic temporary"
        )
        try:
            if target_profile is None:
                published_profile = publish_bound_file_noreplace(
                    parent, temporary, path.name, temporary_profile, f"{path.name} atomic write"
                )
            else:
                require_file_profile(target_profile, f"{path.name} atomic predecessor")
                published_profile, predecessor_profile = exchange_bound_entries(
                    parent,
                    temporary,
                    path.name,
                    temporary_profile,
                    target_profile,
                    f"{path.name} atomic write",
                )
        except BaseException:
            # Retain every temporary or predecessor on ambiguity. Recovery can
            # authenticate it later; rollback never unlinks a raced pathname.
            raise
        durably_authenticate_bound_file(
            parent,
            path.name,
            str(path),
            expected=sha256_bytes(payload),
            mode=mode,
            payload=payload,
        )
        if target_profile is not None:
            unlink_bound_entry(
                parent,
                temporary,
                predecessor_profile,
                f"{path.name} replaced predecessor",
            )
            durably_authenticate_bound_file(
                parent,
                path.name,
                str(path),
                expected=sha256_bytes(payload),
                mode=mode,
                payload=payload,
            )
        durably_authenticate_bound_file(
            parent,
            recovery,
            f"{path.name} retained atomic predecessor recovery",
            expected=predecessor_digest,
            mode=mode,
            payload=predecessor_payload,
        )
        require_bound_directory_descriptor(parent)
        require_active_mutation_lock_bound()


def write_exclusive(path: Path, payload: bytes, mode: int) -> None:
    """Create one exact direct child through a stable parent descriptor."""

    with bound_directory(path.parent, f"{path.parent} directory") as (parent, _):
        write_bound_exclusive(parent, path.name, payload, mode, str(path))


def unlink_durable(path: Path) -> None:
    """Retire only the exact inode authenticated through a stable parent."""

    with bound_directory(path.parent, f"{path.parent} directory") as (parent, _):
        _, _, profile = read_bound_file(parent, path.name, str(path))
        unlink_bound_entry(parent, path.name, profile, str(path))


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


def hardened_git_arguments(repository: Path, arguments: list[str]) -> list[str]:
    """Pin repository mode/worktree and disable local-config execution surfaces."""

    require_directory(repository, "Git repository")
    metadata = repository / ".git"
    retained = list(arguments)
    initializing = bool(retained) and retained[0] == "init"
    if os.path.lexists(metadata):
        require_directory(metadata, "Git metadata directory")
        repository_location = [
            f"--git-dir={metadata}",
            f"--work-tree={repository}",
        ]
    else:
        repository_location = [f"--git-dir={repository}"]
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


def git_run(repository: Path, arguments: list[str], *, capture_output: bool = True,
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
                        str(repository),
                        *hardened_git_arguments(repository, arguments),
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


def exact_option(argv: list[str], option: str) -> str:
    if any(item.startswith(f"{option}=") for item in argv):
        raise ValueError(f"{option} must use one separate exact value")
    matches = [argv[index + 1] for index, item in enumerate(argv[:-1]) if item == option]
    if len(matches) != 1:
        raise ValueError(f"{option} must be supplied exactly once")
    return matches[0]


def initial_source() -> Path:
    """Return the named publisher source for initial descriptor re-execution."""

    source = Path(__file__).absolute()
    if source.is_symlink():
        raise ValueError("publisher source must not be a symbolic link")
    return source.resolve()


def initial_repository(source: Path) -> Path:
    """Return the named repository for initial descriptor re-execution."""

    try:
        return source.parents[2].resolve()
    except IndexError as error:
        raise ValueError(f"publisher source path is invalid: {source}") from error


def require_reexec_source_relationship(source: Path, repository: Path) -> None:
    """Bind private reexec metadata to the publisher's repository-relative path."""

    if source != repository / PUBLISHER_RELATIVE:
        raise ValueError("authenticated publisher source/repository relationship differs")


def inherited_reexec_path(name: str, label: str) -> Path:
    """Read one private reexec path only after authenticating the self descriptor."""

    retained = os.environ.get(name)
    if retained is None:
        raise ValueError(f"authenticated publisher reexecution lacks {label}")
    path = Path(retained)
    if (
        not path.is_absolute()
        or path != Path(os.path.normpath(retained))
        or ".." in path.parts
        or len(path.parts) < 2
    ):
        raise ValueError(f"authenticated publisher {label} is not normalized and absolute")
    require_no_symlinks(path, f"authenticated publisher {label}")
    return path


def reexec_environment(
    descriptor: int, python_descriptor: int, source: Path, repository: Path
) -> dict[str, str]:
    """Return the complete caller-independent authenticated reexec environment."""

    return {
        SELF_DESCRIPTOR_ENV: str(descriptor),
        PYTHON_DESCRIPTOR_ENV: str(python_descriptor),
        SELF_SOURCE_ENV: str(source),
        REPOSITORY_ROOT_ENV: str(repository),
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": TRUSTED_SYSTEM_PATH,
        "PYTHONDONTWRITEBYTECODE": "1",
        "XDG_CONFIG_HOME": "/nonexistent",
    }


def require_authenticated_python_descriptor() -> int:
    """Authenticate the inherited interpreter descriptor against this process."""

    retained = os.environ.get(PYTHON_DESCRIPTOR_ENV)
    if retained is None or re.fullmatch(r"[0-9]+", retained) is None:
        raise ValueError("authenticated publisher Python descriptor marker is invalid")
    descriptor = int(retained)
    profile = os.fstat(descriptor)
    require_system_executable_profile(profile, "authenticated publisher Python interpreter")
    running = os.stat("/proc/self/exe")
    if profile_identity(profile) != profile_identity(running):
        raise ValueError("authenticated publisher Python descriptor is not this interpreter")
    if not sys.flags.isolated:
        raise ValueError("authenticated publisher Python interpreter is not isolated")
    return descriptor


def authenticate_self(argv: list[str]) -> tuple[Path, Path, str]:
    expected = require_sha256(exact_option(argv, "--expected-publisher-sha256"),
                              "expected publisher SHA-256")
    inherited = os.environ.get(SELF_DESCRIPTOR_ENV)
    if inherited is not None:
        require_authenticated_python_descriptor()
        if re.fullmatch(r"[0-9]+", inherited) is None:
            raise ValueError("publisher descriptor marker is invalid")
        try:
            descriptor = int(inherited)
        except ValueError as error:
            raise ValueError("publisher descriptor marker is invalid") from error
        if __file__ != f"/proc/self/fd/{descriptor}":
            raise ValueError("publisher descriptor marker is not attached to this execution")
        require_file_profile(os.fstat(descriptor), "authenticated publisher", mode=0o755)
        os.lseek(descriptor, 0, os.SEEK_SET)
        digest = hashlib.sha256()
        while block := os.read(descriptor, 1024 * 1024):
            digest.update(block)
        if digest.hexdigest() != expected:
            raise ValueError("authenticated publisher descriptor checksum differs")
        source = inherited_reexec_path(SELF_SOURCE_ENV, "source path")
        repository = inherited_reexec_path(REPOSITORY_ROOT_ENV, "repository root")
        require_reexec_source_relationship(source, repository)
        with absolute_descriptor(source, "retained publisher", flags=os.O_RDONLY) as named:
            _, retained, _ = read_descriptor(named, "retained publisher", mode=0o755)
            if retained != expected:
                raise ValueError("retained publisher checksum differs after reexecution")
        return source, repository, expected
    orphaned = [
        name for name in (
            PYTHON_DESCRIPTOR_ENV, SELF_SOURCE_ENV, REPOSITORY_ROOT_ENV
        ) if name in os.environ
    ]
    if orphaned:
        raise ValueError(
            "private publisher reexecution path is forbidden without an authenticated "
            "descriptor"
        )
    source = initial_source()
    repository = initial_repository(source)
    require_reexec_source_relationship(source, repository)
    with absolute_descriptor(source, "retained publisher", flags=os.O_RDONLY) as descriptor:
        _, retained, _ = read_descriptor(descriptor, "retained publisher", mode=0o755)
        if retained != expected:
            raise ValueError("retained publisher checksum differs")
        interpreter = Path(sys.executable).resolve(strict=True)
        with absolute_descriptor(
            interpreter, "authenticated publisher Python interpreter", flags=os.O_RDONLY
        ) as python_descriptor:
            require_system_executable_profile(
                os.fstat(python_descriptor),
                "authenticated publisher Python interpreter",
            )
            if profile_identity(os.fstat(python_descriptor)) != profile_identity(
                os.stat("/proc/self/exe")
            ):
                raise ValueError("selected publisher Python interpreter is not this process")
            os.set_inheritable(descriptor, True)
            os.set_inheritable(python_descriptor, True)
            os.execve(
                f"/proc/self/fd/{python_descriptor}",
                [
                    str(interpreter),
                    "-I",
                    "-B",
                    f"/proc/self/fd/{descriptor}",
                    *argv,
                ],
                reexec_environment(
                    descriptor, python_descriptor, source, repository
                ),
            )
    raise AssertionError("publisher descriptor re-execution unexpectedly returned")


def repository_head(repository: Path) -> str:
    completed = git_run(repository, ["rev-parse", "--verify", "HEAD^{commit}"])
    if completed.returncode:
        raise ValueError("cannot resolve repository HEAD")
    try:
        return require_revision(completed.stdout.decode().strip(), "repository HEAD")
    except UnicodeDecodeError as error:
        raise ValueError("repository HEAD is not ASCII") from error


def committed_tools(repository: Path, expected_head: str,
                    publisher_sha256: str) -> list[dict[str, object]]:
    if repository_head(repository) != expected_head:
        raise ValueError("final source-authority revision is not live repository HEAD")
    retained = []
    for relative, expected_mode in sorted(REQUIRED_TOOLS.items()):
        path = repository / relative
        if git_run(repository, ["ls-files", "--error-unmatch", "--", relative]).returncode:
            raise ValueError(f"required source-authority tool is not tracked: {relative}")
        for arguments in (
            ["diff", "--quiet", "--", relative],
            ["diff", "--cached", "--quiet", "--", relative],
        ):
            if git_run(repository, arguments).returncode:
                raise ValueError(f"required source-authority tool is not committed: {relative}")
        payload, digest = read_file(
            path, f"committed tool {relative}", mode=int(expected_mode, 8)
        )
        committed = git_run(repository, ["show", f"{expected_head}:{relative}"])
        if committed.returncode or committed.stdout != payload:
            raise ValueError(f"committed tool bytes differ at final HEAD: {relative}")
        tree = git_run(repository, ["ls-tree", expected_head, "--", relative])
        expected_git_mode = "100755" if expected_mode == "0755" else "100644"
        if tree.returncode:
            raise ValueError(f"cannot inspect committed tool mode: {relative}")
        try:
            tree_line = tree.stdout.decode("ascii").strip()
        except UnicodeDecodeError as error:
            raise ValueError(f"committed tool mode is not ASCII: {relative}") from error
        if not tree_line.startswith(f"{expected_git_mode} blob ") or not tree_line.endswith(
            f"\t{relative}"
        ):
            raise ValueError(f"committed tool mode differs at final HEAD: {relative}")
        if relative == PUBLISHER_RELATIVE.as_posix() and digest != publisher_sha256:
            raise ValueError("committed publisher digest differs from authenticated self")
        retained.append(
            {
                "path": relative,
                "revision": expected_head,
                "sha256": digest,
                "mode": expected_mode,
            }
        )
    return retained


def committed_tools_at_revision(repository: Path, revision: str,
                                label: str) -> list[dict[str, object]]:
    """Authenticate required tool bytes and modes at one retained revision."""

    revision = require_revision(revision, f"{label} revision")
    retained = []
    for relative, expected_mode in sorted(REQUIRED_TOOLS.items()):
        committed = git_run(repository, ["show", f"{revision}:{relative}"])
        if committed.returncode:
            raise ValueError(f"{label} required tool is absent: {relative}")
        tree = git_run(repository, ["ls-tree", revision, "--", relative])
        expected_git_mode = "100755" if expected_mode == "0755" else "100644"
        if tree.returncode:
            raise ValueError(f"cannot inspect {label} tool mode: {relative}")
        try:
            tree_line = tree.stdout.decode("ascii").strip()
        except UnicodeDecodeError as error:
            raise ValueError(f"{label} tool mode is not ASCII: {relative}") from error
        if not tree_line.startswith(f"{expected_git_mode} blob ") or not tree_line.endswith(
            f"\t{relative}"
        ):
            raise ValueError(f"{label} tool mode differs: {relative}")
        retained.append(
            {
                "path": relative,
                "revision": revision,
                "sha256": sha256_bytes(committed.stdout),
                "mode": expected_mode,
            }
        )
    return retained


def stable_bundle_validation(repository: Path, path: Path, expected_sha256: str,
                             expected_revision: str, expected_name: str,
                             required_revisions: list[str], label: str
                             ) -> dict[str, object]:
    path = normalized_absolute(path, label)
    expected_sha256 = require_sha256(expected_sha256, f"{label} SHA-256")
    expected_revision = require_revision(expected_revision, f"{label} tip")
    revisions = [require_revision(value, f"{label} required revision") for value in required_revisions]
    if not revisions or len(set(revisions)) != len(revisions):
        raise ValueError(f"{label} required revisions are empty or duplicated")
    require_no_symlinks(path, label)
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        payload, observed, descriptor_before = read_descriptor(
            descriptor, label, expected=expected_sha256, mode=0o644
        )
        named_before = os.stat(path, follow_symlinks=False)
        if (
            named_before.st_dev,
            named_before.st_ino,
        ) != (
            descriptor_before.st_dev,
            descriptor_before.st_ino,
        ):
            raise ValueError(f"{label} pathname changed before validation")
        header = payload.split(b"\n\n", 1)[0]
        try:
            lines = header.decode("utf-8").splitlines()
        except UnicodeDecodeError as error:
            raise ValueError(f"{label} header is not UTF-8") from error
        if not lines or lines[0] not in {"# v2 git bundle", "# v3 git bundle"}:
            raise ValueError(f"{label} header version is invalid")
        if any(line.startswith("-") for line in lines[1:]):
            raise ValueError(f"{label} must be self-contained without prerequisites")
        advertised = []
        for line in lines[1:]:
            if line.startswith("@"):
                continue
            match = re.fullmatch(r"([0-9a-f]{40}) (.+)", line)
            if match is None:
                raise ValueError(f"{label} advertised reference is malformed")
            advertised.append((match.group(1), match.group(2)))
        if advertised != [(expected_revision, expected_name)]:
            raise ValueError(f"{label} must advertise exactly {expected_revision} {expected_name}")

        # Every Git operation opens this inherited descriptor path, so it
        # authenticates the exact inode already checksummed above rather than
        # a pathname that can be exchanged between opens.
        descriptor_path = f"/proc/self/fd/{descriptor}"
        verify = git_run(repository, ["bundle", "verify", descriptor_path], pass_fds=(descriptor,))
        if verify.returncode:
            raise ValueError(f"{label} fails git bundle verify")
        heads = git_run(repository, ["bundle", "list-heads", descriptor_path],
                        pass_fds=(descriptor,))
        if heads.returncode:
            raise ValueError(f"{label} heads cannot be listed")
        if heads.stdout.decode("utf-8").splitlines() != [
            f"{expected_revision} {expected_name}"
        ]:
            raise ValueError(f"{label} Git-advertised head differs")
        with tempfile.TemporaryDirectory(prefix="cgl-lf-source-authority-bundle-") as directory:
            isolated = Path(directory) / "repository.git"
            if git_run(repository, ["init", "--bare", str(isolated)]).returncode:
                raise ValueError(f"{label} isolated repository initialization failed")
            unbundle = git_run(
                isolated, ["bundle", "unbundle", descriptor_path], pass_fds=(descriptor,)
            )
            if unbundle.returncode:
                raise ValueError(f"{label} cannot be reconstructed in isolation")
            if git_run(
                isolated,
                ["fsck", "--full", "--strict", "--no-reflogs", expected_revision],
            ).returncode:
                raise ValueError(f"{label} advertised history is incomplete")
            for revision in revisions:
                if git_run(isolated, ["cat-file", "-e", f"{revision}^{{commit}}"]).returncode:
                    raise ValueError(f"{label} omits required revision {revision}")
                if git_run(
                    isolated, ["merge-base", "--is-ancestor", revision, expected_revision]
                ).returncode:
                    raise ValueError(f"{label} does not cover required revision {revision}")
        retained_payload, retained_sha256, descriptor_after = read_descriptor(
            descriptor, label, expected=expected_sha256, mode=0o644
        )
        if retained_sha256 != observed or retained_payload != payload:
            raise ValueError(f"{label} bytes changed during validation")
        require_no_symlinks(path, label)
        named_after = os.stat(path, follow_symlinks=False)
        if (
            named_after.st_dev,
            named_after.st_ino,
        ) != (
            descriptor_after.st_dev,
            descriptor_after.st_ino,
        ):
            raise ValueError(f"{label} pathname changed during validation")
        fields = (
            "st_dev", "st_ino", "st_mode", "st_uid", "st_nlink", "st_size",
            "st_mtime_ns", "st_ctime_ns",
        )
        if any(
            getattr(descriptor_before, field) != getattr(descriptor_after, field)
            for field in fields
        ):
            raise ValueError(f"{label} descriptor changed during validation")
    finally:
        os.close(descriptor)
    return {
        "sha256": observed,
        "complete_history": True,
        "advertised_tip": {"revision": expected_revision, "name": expected_name},
        "required_revisions": revisions,
    }


def binding(value: object, label: str) -> tuple[PurePosixPath, str]:
    retained = require_exact_keys(value, {"path", "sha256"}, label)
    return require_relative(retained["path"], f"{label} path"), require_sha256(
        retained["sha256"], f"{label} SHA-256"
    )


def declared_binding(value: object, path: Path, digest: str, label: str,
                     *, mode: str = "0444") -> None:
    expected = {"path": str(path), "sha256": digest, "mode": mode, "links": 1}
    if value != expected:
        raise ValueError(f"{label} differs from exact published binding")


def historical_f115(root: Path, repository: Path, bindings: object,
                    *, canonical: bool) -> dict[str, object]:
    retained_bindings = require_exact_keys(
        bindings, set(F115_PATHS), "historical F115 bindings"
    )
    loaded: dict[str, tuple[dict[str, object], str, Path]] = {}
    for key, relative in F115_PATHS.items():
        selected_relative, expected = binding(
            retained_bindings[key], f"historical F115 {key}"
        )
        if selected_relative != relative:
            raise ValueError(f"historical F115 {key} path differs")
        if canonical and expected != F115_CANONICAL_SHA256[key]:
            raise ValueError(f"historical F115 {key} digest differs from canonical authority")
        path = root / relative
        value, _, digest = read_json(
            path, f"historical F115 {key}", expected=expected, mode=0o444
        )
        loaded[key] = value, digest, path
    evidence, evidence_sha, evidence_path = loaded["evidence"]
    if (
        evidence.get("schema_version") != 1
        or evidence.get("record_type")
        != "stage-i-source-bundle-recovery-supersession-evidence"
        or evidence.get("checkpoint") != "F-115"
        or evidence.get("execution_epoch") != EXECUTION_EPOCH
    ):
        raise ValueError("historical F115 evidence identity differs")
    implementation = evidence.get("implementation")
    if not isinstance(implementation, dict) or not isinstance(
        implementation.get("source_bundle"), dict
    ):
        raise ValueError("historical F115 source-bundle implementation is missing")
    bundle = implementation["source_bundle"]
    bundle_relative = require_relative(bundle.get("path"), "historical F115 bundle path")
    bundle_sha = require_sha256(bundle.get("sha256"), "historical F115 bundle SHA-256")
    bundle_head = require_revision(bundle.get("head"), "historical F115 bundle head")
    revisions = bundle.get("verified_revisions")
    if (
        bundle.get("complete_history") is not True
        or not isinstance(revisions, list)
        or not revisions
    ):
        raise ValueError("historical F115 bundle declaration differs")
    f115_bundle = root / bundle_relative
    stable_bundle_validation(
        repository,
        f115_bundle,
        bundle_sha,
        bundle_head,
        "HEAD",
        revisions,
        "historical F115 source bundle",
    )
    audit, _, _ = loaded["publication_audit"]
    if (
        audit.get("schema_version") != 1
        or audit.get("record_type")
        != "stage-i-source-bundle-recovery-supersession-publication-audit"
        or audit.get("checkpoint") != "F-115"
        or audit.get("execution_epoch") != EXECUTION_EPOCH
    ):
        raise ValueError("historical F115 publication audit identity differs")
    declared_binding(audit.get("artifact"), evidence_path, evidence_sha, "historical F115 artifact")
    reviews = audit.get("independent_reviews")
    if not isinstance(reviews, dict) or reviews.get(
        "reviews_bind_exact_published_f115_sha256"
    ) != evidence_sha:
        raise ValueError("historical F115 review binding differs")
    agents = set()
    for key, audit_key, kind, decision in (
        ("provenance_review", "provenance_security", "provenance-security",
         "approved-for-publication"),
        ("plasma_review", "plasma_scientific_continuation",
         "plasma-scientific-continuation", "approved"),
    ):
        review, review_sha, review_path = loaded[key]
        declared_binding(reviews.get(audit_key), review_path, review_sha, f"historical F115 {key}")
        published = review.get("published_f115")
        reviewer = review.get("reviewer")
        if (
            review.get("schema_version") != 1
            or review.get("record_type")
            != "stage-i-source-bundle-recovery-supersession-independent-review"
            or review.get("checkpoint") != "F-115"
            or review.get("execution_epoch") != EXECUTION_EPOCH
            or review.get("review_kind") != kind
            or review.get("decision") != decision
            or published != {"path": str(evidence_path), "sha256": evidence_sha}
            or not isinstance(reviewer, dict)
        ):
            raise ValueError(f"historical F115 {key} identity differs")
        agent = require_nonempty(reviewer.get("agent_id"), f"historical F115 {key} reviewer")
        if agent in agents:
            raise ValueError("historical F115 reviewers are not independent")
        agents.add(agent)
    authority = audit.get("authority_and_enforcement")
    if not isinstance(authority, dict) or authority.get("direct_sbatch_authorized") is not False:
        raise ValueError("historical F115 publication audit over-authorizes")
    return {
        "digests": {
            "evidence_sha256": evidence_sha,
            "publication_audit_sha256": loaded["publication_audit"][1],
            "provenance_review_sha256": loaded["provenance_review"][1],
            "plasma_review_sha256": loaded["plasma_review"][1],
        },
        "bundle": {
            "path": bundle_relative.as_posix(),
            "sha256": bundle_sha,
            "head": bundle_head,
            "verified_revisions": revisions,
        },
    }


def historical_f116(root: Path, repository: Path, bindings: object,
                    *, canonical: bool) -> dict[str, object]:
    """Authenticate the immutable four-part F116 authority and nested F115."""

    retained_bindings = require_exact_keys(
        bindings, set(F116_PATHS), "historical F116 bindings"
    )
    loaded: dict[str, tuple[dict[str, object], str, Path]] = {}
    for key, relative in F116_PATHS.items():
        selected_relative, expected = binding(
            retained_bindings[key], f"historical F116 {key}"
        )
        if selected_relative != relative:
            raise ValueError(f"historical F116 {key} path differs")
        if canonical and expected != F116_CANONICAL_SHA256[key]:
            raise ValueError(f"historical F116 {key} digest differs from canonical authority")
        path = root / relative
        value, _, digest = read_json(
            path, f"historical F116 {key}", expected=expected, mode=0o444
        )
        loaded[key] = value, digest, path

    evidence, evidence_sha, evidence_path = loaded["evidence"]
    evidence = require_exact_keys(
        evidence,
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "generated_utc", "scope", "predecessor_authorities", "implementation",
            "source_archive_catalog", "authorization", "validation",
            "publication_requirements",
        },
        "historical F116 evidence",
    )
    if (
        evidence["schema_version"] != 1
        or evidence["record_type"] != "stage-i-current-source-authority-supersession-evidence"
        or evidence["checkpoint"] != "F-116"
        or evidence["execution_epoch"] != EXECUTION_EPOCH
        or evidence["authorization"] != AUTHORIZATION
        or evidence["validation"] != F116_VALIDATION_CLAIMS
        or evidence["publication_requirements"] != PUBLICATION_REQUIREMENTS
    ):
        raise ValueError("historical F116 evidence identity or authority differs")
    generated = require_utc(evidence["generated_utc"], "historical F116 generation timestamp")
    scope = require_exact_keys(
        evidence["scope"],
        {"relationship", "summary", "preserves", "does_not_authorize"},
        "historical F116 scope",
    )
    if (
        scope["relationship"] != "current-source-selection-only-supersession"
        or not require_nonempty(scope["summary"], "historical F116 scope summary")
        or scope["preserves"] != F116_PRESERVES
        or scope["does_not_authorize"] != DOES_NOT_AUTHORIZE
    ):
        raise ValueError("historical F116 scope differs or broadens authority")
    predecessors = require_exact_keys(
        evidence["predecessor_authorities"], {"historical_f115"},
        "historical F116 predecessor authorities",
    )
    f115 = historical_f115(
        root, repository, predecessors["historical_f115"], canonical=canonical
    )
    implementation = require_exact_keys(
        evidence["implementation"],
        {
            "publisher", "committed_tools", "intermediate_36140_bundle",
            "current_source_bundle",
        },
        "historical F116 implementation",
    )
    bridge = parse_bundle_declaration(
        implementation["intermediate_36140_bundle"],
        "historical F116 bridge bundle",
        current=False,
    )
    current = parse_bundle_declaration(
        implementation["current_source_bundle"],
        "historical F116 current source bundle",
        current=True,
    )
    stable_bundle_validation(
        repository,
        root / str(bridge["path"]),
        str(bridge["sha256"]),
        str(bridge["head"]),
        str(bridge["advertised_tip"]["name"]),
        list(bridge["verified_revisions"]),
        "historical F116 bridge bundle",
    )
    stable_bundle_validation(
        repository,
        root / str(current["path"]),
        str(current["sha256"]),
        str(current["head"]),
        "HEAD",
        list(current["verified_revisions"]),
        "historical F116 current source bundle",
    )
    if (
        f115["bundle"]["head"] not in set(bridge["verified_revisions"])
        or not {f115["bundle"]["head"], bridge["head"], current["head"]}.issubset(
            set(current["verified_revisions"])
        )
    ):
        raise ValueError("historical F116 bundles do not cover the nested F115 chain")
    tools = committed_tools_at_revision(
        repository, str(current["head"]), "historical F116"
    )
    publisher = next(
        item for item in tools if item["path"] == PUBLISHER_RELATIVE.as_posix()
    )
    if implementation["committed_tools"] != tools or implementation["publisher"] != publisher:
        raise ValueError("historical F116 committed-tool vector or publisher differs")
    if current["subject"] != revision_subject(
        repository, str(current["head"]), "historical F116 current source"
    ):
        raise ValueError("historical F116 current source subject differs")

    catalogs = require_exact_keys(
        evidence["source_archive_catalog"], {"before", "after"},
        "historical F116 source-archive catalog",
    )
    before = require_exact_keys(
        catalogs["before"],
        {
            "readme_sha256", "sha256sums_sha256", "bridge_listed",
            "final_bundle_listed", "corrupt_c7_listed",
        },
        "historical F116 source-archive catalog before",
    )
    after = require_exact_keys(
        catalogs["after"],
        {
            "readme_sha256", "sha256sums_sha256", "bridge_listed_exactly_once",
            "final_bundle_listed_exactly_once", "corrupt_c7_listed",
            "historical_f115_preserved", "sole_current_source_bundle",
        },
        "historical F116 source-archive catalog after",
    )
    for catalog, label in ((before, "before"), (after, "after")):
        for key in ("readme_sha256", "sha256sums_sha256"):
            require_sha256(catalog[key], f"historical F116 catalog {label} {key}")
    if (
        before["bridge_listed"] is not False
        or before["final_bundle_listed"] is not False
        or before["corrupt_c7_listed"] is not False
        or after["bridge_listed_exactly_once"] is not True
        or after["final_bundle_listed_exactly_once"] is not True
        or after["corrupt_c7_listed"] is not False
        or after["historical_f115_preserved"] is not True
        or after["sole_current_source_bundle"] != current["path"]
    ):
        raise ValueError("historical F116 source-archive catalog policy differs")

    expected_verified = {
        "authorization_broadening": False,
        "bridge_selected_as_current": False,
        "corrupt_c7_excluded": True,
        "current_source_selection_only": True,
        "final_bundle_sha256": current["sha256"],
        "final_head": current["head"],
        "historical_f115_preserved": True,
    }
    review_times = []
    agents = set()
    for key, kind, decision in (
        ("provenance_review", "provenance-security", "approved-for-publication"),
        ("plasma_review", "plasma-scientific-continuation", "approved"),
    ):
        review, _, review_path = loaded[key]
        review = require_exact_keys(
            review,
            {
                "schema_version", "record_type", "checkpoint", "execution_epoch",
                "review_kind", "decision", "reviewed_candidate", "published_f116",
                "reviewer", "reviewed_utc", "findings", "limitations", "verified",
            },
            f"historical F116 {key}",
        )
        reviewed_candidate = require_exact_keys(
            review["reviewed_candidate"], {"path", "sha256"},
            f"historical F116 {key} reviewed candidate",
        )
        normalized_absolute(
            Path(require_nonempty(
                reviewed_candidate["path"], f"historical F116 {key} reviewed candidate path"
            )),
            f"historical F116 {key} reviewed candidate path",
        )
        reviewer = require_exact_keys(
            review["reviewer"], {"agent_id", "identity"},
            f"historical F116 {key} reviewer",
        )
        agent = require_nonempty(
            reviewer["agent_id"], f"historical F116 {key} reviewer agent ID"
        )
        if (
            review["schema_version"] != 1
            or review["record_type"]
            != "stage-i-current-source-authority-supersession-independent-review"
            or review["checkpoint"] != "F-116"
            or review["execution_epoch"] != EXECUTION_EPOCH
            or review["review_kind"] != kind
            or review["decision"] != decision
            or reviewed_candidate["sha256"] != evidence_sha
            or review["published_f116"]
            != {"path": str(evidence_path), "sha256": evidence_sha}
            or review["verified"] != expected_verified
            or agent in agents
            or not require_nonempty(
                reviewer["identity"], f"historical F116 {key} reviewer identity"
            )
            or not isinstance(review["findings"], list)
            or not review["findings"]
            or any(not isinstance(item, str) or not item for item in review["findings"])
            or not isinstance(review["limitations"], list)
            or not review["limitations"]
            or any(not isinstance(item, str) or not item for item in review["limitations"])
            or INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION not in review["limitations"]
        ):
            raise ValueError(f"historical F116 {key} identity or review differs")
        reviewed = require_utc(review["reviewed_utc"], f"historical F116 {key} timestamp")
        if reviewed < generated:
            raise ValueError(f"historical F116 {key} predates evidence")
        agents.add(agent)
        review_times.append(reviewed)

    audit, _, audit_path = loaded["publication_audit"]
    audit = require_exact_keys(
        audit,
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "published_utc", "artifact", "independent_reviews",
            "historical_f115_authority", "source_archive_catalog",
            "authority_and_enforcement", "publication",
        },
        "historical F116 publication audit",
    )
    if (
        audit["schema_version"] != 1
        or audit["record_type"]
        != "stage-i-current-source-authority-supersession-publication-audit"
        or audit["checkpoint"] != "F-116"
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or audit["historical_f115_authority"] != f115["digests"]
        or audit["authority_and_enforcement"] != AUTHORIZATION
        or audit["publication"]
        != "recoverable-forward-transaction-with-publication-audit-commit-marker-under-stage-i-lock"
    ):
        raise ValueError("historical F116 publication audit identity or authority differs")
    published = require_utc(audit["published_utc"], "historical F116 publication timestamp")
    if any(published < reviewed for reviewed in review_times):
        raise ValueError("historical F116 publication audit predates review")
    declared_binding(
        audit["artifact"], evidence_path, evidence_sha, "historical F116 audit artifact"
    )
    reviews = require_exact_keys(
        audit["independent_reviews"],
        {
            "reviews_bind_exact_published_f116_sha256",
            "provenance_security", "plasma_scientific_continuation",
        },
        "historical F116 audit independent reviews",
    )
    if reviews["reviews_bind_exact_published_f116_sha256"] != evidence_sha:
        raise ValueError("historical F116 audit review digest binding differs")
    declared_binding(
        reviews["provenance_security"],
        loaded["provenance_review"][2],
        loaded["provenance_review"][1],
        "historical F116 provenance review",
    )
    declared_binding(
        reviews["plasma_scientific_continuation"],
        loaded["plasma_review"][2],
        loaded["plasma_review"][1],
        "historical F116 plasma review",
    )
    audit_catalog = require_exact_keys(
        audit["source_archive_catalog"],
        {
            "readme", "sha256sums", "bridge_bundle", "current_source_bundle",
            "corrupt_c7_absent_from_active_checksum_ledger",
            "sole_current_source_bundle",
        },
        "historical F116 audit source-archive catalog",
    )
    declared_binding(
        audit_catalog["readme"], root / "source-archives/README.md",
        str(after["readme_sha256"]), "historical F116 catalog README", mode="0644",
    )
    declared_binding(
        audit_catalog["sha256sums"], root / "source-archives/SHA256SUMS",
        str(after["sha256sums_sha256"]), "historical F116 catalog SHA256SUMS", mode="0644",
    )
    if audit_catalog["bridge_bundle"] != {
        "path": str(root / str(bridge["path"])),
        "sha256": bridge["sha256"],
        "mode": "0644",
        "links": 1,
        "head": bridge["head"],
        "role": "retained-non-current-bridge",
        "selected_as_current": False,
    } or audit_catalog["current_source_bundle"] != {
        "path": str(root / str(current["path"])),
        "sha256": current["sha256"],
        "mode": "0644",
        "links": 1,
        "head": current["head"],
        "selected_as_current": True,
    } or (
        audit_catalog["corrupt_c7_absent_from_active_checksum_ledger"] is not True
        or audit_catalog["sole_current_source_bundle"] != str(root / str(current["path"]))
    ):
        raise ValueError("historical F116 audit catalog binding differs")

    return {
        "digests": {
            "evidence_sha256": evidence_sha,
            "publication_audit_sha256": loaded["publication_audit"][1],
            "provenance_review_sha256": loaded["provenance_review"][1],
            "plasma_review_sha256": loaded["plasma_review"][1],
        },
        "historical_f115": f115,
        "bridge": bridge,
        "bundle": current,
        "catalog_after": after,
        "publication_audit_path": audit_path,
    }


def parse_sha256sums(payload: bytes, label: str) -> list[tuple[str, str]]:
    try:
        lines = payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} is not UTF-8") from error
    entries = []
    names = set()
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


def validate_checksum_entries(source_archives: Path, entries: list[tuple[str, str]],
                              *, final_name: str | None = None,
                              final_payload: bytes | None = None) -> None:
    for digest, name in entries:
        if name == final_name and final_payload is not None:
            if sha256_bytes(final_payload) != digest:
                raise ValueError("final bundle checksum-ledger entry differs")
            continue
        path = source_archives / name
        _, observed = read_file(path, f"source archive {name}", expected=digest, mode=0o644)
        if observed != digest:
            raise ValueError(f"source archive checksum differs: {name}")


def catalog_readme_block(*, final_name: str, final_sha256: str,
                         final_revision: str, final_subject: str) -> bytes:
    return (
        f"`{final_name}` records complete history through commit `{final_revision}` "
        f"(`{final_subject}`). It is the sole current Stage I source-selection bundle "
        "after independently reviewed F-118 publication. It preserves the immutable "
        "four-part F-116 authority and every prior active source-archive catalog entry; "
        "it does not itself authorize prepare or submission. "
        f"Its SHA-256 is `{final_sha256}`.\n\n"
    ).encode()


def catalog_payloads(old_readme: bytes, old_sums: bytes, *,
                     bridge_name: str, predecessor_name: str,
                     final_name: str, final_sha256: str, final_revision: str,
                     final_subject: str) -> tuple[bytes, bytes]:
    marker = b"## AthenaK\n\n"
    if old_readme.count(marker) != 1:
        raise ValueError("source-archive README lacks one exact AthenaK insertion marker")
    for name in (bridge_name, predecessor_name):
        if name.encode() not in old_readme:
            raise ValueError(f"source-archive README does not preserve {name}")
    if final_name.encode() in old_readme:
        raise ValueError(f"source-archive README already mentions {final_name}")
    block = catalog_readme_block(
        final_name=final_name,
        final_sha256=final_sha256,
        final_revision=final_revision,
        final_subject=final_subject,
    )
    new_readme = old_readme.replace(marker, marker + block, 1)
    entries = parse_sha256sums(old_sums, "source-archive SHA256SUMS")
    names = {name for _, name in entries}
    if bridge_name not in names or predecessor_name not in names:
        raise ValueError("source-archive checksum ledger omits an F118 predecessor bundle")
    if final_name in names:
        raise ValueError("source-archive checksum ledger already lists the F118 final bundle")
    new_sums = old_sums
    if new_sums and not new_sums.endswith(b"\n"):
        raise ValueError("source-archive SHA256SUMS must end newline")
    new_sums += f"{final_sha256}  {final_name}\n".encode()
    return new_readme, new_sums


def validate_catalog_before(root: Path, f116: dict[str, object],
                            bridge: dict[str, object], predecessor: dict[str, object],
                            final: dict[str, object], final_payload: bytes
                            ) -> dict[str, object]:
    source_archives = root / "source-archives"
    readme_path = source_archives / "README.md"
    sums_path = source_archives / "SHA256SUMS"
    old_readme, old_readme_sha = read_file(readme_path, "source-archive README", mode=0o644)
    old_sums, old_sums_sha = read_file(sums_path, "source-archive SHA256SUMS", mode=0o644)
    f116_catalog_after = require_exact_keys(
        f116["catalog_after"],
        {
            "readme_sha256", "sha256sums_sha256", "bridge_listed_exactly_once",
            "final_bundle_listed_exactly_once", "corrupt_c7_listed",
            "historical_f115_preserved", "sole_current_source_bundle",
        },
        "authenticated F116 source-archive catalog after",
    )
    if (
        old_readme_sha != f116_catalog_after["readme_sha256"]
        or old_sums_sha != f116_catalog_after["sha256sums_sha256"]
    ):
        raise ValueError(
            "live predecessor source-archive catalog differs from authenticated "
            "F116 catalog_after"
        )
    entries = parse_sha256sums(old_sums, "source-archive SHA256SUMS")
    validate_checksum_entries(source_archives, entries)
    names = [name for _, name in entries]
    f115_name = Path(str(f116["historical_f115"]["bundle"]["path"])).name
    predecessor_name = str(predecessor["name"])
    if names.count(f115_name) != 1:
        raise ValueError("historical F115 bundle is not listed exactly once")
    if names.count(str(bridge["name"])) != 1:
        raise ValueError("historical bridge bundle is not listed exactly once")
    if names.count(predecessor_name) != 1:
        raise ValueError("F116 predecessor bundle is not listed exactly once")
    if CORRUPT_C7_NAME in names:
        raise ValueError("corrupt C7 bundle reappeared in active checksum ledger")
    if names.count(str(final["name"])):
        raise ValueError("final F118 bundle is already cataloged")
    new_readme, new_sums = catalog_payloads(
        old_readme,
        old_sums,
        bridge_name=str(bridge["name"]),
        predecessor_name=predecessor_name,
        final_name=str(final["name"]),
        final_sha256=str(final["sha256"]),
        final_revision=str(final["head"]),
        final_subject=str(final["subject"]),
    )
    new_entries = parse_sha256sums(new_sums, "new source-archive SHA256SUMS")
    validate_checksum_entries(
        source_archives, new_entries, final_name=str(final["name"]), final_payload=final_payload
    )
    new_names = [name for _, name in new_entries]
    if (
        new_names.count(str(bridge["name"])) != 1
        or new_names.count(predecessor_name) != 1
        or new_names.count(str(final["name"])) != 1
        or CORRUPT_C7_NAME in new_names
        or new_names.count(f115_name) != 1
        or new_entries[:-1] != entries
    ):
        raise ValueError("new source-archive checksum ledger violates F-118 catalog policy")
    return {
        "old_readme": old_readme,
        "old_sums": old_sums,
        "new_readme": new_readme,
        "new_sums": new_sums,
        "before": {
            "readme_sha256": old_readme_sha,
            "sha256sums_sha256": old_sums_sha,
            "bridge_listed_exactly_once": True,
            "predecessor_current_source_bundle_listed_exactly_once": True,
            "final_bundle_listed": False,
            "corrupt_c7_listed": False,
            "historical_f115_preserved": True,
        },
        "after": {
            "readme_sha256": sha256_bytes(new_readme),
            "sha256sums_sha256": sha256_bytes(new_sums),
            "bridge_listed_exactly_once": True,
            "predecessor_current_source_bundle_listed_exactly_once": True,
            "final_bundle_listed_exactly_once": True,
            "corrupt_c7_listed": False,
            "historical_f115_preserved": True,
            "historical_f116_preserved": True,
            "all_prior_checksum_entries_preserved": True,
            "sole_current_source_bundle": str(final["path"]),
        },
    }


def parse_bundle_declaration(value: object, label: str, *,
                             current: bool) -> dict[str, object]:
    keys = {
        "path", "sha256", "complete_history", "head", "advertised_tip",
        "verified_revisions", "selected_as_current",
    }
    if current:
        keys |= {"candidate_path", "subject"}
    else:
        keys |= {"role"}
    retained = require_exact_keys(value, keys, label)
    relative = require_relative(retained["path"], f"{label} path")
    if relative.parent != PurePosixPath("source-archives"):
        raise ValueError(f"{label} must be a direct child of source-archives")
    name_match = BUNDLE_NAME_RE.fullmatch(relative.name)
    head = require_revision(retained["head"], f"{label} head")
    if name_match is None or name_match.group(1) != head[:9]:
        raise ValueError(f"{label} filename does not bind its head")
    revisions = retained["verified_revisions"]
    if not isinstance(revisions, list) or not revisions:
        raise ValueError(f"{label} verified revisions must be a nonempty list")
    revisions = [require_revision(item, f"{label} verified revision") for item in revisions]
    if len(set(revisions)) != len(revisions):
        raise ValueError(f"{label} verified revisions are duplicated")
    tip = require_exact_keys(
        retained["advertised_tip"], {"revision", "name"}, f"{label} advertised tip"
    )
    if tip["revision"] != head:
        raise ValueError(f"{label} advertised revision differs from head")
    if retained["complete_history"] is not True:
        raise ValueError(f"{label} is not declared complete history")
    if current:
        candidate_path = normalized_absolute(
            Path(require_nonempty(retained["candidate_path"], f"{label} candidate path")),
            f"{label} candidate path",
        )
        if tip["name"] != "HEAD" or retained["selected_as_current"] is not True:
            raise ValueError(f"{label} is not the sole exact HEAD current bundle")
        subject = require_nonempty(retained["subject"], f"{label} subject")
    else:
        candidate_path = None
        subject = None
        if (
            tip["name"] != "refs/heads/feature/cgl-landau-fluid"
            or retained["selected_as_current"] is not False
            or retained["role"] != "retained-non-current-bridge"
        ):
            raise ValueError(f"{label} role or branch-ref declaration differs")
    return {
        **retained,
        "path": relative.as_posix(),
        "name": relative.name,
        "sha256": require_sha256(retained["sha256"], f"{label} SHA-256"),
        "head": head,
        "verified_revisions": revisions,
        "candidate_path": candidate_path,
        "subject": subject,
    }


def parse_predecessor_bundle_declaration(value: object, label: str) -> dict[str, object]:
    retained = require_exact_keys(
        value,
        {
            "path", "sha256", "complete_history", "head", "advertised_tip",
            "verified_revisions", "selected_as_current", "role", "subject",
        },
        label,
    )
    selected = {
        **retained,
        "candidate_path": "/retained-predecessor-not-a-candidate",
        "selected_as_current": True,
    }
    selected.pop("role")
    parsed = parse_bundle_declaration(selected, label, current=True)
    if (
        retained["selected_as_current"] is not False
        or retained["role"] != "retained-non-current-predecessor"
    ):
        raise ValueError(f"{label} role or selection declaration differs")
    return {
        **retained,
        "path": parsed["path"],
        "name": parsed["name"],
        "sha256": parsed["sha256"],
        "head": parsed["head"],
        "verified_revisions": parsed["verified_revisions"],
        "subject": parsed["subject"],
    }


def retained_predecessor_declaration(f116: dict[str, object]) -> dict[str, object]:
    bundle = f116["bundle"]
    if not isinstance(bundle, dict):
        raise ValueError("historical F116 current source bundle is missing")
    return {
        "path": bundle["path"],
        "sha256": bundle["sha256"],
        "complete_history": True,
        "head": bundle["head"],
        "advertised_tip": bundle["advertised_tip"],
        "verified_revisions": bundle["verified_revisions"],
        "selected_as_current": False,
        "role": "retained-non-current-predecessor",
        "subject": bundle["subject"],
    }


def parse_evidence(value: dict[str, object], *, root: Path, repository: Path,
                   evidence_candidate_path: Path, bundle_candidate_path: Path,
                   publisher_sha256: str, canonical: bool,
                   old_catalog: dict[str, object], new_catalog: dict[str, object]
                   ) -> dict[str, object]:
    evidence = require_exact_keys(
        value,
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "generated_utc", "scope", "predecessor_authorities", "implementation",
            "source_archive_catalog", "authorization", "validation",
            "publication_requirements",
        },
        "F118 evidence",
    )
    if (
        evidence["schema_version"] != 1
        or evidence["record_type"] != "stage-i-current-source-authority-supersession-evidence"
        or evidence["checkpoint"] != CHECKPOINT
        or evidence["execution_epoch"] != EXECUTION_EPOCH
    ):
        raise ValueError("F118 evidence identity differs")
    require_utc(evidence["generated_utc"], "F118 evidence generation timestamp")
    scope = require_exact_keys(
        evidence["scope"],
        {"relationship", "summary", "preserves", "does_not_authorize"},
        "F118 scope",
    )
    if (
        scope["relationship"] != "current-source-selection-only-supersession"
        or not require_nonempty(scope["summary"], "F118 scope summary")
        or scope["preserves"] != PRESERVES
        or scope["does_not_authorize"] != DOES_NOT_AUTHORIZE
    ):
        raise ValueError("F118 scope differs or broadens authority")
    predecessors = require_exact_keys(
        evidence["predecessor_authorities"], {"historical_f116"},
        "F118 predecessor authorities",
    )
    implementation = require_exact_keys(
        evidence["implementation"],
        {
            "publisher", "committed_tools", "intermediate_36140_bundle",
            "predecessor_current_source_bundle", "current_source_bundle",
        },
        "F118 implementation",
    )
    bridge = parse_bundle_declaration(
        implementation["intermediate_36140_bundle"], "F118 bridge bundle", current=False
    )
    final = parse_bundle_declaration(
        implementation["current_source_bundle"], "F118 current source bundle", current=True
    )
    predecessor = parse_predecessor_bundle_declaration(
        implementation["predecessor_current_source_bundle"],
        "F118 predecessor current source bundle",
    )
    if final["candidate_path"] != bundle_candidate_path:
        raise ValueError("F118 final bundle candidate path differs")
    if canonical and (
        bridge["head"] != BRIDGE_REVISION
        or bridge["sha256"] != BRIDGE_SHA256
        or bridge["name"] != BRIDGE_NAME
    ):
        raise ValueError("canonical F118 bridge identity differs")
    head = final["head"]
    tools = committed_tools(repository, head, publisher_sha256)
    if implementation["committed_tools"] != tools:
        raise ValueError("F118 committed-tool vector differs")
    publisher = {
        "path": PUBLISHER_RELATIVE.as_posix(),
        "revision": head,
        "sha256": publisher_sha256,
        "mode": REQUIRED_TOOLS[PUBLISHER_RELATIVE.as_posix()],
    }
    if implementation["publisher"] != publisher:
        raise ValueError("F118 publisher binding differs")
    subject = git_run(repository, ["show", "-s", "--format=%s", head])
    try:
        committed_subject = subject.stdout.decode("utf-8").rstrip("\n")
    except UnicodeDecodeError as error:
        raise ValueError("F118 final HEAD subject is not UTF-8") from error
    if subject.returncode or not committed_subject or final["subject"] != committed_subject:
        raise ValueError("F118 final HEAD subject differs")
    required = {head, bridge["head"], predecessor["head"]}
    if canonical:
        required |= PRODUCTION_REQUIRED_REVISIONS
    if not required.issubset(set(final["verified_revisions"])):
        raise ValueError("F118 final bundle verified revisions omit required history")
    if evidence["source_archive_catalog"] != {
        "before": old_catalog,
        "after": new_catalog,
    }:
        raise ValueError("F118 source-archive catalog binding differs")
    if evidence["authorization"] != AUTHORIZATION:
        raise ValueError("F118 evidence broadens source-selection-only authority")
    if evidence["validation"] != VALIDATION_CLAIMS:
        raise ValueError("F118 validation claims differ")
    if evidence["publication_requirements"] != PUBLICATION_REQUIREMENTS:
        raise ValueError("F118 publication requirements differ")
    return {
        "historical_f116_bindings": predecessors["historical_f116"],
        "bridge": bridge,
        "predecessor": predecessor,
        "final": final,
        "head": head,
        "tools": tools,
        "generated_utc": evidence["generated_utc"],
        "evidence_candidate_path": evidence_candidate_path,
    }


def review_verified(final: dict[str, object]) -> dict[str, object]:
    return {
        "authorization_broadening": False,
        "bridge_selected_as_current": False,
        "predecessor_current_source_bundle_selected_as_current": False,
        "corrupt_c7_excluded": True,
        "current_source_selection_only": True,
        "final_bundle_sha256": final["sha256"],
        "final_head": final["head"],
        "historical_f115_preserved": True,
        "historical_f116_preserved": True,
    }


def declared_process_independence_assurance(
    role_agents: dict[str, str], label: str
) -> dict[str, object]:
    """Represent strict declared process independence without identity overclaim."""

    if len(role_agents) < 2:
        raise ValueError(f"{label} must declare at least two distinct process roles")
    retained = {
        require_nonempty(role, f"{label} role"): require_nonempty(
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


def parse_reviews(provenance: dict[str, object], plasma: dict[str, object], *,
                  evidence_candidate_path: Path, evidence_sha256: str,
                  evidence_path: Path, final: dict[str, object],
                  generated_utc: str) -> tuple[str, str]:
    agents = set()
    role_agents = {}
    reviewed_times = []
    for value, label, kind, decision in (
        (provenance, "F118 provenance review", "provenance-security",
         "approved-for-publication"),
        (plasma, "F118 plasma review", "plasma-scientific-continuation", "approved"),
    ):
        review = require_exact_keys(
            value,
            {
                "schema_version", "record_type", "checkpoint", "execution_epoch",
                "review_kind", "decision", "reviewed_candidate", "published_f118",
                "reviewer", "reviewed_utc", "findings", "limitations", "verified",
            },
            label,
        )
        if (
            review["schema_version"] != 1
            or review["record_type"]
            != "stage-i-current-source-authority-supersession-independent-review"
            or review["checkpoint"] != CHECKPOINT
            or review["execution_epoch"] != EXECUTION_EPOCH
            or review["review_kind"] != kind
            or review["decision"] != decision
            or review["reviewed_candidate"]
            != {"path": str(evidence_candidate_path), "sha256": evidence_sha256}
            or review["published_f118"]
            != {"path": str(evidence_path), "sha256": evidence_sha256}
            or review["verified"] != review_verified(final)
        ):
            raise ValueError(f"{label} identity or authorization boundary differs")
        reviewer = require_exact_keys(review["reviewer"], {"agent_id", "identity"},
                                      f"{label} reviewer")
        agent = require_nonempty(reviewer["agent_id"], f"{label} reviewer agent ID")
        require_nonempty(reviewer["identity"], f"{label} reviewer identity")
        if agent in agents:
            raise ValueError("F118 independent reviews do not have distinct reviewers")
        agents.add(agent)
        role_agents[kind] = agent
        if (
            not isinstance(review["findings"], list)
            or not review["findings"]
            or any(not isinstance(item, str) or not item for item in review["findings"])
        ):
            raise ValueError(f"{label} findings must be nonempty strings")
        if (
            not isinstance(review["limitations"], list)
            or not review["limitations"]
            or any(not isinstance(item, str) or not item for item in review["limitations"])
        ):
            raise ValueError(f"{label} limitations must be nonempty strings")
        if INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION not in review["limitations"]:
            raise ValueError(
                f"{label} must state the non-cryptographic reviewer identity limitation"
            )
        reviewed = require_utc(review["reviewed_utc"], f"{label} timestamp")
        if reviewed < require_utc(generated_utc, "F118 evidence generation timestamp"):
            raise ValueError(f"{label} predates F118 evidence")
        reviewed_times.append(reviewed)
    declared_process_independence_assurance(
        role_agents, "F118 independent-review process"
    )
    return reviewed_times[0].isoformat(), reviewed_times[1].isoformat()


def parse_audit(value: dict[str, object], *, root: Path, evidence_sha256: str,
                provenance_sha256: str, plasma_sha256: str,
                final: dict[str, object], bridge: dict[str, object],
                predecessor: dict[str, object], f116: dict[str, object],
                catalog_after: dict[str, object],
                reviewed_times: tuple[str, str]) -> None:
    audit = require_exact_keys(
        value,
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "published_utc", "artifact", "independent_reviews",
            "historical_f116_authority", "source_archive_catalog",
            "authority_and_enforcement", "publication",
        },
        "F118 publication audit",
    )
    if (
        audit["schema_version"] != 1
        or audit["record_type"]
        != "stage-i-current-source-authority-supersession-publication-audit"
        or audit["checkpoint"] != CHECKPOINT
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or audit["publication"]
        != "recoverable-forward-transaction-with-publication-audit-commit-marker-under-stage-i-lock"
        or audit["historical_f116_authority"] != f116["digests"]
        or audit["authority_and_enforcement"] != AUTHORIZATION
    ):
        raise ValueError("F118 publication audit identity or authority differs")
    published = require_utc(audit["published_utc"], "F118 publication timestamp")
    if any(published < require_utc(value, "F118 review timestamp") for value in reviewed_times):
        raise ValueError("F118 publication audit predates independent review")
    evidence_path = root / F118_PATHS["evidence"]
    provenance_path = root / F118_PATHS["provenance_review"]
    plasma_path = root / F118_PATHS["plasma_review"]
    declared_binding(audit["artifact"], evidence_path, evidence_sha256, "F118 audit artifact")
    reviews = require_exact_keys(
        audit["independent_reviews"],
        {
            "reviews_bind_exact_published_f118_sha256",
            "provenance_security", "plasma_scientific_continuation",
        },
        "F118 audit independent reviews",
    )
    if reviews["reviews_bind_exact_published_f118_sha256"] != evidence_sha256:
        raise ValueError("F118 audit review digest binding differs")
    declared_binding(
        reviews["provenance_security"], provenance_path, provenance_sha256,
        "F118 provenance review",
    )
    declared_binding(
        reviews["plasma_scientific_continuation"], plasma_path, plasma_sha256,
        "F118 plasma review",
    )
    catalog = require_exact_keys(
        audit["source_archive_catalog"],
        {
            "readme", "sha256sums", "bridge_bundle",
            "predecessor_current_source_bundle", "current_source_bundle",
            "corrupt_c7_absent_from_active_checksum_ledger",
            "sole_current_source_bundle",
        },
        "F118 audit source-archive catalog",
    )
    declared_binding(
        catalog["readme"], root / "source-archives/README.md",
        str(catalog_after["readme_sha256"]), "F118 catalog README", mode="0644",
    )
    declared_binding(
        catalog["sha256sums"], root / "source-archives/SHA256SUMS",
        str(catalog_after["sha256sums_sha256"]), "F118 catalog SHA256SUMS", mode="0644",
    )
    if catalog["bridge_bundle"] != {
        "path": str(root / str(bridge["path"])),
        "sha256": bridge["sha256"],
        "mode": "0644",
        "links": 1,
        "head": bridge["head"],
        "role": "retained-non-current-bridge",
        "selected_as_current": False,
    }:
        raise ValueError("F118 audit bridge binding differs")
    if catalog["predecessor_current_source_bundle"] != {
        "path": str(root / str(predecessor["path"])),
        "sha256": predecessor["sha256"],
        "mode": "0644",
        "links": 1,
        "head": predecessor["head"],
        "role": "retained-non-current-predecessor",
        "selected_as_current": False,
    }:
        raise ValueError("F118 audit predecessor current bundle binding differs")
    if catalog["current_source_bundle"] != {
        "path": str(root / str(final["path"])),
        "sha256": final["sha256"],
        "mode": "0644",
        "links": 1,
        "head": final["head"],
        "selected_as_current": True,
    }:
        raise ValueError("F118 audit current bundle binding differs")
    if (
        catalog["corrupt_c7_absent_from_active_checksum_ledger"] is not True
        or catalog["sole_current_source_bundle"] != str(root / str(final["path"]))
    ):
        raise ValueError("F118 audit catalog authority differs")


def validate_candidate_payloads(*, root: Path, repository: Path, canonical: bool,
                                publisher_sha256: str, bundle_path: Path,
                                bundle_sha256: str, evidence_path: Path,
                                expected_bundle_candidate_path: Path | None = None,
                                evidence_payload: bytes, evidence_sha256: str,
                                provenance_payload: bytes, provenance_sha256: str,
                                plasma_payload: bytes, plasma_sha256: str,
                                audit_payload: bytes, audit_sha256: str,
                                old_readme: bytes, old_sums: bytes,
                                new_readme: bytes, new_sums: bytes
                                ) -> dict[str, object]:
    try:
        evidence = json.loads(evidence_payload)
        provenance = json.loads(provenance_payload)
        plasma = json.loads(plasma_payload)
        audit = json.loads(audit_payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("F118 transaction payload is not valid JSON") from error
    for value, payload, label in (
        (evidence, evidence_payload, "F118 evidence"),
        (provenance, provenance_payload, "F118 provenance review"),
        (plasma, plasma_payload, "F118 plasma review"),
        (audit, audit_payload, "F118 publication audit"),
    ):
        if not isinstance(value, dict) or canonical_json(value) != payload:
            raise ValueError(f"{label} payload is not stable canonical JSON")
    preliminary = require_exact_keys(
        evidence.get("implementation"),
        {
            "publisher", "committed_tools", "intermediate_36140_bundle",
            "predecessor_current_source_bundle", "current_source_bundle",
        },
        "F118 implementation",
    )
    bridge = parse_bundle_declaration(
        preliminary["intermediate_36140_bundle"], "F118 bridge bundle", current=False
    )
    final = parse_bundle_declaration(
        preliminary["current_source_bundle"], "F118 current source bundle", current=True
    )
    predecessor = parse_predecessor_bundle_declaration(
        preliminary["predecessor_current_source_bundle"],
        "F118 predecessor current source bundle",
    )
    if final["sha256"] != bundle_sha256:
        raise ValueError("F118 evidence final bundle digest differs from selected candidate")
    if (
        expected_bundle_candidate_path is not None
        and final["candidate_path"] != expected_bundle_candidate_path
    ):
        raise ValueError("F118 transaction final-bundle candidate path differs")
    bridge_path = root / str(bridge["path"])
    stable_bundle_validation(
        repository,
        bridge_path,
        str(bridge["sha256"]),
        str(bridge["head"]),
        str(bridge["advertised_tip"]["name"]),
        list(bridge["verified_revisions"]),
        "F118 retained bridge bundle",
    )
    stable_bundle_validation(
        repository,
        root / str(predecessor["path"]),
        str(predecessor["sha256"]),
        str(predecessor["head"]),
        "HEAD",
        list(predecessor["verified_revisions"]),
        "F118 retained predecessor current source bundle",
    )
    stable_bundle_validation(
        repository,
        bundle_path,
        bundle_sha256,
        str(final["head"]),
        "HEAD",
        list(final["verified_revisions"]),
        "F118 final source bundle",
    )
    final_payload, _ = read_file(
        bundle_path, "F118 final source bundle payload",
        expected=bundle_sha256, mode=0o644,
    )
    catalog = validate_catalog_payload_pair(
        root, old_readme, old_sums, new_readme, new_sums, bridge, predecessor, final,
        final_payload,
    )
    parsed = parse_evidence(
        evidence,
        root=root,
        repository=repository,
        evidence_candidate_path=evidence_path,
        bundle_candidate_path=Path(str(final["candidate_path"])),
        publisher_sha256=publisher_sha256,
        canonical=canonical,
        old_catalog=catalog["before"],
        new_catalog=catalog["after"],
    )
    f116 = historical_f116(
        root, repository, parsed["historical_f116_bindings"], canonical=canonical
    )
    expected_predecessor = parse_predecessor_bundle_declaration(
        retained_predecessor_declaration(f116),
        "authenticated historical F116 predecessor current source bundle",
    )
    if bridge != f116["bridge"] or predecessor != expected_predecessor:
        raise ValueError("F118 retained predecessor bundle declarations differ from F116")
    if f116["bundle"]["head"] not in set(final["verified_revisions"]):
        raise ValueError("F118 final bundle does not cover historical F116 source")
    f115_name = Path(str(f116["historical_f115"]["bundle"]["path"])).name
    for payload, label in (
        (old_sums, "F118 predecessor checksum ledger"),
        (new_sums, "F118 published checksum ledger"),
    ):
        if [name for _, name in parse_sha256sums(payload, label)].count(f115_name) != 1:
            raise ValueError(f"{label} does not preserve historical F115 exactly once")
    reviewed_times = parse_reviews(
        provenance,
        plasma,
        evidence_candidate_path=evidence_path,
        evidence_sha256=evidence_sha256,
        evidence_path=root / F118_PATHS["evidence"],
        final=final,
        generated_utc=str(parsed["generated_utc"]),
    )
    parse_audit(
        audit,
        root=root,
        evidence_sha256=evidence_sha256,
        provenance_sha256=provenance_sha256,
        plasma_sha256=plasma_sha256,
        final=final,
        bridge=bridge,
        predecessor=predecessor,
        f116=f116,
        catalog_after=catalog["after"],
        reviewed_times=reviewed_times,
    )
    return {
        "evidence": evidence,
        "bridge": bridge,
        "predecessor": predecessor,
        "final": final,
        "f116": f116,
        "catalog": catalog,
        "audit_sha256": audit_sha256,
    }


def validate_catalog_payload_pair(root: Path, old_readme: bytes, old_sums: bytes,
                                  new_readme: bytes, new_sums: bytes,
                                  bridge: dict[str, object],
                                  predecessor: dict[str, object],
                                  final: dict[str, object],
                                  final_payload: bytes) -> dict[str, object]:
    generated_readme, generated_sums = catalog_payloads(
        old_readme,
        old_sums,
        bridge_name=str(bridge["name"]),
        predecessor_name=str(predecessor["name"]),
        final_name=str(final["name"]),
        final_sha256=str(final["sha256"]),
        final_revision=str(final["head"]),
        final_subject=str(final["subject"]),
    )
    if generated_readme != new_readme or generated_sums != new_sums:
        raise ValueError("F118 transaction catalog payloads are not deterministic")
    entries = parse_sha256sums(old_sums, "old source-archive SHA256SUMS")
    validate_checksum_entries(root / "source-archives", entries)
    old_names = [name for _, name in entries]
    if (
        old_names.count(str(bridge["name"])) != 1
        or old_names.count(str(predecessor["name"])) != 1
        or str(final["name"]) in old_names
        or CORRUPT_C7_NAME in old_names
    ):
        raise ValueError("old source-archive catalog violates F118 predecessor policy")
    new_entries = parse_sha256sums(new_sums, "new source-archive SHA256SUMS")
    validate_checksum_entries(
        root / "source-archives",
        new_entries,
        final_name=str(final["name"]),
        final_payload=final_payload,
    )
    new_names = [name for _, name in new_entries]
    if (
        new_names.count(str(bridge["name"])) != 1
        or new_names.count(str(predecessor["name"])) != 1
        or new_names.count(str(final["name"])) != 1
        or CORRUPT_C7_NAME in new_names
        or new_entries[:-1] != entries
    ):
        raise ValueError("new source-archive catalog violates F118 policy")
    return {
        "before": {
            "readme_sha256": sha256_bytes(old_readme),
            "sha256sums_sha256": sha256_bytes(old_sums),
            "bridge_listed_exactly_once": True,
            "predecessor_current_source_bundle_listed_exactly_once": True,
            "final_bundle_listed": False,
            "corrupt_c7_listed": False,
            "historical_f115_preserved": True,
        },
        "after": {
            "readme_sha256": sha256_bytes(new_readme),
            "sha256sums_sha256": sha256_bytes(new_sums),
            "bridge_listed_exactly_once": True,
            "predecessor_current_source_bundle_listed_exactly_once": True,
            "final_bundle_listed_exactly_once": True,
            "corrupt_c7_listed": False,
            "historical_f115_preserved": True,
            "historical_f116_preserved": True,
            "all_prior_checksum_entries_preserved": True,
            "sole_current_source_bundle": str(final["path"]),
        },
    }


def root_layout(root: Path) -> dict[str, Path]:
    return {
        "root": root,
        "accounting": root / "accounting",
        "source_archives": root / "source-archives",
        "lock": root / f".mks24_stage_i_{EXECUTION_EPOCH_SLUG}.lock",
        "transactions": root / "accounting" / TRANSACTION_ROOT_NAME,
        "readme": root / "source-archives/README.md",
        "sha256sums": root / "source-archives/SHA256SUMS",
        **{f"f118_{key}": root / relative for key, relative in F118_PATHS.items()},
    }


def validate_root(args: argparse.Namespace) -> tuple[Path, bool, dict[str, Path]]:
    root = normalized_absolute(args.root, "campaign root")
    canonical = root == DEFAULT_ROOT
    if canonical and args.allow_local_root:
        raise ValueError("--allow-local-root is forbidden for canonical execution")
    if not canonical and not args.allow_local_root:
        raise ValueError("--allow-local-root is required outside the canonical root")
    layout = root_layout(root)
    for key in ("root", "accounting", "source_archives"):
        require_directory(layout[key], key.replace("_", " "))
    # Promotion and recovery classify catalog crash states only after the Stage I
    # lock is bound.  Read-only actions continue to require complete catalogs.
    if args.action not in {"promote", "recover"}:
        _, _ = read_file(layout["readme"], "source-archive README", mode=0o644)
        _, _ = read_file(layout["sha256sums"], "source-archive SHA256SUMS", mode=0o644)
    return root, canonical, layout


def require_canonical_repository(repository: Path, canonical: bool) -> None:
    """Pin canonical publication and verification to the retained source repository."""

    if canonical and repository != CANONICAL_REPOSITORY_ROOT:
        raise ValueError(
            "canonical use requires the retained repository root "
            f"{CANONICAL_REPOSITORY_ROOT}: {repository}"
        )


def require_canonical_public_namespace_profile(
    path: Path, profile: os.stat_result
) -> None:
    """Require one exact trusted canonical public-namespace directory profile."""

    expected = next(
        (
            (mode, uid, gid)
            for candidate, mode, uid, gid in CANONICAL_PUBLIC_NAMESPACE_PROFILES
            if candidate == path
        ),
        None,
    )
    if expected is None:
        raise ValueError(f"canonical public namespace path is not trusted: {path}")
    mode, uid, gid = expected
    retained_mode = stat.S_IMODE(profile.st_mode)
    if (
        not stat.S_ISDIR(profile.st_mode)
        or retained_mode != mode
        or profile.st_uid != uid
        or profile.st_gid != gid
        or retained_mode & stat.S_IWOTH
    ):
        raise ValueError(f"canonical public namespace profile is untrusted: {path}")
    if retained_mode & stat.S_IWGRP and path != CANONICAL_TRUSTED_PROJECT_BOUNDARY:
        raise ValueError(f"canonical public namespace profile is group-writable: {path}")


@contextmanager
def bound_canonical_public_namespace(root: Path, canonical: bool) -> Iterator[None]:
    """Bind the exact canonical project namespace across one complete action.

    The root-owned setgid project boundary is intentionally group-writable for
    trusted project GID 31114. Every descendant public boundary is exact,
    user-owned, and non-writable by group or world.
    """

    global _ACTIVE_CANONICAL_PUBLIC_NAMESPACE

    if not canonical:
        yield
        return
    if root != DEFAULT_ROOT:
        raise ValueError("canonical public namespace binding requires the canonical root")
    if _ACTIVE_CANONICAL_PUBLIC_NAMESPACE is not None:
        raise ValueError("nested canonical public namespace bindings are forbidden")
    bindings: list[tuple[Path, int, os.stat_result]] = []
    with ExitStack() as stack:
        for path, _, _, _ in CANONICAL_PUBLIC_NAMESPACE_PROFILES:
            descriptor = stack.enter_context(
                absolute_descriptor(
                    path,
                    f"canonical public namespace {path}",
                    flags=os.O_RDONLY | os.O_DIRECTORY,
                )
            )
            profile = os.fstat(descriptor)
            require_canonical_public_namespace_profile(path, profile)
            named = os.stat(path, follow_symlinks=False)
            require_canonical_public_namespace_profile(path, named)
            if profile_identity(named) != profile_identity(profile):
                raise ValueError(f"canonical public namespace path changed: {path}")
            bindings.append((path, descriptor, profile))
        guard = CanonicalPublicNamespaceGuard(bindings)
        guard.assert_bound()
        _ACTIVE_CANONICAL_PUBLIC_NAMESPACE = guard
        try:
            yield
        finally:
            try:
                guard.assert_bound()
            finally:
                _ACTIVE_CANONICAL_PUBLIC_NAMESPACE = None


@contextmanager
def stage_i_lock(layout: dict[str, Path]) -> Iterator[None]:
    global _ACTIVE_MUTATION_LOCK

    path = layout["lock"]
    require_no_symlinks(path, "Stage I lock")
    if _ACTIVE_MUTATION_LOCK is not None:
        raise ValueError("nested source-authority mutation locks are forbidden")
    with bound_directory(path.parent, "Stage I lock parent") as (parent, _):
        descriptor = os.open(path.name, os.O_RDWR | os.O_NOFOLLOW, dir_fd=parent)
        acquired = False
        try:
            profile = os.fstat(descriptor)
            require_file_profile(profile, "Stage I lock", mode=0o644)
            if profile.st_size != 0:
                raise ValueError("Stage I lock must be empty")
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                acquired = True
            except OSError as error:
                if error.errno not in (errno.EACCES, errno.EAGAIN):
                    raise
                raise ValueError(f"another Stage I mutation holds {path}") from error
            guard = MutationLockGuard(parent, path.name, descriptor, profile, "Stage I lock")
            guard.assert_bound()
            _ACTIVE_MUTATION_LOCK = guard
            try:
                yield
            finally:
                try:
                    guard.assert_bound()
                finally:
                    _ACTIVE_MUTATION_LOCK = None
        finally:
            if acquired:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
            os.close(descriptor)


def require_other_transactions_empty(layout: dict[str, Path]) -> None:
    for name in OTHER_TRANSACTION_NAMES:
        path = layout["accounting"] / name
        if os.path.lexists(path):
            with bound_directory(path, f"{name} directory") as (descriptor, _):
                if os.listdir(descriptor):
                    raise ValueError(f"{name} directory is not empty")


def transaction_root_entry_bindings(
    parent: int, label: str
) -> tuple[tuple[str, tuple[object, ...]], ...]:
    """Bind exact transaction-root names, inode identities, and profiles."""

    require_active_mutation_lock_bound()
    require_bound_directory_descriptor(parent)
    bindings: list[tuple[str, tuple[object, ...]]] = []
    for name in sorted(os.listdir(parent)):
        profile = os.stat(name, dir_fd=parent, follow_symlinks=False)
        require_bound_entry_profile(parent, name, profile, f"{label} entry {name}")
        if stat.S_ISREG(profile.st_mode):
            binding: tuple[object, ...] = (
                "file",
                *profile_identity(profile),
                *file_profile_binding(profile),
            )
        elif stat.S_ISDIR(profile.st_mode):
            binding = (
                "directory",
                *profile_identity(profile),
                *directory_security_identity(profile),
            )
        else:
            raise ValueError(f"{label} entry {name} has an unsupported type")
        bindings.append((name, binding))
    require_bound_directory_descriptor(parent)
    require_active_mutation_lock_bound()
    return tuple(bindings)


def stable_transaction_directory_scan(
    descriptor: int,
) -> tuple[
    tuple[tuple[str, tuple[object, ...]], ...],
    tuple[tuple[str, tuple[int, int], tuple[int, int, int]], ...],
] | None:
    """Return one internally stable identity-bearing transaction-root scan."""

    label = "source-authority transaction root"
    before = transaction_root_entry_bindings(descriptor, label)
    retained = []
    for name in non_forensic_entry_names(descriptor, label):
        with bound_child_directory(
            descriptor, name, "source-authority transaction directory", mode=0o700
        ) as (_, profile):
            retained.append(
                (name, profile_identity(profile), directory_security_identity(profile))
            )
    after = transaction_root_entry_bindings(descriptor, label)
    if before != after:
        return None
    require_bound_directory_descriptor(descriptor)
    require_active_mutation_lock_bound()
    return after, tuple(retained)


def transaction_directories(layout: dict[str, Path]) -> list[Path]:
    path = layout["transactions"]
    if not os.path.lexists(path):
        return []
    with bound_directory(path, "source-authority transaction root") as (descriptor, _):
        previous = None
        consecutive = 0
        for _ in range(8):
            current = stable_transaction_directory_scan(descriptor)
            if current is None:
                previous = None
                consecutive = 0
                continue
            if current == previous:
                consecutive += 1
            else:
                previous = current
                consecutive = 1
            if consecutive >= 2:
                require_bound_directory_descriptor(descriptor)
                require_active_mutation_lock_bound()
                return [path / name for name, _, _ in current[1]]
        raise ValueError("source-authority transaction root did not stabilize")


def non_forensic_entry_names(parent: int, label: str) -> list[str]:
    """Return active entries while validating retained forensic siblings."""

    active = []
    regular_forensics: list[tuple[str, os.stat_result]] = []
    for name in sorted(os.listdir(parent)):
        if FORENSIC_ENTRY_RE.fullmatch(name) is None:
            active.append(name)
            continue
        profile = os.stat(name, dir_fd=parent, follow_symlinks=False)
        forensic_label = f"{label} forensic entry {name}"
        if stat.S_ISREG(profile.st_mode):
            if profile.st_nlink not in {1, 2}:
                raise ValueError(f"{forensic_label} has an invalid link profile")
            require_file_profile(profile, forensic_label, links=profile.st_nlink)
            regular_forensics.append((name, profile))
        elif stat.S_ISDIR(profile.st_mode):
            require_directory_profile(profile, forensic_label)
        else:
            raise ValueError(f"{label} forensic entry {name} has an invalid type")
        require_bound_entry_profile(parent, name, profile, forensic_label)
    for name, profile in regular_forensics:
        if profile.st_nlink == 2 and sum(
            profile_identity(other) == profile_identity(profile)
            for _, other in regular_forensics
        ) != 2:
            raise ValueError(
                f"{label} forensic entry {name} has an external retained hardlink"
            )
    return active


def require_only_forensic_entries(parent: int, label: str) -> None:
    """Require one retired container to hold only forensic siblings."""

    active = non_forensic_entry_names(parent, label)
    if active:
        raise ValueError(f"{label} contains unexpected active entries: {active}")


def candidate_file(path: Path, expected: str, label: str) -> tuple[bytes, str]:
    return read_file(normalized_absolute(path, label), label, expected=expected)


def require_external_draft_output(path: Path, root: Path, label: str) -> Path:
    retained = normalized_absolute(path, label)
    require_no_symlinks(retained, label)
    try:
        retained.relative_to(root)
    except ValueError:
        pass
    else:
        raise ValueError(f"{label} must remain outside the campaign root")
    if retained.exists() or retained.is_symlink():
        raise ValueError(f"{label} already exists")
    require_directory(retained.parent, f"{label} parent")
    return retained


def require_draftable_state(root: Path, repository: Path,
                            layout: dict[str, Path]) -> None:
    require_other_transactions_empty(layout)
    for key in F118_PATHS:
        if os.path.lexists(layout[f"f118_{key}"]):
            raise ValueError(f"F118 target already exists: {layout[f'f118_{key}']}")
    require_inert_stale_transactions(root, repository, layout)


def current_f116_bindings(root: Path, canonical: bool) -> dict[str, object]:
    retained = {}
    for key, relative in F116_PATHS.items():
        _, digest = read_file(root / relative, f"historical F116 {key}", mode=0o444)
        if canonical and digest != F116_CANONICAL_SHA256[key]:
            raise ValueError(f"historical F116 {key} digest differs from canonical authority")
        retained[key] = {"path": relative.as_posix(), "sha256": digest}
    return retained


def revision_subject(repository: Path, revision: str, label: str) -> str:
    completed = git_run(repository, ["show", "-s", "--format=%s", revision])
    if completed.returncode:
        raise ValueError(f"cannot inspect {label} subject")
    try:
        subject = completed.stdout.decode("utf-8").rstrip("\n")
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} subject is not UTF-8") from error
    return require_nonempty(subject, f"{label} subject")


def publication_binding(path: Path, digest: str, mode: str) -> dict[str, object]:
    return {"path": str(path), "sha256": digest, "mode": mode, "links": 1}


def draft_evidence(args: argparse.Namespace, root: Path, repository: Path,
                   canonical: bool, publisher_sha256: str,
                   layout: dict[str, Path]) -> Path:
    require_draftable_state(root, repository, layout)
    output = require_external_draft_output(
        args.evidence_candidate_output, root, "F118 evidence candidate output"
    )
    bundle_path = normalized_absolute(args.bundle_candidate, "final bundle candidate")
    bridge_path = normalized_absolute(args.bridge_bundle, "retained bridge bundle")
    if bridge_path.parent != layout["source_archives"]:
        raise ValueError("retained bridge bundle must be a direct source-archive child")
    head = repository_head(repository)
    bridge_head = require_revision(args.bridge_head, "retained bridge bundle head")
    bridge_sha256 = require_sha256(
        args.expected_bridge_sha256, "retained bridge bundle SHA-256"
    )
    bundle_sha256 = require_sha256(args.expected_bundle_sha256, "final bundle SHA-256")
    if canonical and (
        bridge_path.name != BRIDGE_NAME
        or bridge_head != BRIDGE_REVISION
        or bridge_sha256 != BRIDGE_SHA256
    ):
        raise ValueError("canonical retained bridge identity differs")
    final_target = layout["source_archives"] / bundle_path.name
    if final_target.exists():
        raise ValueError(f"F118 target already exists: {final_target}")

    f116_bindings = current_f116_bindings(root, canonical)
    f116 = historical_f116(root, repository, f116_bindings, canonical=canonical)
    bridge_revisions = set(f116["bridge"]["verified_revisions"])
    final_revisions = set(f116["bundle"]["verified_revisions"]) | {head}
    stable_bundle_validation(
        repository,
        bridge_path,
        bridge_sha256,
        bridge_head,
        "refs/heads/feature/cgl-landau-fluid",
        sorted(bridge_revisions),
        "F118 retained bridge bundle",
    )
    stable_bundle_validation(
        repository,
        bundle_path,
        bundle_sha256,
        head,
        "HEAD",
        sorted(final_revisions),
        "F118 final source bundle candidate",
    )
    final_payload, _ = read_file(
        bundle_path, "F118 final source bundle candidate payload",
        expected=bundle_sha256, mode=0o644,
    )
    bridge = {
        "path": bridge_path.relative_to(root).as_posix(),
        "sha256": bridge_sha256,
        "complete_history": True,
        "head": bridge_head,
        "advertised_tip": {
            "revision": bridge_head,
            "name": "refs/heads/feature/cgl-landau-fluid",
        },
        "verified_revisions": sorted(bridge_revisions),
        "selected_as_current": False,
        "role": "retained-non-current-bridge",
    }
    predecessor = retained_predecessor_declaration(f116)
    final = {
        "candidate_path": str(bundle_path),
        "path": (PurePosixPath("source-archives") / bundle_path.name).as_posix(),
        "sha256": bundle_sha256,
        "complete_history": True,
        "head": head,
        "advertised_tip": {"revision": head, "name": "HEAD"},
        "verified_revisions": sorted(final_revisions),
        "selected_as_current": True,
        "subject": revision_subject(repository, head, "F118 final HEAD"),
    }
    parsed_bridge = parse_bundle_declaration(
        bridge, "F118 retained bridge bundle", current=False
    )
    parsed_predecessor = parse_predecessor_bundle_declaration(
        predecessor, "F118 retained predecessor current source bundle"
    )
    parsed_final = parse_bundle_declaration(
        final, "F118 current source bundle", current=True
    )
    if parsed_bridge != f116["bridge"] or parsed_predecessor != parse_predecessor_bundle_declaration(
        retained_predecessor_declaration(f116),
        "authenticated historical F116 predecessor current source bundle",
    ):
        raise ValueError("F118 retained predecessor bundle arguments differ from F116")
    catalog = validate_catalog_before(
        root, f116, parsed_bridge, parsed_predecessor, parsed_final, final_payload
    )
    generated_utc = require_utc(
        args.generated_utc, "F118 evidence generation timestamp"
    ).isoformat()
    evidence = {
        "schema_version": 1,
        "record_type": "stage-i-current-source-authority-supersession-evidence",
        "checkpoint": CHECKPOINT,
        "execution_epoch": EXECUTION_EPOCH,
        "generated_utc": generated_utc,
        "scope": {
            "relationship": "current-source-selection-only-supersession",
            "summary": "Select the exact committed final tooling source without execution authority.",
            "preserves": PRESERVES,
            "does_not_authorize": DOES_NOT_AUTHORIZE,
        },
        "predecessor_authorities": {"historical_f116": f116_bindings},
        "implementation": {
            "publisher": {
                "path": PUBLISHER_RELATIVE.as_posix(),
                "revision": head,
                "sha256": publisher_sha256,
                "mode": REQUIRED_TOOLS[PUBLISHER_RELATIVE.as_posix()],
            },
            "committed_tools": committed_tools(repository, head, publisher_sha256),
            "intermediate_36140_bundle": bridge,
            "predecessor_current_source_bundle": predecessor,
            "current_source_bundle": final,
        },
        "source_archive_catalog": {
            "before": catalog["before"],
            "after": catalog["after"],
        },
        "authorization": AUTHORIZATION,
        "validation": VALIDATION_CLAIMS,
        "publication_requirements": PUBLICATION_REQUIREMENTS,
    }
    parse_evidence(
        evidence,
        root=root,
        repository=repository,
        evidence_candidate_path=output,
        bundle_candidate_path=bundle_path,
        publisher_sha256=publisher_sha256,
        canonical=canonical,
        old_catalog=catalog["before"],
        new_catalog=catalog["after"],
    )
    payload = canonical_json(evidence)
    write_exclusive(output, payload, 0o444)
    read_file(output, "drafted F118 evidence candidate", expected=sha256_bytes(payload), mode=0o444)
    return output


def reviewed_draft_context(args: argparse.Namespace, root: Path, repository: Path,
                           canonical: bool, publisher_sha256: str,
                           layout: dict[str, Path]) -> dict[str, object]:
    require_draftable_state(root, repository, layout)
    bundle_path = normalized_absolute(args.bundle_candidate, "final bundle candidate")
    evidence_path = normalized_absolute(args.evidence_candidate, "F118 evidence candidate")
    evidence_payload, evidence_sha256 = candidate_file(
        evidence_path, args.expected_evidence_sha256, "F118 evidence candidate"
    )
    evidence = parse_json_payload(evidence_payload, "F118 evidence candidate")
    implementation = evidence.get("implementation")
    if not isinstance(implementation, dict):
        raise ValueError("F118 evidence implementation is missing")
    bridge = parse_bundle_declaration(
        implementation.get("intermediate_36140_bundle"), "F118 bridge bundle", current=False
    )
    final = parse_bundle_declaration(
        implementation.get("current_source_bundle"), "F118 current source bundle", current=True
    )
    predecessor = parse_predecessor_bundle_declaration(
        implementation.get("predecessor_current_source_bundle"),
        "F118 predecessor current source bundle",
    )
    bundle_sha256 = require_sha256(args.expected_bundle_sha256, "final bundle SHA-256")
    if final["candidate_path"] != bundle_path or final["sha256"] != bundle_sha256:
        raise ValueError("F118 evidence does not bind the selected final bundle candidate")
    final_target = root / str(final["path"])
    if final_target.exists():
        raise ValueError(f"F118 target already exists: {final_target}")
    predecessors = require_exact_keys(
        evidence.get("predecessor_authorities"), {"historical_f116"},
        "F118 predecessor authorities",
    )
    f116 = historical_f116(
        root, repository, predecessors["historical_f116"], canonical=canonical
    )
    if bridge != f116["bridge"] or predecessor != parse_predecessor_bundle_declaration(
        retained_predecessor_declaration(f116),
        "authenticated historical F116 predecessor current source bundle",
    ):
        raise ValueError("F118 retained predecessor bundle declarations differ from F116")
    stable_bundle_validation(
        repository,
        root / str(bridge["path"]),
        str(bridge["sha256"]),
        str(bridge["head"]),
        str(bridge["advertised_tip"]["name"]),
        list(bridge["verified_revisions"]),
        "F118 retained bridge bundle",
    )
    stable_bundle_validation(
        repository,
        root / str(predecessor["path"]),
        str(predecessor["sha256"]),
        str(predecessor["head"]),
        "HEAD",
        list(predecessor["verified_revisions"]),
        "F118 retained predecessor current source bundle",
    )
    stable_bundle_validation(
        repository,
        bundle_path,
        bundle_sha256,
        str(final["head"]),
        "HEAD",
        list(final["verified_revisions"]),
        "F118 final source bundle candidate",
    )
    final_payload, _ = read_file(
        bundle_path, "F118 final source bundle candidate payload",
        expected=bundle_sha256, mode=0o644,
    )
    catalog = validate_catalog_before(
        root, f116, bridge, predecessor, final, final_payload
    )
    parsed = parse_evidence(
        evidence,
        root=root,
        repository=repository,
        evidence_candidate_path=evidence_path,
        bundle_candidate_path=bundle_path,
        publisher_sha256=publisher_sha256,
        canonical=canonical,
        old_catalog=catalog["before"],
        new_catalog=catalog["after"],
    )
    reviews = {}
    review_sha256 = {}
    for key, path, expected in (
        ("provenance", args.provenance_review_candidate,
         args.expected_provenance_review_sha256),
        ("plasma", args.plasma_review_candidate, args.expected_plasma_review_sha256),
    ):
        payload, digest = candidate_file(
            path, expected, f"F118 {key} review candidate"
        )
        reviews[key] = parse_json_payload(payload, f"F118 {key} review candidate")
        review_sha256[key] = digest
    reviewed_times = parse_reviews(
        reviews["provenance"],
        reviews["plasma"],
        evidence_candidate_path=evidence_path,
        evidence_sha256=evidence_sha256,
        evidence_path=layout["f118_evidence"],
        final=final,
        generated_utc=str(parsed["generated_utc"]),
    )
    return {
        "evidence_sha256": evidence_sha256,
        "review_sha256": review_sha256,
        "reviewed_times": reviewed_times,
        "f116": f116,
        "bridge": bridge,
        "predecessor": predecessor,
        "final": final,
        "catalog": catalog,
    }


def draft_audit(args: argparse.Namespace, root: Path, repository: Path,
                canonical: bool, publisher_sha256: str,
                layout: dict[str, Path]) -> Path:
    output = require_external_draft_output(
        args.audit_candidate_output, root, "F118 audit candidate output"
    )
    context = reviewed_draft_context(
        args, root, repository, canonical, publisher_sha256, layout
    )
    published_utc = require_utc(
        args.published_utc, "F118 publication timestamp"
    ).isoformat()
    final = context["final"]
    bridge = context["bridge"]
    predecessor = context["predecessor"]
    f116 = context["f116"]
    catalog = context["catalog"]
    review_sha256 = context["review_sha256"]
    if not all(
        isinstance(value, dict)
        for value in (final, bridge, predecessor, f116, catalog, review_sha256)
    ):
        raise ValueError("F118 audit drafting context differs")
    evidence_sha256 = str(context["evidence_sha256"])
    audit = {
        "schema_version": 1,
        "record_type": "stage-i-current-source-authority-supersession-publication-audit",
        "checkpoint": CHECKPOINT,
        "execution_epoch": EXECUTION_EPOCH,
        "published_utc": published_utc,
        "artifact": publication_binding(
            layout["f118_evidence"], evidence_sha256, "0444"
        ),
        "independent_reviews": {
            "reviews_bind_exact_published_f118_sha256": evidence_sha256,
            "provenance_security": publication_binding(
                layout["f118_provenance_review"],
                str(review_sha256["provenance"]),
                "0444",
            ),
            "plasma_scientific_continuation": publication_binding(
                layout["f118_plasma_review"],
                str(review_sha256["plasma"]),
                "0444",
            ),
        },
        "historical_f116_authority": f116["digests"],
        "source_archive_catalog": {
            "readme": publication_binding(
                layout["readme"], str(catalog["after"]["readme_sha256"]), "0644"
            ),
            "sha256sums": publication_binding(
                layout["sha256sums"], str(catalog["after"]["sha256sums_sha256"]), "0644"
            ),
            "bridge_bundle": {
                "path": str(root / str(bridge["path"])),
                "sha256": bridge["sha256"],
                "mode": "0644",
                "links": 1,
                "head": bridge["head"],
                "role": "retained-non-current-bridge",
                "selected_as_current": False,
            },
            "predecessor_current_source_bundle": {
                "path": str(root / str(predecessor["path"])),
                "sha256": predecessor["sha256"],
                "mode": "0644",
                "links": 1,
                "head": predecessor["head"],
                "role": "retained-non-current-predecessor",
                "selected_as_current": False,
            },
            "current_source_bundle": {
                "path": str(root / str(final["path"])),
                "sha256": final["sha256"],
                "mode": "0644",
                "links": 1,
                "head": final["head"],
                "selected_as_current": True,
            },
            "corrupt_c7_absent_from_active_checksum_ledger": True,
            "sole_current_source_bundle": str(root / str(final["path"])),
        },
        "authority_and_enforcement": AUTHORIZATION,
        "publication": (
            "recoverable-forward-transaction-with-publication-audit-commit-marker-"
            "under-stage-i-lock"
        ),
    }
    parse_audit(
        audit,
        root=root,
        evidence_sha256=evidence_sha256,
        provenance_sha256=str(review_sha256["provenance"]),
        plasma_sha256=str(review_sha256["plasma"]),
        final=final,
        bridge=bridge,
        predecessor=predecessor,
        f116=f116,
        catalog_after=catalog["after"],
        reviewed_times=context["reviewed_times"],
    )
    payload = canonical_json(audit)
    write_exclusive(output, payload, 0o444)
    read_file(output, "drafted F118 audit candidate", expected=sha256_bytes(payload), mode=0o444)
    return output


def parse_json_payload(payload: bytes, label: str) -> dict[str, object]:
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON") from error
    if not isinstance(value, dict) or canonical_json(value) != payload:
        raise ValueError(f"{label} must be stable canonical JSON")
    return value


def initial_plan(args: argparse.Namespace, root: Path, repository: Path,
                 canonical: bool, publisher_sha256: str,
                 layout: dict[str, Path]) -> dict[str, object]:
    require_other_transactions_empty(layout)
    require_inert_stale_transactions(root, repository, layout)
    candidate_paths = {
        "bundle": normalized_absolute(args.bundle_candidate, "final bundle candidate"),
        "evidence": normalized_absolute(args.evidence_candidate, "F118 evidence candidate"),
        "provenance_review": normalized_absolute(
            args.provenance_review_candidate, "F118 provenance review candidate"
        ),
        "plasma_review": normalized_absolute(
            args.plasma_review_candidate, "F118 plasma review candidate"
        ),
        "audit": normalized_absolute(args.audit_candidate, "F118 audit candidate"),
    }
    expected = {
        "bundle": args.expected_bundle_sha256,
        "evidence": args.expected_evidence_sha256,
        "provenance_review": args.expected_provenance_review_sha256,
        "plasma_review": args.expected_plasma_review_sha256,
        "audit": args.expected_audit_sha256,
    }
    payloads = {}
    for key, path in candidate_paths.items():
        payloads[key], observed = candidate_file(path, expected[key], f"F118 {key} candidate")
        if observed != expected[key]:
            raise ValueError(f"F118 {key} candidate checksum differs")
    evidence = parse_json_payload(payloads["evidence"], "F118 evidence candidate")
    implementation = evidence.get("implementation")
    if not isinstance(implementation, dict):
        raise ValueError("F118 evidence implementation is missing")
    bridge = parse_bundle_declaration(
        implementation.get("intermediate_36140_bundle"), "F118 bridge bundle", current=False
    )
    final = parse_bundle_declaration(
        implementation.get("current_source_bundle"), "F118 current source bundle", current=True
    )
    predecessor = parse_predecessor_bundle_declaration(
        implementation.get("predecessor_current_source_bundle"),
        "F118 predecessor current source bundle",
    )
    if final["candidate_path"] != candidate_paths["bundle"]:
        raise ValueError("F118 evidence does not bind the selected final bundle candidate")
    if final["sha256"] != expected["bundle"]:
        raise ValueError("F118 evidence final bundle SHA-256 differs")
    if canonical and (
        bridge["head"] != BRIDGE_REVISION
        or bridge["sha256"] != BRIDGE_SHA256
        or bridge["name"] != BRIDGE_NAME
    ):
        raise ValueError("canonical F118 bridge identity differs")
    final_target = root / str(final["path"])
    for path in (
        final_target,
        layout["f118_evidence"],
        layout["f118_provenance_review"],
        layout["f118_plasma_review"],
        layout["f118_publication_audit"],
    ):
        if path.exists():
            raise ValueError(f"F118 target already exists: {path}")
    f116_bindings = require_exact_keys(
        evidence.get("predecessor_authorities"), {"historical_f116"},
        "F118 predecessor authorities",
    )["historical_f116"]
    f116 = historical_f116(root, repository, f116_bindings, canonical=canonical)
    if bridge != f116["bridge"] or predecessor != parse_predecessor_bundle_declaration(
        retained_predecessor_declaration(f116),
        "authenticated historical F116 predecessor current source bundle",
    ):
        raise ValueError("F118 retained predecessor bundle declarations differ from F116")
    stable_bundle_validation(
        repository,
        root / str(bridge["path"]),
        str(bridge["sha256"]),
        str(bridge["head"]),
        str(bridge["advertised_tip"]["name"]),
        list(bridge["verified_revisions"]),
        "F118 retained bridge bundle",
    )
    stable_bundle_validation(
        repository,
        root / str(predecessor["path"]),
        str(predecessor["sha256"]),
        str(predecessor["head"]),
        "HEAD",
        list(predecessor["verified_revisions"]),
        "F118 retained predecessor current source bundle",
    )
    stable_bundle_validation(
        repository,
        candidate_paths["bundle"],
        expected["bundle"],
        str(final["head"]),
        "HEAD",
        list(final["verified_revisions"]),
        "F118 final source bundle candidate",
    )
    catalog = validate_catalog_before(
        root, f116, bridge, predecessor, final, payloads["bundle"]
    )
    parsed = parse_evidence(
        evidence,
        root=root,
        repository=repository,
        evidence_candidate_path=candidate_paths["evidence"],
        bundle_candidate_path=candidate_paths["bundle"],
        publisher_sha256=publisher_sha256,
        canonical=canonical,
        old_catalog=catalog["before"],
        new_catalog=catalog["after"],
    )
    if parsed["historical_f116_bindings"] != f116_bindings:
        raise ValueError("F118 historical F116 authority changed during validation")
    provenance = parse_json_payload(payloads["provenance_review"], "F118 provenance review")
    plasma = parse_json_payload(payloads["plasma_review"], "F118 plasma review")
    reviewed_times = parse_reviews(
        provenance,
        plasma,
        evidence_candidate_path=candidate_paths["evidence"],
        evidence_sha256=expected["evidence"],
        evidence_path=layout["f118_evidence"],
        final=final,
        generated_utc=str(parsed["generated_utc"]),
    )
    audit = parse_json_payload(payloads["audit"], "F118 publication audit")
    parse_audit(
        audit,
        root=root,
        evidence_sha256=expected["evidence"],
        provenance_sha256=expected["provenance_review"],
        plasma_sha256=expected["plasma_review"],
        final=final,
        bridge=bridge,
        predecessor=predecessor,
        f116=f116,
        catalog_after=catalog["after"],
        reviewed_times=reviewed_times,
    )
    payloads.update(
        {
            "readme_before": catalog["old_readme"],
            "sha256sums_before": catalog["old_sums"],
            "readme_after": catalog["new_readme"],
            "sha256sums_after": catalog["new_sums"],
        }
    )
    return {
        "candidate_paths": {key: str(value) for key, value in candidate_paths.items()},
        "expected": expected,
        "payloads": payloads,
        "final": final,
        "bridge": bridge,
        "head": final["head"],
        "catalog_before": catalog["before"],
        "catalog_after": catalog["after"],
    }


TRANSACTION_PAYLOADS = {
    "bundle": ("final.bundle", "0644"),
    "evidence": ("evidence.json", "0444"),
    "provenance_review": ("provenance_review.json", "0444"),
    "plasma_review": ("plasma_review.json", "0444"),
    "audit": ("audit.json", "0444"),
    "readme_before": ("README.before", "0644"),
    "sha256sums_before": ("SHA256SUMS.before", "0644"),
    "readme_after": ("README.after", "0644"),
    "sha256sums_after": ("SHA256SUMS.after", "0644"),
}


def purge_forensic_only_prepublication_transaction_root(
    layout: dict[str, Path],
) -> bool:
    """Retire one authenticated forensic-only root before F118 is visible."""

    transactions = layout["transactions"]
    if not os.path.lexists(transactions):
        return False
    visible_targets = [
        layout[f"f118_{key}"] for key in F118_PATHS if layout[f"f118_{key}"].exists()
    ]
    if visible_targets:
        return False
    with bound_directory(
        transactions, "source-authority pre-publication forensic transaction root"
    ) as (descriptor, profile):
        if non_forensic_entry_names(
            descriptor, "source-authority pre-publication forensic transaction root"
        ):
            return False
        require_only_forensic_entries(
            descriptor, "source-authority pre-publication forensic transaction root"
        )
        contents = directory_content_bindings(
            descriptor, "source-authority pre-publication forensic transaction root"
        )
    rmdir_bound_path(
        transactions,
        profile,
        "source-authority pre-publication forensic transaction root",
        expected_contents=contents,
    )
    return True


def purge_abandoned_staging_transaction(layout: dict[str, Path],
                                        args: argparse.Namespace | None = None) -> bool:
    """Fail closed: retained staging requires explicit operator disposition."""

    raise ValueError(
        "source-authority staging cleanup is unsupported; retain recovery debris"
    )


def create_transaction(args: argparse.Namespace, layout: dict[str, Path],
                       plan: dict[str, object],
                       publisher_sha256: str) -> Path:
    transactions = layout["transactions"]
    transaction_id = f"{utc_now().replace(':', '')}-{uuid.uuid4().hex}"
    staging = transactions / f"{transaction_id}{STAGING_TRANSACTION_SUFFIX}"
    with bound_directory(
        layout["accounting"], "accounting directory"
    ) as (accounting_descriptor, _):
        try:
            os.stat(transactions.name, dir_fd=accounting_descriptor, follow_symlinks=False)
        except FileNotFoundError:
            mkdir_bound_exclusive(
                accounting_descriptor,
                transactions.name,
                0o755,
                "source-authority transaction root",
            )
    with bound_directory(
        transactions, "source-authority transaction root"
    ) as (transactions_descriptor, _):
        before = set(os.listdir(transactions_descriptor))
        mkdir_bound_exclusive(
            transactions_descriptor,
            staging.name,
            0o700,
            "source-authority staging directory",
        )
        simulation(args, "during-staging-after-directory")
        payload_bindings = {}
        with bound_child_directory(
            transactions_descriptor,
            staging.name,
            "source-authority staging directory",
            mode=0o700,
        ) as (staging_descriptor, _):
            for index, (key, (name, mode_text)) in enumerate(TRANSACTION_PAYLOADS.items()):
                payload = plan["payloads"][key]
                if not isinstance(payload, bytes):
                    raise ValueError(f"source-authority staged payload is not bytes: {key}")
                write_bound_exclusive(
                    staging_descriptor,
                    name,
                    payload,
                    int(mode_text, 8),
                    f"source-authority staged payload {key}",
                )
                payload_bindings[key] = {
                    "name": name,
                    "sha256": sha256_bytes(payload),
                    "mode": mode_text,
                }
                if index == 0:
                    simulation(args, "during-staging-after-first-payload")
            journal = {
                "schema_version": 1,
                "record_type": "stage-i-current-source-authority-publication-transaction",
                "transaction_id": transaction_id,
                "execution_epoch": EXECUTION_EPOCH,
                "checkpoint": CHECKPOINT,
                "state": "staged",
                "created_utc": utc_now(),
                "publisher": {
                    "revision": plan["head"],
                    "sha256": publisher_sha256,
                },
                "candidate_paths": plan["candidate_paths"],
                "expected": plan["expected"],
                "payloads": payload_bindings,
                "targets": {
                    "bundle": plan["final"]["path"],
                    **{key: relative.as_posix() for key, relative in F118_PATHS.items()},
                    "readme": "source-archives/README.md",
                    "sha256sums": "source-archives/SHA256SUMS",
                },
                "catalog_before": plan["catalog_before"],
                "catalog_after": plan["catalog_after"],
            }
            simulation(args, "during-staging-before-journal")
            write_bound_exclusive(
                staging_descriptor,
                "journal.json",
                canonical_json(journal),
                0o600,
                "source-authority staged journal",
            )
            simulation(args, "during-staging-after-journal")
        if set(os.listdir(transactions_descriptor)) != before | {staging.name}:
            raise ValueError("source-authority transaction root changed during staging")
    return staging


def logical_transaction_id(transaction: Path) -> str:
    """Return the journal identity encoded by a retained transaction pathname."""

    name = transaction.name
    if name.endswith(STAGING_TRANSACTION_SUFFIX):
        name = name.removesuffix(STAGING_TRANSACTION_SUFFIX)
    elif name.endswith(RETIRED_TRANSACTION_SUFFIX):
        name = name.removesuffix(RETIRED_TRANSACTION_SUFFIX)
    if not is_publisher_transaction_id(name):
        raise ValueError("source-authority transaction pathname is malformed")
    return name


def matching_transaction_directories(
    root: Path,
    repository: Path,
    layout: dict[str, Path],
    expected_audit_sha256: str,
) -> list[Path]:
    """Select exact recovery inputs while requiring nonmatching state to be inert."""

    expected_audit = require_sha256(expected_audit_sha256, "expected audit SHA-256")
    return [
        transaction
        for transaction, journal, _ in classified_prepublication_transactions(
            root, repository, layout, expected_audit
        )
        if transaction_expected_audit(journal) == expected_audit
    ]


def recover_atomic_journal(transaction: Path) -> str:
    """Recover one journal while retaining an authenticated deterministic copy."""

    with bound_directory(
        transaction, "source-authority transaction directory", mode=0o700
    ) as (descriptor, _):
        entries = set(os.listdir(descriptor))
        temporaries = sorted(
            name
            for name in entries
            if JOURNAL_TEMP_RE.fullmatch(name) or name in JOURNAL_RECOVERY_NAMES
        )
        journal_record: tuple[bytes, str, os.stat_result] | None = None
        if "journal.json" in entries:
            try:
                _, journal_payload, journal_digest, journal_profile = read_bound_json(
                    descriptor, "journal.json", "source-authority journal", mode=0o600
                )
            except ValueError:
                pass
            else:
                journal_record = (journal_payload, journal_digest, journal_profile)

        valid: list[tuple[str, bytes, str, os.stat_result]] = []
        observed: dict[str, os.stat_result] = {}
        for temporary in temporaries:
            profile = os.stat(temporary, dir_fd=descriptor, follow_symlinks=False)
            observed[temporary] = profile
            require_bound_entry_profile(
                descriptor,
                temporary,
                profile,
                "source-authority recoverable journal temporary",
            )
            require_file_profile(profile, "source-authority recoverable journal temporary")
            try:
                _, payload, digest, recovered = read_bound_json(
                    descriptor,
                    temporary,
                    "source-authority recoverable journal temporary",
                    mode=0o600,
                )
                require_bound_entry_profile(
                    descriptor,
                    temporary,
                    recovered,
                    "source-authority recoverable journal temporary",
                )
                if file_profile_binding(recovered) != file_profile_binding(profile):
                    raise ValueError(
                        "source-authority recoverable journal temporary profile changed"
                    )
            except ValueError:
                continue
            valid.append((temporary, payload, digest, recovered))

        if journal_record is not None:
            journal_payload, journal_digest, _ = journal_record
            retained, _ = prepare_atomic_recovery_copy(
                descriptor,
                JOURNAL_RECOVERY_NAMES,
                journal_payload,
                journal_digest,
                0o600,
                "source-authority journal",
            )
            durably_authenticate_bound_file(
                descriptor,
                "journal.json",
                "source-authority journal",
                expected=journal_digest,
                mode=0o600,
                payload=journal_payload,
            )
        elif not valid:
            for temporary in temporaries:
                profile = observed[temporary]
                retire_bound_recovery_file(
                    descriptor,
                    temporary,
                    profile,
                    "source-authority incomplete recoverable journal temporary",
                )
            raise ValueError(
                "source-authority transaction lacks a complete recoverable journal"
            )
        else:
            distinct_payloads = {payload for _, payload, _, _ in valid}
            if len(distinct_payloads) != 1:
                raise ValueError(
                    "source-authority transaction has ambiguous recoverable journals"
                )
            _, payload, digest, _ = valid[0]
            retained, retained_profile = prepare_atomic_recovery_copy(
                descriptor,
                JOURNAL_RECOVERY_NAMES,
                payload,
                digest,
                0o600,
                "source-authority recoverable journal",
            )
            publication_source: tuple[str, os.stat_result] | None = None
            for name, candidate_payload, candidate_digest, candidate_profile in valid:
                if (
                    name != retained
                    and candidate_payload == payload
                    and candidate_digest == digest
                ):
                    try:
                        candidate_profile = authenticate_bound_file_content(
                            descriptor,
                            name,
                            candidate_profile,
                            "source-authority recoverable journal publication source",
                            (digest, 0o600, payload, 1),
                        )
                    except ValueError:
                        continue
                    publication_source = name, candidate_profile
                    break
            if publication_source is None:
                source_name = next(name for name in JOURNAL_RECOVERY_NAMES if name != retained)
                source_profile = write_bound_exclusive(
                    descriptor,
                    source_name,
                    payload,
                    0o600,
                    "source-authority recoverable journal publication source",
                )
                publication_source = source_name, source_profile
            temporary, profile = publication_source
            if "journal.json" in entries:
                journal_profile = os.stat(
                    "journal.json", dir_fd=descriptor, follow_symlinks=False
                )
                exchange_bound_entries(
                    descriptor,
                    temporary,
                    "journal.json",
                    profile,
                    journal_profile,
                    "source-authority recoverable journal",
                    source_content=(digest, 0o600, payload, 1),
                )
            else:
                publish_bound_file_noreplace(
                    descriptor,
                    temporary,
                    "journal.json",
                    profile,
                    "source-authority recoverable journal temporary",
                    content=(digest, 0o600, payload, 1),
                )
            durably_authenticate_bound_file(
                descriptor,
                retained,
                "source-authority retained journal recovery",
                expected=digest,
                mode=0o600,
                payload=payload,
            )
            durably_authenticate_bound_file(
                descriptor,
                "journal.json",
                "source-authority recovered journal",
                expected=digest,
                mode=0o600,
                payload=payload,
            )

        retained_payload, retained_digest, retained_profile = read_bound_file(
            descriptor,
            retained,
            "source-authority retained journal recovery",
            mode=0o600,
        )
        for orphan in sorted(os.listdir(descriptor)):
            if orphan == retained or (
                JOURNAL_TEMP_RE.fullmatch(orphan) is None
                and orphan not in JOURNAL_RECOVERY_NAMES
            ):
                continue
            orphan_profile = os.stat(orphan, dir_fd=descriptor, follow_symlinks=False)
            retire_bound_recovery_file(
                descriptor,
                orphan,
                orphan_profile,
                "source-authority orphan journal temporary",
            )
            durably_authenticate_bound_file(
                descriptor,
                "journal.json",
                "source-authority recovered journal",
                expected=retained_digest,
                mode=0o600,
                payload=retained_payload,
            )
            retained_profile = durably_authenticate_bound_file(
                descriptor,
                retained,
                "source-authority retained journal recovery",
                expected=retained_digest,
                mode=0o600,
                payload=retained_payload,
            )
        return retained


def read_complete_transaction(
    transaction: Path,
) -> tuple[dict[str, object], dict[str, bytes]]:
    """Read-only authenticate one complete retained source-authority transaction."""

    with bound_directory(
        transaction, "source-authority transaction directory", mode=0o700
    ) as (descriptor, _):
        entries = set(os.listdir(descriptor))
        if "journal.json" not in entries:
            payload_names = {name for name, _ in TRANSACTION_PAYLOADS.values()}
            if not transaction.name.endswith(STAGING_TRANSACTION_SUFFIX):
                raise ValueError("journal-free source-authority transaction is not staging")
            unexpected = entries - payload_names
            if unexpected:
                raise ValueError(
                    "incomplete source-authority staging contains unexpected entries"
                )
            for name in entries:
                try:
                    profile = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
                except OSError as error:
                    raise ValueError(
                        f"incomplete source-authority staging entry {name} is unavailable"
                    ) from error
                require_file_profile(
                    profile,
                    f"incomplete source-authority staging entry {name}",
                    mode=stat.S_IMODE(profile.st_mode),
                )
            raise ValueError(
                "incomplete source-authority staging requires operator disposition"
            )
        try:
            journal_profile = os.stat(
                "journal.json", dir_fd=descriptor, follow_symlinks=False
            )
        except OSError as error:
            raise ValueError("source-authority journal is unavailable") from error
        if stat.S_IMODE(journal_profile.st_mode) == 0o000:
            raise ValueError(
                "incomplete or downgraded source-authority journal requires "
                "operator disposition"
            )
        journal, journal_payload, journal_digest, _ = read_bound_json(
            descriptor, "journal.json", "source-authority journal", mode=0o600
        )
        journal = require_exact_keys(
            journal,
            {
                "schema_version", "record_type", "transaction_id", "execution_epoch",
                "checkpoint", "state", "created_utc", "publisher", "candidate_paths",
                "expected", "payloads", "targets", "catalog_before", "catalog_after",
            },
            "source-authority journal",
        )
        if (
            journal["schema_version"] != 1
            or journal["record_type"]
            != "stage-i-current-source-authority-publication-transaction"
            or journal["transaction_id"] != logical_transaction_id(transaction)
            or journal["execution_epoch"] != EXECUTION_EPOCH
            or journal["checkpoint"] != CHECKPOINT
            or journal["state"] not in {"staged", "installing", "committed"}
        ):
            raise ValueError("source-authority journal identity differs")
        require_utc(journal["created_utc"], "source-authority journal creation timestamp")
        expected = require_exact_keys(
            journal["expected"],
            {"bundle", "evidence", "provenance_review", "plasma_review", "audit"},
            "source-authority journal expected digests",
        )
        for key, digest in expected.items():
            require_sha256(digest, f"source-authority journal expected {key} SHA-256")
        publisher = require_exact_keys(
            journal["publisher"],
            {"revision", "sha256"},
            "source-authority journal publisher",
        )
        require_revision(publisher["revision"], "source-authority journal publisher revision")
        require_sha256(publisher["sha256"], "source-authority journal publisher SHA-256")
        candidate_paths = require_exact_keys(
            journal["candidate_paths"],
            {"bundle", "evidence", "provenance_review", "plasma_review", "audit"},
            "source-authority journal candidate paths",
        )
        for key, value in candidate_paths.items():
            normalized_absolute(
                Path(require_nonempty(value, f"source-authority journal {key} candidate path")),
                f"source-authority journal {key} candidate path",
            )
        targets = require_exact_keys(
            journal["targets"],
            {
                "bundle", "evidence", "provenance_review", "plasma_review",
                "publication_audit", "readme", "sha256sums",
            },
            "source-authority journal targets",
        )
        for key, value in targets.items():
            require_relative(value, f"source-authority journal {key} target")
        catalog_keys = {
            "catalog_before": {
                "readme_sha256", "sha256sums_sha256",
                "bridge_listed_exactly_once",
                "predecessor_current_source_bundle_listed_exactly_once",
                "final_bundle_listed", "corrupt_c7_listed",
                "historical_f115_preserved",
            },
            "catalog_after": {
                "readme_sha256", "sha256sums_sha256",
                "bridge_listed_exactly_once",
                "predecessor_current_source_bundle_listed_exactly_once",
                "final_bundle_listed_exactly_once", "corrupt_c7_listed",
                "historical_f115_preserved", "historical_f116_preserved",
                "all_prior_checksum_entries_preserved",
                "sole_current_source_bundle",
            },
        }
        catalogs = {}
        for generation, keys in catalog_keys.items():
            catalog = require_exact_keys(
                journal[generation],
                keys,
                f"source-authority journal {generation}",
            )
            for key in ("readme_sha256", "sha256sums_sha256"):
                require_sha256(
                    catalog[key], f"source-authority journal {generation} {key} SHA-256"
                )
            catalogs[generation] = catalog
        if (
            catalogs["catalog_before"]["bridge_listed_exactly_once"] is not True
            or catalogs["catalog_before"][
                "predecessor_current_source_bundle_listed_exactly_once"
            ] is not True
            or catalogs["catalog_before"]["final_bundle_listed"] is not False
            or catalogs["catalog_before"]["corrupt_c7_listed"] is not False
            or catalogs["catalog_before"]["historical_f115_preserved"] is not True
            or catalogs["catalog_after"]["bridge_listed_exactly_once"] is not True
            or catalogs["catalog_after"][
                "predecessor_current_source_bundle_listed_exactly_once"
            ] is not True
            or catalogs["catalog_after"]["final_bundle_listed_exactly_once"] is not True
            or catalogs["catalog_after"]["corrupt_c7_listed"] is not False
            or catalogs["catalog_after"]["historical_f115_preserved"] is not True
            or catalogs["catalog_after"]["historical_f116_preserved"] is not True
            or catalogs["catalog_after"]["all_prior_checksum_entries_preserved"] is not True
            or catalogs["catalog_after"]["sole_current_source_bundle"]
            != targets["bundle"]
        ):
            raise ValueError("source-authority journal catalog policy bindings differ")
        payload_bindings = require_exact_keys(
            journal["payloads"], set(TRANSACTION_PAYLOADS),
            "source-authority journal payloads",
        )
        recoveries = entries & set(JOURNAL_RECOVERY_NAMES)
        temporaries = {name for name in entries if JOURNAL_TEMP_RE.fullmatch(name)}
        for temporary in temporaries:
            try:
                profile = os.stat(temporary, dir_fd=descriptor, follow_symlinks=False)
            except OSError as error:
                raise ValueError(
                    f"source-authority retained journal temporary {temporary} is unavailable"
                ) from error
            require_file_profile(
                profile,
                f"source-authority retained journal temporary {temporary}",
                mode=0o600,
            )
        if len(recoveries) > 1 or temporaries:
            raise ValueError(
                "source-authority transaction contains legacy mutable-journal state"
            )
        for recovery in recoveries:
            recovered, _, _ = read_bound_file(
                descriptor,
                recovery,
                "source-authority retained immutable journal copy",
                expected=journal_digest,
                mode=0o600,
            )
            if recovered != journal_payload:
                raise ValueError(
                    "source-authority retained immutable journal copy differs"
                )
        expected_entries = {"journal.json", *recoveries}
        payloads = {}
        for key, (name, mode_text) in TRANSACTION_PAYLOADS.items():
            retained = require_exact_keys(
                payload_bindings[key], {"name", "sha256", "mode"},
                f"source-authority journal payload {key}",
            )
            if retained != {
                "name": name,
                "sha256": retained["sha256"],
                "mode": mode_text,
            }:
                raise ValueError(f"source-authority journal payload {key} binding differs")
            digest = require_sha256(
                retained["sha256"], f"source-authority payload {key} SHA-256"
            )
            payloads[key], observed, _ = read_bound_file(
                descriptor,
                name,
                f"source-authority payload {key}",
                expected=digest,
                mode=int(mode_text, 8),
            )
            if observed != digest:
                raise ValueError(f"source-authority payload {key} checksum differs")
            expected_entries.add(name)
        if set(os.listdir(descriptor)) != expected_entries:
            raise ValueError("source-authority transaction contains unexpected entries")
    for key in ("bundle", "evidence", "provenance_review", "plasma_review", "audit"):
        if sha256_bytes(payloads[key]) != expected[key]:
            raise ValueError(f"source-authority transaction {key} differs from operator binding")
    for generation, payload_generation in (
        ("catalog_before", "before"),
        ("catalog_after", "after"),
    ):
        catalog = journal[generation]
        if not isinstance(catalog, dict):
            raise ValueError(f"source-authority journal {generation} is not an object")
        for catalog_key, payload_key in (
            ("readme_sha256", f"readme_{payload_generation}"),
            ("sha256sums_sha256", f"sha256sums_{payload_generation}"),
        ):
            if sha256_bytes(payloads[payload_key]) != catalog[catalog_key]:
                raise ValueError(
                    "source-authority journal catalog bindings differ from payloads: "
                    f"{generation} {catalog_key}"
                )
    return journal, payloads


def require_committed_ancestor_publisher(
    repository: Path, journal: dict[str, object]
) -> None:
    """Bind one stale journal publisher to exact bytes in a reachable commit."""

    publisher = require_exact_keys(
        journal["publisher"],
        {"revision", "sha256"},
        "source-authority journal publisher",
    )
    revision = require_revision(
        publisher["revision"], "source-authority journal publisher revision"
    )
    digest = require_sha256(
        publisher["sha256"], "source-authority journal publisher SHA-256"
    )
    if git_run(repository, ["cat-file", "-e", f"{revision}^{{commit}}"]).returncode:
        raise ValueError("source-authority stale publisher revision is not committed")
    head = repository_head(repository)
    if git_run(
        repository, ["merge-base", "--is-ancestor", revision, head]
    ).returncode:
        raise ValueError("source-authority stale publisher revision is not an ancestor")
    committed = git_run(
        repository, ["show", f"{revision}:{PUBLISHER_RELATIVE.as_posix()}"]
    )
    if committed.returncode or sha256_bytes(committed.stdout) != digest:
        raise ValueError(
            "source-authority stale publisher digest differs from committed revision"
        )
    tree = git_run(
        repository, ["ls-tree", revision, "--", PUBLISHER_RELATIVE.as_posix()]
    )
    try:
        tree_line = tree.stdout.decode("ascii").strip()
    except UnicodeDecodeError as error:
        raise ValueError("source-authority stale publisher tree entry is not ASCII") from error
    if (
        tree.returncode
        or not tree_line.startswith("100755 blob ")
        or not tree_line.endswith(f"\t{PUBLISHER_RELATIVE.as_posix()}")
    ):
        raise ValueError("source-authority stale publisher is not a committed executable")


def require_names_absent(parent_path: Path, names: set[str], label: str) -> None:
    """Require exact direct-child names to be absent, including dangling symlinks."""

    with bound_directory(parent_path, f"{label} parent") as (parent, _):
        retained = []
        for name in sorted(names):
            name = require_entry_name(name, f"{label} retained name")
            try:
                os.stat(name, dir_fd=parent, follow_symlinks=False)
            except FileNotFoundError:
                continue
            except OSError as error:
                raise ValueError(f"{label} retained name {name} is unavailable") from error
            retained.append(name)
        if retained:
            raise ValueError(f"{label} retains forbidden names: {retained}")


def direct_child_identity_bindings(
    parent: int, label: str
) -> tuple[tuple[str, tuple[int, ...]], ...] | None:
    """Return one identity-bearing direct-child snapshot, or None after a race."""

    require_bound_directory_descriptor(parent)
    try:
        names = sorted(os.listdir(parent))
    except OSError as error:
        raise ValueError(f"{label} namespace is unavailable") from error
    bindings = []
    for name in names:
        try:
            profile = os.stat(name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            return None
        except OSError as error:
            raise ValueError(f"{label} entry {name} is unavailable") from error
        bindings.append(
            (
                name,
                (
                    profile.st_dev,
                    profile.st_ino,
                    profile.st_mode,
                    profile.st_uid,
                    profile.st_gid,
                    profile.st_nlink,
                    profile.st_size,
                    profile.st_mtime_ns,
                    profile.st_ctime_ns,
                ),
            )
        )
    require_bound_directory_descriptor(parent)
    return tuple(bindings)


def recovery_transaction_id(name: str, prefixes: tuple[str, ...]) -> str | None:
    """Extract an exact publisher-generatable transaction ID from a recovery name."""

    if not name.endswith(".tmp"):
        return None
    for prefix in prefixes:
        if name.startswith(prefix):
            transaction_id = name[len(prefix):-len(".tmp")]
            if is_publisher_transaction_id(transaction_id):
                return transaction_id
    return None


def require_transaction_recovery_names_absent(
    parent_path: Path, target: str, label: str, *, catalog: bool = False
) -> None:
    """Reject legacy recovery names for every valid source-authority transaction ID."""

    target = require_entry_name(target, label)
    prefixes = [f".{target}.", f".{target}.single-link."]
    if catalog:
        prefixes.append(f".{target}.predecessor.")
    with bound_directory(parent_path, f"{label} parent") as (parent, _):
        previous = None
        consecutive = 0
        for _ in range(8):
            current = direct_child_identity_bindings(parent, label)
            if current is None:
                previous = None
                consecutive = 0
                continue
            retained = [
                name for name, _ in current
                if recovery_transaction_id(name, tuple(prefixes)) is not None
            ]
            if retained:
                raise ValueError(f"{label} retains forbidden names: {retained}")
            if current == previous:
                consecutive += 1
            else:
                previous = current
                consecutive = 1
            if consecutive >= 2:
                require_bound_directory_descriptor(parent)
                return
        raise ValueError(f"{label} namespace did not stabilize")


def inert_stale_targets(
    root: Path,
    layout: dict[str, Path],
    journal: dict[str, object],
    payloads: dict[str, bytes],
) -> tuple[tuple[Path, str], ...]:
    """Return exact public targets after binding a stale journal to its evidence."""

    targets = require_exact_keys(
        journal["targets"],
        {
            "bundle", "evidence", "provenance_review", "plasma_review",
            "publication_audit", "readme", "sha256sums",
        },
        "source-authority journal targets",
    )
    for key, relative in F118_PATHS.items():
        if targets[key] != relative.as_posix():
            raise ValueError(f"source-authority stale F118 {key} target differs")
    if (
        targets["readme"] != "source-archives/README.md"
        or targets["sha256sums"] != "source-archives/SHA256SUMS"
    ):
        raise ValueError("source-authority stale catalog targets differ")

    evidence = require_exact_keys(
        parse_json_payload(payloads["evidence"], "source-authority stale evidence"),
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "generated_utc", "scope", "predecessor_authorities", "implementation",
            "source_archive_catalog", "authorization", "validation",
            "publication_requirements",
        },
        "source-authority stale evidence",
    )
    if (
        evidence["schema_version"] != 1
        or evidence["record_type"]
        != "stage-i-current-source-authority-supersession-evidence"
        or evidence["checkpoint"] != CHECKPOINT
        or evidence["execution_epoch"] != EXECUTION_EPOCH
    ):
        raise ValueError("source-authority stale evidence identity differs")
    implementation = require_exact_keys(
        evidence["implementation"],
        {
            "publisher", "committed_tools", "intermediate_36140_bundle",
            "predecessor_current_source_bundle", "current_source_bundle",
        },
        "source-authority stale evidence implementation",
    )
    publisher = require_exact_keys(
        journal["publisher"],
        {"revision", "sha256"},
        "source-authority journal publisher",
    )
    if implementation["publisher"] != {
        "path": PUBLISHER_RELATIVE.as_posix(),
        "revision": publisher["revision"],
        "sha256": publisher["sha256"],
        "mode": REQUIRED_TOOLS[PUBLISHER_RELATIVE.as_posix()],
    }:
        raise ValueError("source-authority stale evidence publisher binding differs")
    final = parse_bundle_declaration(
        implementation["current_source_bundle"],
        "source-authority stale current source bundle",
        current=True,
    )
    if targets["bundle"] != final["path"]:
        raise ValueError("source-authority stale final bundle target differs")
    catalog = require_exact_keys(
        evidence["source_archive_catalog"],
        {"before", "after"},
        "source-authority stale evidence catalog",
    )
    if (
        catalog["before"] != journal["catalog_before"]
        or catalog["after"] != journal["catalog_after"]
    ):
        raise ValueError("source-authority stale evidence catalog binding differs")

    bundle = root / require_relative(targets["bundle"], "source-authority stale bundle target")
    if bundle.parent != layout["source_archives"]:
        raise ValueError("source-authority stale bundle target is not a source-archive child")
    return (
        (bundle, "source-authority stale final source bundle"),
        (layout["f118_evidence"], "source-authority stale F118 evidence"),
        (
            layout["f118_provenance_review"],
            "source-authority stale F118 provenance review",
        ),
        (layout["f118_plasma_review"], "source-authority stale F118 plasma review"),
        (layout["f118_publication_audit"], "source-authority stale F118 publication audit"),
    )


def require_inert_stale_transaction(
    root: Path,
    repository: Path,
    layout: dict[str, Path],
    transaction: Path,
    journal: dict[str, object],
    payloads: dict[str, bytes],
) -> None:
    """Require one non-authoritative prior attempt to have no publication state."""

    if not transaction.name.endswith(STAGING_TRANSACTION_SUFFIX):
        raise ValueError("source-authority stale transaction must retain exact .staging name")
    transaction_id = logical_transaction_id(transaction)
    if journal["state"] != "staged":
        raise ValueError("source-authority stale transaction state must be staged")
    require_committed_ancestor_publisher(repository, journal)
    public_targets = inert_stale_targets(root, layout, journal, payloads)
    for target, label in public_targets:
        require_names_absent(target.parent, {target.name}, label)

    catalogs = (
        (
            layout["readme"],
            "readme_sha256",
            "readme_before",
            "source-authority stale source-archive README",
        ),
        (
            layout["sha256sums"],
            "sha256sums_sha256",
            "sha256sums_before",
            "source-authority stale source-archive SHA256SUMS",
        ),
    )
    for target, digest_key, payload_key, label in catalogs:
        retained, digest = read_file(
            target,
            label,
            expected=str(journal["catalog_before"][digest_key]),
            mode=0o644,
        )
        if retained != payloads[payload_key] or digest != journal["catalog_before"][digest_key]:
            raise ValueError(f"{label} differs from exact catalog_before")

    for target, label in (*public_targets, *((item[0], item[3]) for item in catalogs)):
        names = set(private_publication_names(target.name, label))
        names.update(
            {
                f".{target.name}.{transaction_id}.tmp",
                f".{target.name}.single-link.{transaction_id}.tmp",
            }
        )
        if target in {layout["readme"], layout["sha256sums"]}:
            names.add(f".{target.name}.predecessor.{transaction_id}.tmp")
            names.update(atomic_recovery_names(target))
        require_names_absent(target.parent, names, f"{label} recovery namespace")
        require_transaction_recovery_names_absent(
            target.parent,
            target.name,
            f"{label} transaction recovery namespace",
            catalog=target in {layout["readme"], layout["sha256sums"]},
        )


def authenticated_complete_transactions(
    layout: dict[str, Path],
) -> Iterator[tuple[Path, dict[str, object], dict[str, bytes]]]:
    """Read-only authenticate every retained pre-publication transaction."""

    try:
        for transaction in transaction_directories(layout):
            journal, payloads = read_complete_transaction(transaction)
            yield transaction, journal, payloads
    except ValueError as error:
        raise ValueError(
            "incomplete source-authority staging or unsafe retained transaction "
            "requires operator disposition; only complete authenticated prior-attempt "
            f"debris may coexist before publication: {error}"
        ) from error


def require_inert_stale_transactions(
    root: Path, repository: Path, layout: dict[str, Path]
) -> None:
    """Require every retained drafting/prepublication transaction to be inert."""

    for transaction, journal, payloads in authenticated_complete_transactions(layout):
        require_inert_stale_transaction(
            root, repository, layout, transaction, journal, payloads
        )


def transaction_expected_audit(journal: dict[str, object]) -> str:
    """Return the exact audit digest from one authenticated transaction journal."""

    expected = require_exact_keys(
        journal["expected"],
        {"bundle", "evidence", "provenance_review", "plasma_review", "audit"},
        "source-authority journal expected digests",
    )
    return require_sha256(
        expected["audit"], "source-authority retained journal audit SHA-256"
    )


def classified_prepublication_transactions(
    root: Path,
    repository: Path,
    layout: dict[str, Path],
    expected_audit_sha256: str,
) -> Iterator[tuple[Path, dict[str, object], dict[str, bytes]]]:
    """Authenticate all transactions and require every nonmatching one to be inert."""

    expected_audit = require_sha256(expected_audit_sha256, "expected audit SHA-256")
    for transaction, journal, payloads in authenticated_complete_transactions(layout):
        if transaction_expected_audit(journal) != expected_audit:
            require_inert_stale_transaction(
                root, repository, layout, transaction, journal, payloads
            )
        yield transaction, journal, payloads


def load_transaction(root: Path, repository: Path, layout: dict[str, Path],
                     expected_audit_sha256: str
                     ) -> tuple[Path, dict[str, object], dict[str, bytes]]:
    expected_audit = require_sha256(expected_audit_sha256, "expected audit SHA-256")
    retained = None
    matches = 0
    for transaction, journal, payloads in classified_prepublication_transactions(
        root, repository, layout, expected_audit
    ):
        if transaction_expected_audit(journal) != expected_audit:
            continue
        matches += 1
        if matches == 1:
            retained = transaction, journal, payloads
    if matches != 1 or retained is None:
        raise ValueError("source-authority recovery requires exactly one transaction")
    transaction, journal, payloads = retained
    if transaction.name.endswith(RETIRED_TRANSACTION_SUFFIX):
        raise ValueError("source-authority committed transaction retirement requires recovery")
    return transaction, journal, payloads


def update_journal(transaction: Path, journal: dict[str, object], state: str) -> None:
    raise ValueError(
        "source-authority journals are immutable; publication state is the public audit"
    )


def private_publication_names(target: str, label: str) -> tuple[str, ...]:
    """Return the bounded deterministic private slots for one public name."""

    target = require_entry_name(target, label)
    return tuple(
        require_entry_name(
            f".{target}{PRIVATE_PUBLICATION_SUFFIX}.{slot}",
            f"{label} private publication slot",
        )
        for slot in range(PRIVATE_PUBLICATION_SLOT_COUNT)
    )


def unlink_bound_name_lustre(
    parent: int, name: str, expected: os.stat_result, label: str
) -> None:
    """Unlink one exact bound name without renameat2 or replacement semantics."""

    require_active_mutation_lock_bound()
    require_bound_directory_descriptor(parent)
    name = require_entry_name(name, label)
    require_bound_entry_profile(parent, name, expected, label)
    try:
        os.unlink(name, dir_fd=parent)
    except BaseException as error:
        fsync_bound_namespace(parent)
        try:
            retained = os.stat(name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            return
        require_bound_entry_profile(parent, name, expected, label)
        if file_profile_binding(retained) != file_profile_binding(expected):
            raise ValueError(f"{label} changed during unlink") from error
        raise ValueError(f"{label} unlink failed with its exact name retained") from error
    fsync_bound_namespace(parent)
    require_bound_entry_absent(parent, name, label)
    require_bound_directory_descriptor(parent)
    require_active_mutation_lock_bound()


def prepare_private_publication(
    parent: int,
    target: str,
    payload: bytes,
    expected: str,
    mode: int,
    label: str,
) -> tuple[str, os.stat_result]:
    """Return a finalized private inode, recycling only authenticated incomplete slots."""

    absent: list[str] = []
    incomplete: list[tuple[str, os.stat_result]] = []
    complete: tuple[str, os.stat_result] | None = None
    for name in private_publication_names(target, label):
        try:
            observed = os.stat(name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            absent.append(name)
            continue
        retained_mode = stat.S_IMODE(observed.st_mode)
        if observed.st_nlink == 2:
            try:
                target_profile = os.stat(target, dir_fd=parent, follow_symlinks=False)
            except FileNotFoundError as error:
                raise ValueError(
                    f"{label} private publication slot has an external hardlink"
                ) from error
            if profile_identity(target_profile) != profile_identity(observed):
                raise ValueError(f"{label} private publication slot has an external hardlink")
            retained, _, linked = read_bound_file(
                parent, name, f"{label} linked private publication", expected=expected,
                mode=mode, links=2,
            )
            if retained != payload:
                raise ValueError(f"{label} linked private publication payload differs")
            return name, linked
        require_file_profile(
            observed,
            f"{label} private publication slot",
            mode=retained_mode,
        )
        if retained_mode in {0o000, 0o600}:
            # A crashed private writer is bounded single-link debris.  It may
            # be unlinked and recreated, but its inode is never rewritten.
            incomplete.append((name, observed))
            continue
        if retained_mode != mode:
            raise ValueError(
                f"{label} private publication slot mode is {retained_mode:04o}, "
                f"expected {mode:04o}"
            )
        retained, _, authenticated = read_bound_file(
            parent, name, f"{label} private publication", expected=expected, mode=mode
        )
        if retained != payload:
            raise ValueError(f"{label} private publication payload differs")
        if complete is not None:
            raise ValueError(f"{label} has ambiguous complete private publications")
        complete = name, authenticated
    if complete is not None:
        return complete
    if absent:
        name = absent[0]
    elif incomplete:
        name, observed = incomplete[0]
        unlink_bound_name_lustre(
            parent,
            name,
            observed,
            f"{label} incomplete private publication",
        )
    else:
        raise ValueError(f"{label} bounded private publication slots are exhausted")
    profile = write_bound_exclusive(
        parent, name, payload, mode, f"{label} private publication"
    )
    return name, profile


def finish_private_publication(
    parent: int,
    private: str,
    target: str,
    payload: bytes,
    expected: str,
    mode: int,
    label: str,
) -> os.stat_result:
    """Publish and finish one finalized private inode using Lustre hardlinks."""

    private = require_entry_name(private, f"{label} private publication")
    target = require_entry_name(target, label)
    try:
        target_profile = os.stat(target, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError:
        target_profile = None
    if target_profile is None:
        private_profile = durably_authenticate_bound_file(
            parent,
            private,
            f"{label} private publication",
            expected=expected,
            mode=mode,
            payload=payload,
        )
        require_bound_entry_absent(parent, target, label)
        require_active_mutation_lock_bound()
        require_bound_directory_descriptor(parent)
        try:
            os.link(
                private,
                target,
                src_dir_fd=parent,
                dst_dir_fd=parent,
                follow_symlinks=False,
            )
        except FileExistsError as error:
            fsync_bound_namespace(parent)
            raise ValueError(f"{label} target already exists") from error
        except BaseException as error:
            fsync_bound_namespace(parent)
            raise ValueError(f"{label} no-replace hardlink publication is ambiguous") from error
        fsync_bound_namespace(parent)
        target_profile = os.stat(target, dir_fd=parent, follow_symlinks=False)
        private_profile = os.stat(private, dir_fd=parent, follow_symlinks=False)
        if profile_identity(target_profile) != profile_identity(private_profile):
            raise ValueError(f"{label} hardlink publication inode identity differs")
    if target_profile.st_nlink == 1:
        return durably_authenticate_bound_file(
            parent, target, label, expected=expected, mode=mode, payload=payload
        )
    if target_profile.st_nlink != 2:
        raise ValueError(f"{label} hardlink publication has an unsafe link profile")
    private_profile = os.stat(private, dir_fd=parent, follow_symlinks=False)
    if profile_identity(target_profile) != profile_identity(private_profile):
        raise ValueError(f"{label} public target and private publication are unrelated")
    target_payload, _, target_profile = read_bound_file(
        parent, target, label, expected=expected, mode=mode, links=2
    )
    private_payload, _, private_profile = read_bound_file(
        parent,
        private,
        f"{label} linked private publication",
        expected=expected,
        mode=mode,
        links=2,
    )
    if target_payload != payload or private_payload != payload:
        raise ValueError(f"{label} hardlink publication payload differs")
    # The public name is not final until the private name is durably removed.
    # No bytes or permissions are changed after the public name first appears.
    unlink_bound_name_lustre(
        parent, private, private_profile, f"{label} linked private publication"
    )
    return durably_authenticate_bound_file(
        parent, target, label, expected=expected, mode=mode, payload=payload
    )


def remove_legacy_partial_publication(
    parent: int, target: str, observed: os.stat_result, label: str
) -> None:
    """Remove one pre-redesign incomplete public inode under the Stage I lock."""

    retained_mode = stat.S_IMODE(observed.st_mode)
    require_file_profile(observed, label, mode=retained_mode)
    if retained_mode not in {0o000, 0o600}:
        raise ValueError(f"{label} cannot recover from mode {retained_mode:04o}")
    unlink_bound_name_lustre(parent, target, observed, f"{label} legacy partial")


def cleanup_incomplete_private_publication_slots(target: Path, label: str) -> None:
    """Remove only authenticated single-link remnants before the audit commits."""

    with bound_directory(target.parent, f"{label} parent") as (parent, _):
        for name in private_publication_names(target.name, label):
            try:
                observed = os.stat(name, dir_fd=parent, follow_symlinks=False)
            except FileNotFoundError:
                continue
            retained_mode = stat.S_IMODE(observed.st_mode)
            if retained_mode not in {0o000, 0o600}:
                continue
            require_file_profile(
                observed,
                f"{label} incomplete private publication",
                mode=retained_mode,
            )
            unlink_bound_name_lustre(
                parent,
                name,
                observed,
                f"{label} incomplete private publication",
            )


def require_private_publication_slots_absent(target: Path, label: str) -> None:
    """Require that no deterministic private publication names remain."""

    with bound_directory(target.parent, f"{label} parent") as (parent, _):
        retained = []
        for name in private_publication_names(target.name, label):
            try:
                os.stat(name, dir_fd=parent, follow_symlinks=False)
            except FileNotFoundError:
                continue
            retained.append(name)
        if retained:
            raise ValueError(f"{label} retains private publication slots: {retained}")


def ensure_direct_final_file(
    payload: bytes, target: Path, expected: str, mode: int, label: str
) -> None:
    """Publish exact final bytes without ever rewriting a public inode."""

    expected = require_sha256(expected, f"{label} SHA-256")
    if sha256_bytes(payload) != expected:
        raise ValueError(f"{label} payload checksum differs")
    with bound_directory(target.parent, f"{label} parent") as (parent, _):
        try:
            observed = os.stat(target.name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            observed = None
        else:
            try:
                retained, digest, observed = read_bound_file(
                    parent, target.name, label, expected=expected, mode=mode
                )
            except ValueError:
                if observed.st_nlink == 2:
                    private, _ = prepare_private_publication(
                        parent, target.name, payload, expected, mode, label
                    )
                    finish_private_publication(
                        parent, private, target.name, payload, expected, mode, label
                    )
                    return
                remove_legacy_partial_publication(parent, target.name, observed, label)
            else:
                if retained != payload or digest != expected:
                    raise ValueError(f"{label} published bytes differ")
                return
        private, _ = prepare_private_publication(
            parent, target.name, payload, expected, mode, label
        )
        finish_private_publication(
            parent, private, target.name, payload, expected, mode, label
        )


def ensure_direct_catalog(
    target: Path,
    payload: bytes,
    old_sha256: str,
    new_sha256: str,
    label: str,
) -> None:
    """Publish one reviewed catalog generation without rewriting public bytes."""

    old_sha256 = require_sha256(old_sha256, f"{label} predecessor SHA-256")
    new_sha256 = require_sha256(new_sha256, f"{label} F118 SHA-256")
    if sha256_bytes(payload) != new_sha256:
        raise ValueError(f"{label} F118 payload checksum differs")
    with bound_directory(target.parent, f"{label} parent") as (parent, _):
        try:
            observed = os.stat(target.name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            observed = None
        if observed is not None:
            retained_mode = stat.S_IMODE(observed.st_mode)
            if observed.st_nlink == 2:
                private, _ = prepare_private_publication(
                    parent, target.name, payload, new_sha256, 0o644, label
                )
                finish_private_publication(
                    parent, private, target.name, payload, new_sha256, 0o644, label
                )
                return
            if retained_mode in {0o000, 0o600}:
                remove_legacy_partial_publication(parent, target.name, observed, label)
            else:
                retained, digest, observed = read_bound_file(
                    parent, target.name, label, mode=0o644
                )
                if digest == new_sha256 and retained == payload:
                    return
                if digest != old_sha256:
                    raise ValueError(
                        f"{label} is neither reviewed predecessor nor F118 catalog"
                    )
                # Prepare the complete successor before making the predecessor
                # name absent.  A crash can therefore always continue forward.
                private, _ = prepare_private_publication(
                    parent, target.name, payload, new_sha256, 0o644, label
                )
                unlink_bound_name_lustre(
                    parent, target.name, observed, f"{label} reviewed predecessor"
                )
                finish_private_publication(
                    parent, private, target.name, payload, new_sha256, 0o644, label
                )
                return
        private, _ = prepare_private_publication(
            parent, target.name, payload, new_sha256, 0o644, label
        )
        finish_private_publication(
            parent, private, target.name, payload, new_sha256, 0o644, label
        )


def require_recoverable_catalog_state(
    target: Path, old_sha256: str, new_sha256: str, label: str
) -> None:
    """Require one catalog state that the locked forward transaction can finish."""

    old_sha256 = require_sha256(old_sha256, f"{label} predecessor SHA-256")
    new_sha256 = require_sha256(new_sha256, f"{label} F118 SHA-256")
    with bound_directory(target.parent, f"{label} parent") as (parent, _):
        try:
            observed = os.stat(target.name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            return
        retained_mode = stat.S_IMODE(observed.st_mode)
        if observed.st_nlink == 1 and retained_mode in {0o000, 0o600}:
            require_file_profile(observed, f"{label} legacy partial", mode=retained_mode)
            return
        if observed.st_nlink == 1:
            _, digest, _ = read_bound_file(parent, target.name, label, mode=0o644)
            if digest not in {old_sha256, new_sha256}:
                raise ValueError(f"{label} differs from both reviewed transaction states")
            return
        if observed.st_nlink != 2:
            raise ValueError(f"{label} has an unrecoverable link profile")
        target_payload, _, target_profile = read_bound_file(
            parent, target.name, label, expected=new_sha256, mode=0o644, links=2
        )
        for private in private_publication_names(target.name, label):
            try:
                private_profile = os.stat(private, dir_fd=parent, follow_symlinks=False)
            except FileNotFoundError:
                continue
            if profile_identity(private_profile) != profile_identity(target_profile):
                continue
            private_payload, _, _ = read_bound_file(
                parent,
                private,
                f"{label} linked private publication",
                expected=new_sha256,
                mode=0o644,
                links=2,
            )
            if private_payload != target_payload:
                raise ValueError(f"{label} linked private publication payload differs")
            return
        raise ValueError(f"{label} two-link target lacks its deterministic private name")


def publication_audit_committed(
    layout: dict[str, Path], expected_audit_sha256: str
) -> bool:
    """Return whether the exact single-link 0444 F118 commit marker is visible."""

    expected = require_sha256(expected_audit_sha256, "expected audit SHA-256")
    target = layout["f118_publication_audit"]
    with bound_directory(target.parent, "F118 publication-audit parent") as (parent, _):
        try:
            observed = os.stat(target.name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            return False
        retained_mode = stat.S_IMODE(observed.st_mode)
        if observed.st_nlink == 1 and retained_mode in {0o000, 0o600}:
            require_file_profile(
                observed, "F118 publication audit", mode=retained_mode
            )
            return False
        if observed.st_nlink == 2 and retained_mode == 0o444:
            audit_payload, _, audit_profile = read_bound_file(
                parent,
                target.name,
                "F118 linked publication audit",
                expected=expected,
                mode=0o444,
                links=2,
            )
            for private in private_publication_names(
                target.name, "F118 publication audit"
            ):
                try:
                    private_profile = os.stat(
                        private, dir_fd=parent, follow_symlinks=False
                    )
                except FileNotFoundError:
                    continue
                if profile_identity(private_profile) != profile_identity(audit_profile):
                    continue
                retained, _, _ = read_bound_file(
                    parent,
                    private,
                    "F118 linked private publication audit",
                    expected=expected,
                    mode=0o444,
                    links=2,
                )
                if retained != audit_payload:
                    raise ValueError("F118 linked publication audit payload differs")
                return False
            raise ValueError(
                "F118 linked publication audit lacks its deterministic private name"
            )
        if retained_mode != 0o444:
            raise ValueError(
                f"F118 publication audit mode is {retained_mode:04o}, expected 0444"
            )
        durably_authenticate_bound_file(
            parent,
            target.name,
            "F118 publication audit",
            expected=expected,
            mode=0o444,
            payload=None,
        )
        return True


def ensure_single_link_recovery_copy(
    parent: int,
    name: str,
    payload: bytes,
    expected: str,
    mode: int,
    observed: os.stat_result | None,
    label: str,
) -> os.stat_result:
    """Reuse one complete recovery copy or forensically retire an incomplete one."""

    recovered = authenticate_or_retire_recovery_file(
        parent,
        name,
        observed,
        label,
        expected=expected,
        mode=mode,
        payload=payload,
    )
    if recovered is not None:
        return recovered
    return write_bound_exclusive(parent, name, payload, mode, label)


def ensure_installed(payload: bytes, target: Path, expected: str, mode: int,
                     transaction_id: str, label: str) -> None:
    temporary = f".{target.name}.{transaction_id}.tmp"
    single_link_recovery = f".{target.name}.single-link.{transaction_id}.tmp"
    with bound_directory(target.parent, f"{label} parent") as (parent, _):
        try:
            target_profile = os.stat(target.name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            target_profile = None
        try:
            temporary_profile = os.stat(temporary, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            temporary_profile = None
        try:
            recovery_profile = os.stat(
                single_link_recovery, dir_fd=parent, follow_symlinks=False
            )
        except FileNotFoundError:
            recovery_profile = None
        if target_profile is not None:
            if target_profile.st_nlink == 1:
                try:
                    target_profile = durably_authenticate_bound_file(
                        parent,
                        target.name,
                        label,
                        expected=expected,
                        mode=mode,
                        payload=payload,
                    )
                except ValueError:
                    if temporary_profile is None and recovery_profile is None:
                        raise
                    temporary_profile = ensure_single_link_recovery_copy(
                        parent,
                        temporary,
                        payload,
                        expected,
                        mode,
                        temporary_profile,
                        f"{label} repair temporary",
                    )
                    recovery_profile = ensure_single_link_recovery_copy(
                        parent,
                        single_link_recovery,
                        payload,
                        expected,
                        mode,
                        recovery_profile,
                        f"{label} repair recovery",
                    )
                    target_profile = os.stat(
                        target.name, dir_fd=parent, follow_symlinks=False
                    )
                    _, forged_profile = exchange_bound_entries(
                        parent,
                        single_link_recovery,
                        target.name,
                        recovery_profile,
                        target_profile,
                        f"{label} public target repair",
                        source_content=(expected, mode, payload, 1),
                    )
                    durably_authenticate_bound_file(
                        parent,
                        target.name,
                        label,
                        expected=expected,
                        mode=mode,
                        payload=payload,
                    )
                    retire_bound_recovery_file(
                        parent,
                        single_link_recovery,
                        forged_profile,
                        f"{label} forged public predecessor",
                        links=forged_profile.st_nlink,
                    )
                    durably_authenticate_bound_file(
                        parent,
                        target.name,
                        label,
                        expected=expected,
                        mode=mode,
                        payload=payload,
                    )
                    unlink_bound_entry(
                        parent,
                        temporary,
                        temporary_profile,
                        f"{label} repair temporary",
                    )
                    try:
                        durably_authenticate_bound_file(
                            parent,
                            target.name,
                            label,
                            expected=expected,
                            mode=mode,
                            payload=payload,
                        )
                    except ValueError:
                        ensure_single_link_recovery_copy(
                            parent,
                            single_link_recovery,
                            payload,
                            expected,
                            mode,
                            None,
                            f"{label} retained repair recovery",
                        )
                        raise
                    return
                retained_entries = [
                    [temporary, temporary_profile, f"{label} retained temporary"],
                    [
                        single_link_recovery,
                        recovery_profile,
                        f"{label} single-link recovery predecessor",
                    ],
                ]
                for index, (name, profile, retained_label) in enumerate(retained_entries):
                    if profile is None:
                        continue
                    retained_profile = authenticate_or_retire_recovery_file(
                        parent,
                        name,
                        profile,
                        retained_label,
                        expected=expected,
                        mode=mode,
                        payload=payload,
                        links=profile.st_nlink,
                    )
                    retired_profile = profile
                    if retained_profile is not None:
                        unlink_bound_entry(parent, name, retained_profile, retained_label)
                        retired_profile = retained_profile
                    for later in retained_entries[index + 1:]:
                        later_name, later_profile, later_label = later
                        if (
                            later_profile is not None
                            and profile_identity(later_profile)
                            == profile_identity(retired_profile)
                        ):
                            later[1] = require_renamed_bound_entry_profile(
                                parent, later_name, later_profile, later_label
                            )
                try:
                    durably_authenticate_bound_file(
                        parent,
                        target.name,
                        label,
                        expected=expected,
                        mode=mode,
                        payload=payload,
                    )
                except ValueError:
                    ensure_single_link_recovery_copy(
                        parent,
                        single_link_recovery,
                        payload,
                        expected,
                        mode,
                        None,
                        f"{label} retained recovery",
                    )
                    raise
                return
            if target_profile.st_nlink == 2 and temporary_profile is not None:
                if profile_identity(target_profile) != profile_identity(temporary_profile):
                    raise ValueError(f"{label} target and retained temporary are unrelated")
                target_payload, _, target_profile = read_bound_file(
                    parent, target.name, label, expected=expected, mode=mode, links=2
                )
                temporary_payload, _, temporary_profile = read_bound_file(
                    parent, temporary, f"{label} retained temporary",
                    expected=expected, mode=mode, links=2,
                )
                if target_payload != payload or temporary_payload != payload:
                    raise ValueError(f"{label} target or retained temporary payload differs")
                recovery_profile = ensure_single_link_recovery_copy(
                    parent,
                    single_link_recovery,
                    payload,
                    expected,
                    mode,
                    recovery_profile,
                    f"{label} single-link recovery",
                )
                _, predecessor_profile = exchange_bound_entries(
                    parent,
                    single_link_recovery,
                    target.name,
                    recovery_profile,
                    target_profile,
                    f"{label} single-link recovery",
                    source_content=(expected, mode, payload, 1),
                    target_content=(expected, mode, payload, 2),
                )
                _, _, temporary_profile = read_bound_file(
                    parent,
                    temporary,
                    f"{label} retained temporary",
                    expected=expected,
                    mode=mode,
                    links=2,
                )
                unlink_bound_entry(
                    parent, temporary, temporary_profile, f"{label} retained temporary"
                )
                _, _, predecessor_profile = read_bound_file(
                    parent,
                    single_link_recovery,
                    f"{label} single-link recovery predecessor",
                    expected=expected,
                    mode=mode,
                    links=2,
                )
                unlink_bound_entry(
                    parent,
                    single_link_recovery,
                    predecessor_profile,
                    f"{label} single-link recovery predecessor",
                )
                durably_authenticate_bound_file(
                    parent,
                    target.name,
                    label,
                    expected=expected,
                    mode=mode,
                    payload=payload,
                )
                return
            raise ValueError(f"{label} target has an unrecoverable link profile")
        if temporary_profile is not None:
            temporary_profile = authenticate_or_retire_recovery_file(
                parent,
                temporary,
                temporary_profile,
                f"{label} retained temporary",
                expected=expected,
                mode=mode,
                payload=payload,
            )
        if temporary_profile is None:
            temporary_profile = write_bound_exclusive(
                parent, temporary, payload, mode, f"{label} temporary"
            )
        recovery_profile = ensure_single_link_recovery_copy(
            parent,
            single_link_recovery,
            payload,
            expected,
            mode,
            recovery_profile,
            f"{label} publication recovery",
        )
        publish_bound_file_noreplace(
            parent,
            temporary,
            target.name,
            temporary_profile,
            label,
            content=(expected, mode, payload, 1),
        )
        durably_authenticate_bound_file(
            parent,
            target.name,
            label,
            expected=expected,
            mode=mode,
            payload=payload,
        )
        unlink_bound_entry(
            parent,
            single_link_recovery,
            recovery_profile,
            f"{label} publication recovery",
        )
        try:
            durably_authenticate_bound_file(
                parent,
                target.name,
                label,
                expected=expected,
                mode=mode,
                payload=payload,
            )
        except ValueError:
            ensure_single_link_recovery_copy(
                parent,
                single_link_recovery,
                payload,
                expected,
                mode,
                None,
                f"{label} retained publication recovery",
            )
            raise


def restore_catalog_predecessor(
    parent: int,
    target: str,
    recovery: str,
    old_payload: bytes,
    old_sha256: str,
    label: str,
) -> None:
    """Restore the reviewed catalog predecessor after an ambiguous exchange."""

    try:
        durably_authenticate_bound_file(
            parent,
            target,
            f"{label} reviewed predecessor",
            expected=old_sha256,
            mode=0o644,
            payload=old_payload,
        )
        return
    except ValueError:
        pass
    try:
        recovery_profile = os.stat(recovery, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError:
        recovery_profile = None
    recovery_profile = ensure_single_link_recovery_copy(
        parent,
        recovery,
        old_payload,
        old_sha256,
        0o644,
        recovery_profile,
        f"{label} predecessor recovery",
    )
    target_profile = os.stat(target, dir_fd=parent, follow_symlinks=False)
    exchange_bound_entries(
        parent,
        recovery,
        target,
        recovery_profile,
        target_profile,
        f"{label} predecessor rollback",
        source_content=(old_sha256, 0o644, old_payload, 1),
    )
    durably_authenticate_bound_file(
        parent,
        target,
        f"{label} restored predecessor",
        expected=old_sha256,
        mode=0o644,
        payload=old_payload,
    )


def ensure_catalog(target: Path, payload: bytes, old_sha256: str, new_sha256: str,
                   transaction_id: str, label: str) -> None:
    temporary = f".{target.name}.{transaction_id}.tmp"
    predecessor_recovery = f".{target.name}.predecessor.{transaction_id}.tmp"
    with bound_directory(target.parent, f"{label} parent") as (parent, _):
        target_payload, observed, target_profile = read_bound_file(
            parent, target.name, label, mode=0o644
        )
        if observed == new_sha256:
            old_payload: bytes | None = None
            try:
                temporary_profile = os.stat(
                    temporary, dir_fd=parent, follow_symlinks=False
                )
            except FileNotFoundError:
                temporary_profile = None
            try:
                recovery_profile = os.stat(
                    predecessor_recovery, dir_fd=parent, follow_symlinks=False
                )
            except FileNotFoundError:
                recovery_profile = None
            for name, profile, retained_label in (
                (temporary, temporary_profile, f"{label} retained temporary"),
                (
                    predecessor_recovery,
                    recovery_profile,
                    f"{label} retained predecessor recovery",
                ),
            ):
                if profile is None:
                    continue
                retained_profile = authenticate_or_retire_recovery_file(
                    parent,
                    name,
                    profile,
                    retained_label,
                    expected=old_sha256,
                    mode=0o644,
                    payload=None,
                )
                if retained_profile is None:
                    continue
                retained_payload, _, retained_profile = read_bound_file(
                    parent, name, retained_label, expected=old_sha256, mode=0o644
                )
                if file_profile_binding(retained_profile) != file_profile_binding(profile):
                    raise ValueError(f"{retained_label} profile changed before cleanup")
                if old_payload is None:
                    old_payload = retained_payload
                elif old_payload != retained_payload:
                    raise ValueError(f"{label} retained predecessors differ")
                unlink_bound_entry(parent, name, retained_profile, retained_label)
                try:
                    durably_authenticate_bound_file(
                        parent,
                        target.name,
                        label,
                        expected=new_sha256,
                        mode=0o644,
                        payload=payload,
                    )
                except ValueError:
                    if old_payload is None:
                        raise
                    restore_catalog_predecessor(
                        parent,
                        target.name,
                        predecessor_recovery,
                        old_payload,
                        old_sha256,
                        label,
                    )
                    raise
            durably_authenticate_bound_file(
                parent,
                target.name,
                label,
                expected=new_sha256,
                mode=0o644,
                payload=payload,
            )
            return
        if observed != old_sha256:
            try:
                recovery_profile = os.stat(
                    predecessor_recovery, dir_fd=parent, follow_symlinks=False
                )
            except FileNotFoundError:
                recovery_profile = None
            if recovery_profile is None:
                raise ValueError(f"{label} is neither reviewed predecessor nor F118 catalog")
            recovery_profile = authenticate_or_retire_recovery_file(
                parent,
                predecessor_recovery,
                recovery_profile,
                f"{label} predecessor recovery",
                expected=old_sha256,
                mode=0o644,
                payload=None,
            )
            if recovery_profile is None:
                raise ValueError(f"{label} lacks its reviewed predecessor recovery")
            authenticated_recovery = recovery_profile
            old_payload, _, recovery_profile = read_bound_file(
                parent,
                predecessor_recovery,
                f"{label} predecessor recovery",
                expected=old_sha256,
                mode=0o644,
            )
            if file_profile_binding(recovery_profile) != file_profile_binding(
                authenticated_recovery
            ):
                raise ValueError(f"{label} predecessor recovery profile changed")
            restore_catalog_predecessor(
                parent,
                target.name,
                predecessor_recovery,
                old_payload,
                old_sha256,
                label,
            )
            raise ValueError(f"{label} forged public target was rolled back")
        old_payload = target_payload
        try:
            recovery_profile = os.stat(
                predecessor_recovery, dir_fd=parent, follow_symlinks=False
            )
        except FileNotFoundError:
            recovery_profile = None
        recovery_profile = ensure_single_link_recovery_copy(
            parent,
            predecessor_recovery,
            old_payload,
            old_sha256,
            0o644,
            recovery_profile,
            f"{label} predecessor recovery",
        )
        try:
            temporary_profile = os.stat(
                temporary, dir_fd=parent, follow_symlinks=False
            )
        except FileNotFoundError:
            temporary_profile = None
        temporary_profile = authenticate_or_retire_recovery_file(
            parent,
            temporary,
            temporary_profile,
            f"{label} temporary",
            expected=new_sha256,
            mode=0o644,
            payload=payload,
        )
        if temporary_profile is None:
            temporary_profile = write_bound_exclusive(
                parent, temporary, payload, 0o644, f"{label} temporary"
            )
        try:
            _, predecessor_profile = exchange_bound_entries(
                parent,
                temporary,
                target.name,
                temporary_profile,
                target_profile,
                label,
                source_content=(new_sha256, 0o644, payload, 1),
                target_content=(old_sha256, 0o644, old_payload, 1),
            )
            durably_authenticate_bound_file(
                parent,
                target.name,
                label,
                expected=new_sha256,
                mode=0o644,
                payload=payload,
            )
            _, _, predecessor_profile = read_bound_file(
                parent,
                temporary,
                f"{label} predecessor",
                expected=old_sha256,
                mode=0o644,
            )
            unlink_bound_entry(parent, temporary, predecessor_profile, f"{label} predecessor")
            durably_authenticate_bound_file(
                parent,
                target.name,
                label,
                expected=new_sha256,
                mode=0o644,
                payload=payload,
            )
            recovery_profile = os.stat(
                predecessor_recovery, dir_fd=parent, follow_symlinks=False
            )
            recovery_profile = authenticate_bound_file_content(
                parent,
                predecessor_recovery,
                recovery_profile,
                f"{label} predecessor recovery",
                (old_sha256, 0o644, old_payload, 1),
            )
            unlink_bound_entry(
                parent,
                predecessor_recovery,
                recovery_profile,
                f"{label} predecessor recovery",
            )
            durably_authenticate_bound_file(
                parent,
                target.name,
                label,
                expected=new_sha256,
                mode=0o644,
                payload=payload,
            )
        except BaseException:
            restore_catalog_predecessor(
                parent,
                target.name,
                predecessor_recovery,
                old_payload,
                old_sha256,
                label,
            )
            raise


def simulation(args: argparse.Namespace, point: str) -> None:
    if args.simulate_interruption == point:
        raise ValueError(f"simulated interruption {point}")


def validate_transaction_context(root: Path, repository: Path, canonical: bool,
                                 publisher_sha256: str, layout: dict[str, Path],
                                 transaction: Path, journal: dict[str, object],
                                 payloads: dict[str, bytes]) -> dict[str, object]:
    publisher = require_exact_keys(
        journal["publisher"], {"revision", "sha256"}, "source-authority journal publisher"
    )
    if publisher["sha256"] != publisher_sha256:
        raise ValueError("source-authority journal publisher digest differs")
    head = require_revision(publisher["revision"], "source-authority journal publisher revision")
    if repository_head(repository) != head:
        raise ValueError("live repository HEAD moved after F118 transaction staging")
    targets = require_exact_keys(
        journal["targets"],
        {
            "bundle", "evidence", "provenance_review", "plasma_review",
            "publication_audit", "readme", "sha256sums",
        },
        "source-authority journal targets",
    )
    for key, relative in F118_PATHS.items():
        if targets[key] != relative.as_posix():
            raise ValueError(f"source-authority journal F118 {key} target differs")
    if targets["readme"] != "source-archives/README.md" or targets[
        "sha256sums"
    ] != "source-archives/SHA256SUMS":
        raise ValueError("source-authority journal catalog targets differ")
    candidate_paths = require_exact_keys(
        journal["candidate_paths"],
        {"bundle", "evidence", "provenance_review", "plasma_review", "audit"},
        "source-authority journal candidate paths",
    )
    normalized_candidates = {
        key: normalized_absolute(
            Path(require_nonempty(value, f"source-authority {key} candidate path")),
            f"source-authority {key} candidate path",
        )
        for key, value in candidate_paths.items()
    }
    expected = journal["expected"]
    if not isinstance(expected, dict):
        raise ValueError("source-authority journal expected digests must be an object")
    context = validate_candidate_payloads(
        root=root,
        repository=repository,
        canonical=canonical,
        publisher_sha256=publisher_sha256,
        bundle_path=transaction / TRANSACTION_PAYLOADS["bundle"][0],
        bundle_sha256=str(expected["bundle"]),
        evidence_path=normalized_candidates["evidence"],
        expected_bundle_candidate_path=normalized_candidates["bundle"],
        evidence_payload=payloads["evidence"],
        evidence_sha256=str(expected["evidence"]),
        provenance_payload=payloads["provenance_review"],
        provenance_sha256=str(expected["provenance_review"]),
        plasma_payload=payloads["plasma_review"],
        plasma_sha256=str(expected["plasma_review"]),
        audit_payload=payloads["audit"],
        audit_sha256=str(expected["audit"]),
        old_readme=payloads["readme_before"],
        old_sums=payloads["sha256sums_before"],
        new_readme=payloads["readme_after"],
        new_sums=payloads["sha256sums_after"],
    )
    if targets["bundle"] != context["final"]["path"]:
        raise ValueError("source-authority journal final bundle target differs")
    if (
        journal["catalog_before"] != context["catalog"]["before"]
        or journal["catalog_after"] != context["catalog"]["after"]
    ):
        raise ValueError("source-authority journal catalog bindings differ from payloads")
    for target, old_key, new_key, label in (
        (layout["readme"], "readme_sha256", "readme_sha256", "source-archive README"),
        (layout["sha256sums"], "sha256sums_sha256", "sha256sums_sha256",
         "source-archive SHA256SUMS"),
    ):
        require_recoverable_catalog_state(
            target,
            str(journal["catalog_before"][old_key]),
            str(journal["catalog_after"][new_key]),
            label,
        )
    return context


def f118_public_authority_declarations(
    layout: dict[str, Path],
    final_target: Path,
    expected: dict[str, object],
    catalog_after: dict[str, object],
    *,
    include_audit: bool,
) -> list[tuple[Path, str, int, str]]:
    """Return the exact aggregate public F118 generation to retain."""

    declarations = [
        (
            normalized_absolute(final_target, "F118 final source bundle"),
            require_sha256(expected.get("bundle"), "F118 final source bundle SHA-256"),
            0o644,
            "F118 final source bundle",
        ),
        (
            layout["f118_evidence"],
            require_sha256(expected.get("evidence"), "F118 evidence SHA-256"),
            0o444,
            "F118 evidence",
        ),
        (
            layout["f118_provenance_review"],
            require_sha256(
                expected.get("provenance_review"), "F118 provenance review SHA-256"
            ),
            0o444,
            "F118 provenance review",
        ),
        (
            layout["f118_plasma_review"],
            require_sha256(
                expected.get("plasma_review"), "F118 plasma review SHA-256"
            ),
            0o444,
            "F118 plasma review",
        ),
        (
            layout["readme"],
            require_sha256(
                catalog_after.get("readme_sha256"), "F118 source-archive README SHA-256"
            ),
            0o644,
            "source-archive README",
        ),
        (
            layout["sha256sums"],
            require_sha256(
                catalog_after.get("sha256sums_sha256"),
                "F118 source-archive SHA256SUMS SHA-256",
            ),
            0o644,
            "source-archive SHA256SUMS",
        ),
    ]
    if include_audit:
        declarations.append(
            (
                layout["f118_publication_audit"],
                require_sha256(expected.get("audit"), "F118 publication audit SHA-256"),
                0o444,
                "F118 publication audit",
            )
        )
    return declarations


def continue_transaction(args: argparse.Namespace, root: Path, repository: Path,
                         canonical: bool, publisher_sha256: str,
                         layout: dict[str, Path], transaction: Path,
                         journal: dict[str, object], payloads: dict[str, bytes]) -> None:
    require_other_transactions_empty(layout)
    context = validate_transaction_context(
        root, repository, canonical, publisher_sha256, layout, transaction, journal, payloads
    )
    expected = journal["expected"]
    if not isinstance(expected, dict):
        raise ValueError("source-authority journal expected digests must be an object")
    final_target = root / str(context["final"]["path"])
    precommit_targets = (
        (final_target, "F118 final source bundle"),
        (layout["f118_evidence"], "F118 evidence"),
        (layout["f118_provenance_review"], "F118 provenance review"),
        (layout["f118_plasma_review"], "F118 plasma review"),
        (layout["readme"], "source-archive README"),
        (layout["sha256sums"], "source-archive SHA256SUMS"),
        (layout["f118_publication_audit"], "F118 publication audit"),
    )
    if publication_audit_committed(layout, str(expected["audit"])):
        # Once the authority commit marker exists, recovery may only verify and
        # authenticate an already complete generation; it must never repair state
        # underneath a visible authority marker.
        verify_promoted(
            root, repository, canonical, publisher_sha256, layout, str(expected["audit"]),
            allow_transaction=True,
        )
        return
    ensure_direct_final_file(
        payloads["bundle"], final_target, str(expected["bundle"]),
        0o644, "F118 final source bundle",
    )
    simulation(args, "after-bundle")
    for key, label in (
        ("evidence", "F118 evidence"),
        ("provenance_review", "F118 provenance review"),
        ("plasma_review", "F118 plasma review"),
    ):
        ensure_direct_final_file(
            payloads[key], layout[f"f118_{key}"], str(expected[key]), 0o444,
            label,
        )
    simulation(args, "after-artifacts")
    ensure_direct_catalog(
        layout["readme"], payloads["readme_after"],
        str(journal["catalog_before"]["readme_sha256"]),
        str(journal["catalog_after"]["readme_sha256"]),
        "source-archive README",
    )
    simulation(args, "after-readme")
    ensure_direct_catalog(
        layout["sha256sums"], payloads["sha256sums_after"],
        str(journal["catalog_before"]["sha256sums_sha256"]),
        str(journal["catalog_after"]["sha256sums_sha256"]),
        "source-archive SHA256SUMS",
    )
    simulation(args, "after-catalogs")
    if publication_audit_committed(layout, str(expected["audit"])):
        verify_promoted(
            root, repository, canonical, publisher_sha256, layout, str(expected["audit"]),
            allow_transaction=True,
        )
        return
    # Remove bounded writer remnants before the sole commit marker appears.
    # Once the exact audit commits, this transaction performs no more mutations.
    for target, label in precommit_targets:
        cleanup_incomplete_private_publication_slots(target, label)
    for target, label in precommit_targets[:-1]:
        require_private_publication_slots_absent(target, label)
    # Bind the complete precommit public generation through the audit-last
    # mutation. Any drift leaves the visible audit non-authorizing and fails
    # closed; recovery never repairs beneath that marker.
    precommit_declarations = f118_public_authority_declarations(
        layout,
        final_target,
        expected,
        journal["catalog_after"],
        include_audit=False,
    )
    with f118_public_authority_closure(precommit_declarations):
        # This immutable audit is the sole authority commit marker and is always last.
        ensure_direct_final_file(
            payloads["audit"], layout["f118_publication_audit"], str(expected["audit"]),
            0o444, "F118 publication audit",
        )
        require_private_publication_slots_absent(
            layout["f118_publication_audit"], "F118 publication audit"
        )
    simulation(args, "after-audit")
    verify_promoted(
        root, repository, canonical, publisher_sha256, layout, str(expected["audit"]),
        allow_transaction=True,
    )


def purge_retired_transaction(layout: dict[str, Path], transaction: Path) -> None:
    raise ValueError(
        "source-authority transaction retirement is unsupported; retain recovery debris"
    )


def cleanup_transaction(args: argparse.Namespace, layout: dict[str, Path],
                        transaction: Path) -> None:
    raise ValueError(
        "source-authority transaction cleanup is unsupported; retain recovery debris"
    )


def classify_non_authoritative_recovery_debris(layout: dict[str, Path]) -> None:
    """Validate only the container profile of post-commit recovery debris."""

    transactions = layout["transactions"]
    if not os.path.lexists(transactions):
        return
    with bound_directory(
        transactions, "non-authoritative source-authority recovery-debris root"
    ) as (descriptor, _):
        for name in non_forensic_entry_names(
            descriptor, "non-authoritative source-authority recovery-debris root"
        ):
            logical_transaction_id(transactions / name)
            with bound_child_directory(
                descriptor,
                name,
                "non-authoritative source-authority recovery debris",
                mode=0o700,
            ):
                pass


def promoted_payloads(layout: dict[str, Path], expected_audit_sha256: str
                      ) -> tuple[dict[str, bytes], dict[str, str]]:
    audit, audit_payload, audit_sha = read_json(
        layout["f118_publication_audit"], "F118 publication audit",
        expected=expected_audit_sha256, mode=0o444,
    )
    artifact = audit.get("artifact")
    reviews = audit.get("independent_reviews")
    if not isinstance(artifact, dict) or not isinstance(reviews, dict):
        raise ValueError("F118 publication audit lacks exact publication bindings")
    expected = {
        "audit": audit_sha,
        "evidence": require_sha256(artifact.get("sha256"), "F118 published evidence SHA-256"),
        "provenance_review": require_sha256(
            reviews.get("provenance_security", {}).get("sha256")
            if isinstance(reviews.get("provenance_security"), dict) else None,
            "F118 published provenance review SHA-256",
        ),
        "plasma_review": require_sha256(
            reviews.get("plasma_scientific_continuation", {}).get("sha256")
            if isinstance(reviews.get("plasma_scientific_continuation"), dict) else None,
            "F118 published plasma review SHA-256",
        ),
    }
    payloads = {"audit": audit_payload}
    for key in ("evidence", "provenance_review", "plasma_review"):
        payloads[key], _ = read_file(
            layout[f"f118_{key}"], f"published F118 {key}",
            expected=expected[key], mode=0o444,
        )
    return payloads, expected


def promoted_f118_public_authority_lease(
    layout: dict[str, Path], expected_audit_sha256: str
) -> F118PublicAuthorityLease:
    """Acquire the exact seven-member public F118 generation after verification."""

    payloads, expected = promoted_payloads(layout, expected_audit_sha256)
    evidence = parse_json_payload(payloads["evidence"], "published F118 evidence")
    implementation = evidence.get("implementation")
    catalog_after = evidence.get("source_archive_catalog", {}).get("after")
    if not isinstance(implementation, dict) or not isinstance(catalog_after, dict):
        raise ValueError("published F118 authority lease binding is missing")
    final = parse_bundle_declaration(
        implementation.get("current_source_bundle"),
        "published F118 current bundle",
        current=True,
    )
    expected["bundle"] = str(final["sha256"])
    return F118PublicAuthorityLease(
        f118_public_authority_declarations(
            layout,
            layout["root"] / str(final["path"]),
            expected,
            catalog_after,
            include_audit=True,
        )
    )


def verify_promoted(root: Path, repository: Path, canonical: bool, publisher_sha256: str,
                    layout: dict[str, Path], expected_audit_sha256: str, *,
                    allow_transaction: bool = False) -> None:
    del allow_transaction
    require_other_transactions_empty(layout)
    payloads, expected = promoted_payloads(layout, expected_audit_sha256)
    evidence = parse_json_payload(payloads["evidence"], "published F118 evidence")
    implementation = evidence.get("implementation")
    if not isinstance(implementation, dict):
        raise ValueError("published F118 implementation is missing")
    final = parse_bundle_declaration(
        implementation.get("current_source_bundle"), "published F118 current bundle",
        current=True,
    )
    old_catalog = evidence.get("source_archive_catalog", {}).get("before")
    new_catalog = evidence.get("source_archive_catalog", {}).get("after")
    if not isinstance(old_catalog, dict) or not isinstance(new_catalog, dict):
        raise ValueError("published F118 catalog binding is missing")
    readme, readme_sha = read_file(layout["readme"], "published source-archive README", mode=0o644)
    sums, sums_sha = read_file(
        layout["sha256sums"], "published source-archive SHA256SUMS", mode=0o644
    )
    if (
        readme_sha != new_catalog.get("readme_sha256")
        or sums_sha != new_catalog.get("sha256sums_sha256")
    ):
        raise ValueError("published source-archive catalogs differ from F118 authority")
    # Reconstruct the predecessor catalogs by removing the deterministic F118 additions.
    bridge = parse_bundle_declaration(
        implementation.get("intermediate_36140_bundle"), "published F118 bridge bundle",
        current=False,
    )
    predecessor = parse_predecessor_bundle_declaration(
        implementation.get("predecessor_current_source_bundle"),
        "published F118 predecessor current source bundle",
    )
    inserted = catalog_readme_block(
        final_name=str(final["name"]),
        final_sha256=str(final["sha256"]),
        final_revision=str(final["head"]),
        final_subject=str(final["subject"]),
    )
    if readme.count(inserted) != 1:
        raise ValueError("published README does not contain one exact F118 authority block")
    old_readme = readme.replace(inserted, b"", 1)
    appended = f"{final['sha256']}  {final['name']}\n".encode()
    if not sums.endswith(appended):
        raise ValueError("published SHA256SUMS lacks exact terminal F118 entry")
    old_sums = sums[:-len(appended)]
    if (
        sha256_bytes(old_readme) != old_catalog.get("readme_sha256")
        or sha256_bytes(old_sums) != old_catalog.get("sha256sums_sha256")
    ):
        raise ValueError("published F118 predecessor catalogs cannot be reconstructed")
    provenance_preview = parse_json_payload(
        payloads["provenance_review"], "published F118 provenance review"
    )
    reviewed_candidate = require_exact_keys(
        provenance_preview.get("reviewed_candidate"), {"path", "sha256"},
        "published F118 reviewed candidate",
    )
    if reviewed_candidate["sha256"] != expected["evidence"]:
        raise ValueError("published F118 reviewed-candidate digest differs")
    reviewed_evidence_path = normalized_absolute(
        Path(require_nonempty(reviewed_candidate["path"], "published F118 reviewed candidate path")),
        "published F118 reviewed candidate path",
    )
    validate_candidate_payloads(
        root=root,
        repository=repository,
        canonical=canonical,
        publisher_sha256=publisher_sha256,
        bundle_path=root / str(final["path"]),
        bundle_sha256=str(final["sha256"]),
        evidence_path=reviewed_evidence_path,
        evidence_payload=payloads["evidence"],
        evidence_sha256=expected["evidence"],
        provenance_payload=payloads["provenance_review"],
        provenance_sha256=expected["provenance_review"],
        plasma_payload=payloads["plasma_review"],
        plasma_sha256=expected["plasma_review"],
        audit_payload=payloads["audit"],
        audit_sha256=expected["audit"],
        old_readme=old_readme,
        old_sums=old_sums,
        new_readme=readme,
        new_sums=sums,
    )
    expected["bundle"] = str(final["sha256"])
    declarations = f118_public_authority_declarations(
        layout,
        root / str(final["path"]),
        expected,
        new_catalog,
        include_audit=True,
    )
    # Once the exact audit commits F118, transaction contents have no authority
    # over the published generation. Re-establish the complete public closure
    # before classifying debris and retain it through successful return.
    with f118_public_authority_closure(declarations):
        classify_non_authoritative_recovery_debris(layout)


def promote(args: argparse.Namespace, root: Path, repository: Path, canonical: bool,
            publisher_sha256: str, layout: dict[str, Path]) -> None:
    if publication_audit_committed(layout, args.expected_audit_sha256):
        verify_promoted(
            root, repository, canonical, publisher_sha256, layout,
            args.expected_audit_sha256, allow_transaction=True,
        )
        return
    retained = matching_transaction_directories(
        root, repository, layout, args.expected_audit_sha256
    )
    if retained:
        transaction, journal, payloads = load_transaction(
            root, repository, layout, args.expected_audit_sha256
        )
        continue_transaction(
            args, root, repository, canonical, publisher_sha256, layout,
            transaction, journal, payloads,
        )
        return
    plan = initial_plan(args, root, repository, canonical, publisher_sha256, layout)
    transaction = create_transaction(args, layout, plan, publisher_sha256)
    simulation(args, "after-staging")
    loaded_transaction, journal, payloads = load_transaction(
        root, repository, layout, args.expected_audit_sha256
    )
    if loaded_transaction != transaction:
        raise ValueError("new source-authority transaction identity changed")
    continue_transaction(
        args, root, repository, canonical, publisher_sha256, layout,
        transaction, journal, payloads,
    )


def recover(args: argparse.Namespace, root: Path, repository: Path, canonical: bool,
            publisher_sha256: str, layout: dict[str, Path]) -> None:
    if publication_audit_committed(layout, args.expected_audit_sha256):
        verify_promoted(
            root, repository, canonical, publisher_sha256, layout,
            args.expected_audit_sha256, allow_transaction=True,
        )
        return
    transaction, journal, payloads = load_transaction(
        root, repository, layout, args.expected_audit_sha256
    )
    continue_transaction(
        args, root, repository, canonical, publisher_sha256, layout,
        transaction, journal, payloads,
    )


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    command.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    command.add_argument("--allow-local-root", action="store_true")
    command.add_argument("--expected-publisher-sha256", required=True)
    actions = command.add_subparsers(dest="action", required=True)
    evidence_parser = actions.add_parser("draft-evidence", allow_abbrev=False)
    evidence_parser.add_argument("--bundle-candidate", type=Path, required=True)
    evidence_parser.add_argument("--expected-bundle-sha256", required=True)
    evidence_parser.add_argument("--bridge-bundle", type=Path, required=True)
    evidence_parser.add_argument("--expected-bridge-sha256", required=True)
    evidence_parser.add_argument("--bridge-head", required=True)
    evidence_parser.add_argument("--evidence-candidate-output", type=Path, required=True)
    evidence_parser.add_argument("--generated-utc", required=True)
    audit_parser = actions.add_parser("draft-audit", allow_abbrev=False)
    audit_parser.add_argument("--bundle-candidate", type=Path, required=True)
    audit_parser.add_argument("--expected-bundle-sha256", required=True)
    audit_parser.add_argument("--evidence-candidate", type=Path, required=True)
    audit_parser.add_argument("--expected-evidence-sha256", required=True)
    audit_parser.add_argument("--provenance-review-candidate", type=Path, required=True)
    audit_parser.add_argument("--expected-provenance-review-sha256", required=True)
    audit_parser.add_argument("--plasma-review-candidate", type=Path, required=True)
    audit_parser.add_argument("--expected-plasma-review-sha256", required=True)
    audit_parser.add_argument("--audit-candidate-output", type=Path, required=True)
    audit_parser.add_argument("--published-utc", required=True)
    promote_parser = actions.add_parser("promote", allow_abbrev=False)
    promote_parser.add_argument("--bundle-candidate", type=Path, required=True)
    promote_parser.add_argument("--expected-bundle-sha256", required=True)
    promote_parser.add_argument("--evidence-candidate", type=Path, required=True)
    promote_parser.add_argument("--expected-evidence-sha256", required=True)
    promote_parser.add_argument("--provenance-review-candidate", type=Path, required=True)
    promote_parser.add_argument("--expected-provenance-review-sha256", required=True)
    promote_parser.add_argument("--plasma-review-candidate", type=Path, required=True)
    promote_parser.add_argument("--expected-plasma-review-sha256", required=True)
    promote_parser.add_argument("--audit-candidate", type=Path, required=True)
    promote_parser.add_argument("--expected-audit-sha256", required=True)
    promote_parser.add_argument(
        "--simulate-interruption",
        choices=(
            "during-staging-after-directory",
            "during-staging-after-first-payload",
            "during-staging-before-journal",
            "during-staging-after-journal",
            "after-staging",
            "after-bundle",
            "after-artifacts",
            "after-readme",
            "after-catalogs",
            "after-audit",
        ),
    )
    recover_parser = actions.add_parser("recover", allow_abbrev=False)
    recover_parser.add_argument("--expected-audit-sha256", required=True)
    recover_parser.add_argument(
        "--simulate-interruption",
        choices=(
            "after-bundle", "after-artifacts", "after-readme", "after-catalogs",
            "after-audit",
        ),
    )
    verify_parser = actions.add_parser("verify", allow_abbrev=False)
    verify_parser.add_argument("--expected-audit-sha256", required=True)
    return command


def validate_cli_hashes(args: argparse.Namespace) -> None:
    args.expected_publisher_sha256 = require_sha256(
        args.expected_publisher_sha256, "expected publisher SHA-256"
    )
    for key in (
        "expected_bundle_sha256", "expected_evidence_sha256",
        "expected_provenance_review_sha256", "expected_plasma_review_sha256",
        "expected_audit_sha256", "expected_bridge_sha256",
    ):
        if hasattr(args, key):
            setattr(args, key, require_sha256(getattr(args, key), key.replace("_", " ")))
    if getattr(args, "simulate_interruption", None) is not None and not args.allow_local_root:
        raise ValueError("simulation options are restricted to local fixtures")


def reject_duplicate_value_options(argv: list[str]) -> None:
    """Reject argparse's otherwise ambiguous last-value-wins behavior."""

    for option in (
        "--root",
        "--expected-publisher-sha256",
        "--bundle-candidate",
        "--expected-bundle-sha256",
        "--evidence-candidate",
        "--expected-evidence-sha256",
        "--provenance-review-candidate",
        "--expected-provenance-review-sha256",
        "--plasma-review-candidate",
        "--expected-plasma-review-sha256",
        "--audit-candidate",
        "--bridge-bundle",
        "--expected-bridge-sha256",
        "--bridge-head",
        "--evidence-candidate-output",
        "--audit-candidate-output",
        "--generated-utc",
        "--published-utc",
        "--expected-audit-sha256",
        "--simulate-interruption",
    ):
        if any(item.startswith(f"{option}=") for item in argv):
            raise ValueError(f"{option} must use a separate exact value")
        if argv.count(option) > 1:
            raise ValueError(f"{option} must not be supplied more than once")


def main(argv: list[str] | None = None) -> int:
    retained_argv = sys.argv[1:] if argv is None else argv
    source, repository, publisher_sha256 = authenticate_self(retained_argv)
    reject_duplicate_value_options(retained_argv)
    args = parser().parse_args(retained_argv)
    validate_cli_hashes(args)
    if args.expected_publisher_sha256 != publisher_sha256:
        raise ValueError("parsed publisher SHA-256 differs from authenticated self")
    if source != repository / PUBLISHER_RELATIVE:
        raise ValueError("publisher source path differs from the committed repository path")
    requested_root = normalized_absolute(args.root, "campaign root")
    requested_canonical = requested_root == DEFAULT_ROOT
    success_lease = None
    try:
        with bound_canonical_public_namespace(requested_root, requested_canonical):
            root, canonical, layout = validate_root(args)
            require_canonical_repository(repository, canonical)
            with stage_i_lock(layout):
                if args.action == "draft-evidence":
                    result = draft_evidence(
                        args, root, repository, canonical, publisher_sha256, layout
                    )
                elif args.action == "draft-audit":
                    result = draft_audit(
                        args, root, repository, canonical, publisher_sha256, layout
                    )
                elif args.action == "promote":
                    promote(args, root, repository, canonical, publisher_sha256, layout)
                    result = layout["f118_publication_audit"]
                elif args.action == "recover":
                    recover(args, root, repository, canonical, publisher_sha256, layout)
                    result = layout["f118_publication_audit"]
                else:
                    verify_promoted(
                        root, repository, canonical, publisher_sha256, layout,
                        args.expected_audit_sha256,
                    )
                    result = layout["f118_publication_audit"]
                if args.action in {"promote", "recover", "verify"}:
                    success_lease = promoted_f118_public_authority_lease(
                        layout, args.expected_audit_sha256
                    )
        if success_lease is None:
            print(result)
            return 0
        with active_f118_public_authority_lease(success_lease):
            print(result)
            return 0
    finally:
        if success_lease is not None:
            success_lease.close()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        KeyError,
        OSError,
        OverflowError,
        subprocess.SubprocessError,
        TypeError,
        ValueError,
    ) as error:
        print(f"Stage I source-authority publisher failed: {error}", file=sys.stderr)
        raise SystemExit(1)
