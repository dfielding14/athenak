#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Capture authenticated write-sync-remove evidence for the fixed PIC roots."""

from __future__ import annotations

import argparse
from collections.abc import Callable
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import secrets
import stat
import subprocess
import sys
import time
import uuid

from control_plane_common import stable_serialization_anchor


AUTHORIZED_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
AUTHORIZED_PROJECT_HOME_ROOT = Path(
    "/autofs/nccs-svm1_proj/ast207/proj-shared/PIC"
)
EVIDENCE_PARENT_PARTS = ("policy", "storage_preflight_evidence")
ENTRYPOINT_NAME = "capture_storage_preflight_evidence.py"
RUNNER_NAME = "run_control_plane.py"
SCHEMA_NAME = "storage_preflight.schema.json"
COMMON_NAME = "control_plane_common.py"
METHOD = "local_create_write_sync_remove_probe"
PAYLOAD_BYTES = 32
RECORD_TYPE = "frontier_pic_storage_preflight_evidence"
TRUSTED_GIT = "/usr/bin/git"
TRUSTED_GIT_OPTIONS = [
    "-c",
    "core.fsmonitor=false",
    "-c",
    "core.hooksPath=/dev/null",
]
OPERATIONS = [
    "create_exclusive",
    "write_all",
    "fsync_file",
    "read_back_exact",
    "unlink",
    "fsync_parent",
    "verify_absent",
]
SOURCE_AUTHENTICATION_KEYS = {
    "common_sha256",
    "entrypoint_sha256",
    "git_commit",
    "runner_sha256",
    "schema_sha256",
    "tracked_clean_head_blobs",
}
ARTIFACT_KEYS = {
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
PROBE_KEYS = {
    "operations",
    "path",
    "payload_bytes",
    "payload_sha256",
    "role",
    "st_dev",
    "st_ino",
    "status",
}
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
GIT_COMMIT_PATTERN = re.compile(r"[0-9a-f]{40}")
EVIDENCE_FILENAME_PATTERN = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}\.json"
)
DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
FILE_READ_FLAGS = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK
FILE_CREATE_FLAGS = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW
ORION_ROLE = "orion_simulation_root"
PROJECT_HOME_ROLE = "project_home_mirror_root"
RECOVERY_ROLES = {ORION_ROLE, PROJECT_HOME_ROLE}
RECOVERY_STAGING_SUFFIX = ".recovery-staging"
COMMITTED_LINK_COUNT_TIMEOUT_SECONDS = 2.0
COMMITTED_LINK_COUNT_INITIAL_RETRY_SECONDS = 0.01
COMMITTED_LINK_COUNT_MAX_RETRY_SECONDS = 0.25
STABLE_FIELDS = (
    "st_dev",
    "st_ino",
    "st_mode",
    "st_size",
    "st_mtime_ns",
    "st_ctime_ns",
)


Clock = Callable[[], datetime]
SourceAuthenticator = Callable[[], dict[str, object]]
TokenBytes = Callable[[int], bytes]


def _canonical_json_bytes(value: object) -> bytes:
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


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _timestamp(now: datetime) -> str:
    if now.tzinfo is None:
        raise ValueError("Injected clock must return a timezone-aware datetime")
    return now.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def _canonical_probe_id(value: str | None) -> str:
    probe_id = str(uuid.uuid4()) if value is None else value
    try:
        parsed = uuid.UUID(probe_id)
    except ValueError as error:
        raise ValueError("Probe ID must be one canonical UUID") from error
    if str(parsed) != probe_id:
        raise ValueError("Probe ID must be one canonical UUID")
    return probe_id


def _expected_git_commit(value: str) -> str:
    if not isinstance(value, str) or GIT_COMMIT_PATTERN.fullmatch(value) is None:
        raise ValueError("Expected source Git commit is malformed")
    return value


def _expected_existing_role(value: str) -> str:
    if value not in RECOVERY_ROLES:
        raise ValueError("Expected existing evidence role is malformed")
    return value


def _absolute(path: Path) -> Path:
    return Path(os.path.abspath(path))


def _validate_component(name: str) -> None:
    if not name or "/" in name or name in {".", ".."} or Path(name).name != name:
        raise ValueError(f"Unsafe path component: {name!r}")


def _open_absolute_directory(path: Path, *, label: str) -> tuple[Path, int]:
    lexical = _absolute(path)
    descriptor = os.open("/", DIRECTORY_FLAGS)
    try:
        for part in lexical.parts[1:]:
            _validate_component(part)
            child = os.open(part, DIRECTORY_FLAGS, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        metadata = os.fstat(descriptor)
        if not stat.S_ISDIR(metadata.st_mode):
            raise ValueError(f"{label} is not a directory: {lexical}")
        return lexical, descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _require_same_directory(
    path: Path, descriptor: int, *, label: str
) -> tuple[int, int]:
    expected = os.fstat(descriptor)
    lexical, current_descriptor = _open_absolute_directory(path, label=label)
    del lexical
    try:
        current = os.fstat(current_descriptor)
        if (current.st_dev, current.st_ino) != (expected.st_dev, expected.st_ino):
            raise ValueError(f"{label} changed during storage-preflight recovery")
    finally:
        os.close(current_descriptor)
    return expected.st_dev, expected.st_ino


@contextmanager
def _serialization_anchor_lock(orion_root: Path):
    anchor = stable_serialization_anchor(orion_root)
    lexical, descriptor = _open_absolute_directory(
        anchor, label="Storage-preflight serialization anchor"
    )
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        try:
            _require_same_directory(
                lexical,
                descriptor,
                label="Storage-preflight serialization anchor",
            )
            yield
            _require_same_directory(
                lexical,
                descriptor,
                label="Storage-preflight serialization anchor",
            )
        finally:
            try:
                _require_same_directory(
                    lexical,
                    descriptor,
                    label="Storage-preflight serialization anchor",
                )
            finally:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
    finally:
        os.close(descriptor)


def _require_root_identity(
    descriptor: int, expected_identity: tuple[int, int] | None, *, label: str
) -> tuple[int, int]:
    metadata = os.fstat(descriptor)
    identity = (metadata.st_dev, metadata.st_ino)
    if expected_identity is not None and identity != expected_identity:
        raise ValueError(f"{label} changed after its storage probe")
    return identity


def _write_all(descriptor: int, payload: bytes) -> None:
    offset = 0
    while offset < len(payload):
        written = os.write(descriptor, payload[offset:])
        if written <= 0:
            raise OSError("Storage probe write did not make progress")
        offset += written


def _read_all(descriptor: int) -> bytes:
    chunks: list[bytes] = []
    while True:
        chunk = os.read(descriptor, 1024 * 1024)
        if not chunk:
            return b"".join(chunks)
        chunks.append(chunk)


def _same_metadata(before: os.stat_result, after: os.stat_result) -> bool:
    return all(getattr(before, field) == getattr(after, field) for field in STABLE_FIELDS)


def _require_same_regular_entry(
    directory_descriptor: int,
    name: str,
    expected: os.stat_result,
    *,
    label: str,
) -> os.stat_result:
    current = os.stat(name, dir_fd=directory_descriptor, follow_symlinks=False)
    if (
        not stat.S_ISREG(current.st_mode)
        or (current.st_dev, current.st_ino) != (expected.st_dev, expected.st_ino)
    ):
        raise ValueError(f"{label} namespace entry changed")
    return current


def _unlink_same_regular_entry_if_present(
    directory_descriptor: int,
    name: str,
    expected: os.stat_result,
) -> bool:
    try:
        _require_same_regular_entry(
            directory_descriptor,
            name,
            expected,
            label="Storage-preflight cleanup",
        )
    except (FileNotFoundError, ValueError):
        return False
    os.unlink(name, dir_fd=directory_descriptor)
    return True


def _probe_root(
    root: Path,
    *,
    expected_root_identity: tuple[int, int],
    role: str,
    probe_id: str,
    payload: bytes,
) -> tuple[dict[str, object], tuple[int, int]]:
    lexical, root_descriptor = _open_absolute_directory(root, label=f"{role} root")
    filename = (
        f".pic-storage-preflight-{probe_id}-{role}-{secrets.token_hex(16)}"
    )
    descriptor: int | None = None
    created_metadata: os.stat_result | None = None
    try:
        root_identity = _require_root_identity(
            root_descriptor,
            expected_root_identity,
            label=f"{role} root",
        )
        descriptor = os.open(
            filename,
            FILE_CREATE_FLAGS,
            0o600,
            dir_fd=root_descriptor,
        )
        created_metadata = os.fstat(descriptor)
        if not stat.S_ISREG(created_metadata.st_mode):
            raise ValueError(f"{role} storage probe did not create a regular file")
        if created_metadata.st_nlink != 1:
            raise ValueError(f"{role} storage probe file link count is not one")
        _write_all(descriptor, payload)
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"{role} storage probe did not create a regular file")
        if metadata.st_nlink != 1:
            raise ValueError(f"{role} storage probe file link count is not one")
        if (
            metadata.st_dev,
            metadata.st_ino,
        ) != (
            created_metadata.st_dev,
            created_metadata.st_ino,
        ):
            raise ValueError(f"{role} storage probe file changed during write")
        os.fsync(descriptor)
        os.close(descriptor)
        descriptor = None

        descriptor = os.open(filename, FILE_READ_FLAGS, dir_fd=root_descriptor)
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ValueError(f"{role} storage probe readback is not a regular file")
        if (before.st_dev, before.st_ino) != (metadata.st_dev, metadata.st_ino):
            raise ValueError(f"{role} storage probe file changed before readback")
        actual = _read_all(descriptor)
        after = os.fstat(descriptor)
        if not _same_metadata(before, after) or len(actual) != after.st_size:
            raise ValueError(f"{role} storage probe file changed during readback")
        if actual != payload:
            raise ValueError(f"{role} storage probe readback differs from written bytes")
        _require_same_regular_entry(
            root_descriptor,
            filename,
            after,
            label=f"{role} storage probe",
        )
        os.unlink(filename, dir_fd=root_descriptor)
        created_metadata = None
        os.fsync(root_descriptor)
        try:
            os.stat(filename, dir_fd=root_descriptor, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise ValueError(f"{role} storage probe file remains after unlink")
        return (
            {
                "operations": OPERATIONS,
                "path": str(lexical),
                "payload_bytes": len(payload),
                "payload_sha256": _sha256(payload),
                "role": role,
                "st_dev": root_identity[0],
                "st_ino": root_identity[1],
                "status": "passed",
            },
            root_identity,
        )
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if created_metadata is not None:
            try:
                if _unlink_same_regular_entry_if_present(
                    root_descriptor,
                    filename,
                    created_metadata,
                ):
                    os.fsync(root_descriptor)
            except OSError:
                pass
        os.close(root_descriptor)


def _open_evidence_parent(
    root: Path, *, expected_root_identity: tuple[int, int]
) -> tuple[Path, int]:
    lexical, descriptor = _open_absolute_directory(root, label="Evidence root")
    try:
        _require_root_identity(descriptor, expected_root_identity, label="Evidence root")
        for part in EVIDENCE_PARENT_PARTS:
            _validate_component(part)
            try:
                os.mkdir(part, 0o700, dir_fd=descriptor)
            except FileExistsError:
                pass
            else:
                os.fsync(descriptor)
            child = os.open(part, DIRECTORY_FLAGS, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        return lexical, descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _open_existing_evidence_parent(
    root: Path, *, expected_root_identity: tuple[int, int]
) -> tuple[Path, int]:
    lexical, descriptor = _open_absolute_directory(root, label="Evidence root")
    try:
        _require_root_identity(descriptor, expected_root_identity, label="Evidence root")
        for part in EVIDENCE_PARENT_PARTS:
            _validate_component(part)
            child = os.open(part, DIRECTORY_FLAGS, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        return lexical, descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _read_published_evidence(
    root: Path,
    *,
    expected_root_identity: tuple[int, int],
    filename: str,
) -> bytes:
    _, parent_descriptor = _open_existing_evidence_parent(
        root,
        expected_root_identity=expected_root_identity,
    )
    descriptor: int | None = None
    try:
        descriptor = os.open(filename, FILE_READ_FLAGS, dir_fd=parent_descriptor)
        before = os.fstat(descriptor)
        _validate_evidence_metadata(before)
        payload = _read_all(descriptor)
        after = os.fstat(descriptor)
        if not _same_metadata(before, after) or len(payload) != after.st_size:
            raise ValueError(
                "Published storage-preflight evidence changed during readback"
            )
        _require_same_regular_entry(
            parent_descriptor,
            filename,
            after,
            label="Published storage-preflight evidence",
        )
        return payload
    finally:
        if descriptor is not None:
            os.close(descriptor)
        os.close(parent_descriptor)


def _validate_evidence_metadata_without_link_count(metadata: os.stat_result) -> None:
    if not stat.S_ISREG(metadata.st_mode):
        raise ValueError("Published storage-preflight evidence is not a regular file")
    if stat.S_IMODE(metadata.st_mode) != 0o400:
        raise ValueError("Published storage-preflight evidence mode is not exactly 0400")
    if metadata.st_uid != os.getuid():
        raise ValueError("Published storage-preflight evidence owner differs")


def _validate_evidence_metadata(metadata: os.stat_result) -> None:
    _validate_evidence_metadata_without_link_count(metadata)
    if metadata.st_nlink != 1:
        raise ValueError("Published storage-preflight evidence link count is not one")


def _read_optional_published_evidence(
    root: Path,
    *,
    expected_root_identity: tuple[int, int],
    filename: str,
) -> bytes | None:
    try:
        _, parent_descriptor = _open_existing_evidence_parent(
            root,
            expected_root_identity=expected_root_identity,
        )
    except FileNotFoundError:
        return None
    descriptor: int | None = None
    try:
        try:
            descriptor = os.open(filename, FILE_READ_FLAGS, dir_fd=parent_descriptor)
        except FileNotFoundError:
            return None
        before = os.fstat(descriptor)
        _validate_evidence_metadata(before)
        payload = _read_all(descriptor)
        after = os.fstat(descriptor)
        if not _same_metadata(before, after) or len(payload) != after.st_size:
            raise ValueError(
                "Published storage-preflight evidence changed during readback"
            )
        _require_same_regular_entry(
            parent_descriptor,
            filename,
            after,
            label="Published storage-preflight evidence",
        )
        return payload
    finally:
        if descriptor is not None:
            os.close(descriptor)
        os.close(parent_descriptor)


def _expected_sha256(value: str) -> str:
    if not isinstance(value, str) or SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError("Expected storage-preflight evidence SHA-256 is malformed")
    return value


def _scan_root_probe_residue(
    root: Path, *, expected_root_identity: tuple[int, int], role: str
) -> None:
    _, descriptor = _open_absolute_directory(root, label=f"{role} root")
    try:
        _require_root_identity(
            descriptor, expected_root_identity, label=f"{role} root"
        )
        residue = sorted(
            name
            for name in os.listdir(descriptor)
            if name.startswith(".pic-storage-preflight-")
        )
        if residue:
            raise ValueError(
                f"{role} root has unresolved storage-preflight probe residue"
            )
    finally:
        os.close(descriptor)


def _evidence_entries(
    root: Path, *, expected_root_identity: tuple[int, int]
) -> set[str]:
    try:
        _, descriptor = _open_existing_evidence_parent(
            root, expected_root_identity=expected_root_identity
        )
    except FileNotFoundError:
        return set()
    try:
        entries = set(os.listdir(descriptor))
    finally:
        os.close(descriptor)
    for name in entries:
        if EVIDENCE_FILENAME_PATTERN.fullmatch(name) is None:
            raise ValueError(
                "Storage-preflight evidence directory has unresolved recovery residue"
            )
    return entries


def _audit_exact_pair_with_pinned_roots(
    *,
    probe_id: str,
    expected_sha256: str,
    source_authentication: dict[str, object],
    roots: list[tuple[str, Path]],
    identities: dict[str, tuple[int, int]],
) -> tuple[str, dict[str, bytes | None], dict[str, object] | None]:
    filename = f"{probe_id}.json"
    entries: dict[str, set[str]] = {}
    for role, root in roots:
        _scan_root_probe_residue(
            root,
            expected_root_identity=identities[role],
            role=role,
        )
        entries[role] = _evidence_entries(
            root, expected_root_identity=identities[role]
        )

    all_evidence_names = set().union(*(entries[role] for role, _ in roots))
    for other_name in sorted(all_evidence_names - {filename}):
        if any(other_name not in entries[role] for role, _ in roots):
            raise ValueError(
                "Storage-preflight evidence directories have another one-sided pair"
            )
        mirrored = [
            _read_published_evidence(
                root,
                expected_root_identity=identities[role],
                filename=other_name,
            )
            for role, root in roots
        ]
        if mirrored[0] != mirrored[1]:
            raise ValueError(
                "Storage-preflight evidence directories have another divergent pair"
            )

    payloads = {
        role: _read_optional_published_evidence(
            root,
            expected_root_identity=identities[role],
            filename=filename,
        )
        for role, root in roots
    }
    present_roles = [role for role, _ in roots if payloads[role] is not None]
    if not present_roles:
        return "absent_both", payloads, None
    if len(present_roles) == 2 and payloads[roots[0][0]] != payloads[roots[1][0]]:
        raise ValueError("Mirrored storage-preflight recovery evidence differs")
    payload = payloads[present_roles[0]]
    assert payload is not None
    artifact = _validate_recovery_artifact(
        payload,
        probe_id=probe_id,
        expected_sha256=expected_sha256,
        source_authentication=source_authentication,
        roots=roots,
    )
    if len(present_roles) == 2:
        state = "valid_identical_pair"
    elif present_roles[0] == ORION_ROLE:
        state = "valid_orion_only"
    else:
        state = "valid_project_home_only"
    return state, payloads, artifact


def _require_committed_evidence_metadata(
    parent: Path,
    parent_descriptor: int,
    filename: str,
    expected: os.stat_result,
) -> os.stat_result:
    deadline = time.monotonic() + COMMITTED_LINK_COUNT_TIMEOUT_SECONDS
    retry_seconds = COMMITTED_LINK_COUNT_INITIAL_RETRY_SECONDS
    retrying = False
    while True:
        _require_same_directory(parent, parent_descriptor, label="Evidence parent")
        committed = _require_same_regular_entry(
            parent_descriptor,
            filename,
            expected,
            label="Published storage-preflight evidence",
        )
        _validate_evidence_metadata_without_link_count(committed)
        if committed.st_nlink == 1:
            if retrying and time.monotonic() >= deadline:
                raise ValueError(
                    "Published storage-preflight evidence link count did not "
                    "converge to one before the commit deadline"
                )
            return committed
        if committed.st_nlink != 2:
            raise ValueError(
                "Published storage-preflight evidence link count is not one"
            )
        remaining_seconds = deadline - time.monotonic()
        if remaining_seconds <= 0:
            raise ValueError(
                "Published storage-preflight evidence link count did not converge "
                "to one after commit"
            )
        retrying = True
        os.fsync(parent_descriptor)
        time.sleep(min(retry_seconds, remaining_seconds))
        retry_seconds = min(
            retry_seconds * 2,
            COMMITTED_LINK_COUNT_MAX_RETRY_SECONDS,
        )


def _publish_staged_evidence(
    root: Path,
    *,
    expected_root_identity: tuple[int, int],
    filename: str,
    payload: bytes,
) -> Path:
    lexical, parent_descriptor = _open_evidence_parent(
        root,
        expected_root_identity=expected_root_identity,
    )
    parent = lexical.joinpath(*EVIDENCE_PARENT_PARTS)
    staging_name = f".{filename}.{secrets.token_hex(16)}{RECOVERY_STAGING_SUFFIX}"
    descriptor: int | None = None
    try:
        _require_same_directory(parent, parent_descriptor, label="Evidence parent")
        descriptor = os.open(
            staging_name,
            FILE_CREATE_FLAGS,
            0o600,
            dir_fd=parent_descriptor,
        )
        _write_all(descriptor, payload)
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o400)
        os.fsync(descriptor)
        staged = os.fstat(descriptor)
        _validate_evidence_metadata(staged)
        _require_same_directory(parent, parent_descriptor, label="Evidence parent")
        _require_same_regular_entry(
            parent_descriptor,
            staging_name,
            staged,
            label="Staged storage-preflight evidence",
        )
        os.link(
            staging_name,
            filename,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        published = _require_same_regular_entry(
            parent_descriptor,
            filename,
            staged,
            label="Published storage-preflight evidence",
        )
        if published.st_nlink != 2:
            raise ValueError(
                "Published storage-preflight evidence link count differs before commit"
            )
        os.fsync(parent_descriptor)
        _require_same_directory(parent, parent_descriptor, label="Evidence parent")
        _require_same_regular_entry(
            parent_descriptor,
            staging_name,
            staged,
            label="Staged storage-preflight evidence",
        )
        os.unlink(staging_name, dir_fd=parent_descriptor)
        os.fsync(parent_descriptor)
        _require_committed_evidence_metadata(
            parent,
            parent_descriptor,
            filename,
            staged,
        )
        return lexical.joinpath(*EVIDENCE_PARENT_PARTS, filename)
    finally:
        if descriptor is not None:
            os.close(descriptor)
        os.close(parent_descriptor)


def _open_pinned_roots(
    roots: list[tuple[str, Path]],
) -> tuple[dict[str, tuple[int, int]], dict[str, int]]:
    identities: dict[str, tuple[int, int]] = {}
    descriptors: dict[str, int] = {}
    try:
        for role, root in roots:
            _, descriptor = _open_absolute_directory(root, label=f"{role} root")
            descriptors[role] = descriptor
            identities[role] = _require_root_identity(
                descriptor, None, label=f"{role} root"
            )
        if identities[roots[0][0]] == identities[roots[1][0]]:
            raise ValueError("Orion and Project Home roots must be distinct directories")
        return identities, descriptors
    except BaseException:
        for descriptor in descriptors.values():
            os.close(descriptor)
        raise


def _require_pinned_roots(
    roots: list[tuple[str, Path]], descriptors: dict[str, int]
) -> None:
    for role, root in roots:
        _require_same_directory(root, descriptors[role], label=f"{role} root")


def _git(*arguments: str) -> list[str]:
    return [TRUSTED_GIT, *TRUSTED_GIT_OPTIONS, *arguments]


def _git_environment() -> dict[str, str]:
    return {
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1",
        "HOME": "/",
        "LANG": "C",
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
    }


def _authenticated_source_record() -> dict[str, object]:
    script_dir = _absolute(Path(__file__)).parent
    environment = _git_environment()
    repository = Path(
        subprocess.check_output(
            _git("-C", str(script_dir), "rev-parse", "--show-toplevel"),
            text=True,
            env=environment,
        ).strip()
    )
    names = [ENTRYPOINT_NAME, RUNNER_NAME, SCHEMA_NAME, COMMON_NAME]
    paths = [str((script_dir / name).relative_to(repository)) for name in names]
    subprocess.run(
        _git("-C", str(repository), "ls-files", "--error-unmatch", "--", *paths),
        check=True,
        stdout=subprocess.DEVNULL,
        env=environment,
    )
    head = subprocess.check_output(
        _git("-C", str(repository), "rev-parse", "HEAD"),
        text=True,
        env=environment,
    ).strip()
    if re.fullmatch(r"[0-9a-f]{40}", head) is None:
        raise ValueError("Authenticated source Git commit is malformed")
    status_command = _git(
        "-C",
        str(repository),
        "status",
        "--porcelain=v1",
        "--untracked-files=all",
        "--",
        *paths,
    )
    if subprocess.check_output(status_command, text=True, env=environment):
        raise ValueError(
            "Storage-preflight source files must be clean tracked HEAD blobs"
        )
    sources = {
        name: subprocess.check_output(
            _git("-C", str(repository), "show", f"{head}:{path}"),
            env=environment,
        )
        for name, path in zip(names, paths)
    }
    if (
        subprocess.check_output(
            _git("-C", str(repository), "rev-parse", "HEAD"),
            text=True,
            env=environment,
        ).strip()
        != head
        or subprocess.check_output(status_command, text=True, env=environment)
    ):
        raise ValueError("Storage-preflight source files changed during authentication")
    return {
        "common_sha256": _sha256(sources[COMMON_NAME]),
        "entrypoint_sha256": _sha256(sources[ENTRYPOINT_NAME]),
        "git_commit": head,
        "runner_sha256": _sha256(sources[RUNNER_NAME]),
        "schema_sha256": _sha256(sources[SCHEMA_NAME]),
        "tracked_clean_head_blobs": True,
    }


def _validated_source_record(
    authenticate_source: SourceAuthenticator,
) -> dict[str, object]:
    record = authenticate_source()
    if (
        set(record) != SOURCE_AUTHENTICATION_KEYS
        or not isinstance(record["git_commit"], str)
        or re.fullmatch(r"[0-9a-f]{40}", record["git_commit"]) is None
        or record["tracked_clean_head_blobs"] is not True
    ):
        raise ValueError("Storage-preflight source authentication record is malformed")
    for key in [
        "common_sha256",
        "entrypoint_sha256",
        "runner_sha256",
        "schema_sha256",
    ]:
        if (
            not isinstance(record[key], str)
            or re.fullmatch(r"[0-9a-f]{64}", record[key]) is None
        ):
            raise ValueError(
                "Storage-preflight source authentication digest is malformed"
            )
    return dict(record)


def _require_expected_source_record(
    authenticate_source: SourceAuthenticator,
    *,
    expected_git_commit: str,
) -> dict[str, object]:
    record = _validated_source_record(authenticate_source)
    if record["git_commit"] != _expected_git_commit(expected_git_commit):
        raise ValueError(
            "Storage-preflight source Git commit differs from expected Git commit"
        )
    return record


def _policy_fragment(
    artifact: dict[str, object],
    *,
    artifact_sha256: str,
    orion_root: Path,
    project_home_root: Path,
) -> dict[str, object]:
    probe_id = str(artifact["probe_id"])
    publication = artifact["publication"]
    assert isinstance(publication, dict)
    return {
        "last_preflight_utc": str(artifact["completed_utc"]),
        "orion_simulation_root_preflight": {
            "method": METHOD,
            "path": str(_absolute(orion_root)),
            "status": "passed",
        },
        "project_home_preflight": {
            "method": METHOD,
            "path": str(_absolute(project_home_root)),
            "status": "passed",
        },
        "storage_preflight_evidence": {
            "orion_path": str(publication["orion_path"]),
            "probe_id": probe_id,
            "project_home_path": str(publication["project_home_path"]),
            "sha256": artifact_sha256,
        },
    }


def _validate_recovery_artifact(
    payload: bytes,
    *,
    probe_id: str,
    expected_sha256: str,
    source_authentication: dict[str, object],
    roots: list[tuple[str, Path]],
) -> dict[str, object]:
    _expected_sha256(expected_sha256)
    if _sha256(payload) != expected_sha256:
        raise ValueError("Storage-preflight recovery evidence digest differs")
    try:
        artifact = json.loads(payload.decode("utf-8"))
        canonical = _canonical_json_bytes(artifact)
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError) as error:
        raise ValueError(
            "Storage-preflight recovery evidence is not canonical JSON"
        ) from error
    if canonical != payload or not isinstance(artifact, dict):
        raise ValueError("Storage-preflight recovery evidence is not canonical JSON")
    timestamps = []
    for key in ["started_utc", "completed_utc"]:
        value = artifact.get(key)
        if not isinstance(value, str) or not value.endswith("Z"):
            raise ValueError("Storage-preflight recovery evidence timestamp differs")
        try:
            parsed = datetime.fromisoformat(value[:-1] + "+00:00")
        except ValueError as error:
            raise ValueError(
                "Storage-preflight recovery evidence timestamp differs"
            ) from error
        if _timestamp(parsed) != value:
            raise ValueError("Storage-preflight recovery evidence timestamp differs")
        timestamps.append(parsed)
    if timestamps[1] < timestamps[0]:
        raise ValueError("Storage-preflight recovery evidence timestamps are reversed")
    expected_paths = {
        role: root.joinpath(*EVIDENCE_PARENT_PARTS, f"{probe_id}.json")
        for role, root in roots
    }
    if (
        set(artifact) != ARTIFACT_KEYS
        or type(artifact.get("schema_version")) is not int
        or artifact.get("schema_version") != 2
        or artifact.get("record_type") != RECORD_TYPE
        or artifact.get("probe_id") != probe_id
        or artifact.get("method") != METHOD
        or artifact.get("status") != "passed"
        or artifact.get("source_authentication") != source_authentication
        or artifact.get("publication")
        != {
            "orion_path": str(expected_paths["orion_simulation_root"]),
            "project_home_path": str(
                expected_paths["project_home_mirror_root"]
            ),
        }
    ):
        raise ValueError("Storage-preflight recovery evidence contract differs")
    probes = artifact.get("probes")
    if not isinstance(probes, list) or len(probes) != len(roots):
        raise ValueError("Storage-preflight recovery evidence probes differ")
    for probe, (role, root) in zip(probes, roots):
        if (
            not isinstance(probe, dict)
            or set(probe) != PROBE_KEYS
            or probe.get("role") != role
            or probe.get("path") != str(root)
            or type(probe.get("st_dev")) is not int
            or probe.get("st_dev") < 0
            or type(probe.get("st_ino")) is not int
            or probe.get("st_ino") <= 0
            or type(probe.get("payload_bytes")) is not int
            or probe.get("payload_bytes") != PAYLOAD_BYTES
            or not isinstance(probe.get("payload_sha256"), str)
            or SHA256_PATTERN.fullmatch(str(probe["payload_sha256"])) is None
            or probe.get("operations") != OPERATIONS
            or probe.get("status") != "passed"
        ):
            raise ValueError("Storage-preflight recovery evidence probes differ")
    return artifact


def capture_storage_preflight_evidence(
    *,
    expected_git_commit: str,
    orion_root: Path = AUTHORIZED_PIC_ROOT,
    project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    probe_id: str | None = None,
    now: Clock = _utc_now,
    token_bytes: TokenBytes = secrets.token_bytes,
    authenticate_source: SourceAuthenticator = _authenticated_source_record,
) -> dict[str, object]:
    """Probe both roots and publish one byte-identical immutable evidence pair."""
    canonical_probe_id = _canonical_probe_id(probe_id)
    source_authentication = _require_expected_source_record(
        authenticate_source,
        expected_git_commit=expected_git_commit,
    )
    started_utc = _timestamp(now())
    roots = [
        (ORION_ROLE, _absolute(orion_root)),
        (PROJECT_HOME_ROLE, _absolute(project_home_root)),
    ]
    filename = f"{canonical_probe_id}.json"
    evidence_paths = {
        role: root.joinpath(*EVIDENCE_PARENT_PARTS, filename)
        for role, root in roots
    }
    with _serialization_anchor_lock(roots[0][1]):
        identities, descriptors = _open_pinned_roots(roots)
        try:
            current_source_authentication = _require_expected_source_record(
                authenticate_source,
                expected_git_commit=expected_git_commit,
            )
            if current_source_authentication != source_authentication:
                raise ValueError("Storage-preflight source changed before probing roots")
            state, _, _ = _audit_exact_pair_with_pinned_roots(
                probe_id=canonical_probe_id,
                expected_sha256="0" * 64,
                source_authentication=source_authentication,
                roots=roots,
                identities=identities,
            )
            if state != "absent_both":
                raise ValueError("Storage-preflight probe ID evidence already exists")

            probes: list[dict[str, object]] = []
            for role, root in roots:
                payload = token_bytes(PAYLOAD_BYTES)
                if not isinstance(payload, bytes) or len(payload) != PAYLOAD_BYTES:
                    raise ValueError(
                        "Injected token generator returned malformed probe bytes"
                    )
                probe, identity = _probe_root(
                    root,
                    expected_root_identity=identities[role],
                    role=role,
                    probe_id=canonical_probe_id,
                    payload=payload,
                )
                if identity != identities[role]:
                    raise ValueError(f"{role} root changed during its storage probe")
                probes.append(probe)
            _require_pinned_roots(roots, descriptors)
            try:
                current_source_authentication = _require_expected_source_record(
                    authenticate_source,
                    expected_git_commit=expected_git_commit,
                )
            except ValueError as error:
                raise ValueError(
                    "Storage-preflight source changed after probing roots"
                ) from error
            if current_source_authentication != source_authentication:
                raise ValueError("Storage-preflight source changed after probing roots")

            completed_utc = _timestamp(now())
            artifact = {
                "completed_utc": completed_utc,
                "method": METHOD,
                "probes": probes,
                "probe_id": canonical_probe_id,
                "publication": {
                    "orion_path": str(evidence_paths[ORION_ROLE]),
                    "project_home_path": str(evidence_paths[PROJECT_HOME_ROLE]),
                },
                "record_type": RECORD_TYPE,
                "schema_version": 2,
                "source_authentication": source_authentication,
                "started_utc": started_utc,
                "status": "passed",
            }
            artifact_payload = _canonical_json_bytes(artifact)
            artifact_sha256 = _sha256(artifact_payload)
            state, _, _ = _audit_exact_pair_with_pinned_roots(
                probe_id=canonical_probe_id,
                expected_sha256=artifact_sha256,
                source_authentication=source_authentication,
                roots=roots,
                identities=identities,
            )
            if state != "absent_both":
                raise ValueError("Storage-preflight probe ID evidence already exists")
            for role, root in roots:
                published_path = _publish_staged_evidence(
                    root,
                    expected_root_identity=identities[role],
                    filename=filename,
                    payload=artifact_payload,
                )
                if published_path != evidence_paths[role]:
                    raise ValueError(
                        "Published storage-preflight path differs from fixed path"
                    )
                _require_pinned_roots(roots, descriptors)
            state, _, final_artifact = _audit_exact_pair_with_pinned_roots(
                probe_id=canonical_probe_id,
                expected_sha256=artifact_sha256,
                source_authentication=source_authentication,
                roots=roots,
                identities=identities,
            )
            if state != "valid_identical_pair" or final_artifact != artifact:
                raise ValueError("Published storage-preflight evidence bytes differ")
            try:
                current_source_authentication = _require_expected_source_record(
                    authenticate_source,
                    expected_git_commit=expected_git_commit,
                )
            except ValueError as error:
                raise ValueError(
                    "Storage-preflight source changed during evidence publication"
                ) from error
            if current_source_authentication != source_authentication:
                raise ValueError(
                    "Storage-preflight source changed during evidence publication"
                )
            _require_pinned_roots(roots, descriptors)
        finally:
            for descriptor in descriptors.values():
                os.close(descriptor)
    return _policy_fragment(
        artifact,
        artifact_sha256=artifact_sha256,
        orion_root=roots[0][1],
        project_home_root=roots[1][1],
    )


def audit_storage_preflight_evidence_pair(
    *,
    probe_id: str,
    expected_sha256: str,
    expected_git_commit: str,
    orion_root: Path = AUTHORIZED_PIC_ROOT,
    project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authenticate_source: SourceAuthenticator = _authenticated_source_record,
) -> dict[str, object]:
    """Read and classify one exact evidence pair without mutating either root."""
    canonical_probe_id = _canonical_probe_id(probe_id)
    canonical_sha256 = _expected_sha256(expected_sha256)
    source_authentication = _require_expected_source_record(
        authenticate_source,
        expected_git_commit=expected_git_commit,
    )
    roots = [
        (ORION_ROLE, _absolute(orion_root)),
        (PROJECT_HOME_ROLE, _absolute(project_home_root)),
    ]
    with _serialization_anchor_lock(roots[0][1]):
        identities, descriptors = _open_pinned_roots(roots)
        try:
            state, _, _ = _audit_exact_pair_with_pinned_roots(
                probe_id=canonical_probe_id,
                expected_sha256=canonical_sha256,
                source_authentication=source_authentication,
                roots=roots,
                identities=identities,
            )
            _require_pinned_roots(roots, descriptors)
        finally:
            for descriptor in descriptors.values():
                os.close(descriptor)
    current_source_authentication = _require_expected_source_record(
        authenticate_source,
        expected_git_commit=expected_git_commit,
    )
    if current_source_authentication != source_authentication:
        raise ValueError("Storage-preflight source changed during evidence audit")
    return {
        "expected_sha256": canonical_sha256,
        "probe_id": canonical_probe_id,
        "state": state,
    }


def recover_storage_preflight_evidence_pair(
    *,
    probe_id: str,
    expected_sha256: str,
    expected_existing_role: str,
    expected_git_commit: str,
    orion_root: Path = AUTHORIZED_PIC_ROOT,
    project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authenticate_source: SourceAuthenticator = _authenticated_source_record,
) -> dict[str, object]:
    """Complete one exact interrupted evidence pair without deleting evidence."""
    canonical_probe_id = _canonical_probe_id(probe_id)
    canonical_sha256 = _expected_sha256(expected_sha256)
    canonical_existing_role = _expected_existing_role(expected_existing_role)
    source_authentication = _require_expected_source_record(
        authenticate_source,
        expected_git_commit=expected_git_commit,
    )
    roots = [
        (ORION_ROLE, _absolute(orion_root)),
        (PROJECT_HOME_ROLE, _absolute(project_home_root)),
    ]
    filename = f"{canonical_probe_id}.json"
    with _serialization_anchor_lock(roots[0][1]):
        identities, descriptors = _open_pinned_roots(roots)
        try:
            state, payloads, artifact = _audit_exact_pair_with_pinned_roots(
                probe_id=canonical_probe_id,
                expected_sha256=canonical_sha256,
                source_authentication=source_authentication,
                roots=roots,
                identities=identities,
            )
            if state == "absent_both":
                raise ValueError("Storage-preflight recovery found no published evidence")
            if state not in {
                "valid_identical_pair",
                f"valid_{'orion' if canonical_existing_role == ORION_ROLE else 'project_home'}_only",
            }:
                raise ValueError(
                    "Storage-preflight recovery existing role differs from expected role"
                )
            assert artifact is not None
            payload = next(
                item for item in payloads.values() if item is not None
            )
            try:
                current_source_authentication = _require_expected_source_record(
                    authenticate_source,
                    expected_git_commit=expected_git_commit,
                )
            except ValueError as error:
                raise ValueError(
                    "Storage-preflight source changed before evidence recovery publication"
                ) from error
            if current_source_authentication != source_authentication:
                raise ValueError(
                    "Storage-preflight source changed before evidence recovery publication"
                )
            _require_pinned_roots(roots, descriptors)
            if state != "valid_identical_pair":
                missing_role, missing_root = next(
                    (role, root)
                    for role, root in roots
                    if payloads[role] is None
                )
                _publish_staged_evidence(
                    missing_root,
                    expected_root_identity=identities[missing_role],
                    filename=filename,
                    payload=payload,
                )
            state, _, final_artifact = _audit_exact_pair_with_pinned_roots(
                probe_id=canonical_probe_id,
                expected_sha256=canonical_sha256,
                source_authentication=source_authentication,
                roots=roots,
                identities=identities,
            )
            if state != "valid_identical_pair" or final_artifact != artifact:
                raise ValueError("Recovered storage-preflight evidence bytes differ")
            _require_pinned_roots(roots, descriptors)
        finally:
            for descriptor in descriptors.values():
                os.close(descriptor)
    try:
        current_source_authentication = _require_expected_source_record(
            authenticate_source,
            expected_git_commit=expected_git_commit,
        )
    except ValueError as error:
        raise ValueError(
            "Storage-preflight source changed during evidence recovery"
        ) from error
    if current_source_authentication != source_authentication:
        raise ValueError("Storage-preflight source changed during evidence recovery")
    return _policy_fragment(
        artifact,
        artifact_sha256=canonical_sha256,
        orion_root=roots[0][1],
        project_home_root=roots[1][1],
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Capture authenticated storage preflight evidence for fixed PIC roots."
        )
    )
    parser.add_argument("--expected-git-commit", required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--audit-exact-pair", action="store_true")
    mode.add_argument("--recover-exact-pair", action="store_true")
    parser.add_argument("--probe-id")
    parser.add_argument("--expected-evidence-sha256")
    parser.add_argument("--expected-existing-role", choices=sorted(RECOVERY_ROLES))
    return parser


def main() -> None:
    if not getattr(sys, "_pic_control_plane_bootstrapped", False):
        raise SystemExit("Run source control-plane tools through run_control_plane.py")
    arguments = _parser().parse_args()
    if arguments.audit_exact_pair:
        if (
            arguments.probe_id is None
            or arguments.expected_evidence_sha256 is None
            or arguments.expected_existing_role is not None
        ):
            raise ValueError(
                "Exact storage-preflight audit requires probe ID and evidence digest only"
            )
        result = audit_storage_preflight_evidence_pair(
            probe_id=arguments.probe_id,
            expected_sha256=arguments.expected_evidence_sha256,
            expected_git_commit=arguments.expected_git_commit,
        )
    elif arguments.recover_exact_pair:
        if arguments.probe_id is None or arguments.expected_evidence_sha256 is None:
            raise ValueError(
                "Exact storage-preflight recovery requires probe ID, evidence digest, "
                "and expected existing role"
            )
        if arguments.expected_existing_role is None:
            raise ValueError(
                "Exact storage-preflight recovery requires probe ID, evidence digest, "
                "and expected existing role"
            )
        result = recover_storage_preflight_evidence_pair(
            probe_id=arguments.probe_id,
            expected_sha256=arguments.expected_evidence_sha256,
            expected_existing_role=arguments.expected_existing_role,
            expected_git_commit=arguments.expected_git_commit,
        )
    else:
        if (
            arguments.probe_id is not None
            or arguments.expected_evidence_sha256 is not None
            or arguments.expected_existing_role is not None
        ):
            raise ValueError(
                "Capture does not accept recovery-only probe or evidence bindings"
            )
        result = capture_storage_preflight_evidence(
            expected_git_commit=arguments.expected_git_commit,
        )
    print(_canonical_json_bytes(result).decode(), end="")


if __name__ == "__main__":
    main()
