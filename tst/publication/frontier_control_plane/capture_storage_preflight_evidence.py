#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Capture authenticated write-sync-remove evidence for the fixed PIC roots."""

from __future__ import annotations

import argparse
from collections.abc import Callable
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import secrets
import stat
import subprocess
import sys
import uuid


AUTHORIZED_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
AUTHORIZED_PROJECT_HOME_ROOT = Path("/ccs/proj/ast207/proj-shared/PIC")
EVIDENCE_PARENT_PARTS = ("policy", "storage_preflight_evidence")
ENTRYPOINT_NAME = "capture_storage_preflight_evidence.py"
RUNNER_NAME = "run_control_plane.py"
SCHEMA_NAME = "storage_preflight.schema.json"
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
    "entrypoint_sha256",
    "git_commit",
    "runner_sha256",
    "schema_sha256",
    "tracked_clean_head_blobs",
}
DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
FILE_READ_FLAGS = os.O_RDONLY | os.O_NOFOLLOW
FILE_CREATE_FLAGS = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW
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


def _probe_root(
    root: Path,
    *,
    role: str,
    probe_id: str,
    payload: bytes,
) -> tuple[dict[str, object], tuple[int, int]]:
    lexical, root_descriptor = _open_absolute_directory(root, label=f"{role} root")
    filename = f".pic-storage-preflight-{probe_id}-{role}"
    descriptor: int | None = None
    created = False
    try:
        root_identity = _require_root_identity(
            root_descriptor, None, label=f"{role} root"
        )
        descriptor = os.open(
            filename,
            FILE_CREATE_FLAGS,
            0o600,
            dir_fd=root_descriptor,
        )
        created = True
        _write_all(descriptor, payload)
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"{role} storage probe did not create a regular file")
        if metadata.st_nlink != 1:
            raise ValueError(f"{role} storage probe file link count is not one")
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
        os.close(descriptor)
        descriptor = None

        os.unlink(filename, dir_fd=root_descriptor)
        created = False
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
        if created:
            try:
                os.unlink(filename, dir_fd=root_descriptor)
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


def _publish_evidence(
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
    descriptor: int | None = None
    created = False
    try:
        descriptor = os.open(
            filename,
            FILE_CREATE_FLAGS,
            0o600,
            dir_fd=parent_descriptor,
        )
        created = True
        _write_all(descriptor, payload)
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o400)
        os.fsync(descriptor)
        os.close(descriptor)
        descriptor = None
        os.fsync(parent_descriptor)
        return lexical.joinpath(*EVIDENCE_PARENT_PARTS, filename)
    except BaseException:
        if descriptor is not None:
            os.close(descriptor)
            descriptor = None
        if created:
            try:
                os.unlink(filename, dir_fd=parent_descriptor)
                os.fsync(parent_descriptor)
            except OSError:
                pass
        raise
    finally:
        if descriptor is not None:
            os.close(descriptor)
        os.close(parent_descriptor)


def _read_published_evidence(
    root: Path,
    *,
    expected_root_identity: tuple[int, int],
    filename: str,
) -> bytes:
    _, parent_descriptor = _open_evidence_parent(
        root,
        expected_root_identity=expected_root_identity,
    )
    descriptor: int | None = None
    try:
        descriptor = os.open(filename, FILE_READ_FLAGS, dir_fd=parent_descriptor)
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ValueError("Published storage-preflight evidence is not a regular file")
        if before.st_mode & 0o222:
            raise ValueError("Published storage-preflight evidence is not read-only")
        payload = _read_all(descriptor)
        after = os.fstat(descriptor)
        if not _same_metadata(before, after) or len(payload) != after.st_size:
            raise ValueError(
                "Published storage-preflight evidence changed during readback"
            )
        return payload
    finally:
        if descriptor is not None:
            os.close(descriptor)
        os.close(parent_descriptor)


def _remove_published_evidence(
    root: Path,
    *,
    expected_root_identity: tuple[int, int],
    filename: str,
    expected_payload: bytes,
) -> None:
    _, parent_descriptor = _open_evidence_parent(
        root,
        expected_root_identity=expected_root_identity,
    )
    try:
        descriptor = os.open(filename, FILE_READ_FLAGS, dir_fd=parent_descriptor)
        try:
            metadata = os.fstat(descriptor)
            if not stat.S_ISREG(metadata.st_mode) or metadata.st_mode & 0o222:
                raise ValueError(
                    "Refusing to remove unexpected storage-preflight artifact"
                )
            if _read_all(descriptor) != expected_payload:
                raise ValueError(
                    "Refusing to remove divergent storage-preflight artifact"
                )
        finally:
            os.close(descriptor)
        os.unlink(filename, dir_fd=parent_descriptor)
        os.fsync(parent_descriptor)
    finally:
        os.close(parent_descriptor)


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
    names = [ENTRYPOINT_NAME, RUNNER_NAME, SCHEMA_NAME]
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
    for key in ["entrypoint_sha256", "runner_sha256", "schema_sha256"]:
        if (
            not isinstance(record[key], str)
            or re.fullmatch(r"[0-9a-f]{64}", record[key]) is None
        ):
            raise ValueError(
                "Storage-preflight source authentication digest is malformed"
            )
    return dict(record)


def capture_storage_preflight_evidence(
    *,
    orion_root: Path = AUTHORIZED_PIC_ROOT,
    project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    probe_id: str | None = None,
    now: Clock = _utc_now,
    token_bytes: TokenBytes = secrets.token_bytes,
    authenticate_source: SourceAuthenticator = _authenticated_source_record,
) -> dict[str, object]:
    """Probe both roots and publish one byte-identical immutable evidence pair."""
    canonical_probe_id = _canonical_probe_id(probe_id)
    source_authentication = _validated_source_record(authenticate_source)
    started_utc = _timestamp(now())
    roots = [
        ("orion_simulation_root", _absolute(orion_root)),
        ("project_home_mirror_root", _absolute(project_home_root)),
    ]
    probes: list[dict[str, object]] = []
    identities: dict[str, tuple[int, int]] = {}
    for role, root in roots:
        payload = token_bytes(PAYLOAD_BYTES)
        if not isinstance(payload, bytes) or len(payload) != PAYLOAD_BYTES:
            raise ValueError("Injected token generator returned malformed probe bytes")
        probe, identity = _probe_root(
            root,
            role=role,
            probe_id=canonical_probe_id,
            payload=payload,
        )
        probes.append(probe)
        identities[role] = identity
    if identities[roots[0][0]] == identities[roots[1][0]]:
        raise ValueError("Orion and Project Home roots must be distinct directories")
    if _validated_source_record(authenticate_source) != source_authentication:
        raise ValueError("Storage-preflight source changed after probing roots")

    completed_utc = _timestamp(now())
    filename = f"{canonical_probe_id}.json"
    evidence_paths = {
        role: root.joinpath(*EVIDENCE_PARENT_PARTS, filename)
        for role, root in roots
    }
    artifact = {
        "completed_utc": completed_utc,
        "method": METHOD,
        "probes": probes,
        "probe_id": canonical_probe_id,
        "publication": {
            "orion_path": str(evidence_paths["orion_simulation_root"]),
            "project_home_path": str(evidence_paths["project_home_mirror_root"]),
        },
        "record_type": RECORD_TYPE,
        "schema_version": 1,
        "source_authentication": source_authentication,
        "started_utc": started_utc,
        "status": "passed",
    }
    artifact_payload = _canonical_json_bytes(artifact)
    artifact_sha256 = _sha256(artifact_payload)
    published: list[tuple[str, Path]] = []
    try:
        for role, root in roots:
            published_path = _publish_evidence(
                root,
                expected_root_identity=identities[role],
                filename=filename,
                payload=artifact_payload,
            )
            if published_path != evidence_paths[role]:
                raise ValueError(
                    "Published storage-preflight path differs from fixed path"
                )
            published.append((role, root))
        for role, root in roots:
            if (
                _read_published_evidence(
                    root,
                    expected_root_identity=identities[role],
                    filename=filename,
                )
                != artifact_payload
            ):
                raise ValueError("Published storage-preflight evidence bytes differ")
        if _validated_source_record(authenticate_source) != source_authentication:
            raise ValueError(
                "Storage-preflight source changed during evidence publication"
            )
    except BaseException:
        for role, root in reversed(published):
            try:
                _remove_published_evidence(
                    root,
                    expected_root_identity=identities[role],
                    filename=filename,
                    expected_payload=artifact_payload,
                )
            except (OSError, ValueError):
                pass
        raise
    return {
        "last_preflight_utc": completed_utc,
        "orion_simulation_root_preflight": {
            "method": METHOD,
            "path": str(roots[0][1]),
            "status": "passed",
        },
        "project_home_preflight": {
            "method": METHOD,
            "path": str(roots[1][1]),
            "status": "passed",
        },
        "storage_preflight_evidence": {
            "orion_path": str(evidence_paths["orion_simulation_root"]),
            "probe_id": canonical_probe_id,
            "project_home_path": str(evidence_paths["project_home_mirror_root"]),
            "sha256": artifact_sha256,
        },
    }


def _parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        description=(
            "Capture authenticated storage preflight evidence for fixed PIC roots."
        )
    )


def main() -> None:
    if not getattr(sys, "_pic_control_plane_bootstrapped", False):
        raise SystemExit("Run source control-plane tools through run_control_plane.py")
    _parser().parse_args()
    print(_canonical_json_bytes(capture_storage_preflight_evidence()).decode(), end="")


if __name__ == "__main__":
    main()
