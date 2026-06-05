#!/usr/bin/env python3
"""Publish the sealed four-case Q-011 pressure-pilot raw bundle.

Only analyzer-declared raw products enter the aggregate tree.  Trusted-
trampoline inventories, runtime metadata, stderr, error histories and emitted
case descriptors remain in their descriptor-pinned source trees.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import ctypes
import errno
import fcntl
import hashlib
import io
import json
import os
from pathlib import Path
from pathlib import PurePosixPath
import re
import stat
import subprocess
import tarfile
from typing import Any, Mapping, Sequence
import uuid

if __package__:
    from . import analyze_q011_section54_pressure_pilot as pilot
    from . import analyze_q011_section54_pressure_pilot_case as case_verifier
else:
    import analyze_q011_section54_pressure_pilot as pilot
    import analyze_q011_section54_pressure_pilot_case as case_verifier


_ARTIFACT_HELPERS = case_verifier._ARTIFACT_HELPERS
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_FILE_FLAGS = (
    os.O_RDONLY
    | getattr(os, "O_NOFOLLOW", 0)
    | getattr(os, "O_NONBLOCK", 0)
)
MAX_RETAINED_FILE_BYTES = 128 * 1024 * 1024
MAX_JSON_BYTES = 8 * 1024 * 1024
MAX_DIRECTORY_ENTRIES = 256
MAX_TREE_ENTRIES = 1024
MAX_TREE_DEPTH = 16
_READ_CHUNK_BYTES = 1024 * 1024
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT_PATTERN = re.compile(r"[0-9a-f]{40}")
_RENAME_NOREPLACE = 1
_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH
_REGISTERED_EXECUTION_PREREGISTRATION_PATH = (
    pilot.REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_pressure_pilot_registered_execution_preregistration_2026-06-02.json"
)
AUTHORIZED_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
AUTHORIZED_PUBLICATION_ROOT = AUTHORIZED_PIC_ROOT / "publication"
TRUSTED_SOURCE_REPOSITORY = Path("/ccs/home/dfielding/athenak-pic")
GIT_EXECUTABLE = Path("/usr/bin/git")
EXECUTING_SOURCE_ROOT = Path(__file__).resolve().parents[2]
WORKER_SOURCE_SNAPSHOT_ROOT_ENV = "PIC_PRESSURE_PUBLICATION_SOURCE_SNAPSHOT_ROOT"
_PUBLICATION_GUARD_SUFFIX = ".publication-invalid"
_PUBLICATION_GUARD_PAYLOAD = b"receipt publication is not authoritative\n"
_PUBLICATION_SEAL_SUFFIX = ".publication-success"
PUBLICATION_ACCEPTANCE_DIRECTORY = "publication_acceptance"
CONSUMPTION_RULE = (
    "receipt_plus_inode_bound_success_seal_are_required_acceptance_markers_"
    "artifacts_without_both_verified_markers_are_unpublished"
)


class PressurePilotPublicationError(ValueError):
    """Raised when aggregate raw-bundle publication fails closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PressurePilotPublicationError(message)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical_json_bytes(value: Mapping[str, object]) -> bytes:
    try:
        payload = (
            json.dumps(dict(value), indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (RecursionError, TypeError, ValueError) as error:
        raise PressurePilotPublicationError("pressure-pilot manifest is not canonical JSON") from error
    _require(
        len(payload) <= MAX_JSON_BYTES,
        f"pressure-pilot JSON exceeds the {MAX_JSON_BYTES}-byte size limit",
    )
    return payload


def _fsync_descriptor(descriptor: int) -> None:
    os.fsync(descriptor)


def _require_same_directory(path: Path, descriptor: int, label: str) -> None:
    expected = os.fstat(descriptor)
    actual = os.stat(path, follow_symlinks=False)
    _require(
        stat.S_ISDIR(actual.st_mode)
        and (expected.st_dev, expected.st_ino) == (actual.st_dev, actual.st_ino),
        f"{label} changed during publication",
    )


def _directory_identity(descriptor: int) -> dict[str, int]:
    metadata = os.fstat(descriptor)
    _require(stat.S_ISDIR(metadata.st_mode), "retained publication root is not a directory")
    return {"device": metadata.st_dev, "inode": metadata.st_ino}


def _require_directory_identity(value: object, descriptor: int, label: str) -> None:
    _require(
        type(value) is dict
        and set(value) == {"device", "inode"}
        and type(value["device"]) is int
        and value["device"] >= 0
        and type(value["inode"]) is int
        and value["inode"] >= 0
        and value == _directory_identity(descriptor),
        f"{label} identity drifted",
    )


def _require_same_directory_at(
    parent_descriptor: int, name: str, descriptor: int, label: str
) -> None:
    try:
        actual = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    except OSError as error:
        raise PressurePilotPublicationError(f"{label} is unavailable") from error
    expected = os.fstat(descriptor)
    _require(
        stat.S_ISDIR(actual.st_mode)
        and (expected.st_dev, expected.st_ino) == (actual.st_dev, actual.st_ino),
        f"{label} changed during publication",
    )


def _require_absent_at(parent_descriptor: int, name: str, label: str) -> None:
    try:
        os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        return
    except OSError as error:
        raise PressurePilotPublicationError(f"{label} cannot be inspected") from error
    raise PressurePilotPublicationError(f"{label} already exists")


def _file_identity_at(parent_descriptor: int, name: str, label: str) -> tuple[int, int]:
    metadata = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    _require(stat.S_ISREG(metadata.st_mode), f"{label} is not a regular file")
    return metadata.st_dev, metadata.st_ino


def _require_same_file_at(
    parent_descriptor: int, name: str, identity: tuple[int, int], label: str
) -> None:
    try:
        actual = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    except OSError as error:
        raise PressurePilotPublicationError(f"{label} is unavailable") from error
    _require(
        stat.S_ISREG(actual.st_mode) and (actual.st_dev, actual.st_ino) == identity,
        f"{label} changed during publication",
    )


def _canonical_existing_directory(path: Path, label: str) -> Path:
    absolute = Path(os.path.abspath(path))
    try:
        resolved = absolute.resolve(strict=True)
    except OSError as error:
        raise PressurePilotPublicationError(f"{label} is unavailable") from error
    _require(resolved == absolute and resolved.is_dir(), f"{label} is not canonical")
    return resolved


def _publication_root(authorized_pic_root: Path) -> tuple[Path, Path]:
    pic_root = _canonical_existing_directory(authorized_pic_root, "authorized PIC root")
    publication_root = _canonical_existing_directory(
        pic_root / "publication", "authorized PIC publication root"
    )
    return pic_root, publication_root


def _publication_acceptance_root(pic_root: Path) -> Path:
    return _canonical_existing_directory(
        pic_root / PUBLICATION_ACCEPTANCE_DIRECTORY,
        "authorized PIC publication acceptance root",
    )


def _publication_transaction_anchor(pic_root: Path) -> Path:
    production_root = Path(os.path.abspath(AUTHORIZED_PIC_ROOT))
    anchor = production_root.parents[2] if pic_root == production_root else pic_root.parent
    return _canonical_existing_directory(anchor, "stable publication transaction anchor")


def _is_production_pic_root(pic_root: Path) -> bool:
    return pic_root == Path(os.path.abspath(AUTHORIZED_PIC_ROOT))


def _direct_publication_target(path: str | Path, publication_root: Path, label: str) -> Path:
    target = Path(os.path.abspath(path))
    _require(
        target.parent == publication_root and target.name not in {"", ".", ".."},
        f"{label} is outside the authorized PIC publication root",
    )
    return target


def _authorized_raw_case_root(path: str | Path, pic_root: Path) -> Path:
    runs_root = _canonical_existing_directory(pic_root / "runs", "authorized PIC runs root")
    raw_root = _canonical_existing_directory(Path(path), "pressure-pilot raw-case root")
    _require(
        raw_root != runs_root and runs_root in raw_root.parents,
        "pressure-pilot raw-case root is outside the authorized PIC runs root",
    )
    return raw_root


def _write_exclusive_at(parent_descriptor: int, name: str, payload: bytes) -> None:
    _require("/" not in name, "pressure-pilot descriptor-relative write received a nested path")
    descriptor = os.open(
        name,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o600,
        dir_fd=parent_descriptor,
    )
    try:
        with os.fdopen(descriptor, "wb", closefd=False) as stream:
            stream.write(payload)
            stream.flush()
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _publication_guard_name(receipt_name: str) -> str:
    _require(
        "/" not in receipt_name and receipt_name not in {"", ".", ".."},
        "pressure-pilot receipt guard received an invalid name",
    )
    return f".{receipt_name}{_PUBLICATION_GUARD_SUFFIX}"


def _require_publication_guard_absent_at(
    parent_descriptor: int, receipt_name: str, label: str
) -> None:
    _require_absent_at(
        parent_descriptor,
        _publication_guard_name(receipt_name),
        f"{label} fail-closed guard",
    )


def _validated_publication_guard_identity_at(
    parent_descriptor: int, receipt_name: str
) -> tuple[int, int]:
    guard_name = _publication_guard_name(receipt_name)
    identity = _file_identity_at(
        parent_descriptor, guard_name, "receipt publication guard"
    )
    _require(
        _read_stable_readonly_regular_at(
            parent_descriptor,
            guard_name,
            "receipt publication guard",
            max_bytes=len(_PUBLICATION_GUARD_PAYLOAD),
        )
        == _PUBLICATION_GUARD_PAYLOAD,
        "receipt publication guard payload drifted",
    )
    _require_same_file_at(
        parent_descriptor, guard_name, identity, "receipt publication guard"
    )
    return identity


def _arm_publication_guard_at(
    parent_descriptor: int, receipt_name: str
) -> tuple[int, int]:
    guard_name = _publication_guard_name(receipt_name)
    _require_absent_at(parent_descriptor, guard_name, "receipt publication guard")
    _write_exclusive_at(parent_descriptor, guard_name, _PUBLICATION_GUARD_PAYLOAD)
    _fsync_descriptor(parent_descriptor)
    return _validated_publication_guard_identity_at(parent_descriptor, receipt_name)


def _ensure_publication_guard_at(
    parent_descriptor: int, receipt_name: str
) -> tuple[int, int]:
    guard_name = _publication_guard_name(receipt_name)
    try:
        os.stat(guard_name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        _write_exclusive_at(parent_descriptor, guard_name, _PUBLICATION_GUARD_PAYLOAD)
    _fsync_descriptor(parent_descriptor)
    return _validated_publication_guard_identity_at(parent_descriptor, receipt_name)


def _disarm_publication_guard_at(
    parent_descriptor: int, receipt_name: str, identity: tuple[int, int]
) -> None:
    guard_name = _publication_guard_name(receipt_name)
    _require_same_file_at(
        parent_descriptor, guard_name, identity, "receipt publication guard"
    )
    os.unlink(guard_name, dir_fd=parent_descriptor)
    _fsync_descriptor(parent_descriptor)
    _require_absent_at(parent_descriptor, guard_name, "receipt publication guard")


def _publication_seal_name(receipt_name: str) -> str:
    _require(
        "/" not in receipt_name and receipt_name not in {"", ".", ".."},
        "pressure-pilot receipt seal received an invalid name",
    )
    return f".{receipt_name}{_PUBLICATION_SEAL_SUFFIX}"


def _publication_seal_payload(
    publication_descriptor: int,
    receipt_name: str,
    receipt_sha256: str,
    receipt_identity: tuple[int, int],
) -> bytes:
    _require(
        _SHA256_PATTERN.fullmatch(receipt_sha256) is not None,
        "receipt publication seal SHA-256 is malformed",
    )
    return _canonical_json_bytes(
        {
            "schema_version": 1,
            "record_type": "q011_receipt_inode_bound_publication_success_seal",
            "publication_root_identity": _directory_identity(publication_descriptor),
            "receipt_name": receipt_name,
            "receipt_sha256": receipt_sha256,
            "receipt_identity": {
                "device": receipt_identity[0],
                "inode": receipt_identity[1],
            },
        }
    )


def _require_publication_seal_absent_at(
    parent_descriptor: int, receipt_name: str, label: str
) -> None:
    _require_absent_at(
        parent_descriptor,
        _publication_seal_name(receipt_name),
        f"{label} durable success seal",
    )


def _publish_publication_seal_at(
    acceptance_descriptor: int,
    publication_descriptor: int,
    receipt_name: str,
    receipt_payload: bytes,
    receipt_identity: tuple[int, int],
) -> tuple[int, int]:
    _require_same_file_at(
        publication_descriptor, receipt_name, receipt_identity, "canonical published receipt"
    )
    seal_name = _publication_seal_name(receipt_name)
    staging_name = f".{seal_name}.staging-{uuid.uuid4()}"
    _require_absent_at(acceptance_descriptor, seal_name, "receipt durable success seal")
    seal_payload = _publication_seal_payload(
        publication_descriptor,
        receipt_name,
        _sha256(receipt_payload),
        receipt_identity,
    )
    _write_exclusive_at(
        acceptance_descriptor,
        staging_name,
        seal_payload,
    )
    staging_identity = _file_identity_at(
        acceptance_descriptor, staging_name, "receipt staged durable success seal"
    )
    _fsync_descriptor(acceptance_descriptor)
    _require_same_file_at(
        publication_descriptor,
        receipt_name,
        receipt_identity,
        "canonical published receipt",
    )
    try:
        _rename_no_replace_at(acceptance_descriptor, staging_name, seal_name)
    except BaseException:
        # A rename wrapper can raise after the kernel committed the marker.
        # Reconcile the exact canonical seal so callers never report failure
        # after making this receipt externally consumable.
        _require_absent_at(
            acceptance_descriptor,
            staging_name,
            "receipt staged durable success seal",
        )
        _require_same_file_at(
            acceptance_descriptor,
            seal_name,
            staging_identity,
            "receipt durable success seal",
        )
        _require_same_file_at(
            publication_descriptor,
            receipt_name,
            receipt_identity,
            "canonical published receipt",
        )
        canonical_payload = _read_stable_readonly_regular_at(
            acceptance_descriptor,
            seal_name,
            "receipt durable success seal",
            max_bytes=MAX_JSON_BYTES,
        )
        _require(
            canonical_payload == seal_payload,
            "receipt durable success seal drifted during commit reconciliation",
        )
        _fsync_descriptor(acceptance_descriptor)
        return staging_identity
    _require_absent_at(
        acceptance_descriptor,
        staging_name,
        "receipt staged durable success seal",
    )
    _require_same_file_at(
        acceptance_descriptor,
        seal_name,
        staging_identity,
        "receipt durable success seal",
    )
    _require(
        _read_stable_readonly_regular_at(
            acceptance_descriptor,
            seal_name,
            "receipt durable success seal",
            max_bytes=MAX_JSON_BYTES,
        )
        == seal_payload,
        "receipt durable success seal drifted after commit",
    )
    _fsync_descriptor(acceptance_descriptor)
    return staging_identity


def _require_publication_seal_at(
    acceptance_descriptor: int,
    publication_descriptor: int,
    receipt_name: str,
    receipt_payload: bytes,
    receipt_identity: tuple[int, int],
    label: str,
) -> None:
    _require_same_file_at(publication_descriptor, receipt_name, receipt_identity, label)
    seal_name = _publication_seal_name(receipt_name)
    try:
        payload = _read_stable_readonly_regular_at(
            acceptance_descriptor,
            seal_name,
            f"{label} durable success seal",
            max_bytes=MAX_JSON_BYTES,
        )
    except OSError as error:
        raise PressurePilotPublicationError(
            f"{label} durable success seal is unavailable"
        ) from error
    _require(
        payload
        == _publication_seal_payload(
            publication_descriptor,
            receipt_name,
            _sha256(receipt_payload),
            receipt_identity,
        ),
        f"{label} durable success seal drifted",
    )
    _require_same_file_at(publication_descriptor, receipt_name, receipt_identity, label)


def _write_exclusive_below(root_descriptor: int, relative: str, payload: bytes) -> None:
    path = PurePosixPath(relative)
    _require(
        not path.is_absolute()
        and path.parts
        and all(part not in {"", ".", ".."} for part in path.parts),
        "pressure-pilot bundle member path is not canonical",
    )
    descriptor = os.dup(root_descriptor)
    try:
        for part in path.parts[:-1]:
            try:
                os.mkdir(part, mode=0o700, dir_fd=descriptor)
            except FileExistsError:
                pass
            child = os.open(part, _DIRECTORY_FLAGS, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        _write_exclusive_at(descriptor, path.parts[-1], payload)
    finally:
        os.close(descriptor)


def _freeze_anchored_tree(root_descriptor: int) -> None:
    """Remove write bits recursively without reopening the staging pathname."""

    entry_count = [0]

    def freeze(directory_descriptor: int, depth: int = 0) -> None:
        _require(
            depth <= MAX_TREE_DEPTH,
            f"pressure-pilot tree exceeds the {MAX_TREE_DEPTH}-level depth limit",
        )
        for name in _bounded_directory_names(
            directory_descriptor,
            "pressure-pilot staging directory",
        ):
            _require(
                entry_count[0] < MAX_TREE_ENTRIES,
                f"pressure-pilot tree exceeds the {MAX_TREE_ENTRIES}-entry limit",
            )
            entry_count[0] += 1
            observed = os.stat(name, dir_fd=directory_descriptor, follow_symlinks=False)
            if stat.S_ISDIR(observed.st_mode):
                child = os.open(name, _DIRECTORY_FLAGS, dir_fd=directory_descriptor)
                try:
                    opened = os.fstat(child)
                    _require(
                        (observed.st_dev, observed.st_ino)
                        == (opened.st_dev, opened.st_ino),
                        "pressure-pilot tree changed during anchored freeze",
                    )
                    freeze(child, depth + 1)
                    os.fchmod(child, opened.st_mode & ~_WRITE_BITS)
                    os.fsync(child)
                finally:
                    os.close(child)
                continue
            _require(
                stat.S_ISREG(observed.st_mode),
                "pressure-pilot tree has unsupported entry during freeze",
            )
            descriptor = os.open(name, _FILE_FLAGS, dir_fd=directory_descriptor)
            try:
                opened = os.fstat(descriptor)
                _require(
                    stat.S_ISREG(opened.st_mode)
                    and (observed.st_dev, observed.st_ino)
                    == (opened.st_dev, opened.st_ino),
                    "pressure-pilot tree changed during anchored freeze",
                )
                os.fchmod(descriptor, opened.st_mode & ~_WRITE_BITS)
                os.fsync(descriptor)
            finally:
                os.close(descriptor)

    freeze(root_descriptor)
    os.fchmod(root_descriptor, os.fstat(root_descriptor).st_mode & ~_WRITE_BITS)
    os.fsync(root_descriptor)


def _require_same_account_isolated_parent(parent_descriptor: int) -> None:
    metadata = os.fstat(parent_descriptor)
    _require(
        stat.S_ISDIR(metadata.st_mode)
        and metadata.st_uid == os.geteuid()
        and not metadata.st_mode & (stat.S_IWGRP | stat.S_IWOTH),
        "pressure-pilot publication requires a same-account isolated parent",
    )


def _lock_publication_transaction(
    anchor_descriptor: int, acceptance_descriptor: int
) -> None:
    try:
        fcntl.flock(anchor_descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as error:
        raise PressurePilotPublicationError(
            "pressure-pilot stable publication transaction lock is unavailable"
        ) from error
    try:
        fcntl.flock(acceptance_descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as error:
        fcntl.flock(anchor_descriptor, fcntl.LOCK_UN)
        raise PressurePilotPublicationError(
            "pressure-pilot acceptance-root publication transaction lock is unavailable"
        ) from error


def _close_descriptors(descriptors: Sequence[int | None]) -> OSError | None:
    first_error: OSError | None = None
    for descriptor in descriptors:
        if descriptor is None:
            continue
        try:
            os.close(descriptor)
        except OSError as error:
            if first_error is None:
                first_error = error
    return first_error


def _require_linked_commit_complete_at(
    parent_descriptor: int,
    source_name: str,
    destination_name: str,
    identity: tuple[int, int],
) -> None:
    _require_absent_at(
        parent_descriptor, source_name, "pressure-pilot linked fallback source"
    )
    _require_same_file_at(
        parent_descriptor,
        destination_name,
        identity,
        "pressure-pilot linked fallback destination",
    )
    destination = os.stat(
        destination_name, dir_fd=parent_descriptor, follow_symlinks=False
    )
    _require(
        destination.st_nlink == 1,
        "pressure-pilot linked fallback destination has an invalid link count",
    )


def _link_no_replace_file_at(
    parent_descriptor: int,
    source_name: str,
    destination_name: str,
    identity: tuple[int, int],
) -> None:
    source = os.stat(source_name, dir_fd=parent_descriptor, follow_symlinks=False)
    _require(
        stat.S_ISREG(source.st_mode)
        and (source.st_dev, source.st_ino) == identity
        and source.st_nlink == 1,
        "pressure-pilot linked fallback source is not singly linked",
    )
    try:
        os.link(
            source_name,
            destination_name,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
    except OSError as error:
        if error.errno == errno.EEXIST:
            raise PressurePilotPublicationError(
                f"pressure-pilot output collided during rename: {destination_name}"
            ) from error
        try:
            _require_same_file_at(
                parent_descriptor,
                source_name,
                identity,
                "pressure-pilot linked fallback source",
            )
            _require_same_file_at(
                parent_descriptor,
                destination_name,
                identity,
                "pressure-pilot linked fallback destination",
            )
        except (OSError, PressurePilotPublicationError):
            raise error
    except BaseException as error:
        try:
            _require_same_file_at(
                parent_descriptor,
                source_name,
                identity,
                "pressure-pilot linked fallback source",
            )
            _require_same_file_at(
                parent_descriptor,
                destination_name,
                identity,
                "pressure-pilot linked fallback destination",
            )
        except (OSError, PressurePilotPublicationError):
            raise error

    try:
        _require_same_file_at(
            parent_descriptor,
            source_name,
            identity,
            "pressure-pilot linked fallback source",
        )
        _require_same_file_at(
            parent_descriptor,
            destination_name,
            identity,
            "pressure-pilot linked fallback destination",
        )
        linked = os.stat(
            destination_name, dir_fd=parent_descriptor, follow_symlinks=False
        )
        _require(
            linked.st_nlink == 2,
            "pressure-pilot linked fallback commit has an invalid link count",
        )
    except (OSError, PressurePilotPublicationError):
        raise

    try:
        os.unlink(source_name, dir_fd=parent_descriptor)
        _require_linked_commit_complete_at(
            parent_descriptor, source_name, destination_name, identity
        )
    except BaseException as error:
        try:
            _require_linked_commit_complete_at(
                parent_descriptor, source_name, destination_name, identity
            )
        except BaseException:
            raise error


def _rename_no_replace_at(parent_descriptor: int, source_name: str, destination_name: str) -> None:
    _require(
        all(
            "/" not in name and name not in {"", ".", ".."}
            for name in (source_name, destination_name)
        ),
        "pressure-pilot descriptor-relative rename received an invalid name",
    )
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    error_number = errno.ENOSYS
    if renameat2 is not None:
        renameat2.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renameat2.restype = ctypes.c_int
        if renameat2(
            parent_descriptor,
            os.fsencode(source_name),
            parent_descriptor,
            os.fsencode(destination_name),
            _RENAME_NOREPLACE,
        ) == 0:
            return
        error_number = ctypes.get_errno()
    if error_number == errno.EEXIST:
        raise PressurePilotPublicationError(
            f"pressure-pilot output collided during rename: {destination_name}"
        )
    unsupported = {
        errno.EINVAL,
        errno.ENOSYS,
        getattr(errno, "ENOTSUP", errno.EINVAL),
        getattr(errno, "EOPNOTSUPP", errno.EINVAL),
    }
    if error_number not in unsupported:
        raise OSError(error_number, os.strerror(error_number), destination_name)
    _require_same_account_isolated_parent(parent_descriptor)
    source = os.stat(source_name, dir_fd=parent_descriptor, follow_symlinks=False)
    identity = source.st_dev, source.st_ino
    if stat.S_ISREG(source.st_mode):
        _link_no_replace_file_at(
            parent_descriptor, source_name, destination_name, identity
        )
        return
    raise PressurePilotPublicationError(
        "pressure-pilot directory publication requires exclusive final-name creation"
    )


def _bounded_directory_names(directory_fd: int, label: str) -> list[str]:
    names: list[str] = []
    with os.scandir(directory_fd) as entries:
        for entry in entries:
            _require(
                len(names) < MAX_DIRECTORY_ENTRIES,
                f"{label} exceeds the {MAX_DIRECTORY_ENTRIES}-entry directory limit",
            )
            names.append(entry.name)
    return sorted(names)


def _read_bounded_descriptor(descriptor: int, label: str, max_bytes: int) -> bytes:
    _require(
        os.fstat(descriptor).st_size <= max_bytes,
        f"{label} exceeds the {max_bytes}-byte size limit",
    )
    payload = bytearray()
    while True:
        chunk = os.read(
            descriptor,
            min(_READ_CHUNK_BYTES, max_bytes + 1 - len(payload)),
        )
        if not chunk:
            return bytes(payload)
        payload.extend(chunk)
        _require(
            len(payload) <= max_bytes,
            f"{label} exceeds the {max_bytes}-byte size limit",
        )


def _read_stable_readonly_regular(
    path: Path,
    label: str,
    *,
    max_bytes: int = MAX_RETAINED_FILE_BYTES,
) -> bytes:
    descriptor = os.open(path, _FILE_FLAGS)
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode) and not before.st_mode & 0o222,
            f"{label} is not a read-only regular file",
        )
        payload = _read_bounded_descriptor(descriptor, label, max_bytes)
        after = os.fstat(descriptor)
        current = os.stat(path, follow_symlinks=False)
        _require(
            (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
            == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and len(payload) == after.st_size,
            f"{label} changed while reading",
        )
        return payload
    finally:
        os.close(descriptor)


def _read_stable_readonly_regular_at(
    parent_descriptor: int,
    name: str,
    label: str,
    *,
    max_bytes: int = MAX_RETAINED_FILE_BYTES,
) -> bytes:
    _require("/" not in name, f"{label}: descriptor-relative read received a nested path")
    descriptor = os.open(name, _FILE_FLAGS, dir_fd=parent_descriptor)
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and not before.st_mode & 0o222
            and before.st_nlink == 1,
            f"{label} is not a singly linked read-only regular file",
        )
        payload = _read_bounded_descriptor(descriptor, label, max_bytes)
        after = os.fstat(descriptor)
        current = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
        _require(
            (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
            == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and after.st_nlink == current.st_nlink == 1
            and len(payload) == after.st_size,
            f"{label} changed while reading",
        )
        return payload
    finally:
        os.close(descriptor)


def _verify_published_case_descriptor_bounded(
    tree: object,
    case_id: str,
    expected_descriptor_sha256: str,
) -> dict[str, object]:
    """Recompute a raw-case descriptor without the producer's unbounded reader."""
    try:
        _require(
            _SHA256_PATTERN.fullmatch(expected_descriptor_sha256) is not None,
            "expected raw-case descriptor SHA-256 is malformed",
        )
        case_verifier.load_inventory(tree)
        payload = _read_stable_readonly_regular_at(
            tree.analysis_fd,
            "analysis.json",
            "published raw-case descriptor",
            max_bytes=MAX_JSON_BYTES,
        )
        _require(
            _sha256(payload) == expected_descriptor_sha256,
            "raw-case descriptor SHA-256 drifted",
        )
        decoded = case_verifier._decode_json(payload, "raw-case descriptor")
        _require(type(decoded) is dict, "raw-case descriptor must be an object")
        _require(
            payload == case_verifier.canonical_json_bytes(decoded),
            "raw-case descriptor is not canonical JSON",
        )
        recomputed = case_verifier._analyze_tree(tree, case_id)
        case_verifier._strict_equal(decoded, recomputed, "raw-case descriptor")
        tree.require_tree_closure()
        tree.require_analysis_identity()
        _require(
            _read_stable_readonly_regular_at(
                tree.analysis_fd,
                "analysis.json",
                "published raw-case descriptor",
                max_bytes=MAX_JSON_BYTES,
            )
            == payload,
            "raw-case descriptor changed after verification",
        )
        return recomputed
    except RecursionError as error:
        raise PressurePilotPublicationError(
            "raw-case descriptor exceeds the supported nesting depth"
        ) from error


def _open_absolute_directory(path: Path) -> int:
    absolute = Path(os.path.abspath(path))
    descriptor = os.open("/", _DIRECTORY_FLAGS)
    try:
        for part in absolute.parts[1:]:
            child = os.open(part, _DIRECTORY_FLAGS, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


class ImmutablePressurePilotBundle:
    """Retain one no-follow root descriptor while verifying aggregate closure."""

    def __init__(
        self,
        root: str | Path,
        *,
        inherited_root_fd: int | None = None,
        parent_fd: int | None = None,
        entry_name: str | None = None,
    ) -> None:
        self.root = Path(os.path.abspath(root))
        self._inherited_root_fd = inherited_root_fd
        self._parent_fd = parent_fd
        self._entry_name = entry_name
        self._root_fd: int | None = None
        self._directory_identities: dict[str, tuple[int, int]] | None = None
        self._file_identities: dict[str, tuple[int, int]] | None = None

    @property
    def root_fd(self) -> int:
        if self._root_fd is None:
            raise PressurePilotPublicationError("pressure-pilot bundle verifier is not open")
        return self._root_fd

    def __enter__(self) -> "ImmutablePressurePilotBundle":
        self._root_fd = (
            _open_absolute_directory(self.root)
            if self._inherited_root_fd is None
            else os.dup(self._inherited_root_fd)
        )
        try:
            self.require_path_identity()
            _require(
                stat.S_ISDIR(os.fstat(self.root_fd).st_mode)
                and not os.fstat(self.root_fd).st_mode & 0o222,
                "pressure-pilot bundle root is not read-only",
            )
            return self
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def __exit__(self, *_: object) -> None:
        if self._root_fd is not None:
            os.close(self._root_fd)
            self._root_fd = None
        self._directory_identities = None
        self._file_identities = None

    def require_path_identity(self) -> None:
        if self._parent_fd is not None and self._entry_name is not None:
            _require_same_directory_at(
                self._parent_fd,
                self._entry_name,
                self.root_fd,
                "pressure-pilot bundle root",
            )
            return
        actual_fd = _open_absolute_directory(self.root)
        try:
            expected = os.fstat(self.root_fd)
            actual = os.fstat(actual_fd)
            _require(
                (expected.st_dev, expected.st_ino) == (actual.st_dev, actual.st_ino),
                "pressure-pilot bundle root path changed during verification",
            )
        finally:
            os.close(actual_fd)

    def file_paths(self) -> set[str]:
        _require(
            self._file_identities is not None,
            "pressure-pilot bundle verifier has not established file closure",
        )
        return set(self._file_identities)

    def _scan(
        self,
        directory_fd: int,
        prefix: tuple[str, ...] = (),
        *,
        _entry_count: list[int] | None = None,
        _depth: int = 0,
    ) -> tuple[list[str], dict[str, tuple[int, int]], dict[str, tuple[int, int]]]:
        _require(
            _depth <= MAX_TREE_DEPTH,
            f"pressure-pilot tree exceeds the {MAX_TREE_DEPTH}-level depth limit",
        )
        if _entry_count is None:
            _entry_count = [0]
        paths = []
        directories = {}
        files = {}
        label = "/".join(prefix) or "."
        names = _bounded_directory_names(
            directory_fd,
            f"pressure-pilot directory {label}",
        )
        for name in names:
            _require(
                _entry_count[0] < MAX_TREE_ENTRIES,
                f"pressure-pilot tree exceeds the {MAX_TREE_ENTRIES}-entry limit",
            )
            _entry_count[0] += 1
            relative = "/".join((*prefix, name))
            metadata = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            if stat.S_ISREG(metadata.st_mode):
                _require(not metadata.st_mode & 0o222, f"pressure-pilot member is writable: {relative}")
                descriptor = os.open(name, _FILE_FLAGS, dir_fd=directory_fd)
                try:
                    opened = os.fstat(descriptor)
                    _require(
                        stat.S_ISREG(opened.st_mode)
                        and not opened.st_mode & 0o222
                        and (metadata.st_dev, metadata.st_ino)
                        == (opened.st_dev, opened.st_ino),
                        f"pressure-pilot member changed during verification: {relative}",
                    )
                    _require(
                        opened.st_size <= MAX_RETAINED_FILE_BYTES,
                        f"pressure-pilot member exceeds the {MAX_RETAINED_FILE_BYTES}-byte "
                        f"size limit: {relative}",
                    )
                    files[relative] = (opened.st_dev, opened.st_ino)
                finally:
                    os.close(descriptor)
                paths.append(relative)
            elif stat.S_ISDIR(metadata.st_mode):
                _require(not metadata.st_mode & 0o222, f"pressure-pilot directory is writable: {relative}")
                child_fd = os.open(name, _DIRECTORY_FLAGS, dir_fd=directory_fd)
                try:
                    opened = os.fstat(child_fd)
                    _require(
                        stat.S_ISDIR(opened.st_mode)
                        and not opened.st_mode & 0o222
                        and (metadata.st_dev, metadata.st_ino)
                        == (opened.st_dev, opened.st_ino),
                        f"pressure-pilot directory changed during verification: {relative}",
                    )
                    directories[relative] = (opened.st_dev, opened.st_ino)
                    child_paths, child_directories, child_files = self._scan(
                        child_fd,
                        (*prefix, name),
                        _entry_count=_entry_count,
                        _depth=_depth + 1,
                    )
                    _require(bool(child_paths), f"pressure-pilot tree contains empty directory: {relative}")
                    paths.extend(child_paths)
                    directories.update(child_directories)
                    files.update(child_files)
                    after = os.fstat(child_fd)
                    entry = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
                    _require(
                        (entry.st_dev, entry.st_ino) == (after.st_dev, after.st_ino),
                        f"pressure-pilot directory changed during verification: {relative}",
                    )
                finally:
                    os.close(child_fd)
            else:
                raise PressurePilotPublicationError(f"pressure-pilot tree has unsupported entry: {relative}")
        return paths, directories, files

    def read(
        self,
        relative: str,
        *,
        max_bytes: int = MAX_RETAINED_FILE_BYTES,
    ) -> bytes:
        if relative == pilot.MANIFEST_NAME or relative.endswith(".rst.manifest"):
            max_bytes = min(max_bytes, MAX_JSON_BYTES)
        self.require_path_identity()
        payload = _ARTIFACT_HELPERS._read_at(
            self.root_fd,
            relative,
            self._directory_identities,
            self._file_identities,
            max_bytes=max_bytes,
        )
        self.require_path_identity()
        return payload

    def verify(self, expected_manifest_sha256: str) -> dict[str, object]:
        _require(
            _SHA256_PATTERN.fullmatch(expected_manifest_sha256) is not None,
            "expected pressure-pilot manifest SHA-256 is malformed",
        )
        manifest_payload = self.read(pilot.MANIFEST_NAME, max_bytes=MAX_JSON_BYTES)
        _require(_sha256(manifest_payload) == expected_manifest_sha256, "pressure-pilot manifest SHA-256 drifted")
        try:
            manifest = pilot._manifest_schema(manifest_payload)
        except RecursionError as error:
            raise PressurePilotPublicationError(
                "pressure-pilot manifest exceeds the supported JSON nesting depth"
            ) from error
        pilot._validate_case_identity(manifest["cases"])
        expected = {pilot.MANIFEST_NAME: (len(manifest_payload), _sha256(manifest_payload))}
        for case in manifest["cases"]:
            expected[case["stdout"]["path"]] = (
                None,
                case["stdout"]["sha256"],
            )
            for snapshot in case["snapshots"]:
                for name in ("mhd_w_bcc", "bmag", "prtcl_jx", "j2", "prtcl_all"):
                    binding = snapshot[name]
                    expected[binding["path"]] = (None, binding["sha256"])
            restart = case["terminal_restart"]
            for binding in (restart["manifest"], restart["manifest_complete"]):
                expected[binding["path"]] = (None, binding["sha256"])
            for member in restart["members"]:
                for binding in (member["artifact"], member["complete"]):
                    expected[binding["path"]] = (None, binding["sha256"])
        _require(set(expected) == pilot._declared_paths(manifest), "pressure-pilot manifest member closure drifted")
        paths, directories, files = self._scan(self.root_fd)
        _require(paths == sorted(expected), "pressure-pilot tree closure drifted")
        if self._directory_identities is None:
            self._directory_identities = directories
            self._file_identities = files
        else:
            _require(directories == self._directory_identities, "pressure-pilot directory identity changed")
            _require(files == self._file_identities, "pressure-pilot member identity changed")
        for relative, (expected_size, expected_digest) in expected.items():
            payload = self.read(relative)
            _require(
                (expected_size is None or len(payload) == expected_size)
                and _sha256(payload) == expected_digest,
                f"pressure-pilot member checksum drifted: {relative}",
            )
        self.require_path_identity()
        return manifest


def _verify_pressure_pilot_bundle_at(
    root: Path,
    parent_descriptor: int,
    entry_name: str,
    root_descriptor: int,
    expected_manifest_sha256: str,
    *,
    authorized_publication_root: Path,
) -> dict[str, Any]:
    """Verify one aggregate tree through the caller's retained root descriptor."""
    try:
        with ImmutablePressurePilotBundle(
            root,
            inherited_root_fd=root_descriptor,
            parent_fd=parent_descriptor,
            entry_name=entry_name,
        ) as bundle:
            bundle.verify(expected_manifest_sha256)
            result = pilot.analyze_pressure_pilot_bundle(
                root,
                expected_manifest_sha256,
                authorized_publication_root=authorized_publication_root,
                member_reader=bundle.read,
                actual_files=bundle.file_paths(),
            )
            bundle.verify(expected_manifest_sha256)
            return result
    except RecursionError as error:
        raise PressurePilotPublicationError(
            "pressure-pilot aggregate evidence exceeds the supported nesting depth"
        ) from error


def verify_published_pressure_pilot_bundle(
    output_path: str | Path,
    expected_manifest_sha256: str,
    *,
    authorized_publication_root: Path = AUTHORIZED_PUBLICATION_ROOT,
) -> dict[str, Any]:
    """Verify a nonproduction fixture bundle without accepting it for consumption."""
    publication_root = _canonical_existing_directory(
        authorized_publication_root, "authorized PIC publication root"
    )
    _require(
        publication_root != Path(os.path.abspath(AUTHORIZED_PUBLICATION_ROOT)),
        "production pressure-pilot bundle verification requires receipt consumption",
    )
    root = _direct_publication_target(
        output_path, publication_root, "pressure-pilot retained bundle"
    )
    parent_descriptor = _open_absolute_directory(publication_root)
    root_descriptor: int | None = None
    try:
        _require_same_directory(
            publication_root, parent_descriptor, "authorized PIC publication root"
        )
        root_descriptor = os.open(root.name, _DIRECTORY_FLAGS, dir_fd=parent_descriptor)
        _require_same_directory_at(
            parent_descriptor, root.name, root_descriptor, "pressure-pilot retained bundle"
        )
        return _verify_pressure_pilot_bundle_at(
            root,
            parent_descriptor,
            root.name,
            root_descriptor,
            expected_manifest_sha256,
            authorized_publication_root=publication_root,
        )
    finally:
        if root_descriptor is not None:
            os.close(root_descriptor)
        os.close(parent_descriptor)


def _static_source_bindings() -> dict[str, object]:
    source_payloads: dict[str, bytes] = {}
    successor = pilot._load_postrun_source_authorization_successor(
        source_payloads=source_payloads
    )
    successor_payload = source_payloads["postrun_aggregate_source_authorization"]
    registered_execution_payload = pilot._regular_bytes(
        _REGISTERED_EXECUTION_PREREGISTRATION_PATH,
        "pressure-pilot registered-execution preregistration",
    )
    _require(
        _sha256(registered_execution_payload)
        == pilot.REGISTERED_EXECUTION_PREREGISTRATION_SHA256,
        "pressure-pilot registered-execution preregistration SHA-256 drifted",
    )
    return {
        "postrun_aggregate_source_authorization": {
            "path": pilot.PREREGISTRATION_PATH.relative_to(pilot.REPO_ROOT).as_posix(),
            "sha256": _sha256(successor_payload),
        },
        "registered_execution_preregistration": {
            "path": _REGISTERED_EXECUTION_PREREGISTRATION_PATH.relative_to(
                pilot.REPO_ROOT
            ).as_posix(),
            "sha256": _sha256(registered_execution_payload),
        },
        "historical_v2_execution_preregistration": successor[
            "historical_v2_execution_preregistration"
        ],
        "reviewed_source_closure": successor["source_closure"],
    }


def _source_closure_sha256(source_bindings: Mapping[str, object]) -> str:
    return _sha256(
        _canonical_json_bytes(
            {
                "postrun_aggregate_source_authorization": source_bindings[
                    "postrun_aggregate_source_authorization"
                ],
                "reviewed_source_closure": source_bindings["reviewed_source_closure"],
            }
        )
    )


def _archive_member_payload(
    archive: tarfile.TarFile, relative: str, label: str
) -> bytes:
    members = [member for member in archive.getmembers() if member.name == relative]
    _require(
        len(members) == 1 and members[0].isfile(),
        f"{label} is absent or not one regular archive member",
    )
    _require(
        members[0].size <= MAX_RETAINED_FILE_BYTES,
        f"{label} exceeds the {MAX_RETAINED_FILE_BYTES}-byte size limit",
    )
    stream = archive.extractfile(members[0])
    _require(stream is not None, f"{label} cannot be read from source archive")
    payload = stream.read(MAX_RETAINED_FILE_BYTES + 1)
    _require(
        len(payload) == members[0].size and len(payload) <= MAX_RETAINED_FILE_BYTES,
        f"{label} changed while reading or exceeds its size limit",
    )
    return payload


def _trusted_source_archive(commitish: str = "HEAD") -> tuple[str, bytes]:
    try:
        repository = TRUSTED_SOURCE_REPOSITORY.resolve(strict=True)
    except OSError as error:
        raise PressurePilotPublicationError(
            "trusted pressure-pilot source repository is unavailable"
        ) from error
    _require(repository.is_dir(), "trusted pressure-pilot source repository is invalid")
    _require(GIT_EXECUTABLE.is_file(), "trusted Git executable is unavailable")
    commit_result = subprocess.run(
        [
            str(GIT_EXECUTABLE),
            "-C",
            str(repository),
            "rev-parse",
            "--verify",
            f"{commitish}^{{commit}}",
        ],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    _require(
        commit_result.returncode == 0,
        "pressure-pilot trusted source commit is unavailable",
    )
    commit = commit_result.stdout.decode("ascii").strip()
    _require(
        _GIT_COMMIT_PATTERN.fullmatch(commit) is not None,
        "pressure-pilot trusted source commit is malformed",
    )
    archive_result = subprocess.run(
        [
            str(GIT_EXECUTABLE),
            "-C",
            str(repository),
            "archive",
            "--format=tar",
            commit,
        ],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    _require(
        archive_result.returncode == 0,
        "pressure-pilot trusted source archive cannot be reconstructed",
    )
    return commit, archive_result.stdout


def _verify_source_archive_payload(
    source_bindings: Mapping[str, object],
    archive_payload: bytes,
    *,
    expected_commit: str,
) -> None:
    try:
        with tarfile.open(fileobj=io.BytesIO(archive_payload), mode="r:") as archive:
            commit = archive.pax_headers.get("comment")
            _require(
                commit == expected_commit,
                "pressure-pilot worker source archive Git commit drifted",
            )
            required = {
                record["path"]: record["sha256"]
                for record in source_bindings["reviewed_source_closure"]
            }
            authorization = source_bindings["postrun_aggregate_source_authorization"]
            required[authorization["path"]] = authorization["sha256"]
            for relative, digest in required.items():
                payload = _archive_member_payload(
                    archive, relative, f"reviewed source archive member {relative}"
                )
                _require(
                    _sha256(payload) == digest,
                    f"reviewed source archive member SHA-256 drifted: {relative}",
                )
    except tarfile.TarError as error:
        raise PressurePilotPublicationError(
            "pressure-pilot worker source archive is not a readable tar file"
        ) from error


def _runtime_source_archive_binding(
    source_bindings: Mapping[str, object],
    *,
    require_verified_archive: bool,
) -> dict[str, object]:
    archive_path_text = os.environ.get("PIC_PRESSURE_PUBLICATION_SOURCE_ARCHIVE_PATH")
    if archive_path_text is None:
        _require(
            not require_verified_archive,
            "pressure-pilot production publication requires a verified worker source archive",
        )
        return {
            "execution_mode": "direct_api_nonproduction_only",
            "git_commit": None,
            "archive_sha256": None,
            "verified_source_closure_sha256": None,
        }
    archive_path = Path(archive_path_text)
    _require(
        archive_path.is_absolute(),
        "pressure-pilot worker source archive path must be absolute",
    )
    archive_payload = _read_stable_readonly_regular(
        archive_path, "pressure-pilot worker source archive"
    )
    snapshot_root_text = os.environ.get(WORKER_SOURCE_SNAPSHOT_ROOT_ENV)
    _require(
        snapshot_root_text is not None,
        "pressure-pilot production publication requires an extracted worker source snapshot",
    )
    snapshot_root = Path(snapshot_root_text)
    _require(
        snapshot_root.is_absolute(),
        "pressure-pilot worker source snapshot path must be absolute",
    )
    try:
        snapshot_root = snapshot_root.resolve(strict=True)
        executing_source_root = EXECUTING_SOURCE_ROOT.resolve(strict=True)
        trusted_repository = TRUSTED_SOURCE_REPOSITORY.resolve(strict=True)
    except OSError as error:
        raise PressurePilotPublicationError(
            "pressure-pilot worker source snapshot is unavailable"
        ) from error
    _require(
        snapshot_root.is_dir() and snapshot_root == executing_source_root,
        "pressure-pilot production API is not executing from its worker source snapshot",
    )
    _require(
        snapshot_root != trusted_repository,
        "pressure-pilot production API must not execute from the trusted checkout",
    )
    commit, trusted_archive_payload = _trusted_source_archive()
    _require(
        archive_payload == trusted_archive_payload,
        "pressure-pilot worker source archive differs from trusted repository HEAD archive",
    )
    _verify_source_archive_payload(
        source_bindings, archive_payload, expected_commit=commit
    )
    try:
        with tarfile.open(fileobj=io.BytesIO(archive_payload), mode="r:") as archive:
            required = {
                record["path"] for record in source_bindings["reviewed_source_closure"]
            }
            required.add(source_bindings["postrun_aggregate_source_authorization"]["path"])
            for relative in sorted(required):
                archived = _archive_member_payload(
                    archive, relative, f"worker snapshot archive member {relative}"
                )
                snapshot = _read_stable_readonly_regular(
                    snapshot_root / relative,
                    f"worker snapshot source member {relative}",
                )
                _require(
                    snapshot == archived,
                    f"worker snapshot source member drifted: {relative}",
                )
    except tarfile.TarError as error:
        raise PressurePilotPublicationError(
            "pressure-pilot worker source archive is not a readable tar file"
        ) from error
    return {
        "execution_mode": "worker_extracted_git_archive_head_verified",
        "git_commit": commit,
        "archive_sha256": _sha256(archive_payload),
        "verified_source_closure_sha256": _source_closure_sha256(source_bindings),
    }


def _source_bindings(*, require_verified_archive: bool = False) -> dict[str, object]:
    static = _static_source_bindings()
    return {
        **static,
        "runtime_source_archive": _runtime_source_archive_binding(
            static, require_verified_archive=require_verified_archive
        ),
    }


def _validate_retained_source_bindings(
    value: object, *, require_verified_archive: bool = False
) -> None:
    _require(type(value) is dict, "pressure-pilot receipt source bindings are malformed")
    expected = _static_source_bindings()
    _require(
        set(value) == {*expected, "runtime_source_archive"},
        "pressure-pilot receipt source bindings drifted",
    )
    for key, binding in expected.items():
        _require(
            value[key] == binding,
            f"pressure-pilot receipt source binding drifted: {key}",
        )
    archive = value["runtime_source_archive"]
    _require(
        type(archive) is dict
        and set(archive)
        == {
            "execution_mode",
            "git_commit",
            "archive_sha256",
            "verified_source_closure_sha256",
        },
        "pressure-pilot runtime source archive binding is malformed",
    )
    if archive["execution_mode"] == "direct_api_nonproduction_only":
        _require(
            not require_verified_archive
            and archive["git_commit"] is None
            and archive["archive_sha256"] is None
            and archive["verified_source_closure_sha256"] is None,
            "pressure-pilot direct API archive binding drifted",
        )
        return
    _require(
        archive["execution_mode"] == "worker_extracted_git_archive_head_verified"
        and type(archive["git_commit"]) is str
        and _GIT_COMMIT_PATTERN.fullmatch(archive["git_commit"]) is not None
        and type(archive["archive_sha256"]) is str
        and _SHA256_PATTERN.fullmatch(archive["archive_sha256"]) is not None,
        "pressure-pilot runtime source archive binding is malformed",
    )
    _require(
        archive["verified_source_closure_sha256"] == _source_closure_sha256(expected),
        "pressure-pilot retained verified source closure SHA-256 drifted",
    )
    commit, archive_payload = _trusted_source_archive(archive["git_commit"])
    _require(
        commit == archive["git_commit"]
        and _sha256(archive_payload) == archive["archive_sha256"],
        "pressure-pilot retained trusted source archive binding drifted",
    )
    _verify_source_archive_payload(expected, archive_payload, expected_commit=commit)


def _manifest(case_descriptors: Mapping[str, Mapping[str, object]]) -> dict[str, object]:
    policy, _approved_snapshot_metadata, preregistration_payload = (
        pilot._load_policy_and_snapshot_metadata()
    )
    registered_execution_payload = pilot._regular_bytes(
        _REGISTERED_EXECUTION_PREREGISTRATION_PATH,
        "pressure-pilot registered-execution preregistration",
    )
    _require(
        _sha256(registered_execution_payload)
        == pilot.REGISTERED_EXECUTION_PREREGISTRATION_SHA256,
        "pressure-pilot registered-execution preregistration SHA-256 drifted",
    )
    return {
        "schema_version": 1,
        "record_type": "q011_section54_pressure_pilot_bundle_manifest",
        "evidence_class": pilot.EVIDENCE_CLASS,
        "qualification_effect": pilot.QUALIFICATION_EFFECT,
        "active_deck_binding": policy["active_deck_binding"],
        "preregistration_binding": {
            "path": pilot.PREREGISTRATION_PATH.relative_to(pilot.REPO_ROOT).as_posix(),
            "sha256": _sha256(preregistration_payload),
        },
        "registered_execution_preregistration_binding": {
            "path": _REGISTERED_EXECUTION_PREREGISTRATION_PATH.relative_to(
                pilot.REPO_ROOT
            ).as_posix(),
            "sha256": _sha256(registered_execution_payload),
        },
        "cases": [
            case_descriptors[case_id]["manifest_case"]
            for case_id in case_verifier.CASE_IDS
        ],
    }


def publish_pressure_pilot_bundle(
    output_path: str | Path,
    *,
    receipt_path: str | Path,
    analysis_result_path: str | Path,
    case_artifact_dirs: Mapping[str, str | Path],
    case_descriptor_sha256: Mapping[str, str],
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> dict[str, str]:
    """Copy four descriptor-verified raw cases into one immutable aggregate tree."""
    pic_root, parent = _publication_root(authorized_pic_root)
    acceptance_root = _publication_acceptance_root(pic_root)
    transaction_anchor = _publication_transaction_anchor(pic_root)
    target = _direct_publication_target(output_path, parent, "pressure-pilot output")
    receipt_target = _direct_publication_target(
        receipt_path, parent, "pressure-pilot receipt"
    )
    result_target = _direct_publication_target(
        analysis_result_path, parent, "pressure-pilot aggregate analysis"
    )
    _require(set(case_artifact_dirs) == set(case_verifier.CASE_IDS), "pressure-pilot source tree set drifted")
    _require(set(case_descriptor_sha256) == set(case_verifier.CASE_IDS), "pressure-pilot descriptor set drifted")
    case_artifact_dirs = {
        case_id: _authorized_raw_case_root(case_artifact_dirs[case_id], pic_root)
        for case_id in case_verifier.CASE_IDS
    }
    with ExitStack() as preflight_stack:
        preflight_descriptors = {}
        for case_id in case_verifier.CASE_IDS:
            tree = preflight_stack.enter_context(
                case_verifier.StructuredArtifactTree(Path(case_artifact_dirs[case_id]))
            )
            preflight_descriptors[case_id] = _verify_published_case_descriptor_bounded(
                tree, case_id, case_descriptor_sha256[case_id]
            )
            tree.require_tree_closure()
        preflight_manifest = _manifest(preflight_descriptors)
        pilot._validate_case_identity(preflight_manifest["cases"])
        preflight_expected_targets = pilot._declared_paths(preflight_manifest) - {
            pilot.MANIFEST_NAME
        }
        preflight_declared_targets = {
            str(member["path"])
            for descriptor in preflight_descriptors.values()
            for member in descriptor["bundle_members"]
        }
        _require(
            preflight_declared_targets == preflight_expected_targets,
            "verified raw-case bundle member closure drifted",
        )
    _source_bindings(require_verified_archive=_is_production_pic_root(pic_root))
    build_target = target
    receipt_staging = parent / f".{receipt_target.name}.staging-{uuid.uuid4()}"
    result_staging = parent / f".{result_target.name}.staging-{uuid.uuid4()}"
    publication_descriptor = _open_absolute_directory(parent)
    acceptance_descriptor = _open_absolute_directory(acceptance_root)
    transaction_descriptor = _open_absolute_directory(transaction_anchor)
    staging_descriptor: int | None = None
    renamed = False
    receipt_renamed = False
    result_renamed = False
    guard_armed = False
    seal_committed = False
    seal_identity: tuple[int, int] | None = None
    receipt_identity: tuple[int, int] | None = None
    guard_identity: tuple[int, int] | None = None
    try:
        _lock_publication_transaction(transaction_descriptor, acceptance_descriptor)
        _require_same_directory(
            transaction_anchor,
            transaction_descriptor,
            "stable publication transaction anchor",
        )
        _require_same_directory(
            parent, publication_descriptor, "authorized PIC publication root"
        )
        _require_same_directory(
            acceptance_root,
            acceptance_descriptor,
            "authorized PIC publication acceptance root",
        )
        _require_absent_at(publication_descriptor, target.name, "pressure-pilot output")
        _require_absent_at(
            publication_descriptor, receipt_target.name, "pressure-pilot receipt"
        )
        _require_absent_at(
            publication_descriptor, result_target.name, "pressure-pilot result"
        )
        _require_publication_guard_absent_at(
            publication_descriptor,
            receipt_target.name,
            "pressure-pilot receipt",
        )
        _require_publication_seal_absent_at(
            acceptance_descriptor,
            receipt_target.name,
            "pressure-pilot receipt",
        )
        guard_identity = _arm_publication_guard_at(
            publication_descriptor, receipt_target.name
        )
        guard_armed = True
        _require_same_account_isolated_parent(publication_descriptor)
        os.mkdir(build_target.name, mode=0o700, dir_fd=publication_descriptor)
        renamed = True
        staging_descriptor = os.open(
            build_target.name, _DIRECTORY_FLAGS, dir_fd=publication_descriptor
        )
        _require_same_directory_at(
            publication_descriptor,
            build_target.name,
            staging_descriptor,
            "pressure-pilot public build tree",
        )
        with ExitStack() as stack:
            trees = {}
            descriptors = {}
            inventories = {}
            for case_id in case_verifier.CASE_IDS:
                tree = stack.enter_context(
                    case_verifier.StructuredArtifactTree(Path(case_artifact_dirs[case_id]))
                )
                descriptor = _verify_published_case_descriptor_bounded(
                    tree, case_id, case_descriptor_sha256[case_id]
                )
                trees[case_id] = tree
                descriptors[case_id] = descriptor
                inventories[case_id] = case_verifier.load_inventory(tree)
            manifest = _manifest(descriptors)
            pilot._validate_case_identity(manifest["cases"])
            expected_targets = pilot._declared_paths(manifest) - {pilot.MANIFEST_NAME}
            declared_targets = {
                str(member["path"])
                for descriptor in descriptors.values()
                for member in descriptor["bundle_members"]
            }
            _require(declared_targets == expected_targets, "verified raw-case bundle member closure drifted")

            copied = set()
            for case_id in case_verifier.CASE_IDS:
                tree = trees[case_id]
                inventory = inventories[case_id]
                for member in descriptors[case_id]["bundle_members"]:
                    source_path = str(member["source_path"])
                    target_path = str(member["path"])
                    _require(target_path not in copied, f"duplicate pressure-pilot target member: {target_path}")
                    payload = case_verifier.read_inventory_bytes(tree, inventory, source_path)
                    _require(
                        len(payload) == member["size"] and _sha256(payload) == member["sha256"],
                        f"verified raw-case member changed before copy: {source_path}",
                    )
                    _write_exclusive_below(staging_descriptor, target_path, payload)
                    copied.add(target_path)
            _require(copied == expected_targets, "pressure-pilot copied member closure drifted")
            manifest_payload = _canonical_json_bytes(manifest)
            manifest_sha256 = _sha256(manifest_payload)
            _write_exclusive_at(staging_descriptor, pilot.MANIFEST_NAME, manifest_payload)
            _require_same_directory_at(
                publication_descriptor,
                build_target.name,
                staging_descriptor,
                "pressure-pilot public build tree",
            )
            _freeze_anchored_tree(staging_descriptor)
            _require_same_directory_at(
                publication_descriptor,
                build_target.name,
                staging_descriptor,
                "pressure-pilot public build tree",
            )
            result = _verify_pressure_pilot_bundle_at(
                build_target,
                publication_descriptor,
                build_target.name,
                staging_descriptor,
                manifest_sha256,
                authorized_publication_root=parent,
            )
            result_payload = _canonical_json_bytes(result)
            result_sha256 = _sha256(result_payload)
            _write_exclusive_at(publication_descriptor, result_staging.name, result_payload)
            receipt = {
                "schema_version": 1,
                "record_type": "q011_section54_pressure_pilot_bundle_publication_receipt",
                "evidence_class": pilot.EVIDENCE_CLASS,
                "qualification_effect": pilot.QUALIFICATION_EFFECT,
                "consumption_rule": CONSUMPTION_RULE,
                "publication_root_identity": _directory_identity(publication_descriptor),
                "aggregate_bundle": {
                    "path": str(target),
                    "manifest_sha256": manifest_sha256,
                },
                "aggregate_analysis": {
                    "path": str(result_target),
                    "sha256": result_sha256,
                },
                "source_bindings": _source_bindings(
                    require_verified_archive=_is_production_pic_root(pic_root)
                ),
                "raw_cases": [
                    {
                        "case_id": case_id,
                        "artifact_dir": str(Path(case_artifact_dirs[case_id]).resolve()),
                        "descriptor_path": case_verifier.CASE_DESCRIPTOR_PATH,
                        "descriptor_sha256": case_descriptor_sha256[case_id],
                        "artifact_inventory_sha256": descriptors[case_id][
                            "artifact_inventory_sha256"
                        ],
                        "runtime_artifacts": descriptors[case_id]["runtime_artifacts"],
                    }
                    for case_id in case_verifier.CASE_IDS
                ],
            }
            receipt_payload = _canonical_json_bytes(receipt)
            receipt_sha256 = _sha256(receipt_payload)
            _write_exclusive_at(publication_descriptor, receipt_staging.name, receipt_payload)
            receipt_identity = _file_identity_at(
                publication_descriptor, receipt_staging.name, "pressure-pilot staged receipt"
            )
            for tree in trees.values():
                tree.require_tree_closure()
            _require_same_directory(parent, publication_descriptor, "authorized PIC publication root")
            _require_same_directory_at(
                publication_descriptor,
                target.name,
                staging_descriptor,
                "pressure-pilot public build tree",
            )
            _fsync_descriptor(publication_descriptor)
            verified_result = _verify_pressure_pilot_bundle_at(
                target,
                publication_descriptor,
                target.name,
                staging_descriptor,
                manifest_sha256,
                authorized_publication_root=parent,
            )
            _require(
                _canonical_json_bytes(verified_result) == result_payload,
                "pressure-pilot aggregate analysis drifted after publication",
            )
            _require_same_directory(parent, publication_descriptor, "authorized PIC publication root")
            _require_absent_at(
                publication_descriptor, result_target.name, "pressure-pilot result"
            )
            _rename_no_replace_at(
                publication_descriptor, result_staging.name, result_target.name
            )
            result_renamed = True
            _fsync_descriptor(publication_descriptor)
            _require_same_directory(parent, publication_descriptor, "authorized PIC publication root")
            _require_absent_at(
                publication_descriptor, receipt_target.name, "pressure-pilot receipt"
            )
            _rename_no_replace_at(
                publication_descriptor, receipt_staging.name, receipt_target.name
            )
            receipt_renamed = True
            _fsync_descriptor(publication_descriptor)
        verified_result = _verify_pressure_pilot_bundle_at(
            target,
            publication_descriptor,
            target.name,
            staging_descriptor,
            manifest_sha256,
            authorized_publication_root=parent,
        )
        _require(
            _canonical_json_bytes(verified_result) == result_payload,
            "pressure-pilot aggregate analysis drifted after receipt publication",
        )
        _require(
            _read_stable_readonly_regular_at(
                publication_descriptor,
                result_target.name,
                "pressure-pilot aggregate analysis result",
                max_bytes=MAX_JSON_BYTES,
            )
            == result_payload,
            "pressure-pilot aggregate analysis result drifted after publication",
        )
        _require_same_directory(
            parent, publication_descriptor, "authorized PIC publication root"
        )
        _require(
            receipt_identity is not None,
            "pressure-pilot published receipt identity is absent",
        )
        _require_same_file_at(
            publication_descriptor,
            receipt_target.name,
            receipt_identity,
            "canonical pressure-pilot receipt",
        )
        _verify_published_pressure_pilot_receipt(
            receipt_target,
            authorized_pic_root=authorized_pic_root,
            allow_publication_guard=True,
            require_publication_seal=False,
        )
        _require_same_directory(
            acceptance_root,
            acceptance_descriptor,
            "authorized PIC publication acceptance root",
        )
        publication_result = {
            "path": str(target),
            "manifest_sha256": manifest_sha256,
            "analysis_result_path": str(result_target),
            "analysis_result_sha256": result_sha256,
            "receipt_path": str(receipt_target),
            "receipt_sha256": receipt_sha256,
        }
        seal_identity = _publish_publication_seal_at(
            acceptance_descriptor,
            publication_descriptor,
            receipt_target.name,
            receipt_payload,
            receipt_identity,
        )
        seal_committed = True
        _require(
            guard_identity is not None,
            "pressure-pilot publication guard identity is absent",
        )
        _disarm_publication_guard_at(
            publication_descriptor, receipt_target.name, guard_identity
        )
        guard_armed = False
        guard_identity = None
        _require_same_directory(
            transaction_anchor,
            transaction_descriptor,
            "stable publication transaction anchor",
        )
        _require_same_directory(
            parent, publication_descriptor, "authorized PIC publication root"
        )
        _require_same_directory(
            acceptance_root,
            acceptance_descriptor,
            "authorized PIC publication acceptance root",
        )
        _require_same_file_at(
            publication_descriptor,
            receipt_target.name,
            receipt_identity,
            "canonical pressure-pilot receipt",
        )
        _require_same_file_at(
            acceptance_descriptor,
            _publication_seal_name(receipt_target.name),
            seal_identity,
            "canonical pressure-pilot receipt durable success seal",
        )
        _require_publication_guard_absent_at(
            publication_descriptor,
            receipt_target.name,
            "canonical pressure-pilot receipt",
        )
        return publication_result
    except BaseException as publication_error:
        if (
            seal_committed
            and receipt_identity is not None
            and seal_identity is not None
        ):
            try:
                _fsync_descriptor(publication_descriptor)
                _require_same_directory(
                    transaction_anchor,
                    transaction_descriptor,
                    "stable publication transaction anchor",
                )
                _require_same_directory(
                    parent, publication_descriptor, "authorized PIC publication root"
                )
                _require_same_directory(
                    acceptance_root,
                    acceptance_descriptor,
                    "authorized PIC publication acceptance root",
                )
                _require_same_file_at(
                    publication_descriptor,
                    receipt_target.name,
                    receipt_identity,
                    "canonical pressure-pilot receipt",
                )
                _require_same_file_at(
                    acceptance_descriptor,
                    _publication_seal_name(receipt_target.name),
                    seal_identity,
                    "canonical pressure-pilot receipt durable success seal",
                )
                _require_publication_guard_absent_at(
                    publication_descriptor,
                    receipt_target.name,
                    "canonical pressure-pilot receipt",
                )
            except BaseException:
                pass
            else:
                guard_armed = False
                return publication_result
        if guard_armed or renamed or result_renamed or receipt_renamed:
            try:
                _require_same_directory(
                    transaction_anchor,
                    transaction_descriptor,
                    "stable publication transaction anchor",
                )
                guard_identity = _ensure_publication_guard_at(
                    publication_descriptor, receipt_target.name
                )
                guard_armed = True
            except BaseException as error:
                raise PressurePilotPublicationError(
                    "pressure-pilot public artifacts require reviewed reconciliation "
                    "and the fail-closed guard could not be assured"
                ) from error
            raise PressurePilotPublicationError(
                "pressure-pilot public artifacts retained under fail-closed guard; "
                f"reviewed reconciliation required: {publication_error}"
            ) from publication_error
        raise
    finally:
        close_error = _close_descriptors(
            (
                staging_descriptor,
                acceptance_descriptor,
                publication_descriptor,
                transaction_descriptor,
            )
        )
        if close_error is not None and not seal_committed:
            raise close_error


def _verify_published_pressure_pilot_receipt(
    receipt_path: str | Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    allow_publication_guard: bool = False,
    require_publication_seal: bool = True,
) -> dict[str, object]:
    """Re-audit one retained publication receipt and every artifact it binds."""
    pic_root, publication_root = _publication_root(authorized_pic_root)
    acceptance_root = _publication_acceptance_root(pic_root)
    receipt_target = _direct_publication_target(
        receipt_path, publication_root, "pressure-pilot receipt"
    )
    publication_descriptor = _open_absolute_directory(publication_root)
    acceptance_descriptor = _open_absolute_directory(acceptance_root)
    bundle_descriptor: int | None = None
    try:
        _require_same_directory(
            publication_root, publication_descriptor, "authorized PIC publication root"
        )
        _require_same_directory(
            acceptance_root,
            acceptance_descriptor,
            "authorized PIC publication acceptance root",
        )
        if not allow_publication_guard:
            _require_publication_guard_absent_at(
                publication_descriptor,
                receipt_target.name,
                "pressure-pilot receipt",
            )
        receipt_identity = _file_identity_at(
            publication_descriptor, receipt_target.name, "pressure-pilot receipt"
        )
        receipt_payload = _read_stable_readonly_regular_at(
            publication_descriptor,
            receipt_target.name,
            "pressure-pilot receipt",
            max_bytes=MAX_JSON_BYTES,
        )
        if require_publication_seal:
            _require_publication_seal_at(
                acceptance_descriptor,
                publication_descriptor,
                receipt_target.name,
                receipt_payload,
                receipt_identity,
                "pressure-pilot receipt",
            )
        try:
            receipt = pilot._decode_json(receipt_payload, "pressure-pilot receipt")
        except RecursionError as error:
            raise PressurePilotPublicationError(
                "pressure-pilot receipt exceeds the supported JSON nesting depth"
            ) from error
        _require(isinstance(receipt, dict), "pressure-pilot receipt is malformed")
        _require(
            set(receipt)
            == {
                "schema_version",
                "record_type",
                "evidence_class",
                "qualification_effect",
                "consumption_rule",
                "publication_root_identity",
                "aggregate_bundle",
                "aggregate_analysis",
                "source_bindings",
                "raw_cases",
            },
            "pressure-pilot receipt schema drifted",
        )
        _require(
            type(receipt["schema_version"]) is int
            and receipt["schema_version"] == 1
            and receipt["record_type"]
            == "q011_section54_pressure_pilot_bundle_publication_receipt"
            and receipt["evidence_class"] == pilot.EVIDENCE_CLASS
            and receipt["qualification_effect"] == pilot.QUALIFICATION_EFFECT
            and receipt["consumption_rule"] == CONSUMPTION_RULE,
            "pressure-pilot receipt identity drifted",
        )
        _require_directory_identity(
            receipt["publication_root_identity"],
            publication_descriptor,
            "pressure-pilot receipt publication root",
        )
        _validate_retained_source_bindings(
            receipt["source_bindings"],
            require_verified_archive=_is_production_pic_root(pic_root),
        )
        aggregate_bundle = receipt["aggregate_bundle"]
        aggregate_analysis = receipt["aggregate_analysis"]
        _require(
            isinstance(aggregate_bundle, dict)
            and set(aggregate_bundle) == {"path", "manifest_sha256"}
            and isinstance(aggregate_analysis, dict)
            and set(aggregate_analysis) == {"path", "sha256"},
            "pressure-pilot retained aggregate bindings are malformed",
        )
        bundle_root = _direct_publication_target(
            str(aggregate_bundle["path"]),
            publication_root,
            "pressure-pilot retained bundle",
        )
        analysis_result = _direct_publication_target(
            str(aggregate_analysis["path"]),
            publication_root,
            "pressure-pilot retained analysis",
        )
        manifest_sha256 = str(aggregate_bundle["manifest_sha256"])
        analysis_sha256 = str(aggregate_analysis["sha256"])
        _require(
            _SHA256_PATTERN.fullmatch(manifest_sha256) is not None
            and _SHA256_PATTERN.fullmatch(analysis_sha256) is not None,
            "pressure-pilot retained aggregate digest is malformed",
        )
        raw_cases = receipt["raw_cases"]
        _require(
            isinstance(raw_cases, list)
            and [case.get("case_id") for case in raw_cases if isinstance(case, dict)]
            == list(case_verifier.CASE_IDS),
            "pressure-pilot retained raw-case ordering drifted",
        )
        with ExitStack() as stack:
            for case in raw_cases:
                _require(
                    isinstance(case, dict)
                    and set(case)
                    == {
                        "case_id",
                        "artifact_dir",
                        "descriptor_path",
                        "descriptor_sha256",
                        "artifact_inventory_sha256",
                        "runtime_artifacts",
                    },
                    "pressure-pilot retained raw-case binding is malformed",
                )
                case_id = str(case["case_id"])
                raw_root = _authorized_raw_case_root(str(case["artifact_dir"]), pic_root)
                tree = stack.enter_context(case_verifier.StructuredArtifactTree(raw_root))
                descriptor = _verify_published_case_descriptor_bounded(
                    tree, case_id, str(case["descriptor_sha256"])
                )
                _require(
                    case["descriptor_path"] == case_verifier.CASE_DESCRIPTOR_PATH
                    and case["artifact_inventory_sha256"]
                    == descriptor["artifact_inventory_sha256"]
                    and case["runtime_artifacts"] == descriptor["runtime_artifacts"],
                    "pressure-pilot retained raw-case provenance drifted",
                )
        bundle_descriptor = os.open(
            bundle_root.name, _DIRECTORY_FLAGS, dir_fd=publication_descriptor
        )
        result = _verify_pressure_pilot_bundle_at(
            bundle_root,
            publication_descriptor,
            bundle_root.name,
            bundle_descriptor,
            manifest_sha256,
            authorized_publication_root=publication_root,
        )
        expected_result_payload = _canonical_json_bytes(result)
        actual_result_payload = _read_stable_readonly_regular_at(
            publication_descriptor,
            analysis_result.name,
            "pressure-pilot retained aggregate analysis",
            max_bytes=MAX_JSON_BYTES,
        )
        _require(
            actual_result_payload == expected_result_payload
            and _sha256(actual_result_payload) == analysis_sha256,
            "pressure-pilot retained aggregate analysis drifted",
        )
        _require_same_directory(
            publication_root, publication_descriptor, "authorized PIC publication root"
        )
        _require_same_directory(
            acceptance_root,
            acceptance_descriptor,
            "authorized PIC publication acceptance root",
        )
        if not allow_publication_guard:
            _require_publication_guard_absent_at(
                publication_descriptor,
                receipt_target.name,
                "pressure-pilot receipt",
            )
        if require_publication_seal:
            _require_publication_seal_at(
                acceptance_descriptor,
                publication_descriptor,
                receipt_target.name,
                receipt_payload,
                receipt_identity,
                "pressure-pilot receipt",
            )
        return {
            "receipt_sha256": _sha256(receipt_payload),
            "manifest_sha256": manifest_sha256,
            "analysis_result_sha256": analysis_sha256,
            "status": result["status"],
        }
    finally:
        if bundle_descriptor is not None:
            os.close(bundle_descriptor)
        os.close(acceptance_descriptor)
        os.close(publication_descriptor)


def verify_published_pressure_pilot_receipt(
    receipt_path: str | Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    """Re-audit one externally consumable retained publication receipt."""
    return _verify_published_pressure_pilot_receipt(
        receipt_path,
        authorized_pic_root=authorized_pic_root,
        allow_publication_guard=False,
        require_publication_seal=True,
    )


def consume_published_pressure_pilot_bundle(
    receipt_path: str | Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    """Accept one aggregate bundle only through its durable verified receipt."""
    return verify_published_pressure_pilot_receipt(
        receipt_path, authorized_pic_root=authorized_pic_root
    )


def _assignments(values: Sequence[str], label: str) -> dict[str, str]:
    result = {}
    for value in values:
        case_id, separator, item = value.partition("=")
        _require(bool(separator and case_id and item), f"{label} must use CASE_ID=VALUE")
        _require(case_id not in result, f"{label} repeats case ID: {case_id}")
        result[case_id] = item
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_path", nargs="?", type=Path)
    parser.add_argument("--receipt-path", type=Path)
    parser.add_argument("--analysis-result-path", type=Path)
    parser.add_argument("--verify-published-receipt", type=Path)
    parser.add_argument(
        "--case-artifact-dir",
        action="append",
        default=[],
        metavar="CASE_ID=PATH",
    )
    parser.add_argument(
        "--case-descriptor-sha256",
        action="append",
        default=[],
        metavar="CASE_ID=SHA256",
    )
    args = parser.parse_args(argv)
    if args.verify_published_receipt is not None:
        _require(
            args.output_path is None
            and args.receipt_path is None
            and args.analysis_result_path is None
            and not args.case_artifact_dir
            and not args.case_descriptor_sha256,
            "--verify-published-receipt cannot be combined with publication arguments",
        )
        verification = verify_published_pressure_pilot_receipt(
            args.verify_published_receipt
        )
        print(verification["receipt_sha256"])
        return 0
    _require(args.output_path is not None, "pressure-pilot output path is required")
    _require(args.receipt_path is not None, "--receipt-path is required")
    _require(args.analysis_result_path is not None, "--analysis-result-path is required")
    receipt = publish_pressure_pilot_bundle(
        args.output_path,
        receipt_path=args.receipt_path,
        analysis_result_path=args.analysis_result_path,
        case_artifact_dirs=_assignments(args.case_artifact_dir, "--case-artifact-dir"),
        case_descriptor_sha256=_assignments(
            args.case_descriptor_sha256, "--case-descriptor-sha256"
        ),
    )
    print(receipt["manifest_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
