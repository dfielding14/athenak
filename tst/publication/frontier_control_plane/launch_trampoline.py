#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Re-verify one reserved immutable snapshot before executing its job contract."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
import base64
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
from typing import Callable, Iterator
import uuid

from control_plane_common import AUTHORIZED_PIC_ROOT, AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import atomic_write_bytes_at
from control_plane_common import durable_mkdir_parents
from control_plane_common import PinnedDirectoryAncestry
from control_plane_common import launch_contract_sha256
from control_plane_common import read_json_bytes, record_for_role, require_ledger_paths
from control_plane_common import require_same_directory, validate_launch_contract
from control_plane_common import stable_serialization_anchor
from control_plane_common import utc_datetime
from control_plane_common import verify_snapshot_files
from ledger import append_primary_event_locked, latest_reservations, ledger_lock
from ledger import repair_mirrored_state_locked, transition_payload
from ledger import validate_mirrored_state
from validate_and_reserve_frontier_job import _require_run_artifact_dir
from validate_and_reserve_frontier_job import executable_reservation_bound_manifest
from validate_and_reserve_frontier_job import reservation_bound_manifest


SRUN = "/usr/bin/srun"
TRUSTED_PYTHON = "/opt/cray/pe/python/3.11.7/bin/python3"
_INPUT_DECK_FD_TOKEN = "__PIC_INPUT_DECK_FD__"
_DIRECTORY_OPEN_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_NEW_ARTIFACT_OPEN_FLAGS = (
    os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
)
_NEW_READ_WRITE_ARTIFACT_OPEN_FLAGS = (
    os.O_RDWR | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
)
_TIMEOUT_MARGIN_KEYS = {
    "athena_walltime_seconds",
    "scheduler_walltime_seconds",
    "environment_profile_sha256",
    "measured_utc",
    "expires_utc",
}
TRAMPOLINE_COMPLETION_NAMESPACE = Path("ledger/trampoline_completion_receipts")
TRAMPOLINE_COMPLETION_NAME = "trampoline_completion_receipt.json"
TRAMPOLINE_COMPLETION_RECORD_TYPE = "trusted_trampoline_completion_receipt"
TRAMPOLINE_COMPLETION_LEDGER_RECORD_TYPE = (
    "trusted_trampoline_completion_ledger_binding"
)
Q043_REGISTERED_CAMPAIGN = "q043_registered_execution_raw_oracle_successor_v1"
Q023_REGISTERED_CAMPAIGN = "q023_paper_bell_linear_joverc_registered_successor_v1"
Q019_REGISTERED_CAMPAIGN = "q019_nonlinear_bell_registered_successor_v1"
REGISTERED_WRAPPER_CAMPAIGNS = {
    Q043_REGISTERED_CAMPAIGN: ("Q043", r"q043-current-oracle-[a-z0-9_-]+"),
    Q023_REGISTERED_CAMPAIGN: ("Q023", r"[a-z0-9][a-z0-9_-]{0,127}"),
    Q019_REGISTERED_CAMPAIGN: ("Q019", r"q019-[a-z0-9][a-z0-9_-]{0,123}"),
}
_TASK_LOCAL_EXEC = r"""
import hashlib
import os
import re
import socket
import stat
import subprocess
import sys

ROOT, EXECUTABLE, EXECUTABLE_SHA256, INPUT_DECK, INPUT_DECK_SHA256, *ATHENA_ARGS = sys.argv[1:]
INPUT_TOKEN = "__PIC_INPUT_DECK_FD__"
DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
FILE_FLAGS = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)

def open_parent(path):
    root = os.path.abspath(ROOT)
    path = os.path.abspath(path)
    if os.path.commonpath([root, path]) != root or path == root:
        raise SystemExit("PIC task input is outside the authorized root")
    relative = os.path.relpath(os.path.dirname(path), root)
    descriptor = os.open(root, DIRECTORY_FLAGS)
    try:
        if relative != ".":
            for part in relative.split(os.sep):
                child = os.open(part, DIRECTORY_FLAGS, dir_fd=descriptor)
                os.close(descriptor)
                descriptor = child
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise

def require_same_parent(path, expected):
    actual = open_parent(path)
    try:
        expected_stat = os.fstat(expected)
        actual_stat = os.fstat(actual)
        if (expected_stat.st_dev, expected_stat.st_ino) != (actual_stat.st_dev, actual_stat.st_ino):
            raise SystemExit("PIC task input parent changed before execution")
    finally:
        os.close(actual)

def open_verified(path, expected_sha256, *, executable):
    if re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is None:
        raise SystemExit("PIC task input checksum is malformed")
    parent = open_parent(path)
    try:
        descriptor = os.open(os.path.basename(path), FILE_FLAGS, dir_fd=parent)
    except BaseException:
        os.close(parent)
        raise
    metadata = os.fstat(descriptor)
    if not stat.S_ISREG(metadata.st_mode) or metadata.st_mode & 0o222:
        raise SystemExit("PIC task input is not a read-only regular file")
    if executable and not metadata.st_mode & 0o111:
        raise SystemExit("PIC task executable is not executable")
    digest = hashlib.sha256()
    while True:
        data = os.read(descriptor, 1024 * 1024)
        if not data:
            break
        digest.update(data)
    if digest.hexdigest() != expected_sha256:
        raise SystemExit("PIC task input checksum mismatch")
    os.lseek(descriptor, 0, os.SEEK_SET)
    require_same_parent(path, parent)
    return parent, descriptor

executable_parent, executable_fd = open_verified(
    EXECUTABLE, EXECUTABLE_SHA256, executable=True
)
input_parent, input_fd = open_verified(INPUT_DECK, INPUT_DECK_SHA256, executable=False)
require_same_parent(EXECUTABLE, executable_parent)
require_same_parent(INPUT_DECK, input_parent)
rank = os.environ.get("SLURM_PROCID", "")
rocr_visible_devices = os.environ.get("ROCR_VISIBLE_DEVICES", "")
if re.fullmatch(r"[0-9]+", rank) is None:
    raise SystemExit("PIC task has no numeric Slurm rank")
if re.fullmatch(r"[0-9]+", rocr_visible_devices) is None:
    raise SystemExit("PIC task has no numeric ROCR_VISIBLE_DEVICES binding")
linked = subprocess.check_output(
    ["/usr/bin/ldd", f"/proc/self/fd/{executable_fd}"],
    text=True,
    pass_fds=(executable_fd,),
)
for library in ("libamdhip64", "libmpi_amd", "libmpi_gtl_hsa"):
    if re.search(rf"^\s*{library}[.]so(?:[.][0-9]+)*\s+=>\s+(?!not found\b)\S+", linked, re.MULTILINE) is None:
        raise SystemExit(f"PIC task executable is not linked against {library}")
print(
    "PIC trusted GPU launch: "
    f"rank={rank} host={socket.gethostname()} "
    f"ROCR_VISIBLE_DEVICES={rocr_visible_devices} "
    "linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa",
    flush=True,
)
if ATHENA_ARGS.count(INPUT_TOKEN) != 1:
    raise SystemExit("PIC task input-deck binding is malformed")
os.set_inheritable(input_fd, True)
ATHENA_ARGS = [
    f"/proc/self/fd/{input_fd}" if argument == INPUT_TOKEN else argument
    for argument in ATHENA_ARGS
]
os.execve(executable_fd, [EXECUTABLE, *ATHENA_ARGS], os.environ)
"""


def _canonical_sha256(value: object) -> str:
    payload = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


@contextmanager
def _deterministic_artifact_umask() -> Iterator[None]:
    inherited_umask = os.umask(0o022)
    try:
        yield
    finally:
        os.umask(inherited_umask)


class _PinnedSnapshot:
    """Retain one verified launch-node snapshot and its lexical parent."""

    def __init__(
        self,
        path: Path,
        *,
        root: Path,
        expected_sha256: str,
        require_executable: bool = False,
    ) -> None:
        self.path = Path(os.path.abspath(path))
        self.root = Path(os.path.abspath(root))
        self.expected_sha256 = expected_sha256
        self.require_executable = require_executable
        self.parent_ancestry: PinnedDirectoryAncestry | None = None
        self.parent_descriptor: int | None = None
        self.descriptor: int | None = None

    def __enter__(self) -> "_PinnedSnapshot":
        self.parent_ancestry = PinnedDirectoryAncestry(
            self.path.parent, root=self.root
        )
        self.parent_descriptor = self.parent_ancestry.descriptor
        try:
            self.descriptor = os.open(
                self.path.name,
                os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=self.parent_descriptor,
            )
            metadata = os.fstat(self.descriptor)
            if (
                not stat.S_ISREG(metadata.st_mode)
                or metadata.st_nlink != 1
                or metadata.st_mode & 0o222
            ):
                raise ValueError(f"Snapshot is not a read-only regular file: {self.path}")
            if self.require_executable and not metadata.st_mode & 0o111:
                raise ValueError(f"Snapshot executable is not executable: {self.path}")
            digest = hashlib.sha256()
            with os.fdopen(self.descriptor, "rb", closefd=False) as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
            if digest.hexdigest() != self.expected_sha256:
                raise ValueError(f"Snapshot checksum mismatch: {self.path}")
            self.require_lexical_parent()
            return self
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def require_lexical_parent(self) -> None:
        if self.parent_descriptor is None or self.parent_ancestry is None:
            raise ValueError("Pinned snapshot parent is not open")
        require_same_directory(self.path.parent, self.parent_descriptor, root=self.root)
        self.parent_ancestry.require_same()

    def read_bytes(self) -> bytes:
        if self.descriptor is None:
            raise ValueError("Pinned snapshot is not open")
        self.require_lexical_parent()
        os.lseek(self.descriptor, 0, os.SEEK_SET)
        with os.fdopen(self.descriptor, "rb", closefd=False) as stream:
            data = stream.read()
        if hashlib.sha256(data).hexdigest() != self.expected_sha256:
            raise ValueError(f"Snapshot checksum mismatch: {self.path}")
        self.require_lexical_parent()
        return data

    def __exit__(self, *_: object) -> None:
        if self.descriptor is not None:
            os.close(self.descriptor)
            self.descriptor = None
        if self.parent_ancestry is not None:
            self.parent_ancestry.close()
            self.parent_ancestry = None
            self.parent_descriptor = None


def _strict_json_equal(left: object, right: object) -> bool:
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return set(left) == set(right) and all(
            _strict_json_equal(left[key], right[key]) for key in left
        )
    if isinstance(left, list):
        return len(left) == len(right) and all(
            _strict_json_equal(left_value, right_value)
            for left_value, right_value in zip(left, right)
        )
    return left == right


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _trusted_athena_timeout_arguments(
    manifest: dict[str, object], *, root: Path
) -> tuple[str, str]:
    record = record_for_role(manifest, "timeout-margin")
    with _PinnedSnapshot(
        Path(str(record["path"])),
        root=root,
        expected_sha256=str(record["sha256"]),
    ) as timeout_margin:
        margin = read_json_bytes(
            timeout_margin.read_bytes(), label="Timeout-margin snapshot"
        )
    if (
        set(margin) != _TIMEOUT_MARGIN_KEYS
        or not _strict_json_equal(margin, manifest.get("timeout_margin"))
    ):
        raise ValueError("Timeout-margin snapshot differs from manifest record")
    scheduler = margin.get("scheduler_walltime_seconds")
    athena = margin.get("athena_walltime_seconds")
    if type(scheduler) is not int or type(athena) is not int:
        raise ValueError("Timeout-margin walltimes must be exact integers")
    if not 0 < athena < scheduler:
        raise ValueError("Athena timeout must be positive and below Slurm walltime")
    profile = record_for_role(manifest, "environment-profile")
    if margin.get("environment_profile_sha256") != profile.get("sha256"):
        raise ValueError("Timeout-margin artifact does not match environment profile")
    measured = utc_datetime(margin.get("measured_utc"), field="measured_utc")
    expires = utc_datetime(margin.get("expires_utc"), field="expires_utc")
    if not measured <= _utc_now() < expires:
        raise ValueError("Timeout-margin artifact is stale or not yet valid")
    hours, remainder = divmod(athena, 60 * 60)
    minutes, seconds = divmod(remainder, 60)
    return "-t", f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def _require_no_action_timeout_override(value: object) -> None:
    if not isinstance(value, dict) or not isinstance(value.get("actions"), list):
        return
    for action in value["actions"]:
        if not isinstance(action, dict) or not isinstance(action.get("arguments"), list):
            continue
        for argument in action["arguments"]:
            if not isinstance(argument, dict):
                continue
            literal = argument.get("literal")
            if isinstance(literal, str) and (
                literal == "-t" or literal.startswith("-t=")
            ):
                raise ValueError("Launch-action timeout override is not authorized")


def _artifact_path(artifact_dir: Path, relative: object) -> Path:
    text = str(relative)
    path = PurePosixPath(text)
    if (
        not text
        or path.is_absolute()
        or not path.parts
        or path.parts != tuple(part for part in path.parts if part not in {"", ".", ".."})
    ):
        raise ValueError("Artifact path must be a non-empty relative path")
    return artifact_dir.joinpath(*path.parts)


def _artifact_relative_parts(artifact_dir: Path, path: Path) -> tuple[str, ...]:
    try:
        parts = path.relative_to(artifact_dir).parts
    except ValueError as error:
        raise ValueError(f"Artifact path is outside launch artifact directory: {path}") from error
    if any(part in {"", ".", ".."} for part in parts):
        raise ValueError(f"Artifact path contains an unsafe component: {path}")
    return parts


def _directory_identity(metadata: os.stat_result) -> tuple[int, int]:
    return metadata.st_dev, metadata.st_ino


def _open_created_directory_at(
    parent_fd: int, name: str, *, mode: int, label: str
) -> int:
    # POSIX mkdirat does not return a descriptor. Same-UID process isolation is
    # therefore an operational prerequisite until the first no-follow open.
    try:
        os.mkdir(name, mode=mode, dir_fd=parent_fd)
    except FileExistsError as error:
        raise ValueError(f"{label} already exists: {name}") from error
    created = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    descriptor = os.open(name, _DIRECTORY_OPEN_FLAGS, dir_fd=parent_fd)
    try:
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISDIR(created.st_mode)
            or not stat.S_ISDIR(opened.st_mode)
            or _directory_identity(created) != _directory_identity(opened)
        ):
            raise ValueError(f"{label} changed between creation and open: {name}")
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _open_artifact_directory(
    artifact_dir_fd: int,
    artifact_dir: Path,
    directory: Path,
    *,
    create: bool,
    directory_identities: dict[str, tuple[int, int]] | None = None,
) -> int:
    descriptor = os.dup(artifact_dir_fd)
    try:
        for index, part in enumerate(_artifact_relative_parts(artifact_dir, directory)):
            created = False
            try:
                child_descriptor = os.open(part, _DIRECTORY_OPEN_FLAGS, dir_fd=descriptor)
            except FileNotFoundError:
                if not create:
                    raise
                child_descriptor = _open_created_directory_at(
                    descriptor,
                    part,
                    mode=0o755,
                    label="Launch artifact directory",
                )
                created = True
            try:
                child = os.fstat(child_descriptor)
                relative = "/".join(
                    _artifact_relative_parts(artifact_dir, directory)[: index + 1]
                )
                identity = _directory_identity(child)
                if directory_identities is not None:
                    expected = directory_identities.setdefault(relative, identity)
                    if expected != identity:
                        raise ValueError(
                            f"Launch artifact directory changed during execution: {relative}"
                        )
                if created:
                    os.fsync(child_descriptor)
                    os.fsync(descriptor)
            except BaseException:
                os.close(child_descriptor)
                raise
            os.close(descriptor)
            descriptor = child_descriptor
    except BaseException:
        os.close(descriptor)
        raise
    return descriptor


def _mkdir_artifact_directory(
    artifact_dir_fd: int,
    artifact_dir: Path,
    directory: Path,
    directory_identities: dict[str, tuple[int, int]],
) -> None:
    descriptor = _open_artifact_directory(
        artifact_dir_fd,
        artifact_dir,
        directory,
        create=True,
        directory_identities=directory_identities,
    )
    os.close(descriptor)


def _require_retained_artifact_directories_at(
    artifact_dir_fd: int,
    artifact_dir: Path,
    directory_identities: dict[str, tuple[int, int]],
) -> None:
    for relative in sorted(directory_identities):
        descriptor = _open_artifact_directory(
            artifact_dir_fd,
            artifact_dir,
            artifact_dir.joinpath(*PurePosixPath(relative).parts),
            create=False,
            directory_identities=directory_identities,
        )
        os.close(descriptor)


def _capture_artifact_directory_identities_at(
    directory_fd: int,
    directory_identities: dict[str, tuple[int, int]],
    prefix: tuple[str, ...] = (),
) -> None:
    for name in sorted(os.listdir(directory_fd)):
        if not prefix and name == "analysis":
            continue
        metadata = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if not stat.S_ISDIR(metadata.st_mode):
            continue
        relative = "/".join((*prefix, name))
        child_fd = os.open(name, _DIRECTORY_OPEN_FLAGS, dir_fd=directory_fd)
        try:
            child = os.fstat(child_fd)
            identity = _directory_identity(child)
            if _directory_identity(metadata) != identity:
                raise ValueError(
                    f"Launch artifact directory changed during capture: {relative}"
                )
            expected = directory_identities.setdefault(relative, identity)
            if expected != identity:
                raise ValueError(
                    f"Launch artifact directory changed during execution: {relative}"
                )
            _capture_artifact_directory_identities_at(
                child_fd,
                directory_identities,
                (*prefix, name),
            )
        finally:
            os.close(child_fd)


def _open_artifact_file(
    artifact_dir_fd: int,
    artifact_dir: Path,
    path: Path,
    flags: int,
    *,
    create_parent: bool,
    mode: int = 0o600,
) -> int:
    parts = _artifact_relative_parts(artifact_dir, path)
    if not parts:
        raise ValueError("Artifact file path must name one file")
    parent_descriptor = _open_artifact_directory(
        artifact_dir_fd, artifact_dir, path.parent, create=create_parent
    )
    try:
        descriptor = os.open(parts[-1], flags, mode, dir_fd=parent_descriptor)
    finally:
        os.close(parent_descriptor)
    try:
        regular = stat.S_ISREG(os.fstat(descriptor).st_mode)
    except BaseException:
        os.close(descriptor)
        raise
    if not regular:
        os.close(descriptor)
        raise ValueError(f"Launch artifact is not a regular file: {path}")
    return descriptor


def _open_new_artifact(artifact_dir_fd: int, artifact_dir: Path, path: Path) -> int:
    return _open_artifact_file(
        artifact_dir_fd,
        artifact_dir,
        path,
        _NEW_ARTIFACT_OPEN_FLAGS,
        create_parent=True,
    )


def _read_artifact_bytes(artifact_dir_fd: int, artifact_dir: Path, path: Path) -> bytes:
    descriptor = _open_artifact_file(
        artifact_dir_fd,
        artifact_dir,
        path,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        create_parent=False,
    )
    try:
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            return stream.read()
    finally:
        os.close(descriptor)


def _artifact_size(artifact_dir_fd: int, artifact_dir: Path, path: Path) -> int:
    descriptor = _open_artifact_file(
        artifact_dir_fd,
        artifact_dir,
        path,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        create_parent=False,
    )
    try:
        return os.fstat(descriptor).st_size
    finally:
        os.close(descriptor)


def _write_new_text_artifact(
    artifact_dir_fd: int, artifact_dir: Path, path: Path, text: str
) -> None:
    descriptor = _open_new_artifact(artifact_dir_fd, artifact_dir, path)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", closefd=False) as stream:
            stream.write(text)
    finally:
        os.close(descriptor)


def _freeze_artifact_file_at(
    directory_fd: int,
    name: str,
    relative: str,
    *,
    expected_metadata: os.stat_result | None = None,
    file_identities: dict[str, tuple[int, int]] | None = None,
) -> dict[str, object]:
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_fd,
    )
    try:
        return _freeze_open_artifact_file_at(
            directory_fd,
            name,
            relative,
            descriptor,
            expected_metadata=expected_metadata,
            file_identities=file_identities,
        )
    finally:
        os.close(descriptor)


def _freeze_open_artifact_file_at(
    directory_fd: int,
    name: str,
    relative: str,
    descriptor: int,
    *,
    expected_metadata: os.stat_result | None = None,
    file_identities: dict[str, tuple[int, int]] | None = None,
) -> dict[str, object]:
    try:
        initial = os.fstat(descriptor)
        if (
            expected_metadata is not None
            and _directory_identity(expected_metadata) != _directory_identity(initial)
        ):
            raise ValueError(f"Launch artifact changed before freezing: {relative}")
        if not stat.S_ISREG(initial.st_mode):
            raise ValueError(f"Launch artifact is not a regular file: {relative}")
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
        os.lseek(descriptor, 0, os.SEEK_SET)
        before = os.fstat(descriptor)
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            data = stream.read()
        after_read = os.fstat(descriptor)
        stable = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns")
        if any(getattr(before, field) != getattr(after_read, field) for field in stable):
            raise ValueError(f"Launch artifact changed while freezing: {relative}")
        after = os.fstat(descriptor)
        entry = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if (
            (entry.st_dev, entry.st_ino) != (after.st_dev, after.st_ino)
            or len(data) != after.st_size
            or after.st_mode & 0o222
        ):
            raise ValueError(f"Launch artifact changed while freezing: {relative}")
        if file_identities is not None:
            identity = (after.st_dev, after.st_ino)
            expected = file_identities.setdefault(relative, identity)
            if expected != identity:
                raise ValueError(f"Launch artifact changed while freezing: {relative}")
        return {
            "path": relative,
            "sha256": hashlib.sha256(data).hexdigest(),
            "size": len(data),
        }
    except OSError as error:
        raise ValueError(f"Launch artifact cannot be frozen: {relative}") from error


def _freeze_artifact_tree_at(
    directory_fd: int,
    prefix: tuple[str, ...] = (),
    directory_identities: dict[str, tuple[int, int]] | None = None,
    file_identities: dict[str, tuple[int, int]] | None = None,
) -> list[dict[str, object]]:
    records = []
    for name in sorted(os.listdir(directory_fd)):
        if not prefix and name == "analysis":
            continue
        relative = "/".join((*prefix, name))
        metadata = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if stat.S_ISREG(metadata.st_mode):
            records.append(
                _freeze_artifact_file_at(
                    directory_fd,
                    name,
                    relative,
                    expected_metadata=metadata,
                    file_identities=file_identities,
                )
            )
        elif stat.S_ISDIR(metadata.st_mode):
            child_fd = os.open(name, _DIRECTORY_OPEN_FLAGS, dir_fd=directory_fd)
            try:
                child = os.fstat(child_fd)
                if _directory_identity(metadata) != _directory_identity(child):
                    raise ValueError(
                        f"Launch artifact directory changed before freezing: {relative}"
                    )
                child_records = _freeze_artifact_tree_at(
                    child_fd,
                    (*prefix, name),
                    directory_identities,
                    file_identities,
                )
                if not child_records:
                    raise ValueError(
                        f"Launch artifact tree contains an empty directory: {relative}"
                    )
                records.extend(child_records)
                os.fchmod(child_fd, 0o555)
                os.fsync(child_fd)
                after = os.fstat(child_fd)
                entry = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
                if (entry.st_dev, entry.st_ino) != (after.st_dev, after.st_ino):
                    raise ValueError(
                        f"Launch artifact directory changed while freezing: {relative}"
                    )
                if directory_identities is not None:
                    identity = (after.st_dev, after.st_ino)
                    expected = directory_identities.setdefault(relative, identity)
                    if expected != identity:
                        raise ValueError(
                            f"Launch artifact directory changed before freezing: {relative}"
                        )
            finally:
                os.close(child_fd)
        else:
            raise ValueError(f"Launch artifact tree contains an unsupported entry: {relative}")
    return records


def _verify_frozen_artifact_file_at(
    directory_fd: int,
    name: str,
    relative: str,
    expected: dict[str, object],
    *,
    expected_identity: tuple[int, int] | None = None,
) -> None:
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_fd,
    )
    try:
        _verify_open_frozen_artifact_file_at(
            directory_fd,
            name,
            relative,
            expected,
            descriptor,
            expected_identity=expected_identity,
        )
    finally:
        os.close(descriptor)


def _verify_open_frozen_artifact_file_at(
    directory_fd: int,
    name: str,
    relative: str,
    expected: dict[str, object],
    descriptor: int,
    *,
    expected_identity: tuple[int, int] | None = None,
) -> None:
    os.lseek(descriptor, 0, os.SEEK_SET)
    before = os.fstat(descriptor)
    if not stat.S_ISREG(before.st_mode) or before.st_mode & 0o222:
        raise ValueError(f"Launch artifact is not frozen: {relative}")
    with os.fdopen(descriptor, "rb", closefd=False) as stream:
        data = stream.read()
    after = os.fstat(descriptor)
    entry = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    stable = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns")
    if (
        any(getattr(before, field) != getattr(after, field) for field in stable)
        or (
            expected_identity is not None
            and expected_identity != (after.st_dev, after.st_ino)
        )
        or (entry.st_dev, entry.st_ino) != (after.st_dev, after.st_ino)
        or len(data) != after.st_size
        or expected.get("path") != relative
        or expected.get("size") != len(data)
        or expected.get("sha256") != hashlib.sha256(data).hexdigest()
    ):
        raise ValueError(f"Launch artifact changed after freezing: {relative}")


def _verify_frozen_artifact_tree_at(
    directory_fd: int,
    records: list[dict[str, object]],
    directory_identities: dict[str, tuple[int, int]],
    analysis_fd: int,
    file_identities: dict[str, tuple[int, int]] | None = None,
    prefix: tuple[str, ...] = (),
) -> None:
    expected = {str(record["path"]): record for record in records}
    all_records = dict(expected)
    expected_directories = dict(directory_identities)
    expected_files = None if file_identities is None else dict(file_identities)
    retained_descriptors: list[int] = []
    retained_directory_entries: list[tuple[int, str, str, int]] = []
    retained_directory_names: list[tuple[int, tuple[str, ...], list[str]]] = []
    retained_file_entries: list[
        tuple[int, str, str, dict[str, object], tuple[int, int] | None, int]
    ] = []

    def verify_tree(current_fd: int, current_prefix: tuple[str, ...]) -> None:
        before = os.fstat(current_fd)
        if not stat.S_ISDIR(before.st_mode) or before.st_mode & 0o222:
            relative = "/".join(current_prefix) or "."
            raise ValueError(f"Launch artifact directory is not frozen: {relative}")
        names = sorted(os.listdir(current_fd))
        for name in names:
            if not current_prefix and name == "analysis":
                analysis = os.fstat(analysis_fd)
                entry = os.stat(name, dir_fd=current_fd, follow_symlinks=False)
                if (
                    not stat.S_ISDIR(analysis.st_mode)
                    or stat.S_IMODE(analysis.st_mode) != 0o700
                    or stat.S_IMODE(entry.st_mode) != 0o700
                    or os.listdir(analysis_fd)
                    or (entry.st_dev, entry.st_ino)
                    != (analysis.st_dev, analysis.st_ino)
                ):
                    raise ValueError("Launch artifact analysis directory changed while freezing")
                continue
            if not current_prefix and name == "artifact_inventory.json":
                continue
            relative = "/".join((*current_prefix, name))
            metadata = os.stat(name, dir_fd=current_fd, follow_symlinks=False)
            if stat.S_ISREG(metadata.st_mode):
                record = expected.pop(relative, None)
                if record is None:
                    raise ValueError(f"Launch artifact tree gained an unlisted file: {relative}")
                expected_identity = (
                    None if expected_files is None else expected_files.pop(relative, None)
                )
                if file_identities is not None and expected_identity is None:
                    raise ValueError(f"Launch artifact tree gained an unbound file: {relative}")
                descriptor = os.open(
                    name,
                    os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=current_fd,
                )
                retained_descriptors.append(descriptor)
                _verify_open_frozen_artifact_file_at(
                    current_fd,
                    name,
                    relative,
                    record,
                    descriptor,
                    expected_identity=expected_identity,
                )
                retained_file_entries.append(
                    (
                        current_fd,
                        name,
                        relative,
                        record,
                        expected_identity,
                        descriptor,
                    )
                )
            elif stat.S_ISDIR(metadata.st_mode):
                child_fd = os.open(name, _DIRECTORY_OPEN_FLAGS, dir_fd=current_fd)
                retained_descriptors.append(child_fd)
                child = os.fstat(child_fd)
                expected_identity = expected_directories.pop(relative, None)
                if expected_identity != (child.st_dev, child.st_ino):
                    raise ValueError(
                        f"Launch artifact directory changed after freezing: {relative}"
                    )
                retained_directory_entries.append((current_fd, name, relative, child_fd))
                verify_tree(child_fd, (*current_prefix, name))
                after = os.fstat(child_fd)
                entry = os.stat(name, dir_fd=current_fd, follow_symlinks=False)
                if (entry.st_dev, entry.st_ino) != (after.st_dev, after.st_ino):
                    raise ValueError(
                        f"Launch artifact directory changed after freezing: {relative}"
                    )
            else:
                raise ValueError(f"Launch artifact tree contains an unsupported entry: {relative}")
        after = os.fstat(current_fd)
        if (
            (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino)
            or sorted(os.listdir(current_fd)) != names
        ):
            relative = "/".join(current_prefix) or "."
            raise ValueError(f"Launch artifact directory changed after freezing: {relative}")
        retained_directory_names.append((current_fd, current_prefix, names))

    try:
        verify_tree(directory_fd, prefix)
        if expected:
            raise ValueError("Launch artifact tree lost a frozen file")
        if expected_directories:
            raise ValueError("Launch artifact tree lost a frozen directory")
        if expected_files:
            raise ValueError("Launch artifact tree lost a frozen file identity")
        def verify_retained_directories() -> None:
            for current_fd, current_prefix, names in retained_directory_names:
                metadata = os.fstat(current_fd)
                if (
                    not stat.S_ISDIR(metadata.st_mode)
                    or metadata.st_mode & 0o222
                    or sorted(os.listdir(current_fd)) != names
                ):
                    relative = "/".join(current_prefix) or "."
                    raise ValueError(
                        f"Launch artifact directory changed after freezing: {relative}"
                    )
            for parent_fd, name, relative, descriptor in retained_directory_entries:
                metadata = os.fstat(descriptor)
                entry = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
                if (
                    not stat.S_ISDIR(metadata.st_mode)
                    or metadata.st_mode & 0o222
                    or (entry.st_dev, entry.st_ino) != (metadata.st_dev, metadata.st_ino)
                ):
                    raise ValueError(
                        f"Launch artifact directory changed after freezing: {relative}"
                    )

        verify_retained_directories()
        for parent_fd, name, relative, record, identity, descriptor in retained_file_entries:
            _verify_open_frozen_artifact_file_at(
                parent_fd,
                name,
                relative,
                all_records[relative],
                descriptor,
                expected_identity=identity,
            )
        verify_retained_directories()
        analysis = os.fstat(analysis_fd)
        entry = os.stat("analysis", dir_fd=directory_fd, follow_symlinks=False)
        if (
            not stat.S_ISDIR(analysis.st_mode)
            or stat.S_IMODE(analysis.st_mode) != 0o700
            or stat.S_IMODE(entry.st_mode) != 0o700
            or os.listdir(analysis_fd)
            or (entry.st_dev, entry.st_ino) != (analysis.st_dev, analysis.st_ino)
        ):
            raise ValueError("Launch artifact analysis directory changed while freezing")
    finally:
        for descriptor in reversed(retained_descriptors):
            os.close(descriptor)


def _publish_empty_analysis_directory_at(artifact_dir_fd: int) -> int:
    try:
        os.stat("analysis", dir_fd=artifact_dir_fd, follow_symlinks=False)
    except FileNotFoundError:
        pass
    else:
        raise ValueError("Launch artifact analysis directory already exists")
    staging_name = f".analysis.staging-{uuid.uuid4()}"
    analysis_fd: int | None = None
    try:
        analysis_fd = _open_created_directory_at(
            artifact_dir_fd,
            staging_name,
            mode=0o700,
            label="Launch artifact staged analysis directory",
        )
        os.fchmod(analysis_fd, 0o700)
        os.fsync(analysis_fd)
        before = os.fstat(analysis_fd)
        entry = os.stat(staging_name, dir_fd=artifact_dir_fd, follow_symlinks=False)
        if (
            not stat.S_ISDIR(before.st_mode)
            or stat.S_IMODE(before.st_mode) != 0o700
            or os.listdir(analysis_fd)
            or (entry.st_dev, entry.st_ino) != (before.st_dev, before.st_ino)
        ):
            raise ValueError("Launch artifact staged analysis directory changed")
        os.rename(
            staging_name,
            "analysis",
            src_dir_fd=artifact_dir_fd,
            dst_dir_fd=artifact_dir_fd,
        )
        os.fsync(artifact_dir_fd)
        after = os.fstat(analysis_fd)
        entry = os.stat("analysis", dir_fd=artifact_dir_fd, follow_symlinks=False)
        if (
            not stat.S_ISDIR(after.st_mode)
            or stat.S_IMODE(after.st_mode) != 0o700
            or os.listdir(analysis_fd)
            or (entry.st_dev, entry.st_ino) != (after.st_dev, after.st_ino)
        ):
            raise ValueError("Launch artifact analysis directory changed during publication")
        return analysis_fd
    except BaseException:
        if analysis_fd is not None:
            os.close(analysis_fd)
        raise


def _publish_frozen_artifact_inventory_at(
    artifact_dir_fd: int,
    artifact_dir: Path,
    *,
    directory_identities: dict[str, tuple[int, int]] | None = None,
) -> None:
    try:
        analysis_fd = _publish_empty_analysis_directory_at(artifact_dir_fd)
    except FileExistsError as error:
        raise ValueError("Launch artifact analysis directory already exists") from error
    try:
        frozen_directory_identities = (
            {} if directory_identities is None else dict(directory_identities)
        )
        frozen_file_identities: dict[str, tuple[int, int]] = {}
        records = _freeze_artifact_tree_at(
            artifact_dir_fd,
            directory_identities=frozen_directory_identities,
            file_identities=frozen_file_identities,
        )
        inventory_descriptor = _open_artifact_file(
            artifact_dir_fd,
            artifact_dir,
            artifact_dir / "artifact_inventory.json",
            _NEW_READ_WRITE_ARTIFACT_OPEN_FLAGS,
            create_parent=False,
        )
        try:
            inventory_text = json.dumps(
                {"schema_version": 1, "files": records},
                indent=2,
                sort_keys=True,
            ) + "\n"
            with os.fdopen(
                inventory_descriptor, "w", encoding="utf-8", closefd=False
            ) as stream:
                stream.write(inventory_text)
            inventory_record = _freeze_open_artifact_file_at(
                artifact_dir_fd,
                "artifact_inventory.json",
                "artifact_inventory.json",
                inventory_descriptor,
            )
            os.fchmod(artifact_dir_fd, 0o555)
            os.fsync(artifact_dir_fd)
            _verify_open_frozen_artifact_file_at(
                artifact_dir_fd,
                "artifact_inventory.json",
                "artifact_inventory.json",
                inventory_record,
                inventory_descriptor,
            )
            _verify_frozen_artifact_tree_at(
                artifact_dir_fd,
                records,
                frozen_directory_identities,
                analysis_fd,
                frozen_file_identities,
            )
            _verify_open_frozen_artifact_file_at(
                artifact_dir_fd,
                "artifact_inventory.json",
                "artifact_inventory.json",
                inventory_record,
                inventory_descriptor,
            )
        finally:
            os.close(inventory_descriptor)
    finally:
        os.close(analysis_fd)


def _publish_frozen_artifact_inventory(
    artifact_dir_fd: int,
    artifact_dir: Path,
    *,
    directory_identities: dict[str, tuple[int, int]] | None = None,
) -> None:
    with _deterministic_artifact_umask():
        _publish_frozen_artifact_inventory_at(
            artifact_dir_fd,
            artifact_dir,
            directory_identities=directory_identities,
        )


def _open_frozen_artifact_member_at(
    artifact_dir_fd: int, relative: str
) -> tuple[int, int, str]:
    parts = PurePosixPath(relative).parts
    if (
        not relative
        or PurePosixPath(relative).is_absolute()
        or PurePosixPath(relative).as_posix() != relative
        or any(part in {"", ".", ".."} for part in parts)
    ):
        raise ValueError("Trampoline completion artifact path is unsafe")
    parent_fd = os.dup(artifact_dir_fd)
    try:
        for part in parts[:-1]:
            child_fd = os.open(part, _DIRECTORY_OPEN_FLAGS, dir_fd=parent_fd)
            child = os.fstat(child_fd)
            entry = os.stat(part, dir_fd=parent_fd, follow_symlinks=False)
            if (
                not stat.S_ISDIR(child.st_mode)
                or child.st_mode & 0o222
                or (entry.st_dev, entry.st_ino) != (child.st_dev, child.st_ino)
            ):
                os.close(child_fd)
                raise ValueError(
                    f"Trampoline completion artifact directory changed: {relative}"
                )
            os.close(parent_fd)
            parent_fd = child_fd
        descriptor = os.open(
            parts[-1],
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_fd,
        )
        return parent_fd, descriptor, parts[-1]
    except BaseException:
        os.close(parent_fd)
        raise


def _read_retained_completion_artifact(
    parent_fd: int,
    descriptor: int,
    name: str,
    *,
    relative: str,
    expected_sha256: str | None = None,
    expected_size: int | None = None,
) -> tuple[bytes, dict[str, int]]:
    os.lseek(descriptor, 0, os.SEEK_SET)
    before = os.fstat(descriptor)
    if (
        not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or before.st_mode & 0o222
    ):
        raise ValueError(f"Trampoline completion artifact is not immutable: {relative}")
    payload = bytearray()
    while chunk := os.read(descriptor, 1024 * 1024):
        payload.extend(chunk)
    after = os.fstat(descriptor)
    entry = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    stable = (
        "st_dev",
        "st_ino",
        "st_mode",
        "st_nlink",
        "st_size",
        "st_mtime_ns",
        "st_ctime_ns",
    )
    digest = hashlib.sha256(payload).hexdigest()
    if (
        any(getattr(before, field) != getattr(after, field) for field in stable)
        or (entry.st_dev, entry.st_ino) != (after.st_dev, after.st_ino)
        or len(payload) != after.st_size
        or (expected_sha256 is not None and digest != expected_sha256)
        or (expected_size is not None and len(payload) != expected_size)
    ):
        raise ValueError(f"Trampoline completion artifact changed: {relative}")
    return bytes(payload), {"device": after.st_dev, "inode": after.st_ino}


def _completion_receipt_paths(
    submission_id: str,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> tuple[Path, Path]:
    relative = (
        TRAMPOLINE_COMPLETION_NAMESPACE
        / submission_id
        / TRAMPOLINE_COMPLETION_NAME
    )
    return authorized_pic_root / relative, authorized_project_home_root / relative


def _publish_exact_completion_receipt_at(
    path: Path,
    payload: bytes,
    *,
    ancestry: PinnedDirectoryAncestry,
) -> tuple[int, int]:
    metadata = os.fstat(ancestry.descriptor)
    if (
        not stat.S_ISDIR(metadata.st_mode)
        or stat.S_IMODE(metadata.st_mode) not in {0o500, 0o700}
    ):
        raise ValueError("Trampoline completion receipt directory is not private")
    try:
        descriptor = os.open(
            path.name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=ancestry.descriptor,
        )
    except FileNotFoundError:
        atomic_write_bytes_at(
            ancestry.descriptor,
            path.name,
            payload,
            mode=0o444,
            replace=False,
            post_publish_check=ancestry.require_same,
        )
        descriptor = os.open(
            path.name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=ancestry.descriptor,
        )
    try:
        observed, identity = _read_retained_completion_artifact(
            ancestry.descriptor,
            descriptor,
            path.name,
            relative=str(path),
            expected_sha256=hashlib.sha256(payload).hexdigest(),
            expected_size=len(payload),
        )
        if observed != payload:
            raise ValueError("Trampoline completion receipt retry bytes differ")
    finally:
        os.close(descriptor)
    os.fchmod(ancestry.descriptor, 0o500)
    os.fsync(ancestry.descriptor)
    ancestry.require_same()
    if os.fstat(ancestry.descriptor).st_mode & 0o222:
        raise ValueError("Trampoline completion receipt directory remains mutable")
    return identity["device"], identity["inode"]


def _prepare_completion_receipt_ancestry(
    path: Path, *, root: Path
) -> PinnedDirectoryAncestry:
    durable_mkdir_parents(path.parent, mode=0o700, root=root)
    ancestry = PinnedDirectoryAncestry(
        path.parent, root=Path(os.path.abspath(root))
    )
    try:
        metadata = os.fstat(ancestry.descriptor)
        if (
            not stat.S_ISDIR(metadata.st_mode)
            or stat.S_IMODE(metadata.st_mode) not in {0o500, 0o700}
        ):
            raise ValueError("Trampoline completion receipt directory is not private")
        ancestry.require_same()
        return ancestry
    except BaseException:
        ancestry.close()
        raise


def _publish_trampoline_completion_receipt_at(
    artifact_dir_fd: int,
    artifact_dir: Path,
    manifest: dict[str, object],
    *,
    append_anchor: Callable[[dict[str, object]], dict[str, object]],
    manifest_path: Path,
    manifest_sha256: str,
    reservation_id: str,
    submission_id: str,
    slurm_job_id: str,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> dict[str, object]:
    """Anchor then publish paired exact receipts while the run tree is pinned."""
    trusted_anchor = stable_serialization_anchor(authorized_pic_root)
    require_same_directory(artifact_dir, artifact_dir_fd, root=trusted_anchor)
    root_metadata = os.fstat(artifact_dir_fd)
    if (
        not stat.S_ISDIR(root_metadata.st_mode)
        or root_metadata.st_mode & 0o222
    ):
        raise ValueError("Trampoline completion run root is not frozen")

    inventory_parent_fd, inventory_fd, inventory_name = _open_frozen_artifact_member_at(
        artifact_dir_fd, "artifact_inventory.json"
    )
    retained: list[tuple[int, int, str, str, dict[str, object]]] = []
    try:
        inventory_payload, inventory_identity = _read_retained_completion_artifact(
            inventory_parent_fd,
            inventory_fd,
            inventory_name,
            relative="artifact_inventory.json",
        )
        inventory = read_json_bytes(
            inventory_payload, label="trampoline completion artifact inventory"
        )
        records = inventory.get("files")
        if (
            set(inventory) != {"schema_version", "files"}
            or inventory.get("schema_version") != 1
            or not isinstance(records, list)
        ):
            raise ValueError("Trampoline completion artifact inventory is malformed")
        completion_records = []
        for raw in records:
            if (
                not isinstance(raw, dict)
                or set(raw) != {"path", "sha256", "size"}
                or not isinstance(raw["path"], str)
                or not isinstance(raw["sha256"], str)
                or re.fullmatch(r"[0-9a-f]{64}", raw["sha256"]) is None
                or type(raw["size"]) is not int
                or raw["size"] < 0
            ):
                raise ValueError("Trampoline completion artifact record is malformed")
            parent_fd, descriptor, name = _open_frozen_artifact_member_at(
                artifact_dir_fd, raw["path"]
            )
            retained.append((parent_fd, descriptor, name, raw["path"], raw))
            _, identity = _read_retained_completion_artifact(
                parent_fd,
                descriptor,
                name,
                relative=raw["path"],
                expected_sha256=raw["sha256"],
                expected_size=raw["size"],
            )
            completion_records.append(
                {
                    "path": raw["path"],
                    "sha256": raw["sha256"],
                    "byte_count": raw["size"],
                    "filesystem_identity": identity,
                }
            )
        contract = validate_launch_contract(manifest.get("launch_contract"))
        mandatory_stdio = sorted(
            {
                str(action[key])
                for action in contract["actions"]
                for key in ("stdout_artifact", "stderr_artifact")
            }
        )
        indexed = {record["path"]: record for record in completion_records}
        if any(path not in indexed for path in mandatory_stdio):
            raise ValueError("Trampoline completion receipt lacks mandatory stdout/stderr")
        orion_path, project_home_path = _completion_receipt_paths(
            submission_id,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        with _prepare_completion_receipt_ancestry(
            orion_path, root=authorized_pic_root
        ) as orion_ancestry, _prepare_completion_receipt_ancestry(
            project_home_path, root=authorized_project_home_root
        ) as project_home_ancestry:
            orion_parent = os.fstat(orion_ancestry.descriptor)
            project_home_parent = os.fstat(project_home_ancestry.descriptor)
            receipt = {
                "schema_version": 1,
                "record_type": TRAMPOLINE_COMPLETION_RECORD_TYPE,
                "receipt_role": "paired_immutable_pre_reconciliation_execution_anchor",
                "authority": {
                    "launch_authorized": False,
                    "scientific_claim_authorized": False,
                    "publication_authorized": False,
                },
                "paired_paths": {
                    "orion": str(orion_path),
                    "project_home": str(project_home_path),
                },
                "receipt_parent_identities": {
                    "orion": {
                        "device": orion_parent.st_dev,
                        "inode": orion_parent.st_ino,
                    },
                    "project_home": {
                        "device": project_home_parent.st_dev,
                        "inode": project_home_parent.st_ino,
                    },
                },
                "manifest_path": str(manifest_path),
                "manifest_sha256": manifest_sha256,
                "reservation_id": reservation_id,
                "submission_id": submission_id,
                "slurm_job_id": slurm_job_id,
                "control_plane_version": manifest["control_plane_version"],
                "campaign": manifest["campaign"],
                "test_id": manifest["test_id"],
                "registered_science_authorization_id": manifest.get(
                    "registered_science_authorization_id"
                ),
                "artifact_dir": str(artifact_dir),
                "artifact_root_identity": {
                    "device": root_metadata.st_dev,
                    "inode": root_metadata.st_ino,
                },
                "artifact_inventory": {
                    "path": str(artifact_dir / "artifact_inventory.json"),
                    "sha256": hashlib.sha256(inventory_payload).hexdigest(),
                    "byte_count": len(inventory_payload),
                    "filesystem_identity": inventory_identity,
                    "payload_base64": base64.b64encode(inventory_payload).decode("ascii"),
                },
                "artifact_records": completion_records,
                "mandatory_stdout_stderr": [
                    indexed[path] for path in mandatory_stdio
                ],
                "execution_binding": {
                    "launch_contract_sha256": launch_contract_sha256(contract),
                    "job_script_sha256": record_for_role(manifest, "job-script")[
                        "sha256"
                    ],
                    "executable_sha256": record_for_role(manifest, "executable")[
                        "sha256"
                    ],
                    "input_deck_sha256": record_for_role(manifest, "input-deck")[
                        "sha256"
                    ],
                },
            }
            payload = (
                json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
            ).encode("utf-8")
            receipt_sha256 = hashlib.sha256(payload).hexdigest()
            ledger_binding = {
                "schema_version": 1,
                "record_type": TRAMPOLINE_COMPLETION_LEDGER_RECORD_TYPE,
                "authority": receipt["authority"],
                "receipt_sha256": receipt_sha256,
                "receipt_byte_count": len(payload),
                "paired_receipts": {
                    "orion": {
                        "path": str(orion_path),
                        "parent_identity": receipt["receipt_parent_identities"]["orion"],
                    },
                    "project_home": {
                        "path": str(project_home_path),
                        "parent_identity": receipt["receipt_parent_identities"][
                            "project_home"
                        ],
                    },
                },
                "artifact_root_identity": receipt["artifact_root_identity"],
                "artifact_inventory": {
                    key: receipt["artifact_inventory"][key]
                    for key in ("sha256", "byte_count", "filesystem_identity")
                },
                "artifact_records_sha256": _canonical_sha256(completion_records),
                "mandatory_stdout_stderr_sha256": _canonical_sha256(
                    receipt["mandatory_stdout_stderr"]
                ),
            }
            anchored_event = append_anchor(ledger_binding)
            if anchored_event.get("trampoline_completion") != ledger_binding:
                raise ValueError(
                    "Trampoline completion differs from pre-publication ledger anchor"
                )
            orion_identity = _publish_exact_completion_receipt_at(
                orion_path, payload, ancestry=orion_ancestry
            )
            project_home_identity = _publish_exact_completion_receipt_at(
                project_home_path, payload, ancestry=project_home_ancestry
            )
            if orion_identity == project_home_identity:
                raise ValueError(
                    "Trampoline completion receipts reuse one filesystem object"
                )
            if (
                _publish_exact_completion_receipt_at(
                    orion_path, payload, ancestry=orion_ancestry
                )
                != orion_identity
                or _publish_exact_completion_receipt_at(
                    project_home_path, payload, ancestry=project_home_ancestry
                )
                != project_home_identity
            ):
                raise ValueError("Trampoline completion receipt identity changed")
        require_same_directory(artifact_dir, artifact_dir_fd, root=trusted_anchor)
        if _read_retained_completion_artifact(
            inventory_parent_fd,
            inventory_fd,
            inventory_name,
            relative="artifact_inventory.json",
            expected_sha256=receipt["artifact_inventory"]["sha256"],
            expected_size=receipt["artifact_inventory"]["byte_count"],
        )[0] != inventory_payload:
            raise ValueError("Trampoline completion inventory changed after publication")
        for parent_fd, descriptor, name, relative, raw in retained:
            _read_retained_completion_artifact(
                parent_fd,
                descriptor,
                name,
                relative=relative,
                expected_sha256=str(raw["sha256"]),
                expected_size=int(raw["size"]),
            )
        require_same_directory(artifact_dir, artifact_dir_fd, root=trusted_anchor)
        return {
            "orion_path": str(orion_path),
            "project_home_path": str(project_home_path),
            "sha256": receipt_sha256,
            "byte_count": len(payload),
            "ledger_binding": ledger_binding,
            "anchored_event": anchored_event,
        }
    finally:
        for parent_fd, descriptor, *_ in reversed(retained):
            os.close(descriptor)
            os.close(parent_fd)
        os.close(inventory_fd)
        os.close(inventory_parent_fd)


def _append_trampoline_completion_event(
    reservation: dict[str, object],
    ledger_binding: dict[str, object],
    *,
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    authorized_pic_root: Path,
) -> dict[str, object]:
    """Append or recover the exact completion event before receipt publication."""
    ledger_csv = Path(os.path.abspath(authorized_pic_root)) / "ledger" / "node_hours.csv"
    completion_event = transition_payload(reservation)
    completion_event.update(
        {
            "event_type": "trampoline_completion",
            "trampoline_completion": ledger_binding,
        }
    )
    reservation_id = str(reservation["reservation_id"])
    for attempt in range(2):
        try:
            with ledger_lock(ledger_jsonl, mirror_jsonl):
                try:
                    records = validate_mirrored_state(
                        ledger_jsonl, receipts_jsonl, mirror_jsonl
                    )
                except ValueError:
                    repair_mirrored_state_locked(
                        ledger_jsonl,
                        ledger_csv,
                        receipts_jsonl,
                        mirror_jsonl,
                        mirror_transport="filesystem_copy",
                    )
                    records = validate_mirrored_state(
                        ledger_jsonl, receipts_jsonl, mirror_jsonl
                    )
                latest = latest_reservations(records).get(reservation_id)
                if latest is None:
                    raise ValueError(
                        "Trampoline completion reservation disappeared before append"
                    )
                if latest.get("event_type") == "trampoline_completion":
                    if transition_payload(latest) != transition_payload(completion_event):
                        raise ValueError(
                            "Existing trampoline completion differs from retry"
                        )
                    return latest
                if transition_payload(latest) != transition_payload(reservation):
                    raise ValueError(
                        "Trampoline completion reservation changed before append"
                    )
                return append_primary_event_locked(
                    ledger_jsonl,
                    ledger_csv,
                    receipts_jsonl,
                    mirror_jsonl,
                    completion_event,
                    mirror_transport="filesystem_copy",
                )
        except Exception:
            if attempt:
                raise
    raise AssertionError("unreachable trampoline completion append retry")


def recover_published_trampoline_completion(
    *,
    manifest_path: Path,
    manifest_sha256: str,
    job_script_sha256: str,
    executable_sha256: str,
    reservation_id: str,
    submission_id: str,
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    control_plane_dir: Path = Path(__file__).absolute().parent,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    """Recover a published completion anchor without rerunning a zero-retry job."""
    manifest, reservation = reservation_bound_manifest(
        manifest_path,
        reservation_id,
        ledger_jsonl=ledger_jsonl,
        receipts_jsonl=receipts_jsonl,
        mirror_jsonl=mirror_jsonl,
        control_plane_dir=control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    if (
        reservation.get("event_type") not in {"trampoline_completion", "reconciliation"}
        or not isinstance(reservation.get("trampoline_completion"), dict)
    ):
        raise ValueError(
            "Completion recovery requires an existing pre-publication ledger anchor"
        )
    job_id = str(reservation.get("job_id", ""))
    if (
        re.fullmatch(r"[0-9]+", job_id) is None
        or reservation.get("manifest_sha256") != manifest_sha256
        or reservation.get("submission_id") != submission_id
        or record_for_role(manifest, "job-script").get("sha256") != job_script_sha256
        or record_for_role(manifest, "executable").get("sha256") != executable_sha256
    ):
        raise ValueError("Completion recovery arguments differ from the submitted job")
    artifact_dir = _require_run_artifact_dir(manifest)
    trusted_anchor = stable_serialization_anchor(authorized_pic_root)

    def require_existing_anchor(binding: dict[str, object]) -> dict[str, object]:
        if reservation.get("trampoline_completion") != binding:
            raise ValueError("Existing trampoline completion differs from retry")
        return reservation

    with PinnedDirectoryAncestry(artifact_dir, root=trusted_anchor) as ancestry:
        completion = _publish_trampoline_completion_receipt_at(
            ancestry.descriptor,
            artifact_dir,
            manifest,
            append_anchor=require_existing_anchor,
            manifest_path=manifest_path,
            manifest_sha256=manifest_sha256,
            reservation_id=reservation_id,
            submission_id=submission_id,
            slurm_job_id=job_id,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        appended = completion["anchored_event"]
        if appended.get("trampoline_completion") != completion["ledger_binding"]:
            raise ValueError("Recovered trampoline completion differs from ledger anchor")
        ancestry.require_same()
        return appended


def _profile_environment() -> dict[str, str]:
    return {
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
        "PIC_FRONTIER_PROFILE": "frontier_minimum_supported",
    }


def _registered_trusted_wrapper_evidence_bytes(
    manifest: dict[str, object],
    action: dict[str, object],
    stdout_payload: bytes,
) -> bytes:
    """Return aggregate registered-execution evidence after rank verification."""
    specification = REGISTERED_WRAPPER_CAMPAIGNS.get(manifest.get("campaign"))
    if specification is None:
        return b""
    label, case_pattern = specification
    case_id = manifest.get("test_id")
    resources = action.get("resources")
    if (
        not isinstance(case_id, str)
        or not re.fullmatch(case_pattern, case_id)
        or not isinstance(resources, dict)
        or type(resources.get("tasks")) is not int
        or int(resources["tasks"]) <= 0
    ):
        raise ValueError(f"{label} trusted-wrapper manifest identity is malformed")
    try:
        lines = stdout_payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} trusted-wrapper stdout is not UTF-8") from error
    pattern = re.compile(
        r"^PIC trusted GPU launch: rank=([0-9]+) host=\S+ "
        r"ROCR_VISIBLE_DEVICES=[0-9]+ "
        r"linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa$"
    )
    observed = [
        int(match.group(1))
        for line in lines
        if (match := pattern.fullmatch(line)) is not None
    ]
    expected = list(range(int(resources["tasks"])))
    if sorted(observed) != expected or len(observed) != len(expected):
        raise ValueError(
            f"{label} trusted-wrapper task-rank evidence is incomplete or duplicated"
        )
    return (
        f"{label}_REGISTERED_EXECUTION case_id={case_id} "
        f"mpi_world_size={len(expected)} "
        f"rank_ids={','.join(str(rank) for rank in expected)}\n"
        f"{label}_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0\n"
    ).encode("utf-8")


def _q043_trusted_wrapper_evidence_bytes(
    manifest: dict[str, object],
    action: dict[str, object],
    stdout_payload: bytes,
) -> bytes:
    """Backward-compatible Q043 helper retained for focused tests."""
    return _registered_trusted_wrapper_evidence_bytes(
        manifest, action, stdout_payload
    )


def _create_artifact_directory(
    artifact_dir: Path, *, pic_root: Path
) -> tuple[int, PinnedDirectoryAncestry]:
    trusted_anchor = stable_serialization_anchor(pic_root)
    durable_mkdir_parents(artifact_dir.parent, root=trusted_anchor)
    parent_ancestry = PinnedDirectoryAncestry(
        artifact_dir.parent, root=trusted_anchor
    )
    parent_fd = parent_ancestry.descriptor
    child_fd: int | None = None
    try:
        child_fd = _open_created_directory_at(
            parent_fd,
            artifact_dir.name,
            mode=0o755,
            label="Launch artifact directory",
        )
        os.fsync(child_fd)
        os.fsync(parent_fd)
        parent_ancestry.require_same()
        require_same_directory(artifact_dir, child_fd, root=trusted_anchor)
        return child_fd, parent_ancestry
    except BaseException:
        if child_fd is not None:
            os.close(child_fd)
        parent_ancestry.close()
        raise


_MEMFD_SEALS = (
    fcntl.F_SEAL_SEAL
    | fcntl.F_SEAL_SHRINK
    | fcntl.F_SEAL_GROW
    | fcntl.F_SEAL_WRITE
)


def _captured_control_plane_member_at(
    directory_descriptor: int,
    name: str,
    manifest: dict[str, object],
    *,
    executable: bool,
) -> int:
    """Copy one manifest-bound installed member into a sealed anonymous file."""
    records = manifest.get("control_plane_inventory")
    if not isinstance(records, list):
        raise ValueError("Manifest control-plane inventory is malformed")
    expected = [
        record
        for record in records if isinstance(record, dict)
        and record.get("path") == name
    ]
    if (
        len(expected) != 1
        or set(expected[0]) != {"path", "sha256"}
        or re.fullmatch(r"[0-9a-f]{64}", str(expected[0]["sha256"])) is None
    ):
        raise ValueError(f"Control-plane inventory lacks one exact member: {name}")
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_descriptor,
    )
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_mode & 0o222
            or executable and not before.st_mode & 0o111
        ):
            raise ValueError(
                f"Control-plane captured member is not an eligible read-only file: {name}"
            )
        chunks = []
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
        if (
            (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
            != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
            or hashlib.sha256(payload).hexdigest() != expected[0]["sha256"]
        ):
            raise ValueError(f"Control-plane captured member changed or drifted: {name}")
    finally:
        os.close(descriptor)
    captured = os.memfd_create(
        f"athenak-pic-{name}", flags=os.MFD_CLOEXEC | os.MFD_ALLOW_SEALING
    )
    try:
        view = memoryview(payload)
        while view:
            written = os.write(captured, view)
            if written <= 0:
                raise OSError("short write while capturing control-plane member")
            view = view[written:]
        os.fchmod(captured, 0o555 if executable else 0o444)
        os.lseek(captured, 0, os.SEEK_SET)
        fcntl.fcntl(captured, fcntl.F_ADD_SEALS, _MEMFD_SEALS)
    except BaseException:
        os.close(captured)
        raise
    return captured


def _launch_actions(
    manifest: dict[str, object],
    *,
    reservation: dict[str, object],
    executable: _PinnedSnapshot,
    input_deck: _PinnedSnapshot,
    athena_timeout_arguments: tuple[str, str],
    profile_launcher: str,
    profile_launcher_fd: int,
    profile_environment_fd: int,
    slurm_job_id: str,
    runner: Callable[..., object],
    manifest_path: Path,
    manifest_sha256: str,
    reservation_id: str,
    submission_id: str,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
) -> None:
    raw_contract = manifest.get("launch_contract")
    _require_no_action_timeout_override(raw_contract)
    contract = validate_launch_contract(raw_contract)
    pic_root = Path(str(manifest["pic_root"]))
    trusted_anchor = stable_serialization_anchor(pic_root)
    artifact_dir = _require_run_artifact_dir(manifest)
    artifact_dir_fd, artifact_parent_ancestry = _create_artifact_directory(
        artifact_dir, pic_root=pic_root
    )
    launch_directory_identities: dict[str, tuple[int, int]] = {}

    def require_artifact_directory() -> None:
        require_same_directory(artifact_dir, artifact_dir_fd, root=trusted_anchor)
        artifact_parent_ancestry.require_same()

    try:
        require_artifact_directory()
        _bounded_actions(contract["pre_actions"], manifest, artifact_dir, artifact_dir_fd)
        for action in contract["actions"]:
            require_artifact_directory()
            executable.require_lexical_parent()
            input_deck.require_lexical_parent()
            resources = action["resources"]
            command = [
                profile_launcher,
                SRUN,
                f"--jobid={slurm_job_id}",
                f"-N{resources['nodes']}",
                f"-n{resources['tasks']}",
                f"-c{resources['cpus_per_task']}",
                f"--gpus-per-task={resources['gpus_per_task']}",
                f"--gpu-bind={resources['gpu_bind']}",
                TRUSTED_PYTHON,
                "-I",
                "-c",
                _TASK_LOCAL_EXEC,
                str(pic_root),
                str(executable.path),
                executable.expected_sha256,
                str(input_deck.path),
                input_deck.expected_sha256,
                *athena_timeout_arguments,
            ]
            for argument in action["arguments"]:
                if "literal" in argument:
                    command.append(str(argument["literal"]))
                elif argument.get("snapshot_role") == "input-deck":
                    command.append(_INPUT_DECK_FD_TOKEN)
                else:
                    directory = _artifact_path(artifact_dir, argument["artifact_directory"])
                    _mkdir_artifact_directory(
                        artifact_dir_fd,
                        artifact_dir,
                        directory,
                        launch_directory_identities,
                    )
                    command.append(str(directory))
            stdout_path = _artifact_path(artifact_dir, action["stdout_artifact"])
            stderr_path = _artifact_path(artifact_dir, action["stderr_artifact"])
            allowlist_path = _artifact_path(
                artifact_dir, f"{action['action_id']}.environment.allowlist.txt"
            )
            environment = _profile_environment()
            allowlist_fd: int | None = None
            try:
                _require_retained_artifact_directories_at(
                    artifact_dir_fd, artifact_dir, launch_directory_identities
                )
                allowlist_fd = _open_new_artifact(
                    artifact_dir_fd, artifact_dir, allowlist_path
                )
                environment["PIC_RUNTIME_ALLOWLIST_FD"] = str(allowlist_fd)
                environment["PIC_RUNTIME_ALLOWLIST_DIR_FD"] = str(artifact_dir_fd)
                environment["PIC_FRONTIER_PROFILE_FD"] = str(profile_environment_fd)
                with os.fdopen(
                    _open_new_artifact(artifact_dir_fd, artifact_dir, stdout_path), "wb"
                ) as stdout, os.fdopen(
                    _open_new_artifact(artifact_dir_fd, artifact_dir, stderr_path), "wb"
                ) as stderr:
                    try:
                        runner(
                            command,
                            check=True,
                            stdout=stdout,
                            stderr=stderr,
                            env=environment,
                            pass_fds=(
                                allowlist_fd,
                                artifact_dir_fd,
                                profile_launcher_fd,
                                profile_environment_fd,
                            ),
                        )
                        stdout.flush()
                        os.fsync(stdout.fileno())
                        wrapper_evidence = _registered_trusted_wrapper_evidence_bytes(
                            manifest,
                            action,
                            _read_artifact_bytes(
                                artifact_dir_fd, artifact_dir, stdout_path
                            ),
                        )
                        if wrapper_evidence:
                            stdout.write(wrapper_evidence)
                            stdout.flush()
                            os.fsync(stdout.fileno())
                    finally:
                        require_artifact_directory()
                        _capture_artifact_directory_identities_at(
                            artifact_dir_fd, launch_directory_identities
                        )
                        _require_retained_artifact_directories_at(
                            artifact_dir_fd, artifact_dir, launch_directory_identities
                        )
                executable.require_lexical_parent()
                input_deck.require_lexical_parent()
            finally:
                if allowlist_fd is not None:
                    os.close(allowlist_fd)
        _bounded_actions(contract["post_actions"], manifest, artifact_dir, artifact_dir_fd)
        require_artifact_directory()
        _require_retained_artifact_directories_at(
            artifact_dir_fd, artifact_dir, launch_directory_identities
        )
        _capture_artifact_directory_identities_at(
            artifact_dir_fd, launch_directory_identities
        )
        _publish_frozen_artifact_inventory(
            artifact_dir_fd,
            artifact_dir,
            directory_identities=launch_directory_identities,
        )
        require_artifact_directory()
        completion = _publish_trampoline_completion_receipt_at(
            artifact_dir_fd,
            artifact_dir,
            manifest,
            append_anchor=lambda binding: _append_trampoline_completion_event(
                reservation,
                binding,
                ledger_jsonl=ledger_jsonl,
                receipts_jsonl=receipts_jsonl,
                mirror_jsonl=mirror_jsonl,
                authorized_pic_root=authorized_pic_root,
            ),
            manifest_path=manifest_path,
            manifest_sha256=manifest_sha256,
            reservation_id=reservation_id,
            submission_id=submission_id,
            slurm_job_id=slurm_job_id,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        appended = completion["anchored_event"]
        if appended.get("trampoline_completion") != completion["ledger_binding"]:
            raise ValueError("Trampoline completion ledger anchor differs after append")
        require_artifact_directory()
    finally:
        try:
            require_artifact_directory()
        finally:
            try:
                artifact_parent_ancestry.close()
            finally:
                os.close(artifact_dir_fd)


def _bounded_actions(
    actions: list[dict[str, object]],
    manifest: dict[str, object],
    artifact_dir: Path,
    artifact_dir_fd: int,
) -> None:
    for action in actions:
        kind = str(action["kind"])
        if kind == "snapshot_sha256":
            record = record_for_role(manifest, str(action["snapshot_role"]))
            source = Path(str(record["path"]))
            with _PinnedSnapshot(
                source,
                root=Path(str(manifest["pic_root"])),
                expected_sha256=str(record["sha256"]),
            ):
                digest = str(record["sha256"])
            output = _artifact_path(artifact_dir, action["output_artifact"])
            _write_new_text_artifact(artifact_dir_fd, artifact_dir, output, digest + "\n")
        elif kind == "artifact_sha256":
            source = _artifact_path(artifact_dir, action["artifact"])
            try:
                source_bytes = _read_artifact_bytes(artifact_dir_fd, artifact_dir, source)
            except (FileNotFoundError, NotADirectoryError):
                raise ValueError(f"Missing artifact for checksum action: {source}")
            output = _artifact_path(artifact_dir, action["output_artifact"])
            _write_new_text_artifact(
                artifact_dir_fd,
                artifact_dir,
                output,
                hashlib.sha256(source_bytes).hexdigest() + "\n",
            )
        else:
            source = _artifact_path(artifact_dir, action["artifact"])
            try:
                size = _artifact_size(artifact_dir_fd, artifact_dir, source)
            except (FileNotFoundError, NotADirectoryError):
                size = 0
            if size <= 0:
                raise ValueError(f"Expected a non-empty launch artifact: {source}")


def launch(
    *,
    manifest_path: Path,
    manifest_sha256: str,
    job_script_sha256: str,
    executable_sha256: str,
    reservation_id: str,
    submission_id: str,
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    runner: Callable[..., object] = subprocess.run,
    control_plane_dir: Path = Path(__file__).absolute().parent,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> None:
    slurm_job_id = os.environ.get("SLURM_JOB_ID", "")
    if not re.fullmatch(r"[0-9]+", slurm_job_id):
        raise ValueError("Trampoline requires a live numeric SLURM_JOB_ID")
    trusted_anchor = stable_serialization_anchor(authorized_pic_root)
    require_ledger_paths(
        ledger_jsonl,
        authorized_pic_root.resolve() / "ledger" / "node_hours.csv",
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    manifest, reservation = executable_reservation_bound_manifest(
        manifest_path,
        reservation_id,
        ledger_jsonl=ledger_jsonl,
        receipts_jsonl=receipts_jsonl,
        mirror_jsonl=mirror_jsonl,
        control_plane_dir=control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        executable_job_id=slurm_job_id,
    )
    if reservation.get("manifest_sha256") != manifest_sha256:
        raise ValueError("Scheduled manifest checksum differs from reservation ledger")
    if reservation.get("submission_id") != submission_id:
        raise ValueError("Scheduled submission ID differs from reservation ledger")
    verify_snapshot_files(manifest, root=authorized_pic_root)
    for record in manifest["snapshot_files"]:
        with _PinnedSnapshot(
            Path(str(record["path"])),
            root=trusted_anchor,
            expected_sha256=str(record["sha256"]),
        ):
            pass

    job_script = record_for_role(manifest, "job-script")
    executable = record_for_role(manifest, "executable")
    input_deck = record_for_role(manifest, "input-deck")
    athena_timeout_arguments = _trusted_athena_timeout_arguments(
        manifest, root=trusted_anchor
    )
    if (
        job_script.get("sha256") != job_script_sha256
        or reservation.get("job_script_sha256") != job_script_sha256
    ):
        raise ValueError("Snapshotted job-script digest differs from scheduled binding")
    if (
        executable.get("sha256") != executable_sha256
        or reservation.get("executable_sha256") != executable_sha256
    ):
        raise ValueError("Snapshotted executable digest differs from scheduled binding")
    executable_path = os.path.abspath(str(executable["path"]))
    if manifest.get("job_script_executable_env") != "PIC_EXECUTABLE":
        raise ValueError("Manifest does not declare the PIC_EXECUTABLE launch contract")
    declared = os.environ.get("PIC_EXECUTABLE")
    if declared and os.path.abspath(declared) != executable_path:
        raise ValueError("PIC_EXECUTABLE differs from verified executable snapshot")
    os.environ["PIC_EXECUTABLE"] = executable_path
    control_plane_dir_ancestry = PinnedDirectoryAncestry(
        control_plane_dir, root=trusted_anchor
    )
    control_plane_dir_fd = control_plane_dir_ancestry.descriptor
    profile_launcher_fd: int | None = None
    profile_environment_fd: int | None = None
    try:
        require_same_directory(
            control_plane_dir, control_plane_dir_fd, root=trusted_anchor
        )
        control_plane_dir_ancestry.require_same()
        profile_launcher_fd = _captured_control_plane_member_at(
            control_plane_dir_fd,
            "launch_with_frontier_profile.sh",
            manifest,
            executable=True,
        )
        profile_environment_fd = _captured_control_plane_member_at(
            control_plane_dir_fd,
            "frontier_pic_environment.sh",
            manifest,
            executable=False,
        )
        profile_launcher = f"/proc/self/fd/{profile_launcher_fd}"
        with _PinnedSnapshot(
            Path(executable_path),
            root=trusted_anchor,
            expected_sha256=str(executable["sha256"]),
            require_executable=True,
        ) as pinned_executable, _PinnedSnapshot(
            Path(str(input_deck["path"])),
            root=trusted_anchor,
            expected_sha256=str(input_deck["sha256"]),
        ) as pinned_input_deck:
            with _deterministic_artifact_umask():
                _launch_actions(
                    manifest,
                    reservation=reservation,
                    executable=pinned_executable,
                    input_deck=pinned_input_deck,
                    athena_timeout_arguments=athena_timeout_arguments,
                    profile_launcher=profile_launcher,
                    profile_launcher_fd=profile_launcher_fd,
                    profile_environment_fd=profile_environment_fd,
                    slurm_job_id=slurm_job_id,
                    runner=runner,
                    manifest_path=manifest_path,
                    manifest_sha256=manifest_sha256,
                    reservation_id=reservation_id,
                    submission_id=submission_id,
                    authorized_pic_root=authorized_pic_root,
                    authorized_project_home_root=authorized_project_home_root,
                    ledger_jsonl=ledger_jsonl,
                    receipts_jsonl=receipts_jsonl,
                    mirror_jsonl=mirror_jsonl,
                )
        require_same_directory(
            control_plane_dir, control_plane_dir_fd, root=trusted_anchor
        )
        control_plane_dir_ancestry.require_same()
    finally:
        try:
            if profile_environment_fd is not None:
                os.close(profile_environment_fd)
            if profile_launcher_fd is not None:
                os.close(profile_launcher_fd)
        finally:
            control_plane_dir_ancestry.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recover-published-completion", action="store_true")
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--job-script-sha256", required=True)
    parser.add_argument("--executable-sha256", required=True)
    parser.add_argument("--reservation-id", required=True)
    parser.add_argument("--submission-id", required=True)
    parser.add_argument("--ledger-jsonl", required=True, type=Path)
    parser.add_argument("--receipts-jsonl", required=True, type=Path)
    parser.add_argument("--mirror-jsonl", required=True, type=Path)
    args = parser.parse_args()
    arguments = {
        "manifest_path": args.manifest,
        "manifest_sha256": args.manifest_sha256,
        "job_script_sha256": args.job_script_sha256,
        "executable_sha256": args.executable_sha256,
        "reservation_id": args.reservation_id,
        "submission_id": args.submission_id,
        "ledger_jsonl": args.ledger_jsonl,
        "receipts_jsonl": args.receipts_jsonl,
        "mirror_jsonl": args.mirror_jsonl,
    }
    if args.recover_published_completion:
        recovered = recover_published_trampoline_completion(**arguments)
        print(str(recovered["event_sha256"]))
    else:
        launch(**arguments)


if __name__ == "__main__":
    main()
