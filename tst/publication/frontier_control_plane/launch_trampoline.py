#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Re-verify one reserved immutable snapshot before executing its job contract."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
from contextlib import contextmanager
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
from control_plane_common import durable_mkdir_parents
from control_plane_common import PinnedDirectoryAncestry
from control_plane_common import record_for_role, require_ledger_paths
from control_plane_common import require_same_directory, validate_launch_contract
from control_plane_common import stable_serialization_anchor
from control_plane_common import verify_snapshot_files
from validate_and_reserve_frontier_job import _require_run_artifact_dir
from validate_and_reserve_frontier_job import executable_reservation_bound_manifest


SRUN = "/usr/bin/srun"
TRUSTED_PYTHON = "/opt/cray/pe/python/3.11.7/bin/python3"
_INPUT_DECK_FD_TOKEN = "__PIC_INPUT_DECK_FD__"
_DIRECTORY_OPEN_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_NEW_ARTIFACT_OPEN_FLAGS = (
    os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
)
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
            if not stat.S_ISREG(metadata.st_mode) or metadata.st_mode & 0o222:
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

    def __exit__(self, *_: object) -> None:
        if self.descriptor is not None:
            os.close(self.descriptor)
            self.descriptor = None
        if self.parent_ancestry is not None:
            self.parent_ancestry.close()
            self.parent_ancestry = None
            self.parent_descriptor = None


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
) -> dict[str, object]:
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_fd,
    )
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
        return {
            "path": relative,
            "sha256": hashlib.sha256(data).hexdigest(),
            "size": len(data),
        }
    finally:
        os.close(descriptor)


def _freeze_artifact_tree_at(
    directory_fd: int,
    prefix: tuple[str, ...] = (),
    directory_identities: dict[str, tuple[int, int]] | None = None,
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
) -> None:
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_fd,
    )
    try:
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
            or (entry.st_dev, entry.st_ino) != (after.st_dev, after.st_ino)
            or len(data) != after.st_size
            or expected.get("path") != relative
            or expected.get("size") != len(data)
            or expected.get("sha256") != hashlib.sha256(data).hexdigest()
        ):
            raise ValueError(f"Launch artifact changed after freezing: {relative}")
    finally:
        os.close(descriptor)


def _verify_frozen_artifact_tree_at(
    directory_fd: int,
    records: list[dict[str, object]],
    directory_identities: dict[str, tuple[int, int]],
    analysis_fd: int,
    prefix: tuple[str, ...] = (),
) -> None:
    expected = {str(record["path"]): record for record in records}
    expected_directories = dict(directory_identities)

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
                _verify_frozen_artifact_file_at(current_fd, name, relative, record)
            elif stat.S_ISDIR(metadata.st_mode):
                child_fd = os.open(name, _DIRECTORY_OPEN_FLAGS, dir_fd=current_fd)
                try:
                    child = os.fstat(child_fd)
                    expected_identity = expected_directories.pop(relative, None)
                    if expected_identity != (child.st_dev, child.st_ino):
                        raise ValueError(
                            f"Launch artifact directory changed after freezing: {relative}"
                        )
                    verify_tree(child_fd, (*current_prefix, name))
                    after = os.fstat(child_fd)
                    entry = os.stat(name, dir_fd=current_fd, follow_symlinks=False)
                    if (entry.st_dev, entry.st_ino) != (after.st_dev, after.st_ino):
                        raise ValueError(
                            f"Launch artifact directory changed after freezing: {relative}"
                        )
                finally:
                    os.close(child_fd)
            else:
                raise ValueError(f"Launch artifact tree contains an unsupported entry: {relative}")
        after = os.fstat(current_fd)
        if (
            (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino)
            or sorted(os.listdir(current_fd)) != names
        ):
            relative = "/".join(current_prefix) or "."
            raise ValueError(f"Launch artifact directory changed after freezing: {relative}")

    verify_tree(directory_fd, prefix)
    if expected:
        raise ValueError("Launch artifact tree lost a frozen file")
    if expected_directories:
        raise ValueError("Launch artifact tree lost a frozen directory")


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
        records = _freeze_artifact_tree_at(
            artifact_dir_fd,
            directory_identities=frozen_directory_identities,
        )
        _write_new_text_artifact(
            artifact_dir_fd,
            artifact_dir,
            artifact_dir / "artifact_inventory.json",
            json.dumps(
                {"schema_version": 1, "files": records},
                indent=2,
                sort_keys=True,
            )
            + "\n",
        )
        inventory_record = _freeze_artifact_file_at(
            artifact_dir_fd, "artifact_inventory.json", "artifact_inventory.json"
        )
        os.fchmod(artifact_dir_fd, 0o555)
        os.fsync(artifact_dir_fd)
        _verify_frozen_artifact_file_at(
            artifact_dir_fd,
            "artifact_inventory.json",
            "artifact_inventory.json",
            inventory_record,
        )
        _verify_frozen_artifact_tree_at(
            artifact_dir_fd,
            records,
            frozen_directory_identities,
            analysis_fd,
        )
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


def _profile_environment() -> dict[str, str]:
    return {
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
        "PIC_FRONTIER_PROFILE": "frontier_minimum_supported",
    }


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


def _require_read_only_regular_at(directory_descriptor: int, name: str) -> None:
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_descriptor,
    )
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_mode & 0o222:
            raise ValueError(f"Control-plane launcher is not a read-only regular file: {name}")
    finally:
        os.close(descriptor)


def _launch_actions(
    manifest: dict[str, object],
    *,
    executable: _PinnedSnapshot,
    input_deck: _PinnedSnapshot,
    profile_launcher: str,
    control_plane_dir_fd: int,
    slurm_job_id: str,
    runner: Callable[..., object],
) -> None:
    contract = validate_launch_contract(manifest.get("launch_contract"))
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
                environment["PIC_CONTROL_PLANE_DIR_FD"] = str(control_plane_dir_fd)
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
                                control_plane_dir_fd,
                            ),
                        )
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
    try:
        require_same_directory(
            control_plane_dir, control_plane_dir_fd, root=trusted_anchor
        )
        control_plane_dir_ancestry.require_same()
        _require_read_only_regular_at(
            control_plane_dir_fd, "launch_with_frontier_profile.sh"
        )
        profile_launcher = (
            f"/proc/self/fd/{control_plane_dir_fd}/launch_with_frontier_profile.sh"
        )
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
                    executable=pinned_executable,
                    input_deck=pinned_input_deck,
                    profile_launcher=profile_launcher,
                    control_plane_dir_fd=control_plane_dir_fd,
                    slurm_job_id=slurm_job_id,
                    runner=runner,
                )
        require_same_directory(
            control_plane_dir, control_plane_dir_fd, root=trusted_anchor
        )
        control_plane_dir_ancestry.require_same()
    finally:
        control_plane_dir_ancestry.close()


def main() -> None:
    parser = argparse.ArgumentParser()
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
    launch(
        manifest_path=args.manifest,
        manifest_sha256=args.manifest_sha256,
        job_script_sha256=args.job_script_sha256,
        executable_sha256=args.executable_sha256,
        reservation_id=args.reservation_id,
        submission_id=args.submission_id,
        ledger_jsonl=args.ledger_jsonl,
        receipts_jsonl=args.receipts_jsonl,
        mirror_jsonl=args.mirror_jsonl,
    )


if __name__ == "__main__":
    main()
