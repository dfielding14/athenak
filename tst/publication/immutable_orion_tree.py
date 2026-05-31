#!/usr/bin/env python3
"""Fail-closed helpers for recursively frozen Orion-backed evidence trees."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import stat
import tarfile
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, Iterator


INVENTORY_NAME = "artifact_inventory.sha256"
FREEZE_RECEIPT_NAME = "freeze_receipt.json"
INVENTORY_ALGORITHM = (
    "sha256 of '<file_sha256>  <root-relative-path>\\n' entries ordered "
    "lexically by root-relative path"
)
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH
_RESERVED_METADATA = {INVENTORY_NAME, FREEZE_RECEIPT_NAME}
_MEMFD_SEALS = (
    fcntl.F_SEAL_WRITE
    | fcntl.F_SEAL_GROW
    | fcntl.F_SEAL_SHRINK
    | fcntl.F_SEAL_SEAL
)


def _raise(error_type: type[ValueError], label: str, message: str) -> None:
    raise error_type(f"{label}: {message}")


def require_exact_primitive_types(
    actual: Any,
    expected: Any,
    *,
    error_type: type[ValueError] = ValueError,
    label: str = "immutable-tree retained metadata",
) -> None:
    """Reject JSON primitive aliases before a separate whole-object comparison."""
    if type(actual) is not type(expected):
        _raise(error_type, label, "retained primitive type drifted")
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            _raise(error_type, label, "retained object keys drifted")
        for name in sorted(expected):
            require_exact_primitive_types(
                actual[name],
                expected[name],
                error_type=error_type,
                label=f"{label} {name}",
            )
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            _raise(error_type, label, "retained list length drifted")
        for index, (actual_item, expected_item) in enumerate(zip(actual, expected)):
            require_exact_primitive_types(
                actual_item,
                expected_item,
                error_type=error_type,
                label=f"{label} item {index}",
            )


def _is_sealed_memfd_reference(path: Path) -> bool:
    """Return whether one procfs fd path names a fully sealed regular memfd."""
    if path.parent != Path("/proc/self/fd") or not path.name.isdigit():
        return False
    try:
        fd = os.open(path, os.O_RDONLY)
    except OSError:
        return False
    try:
        status = os.fstat(fd)
        return (
            stat.S_ISREG(status.st_mode)
            and status.st_nlink == 0
            and fcntl.fcntl(fd, fcntl.F_GET_SEALS) & _MEMFD_SEALS == _MEMFD_SEALS
        )
    except OSError:
        return False
    finally:
        os.close(fd)


def is_sealed_snapshot_member(path: str | Path) -> bool:
    """Expose the sealed-member check to retained-artifact analyzers."""
    return _is_sealed_memfd_reference(Path(path))


def _open_self_contained_regular(path: Path, *, error_type: type[ValueError], label: str) -> int:
    """Open a canonical retained file or one fully sealed snapshot memfd."""
    if path.parent == Path("/proc/self/fd") and path.name.isdigit():
        try:
            fd = os.open(path, os.O_RDONLY)
        except OSError as error:
            _raise(error_type, label, f"cannot open sealed snapshot member {path}: {error}")
        try:
            status = os.fstat(fd)
            if (
                stat.S_ISREG(status.st_mode)
                and status.st_nlink == 0
                and fcntl.fcntl(fd, fcntl.F_GET_SEALS) & _MEMFD_SEALS == _MEMFD_SEALS
            ):
                return fd
        except OSError as error:
            os.close(fd)
            _raise(error_type, label, f"cannot inspect sealed snapshot member {path}: {error}")
        os.close(fd)
        _raise(error_type, label, f"snapshot member is not a fully sealed regular memfd: {path}")
    try:
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    except OSError as error:
        _raise(error_type, label, f"cannot open retained regular file {path}: {error}")
    status = os.fstat(fd)
    if not stat.S_ISREG(status.st_mode) or status.st_nlink != 1:
        os.close(fd)
        _raise(error_type, label, f"retained regular file is unsafe: {path}")
    return fd


def _regular_identity(status: os.stat_result) -> tuple[int, ...]:
    """Return the fields that must remain stable across one ordinary-file read."""
    return (
        status.st_dev,
        status.st_ino,
        status.st_mode,
        status.st_nlink,
        status.st_size,
        status.st_mtime_ns,
        status.st_ctime_ns,
    )


def _read_open_regular_bytes(
    fd: int,
    path: Path,
    *,
    error_type: type[ValueError],
    label: str,
) -> tuple[bytes, int]:
    """Read one open regular file while rejecting mutation during the read."""
    os.lseek(fd, 0, os.SEEK_SET)
    before = os.fstat(fd)
    if not stat.S_ISREG(before.st_mode) or before.st_nlink not in (0, 1):
        _raise(error_type, label, f"retained regular file is unsafe: {path}")
    payload = bytearray()
    while chunk := os.read(fd, 1024 * 1024):
        payload.extend(chunk)
    after = os.fstat(fd)
    if _regular_identity(before) != _regular_identity(after):
        _raise(error_type, label, f"retained regular file changed while reading: {path}")
    os.lseek(fd, 0, os.SEEK_SET)
    return bytes(payload), before.st_mode


def _sha256_open_regular(
    fd: int,
    path: Path,
    *,
    error_type: type[ValueError],
    label: str,
) -> str:
    """Hash one already-open regular file and reject mutation during the read."""
    payload, _ = _read_open_regular_bytes(fd, path, error_type=error_type, label=label)
    return hashlib.sha256(payload).hexdigest()


def _sealed_memfd(payload: bytes, mode: int, *, name: str) -> int:
    """Materialize immutable bytes into one fully sealed anonymous regular file."""
    fd = os.memfd_create(name, flags=os.MFD_CLOEXEC | os.MFD_ALLOW_SEALING)
    try:
        view = memoryview(payload)
        while view:
            written = os.write(fd, view)
            if written <= 0:
                raise OSError("short write while materializing sealed regular bytes")
            view = view[written:]
        os.fchmod(fd, stat.S_IMODE(mode) & ~_WRITE_BITS)
        os.lseek(fd, 0, os.SEEK_SET)
        fcntl.fcntl(fd, fcntl.F_ADD_SEALS, _MEMFD_SEALS)
    except Exception:
        os.close(fd)
        raise
    return fd


def _sealed_regular_copy(
    path: Path,
    *,
    error_type: type[ValueError],
    label: str,
) -> tuple[int, str, int]:
    """Copy one mutation-checked regular file into a sealed descriptor."""
    source_fd = _open_self_contained_regular(path, error_type=error_type, label=label)
    try:
        payload, mode = _read_open_regular_bytes(
            source_fd,
            path,
            error_type=error_type,
            label=label,
        )
    finally:
        os.close(source_fd)
    return (
        _sealed_memfd(payload, mode, name=f"athenak-pic-{path.name}"),
        hashlib.sha256(payload).hexdigest(),
        mode,
    )


def _sha256_regular(
    path: Path,
    *,
    error_type: type[ValueError],
    label: str,
) -> str:
    """Hash one self-contained regular file without following its final link."""
    fd = _open_self_contained_regular(path, error_type=error_type, label=label)
    try:
        return _sha256_open_regular(fd, path, error_type=error_type, label=label)
    finally:
        os.close(fd)


def _read_regular_text(
    path: Path,
    *,
    error_type: type[ValueError],
    label: str,
) -> str:
    """Read one sealed regular UTF-8 snapshot member."""
    if not _is_sealed_memfd_reference(path):
        _raise(error_type, label, f"retained metadata file is not a sealed snapshot member: {path}")
    fd = _open_self_contained_regular(path, error_type=error_type, label=label)
    try:
        payload, _ = _read_open_regular_bytes(fd, path, error_type=error_type, label=label)
    finally:
        os.close(fd)
    try:
        return payload.decode("utf-8")
    except UnicodeDecodeError as error:
        _raise(error_type, label, f"retained metadata file is not UTF-8: {path}: {error}")


def _canonical_authorized_root(
    runtime_root: str | Path,
    *,
    authorized_root: Path,
    error_type: type[ValueError],
    label: str,
) -> Path:
    """Require a canonical real directory strictly below the authorized root."""
    candidate = Path(runtime_root)
    if not candidate.is_absolute():
        _raise(error_type, label, "retained tree root must be absolute")
    try:
        authorized = authorized_root.resolve(strict=True)
        resolved = candidate.resolve(strict=True)
        mode = os.lstat(candidate).st_mode
        resolved.relative_to(authorized)
    except (OSError, RuntimeError, ValueError) as error:
        _raise(error_type, label, f"retained tree must remain below {authorized_root}: {error}")
    if resolved != candidate:
        _raise(error_type, label, "retained tree root must be a canonical path without aliases")
    if resolved == authorized:
        _raise(error_type, label, "retained tree must not be the authorized bulk root")
    if not stat.S_ISDIR(mode):
        _raise(error_type, label, "retained tree root must be a real directory")
    return resolved


def authorized_tree_root(
    runtime_root: str | Path,
    *,
    authorized_root: Path,
    error_type: type[ValueError] = ValueError,
    label: str = "immutable-tree",
) -> Path:
    """Return one canonical retained-tree root strictly below the authorized root."""
    return _canonical_authorized_root(
        runtime_root,
        authorized_root=authorized_root,
        error_type=error_type,
        label=label,
    )


def validate_source_archive(
    archive_path: str | Path,
    expected_sha256: str,
    *,
    error_type: type[ValueError] = ValueError,
    label: str = "immutable-tree source archive",
) -> dict[str, int | str | bool]:
    """Independently reject unsafe or cache-bearing retained source archives."""
    path = Path(archive_path)
    fd, measured_sha256, _ = _sealed_regular_copy(path, error_type=error_type, label=label)
    try:
        if measured_sha256 != expected_sha256:
            _raise(error_type, label, "source archive SHA-256 drifted")
        names = set()
        regular_file_count = 0
        try:
            with os.fdopen(os.dup(fd), "rb") as handle:
                with tarfile.open(fileobj=handle, mode="r:*") as archive:
                    members = archive.getmembers()
                    for member in members:
                        name = member.name
                        relative = PurePosixPath(name)
                        if (
                            not name
                            or name.splitlines() != [name]
                            or relative.is_absolute()
                            or any(part in ("", ".", "..") for part in relative.parts)
                            or name != relative.as_posix()
                        ):
                            _raise(
                                error_type,
                                label,
                                f"source archive member name is unsafe: {name!r}",
                            )
                        if name in names:
                            _raise(
                                error_type,
                                label,
                                f"source archive member is duplicated: {name!r}",
                            )
                        names.add(name)
                        if "__pycache__" in relative.parts or relative.suffix in (".pyc", ".pyo"):
                            _raise(
                                error_type,
                                label,
                                f"source archive contains bytecode cache: {name!r}",
                            )
                        if member.isdir():
                            continue
                        if not member.isreg():
                            _raise(
                                error_type,
                                label,
                                f"source archive member is not regular: {name!r}",
                            )
                        regular_file_count += 1
        except (tarfile.TarError, OSError) as error:
            _raise(error_type, label, f"cannot decode retained source archive: {error}")
    finally:
        os.close(fd)
    return {
        "archive_sha256": measured_sha256,
        "member_count": len(names),
        "regular_file_count": regular_file_count,
        "passed": True,
    }


def validate_executable_elf(
    executable_path: str | Path,
    expected_sha256: str,
    *,
    error_type: type[ValueError] = ValueError,
    label: str = "immutable-tree executable",
) -> dict[str, str | bool]:
    """Require one retained regular executable with an ELF identity."""
    path = Path(executable_path)
    fd, measured_sha256, mode = _sealed_regular_copy(path, error_type=error_type, label=label)
    try:
        if measured_sha256 != expected_sha256:
            _raise(error_type, label, "executable SHA-256 drifted")
        if stat.S_IMODE(mode) & 0o111 != 0o111:
            _raise(error_type, label, "retained executable must retain all execute bits")
        if os.read(fd, 4) != b"\x7fELF":
            _raise(error_type, label, "retained executable is not an ELF binary")
    finally:
        os.close(fd)
    return {
        "executable_sha256": measured_sha256,
        "elf_identity": True,
        "executable_mode": oct(stat.S_IMODE(mode)),
    }


def validate_source_archive_dependencies(
    archive_path: str | Path,
    dependencies: Any,
    required_paths: set[str],
    *,
    error_type: type[ValueError] = ValueError,
    label: str = "immutable-tree source archive dependencies",
) -> dict[str, int | bool]:
    """Hash required build dependencies directly from one retained source archive."""
    if not isinstance(dependencies, dict) or set(dependencies) != required_paths:
        _raise(error_type, label, "retained dependency manifest path set drifted")
    for name, digest in dependencies.items():
        relative = PurePosixPath(name)
        if (
            not name
            or name != relative.as_posix()
            or relative.is_absolute()
            or any(part in ("", ".", "..") for part in relative.parts)
            or not isinstance(digest, str)
            or _SHA256_PATTERN.fullmatch(digest) is None
        ):
            _raise(error_type, label, f"retained dependency manifest entry is invalid: {name!r}")
    path = Path(archive_path)
    fd, _, _ = _sealed_regular_copy(path, error_type=error_type, label=label)
    try:
        try:
            with os.fdopen(os.dup(fd), "rb") as handle:
                with tarfile.open(fileobj=handle, mode="r:*") as archive:
                    members = {member.name: member for member in archive.getmembers()}
                    for name, expected_sha256 in dependencies.items():
                        member = members.get(name)
                        if member is None or not member.isreg():
                            _raise(
                                error_type,
                                label,
                                f"required archived dependency is absent: {name}",
                            )
                        extracted = archive.extractfile(member)
                        if extracted is None:
                            _raise(error_type, label, f"cannot read archived dependency: {name}")
                        measured_sha256 = hashlib.sha256(extracted.read()).hexdigest()
                        if measured_sha256 != expected_sha256:
                            _raise(
                                error_type,
                                label,
                                f"archived dependency SHA-256 drifted: {name}",
                            )
        except (tarfile.TarError, OSError) as error:
            _raise(error_type, label, f"cannot decode retained source archive: {error}")
    finally:
        os.close(fd)
    return {"dependency_count": len(dependencies), "passed": True}


def validate_serial_host_build_evidence(
    pinned_root: str | Path,
    *,
    file_resolver: Any | None = None,
    error_type: type[ValueError] = ValueError,
    label: str = "immutable-tree serial-host build evidence",
) -> dict[str, str | bool]:
    """Validate retained serial-host build controls and successful nonempty logs."""
    root = Path(pinned_root)

    def read(relative: str) -> str:
        path = root / relative if file_resolver is None else Path(file_resolver(relative))
        return _read_regular_text(path, error_type=error_type, label=label)

    try:
        profile = json.loads(read("build/build_profile.json"))
        compile_commands = json.loads(read("build/compile_commands.json"))
    except json.JSONDecodeError as error:
        _raise(error_type, label, f"retained build JSON is invalid: {error}")
    expected_profile = {
        "schema_version": 1,
        "profile": "bounded_serial_host_clean_build_provenance_only",
        "qualifying_evidence": False,
        "mpi": False,
        "gpu": False,
        "frontier_scheduler": False,
        "kronos": False,
    }
    require_exact_primitive_types(
        profile,
        expected_profile,
        error_type=error_type,
        label=f"{label} profile",
    )
    if profile != expected_profile:
        _raise(error_type, label, "retained serial-host build profile drifted")
    if not isinstance(compile_commands, list) or not compile_commands:
        _raise(error_type, label, "retained compile commands must be a nonempty list")
    cache = read("build/CMakeCache.txt")
    for expected in (
        "Athena_ENABLE_MPI:BOOL=OFF",
        "CMAKE_BUILD_TYPE:STRING=Release",
        "Kokkos_ENABLE_HIP:BOOL=OFF",
        "Kokkos_ENABLE_SERIAL:BOOL=ON",
    ):
        key = expected.partition("=")[0]
        matches = [line for line in cache.splitlines() if line.partition("=")[0] == key]
        if matches != [expected]:
            _raise(error_type, label, f"retained CMake cache binding drifted for {key!r}")
    config = read("build/config.hpp")
    for expected in ("#define MPI_PARALLEL_ENABLED 0", "#define OPENMP_PARALLEL_ENABLED 0"):
        key = expected.split(maxsplit=2)[1]
        matches = [
            line
            for line in config.splitlines()
            if line.split(maxsplit=2)[:2] == ["#define", key]
        ]
        if matches != [expected]:
            _raise(error_type, label, f"retained config.hpp binding drifted for {key!r}")

    def canonical_command(relative: str) -> list[str]:
        text = read(relative)
        try:
            argv = shlex.split(text)
        except ValueError as error:
            _raise(error_type, label, f"retained command sidecar is invalid: {relative}: {error}")
        if text != f"{shlex.join(argv)}\n":
            _raise(error_type, label, f"retained command sidecar is noncanonical: {relative}")
        return argv

    configure = canonical_command("build/configure_command.txt")
    refresh = canonical_command("build/configure_compile_commands_refresh_command.txt")
    build = canonical_command("build/build_command.txt")
    rebuild = canonical_command("build/verbose_clean_rebuild_command.txt")
    if (
        len(configure) != 6
        or configure[:2] != ["cmake", "-S"]
        or configure[3] != "-B"
        or configure[5] != "-DCMAKE_BUILD_TYPE=Release"
        or not Path(configure[2]).is_absolute()
        or not Path(configure[4]).is_absolute()
    ):
        _raise(error_type, label, "retained configure command drifted")
    if refresh != [*configure, "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"]:
        _raise(error_type, label, "retained compile-command refresh drifted")
    expected_build_prefix = ["cmake", "--build", configure[4], "--target", "athena"]
    if build != [*expected_build_prefix, "--", "-j4"]:
        _raise(error_type, label, "retained initial build command drifted")
    if rebuild != [*expected_build_prefix, "--verbose", "--clean-first", "--", "-j4"]:
        _raise(error_type, label, "retained verbose clean rebuild command drifted")
    for relative in (
        "build/configure.log",
        "build/build.log",
        "build/configure_compile_commands_refresh.log",
        "build/verbose_clean_rebuild.log",
    ):
        if not read(relative):
            _raise(error_type, label, f"retained successful log is empty: {relative}")
    return {
        "passed": True,
        "profile": profile["profile"],
        "build_directory": configure[4],
        "mpi": False,
        "gpu": False,
    }


def _relative_text(
    relative: Path,
    *,
    error_type: type[ValueError],
    label: str,
) -> str:
    text = relative.as_posix()
    if text.splitlines() != [text]:
        _raise(error_type, label, "retained inventory paths must not contain line separators")
    return text


@dataclass(frozen=True)
class _AnchoredEntry:
    relative: str
    entry_type: str
    device: int
    inode: int
    mode: int
    link_count: int
    size: int
    mtime_ns: int
    ctime_ns: int
    sha256: str | None = None
    captured_payload: bytes | None = None

    def stable_signature(self) -> tuple[Any, ...]:
        return (
            self.entry_type,
            self.device,
            self.inode,
            self.mode,
            self.link_count,
            self.size,
            self.mtime_ns,
            self.ctime_ns,
        )

    def freeze_boundary_signature(self) -> tuple[Any, ...]:
        return (
            self.entry_type,
            self.device,
            self.inode,
            self.link_count,
            self.size,
            self.mtime_ns,
        )


@dataclass(frozen=True)
class _AnchoredTreeSnapshot:
    root: _AnchoredEntry
    entries: tuple[_AnchoredEntry, ...]

    def by_relative_path(self) -> dict[str, _AnchoredEntry]:
        return {entry.relative: entry for entry in self.entries}


class VerifiedFrozenTree:
    """Expose verified regular payloads only through lazily materialized sealed memfds."""

    def __init__(
        self,
        staged_root: Path,
        *,
        regular_sha256: dict[str, str],
        regular_sizes: dict[str, int],
        directories: frozenset[str],
        loader: Any,
    ) -> None:
        self.staged_root = staged_root
        self._regular_sha256 = regular_sha256
        self._regular_sizes = regular_sizes
        self._directories = directories
        self._loader = loader
        self._fds: dict[str, int] = {}

    @staticmethod
    def _relative_text(relative: str | Path) -> str:
        candidate = PurePosixPath(relative)
        text = candidate.as_posix()
        if (
            not text
            or text == "."
            or candidate.is_absolute()
            or any(part in ("", ".", "..") for part in candidate.parts)
        ):
            raise ValueError(f"unsafe sealed snapshot member path: {relative!r}")
        return text

    def has_file(self, relative: str | Path) -> bool:
        try:
            text = self._relative_text(relative)
        except ValueError:
            return False
        return text in self._regular_sha256

    def has_directory(self, relative: str | Path = ".") -> bool:
        """Return whether one topology-only snapshot directory exists."""
        candidate = PurePosixPath(relative)
        if candidate.as_posix() == ".":
            return True
        try:
            text = self._relative_text(relative)
        except ValueError:
            return False
        return text in self._directories

    def member_path(self, relative: str | Path) -> Path:
        """Return an immutable procfs handle for one verified regular member."""
        text = self._relative_text(relative)
        expected_sha256 = self._regular_sha256.get(text)
        if expected_sha256 is None:
            raise ValueError(f"sealed snapshot regular member is absent: {text}")
        if text not in self._fds:
            payload, mode = self._loader(text, expected_sha256)
            self._fds[text] = _sealed_memfd(
                payload,
                mode,
                name=f"athenak-pic-{Path(text).name}",
            )
        return Path("/proc/self/fd") / str(self._fds[text])

    def _relative_for_io(self, logical_root: Path, path: Path) -> Path | None:
        """Return the captured-tree label for one logical or staged path."""
        try:
            return path.relative_to(logical_root)
        except ValueError:
            try:
                return path.relative_to(self.staged_root)
            except ValueError:
                return None

    def file_io_path(self, logical_root: Path, path: Path) -> Path:
        """Route one captured regular member to an immutable procfs handle."""
        relative = self._relative_for_io(logical_root, path)
        if relative is None:
            return path
        if relative == Path("."):
            raise ValueError("sealed snapshot regular member is absent: .")
        self._relative_text(relative)
        return self.member_path(relative)

    def directory_io_path(self, logical_root: Path, path: Path) -> Path:
        """Route one captured topology-only directory to its staging path."""
        relative = self._relative_for_io(logical_root, path)
        if relative is None:
            return path
        if relative == Path("."):
            return self.staged_root
        self._relative_text(relative)
        if self.has_directory(relative):
            return self.staged_root / relative
        raise ValueError(f"sealed snapshot directory is absent: {relative.as_posix()}")

    def logical_path_for_io(self, logical_root: Path, path: str | Path) -> Path | None:
        """Map a staged-directory or sealed-fd label back to its retained logical label."""
        candidate = Path(path)
        try:
            relative = candidate.relative_to(self.staged_root)
        except ValueError:
            pass
        else:
            if relative == Path("."):
                return logical_root
            try:
                self._relative_text(relative)
            except ValueError:
                return None
            return logical_root / relative
        for relative, fd in self._fds.items():
            if candidate == Path("/proc/self/fd") / str(fd):
                return logical_root / relative
        return None

    def relative_files(self, subtree: str | Path = ".") -> set[str]:
        """Return regular-file labels below one logical subtree, relative to that subtree."""
        prefix = PurePosixPath(subtree)
        if prefix.as_posix() != ".":
            self._relative_text(subtree)
        return {
            PurePosixPath(relative).relative_to(prefix).as_posix()
            for relative in self._regular_sha256
            if PurePosixPath(relative).is_relative_to(prefix)
            and PurePosixPath(relative) != prefix
        }

    def relative_directories(self, subtree: str | Path = ".") -> set[str]:
        """Return directory labels below one logical subtree, relative to that subtree."""
        prefix = PurePosixPath(subtree)
        if prefix.as_posix() != ".":
            self._relative_text(subtree)
        return {
            PurePosixPath(relative).relative_to(prefix).as_posix()
            for relative in self._directories
            if PurePosixPath(relative).is_relative_to(prefix)
            and PurePosixPath(relative) != prefix
        }

    def immediate_directories(self, subtree: str | Path = ".") -> set[str]:
        """Return direct-child directory names below one logical subtree."""
        return {
            PurePosixPath(relative).parts[0]
            for relative in self.relative_directories(subtree)
            if PurePosixPath(relative).parts
        }

    def total_regular_size(self) -> int:
        """Return the verified aggregate size of regular snapshot members."""
        return sum(self._regular_sizes.values())

    def close(self) -> None:
        """Close every sealed member handle."""
        for fd in self._fds.values():
            os.close(fd)
        self._fds.clear()


def _entry_from_status(
    relative: str,
    entry_type: str,
    status: os.stat_result,
    *,
    sha256: str | None = None,
    captured_payload: bytes | None = None,
) -> _AnchoredEntry:
    return _AnchoredEntry(
        relative=relative,
        entry_type=entry_type,
        device=status.st_dev,
        inode=status.st_ino,
        mode=stat.S_IMODE(status.st_mode),
        link_count=status.st_nlink,
        size=status.st_size,
        mtime_ns=status.st_mtime_ns,
        ctime_ns=status.st_ctime_ns,
        sha256=sha256,
        captured_payload=captured_payload,
    )


def _require_same_opened_entry(
    observed: os.stat_result,
    opened: os.stat_result,
    relative: str,
    *,
    error_type: type[ValueError],
    label: str,
) -> None:
    """Reject a leaf or directory replacement between stat and openat."""
    if (observed.st_dev, observed.st_ino) != (opened.st_dev, opened.st_ino):
        _raise(error_type, label, f"retained entry changed during anchored traversal: {relative}")


def _hash_open_regular(
    fd: int,
    relative: str,
    *,
    capture_payload: bool,
    error_type: type[ValueError],
    label: str,
) -> tuple[str, bytes | None, os.stat_result]:
    """Hash one already-open regular file and reject mutation during the read."""
    before = os.fstat(fd)
    digest = hashlib.sha256()
    captured = bytearray() if capture_payload else None
    while chunk := os.read(fd, 1024 * 1024):
        digest.update(chunk)
        if captured is not None:
            captured.extend(chunk)
    after = os.fstat(fd)
    if _entry_from_status(relative, "file", before).stable_signature() != _entry_from_status(
        relative, "file", after
    ).stable_signature():
        _raise(error_type, label, f"retained regular file changed while hashing: {relative}")
    return digest.hexdigest(), bytes(captured) if captured is not None else None, after


def _scan_anchored_tree(
    root_fd: int,
    *,
    hash_regular: bool,
    capture_paths: frozenset[str] = frozenset(),
    error_type: type[ValueError],
    label: str,
) -> _AnchoredTreeSnapshot:
    """Walk one retained tree through open directory fds without reopening paths."""
    entries: list[_AnchoredEntry] = []
    root_before = os.fstat(root_fd)
    if not stat.S_ISDIR(root_before.st_mode):
        _raise(error_type, label, "retained tree root is no longer a directory")

    def walk(directory_fd: int, parent_parts: tuple[str, ...]) -> None:
        try:
            names = sorted(os.listdir(directory_fd))
        except OSError as error:
            _raise(error_type, label, f"cannot list retained tree through anchored fd: {error}")
        for name in names:
            relative = _relative_text(
                Path(*parent_parts, name), error_type=error_type, label=label
            )
            try:
                observed = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            except OSError as error:
                _raise(error_type, label, f"cannot inspect retained entry {relative}: {error}")
            if stat.S_ISLNK(observed.st_mode):
                _raise(error_type, label, f"retained tree contains symlink: {relative}")
            if stat.S_ISDIR(observed.st_mode):
                try:
                    child_fd = os.open(
                        name,
                        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                        dir_fd=directory_fd,
                    )
                except OSError as error:
                    _raise(
                        error_type,
                        label,
                        f"cannot open retained directory {relative} through anchored fd: {error}",
                    )
                try:
                    opened = os.fstat(child_fd)
                    _require_same_opened_entry(
                        observed, opened, relative, error_type=error_type, label=label
                    )
                    if not stat.S_ISDIR(opened.st_mode):
                        _raise(error_type, label, f"retained entry is not a directory: {relative}")
                    walk(child_fd, (*parent_parts, name))
                    after = os.fstat(child_fd)
                    if _entry_from_status(relative, "directory", opened).stable_signature() != (
                        _entry_from_status(relative, "directory", after).stable_signature()
                    ):
                        _raise(
                            error_type,
                            label,
                            f"retained directory changed during anchored traversal: {relative}",
                        )
                    entries.append(_entry_from_status(relative, "directory", after))
                finally:
                    os.close(child_fd)
                continue
            if not stat.S_ISREG(observed.st_mode):
                _raise(error_type, label, f"retained tree contains unsupported entry: {relative}")
            try:
                fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=directory_fd)
            except OSError as error:
                _raise(
                    error_type,
                    label,
                    f"cannot open retained regular file {relative} through anchored fd: {error}",
                )
            try:
                opened = os.fstat(fd)
                _require_same_opened_entry(
                    observed, opened, relative, error_type=error_type, label=label
                )
                if not stat.S_ISREG(opened.st_mode):
                    _raise(error_type, label, f"retained entry is not a regular file: {relative}")
                if opened.st_nlink != 1:
                    _raise(
                        error_type, label, f"retained tree contains hard-linked file: {relative}"
                    )
                digest = None
                captured_payload = None
                if hash_regular:
                    digest, captured_payload, opened = _hash_open_regular(
                        fd,
                        relative,
                        capture_payload=relative in capture_paths,
                        error_type=error_type,
                        label=label,
                    )
                entries.append(
                    _entry_from_status(
                        relative,
                        "file",
                        opened,
                        sha256=digest,
                        captured_payload=captured_payload,
                    )
                )
            finally:
                os.close(fd)

    walk(root_fd, ())
    root_after = os.fstat(root_fd)
    if _entry_from_status(".", "directory", root_before).stable_signature() != _entry_from_status(
        ".", "directory", root_after
    ).stable_signature():
        _raise(error_type, label, "retained tree root changed during anchored traversal")
    return _AnchoredTreeSnapshot(
        root=_entry_from_status(".", "directory", root_after),
        entries=tuple(sorted(entries, key=lambda entry: entry.relative)),
    )


def _require_same_snapshot(
    before: _AnchoredTreeSnapshot,
    after: _AnchoredTreeSnapshot,
    *,
    error_type: type[ValueError],
    label: str,
    phase: str,
) -> None:
    """Require a byte-hashing pass to retain the same anchored tree structure."""
    before_entries = {
        relative: entry.stable_signature()
        for relative, entry in before.by_relative_path().items()
    }
    after_entries = {
        relative: entry.stable_signature()
        for relative, entry in after.by_relative_path().items()
    }
    if (
        before.root.stable_signature() != after.root.stable_signature()
        or before_entries != after_entries
    ):
        _raise(error_type, label, f"retained tree changed during {phase}")


def _require_freeze_boundary(
    preflight: _AnchoredTreeSnapshot,
    boundary: _AnchoredTreeSnapshot,
    *,
    error_type: type[ValueError],
    label: str,
) -> None:
    """Reject staging replacements while allowing the freezer's chmod operations."""
    preflight_entries = {
        relative: entry.freeze_boundary_signature()
        for relative, entry in preflight.by_relative_path().items()
    }
    boundary_entries = {
        relative: entry.freeze_boundary_signature()
        for relative, entry in boundary.by_relative_path().items()
        if relative != FREEZE_RECEIPT_NAME
    }
    root_identity = lambda entry: (entry.entry_type, entry.device, entry.inode)
    if (
        root_identity(preflight.root) != root_identity(boundary.root)
        or preflight_entries != boundary_entries
    ):
        _raise(error_type, label, "retained staging tree changed before immutable boundary")


def _require_published_tree_matches_staging(
    staged: _AnchoredTreeSnapshot,
    published: _AnchoredTreeSnapshot,
    *,
    error_type: type[ValueError],
    label: str,
) -> None:
    """Reject replacements after hashing while allowing inventory publication."""
    staged_entries = {
        relative: entry.stable_signature()
        for relative, entry in staged.by_relative_path().items()
    }
    published_entries = {
        relative: entry.stable_signature()
        for relative, entry in published.by_relative_path().items()
        if relative != INVENTORY_NAME
    }
    root_identity = lambda entry: (entry.entry_type, entry.device, entry.inode)
    if (
        root_identity(staged.root) != root_identity(published.root)
        or staged_entries != published_entries
    ):
        _raise(error_type, label, "retained staging tree changed before inventory publication")


def _writable_entries(snapshot: _AnchoredTreeSnapshot, *, include_root: bool = True) -> list[str]:
    entries = [snapshot.root, *snapshot.entries] if include_root else list(snapshot.entries)
    return [entry.relative for entry in entries if entry.mode & _WRITE_BITS]


def _require_read_only(
    snapshot: _AnchoredTreeSnapshot,
    *,
    include_root: bool,
    error_type: type[ValueError],
    label: str,
) -> None:
    if _writable_entries(snapshot, include_root=include_root):
        _raise(error_type, label, "retained tree still contains writable entries")


def _require_reserved_metadata_absent(
    root_fd: int,
    *,
    error_type: type[ValueError],
    label: str,
) -> None:
    """Reject an already-populated reserved metadata slot before freezing."""
    for name in sorted(_RESERVED_METADATA):
        try:
            os.stat(name, dir_fd=root_fd, follow_symlinks=False)
        except FileNotFoundError:
            continue
        except OSError as error:
            _raise(error_type, label, f"cannot inspect reserved metadata file {name}: {error}")
        _raise(error_type, label, f"reserved metadata file already exists: {name}")


def _validate_receipt_payload(
    receipt: Any,
    *,
    error_type: type[ValueError],
    label: str,
) -> dict[str, Any]:
    """Require stable recursive freeze-policy semantics before any mutation."""
    if not isinstance(receipt, dict):
        _raise(error_type, label, "retained freeze receipt must be a JSON object")
    expected = {
        "schema_version": 1,
        "inventory_excludes": INVENTORY_NAME,
        "freeze_policy": "remove all owner, group and other write bits recursively",
    }
    if type(receipt.get("schema_version")) is not int:
        _raise(error_type, label, "retained freeze receipt schema_version drifted")
    for name, value in expected.items():
        if receipt.get(name) != value:
            _raise(error_type, label, f"retained freeze receipt {name} drifted")
    for name in ("artifact_role", "qualification_effect"):
        if not isinstance(receipt.get(name), str) or not receipt[name]:
            _raise(error_type, label, f"retained freeze receipt requires nonempty {name}")
    return receipt


def _decode_metadata(
    entry: _AnchoredEntry,
    name: str,
    *,
    error_type: type[ValueError],
    label: str,
) -> str:
    """Decode metadata bytes captured from the same fd used for hashing."""
    if entry.captured_payload is None:
        _raise(error_type, label, f"retained metadata file was not captured: {name}")
    try:
        return entry.captured_payload.decode("utf-8")
    except UnicodeDecodeError as error:
        _raise(error_type, label, f"retained metadata file is not UTF-8: {name}: {error}")


def _write_new_metadata(
    root_fd: int,
    name: str,
    payload: str,
    *,
    error_type: type[ValueError],
    label: str,
) -> None:
    """Create one reserved root metadata file exclusively without following links."""
    try:
        fd = os.open(
            name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
            0o444,
            dir_fd=root_fd,
        )
    except OSError as error:
        _raise(error_type, label, f"cannot create reserved metadata file {name}: {error}")
    try:
        encoded = payload.encode("utf-8")
        offset = 0
        while offset != len(encoded):
            offset += os.write(fd, encoded[offset:])
        os.fsync(fd)
    finally:
        os.close(fd)


def _remove_write_bits_below_root(
    root_fd: int,
    *,
    error_type: type[ValueError],
    label: str,
) -> None:
    """Freeze descendants through anchored fds while retaining root publication access."""

    def freeze(directory_fd: int, parent_parts: tuple[str, ...]) -> None:
        try:
            names = sorted(os.listdir(directory_fd))
        except OSError as error:
            _raise(error_type, label, f"cannot list retained tree through anchored fd: {error}")
        for name in names:
            relative = _relative_text(
                Path(*parent_parts, name), error_type=error_type, label=label
            )
            try:
                observed = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            except OSError as error:
                _raise(error_type, label, f"cannot inspect retained entry {relative}: {error}")
            if stat.S_ISLNK(observed.st_mode):
                _raise(error_type, label, f"retained tree contains symlink: {relative}")
            if stat.S_ISDIR(observed.st_mode):
                try:
                    child_fd = os.open(
                        name,
                        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                        dir_fd=directory_fd,
                    )
                except OSError as error:
                    _raise(
                        error_type,
                        label,
                        f"cannot open retained directory {relative} through anchored fd: {error}",
                    )
                try:
                    opened = os.fstat(child_fd)
                    _require_same_opened_entry(
                        observed, opened, relative, error_type=error_type, label=label
                    )
                    freeze(child_fd, (*parent_parts, name))
                    os.fchmod(child_fd, os.fstat(child_fd).st_mode & ~_WRITE_BITS)
                finally:
                    os.close(child_fd)
                continue
            if not stat.S_ISREG(observed.st_mode):
                _raise(error_type, label, f"retained tree contains unsupported entry: {relative}")
            try:
                fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=directory_fd)
            except OSError as error:
                _raise(
                    error_type,
                    label,
                    f"cannot open retained regular file {relative} through anchored fd: {error}",
                )
            try:
                opened = os.fstat(fd)
                _require_same_opened_entry(
                    observed, opened, relative, error_type=error_type, label=label
                )
                if not stat.S_ISREG(opened.st_mode):
                    _raise(error_type, label, f"retained entry is not a regular file: {relative}")
                if opened.st_nlink != 1:
                    _raise(
                        error_type, label, f"retained tree contains hard-linked file: {relative}"
                    )
                os.fchmod(fd, opened.st_mode & ~_WRITE_BITS)
            finally:
                os.close(fd)

    freeze(root_fd, ())


def _parse_inventory(
    text: str,
    *,
    error_type: type[ValueError],
    label: str,
) -> dict[str, str]:
    expected: dict[str, str] = {}
    ordered_paths = []
    for lineno, line in enumerate(text.splitlines(), 1):
        digest, separator, relative = line.partition("  ")
        if separator != "  " or not _SHA256_PATTERN.fullmatch(digest):
            _raise(error_type, label, f"malformed retained inventory line {lineno}")
        candidate = PurePosixPath(relative)
        if (
            not relative
            or candidate.is_absolute()
            or relative != candidate.as_posix()
            or any(part in ("", ".", "..") for part in candidate.parts)
            or relative == INVENTORY_NAME
        ):
            _raise(error_type, label, f"unsafe retained inventory path on line {lineno}")
        if relative in expected:
            _raise(error_type, label, f"duplicate retained inventory path: {relative}")
        expected[relative] = digest
        ordered_paths.append(relative)
    if ordered_paths != sorted(ordered_paths):
        _raise(error_type, label, "retained inventory paths are not sorted")
    canonical = "".join(f"{expected[path]}  {path}\n" for path in ordered_paths)
    if text != canonical:
        _raise(error_type, label, "retained inventory serialization is noncanonical")
    return expected


def _require_root_binding(
    root: Path,
    root_fd: int,
    *,
    authorized_root: Path,
    error_type: type[ValueError],
    label: str,
) -> None:
    """Require the public root pathname to retain the directory held by root_fd."""
    rebound = _canonical_authorized_root(
        root,
        authorized_root=authorized_root,
        error_type=error_type,
        label=label,
    )
    try:
        path_status = os.stat(rebound, follow_symlinks=False)
        opened_status = os.fstat(root_fd)
    except OSError as error:
        _raise(error_type, label, f"cannot revalidate retained tree root binding: {error}")
    if (
        not stat.S_ISDIR(path_status.st_mode)
        or (path_status.st_dev, path_status.st_ino)
        != (opened_status.st_dev, opened_status.st_ino)
    ):
        _raise(error_type, label, "retained tree root binding changed during anchored operation")


@contextmanager
def _open_anchored_root(
    root: Path,
    *,
    authorized_root: Path,
    error_type: type[ValueError],
    label: str,
) -> Iterator[int]:
    """Hold one root directory object for the complete freeze or verify operation."""
    try:
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    except OSError as error:
        _raise(error_type, label, f"cannot open retained tree root through anchored fd: {error}")
    try:
        _require_root_binding(
            root,
            root_fd,
            authorized_root=authorized_root,
            error_type=error_type,
            label=label,
        )
        yield root_fd
    finally:
        os.close(root_fd)


def _read_anchored_regular_bytes(
    root_fd: int,
    relative: str,
    *,
    expected_sha256: str,
    error_type: type[ValueError],
    label: str,
) -> tuple[bytes, int]:
    """Read one immutable member through root-relative descriptors and recheck its digest."""
    candidate = PurePosixPath(relative)
    if (
        not relative
        or candidate.is_absolute()
        or relative != candidate.as_posix()
        or any(part in ("", ".", "..") for part in candidate.parts)
    ):
        _raise(error_type, label, f"unsafe retained snapshot member path: {relative!r}")
    parent_fd = os.dup(root_fd)
    try:
        for part in candidate.parts[:-1]:
            child_fd = os.open(
                part,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=parent_fd,
            )
            os.close(parent_fd)
            parent_fd = child_fd
        fd = os.open(candidate.parts[-1], os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent_fd)
    except OSError as error:
        _raise(error_type, label, f"cannot open retained snapshot member {relative}: {error}")
    finally:
        os.close(parent_fd)
    try:
        before = os.fstat(fd)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_mode & _WRITE_BITS
        ):
            _raise(error_type, label, f"retained snapshot member is unsafe: {relative}")
        payload = bytearray()
        while chunk := os.read(fd, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(fd)
        identity = lambda value: (  # noqa: E731
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_nlink,
            value.st_size,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
        if identity(before) != identity(after):
            _raise(error_type, label, f"retained snapshot member changed while reading: {relative}")
        measured_sha256 = hashlib.sha256(payload).hexdigest()
        if measured_sha256 != expected_sha256:
            _raise(error_type, label, f"retained snapshot member SHA-256 drifted: {relative}")
        return bytes(payload), before.st_mode
    finally:
        os.close(fd)


@contextmanager
def staged_verified_frozen_tree(
    runtime_root: str | Path,
    expected_inventory_sha256: str,
    *,
    authorized_root: Path,
    error_type: type[ValueError] = ValueError,
    label: str = "immutable-tree",
) -> Iterator[tuple[dict[str, Any], VerifiedFrozenTree]]:
    """Expose descriptor-anchored verified bytes through sealed lazy snapshot members."""
    if not _SHA256_PATTERN.fullmatch(expected_inventory_sha256):
        _raise(error_type, label, "expected inventory SHA-256 must be 64 lowercase hex digits")
    root = _canonical_authorized_root(
        runtime_root,
        authorized_root=authorized_root,
        error_type=error_type,
        label=label,
    )
    with _open_anchored_root(
        root,
        authorized_root=authorized_root,
        error_type=error_type,
        label=label,
    ) as root_fd:
        report = _verify_frozen_tree_anchored(
            root,
            root_fd,
            expected_inventory_sha256,
            authorized_root=authorized_root,
            error_type=error_type,
            label=label,
        )
        inventory_payload, _ = _read_anchored_regular_bytes(
            root_fd,
            INVENTORY_NAME,
            expected_sha256=expected_inventory_sha256,
            error_type=error_type,
            label=label,
        )
        try:
            inventory_text = inventory_payload.decode("utf-8")
        except UnicodeDecodeError as error:
            _raise(error_type, label, f"retained inventory is not UTF-8: {error}")
        expected = _parse_inventory(inventory_text, error_type=error_type, label=label)
        with tempfile.TemporaryDirectory(prefix="athenak-pic-verified-tree-") as directory:
            staged_root = Path(directory) / "snapshot"
            staged_root.mkdir()
            snapshot = _scan_anchored_tree(
                root_fd,
                hash_regular=False,
                error_type=error_type,
                label=label,
            )
            for entry in snapshot.entries:
                if entry.entry_type == "directory" and entry.relative:
                    staged_root.joinpath(*PurePosixPath(entry.relative).parts).mkdir(
                        parents=True,
                        exist_ok=True,
                    )
            records = {INVENTORY_NAME: expected_inventory_sha256, **expected}
            _verify_frozen_tree_anchored(
                root,
                root_fd,
                expected_inventory_sha256,
                staged_snapshot=snapshot,
                authorized_root=authorized_root,
                error_type=error_type,
                label=label,
            )

            def load_member(relative: str, digest: str) -> tuple[bytes, int]:
                return _read_anchored_regular_bytes(
                    root_fd,
                    relative,
                    expected_sha256=digest,
                    error_type=error_type,
                    label=label,
                )

            sealed_snapshot = VerifiedFrozenTree(
                staged_root,
                regular_sha256=records,
                regular_sizes={
                    entry.relative: entry.size
                    for entry in snapshot.entries
                    if entry.entry_type == "file"
                },
                directories=frozenset(
                    entry.relative
                    for entry in snapshot.entries
                    if entry.entry_type == "directory" and entry.relative
                ),
                loader=load_member,
            )
            try:
                yield report, sealed_snapshot
            finally:
                sealed_snapshot.close()


@contextmanager
def staged_verified_legacy_read_only_tree(
    runtime_root: str | Path,
    expected_file_count: int,
    expected_inventory_sha256: str,
    *,
    authorized_root: Path,
    error_type: type[ValueError] = ValueError,
    label: str = "immutable-tree legacy read-only",
) -> Iterator[tuple[dict[str, Any], VerifiedFrozenTree]]:
    """Expose a legacy read-only tree through sealed members after anchored aggregate verification."""
    if expected_file_count < 1:
        _raise(error_type, label, "expected legacy file count must be positive")
    if not _SHA256_PATTERN.fullmatch(expected_inventory_sha256):
        _raise(error_type, label, "expected inventory SHA-256 must be 64 lowercase hex digits")
    root = _canonical_authorized_root(
        runtime_root,
        authorized_root=authorized_root,
        error_type=error_type,
        label=label,
    )
    with _open_anchored_root(
        root,
        authorized_root=authorized_root,
        error_type=error_type,
        label=label,
    ) as root_fd:
        baseline = _scan_anchored_tree(
            root_fd,
            hash_regular=False,
            error_type=error_type,
            label=label,
        )
        _require_read_only(baseline, include_root=True, error_type=error_type, label=label)
        measured = _scan_anchored_tree(
            root_fd,
            hash_regular=True,
            error_type=error_type,
            label=label,
        )
        _require_read_only(measured, include_root=True, error_type=error_type, label=label)
        _require_same_snapshot(
            baseline,
            measured,
            error_type=error_type,
            label=label,
            phase="legacy anchored verification hash pass",
        )
        records = {
            entry.relative: entry.sha256
            for entry in measured.entries
            if entry.entry_type == "file" and entry.sha256 is not None
        }
        if len(records) != expected_file_count:
            _raise(
                error_type,
                label,
                f"legacy retained file count drifted: {len(records)} != {expected_file_count}",
            )
        inventory_payload = "".join(
            f"{digest}  {relative}\n" for relative, digest in sorted(records.items())
        )
        inventory_sha256 = hashlib.sha256(inventory_payload.encode("utf-8")).hexdigest()
        if inventory_sha256 != expected_inventory_sha256:
            _raise(error_type, label, "legacy retained aggregate SHA-256 drifted")
        final = _scan_anchored_tree(
            root_fd,
            hash_regular=False,
            error_type=error_type,
            label=label,
        )
        _require_read_only(final, include_root=True, error_type=error_type, label=label)
        _require_same_snapshot(
            measured,
            final,
            error_type=error_type,
            label=label,
            phase="legacy anchored verification stability pass",
        )
        _require_root_binding(
            root,
            root_fd,
            authorized_root=authorized_root,
            error_type=error_type,
            label=label,
        )
        with tempfile.TemporaryDirectory(prefix="athenak-pic-legacy-tree-") as directory:
            staged_root = Path(directory) / "snapshot"
            staged_root.mkdir()
            for entry in measured.entries:
                if entry.entry_type == "directory" and entry.relative:
                    staged_root.joinpath(*PurePosixPath(entry.relative).parts).mkdir(
                        parents=True,
                        exist_ok=True,
                    )

            def load_member(relative: str, digest: str) -> tuple[bytes, int]:
                return _read_anchored_regular_bytes(
                    root_fd,
                    relative,
                    expected_sha256=digest,
                    error_type=error_type,
                    label=label,
                )

            sealed_snapshot = VerifiedFrozenTree(
                staged_root,
                regular_sha256=records,
                regular_sizes={
                    entry.relative: entry.size
                    for entry in measured.entries
                    if entry.entry_type == "file"
                },
                directories=frozenset(
                    entry.relative
                    for entry in measured.entries
                    if entry.entry_type == "directory" and entry.relative
                ),
                loader=load_member,
            )
            try:
                yield {
                    "root": str(root),
                    "file_count": len(records),
                    "inventory_sha256": inventory_sha256,
                    "recursively_read_only": True,
                    "anchored_sealed_snapshot": True,
                }, sealed_snapshot
            finally:
                sealed_snapshot.close()


def _verify_frozen_tree_anchored(
    root: Path,
    root_fd: int,
    expected_inventory_sha256: str,
    *,
    staged_snapshot: _AnchoredTreeSnapshot | None = None,
    authorized_root: Path,
    error_type: type[ValueError],
    label: str,
) -> dict[str, Any]:
    """Verify one frozen tree while retaining the same root directory object."""
    baseline = _scan_anchored_tree(root_fd, hash_regular=False, error_type=error_type, label=label)
    _require_read_only(baseline, include_root=True, error_type=error_type, label=label)
    if staged_snapshot is not None:
        if INVENTORY_NAME in staged_snapshot.by_relative_path():
            _require_same_snapshot(
                staged_snapshot,
                baseline,
                error_type=error_type,
                label=label,
                phase="sealed snapshot topology capture",
            )
        else:
            _require_published_tree_matches_staging(
                staged_snapshot,
                baseline,
                error_type=error_type,
                label=label,
            )
    measured = _scan_anchored_tree(
        root_fd,
        hash_regular=True,
        capture_paths=frozenset(_RESERVED_METADATA),
        error_type=error_type,
        label=label,
    )
    _require_read_only(measured, include_root=True, error_type=error_type, label=label)
    _require_same_snapshot(
        baseline,
        measured,
        error_type=error_type,
        label=label,
        phase="anchored verification hash pass",
    )
    entries = measured.by_relative_path()
    inventory_entry = entries.get(INVENTORY_NAME)
    if inventory_entry is None or inventory_entry.entry_type != "file":
        _raise(error_type, label, "retained inventory metadata file is absent")
    if inventory_entry.sha256 != expected_inventory_sha256:
        _raise(error_type, label, "retained inventory SHA-256 drifted from anchored digest")
    expected = _parse_inventory(
        _decode_metadata(inventory_entry, INVENTORY_NAME, error_type=error_type, label=label),
        error_type=error_type,
        label=label,
    )
    receipt_entry = entries.get(FREEZE_RECEIPT_NAME)
    if receipt_entry is None or receipt_entry.entry_type != "file":
        _raise(error_type, label, "retained freeze receipt metadata file is absent")
    try:
        receipt = json.loads(
            _decode_metadata(receipt_entry, FREEZE_RECEIPT_NAME, error_type=error_type, label=label)
        )
    except json.JSONDecodeError as error:
        _raise(error_type, label, f"retained freeze receipt is not valid JSON: {error}")
    receipt = _validate_receipt_payload(receipt, error_type=error_type, label=label)
    measured_payloads = {
        relative: entry.sha256
        for relative, entry in entries.items()
        if entry.entry_type == "file" and relative != INVENTORY_NAME
    }
    if set(measured_payloads) != set(expected):
        _raise(error_type, label, "retained tree membership drifted from inventory")
    for relative, digest in expected.items():
        if measured_payloads[relative] != digest:
            _raise(error_type, label, f"retained artifact hash drifted: {relative}")
    final = _scan_anchored_tree(root_fd, hash_regular=False, error_type=error_type, label=label)
    _require_read_only(final, include_root=True, error_type=error_type, label=label)
    _require_same_snapshot(
        measured,
        final,
        error_type=error_type,
        label=label,
        phase="anchored verification stability pass",
    )
    _require_root_binding(
        root,
        root_fd,
        authorized_root=authorized_root,
        error_type=error_type,
        label=label,
    )
    return {
        "root": str(root),
        "inventory_path": INVENTORY_NAME,
        "inventory_excludes": INVENTORY_NAME,
        "inventory_algorithm": INVENTORY_ALGORITHM,
        "inventory_sha256": inventory_entry.sha256,
        "inventoried_file_count": len(expected),
        "writable_entries": [],
        "recursively_read_only": True,
        "freeze_receipt": receipt,
    }


def verify_frozen_tree(
    runtime_root: str | Path,
    expected_inventory_sha256: str,
    *,
    authorized_root: Path,
    error_type: type[ValueError] = ValueError,
    label: str = "immutable-tree",
) -> dict[str, Any]:
    """Verify anchored inventory content, exact membership, hashes, and modes."""
    if not _SHA256_PATTERN.fullmatch(expected_inventory_sha256):
        _raise(error_type, label, "expected inventory SHA-256 must be 64 lowercase hex digits")
    root = _canonical_authorized_root(
        runtime_root,
        authorized_root=authorized_root,
        error_type=error_type,
        label=label,
    )
    with _open_anchored_root(
        root,
        authorized_root=authorized_root,
        error_type=error_type,
        label=label,
    ) as root_fd:
        return _verify_frozen_tree_anchored(
            root,
            root_fd,
            expected_inventory_sha256,
            authorized_root=authorized_root,
            error_type=error_type,
            label=label,
        )


def freeze_tree(
    runtime_root: str | Path,
    receipt: dict[str, Any],
    *,
    authorized_root: Path,
    error_type: type[ValueError] = ValueError,
    label: str = "immutable-tree",
) -> dict[str, Any]:
    """Exclusively write metadata, recursively remove write bits, and verify."""
    root = _canonical_authorized_root(
        runtime_root,
        authorized_root=authorized_root,
        error_type=error_type,
        label=label,
    )
    _validate_receipt_payload(receipt, error_type=error_type, label=label)
    with _open_anchored_root(
        root,
        authorized_root=authorized_root,
        error_type=error_type,
        label=label,
    ) as root_fd:
        preflight = _scan_anchored_tree(
            root_fd, hash_regular=False, error_type=error_type, label=label
        )
        _require_reserved_metadata_absent(root_fd, error_type=error_type, label=label)
        _write_new_metadata(
            root_fd,
            FREEZE_RECEIPT_NAME,
            json.dumps(receipt, indent=2, sort_keys=True) + "\n",
            error_type=error_type,
            label=label,
        )
        _remove_write_bits_below_root(root_fd, error_type=error_type, label=label)
        boundary = _scan_anchored_tree(
            root_fd, hash_regular=False, error_type=error_type, label=label
        )
        _require_freeze_boundary(preflight, boundary, error_type=error_type, label=label)
        _require_read_only(boundary, include_root=False, error_type=error_type, label=label)
        hashed = _scan_anchored_tree(
            root_fd, hash_regular=True, error_type=error_type, label=label
        )
        _require_read_only(hashed, include_root=False, error_type=error_type, label=label)
        _require_same_snapshot(
            boundary,
            hashed,
            error_type=error_type,
            label=label,
            phase="immutable staging hash pass",
        )
        inventory = "".join(
            f"{entry.sha256}  {entry.relative}\n"
            for entry in hashed.entries
            if entry.entry_type == "file" and entry.relative != INVENTORY_NAME
        )
        _write_new_metadata(
            root_fd,
            INVENTORY_NAME,
            inventory,
            error_type=error_type,
            label=label,
        )
        os.fchmod(root_fd, os.fstat(root_fd).st_mode & ~_WRITE_BITS)
        return _verify_frozen_tree_anchored(
            root,
            root_fd,
            hashlib.sha256(inventory.encode("utf-8")).hexdigest(),
            staged_snapshot=hashed,
            authorized_root=authorized_root,
            error_type=error_type,
            label=label,
        )
