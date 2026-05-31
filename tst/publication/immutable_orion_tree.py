#!/usr/bin/env python3
"""Fail-closed helpers for recursively frozen Orion-backed evidence trees."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import stat
import tarfile
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


def _raise(error_type: type[ValueError], label: str, message: str) -> None:
    raise error_type(f"{label}: {message}")


def _sha256_regular(
    path: Path,
    *,
    error_type: type[ValueError],
    label: str,
) -> str:
    """Hash one self-contained regular file without following its final link."""
    try:
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    except OSError as error:
        _raise(error_type, label, f"cannot open retained regular file {path}: {error}")
    try:
        mode = os.fstat(fd).st_mode
        if not stat.S_ISREG(mode):
            _raise(error_type, label, f"retained entry is not a regular file: {path}")
        if os.fstat(fd).st_nlink != 1:
            _raise(error_type, label, f"retained regular file must have one link: {path}")
        digest = hashlib.sha256()
        while chunk := os.read(fd, 1024 * 1024):
            digest.update(chunk)
        return digest.hexdigest()
    finally:
        os.close(fd)


def _read_regular_text(
    path: Path,
    *,
    error_type: type[ValueError],
    label: str,
) -> str:
    """Read one self-contained regular UTF-8 file without following its final link."""
    try:
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    except OSError as error:
        _raise(error_type, label, f"cannot open retained metadata file {path}: {error}")
    try:
        status = os.fstat(fd)
        if not stat.S_ISREG(status.st_mode) or status.st_nlink != 1:
            _raise(error_type, label, f"retained metadata file is unsafe: {path}")
        payload = bytearray()
        while chunk := os.read(fd, 1024 * 1024):
            payload.extend(chunk)
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
    measured_sha256 = _sha256_regular(path, error_type=error_type, label=label)
    if measured_sha256 != expected_sha256:
        _raise(error_type, label, "source archive SHA-256 drifted")
    names = set()
    regular_file_count = 0
    try:
        with tarfile.open(path, "r:*") as archive:
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
                    _raise(error_type, label, f"source archive member name is unsafe: {name!r}")
                if name in names:
                    _raise(error_type, label, f"source archive member is duplicated: {name!r}")
                names.add(name)
                if "__pycache__" in relative.parts or relative.suffix in (".pyc", ".pyo"):
                    _raise(error_type, label, f"source archive contains bytecode cache: {name!r}")
                if member.isdir():
                    continue
                if not member.isreg():
                    _raise(error_type, label, f"source archive member is not regular: {name!r}")
                regular_file_count += 1
    except (tarfile.TarError, OSError) as error:
        _raise(error_type, label, f"cannot decode retained source archive: {error}")
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
    measured_sha256 = _sha256_regular(path, error_type=error_type, label=label)
    if measured_sha256 != expected_sha256:
        _raise(error_type, label, "executable SHA-256 drifted")
    try:
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    except OSError as error:
        _raise(error_type, label, f"cannot open retained executable {path}: {error}")
    try:
        mode = os.fstat(fd).st_mode
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
    try:
        with tarfile.open(Path(archive_path), "r:*") as archive:
            members = {member.name: member for member in archive.getmembers()}
            for name, expected_sha256 in dependencies.items():
                member = members.get(name)
                if member is None or not member.isreg():
                    _raise(error_type, label, f"required archived dependency is absent: {name}")
                extracted = archive.extractfile(member)
                if extracted is None:
                    _raise(error_type, label, f"cannot read archived dependency: {name}")
                measured_sha256 = hashlib.sha256(extracted.read()).hexdigest()
                if measured_sha256 != expected_sha256:
                    _raise(error_type, label, f"archived dependency SHA-256 drifted: {name}")
    except (tarfile.TarError, OSError) as error:
        _raise(error_type, label, f"cannot decode retained source archive: {error}")
    return {"dependency_count": len(dependencies), "passed": True}


def validate_serial_host_build_evidence(
    pinned_root: str | Path,
    *,
    error_type: type[ValueError] = ValueError,
    label: str = "immutable-tree serial-host build evidence",
) -> dict[str, str | bool]:
    """Validate retained serial-host build controls and successful nonempty logs."""
    root = Path(pinned_root)

    def read(relative: str) -> str:
        return _read_regular_text(root / relative, error_type=error_type, label=label)

    try:
        profile = json.loads(read("build/build_profile.json"))
        compile_commands = json.loads(read("build/compile_commands.json"))
    except json.JSONDecodeError as error:
        _raise(error_type, label, f"retained build JSON is invalid: {error}")
    if profile != {
        "schema_version": 1,
        "profile": "bounded_serial_host_clean_build_provenance_only",
        "qualifying_evidence": False,
        "mpi": False,
        "gpu": False,
        "frontier_scheduler": False,
        "kronos": False,
    }:
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
        if expected not in cache:
            _raise(error_type, label, f"retained CMake cache is missing {expected!r}")
    config = read("build/config.hpp")
    for expected in ("#define MPI_PARALLEL_ENABLED 0", "#define OPENMP_PARALLEL_ENABLED 0"):
        if expected not in config:
            _raise(error_type, label, f"retained config.hpp is missing {expected!r}")
    configure = shlex.split(read("build/configure_command.txt"))
    refresh = shlex.split(read("build/configure_compile_commands_refresh_command.txt"))
    build = shlex.split(read("build/build_command.txt"))
    rebuild = shlex.split(read("build/verbose_clean_rebuild_command.txt"))
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
