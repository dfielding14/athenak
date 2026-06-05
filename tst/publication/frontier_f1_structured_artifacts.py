#!/usr/bin/env python3
"""Read one pinned structured F1 artifact tree and publish one immutable result."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import stat


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
TRUSTED_PYTHON = "/opt/cray/pe/python/3.11.7/bin/python3"
REVIEWED_FRONTIER_MPICH_DIAGNOSTIC_SHA256 = (
    "fecd6e9635eb80d47efa65ec0789ad03c2ea1375ac96bfa9c801b3cbb7dc00ac"
)


def validate_frontier_mpich_diagnostic_stderr(
    stderr: str,
    *,
    expected_sha256: str = REVIEWED_FRONTIER_MPICH_DIAGNOSTIC_SHA256,
) -> None:
    """Accept only the exact reviewed Cray MPICH informational transcript."""
    if hashlib.sha256(stderr.encode("utf-8")).hexdigest() != expected_sha256:
        raise ValueError("Athena stderr differs from reviewed Cray MPICH diagnostics")


def _parts(relative: str) -> tuple[str, ...]:
    path = PurePosixPath(relative)
    if (
        not relative
        or path.is_absolute()
        or path.as_posix() != relative
        or not path.parts
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise ValueError(f"Unsafe structured artifact path: {relative!r}")
    if len(path.parts) > MAX_TREE_DEPTH:
        raise ValueError(
            f"Structured artifact path exceeds the {MAX_TREE_DEPTH}-component limit: {relative}"
        )
    return path.parts


def _bounded_directory_names(directory_fd: int, label: str) -> list[str]:
    names: list[str] = []
    with os.scandir(directory_fd) as entries:
        for entry in entries:
            if len(names) >= MAX_DIRECTORY_ENTRIES:
                raise ValueError(
                    f"{label} exceeds the {MAX_DIRECTORY_ENTRIES}-entry directory limit"
                )
            names.append(entry.name)
    return sorted(names)


def _read_bounded_descriptor(descriptor: int, label: str, max_bytes: int) -> bytes:
    metadata = os.fstat(descriptor)
    if metadata.st_size > max_bytes:
        raise ValueError(f"{label} exceeds the {max_bytes}-byte size limit")
    payload = bytearray()
    while True:
        chunk = os.read(
            descriptor,
            min(_READ_CHUNK_BYTES, max_bytes + 1 - len(payload)),
        )
        if not chunk:
            return bytes(payload)
        payload.extend(chunk)
        if len(payload) > max_bytes:
            raise ValueError(f"{label} exceeds the {max_bytes}-byte size limit")


def _read_at(
    root_fd: int,
    relative: str,
    directory_identities: dict[str, tuple[int, int]] | None = None,
    file_identities: dict[str, tuple[int, int]] | None = None,
    *,
    max_bytes: int = MAX_RETAINED_FILE_BYTES,
) -> bytes:
    parts = _parts(relative)
    parent_fd = os.dup(root_fd)
    descriptor: int | None = None
    try:
        for index, part in enumerate(parts[:-1]):
            entry = os.stat(part, dir_fd=parent_fd, follow_symlinks=False)
            child_fd = os.open(part, _DIRECTORY_FLAGS, dir_fd=parent_fd)
            try:
                metadata = os.fstat(child_fd)
                if (
                    not stat.S_ISDIR(metadata.st_mode)
                    or metadata.st_mode & 0o222
                    or (entry.st_dev, entry.st_ino) != (metadata.st_dev, metadata.st_ino)
                ):
                    raise ValueError(
                        f"Structured artifact directory is not read-only: {relative}"
                    )
                child_relative = "/".join(parts[: index + 1])
                if (
                    directory_identities is not None
                    and directory_identities.get(child_relative)
                    != (metadata.st_dev, metadata.st_ino)
                ):
                    raise ValueError(
                        f"Structured artifact directory changed during analysis: {child_relative}"
                    )
            except BaseException:
                os.close(child_fd)
                raise
            os.close(parent_fd)
            parent_fd = child_fd
        entry = os.stat(parts[-1], dir_fd=parent_fd, follow_symlinks=False)
        descriptor = os.open(parts[-1], _FILE_FLAGS, dir_fd=parent_fd)
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_mode & 0o222
            or (entry.st_dev, entry.st_ino) != (before.st_dev, before.st_ino)
        ):
            raise ValueError(f"Structured artifact is not a read-only regular file: {relative}")
        if (
            file_identities is not None
            and file_identities.setdefault(relative, (before.st_dev, before.st_ino))
            != (before.st_dev, before.st_ino)
        ):
            raise ValueError(f"Structured artifact changed during analysis: {relative}")
        data = _read_bounded_descriptor(descriptor, relative, max_bytes)
        after = os.fstat(descriptor)
        entry = os.stat(parts[-1], dir_fd=parent_fd, follow_symlinks=False)
        stable = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns")
        if (
            any(getattr(before, field) != getattr(after, field) for field in stable)
            or (entry.st_dev, entry.st_ino) != (after.st_dev, after.st_ino)
            or len(data) != after.st_size
        ):
            raise ValueError(f"Structured artifact changed while reading: {relative}")
        return data
    finally:
        if descriptor is not None:
            os.close(descriptor)
        os.close(parent_fd)


def _open_absolute_directory(path: Path) -> int:
    absolute = Path(os.path.abspath(path))
    descriptor = os.open("/", _DIRECTORY_FLAGS)
    try:
        for part in absolute.parts[1:]:
            child_fd = os.open(part, _DIRECTORY_FLAGS, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child_fd
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _reject_duplicate_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Structured artifact inventory repeats JSON key: {key}")
        value[key] = item
    return value


def canonical_json_bytes(value: dict[str, object]) -> bytes:
    try:
        payload = (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except RecursionError as error:
        raise ValueError("JSON document exceeds the supported nesting depth") from error
    if len(payload) > MAX_JSON_BYTES:
        raise ValueError(f"JSON document exceeds the {MAX_JSON_BYTES}-byte size limit")
    return payload


def read_only_file_sha256(path: Path) -> str:
    descriptor = os.open(path, _FILE_FLAGS)
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_mode & 0o222:
            raise ValueError(f"Offline analysis source is not read-only: {path}")
        data = _read_bounded_descriptor(descriptor, str(path), MAX_RETAINED_FILE_BYTES)
        after = os.fstat(descriptor)
        if (
            (metadata.st_dev, metadata.st_ino, metadata.st_mode, metadata.st_size)
            != (after.st_dev, after.st_ino, after.st_mode, after.st_size)
            or len(data) != after.st_size
        ):
            raise ValueError(f"Offline analysis source changed while reading: {path}")
        return hashlib.sha256(data).hexdigest()
    finally:
        os.close(descriptor)


def offline_analysis_receipt(
    tree: "StructuredArtifactTree",
    *,
    analyzer_path: Path,
    analyzer_sha256: str,
    support_module_sha256: str,
    result: dict[str, object],
) -> dict[str, object]:
    helper_path = Path(__file__)
    for digest in (analyzer_sha256, support_module_sha256):
        if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
            raise ValueError("Offline analysis source checksum is malformed")
    return {
        "schema_version": 1,
        "runner": {
            "python": TRUSTED_PYTHON,
            "flags": ["-I", "-B"],
        },
        "analyzer": {
            "path": analyzer_path.name,
            "sha256": analyzer_sha256,
        },
        "support_modules": [
            {
                "path": helper_path.name,
                "sha256": support_module_sha256,
            }
        ],
        "artifact_inventory": {
            "path": "artifact_inventory.json",
            "sha256": hashlib.sha256(
                tree.read("artifact_inventory.json", max_bytes=MAX_JSON_BYTES)
            ).hexdigest(),
        },
        "analysis_result": {
            "path": "analysis/analysis.json",
            "sha256": hashlib.sha256(canonical_json_bytes(result)).hexdigest(),
        },
    }


class StructuredArtifactTree:
    """Retain one no-follow artifact-root descriptor for an entire analysis."""

    def __init__(
        self, artifact_dir: Path, *, inherited_root_fd: int | None = None
    ) -> None:
        self.artifact_dir = Path(os.path.abspath(artifact_dir))
        self._inherited_root_fd = inherited_root_fd
        self._root_fd: int | None = None
        self._analysis_fd: int | None = None
        self._directory_identities: dict[str, tuple[int, int]] | None = None
        self._file_identities: dict[str, tuple[int, int]] | None = None
        self._inventory_bytes: bytes | None = None
        self._inventory: dict[str, dict[str, object]] | None = None
        self._analysis_result_bindings: dict[str, tuple[int, bytes]] = {}

    @property
    def root_fd(self) -> int:
        if self._root_fd is None:
            raise ValueError("Structured artifact tree is not open")
        return self._root_fd

    @property
    def analysis_fd(self) -> int:
        if self._analysis_fd is None:
            raise ValueError("Structured artifact analysis directory is not open")
        return self._analysis_fd

    def __enter__(self) -> "StructuredArtifactTree":
        self._root_fd = (
            _open_absolute_directory(self.artifact_dir)
            if self._inherited_root_fd is None
            else os.dup(self._inherited_root_fd)
        )
        try:
            if not stat.S_ISDIR(os.fstat(self.root_fd).st_mode):
                raise ValueError("Structured artifact root descriptor is not a directory")
            self.require_path_identity()
            if os.fstat(self.root_fd).st_mode & 0o222:
                raise ValueError("Structured artifact root is not read-only")
            analysis_entry = os.stat("analysis", dir_fd=self.root_fd, follow_symlinks=False)
            self._analysis_fd = os.open("analysis", _DIRECTORY_FLAGS, dir_fd=self.root_fd)
            analysis = os.fstat(self._analysis_fd)
            if (analysis_entry.st_dev, analysis_entry.st_ino) != (
                analysis.st_dev,
                analysis.st_ino,
            ):
                raise ValueError("Structured artifact analysis directory changed during analysis")
            self.require_analysis_identity()
            return self
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def __exit__(self, *_: object) -> None:
        for descriptor, _ in self._analysis_result_bindings.values():
            os.close(descriptor)
        self._analysis_result_bindings.clear()
        if self._analysis_fd is not None:
            os.close(self._analysis_fd)
            self._analysis_fd = None
        if self._root_fd is not None:
            os.close(self._root_fd)
            self._root_fd = None
        self._directory_identities = None
        self._file_identities = None
        self._inventory_bytes = None
        self._inventory = None

    def require_path_identity(self) -> None:
        if self._inherited_root_fd is not None:
            expected = os.fstat(self.root_fd)
            actual = os.fstat(self._inherited_root_fd)
            if (expected.st_dev, expected.st_ino) != (actual.st_dev, actual.st_ino):
                raise ValueError("Structured artifact root descriptor changed during analysis")
            return
        actual_fd = _open_absolute_directory(self.artifact_dir)
        try:
            expected = os.fstat(self.root_fd)
            actual = os.fstat(actual_fd)
            if (expected.st_dev, expected.st_ino) != (actual.st_dev, actual.st_ino):
                raise ValueError("Structured artifact root path changed during analysis")
        finally:
            os.close(actual_fd)

    def require_analysis_identity(self) -> None:
        expected = os.fstat(self.analysis_fd)
        actual = os.stat("analysis", dir_fd=self.root_fd, follow_symlinks=False)
        if (
            not stat.S_ISDIR(expected.st_mode)
            or not stat.S_ISDIR(actual.st_mode)
            or stat.S_IMODE(expected.st_mode) != 0o700
            or stat.S_IMODE(actual.st_mode) != 0o700
            or (expected.st_dev, expected.st_ino) != (actual.st_dev, actual.st_ino)
        ):
            raise ValueError("Structured artifact analysis directory changed during analysis")
        for name in _bounded_directory_names(
            self.analysis_fd,
            "Structured artifact analysis directory",
        ):
            if name not in {"analysis.json", "offline_analysis_receipt.json"}:
                raise ValueError("Structured artifact analysis directory has an unexpected entry")
            metadata = os.stat(name, dir_fd=self.analysis_fd, follow_symlinks=False)
            if not stat.S_ISREG(metadata.st_mode) or metadata.st_mode & 0o222:
                raise ValueError("Structured artifact analysis result is not read-only")
            if metadata.st_size > MAX_JSON_BYTES:
                raise ValueError(
                    f"Structured artifact analysis result exceeds the {MAX_JSON_BYTES}-byte "
                    "size limit"
                )

    def read(self, relative: str, *, max_bytes: int = MAX_RETAINED_FILE_BYTES) -> bytes:
        self.require_path_identity()
        if self._inventory is not None:
            self.require_tree_closure()
        data = _read_at(
            self.root_fd,
            relative,
            self._directory_identities,
            self._file_identities,
            max_bytes=max_bytes,
        )
        if self._inventory is not None:
            self.require_tree_closure()
        self.require_path_identity()
        return data

    def _tree_files(
        self,
        directory_fd: int,
        prefix: tuple[str, ...] = (),
        directory_identities: dict[str, tuple[int, int]] | None = None,
        file_identities: dict[str, tuple[int, int]] | None = None,
        *,
        _entry_count: list[int] | None = None,
        _depth: int = 0,
    ) -> list[str]:
        if _depth > MAX_TREE_DEPTH:
            raise ValueError(
                f"Structured artifact tree exceeds the {MAX_TREE_DEPTH}-level depth limit"
            )
        if _entry_count is None:
            _entry_count = [0]
        paths = []
        label = "/".join(prefix) or "."
        for name in _bounded_directory_names(
            directory_fd,
            f"Structured artifact directory {label}",
        ):
            if _entry_count[0] >= MAX_TREE_ENTRIES:
                raise ValueError(
                    f"Structured artifact tree exceeds the {MAX_TREE_ENTRIES}-entry limit"
                )
            _entry_count[0] += 1
            if not prefix and name == "artifact_inventory.json":
                continue
            if not prefix and name == "analysis":
                analysis_fd = os.open(name, _DIRECTORY_FLAGS, dir_fd=directory_fd)
                try:
                    if not stat.S_ISDIR(os.fstat(analysis_fd).st_mode):
                        raise ValueError("Structured artifact analysis entry is not a directory")
                finally:
                    os.close(analysis_fd)
                continue
            relative = "/".join((*prefix, name))
            metadata = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            if stat.S_ISREG(metadata.st_mode):
                if metadata.st_mode & 0o222:
                    raise ValueError(f"Structured artifact is not read-only: {relative}")
                descriptor = os.open(name, _FILE_FLAGS, dir_fd=directory_fd)
                try:
                    opened = os.fstat(descriptor)
                    if (
                        not stat.S_ISREG(opened.st_mode)
                        or opened.st_mode & 0o222
                        or (metadata.st_dev, metadata.st_ino)
                        != (opened.st_dev, opened.st_ino)
                    ):
                        raise ValueError(
                            f"Structured artifact changed during analysis: {relative}"
                        )
                    if opened.st_size > MAX_RETAINED_FILE_BYTES:
                        raise ValueError(
                            f"Structured artifact exceeds the {MAX_RETAINED_FILE_BYTES}-byte "
                            f"size limit: {relative}"
                        )
                    if file_identities is not None:
                        identity = (opened.st_dev, opened.st_ino)
                        expected = file_identities.setdefault(relative, identity)
                        if expected != identity:
                            raise ValueError(
                                f"Structured artifact changed during analysis: {relative}"
                            )
                finally:
                    os.close(descriptor)
                paths.append(relative)
            elif stat.S_ISDIR(metadata.st_mode):
                if metadata.st_mode & 0o222:
                    raise ValueError(f"Structured artifact directory is not read-only: {relative}")
                child_fd = os.open(name, _DIRECTORY_FLAGS, dir_fd=directory_fd)
                try:
                    child_metadata = os.fstat(child_fd)
                    if (metadata.st_dev, metadata.st_ino) != (
                        child_metadata.st_dev,
                        child_metadata.st_ino,
                    ):
                        raise ValueError(
                            f"Structured artifact directory changed during analysis: {relative}"
                        )
                    if directory_identities is not None:
                        directory_identities[relative] = (
                            child_metadata.st_dev,
                            child_metadata.st_ino,
                        )
                    child_paths = self._tree_files(
                        child_fd,
                        (*prefix, name),
                        directory_identities,
                        file_identities,
                        _entry_count=_entry_count,
                        _depth=_depth + 1,
                    )
                    if not child_paths:
                        raise ValueError(
                            f"Structured artifact tree contains an empty directory: {relative}"
                        )
                    paths.extend(child_paths)
                    after = os.fstat(child_fd)
                    entry = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
                    if (entry.st_dev, entry.st_ino) != (after.st_dev, after.st_ino):
                        raise ValueError(
                            f"Structured artifact directory changed during analysis: {relative}"
                        )
                finally:
                    os.close(child_fd)
            else:
                raise ValueError(f"Structured artifact tree has unsupported entry: {relative}")
        return paths

    def require_tree_closure(self) -> None:
        if (
            self._inventory is None
            or self._directory_identities is None
            or self._file_identities is None
            or self._inventory_bytes is None
        ):
            raise ValueError("Structured artifact inventory is not loaded")
        self.require_path_identity()
        if os.fstat(self.root_fd).st_mode & 0o222:
            raise ValueError("Structured artifact root is not read-only")
        self.require_analysis_identity()
        identities: dict[str, tuple[int, int]] = {}
        file_identities: dict[str, tuple[int, int]] = {}
        inventory_bytes = _read_at(
            self.root_fd,
            "artifact_inventory.json",
            file_identities=file_identities,
            max_bytes=MAX_JSON_BYTES,
        )
        if inventory_bytes != self._inventory_bytes:
            raise ValueError("Structured artifact inventory changed during analysis")
        if self._tree_files(
            self.root_fd,
            directory_identities=identities,
            file_identities=file_identities,
        ) != list(self._inventory):
            raise ValueError("Structured artifact tree differs from its exact inventory")
        if identities != self._directory_identities:
            raise ValueError("Structured artifact directory identity changed during analysis")
        if file_identities != self._file_identities:
            raise ValueError("Structured artifact file identity changed during analysis")
        for relative, record in self._inventory.items():
            data = _read_at(
                self.root_fd,
                relative,
                identities,
                file_identities,
                max_bytes=min(MAX_RETAINED_FILE_BYTES, int(record["size"])),
            )
            if (
                len(data) != record["size"]
                or hashlib.sha256(data).hexdigest() != record["sha256"]
            ):
                raise ValueError(f"Structured artifact inventory checksum mismatch: {relative}")
        self.require_analysis_identity()
        self.require_path_identity()

    def load_inventory(self) -> dict[str, dict[str, object]]:
        self.require_path_identity()
        file_identities: dict[str, tuple[int, int]] = {}
        inventory_bytes = _read_at(
            self.root_fd,
            "artifact_inventory.json",
            file_identities=file_identities,
            max_bytes=MAX_JSON_BYTES,
        )
        try:
            value = json.loads(
                inventory_bytes.decode("utf-8"),
                object_pairs_hook=_reject_duplicate_keys,
            )
        except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as error:
            raise ValueError("Structured artifact inventory is not UTF-8 JSON") from error
        if (
            not isinstance(value, dict)
            or set(value) != {"schema_version", "files"}
            or type(value.get("schema_version")) is not int
            or value["schema_version"] != 1
            or not isinstance(value.get("files"), list)
        ):
            raise ValueError("Structured artifact inventory is malformed")
        if len(value["files"]) > MAX_TREE_ENTRIES:
            raise ValueError(
                f"Structured artifact inventory exceeds the {MAX_TREE_ENTRIES}-record limit"
            )
        result: dict[str, dict[str, object]] = {}
        for record in value["files"]:
            if (
                not isinstance(record, dict)
                or set(record) != {"path", "sha256", "size"}
                or not isinstance(record["path"], str)
                or not isinstance(record["sha256"], str)
                or len(record["sha256"]) != 64
                or any(character not in "0123456789abcdef" for character in record["sha256"])
                or not isinstance(record["size"], int)
                or isinstance(record["size"], bool)
                or record["size"] < 0
            ):
                raise ValueError("Structured artifact inventory record is malformed")
            _parts(record["path"])
            if record["path"] in result:
                raise ValueError("Structured artifact inventory has duplicate paths")
            result[record["path"]] = record
        if list(result) != sorted(result):
            raise ValueError("Structured artifact inventory paths are not sorted")
        self.require_path_identity()
        directory_identities: dict[str, tuple[int, int]] = {}
        if self._tree_files(
            self.root_fd,
            directory_identities=directory_identities,
            file_identities=file_identities,
        ) != list(result):
            raise ValueError("Structured artifact tree differs from its exact inventory")
        for relative in result:
            data = _read_at(
                self.root_fd,
                relative,
                directory_identities,
                file_identities,
                max_bytes=min(MAX_RETAINED_FILE_BYTES, int(result[relative]["size"])),
            )
            if (
                len(data) != result[relative]["size"]
                or hashlib.sha256(data).hexdigest() != result[relative]["sha256"]
            ):
                raise ValueError(f"Structured artifact inventory checksum mismatch: {relative}")
        self._directory_identities = directory_identities
        self._file_identities = file_identities
        self._inventory_bytes = inventory_bytes
        self._inventory = result
        self.require_tree_closure()
        self.require_path_identity()
        return result

    def read_inventory_bytes(
        self,
        inventory: dict[str, dict[str, object]],
        relative: str,
    ) -> bytes:
        if inventory is not self._inventory:
            raise ValueError("Structured artifact inventory is not the loaded inventory")
        record = inventory.get(relative)
        if record is None:
            raise ValueError(f"Structured artifact inventory omits required path: {relative}")
        data = self.read(
            relative,
            max_bytes=min(MAX_RETAINED_FILE_BYTES, int(record["size"])),
        )
        if (
            len(data) != record["size"]
            or hashlib.sha256(data).hexdigest() != record["sha256"]
        ):
            raise ValueError(f"Structured artifact inventory checksum mismatch: {relative}")
        return data

    def _verify_retained_analysis_results(self) -> None:
        for name, (descriptor, expected_bytes) in self._analysis_result_bindings.items():
            os.lseek(descriptor, 0, os.SEEK_SET)
            before = os.fstat(descriptor)
            data = _read_bounded_descriptor(
                descriptor,
                f"analysis/{name}",
                MAX_JSON_BYTES,
            )
            after = os.fstat(descriptor)
            entry = os.stat(name, dir_fd=self.analysis_fd, follow_symlinks=False)
            stable = (
                "st_dev",
                "st_ino",
                "st_mode",
                "st_size",
                "st_mtime_ns",
                "st_ctime_ns",
            )
            if (
                not stat.S_ISREG(after.st_mode)
                or after.st_mode & 0o222
                or any(getattr(before, field) != getattr(after, field) for field in stable)
                or (entry.st_dev, entry.st_ino) != (after.st_dev, after.st_ino)
                or data != expected_bytes
                or len(data) != after.st_size
            ):
                raise ValueError("Published analysis result changed after publication")

    def write_result_exclusive(self, relative: str, result: dict[str, object]) -> None:
        parts = _parts(relative)
        if parts not in {
            ("analysis", "analysis.json"),
            ("analysis", "offline_analysis_receipt.json"),
        }:
            raise ValueError("Structured analysis result path is not authorized")
        name = parts[-1]
        if (
            name == "offline_analysis_receipt.json"
            and "analysis.json" not in self._analysis_result_bindings
        ):
            raise ValueError("Offline analysis receipt requires a published analysis result")
        self.require_tree_closure()
        self.require_analysis_identity()
        self._verify_retained_analysis_results()
        analysis_fd = os.dup(self.analysis_fd)
        descriptor: int | None = None
        try:
            before = os.fstat(analysis_fd)
            descriptor = os.open(
                name,
                os.O_RDWR | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
                0o600,
                dir_fd=analysis_fd,
            )
            data = canonical_json_bytes(result)
            with os.fdopen(descriptor, "wb", closefd=False) as stream:
                stream.write(data)
                stream.flush()
            os.fsync(descriptor)
            os.fchmod(descriptor, 0o444)
            os.fsync(descriptor)
            expected = os.fstat(descriptor)
            actual = os.stat(name, dir_fd=analysis_fd, follow_symlinks=False)
            analysis_entry = os.stat("analysis", dir_fd=self.root_fd, follow_symlinks=False)
            if (
                not stat.S_ISREG(actual.st_mode)
                or actual.st_mode & 0o222
                or (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino)
                or (analysis_entry.st_dev, analysis_entry.st_ino)
                != (before.st_dev, before.st_ino)
            ):
                raise ValueError("Published analysis result changed before directory sync")
            os.fsync(analysis_fd)
            self._analysis_result_bindings[name] = (descriptor, data)
            descriptor = None
            self.require_analysis_identity()
            self.require_tree_closure()
            self.require_path_identity()
            self._verify_retained_analysis_results()
        finally:
            if descriptor is not None:
                os.close(descriptor)
            os.close(analysis_fd)


def load_inventory(tree: StructuredArtifactTree) -> dict[str, dict[str, object]]:
    return tree.load_inventory()


def read_inventory_bytes(
    tree: StructuredArtifactTree,
    inventory: dict[str, dict[str, object]],
    relative: str,
) -> bytes:
    return tree.read_inventory_bytes(inventory, relative)


def require_inventory_sha256(tree: StructuredArtifactTree, expected: str) -> None:
    if (
        len(expected) != 64
        or any(character not in "0123456789abcdef" for character in expected)
        or hashlib.sha256(
            tree.read("artifact_inventory.json", max_bytes=MAX_JSON_BYTES)
        ).hexdigest()
        != expected
    ):
        raise ValueError("Structured artifact inventory differs from parent binding")


def write_result_exclusive(
    tree: StructuredArtifactTree, relative: str, result: dict[str, object]
) -> None:
    tree.write_result_exclusive(relative, result)
