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
import hashlib
import json
import os
from pathlib import Path
from pathlib import PurePosixPath
import re
import shutil
import stat
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
_FILE_FLAGS = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_RENAME_NOREPLACE = 1
_REGISTERED_EXECUTION_PREREGISTRATION_PATH = (
    pilot.REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_pressure_pilot_registered_execution_preregistration_2026-06-02.json"
)
AUTHORIZED_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")


class PressurePilotPublicationError(ValueError):
    """Raised when aggregate raw-bundle publication fails closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PressurePilotPublicationError(message)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical_json_bytes(value: Mapping[str, object]) -> bytes:
    try:
        return (
            json.dumps(dict(value), indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise PressurePilotPublicationError("pressure-pilot manifest is not canonical JSON") from error


def _lexists(path: Path) -> bool:
    return os.path.lexists(path)


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, _DIRECTORY_FLAGS)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


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


def _write_exclusive(path: Path, payload: bytes) -> None:
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o600,
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


def _freeze_directories(root: Path) -> None:
    directories = [Path(directory) for directory, _, _ in os.walk(root)]
    for directory in reversed(directories):
        os.chmod(directory, 0o555, follow_symlinks=False)
        _fsync_directory(directory)


def _cleanup_staging(root: Path) -> None:
    if not _lexists(root):
        return
    for directory, names, filenames in os.walk(root, topdown=False):
        base = Path(directory)
        for name in filenames:
            os.chmod(base / name, 0o600, follow_symlinks=False)
        for name in names:
            os.chmod(base / name, 0o700, follow_symlinks=False)
        os.chmod(base, 0o700, follow_symlinks=False)
    shutil.rmtree(root)


def _rename_no_replace_at(parent_descriptor: int, source_name: str, destination_name: str) -> None:
    _require(
        "/" not in source_name and "/" not in destination_name,
        "pressure-pilot descriptor-relative rename received a nested path",
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
    raise PressurePilotPublicationError(
        "pressure-pilot publication requires atomic no-replace rename support"
    )


def _read_stable_readonly_regular(path: Path, label: str) -> bytes:
    descriptor = os.open(path, _FILE_FLAGS)
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode) and not before.st_mode & 0o222,
            f"{label} is not a read-only regular file",
        )
        payload = b""
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            payload += chunk
        after = os.fstat(descriptor)
        current = os.stat(path, follow_symlinks=False)
        _require(
            (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
            == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino),
            f"{label} changed while reading",
        )
        return payload
    finally:
        os.close(descriptor)


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

    def __init__(self, root: str | Path) -> None:
        self.root = Path(os.path.abspath(root))
        self._root_fd: int | None = None
        self._directory_identities: dict[str, tuple[int, int]] | None = None
        self._file_identities: dict[str, tuple[int, int]] | None = None

    @property
    def root_fd(self) -> int:
        if self._root_fd is None:
            raise PressurePilotPublicationError("pressure-pilot bundle verifier is not open")
        return self._root_fd

    def __enter__(self) -> "ImmutablePressurePilotBundle":
        self._root_fd = _open_absolute_directory(self.root)
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

    def _scan(
        self,
        directory_fd: int,
        prefix: tuple[str, ...] = (),
    ) -> tuple[list[str], dict[str, tuple[int, int]], dict[str, tuple[int, int]]]:
        paths = []
        directories = {}
        files = {}
        names = sorted(os.listdir(directory_fd))
        for name in names:
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
                        child_fd, (*prefix, name)
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

    def read(self, relative: str) -> bytes:
        self.require_path_identity()
        payload = _ARTIFACT_HELPERS._read_at(
            self.root_fd,
            relative,
            self._directory_identities,
            self._file_identities,
        )
        self.require_path_identity()
        return payload

    def verify(self, expected_manifest_sha256: str) -> dict[str, object]:
        _require(
            _SHA256_PATTERN.fullmatch(expected_manifest_sha256) is not None,
            "expected pressure-pilot manifest SHA-256 is malformed",
        )
        manifest_payload = self.read(pilot.MANIFEST_NAME)
        _require(_sha256(manifest_payload) == expected_manifest_sha256, "pressure-pilot manifest SHA-256 drifted")
        manifest = pilot._manifest_schema(manifest_payload)
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


def verify_published_pressure_pilot_bundle(
    output_path: str | Path,
    expected_manifest_sha256: str,
    *,
    authorized_publication_root: Path = pilot.AUTHORIZED_PUBLICATION_ROOT,
) -> dict[str, Any]:
    """Verify immutable closure, run the strict analyzer, then recheck closure."""
    root = Path(os.path.abspath(output_path))
    with ImmutablePressurePilotBundle(root) as bundle:
        bundle.verify(expected_manifest_sha256)
        result = pilot.analyze_pressure_pilot_bundle(
            root,
            expected_manifest_sha256,
            authorized_publication_root=authorized_publication_root,
        )
        bundle.verify(expected_manifest_sha256)
        return result


def _source_bindings() -> dict[str, dict[str, str]]:
    publisher_source = Path(__file__).resolve()
    analyzer_source = Path(pilot.__file__).resolve()
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
        "registered_execution_preregistration": {
            "path": _REGISTERED_EXECUTION_PREREGISTRATION_PATH.relative_to(
                pilot.REPO_ROOT
            ).as_posix(),
            "sha256": _sha256(registered_execution_payload),
        },
        "publisher": {
            "path": publisher_source.relative_to(pilot.REPO_ROOT).as_posix(),
            "sha256": _sha256(
                pilot._regular_bytes(
                    publisher_source, "pressure-pilot publisher source"
                )
            ),
        },
        "aggregate_analyzer": {
            "path": analyzer_source.relative_to(pilot.REPO_ROOT).as_posix(),
            "sha256": _sha256(
                pilot._regular_bytes(
                    analyzer_source, "pressure-pilot aggregate analyzer source"
                )
            ),
        },
    }


def _manifest(case_descriptors: Mapping[str, Mapping[str, object]]) -> dict[str, object]:
    policy = pilot._load_policy()
    preregistration_payload = pilot._regular_bytes(
        pilot.PREREGISTRATION_PATH, "pressure-pilot preregistration"
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
    _require(not _lexists(target), "pressure-pilot output already exists")
    _require(not _lexists(receipt_target), "pressure-pilot receipt already exists")
    _require(not _lexists(result_target), "pressure-pilot result already exists")

    staging = parent / f".{target.name}.staging-{uuid.uuid4()}"
    receipt_staging = parent / f".{receipt_target.name}.staging-{uuid.uuid4()}"
    result_staging = parent / f".{result_target.name}.staging-{uuid.uuid4()}"
    publication_descriptor = _open_absolute_directory(parent)
    _require_same_directory(parent, publication_descriptor, "authorized PIC publication root")
    os.mkdir(staging.name, mode=0o700, dir_fd=publication_descriptor)
    staging_descriptor = os.open(staging.name, _DIRECTORY_FLAGS, dir_fd=publication_descriptor)
    renamed = False
    receipt_renamed = False
    result_renamed = False
    try:
        with ExitStack() as stack:
            trees = {}
            descriptors = {}
            inventories = {}
            for case_id in case_verifier.CASE_IDS:
                tree = stack.enter_context(
                    case_verifier.StructuredArtifactTree(Path(case_artifact_dirs[case_id]))
                )
                descriptor = case_verifier.verify_published_case_descriptor(
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
            _require_same_directory(staging, staging_descriptor, "pressure-pilot staging tree")
            _freeze_directories(staging)
            _require_same_directory(staging, staging_descriptor, "pressure-pilot staging tree")
            result = verify_published_pressure_pilot_bundle(
                staging, manifest_sha256, authorized_publication_root=parent
            )
            result_payload = _canonical_json_bytes(result)
            result_sha256 = _sha256(result_payload)
            _write_exclusive_at(publication_descriptor, result_staging.name, result_payload)
            receipt = {
                "schema_version": 1,
                "record_type": "q011_section54_pressure_pilot_bundle_publication_receipt",
                "evidence_class": pilot.EVIDENCE_CLASS,
                "qualification_effect": pilot.QUALIFICATION_EFFECT,
                "aggregate_bundle": {
                    "path": str(target),
                    "manifest_sha256": manifest_sha256,
                },
                "aggregate_analysis": {
                    "path": str(result_target),
                    "sha256": result_sha256,
                },
                "source_bindings": _source_bindings(),
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
            for tree in trees.values():
                tree.require_tree_closure()
            _require_same_directory(parent, publication_descriptor, "authorized PIC publication root")
            _require(not _lexists(target), "pressure-pilot output appeared during publication")
            _rename_no_replace_at(publication_descriptor, staging.name, target.name)
            renamed = True
            _fsync_descriptor(publication_descriptor)
            verified_result = verify_published_pressure_pilot_bundle(
                target, manifest_sha256, authorized_publication_root=parent
            )
            _require(
                _canonical_json_bytes(verified_result) == result_payload,
                "pressure-pilot aggregate analysis drifted after publication",
            )
            _require_same_directory(parent, publication_descriptor, "authorized PIC publication root")
            _require(not _lexists(result_target), "pressure-pilot result appeared during publication")
            _rename_no_replace_at(
                publication_descriptor, result_staging.name, result_target.name
            )
            result_renamed = True
            _fsync_descriptor(publication_descriptor)
            _require_same_directory(parent, publication_descriptor, "authorized PIC publication root")
            _require(not _lexists(receipt_target), "pressure-pilot receipt appeared during publication")
            _rename_no_replace_at(
                publication_descriptor, receipt_staging.name, receipt_target.name
            )
            receipt_renamed = True
            _fsync_descriptor(publication_descriptor)
        verified_result = verify_published_pressure_pilot_bundle(
            target, manifest_sha256, authorized_publication_root=parent
        )
        _require(
            _canonical_json_bytes(verified_result) == result_payload,
            "pressure-pilot aggregate analysis drifted after receipt publication",
        )
        _require(
            _read_stable_readonly_regular(
                result_target, "pressure-pilot aggregate analysis result"
            )
            == result_payload,
            "pressure-pilot aggregate analysis result drifted after publication",
        )
        return {
            "path": str(target),
            "manifest_sha256": manifest_sha256,
            "analysis_result_path": str(result_target),
            "analysis_result_sha256": result_sha256,
            "receipt_path": str(receipt_target),
            "receipt_sha256": receipt_sha256,
        }
    except BaseException:
        if not renamed:
            _cleanup_staging(staging)
        if not renamed and not receipt_renamed and _lexists(receipt_staging):
            os.chmod(receipt_staging, 0o600, follow_symlinks=False)
            os.unlink(receipt_staging)
        if not result_renamed and _lexists(result_staging):
            os.chmod(result_staging, 0o600, follow_symlinks=False)
            os.unlink(result_staging)
        raise
    finally:
        os.close(staging_descriptor)
        os.close(publication_descriptor)


def verify_published_pressure_pilot_receipt(
    receipt_path: str | Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    """Re-audit one retained publication receipt and every artifact it binds."""
    pic_root, publication_root = _publication_root(authorized_pic_root)
    receipt_target = _direct_publication_target(
        receipt_path, publication_root, "pressure-pilot receipt"
    )
    receipt_payload = _read_stable_readonly_regular(
        receipt_target, "pressure-pilot receipt"
    )
    receipt = pilot._decode_json(receipt_payload, "pressure-pilot receipt")
    _require(isinstance(receipt, dict), "pressure-pilot receipt is malformed")
    _require(
        set(receipt)
        == {
            "schema_version",
            "record_type",
            "evidence_class",
            "qualification_effect",
            "aggregate_bundle",
            "aggregate_analysis",
            "source_bindings",
            "raw_cases",
        },
        "pressure-pilot receipt schema drifted",
    )
    _require(
        receipt["schema_version"] == 1
        and receipt["record_type"]
        == "q011_section54_pressure_pilot_bundle_publication_receipt"
        and receipt["evidence_class"] == pilot.EVIDENCE_CLASS
        and receipt["qualification_effect"] == pilot.QUALIFICATION_EFFECT,
        "pressure-pilot receipt identity drifted",
    )
    _require(
        receipt["source_bindings"] == _source_bindings(),
        "pressure-pilot receipt source bindings drifted",
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
            descriptor = case_verifier.verify_published_case_descriptor(
                tree, case_id, str(case["descriptor_sha256"])
            )
            _require(
                case["descriptor_path"] == case_verifier.CASE_DESCRIPTOR_PATH
                and case["artifact_inventory_sha256"]
                == descriptor["artifact_inventory_sha256"]
                and case["runtime_artifacts"] == descriptor["runtime_artifacts"],
                "pressure-pilot retained raw-case provenance drifted",
            )
    result = verify_published_pressure_pilot_bundle(
        bundle_root,
        manifest_sha256,
        authorized_publication_root=publication_root,
    )
    expected_result_payload = _canonical_json_bytes(result)
    actual_result_payload = _read_stable_readonly_regular(
        analysis_result, "pressure-pilot retained aggregate analysis"
    )
    _require(
        actual_result_payload == expected_result_payload
        and _sha256(actual_result_payload) == analysis_sha256,
        "pressure-pilot retained aggregate analysis drifted",
    )
    return {
        "receipt_sha256": _sha256(receipt_payload),
        "manifest_sha256": manifest_sha256,
        "analysis_result_sha256": analysis_sha256,
        "status": result["status"],
    }


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
