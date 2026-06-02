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
_AT_FDCWD = -100
_RENAME_NOREPLACE = 1


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


def _renameat2_no_replace(source: Path, destination: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise OSError(errno.ENOSYS, os.strerror(errno.ENOSYS), str(destination))
    renameat2.argtypes = [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    ]
    renameat2.restype = ctypes.c_int
    if renameat2(
        _AT_FDCWD,
        os.fsencode(source),
        _AT_FDCWD,
        os.fsencode(destination),
        _RENAME_NOREPLACE,
    ) == 0:
        return
    error_number = ctypes.get_errno()
    raise OSError(error_number, os.strerror(error_number), str(destination))


def _rename_no_replace(source: Path, destination: Path) -> None:
    try:
        _renameat2_no_replace(source, destination)
        return
    except OSError as error:
        error_number = error.errno
    if error_number == errno.EEXIST:
        raise PressurePilotPublicationError(f"pressure-pilot output collided during rename: {destination}")
    unsupported = {
        errno.EINVAL,
        errno.ENOSYS,
        getattr(errno, "ENOTSUP", errno.EINVAL),
        getattr(errno, "EOPNOTSUPP", errno.EINVAL),
    }
    if error_number not in unsupported:
        raise
    before = os.lstat(source)
    _require(not _lexists(destination), "pressure-pilot output appeared during publication")
    os.rename(source, destination)
    after = os.lstat(destination)
    _require(
        (before.st_dev, before.st_ino) == (after.st_dev, after.st_ino),
        "pressure-pilot tree identity changed during rename",
    )
    _require(not _lexists(source), "pressure-pilot staging tree remains after rename")


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
) -> dict[str, Any]:
    """Verify immutable closure, run the strict analyzer, then recheck closure."""
    root = Path(os.path.abspath(output_path))
    with ImmutablePressurePilotBundle(root) as bundle:
        bundle.verify(expected_manifest_sha256)
        result = pilot.analyze_pressure_pilot_bundle(root, expected_manifest_sha256)
        bundle.verify(expected_manifest_sha256)
        return result


def _manifest(case_descriptors: Mapping[str, Mapping[str, object]]) -> dict[str, object]:
    policy = pilot._load_policy()
    preregistration_payload = pilot._regular_bytes(
        pilot.PREREGISTRATION_PATH, "pressure-pilot preregistration"
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
        "cases": [
            case_descriptors[case_id]["manifest_case"]
            for case_id in case_verifier.CASE_IDS
        ],
    }


def publish_pressure_pilot_bundle(
    output_path: str | Path,
    *,
    case_artifact_dirs: Mapping[str, str | Path],
    case_descriptor_sha256: Mapping[str, str],
) -> dict[str, str]:
    """Copy four descriptor-verified raw cases into one immutable aggregate tree."""
    target = Path(os.path.abspath(output_path))
    parent = target.parent
    _require(parent.is_dir() and not parent.is_symlink(), "pressure-pilot output parent is unavailable")
    _require(set(case_artifact_dirs) == set(case_verifier.CASE_IDS), "pressure-pilot source tree set drifted")
    _require(set(case_descriptor_sha256) == set(case_verifier.CASE_IDS), "pressure-pilot descriptor set drifted")
    _require(not _lexists(target), "pressure-pilot output already exists")

    staging = parent / f".{target.name}.staging-{uuid.uuid4()}"
    staging.mkdir(mode=0o700)
    renamed = False
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
                    destination = staging / target_path
                    destination.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
                    _write_exclusive(destination, payload)
                    copied.add(target_path)
            _require(copied == expected_targets, "pressure-pilot copied member closure drifted")
            manifest_payload = _canonical_json_bytes(manifest)
            manifest_sha256 = _sha256(manifest_payload)
            _write_exclusive(staging / pilot.MANIFEST_NAME, manifest_payload)
            _freeze_directories(staging)
            verify_published_pressure_pilot_bundle(staging, manifest_sha256)
            for tree in trees.values():
                tree.require_tree_closure()
            _require(not _lexists(target), "pressure-pilot output appeared during publication")
            _rename_no_replace(staging, target)
            renamed = True
            _fsync_directory(parent)
        verify_published_pressure_pilot_bundle(target, manifest_sha256)
        return {"path": str(target), "manifest_sha256": manifest_sha256}
    except BaseException:
        if not renamed:
            _cleanup_staging(staging)
        raise


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
    parser.add_argument("output_path", type=Path)
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
    receipt = publish_pressure_pilot_bundle(
        args.output_path,
        case_artifact_dirs=_assignments(args.case_artifact_dir, "--case-artifact-dir"),
        case_descriptor_sha256=_assignments(
            args.case_descriptor_sha256, "--case-descriptor-sha256"
        ),
    )
    print(receipt["manifest_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
