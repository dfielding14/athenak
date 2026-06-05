#!/usr/bin/env python3
"""Consume only the exact published historical Q-011 pressure-pilot evidence."""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import stat
from typing import Any, Iterator

if __package__:
    from . import analyze_q011_section54_pressure_pilot as _pressure_pilot_analyzer
    from .frontier_control_plane import (
        q011_pressure_review_packet_verifier as _pressure_review_packet_verifier,
    )
else:
    import analyze_q011_section54_pressure_pilot as _pressure_pilot_analyzer
    from frontier_control_plane import (
        q011_pressure_review_packet_verifier as _pressure_review_packet_verifier,
    )


AUTHORIZED_PRODUCTION_PIC_ROOT = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/PIC"
)
AUTHORIZED_PRODUCTION_PUBLICATION_ROOT = AUTHORIZED_PRODUCTION_PIC_ROOT / "publication"
AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING = {
    "path": str(
        AUTHORIZED_PRODUCTION_PUBLICATION_ROOT
        / "q011_section54_pressure_pilot_review_packet_receipt.json"
    ),
    "sha256": "3f20d3d26a479aa508439f9d038ec6510643bf407aa081fae22959a57571de5d",
}
AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING = {
    "path": str(
        AUTHORIZED_PRODUCTION_PUBLICATION_ROOT
        / "q011_section54_pressure_pilot_bundle_receipt.json"
    ),
    "sha256": "9117b3dbc7573187b2d080568e69bdbbee0642f2a965aa543273ab3ea3d67be9",
}
AUTHORIZED_PRODUCTION_AGGREGATE_BUNDLE_BINDING = {
    "path": str(
        AUTHORIZED_PRODUCTION_PUBLICATION_ROOT
        / "q011_section54_pressure_pilot_bundle"
    ),
    "manifest_sha256": "7b3fb8de9e4dc8d6d8b2320a2c2aeaadf8bc6f076dab6b3e98ef8a415abdbf55",
}
AUTHORIZED_PRODUCTION_AGGREGATE_ANALYSIS_BINDING = {
    "path": str(
        AUTHORIZED_PRODUCTION_PUBLICATION_ROOT
        / "q011_section54_pressure_pilot_analysis.json"
    ),
    "sha256": "d55b4c2020716df899c86dfe5a9018d48194d63590616ff60541067243daacb7",
}
AGGREGATE_MANIFEST_NAME = "pressure_pilot_manifest.json"
MAX_RETAINED_FILE_BYTES = 128 * 1024 * 1024
MAX_JSON_BYTES = 8 * 1024 * 1024
MAX_DIRECTORY_ENTRIES = 256
MAX_TREE_ENTRIES = 1024
MAX_TREE_DEPTH = 16
_READ_CHUNK_BYTES = 1024 * 1024
_STABLE_STAT_FIELDS = (
    "st_dev",
    "st_ino",
    "st_mode",
    "st_nlink",
    "st_size",
    "st_mtime_ns",
    "st_ctime_ns",
)
_IDENTITY_STAT_FIELDS = ("st_dev", "st_ino", "st_mode")


class HistoricalPressurePilotConsumerError(ValueError):
    """Raised when the exact historical production evidence fails closed."""


def _fail(message: str) -> None:
    raise HistoricalPressurePilotConsumerError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical_json_bytes(value: object) -> bytes:
    try:
        payload = (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (RecursionError, TypeError, ValueError) as exc:
        _fail(f"recomputed aggregate analysis is not canonical JSON: {exc}")
    _require(
        len(payload) <= MAX_JSON_BYTES,
        "recomputed aggregate analysis exceeds the JSON size limit",
    )
    return payload


def _same_json_value(left: object, right: object) -> bool:
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return set(left) == set(right) and all(
            _same_json_value(left[key], right[key]) for key in left
        )
    if type(left) is list:
        return len(left) == len(right) and all(
            _same_json_value(a, b) for a, b in zip(left, right)
        )
    return left == right


def _stable_stat(metadata: os.stat_result) -> tuple[int, ...]:
    return tuple(getattr(metadata, field) for field in _STABLE_STAT_FIELDS)


def _identity_stat(metadata: os.stat_result) -> tuple[int, ...]:
    return tuple(getattr(metadata, field) for field in _IDENTITY_STAT_FIELDS)


def _open_flags(*, directory: bool = False) -> int:
    flags = os.O_RDONLY
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    flags |= getattr(os, "O_NONBLOCK", 0)
    if directory:
        flags |= getattr(os, "O_DIRECTORY", 0)
    return flags


def _canonical_absolute_path(value: Path, label: str) -> Path:
    _require(value.is_absolute(), f"{label} must be absolute")
    lexical = Path(os.path.abspath(os.fspath(value)))
    _require(value == lexical, f"{label} must be lexically canonical")
    return value


def _require_readonly_directory(metadata: os.stat_result, label: str) -> None:
    _require(stat.S_ISDIR(metadata.st_mode), f"{label} must be a directory")
    _require(metadata.st_mode & 0o222 == 0, f"{label} must be read-only")


def _require_readonly_file(
    metadata: os.stat_result,
    label: str,
    max_bytes: int,
) -> None:
    _require(stat.S_ISREG(metadata.st_mode), f"{label} must be a regular file")
    _require(metadata.st_nlink == 1, f"{label} must have exactly one hard link")
    _require(metadata.st_mode & 0o222 == 0, f"{label} must be read-only")
    _require(metadata.st_size <= max_bytes, f"{label} exceeds the size limit")


def _member_size_limit(relative: str) -> int:
    if relative == AGGREGATE_MANIFEST_NAME or relative.endswith((".json", ".manifest")):
        return MAX_JSON_BYTES
    return MAX_RETAINED_FILE_BYTES


@contextlib.contextmanager
def _retained_child_directory(
    parent_fd: int,
    name: str,
    label: str,
    *,
    immutable: bool,
) -> Iterator[tuple[int, os.stat_result]]:
    _require(name not in {"", ".", ".."} and "/" not in name, f"{label} has an unsafe name")
    fd = -1
    try:
        observed = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        _require(stat.S_ISDIR(observed.st_mode), f"{label} must be a directory")
        if immutable:
            _require_readonly_directory(observed, label)
        fd = os.open(name, _open_flags(directory=True), dir_fd=parent_fd)
        opened = os.fstat(fd)
        expected = _stable_stat(observed) if immutable else _identity_stat(observed)
        actual = _stable_stat(opened) if immutable else _identity_stat(opened)
        _require(actual == expected, f"{label} changed while being opened")
        yield fd, opened
        after = os.fstat(fd)
        current = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        after_value = _stable_stat(after) if immutable else _identity_stat(after)
        current_value = _stable_stat(current) if immutable else _identity_stat(current)
        _require(
            after_value == actual and current_value == actual,
            f"{label} changed during consumption",
        )
    except HistoricalPressurePilotConsumerError:
        raise
    except OSError as exc:
        _fail(f"{label} cannot be retained without following links: {exc}")
    finally:
        if fd >= 0:
            os.close(fd)


@contextlib.contextmanager
def _retained_absolute_directory(
    path: Path,
    label: str,
) -> Iterator[tuple[int, os.stat_result]]:
    canonical = _canonical_absolute_path(path, label)
    root_fd = -1
    try:
        root_fd = os.open("/", _open_flags(directory=True))
        root_identity = _identity_stat(os.fstat(root_fd))
        with contextlib.ExitStack() as stack:
            current_fd = root_fd
            current_metadata = os.fstat(root_fd)
            traversed = Path("/")
            for part in canonical.parts[1:]:
                traversed /= part
                current_fd, current_metadata = stack.enter_context(
                    _retained_child_directory(
                        current_fd,
                        part,
                        f"{label} component {str(traversed)!r}",
                        immutable=False,
                    )
                )
            yield current_fd, current_metadata
            _require(
                _identity_stat(os.fstat(root_fd)) == root_identity,
                "filesystem root changed during consumption",
            )
    except HistoricalPressurePilotConsumerError:
        raise
    except OSError as exc:
        _fail(f"{label} cannot be opened by descriptor-relative ancestry: {exc}")
    finally:
        if root_fd >= 0:
            os.close(root_fd)


def _bounded_directory_names(directory_fd: int, relative: str) -> list[str]:
    names: list[str] = []
    try:
        with os.scandir(directory_fd) as entries:
            for entry in entries:
                _require(
                    len(names) < MAX_DIRECTORY_ENTRIES,
                    f"aggregate bundle directory {relative or '.'!r} exceeds the entry limit",
                )
                names.append(entry.name)
    except HistoricalPressurePilotConsumerError:
        raise
    except OSError as exc:
        _fail(f"aggregate bundle directory {relative or '.'!r} cannot be listed: {exc}")
    _require(bool(names), f"aggregate bundle directory {relative or '.'!r} must not be empty")
    return sorted(names)


def _scan_bundle_directory(
    directory_fd: int,
    relative: str,
    directories: dict[str, tuple[int, ...]],
    files: dict[str, tuple[int, ...]],
    entry_count: list[int],
    depth: int,
) -> None:
    _require(depth <= MAX_TREE_DEPTH, "aggregate bundle exceeds the tree-depth limit")
    for name in _bounded_directory_names(directory_fd, relative):
        _require(
            entry_count[0] < MAX_TREE_ENTRIES,
            "aggregate bundle exceeds the tree-entry limit",
        )
        entry_count[0] += 1
        _require(
            name not in {"", ".", ".."} and "/" not in name,
            "aggregate bundle contains an unsafe entry name",
        )
        child_relative = f"{relative}/{name}" if relative else name
        try:
            observed = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        except OSError as exc:
            _fail(f"aggregate bundle member {child_relative!r} cannot be inspected: {exc}")
        if stat.S_ISREG(observed.st_mode):
            limit = _member_size_limit(child_relative)
            _require_readonly_file(observed, f"aggregate bundle member {child_relative!r}", limit)
            fd = -1
            try:
                fd = os.open(name, _open_flags(), dir_fd=directory_fd)
                opened = os.fstat(fd)
                _require(
                    _stable_stat(opened) == _stable_stat(observed),
                    f"aggregate bundle member {child_relative!r} changed while being opened",
                )
            except HistoricalPressurePilotConsumerError:
                raise
            except OSError as exc:
                _fail(
                    f"aggregate bundle member {child_relative!r} cannot be opened "
                    f"without following links: {exc}"
                )
            finally:
                if fd >= 0:
                    os.close(fd)
            files[child_relative] = _stable_stat(observed)
        elif stat.S_ISDIR(observed.st_mode):
            _require_readonly_directory(
                observed,
                f"aggregate bundle directory {child_relative!r}",
            )
            with _retained_child_directory(
                directory_fd,
                name,
                f"aggregate bundle directory {child_relative!r}",
                immutable=True,
            ) as (child_fd, child_metadata):
                directories[child_relative] = _stable_stat(child_metadata)
                _scan_bundle_directory(
                    child_fd,
                    child_relative,
                    directories,
                    files,
                    entry_count,
                    depth + 1,
                )
        else:
            _fail(f"aggregate bundle member {child_relative!r} has an unsupported file type")


def _scan_bundle(root_fd: int, root_metadata: os.stat_result) -> dict[str, object]:
    current = os.fstat(root_fd)
    _require_readonly_directory(current, "aggregate bundle root")
    _require(
        _stable_stat(current) == _stable_stat(root_metadata),
        "aggregate bundle root changed before scan",
    )
    directories: dict[str, tuple[int, ...]] = {}
    files: dict[str, tuple[int, ...]] = {}
    _scan_bundle_directory(root_fd, "", directories, files, [0], 0)
    _require(
        _stable_stat(os.fstat(root_fd)) == _stable_stat(root_metadata),
        "aggregate bundle root changed during scan",
    )
    return {
        "root": _stable_stat(root_metadata),
        "directories": directories,
        "files": files,
    }


def _safe_relative_parts(relative: str) -> tuple[str, ...]:
    path = PurePosixPath(relative)
    _require(
        type(relative) is str
        and bool(relative)
        and not path.is_absolute()
        and path.as_posix() == relative
        and len(path.parts) <= MAX_TREE_DEPTH
        and all(part not in {"", ".", ".."} for part in path.parts),
        f"aggregate bundle member path is unsafe: {relative!r}",
    )
    return path.parts


def _read_file_at(
    parent_fd: int,
    name: str,
    label: str,
    *,
    max_bytes: int,
    expected: tuple[int, ...] | None = None,
) -> tuple[bytes, tuple[int, ...]]:
    _require(name not in {"", ".", ".."} and "/" not in name, f"{label} has an unsafe name")
    fd = -1
    try:
        observed = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        _require_readonly_file(observed, label, max_bytes)
        observed_stable = _stable_stat(observed)
        if expected is not None:
            _require(observed_stable == expected, f"{label} changed before stable read")
        fd = os.open(name, _open_flags(), dir_fd=parent_fd)
        opened = os.fstat(fd)
        _require(
            _stable_stat(opened) == observed_stable,
            f"{label} changed while being opened",
        )
        payload = bytearray()
        while True:
            chunk = os.read(
                fd,
                min(_READ_CHUNK_BYTES, max_bytes + 1 - len(payload)),
            )
            if not chunk:
                break
            payload.extend(chunk)
            _require(len(payload) <= max_bytes, f"{label} exceeds the size limit")
        after = os.fstat(fd)
        current = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        _require(
            _stable_stat(after) == observed_stable
            and _stable_stat(current) == observed_stable
            and len(payload) == after.st_size,
            f"{label} changed during stable read",
        )
        return bytes(payload), observed_stable
    except HistoricalPressurePilotConsumerError:
        raise
    except OSError as exc:
        _fail(f"{label} cannot be read without following links: {exc}")
    finally:
        if fd >= 0:
            os.close(fd)


class _ImmutableBundleReader:
    def __init__(self, root_fd: int, scan: dict[str, object]) -> None:
        self._root_fd = root_fd
        self._directories = scan["directories"]
        self._files = scan["files"]
        _require(type(self._directories) is dict, "internal directory scan is malformed")
        _require(type(self._files) is dict, "internal file scan is malformed")

    def file_paths(self) -> set[str]:
        return set(self._files)

    def read(self, relative: str) -> bytes:
        parts = _safe_relative_parts(relative)
        _require(relative in self._files, f"aggregate bundle member is outside scan closure: {relative}")
        current_fd = os.dup(self._root_fd)
        try:
            traversed: list[str] = []
            for part in parts[:-1]:
                traversed.append(part)
                current_relative = "/".join(traversed)
                expected = self._directories.get(current_relative)
                _require(
                    type(expected) is tuple,
                    f"aggregate bundle parent is outside scan closure: {current_relative}",
                )
                observed = os.stat(part, dir_fd=current_fd, follow_symlinks=False)
                _require(
                    _stable_stat(observed) == expected,
                    f"aggregate bundle directory changed before read: {current_relative}",
                )
                child_fd = -1
                try:
                    child_fd = os.open(
                        part,
                        _open_flags(directory=True),
                        dir_fd=current_fd,
                    )
                    opened = os.fstat(child_fd)
                    _require(
                        _stable_stat(opened) == expected,
                        f"aggregate bundle directory changed while being opened: {current_relative}",
                    )
                    os.close(current_fd)
                    current_fd = child_fd
                    child_fd = -1
                finally:
                    if child_fd >= 0:
                        os.close(child_fd)
            payload, _metadata = _read_file_at(
                current_fd,
                parts[-1],
                f"aggregate bundle member {relative!r}",
                max_bytes=_member_size_limit(relative),
                expected=self._files[relative],
            )
            return payload
        except HistoricalPressurePilotConsumerError:
            raise
        except OSError as exc:
            _fail(f"aggregate bundle member {relative!r} cannot be retained: {exc}")
        finally:
            os.close(current_fd)


def _verify_exact_packet() -> dict[str, object]:
    try:
        verified = (
            _pressure_review_packet_verifier.consume_published_pressure_pilot_review_packet(
                AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING["path"],
                aggregate_receipt_binding=AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING,
                authorized_pic_root=AUTHORIZED_PRODUCTION_PIC_ROOT,
            )
        )
    except _pressure_review_packet_verifier.PressureReviewPacketVerificationError as exc:
        raise HistoricalPressurePilotConsumerError(
            "exact production pressure-review packet failed immutable verification"
        ) from exc
    _require(
        type(verified) is dict
        and set(verified)
        == {
            "receipt_binding",
            "aggregate_receipt_binding",
            "packet_receipt",
            "aggregate_receipt",
            "aggregate_bundle",
            "aggregate_analysis",
            "source_bindings",
            "inventory",
        },
        "pressure-review packet verifier result schema drifted",
    )
    _require(
        _same_json_value(
            verified["receipt_binding"],
            AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING,
        ),
        "pressure-review packet verifier returned a different packet receipt",
    )
    _require(
        _same_json_value(
            verified["aggregate_receipt_binding"],
            AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING,
        ),
        "pressure-review packet verifier returned a different aggregate receipt",
    )
    _require(
        _same_json_value(
            verified["aggregate_bundle"],
            AUTHORIZED_PRODUCTION_AGGREGATE_BUNDLE_BINDING,
        ),
        "pressure-review packet verifier returned a different aggregate bundle",
    )
    _require(
        _same_json_value(
            verified["aggregate_analysis"],
            AUTHORIZED_PRODUCTION_AGGREGATE_ANALYSIS_BINDING,
        ),
        "pressure-review packet verifier returned a different aggregate analysis",
    )
    return verified


def consume_exact_historical_production_pressure_pilot() -> dict[str, str]:
    """Recompute and consume only the exact historical production pilot pair."""

    try:
        first_packet = _verify_exact_packet()
        bundle_path = Path(AUTHORIZED_PRODUCTION_AGGREGATE_BUNDLE_BINDING["path"])
        analysis_path = Path(AUTHORIZED_PRODUCTION_AGGREGATE_ANALYSIS_BINDING["path"])
        _require(
            bundle_path.parent == AUTHORIZED_PRODUCTION_PUBLICATION_ROOT
            and analysis_path.parent == AUTHORIZED_PRODUCTION_PUBLICATION_ROOT,
            "exact production aggregate artifacts are outside the publication root",
        )
        with _retained_absolute_directory(
            AUTHORIZED_PRODUCTION_PIC_ROOT,
            "authorized production PIC root",
        ) as (pic_root_fd, _pic_root_metadata), _retained_child_directory(
            pic_root_fd,
            "publication",
            "authorized production publication root",
            immutable=False,
        ) as (publication_fd, _publication_metadata), _retained_child_directory(
            publication_fd,
            bundle_path.name,
            "exact production aggregate bundle",
            immutable=True,
        ) as (bundle_fd, bundle_metadata):
            first_scan = _scan_bundle(bundle_fd, bundle_metadata)
            reader = _ImmutableBundleReader(bundle_fd, first_scan)
            manifest_payload = reader.read(AGGREGATE_MANIFEST_NAME)
            manifest_sha256 = AUTHORIZED_PRODUCTION_AGGREGATE_BUNDLE_BINDING[
                "manifest_sha256"
            ]
            _require(
                _sha256(manifest_payload) == manifest_sha256,
                "exact production aggregate manifest hash drifted",
            )
            recomputed = (
                _pressure_pilot_analyzer.analyze_exact_historical_production_pressure_pilot_bundle(
                    bundle_path,
                    manifest_sha256,
                    authorized_publication_root=AUTHORIZED_PRODUCTION_PUBLICATION_ROOT,
                    member_reader=reader.read,
                    actual_files=reader.file_paths(),
                )
            )
            _require(type(recomputed) is dict, "recomputed aggregate analysis must be an object")
            status = recomputed.get("status")
            _require(type(status) is str and bool(status), "recomputed aggregate status is malformed")
            recomputed_payload = _canonical_json_bytes(recomputed)
            retained_payload, retained_metadata = _read_file_at(
                publication_fd,
                analysis_path.name,
                "exact production retained aggregate analysis",
                max_bytes=MAX_JSON_BYTES,
            )
            analysis_sha256 = AUTHORIZED_PRODUCTION_AGGREGATE_ANALYSIS_BINDING["sha256"]
            _require(
                _sha256(retained_payload) == analysis_sha256,
                "exact production retained aggregate analysis hash drifted",
            )
            _require(
                recomputed_payload == retained_payload,
                "recomputed aggregate analysis differs from retained analysis bytes",
            )
            second_scan = _scan_bundle(bundle_fd, bundle_metadata)
            _require(
                first_scan == second_scan,
                "exact production aggregate bundle changed during recomputation",
            )
            _require(
                reader.read(AGGREGATE_MANIFEST_NAME) == manifest_payload,
                "exact production aggregate manifest changed during recomputation",
            )
            final_retained_payload, _final_metadata = _read_file_at(
                publication_fd,
                analysis_path.name,
                "exact production retained aggregate analysis",
                max_bytes=MAX_JSON_BYTES,
                expected=retained_metadata,
            )
            _require(
                final_retained_payload == retained_payload,
                "exact production retained aggregate analysis changed during recomputation",
            )
            second_packet = _verify_exact_packet()
            _require(
                _same_json_value(first_packet, second_packet),
                "exact production packet verification changed during recomputation",
            )
            return {
                "packet_receipt_sha256": AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING[
                    "sha256"
                ],
                "aggregate_receipt_sha256": AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING[
                    "sha256"
                ],
                "manifest_sha256": manifest_sha256,
                "analysis_result_sha256": analysis_sha256,
                "status": status,
            }
    except HistoricalPressurePilotConsumerError:
        raise
    except (AttributeError, OSError, RecursionError, TypeError, ValueError) as exc:
        raise HistoricalPressurePilotConsumerError(
            f"exact historical production pressure-pilot consumption failed closed: {exc}"
        ) from exc


__all__ = [
    "HistoricalPressurePilotConsumerError",
    "consume_exact_historical_production_pressure_pilot",
]
