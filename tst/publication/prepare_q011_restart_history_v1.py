#!/usr/bin/env python3
"""Prepare an AthenaK history file for a Q011 restart continuation.

The restart payload may be very large.  Only the bounded parameter dump and
fixed-size mesh header are retained in memory; whole-file hashes are computed
incrementally.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
import struct
import sys
import tempfile
from typing import BinaryIO, Sequence


RECORD_TYPE = "q011_restart_history_preparation_receipt_v1"
SCHEMA_VERSION = 1

_PARAMETER_END = b"<par_end>\n"
_RESTART_MESH_HEADER = struct.Struct("<ii9d19i19iddii")
_ATHENA_PARAMETER_SCAN_LIMIT = 41 * 4096
_IO_CHUNK_SIZE = 1024 * 1024
_MAX_HISTORY_LINE_BYTES = 1024 * 1024
_MAX_MESHBLOCKS = 1 << 20
_TIME_ULPS = 32
_HISTORY_BANNER = "# Athena++ history data"
_REAL = re.compile(
    r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
)
_HEADER_FIELD = re.compile(r"\[([1-9][0-9]*)\]=([^\s]+)")


class PreparationError(ValueError):
    """Raised when restart or history input cannot be normalized safely."""


@dataclass(frozen=True)
class _FileIdentity:
    device: int
    inode: int
    size: int
    mtime_ns: int
    ctime_ns: int


@dataclass(frozen=True)
class RestartMetadata:
    meshblock_count: int
    root_level: int
    time: float
    timestep: float
    cycle: int
    original_rank_count: int
    identity: _FileIdentity


@dataclass(frozen=True)
class _HistoryRow:
    values: tuple[float, ...]
    raw: bytes
    sequence: int
    line_number: int

    @property
    def time(self) -> float:
        return self.values[0]


@dataclass(frozen=True)
class _HistoryStage:
    archive_temp: Path
    rollback_temp: Path
    source_identity: _FileIdentity
    source_mode: int
    source_sha256: str
    source_byte_count: int
    labels: tuple[str, ...]
    input_row_count: int
    segment_header_count: int
    retained_rows: tuple[_HistoryRow, ...]
    dropped_row_count: int
    duplicate_row_count: int


@dataclass(frozen=True)
class _StagedFile:
    path: Path
    sha256: str
    byte_count: int


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PreparationError(message)


def _identity(file_stat: os.stat_result) -> _FileIdentity:
    return _FileIdentity(
        device=file_stat.st_dev,
        inode=file_stat.st_ino,
        size=file_stat.st_size,
        mtime_ns=file_stat.st_mtime_ns,
        ctime_ns=file_stat.st_ctime_ns,
    )


def _path_identity(path: Path, *, label: str) -> _FileIdentity:
    try:
        file_stat = path.stat()
    except OSError as error:
        raise PreparationError(f"{label}: cannot stat {path}: {error}") from error
    _require(stat.S_ISREG(file_stat.st_mode), f"{label}: expected a regular file")
    return _identity(file_stat)


def _open_binary(path: Path) -> BinaryIO:
    """Open through one mockable boundary used by bounded-I/O tests."""

    return path.open("rb")


def _read_exact(stream: BinaryIO, size: int, *, label: str) -> bytes:
    chunks: list[bytes] = []
    remaining = size
    while remaining:
        chunk = stream.read(min(remaining, _IO_CHUNK_SIZE))
        if not chunk:
            raise PreparationError(f"{label}: restart mesh header is truncated")
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def read_restart_metadata(path: Path) -> RestartMetadata:
    """Read and validate the fixed AthenaK mesh header without reading the payload."""

    path = Path(path)
    try:
        with _open_binary(path) as stream:
            before_stat = os.fstat(stream.fileno())
            _require(
                stat.S_ISREG(before_stat.st_mode),
                "restart: expected a regular file",
            )
            buffered = bytearray()
            marker_at = -1
            while marker_at < 0:
                if len(buffered) >= _ATHENA_PARAMETER_SCAN_LIMIT:
                    raise PreparationError(
                        "restart: <par_end> is absent from the bounded parameter dump"
                    )
                chunk = stream.read(
                    min(
                        4096,
                        _ATHENA_PARAMETER_SCAN_LIMIT - len(buffered),
                    )
                )
                if not chunk:
                    raise PreparationError("restart: missing <par_end> marker")
                buffered.extend(chunk)
                marker_at = buffered.find(_PARAMETER_END)

            _require(marker_at > 0, "restart: parameter dump is empty")
            _require(
                buffered[marker_at - 1] == ord("\n"),
                "restart: <par_end> is not a standalone parameter line",
            )
            mesh_offset = marker_at + len(_PARAMETER_END)
            mesh_bytes = bytes(buffered[mesh_offset : mesh_offset + _RESTART_MESH_HEADER.size])
            if len(mesh_bytes) < _RESTART_MESH_HEADER.size:
                mesh_bytes += _read_exact(
                    stream,
                    _RESTART_MESH_HEADER.size - len(mesh_bytes),
                    label="restart",
                )
            after_stat = os.fstat(stream.fileno())
    except PreparationError:
        raise
    except OSError as error:
        raise PreparationError(f"restart: cannot read {path}: {error}") from error

    before = _identity(before_stat)
    after = _identity(after_stat)
    _require(before == after, "restart: file changed while reading its mesh header")
    try:
        values = _RESTART_MESH_HEADER.unpack(mesh_bytes)
    except struct.error as error:
        raise PreparationError("restart: invalid mesh header") from error

    meshblock_count = values[0]
    root_level = values[1]
    time, timestep, cycle, original_rank_count = values[-4:]
    valid = (
        0 < meshblock_count <= _MAX_MESHBLOCKS
        and 0 <= root_level <= 30
        and math.isfinite(time)
        and time >= 0.0
        and math.isfinite(timestep)
        and timestep > 0.0
        and cycle >= 0
        and 0 < original_rank_count <= meshblock_count
    )
    _require(valid, "restart: mesh chronology or cardinality is invalid")
    return RestartMetadata(
        meshblock_count=meshblock_count,
        root_level=root_level,
        time=float(time),
        timestep=float(timestep),
        cycle=cycle,
        original_rank_count=original_rank_count,
        identity=before,
    )


def _hash_file(path: Path, *, label: str) -> tuple[str, int, _FileIdentity]:
    digest = hashlib.sha256()
    byte_count = 0
    try:
        with _open_binary(path) as stream:
            before_stat = os.fstat(stream.fileno())
            _require(
                stat.S_ISREG(before_stat.st_mode),
                f"{label}: expected a regular file",
            )
            while True:
                chunk = stream.read(_IO_CHUNK_SIZE)
                if not chunk:
                    break
                digest.update(chunk)
                byte_count += len(chunk)
            after_stat = os.fstat(stream.fileno())
    except PreparationError:
        raise
    except OSError as error:
        raise PreparationError(f"{label}: cannot hash {path}: {error}") from error

    before = _identity(before_stat)
    after = _identity(after_stat)
    _require(before == after, f"{label}: file changed while hashing")
    _require(byte_count == before.size, f"{label}: byte count changed while hashing")
    return digest.hexdigest(), byte_count, before


def _time_tolerance(left: float, right: float) -> float:
    return _TIME_ULPS * max(math.ulp(left), math.ulp(right), math.ulp(1.0))


def _times_close(left: float, right: float) -> bool:
    return abs(left - right) <= _time_tolerance(left, right)


def _parse_header_line(line: str, *, line_number: int) -> tuple[str, ...]:
    matches = list(_HEADER_FIELD.finditer(line))
    _require(matches, f"history: malformed column header at line {line_number}")
    _require(
        re.fullmatch(r"#\s+", line[: matches[0].start()]) is not None,
        f"history: malformed column header at line {line_number}",
    )
    for previous, current in zip(matches, matches[1:]):
        _require(
            line[previous.end() : current.start()].isspace(),
            f"history: malformed column header at line {line_number}",
        )
    _require(
        not line[matches[-1].end() :].strip(),
        f"history: malformed column header at line {line_number}",
    )
    indices = tuple(int(match.group(1)) for match in matches)
    labels = tuple(match.group(2) for match in matches)
    _require(
        indices == tuple(range(1, len(indices) + 1)),
        f"history: noncontiguous column indices at line {line_number}",
    )
    _require(
        len(labels) >= 2 and labels[:2] == ("time", "dt"),
        f"history: first columns must be time and dt at line {line_number}",
    )
    _require(
        len(set(labels)) == len(labels),
        f"history: duplicate column label at line {line_number}",
    )
    return labels


def _new_temp(target: Path) -> tuple[int, Path]:
    try:
        descriptor, name = tempfile.mkstemp(
            prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
        )
    except OSError as error:
        raise PreparationError(f"cannot stage output beside {target}: {error}") from error
    return descriptor, Path(name)


def _stage_history(history: Path, archive: Path, restart_time: float) -> _HistoryStage:
    archive_fd, archive_temp = _new_temp(archive)
    try:
        rollback_fd, rollback_temp = _new_temp(history)
    except Exception:
        os.close(archive_fd)
        archive_temp.unlink(missing_ok=True)
        raise
    source_digest = hashlib.sha256()
    source_byte_count = 0
    labels: tuple[str, ...] | None = None
    rows: list[_HistoryRow] = []
    input_row_count = 0
    segment_header_count = 0
    expecting_columns = False
    segment_row_count = 0
    segment_last_time: float | None = None

    try:
        with os.fdopen(archive_fd, "wb") as archive_stream, os.fdopen(
            rollback_fd, "wb"
        ) as rollback_stream:
            archive_fd = -1
            rollback_fd = -1
            with _open_binary(history) as source:
                before_stat = os.fstat(source.fileno())
                _require(
                    stat.S_ISREG(before_stat.st_mode),
                    "history: expected a regular file",
                )
                source_mode = stat.S_IMODE(before_stat.st_mode)
                os.fchmod(archive_stream.fileno(), source_mode)
                os.fchmod(rollback_stream.fileno(), source_mode)
                line_number = 0
                while True:
                    raw = source.readline(_MAX_HISTORY_LINE_BYTES + 1)
                    if not raw:
                        break
                    line_number += 1
                    _require(
                        len(raw) <= _MAX_HISTORY_LINE_BYTES,
                        f"history: line {line_number} exceeds the safety limit",
                    )
                    archive_stream.write(raw)
                    rollback_stream.write(raw)
                    source_digest.update(raw)
                    source_byte_count += len(raw)
                    _require(
                        raw.endswith(b"\n"),
                        f"history: unterminated line {line_number}",
                    )
                    try:
                        line = raw[:-1].decode("ascii")
                    except UnicodeDecodeError as error:
                        raise PreparationError(
                            f"history: non-ASCII content at line {line_number}"
                        ) from error
                    _require("\r" not in line, f"history: CRLF at line {line_number}")
                    _require(line.strip(), f"history: blank line at line {line_number}")

                    if line == _HISTORY_BANNER:
                        _require(
                            not expecting_columns,
                            f"history: duplicate banner at line {line_number}",
                        )
                        _require(
                            segment_header_count == 0 or segment_row_count > 0,
                            f"history: empty segment before line {line_number}",
                        )
                        expecting_columns = True
                        segment_row_count = 0
                        segment_last_time = None
                        continue

                    if line.startswith("#"):
                        _require(
                            expecting_columns,
                            f"history: column header lacks banner at line {line_number}",
                        )
                        observed = _parse_header_line(line, line_number=line_number)
                        if labels is None:
                            labels = observed
                        else:
                            _require(
                                observed == labels,
                                f"history: column drift at line {line_number}",
                            )
                        segment_header_count += 1
                        expecting_columns = False
                        continue

                    _require(
                        labels is not None and segment_header_count > 0,
                        f"history: data precedes its header at line {line_number}",
                    )
                    _require(
                        not expecting_columns,
                        f"history: data interrupts a header at line {line_number}",
                    )
                    fields = line.split()
                    _require(
                        len(fields) == len(labels),
                        f"history: column count drift at line {line_number}",
                    )
                    values: list[float] = []
                    for field in fields:
                        _require(
                            _REAL.fullmatch(field) is not None,
                            f"history: malformed numeric field at line {line_number}",
                        )
                        value = float(field)
                        _require(
                            math.isfinite(value),
                            f"history: nonfinite value at line {line_number}",
                        )
                        values.append(value)
                    _require(
                        values[0] >= 0.0 and values[1] > 0.0,
                        f"history: invalid time or timestep at line {line_number}",
                    )
                    if segment_last_time is not None:
                        _require(
                            values[0] > segment_last_time
                            or _times_close(values[0], segment_last_time),
                            f"history: nonmonotonic segment at line {line_number}",
                        )
                    segment_last_time = values[0]
                    input_row_count += 1
                    segment_row_count += 1
                    rows.append(
                        _HistoryRow(
                            values=tuple(values),
                            raw=raw,
                            sequence=input_row_count - 1,
                            line_number=line_number,
                        )
                    )
                after_stat = os.fstat(source.fileno())
            archive_stream.flush()
            os.fsync(archive_stream.fileno())
            rollback_stream.flush()
            os.fsync(rollback_stream.fileno())
    except Exception:
        if archive_fd >= 0:
            os.close(archive_fd)
        if rollback_fd >= 0:
            os.close(rollback_fd)
        archive_temp.unlink(missing_ok=True)
        rollback_temp.unlink(missing_ok=True)
        raise
    try:
        before = _identity(before_stat)
        after = _identity(after_stat)
        _require(before == after, "history: file changed while reading")
        _require(
            source_byte_count == before.size,
            "history: byte count changed while reading",
        )
        _require(labels is not None, "history: missing canonical header")
        _require(not expecting_columns, "history: incomplete trailing header")
        _require(segment_row_count > 0, "history: trailing segment has no rows")
        _require(input_row_count > 0, "history: no numeric rows")

        cutoff_rows: list[_HistoryRow] = []
        dropped_row_count = 0
        for row in rows:
            if row.time <= restart_time or _times_close(row.time, restart_time):
                cutoff_rows.append(row)
            else:
                dropped_row_count += 1
        _require(cutoff_rows, "history: no rows exist at or before restart time")

        ordered = sorted(cutoff_rows, key=lambda row: (row.time, row.sequence))
        retained: list[_HistoryRow] = []
        duplicate_row_count = 0
        cursor = 0
        while cursor < len(ordered):
            group = [ordered[cursor]]
            anchor_time = ordered[cursor].time
            cursor += 1
            while cursor < len(ordered) and _times_close(
                ordered[cursor].time, anchor_time
            ):
                group.append(ordered[cursor])
                cursor += 1
            selected = min(group, key=lambda row: row.sequence)
            for duplicate in group:
                if duplicate is selected:
                    continue
                _require(
                    duplicate.values[2:] == selected.values[2:],
                    "history: duplicate time has conflicting state at "
                    f"lines {selected.line_number} and {duplicate.line_number}",
                )
            retained.append(selected)
            duplicate_row_count += len(group) - 1

        for previous, current in zip(retained, retained[1:]):
            _require(
                current.time > previous.time,
                "history: normalized output would not be strictly monotonic",
            )
        _require(
            retained[-1].time <= restart_time
            or _times_close(retained[-1].time, restart_time),
            "history: normalized output exceeds restart time",
        )
        return _HistoryStage(
            archive_temp=archive_temp,
            rollback_temp=rollback_temp,
            source_identity=before,
            source_mode=source_mode,
            source_sha256=source_digest.hexdigest(),
            source_byte_count=source_byte_count,
            labels=labels,
            input_row_count=input_row_count,
            segment_header_count=segment_header_count,
            retained_rows=tuple(retained),
            dropped_row_count=dropped_row_count,
            duplicate_row_count=duplicate_row_count,
        )
    except Exception:
        archive_temp.unlink(missing_ok=True)
        rollback_temp.unlink(missing_ok=True)
        raise


def _canonical_history_header(labels: Sequence[str]) -> bytes:
    columns = " ".join(
        f"[{index}]={label}" for index, label in enumerate(labels, start=1)
    )
    return f"{_HISTORY_BANNER}\n#  {columns}\n".encode("ascii")


def _stage_normalized_history(history: Path, stage: _HistoryStage) -> _StagedFile:
    descriptor, temp_path = _new_temp(history)
    digest = hashlib.sha256()
    byte_count = 0
    try:
        with os.fdopen(descriptor, "wb") as output:
            descriptor = -1
            os.fchmod(output.fileno(), stage.source_mode)
            chunks = [_canonical_history_header(stage.labels)]
            chunks.extend(row.raw for row in stage.retained_rows)
            for chunk in chunks:
                output.write(chunk)
                digest.update(chunk)
                byte_count += len(chunk)
            output.flush()
            os.fsync(output.fileno())
    except Exception:
        if descriptor >= 0:
            os.close(descriptor)
        temp_path.unlink(missing_ok=True)
        raise
    return _StagedFile(temp_path, digest.hexdigest(), byte_count)


def _canonical_json(value: object) -> bytes:
    try:
        return (
            json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
            + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise PreparationError("receipt: value is not canonical JSON") from error


def _stage_payload(target: Path, payload: bytes, *, mode: int = 0o600) -> _StagedFile:
    descriptor, temp_path = _new_temp(target)
    try:
        with os.fdopen(descriptor, "wb") as output:
            descriptor = -1
            os.fchmod(output.fileno(), mode)
            output.write(payload)
            output.flush()
            os.fsync(output.fileno())
    except Exception:
        if descriptor >= 0:
            os.close(descriptor)
        temp_path.unlink(missing_ok=True)
        raise
    return _StagedFile(temp_path, hashlib.sha256(payload).hexdigest(), len(payload))


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _atomic_replace(source: Path, destination: Path) -> None:
    os.replace(source, destination)


def _resolved(path: Path) -> str:
    return str(path.resolve(strict=False))


def _lexists(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def _validate_paths(
    restart: Path, history: Path, archive: Path, receipt: Path
) -> None:
    resolved = [_resolved(path) for path in (restart, history, archive, receipt)]
    _require(len(set(resolved)) == 4, "restart, history, archive, and receipt must differ")
    _require(not history.is_symlink(), "history: symbolic links are not supported")
    _require(archive.parent.is_dir(), "archive: parent directory does not exist")
    _require(receipt.parent.is_dir(), "receipt: parent directory does not exist")
    _require(not _lexists(archive), "archive: destination already exists")
    _require(not _lexists(receipt), "receipt: destination already exists")


def prepare_restart_history(
    *, restart: Path, history: Path, archive: Path, receipt: Path
) -> dict[str, object]:
    """Validate, archive, and atomically normalize one Athena history file."""

    restart = Path(restart)
    history = Path(history)
    archive = Path(archive)
    receipt = Path(receipt)
    _validate_paths(restart, history, archive, receipt)

    metadata = read_restart_metadata(restart)
    restart_sha256, restart_byte_count, restart_identity = _hash_file(
        restart, label="restart"
    )
    _require(
        metadata.identity == restart_identity,
        "restart: file changed between mesh-header read and hashing",
    )

    history_stage: _HistoryStage | None = None
    normalized_stage: _StagedFile | None = None
    receipt_stage: _StagedFile | None = None
    archive_published = False
    history_published = False
    receipt_published = False
    try:
        history_stage = _stage_history(history, archive, metadata.time)
        normalized_stage = _stage_normalized_history(history, history_stage)
        record: dict[str, object] = {
            "archive": {
                "byte_count": history_stage.source_byte_count,
                "path": _resolved(archive),
                "sha256": history_stage.source_sha256,
            },
            "history": {
                "column_count": len(history_stage.labels),
                "dropped_row_count": history_stage.dropped_row_count,
                "duplicate_row_count": history_stage.duplicate_row_count,
                "input_byte_count": history_stage.source_byte_count,
                "input_row_count": history_stage.input_row_count,
                "input_sha256": history_stage.source_sha256,
                "normalized_byte_count": normalized_stage.byte_count,
                "normalized_sha256": normalized_stage.sha256,
                "path": _resolved(history),
                "retained_row_count": len(history_stage.retained_rows),
                "segment_header_count": history_stage.segment_header_count,
            },
            "record_type": RECORD_TYPE,
            "restart": {
                "byte_count": restart_byte_count,
                "cycle": metadata.cycle,
                "meshblock_count": metadata.meshblock_count,
                "original_rank_count": metadata.original_rank_count,
                "path": _resolved(restart),
                "sha256": restart_sha256,
                "time": metadata.time,
                "timestep": metadata.timestep,
            },
            "schema_version": SCHEMA_VERSION,
        }
        receipt_stage = _stage_payload(receipt, _canonical_json(record))

        _require(
            _path_identity(restart, label="restart") == restart_identity,
            "restart: file changed before commit",
        )
        _require(
            _path_identity(history, label="history") == history_stage.source_identity,
            "history: file changed before commit",
        )
        _require(
            not _lexists(archive), "archive: destination appeared before commit"
        )
        _require(
            not _lexists(receipt), "receipt: destination appeared before commit"
        )

        _atomic_replace(history_stage.archive_temp, archive)
        archive_published = True
        _fsync_directory(archive.parent)
        _atomic_replace(normalized_stage.path, history)
        history_published = True
        _fsync_directory(history.parent)
        _atomic_replace(receipt_stage.path, receipt)
        receipt_published = True
        _fsync_directory(receipt.parent)
        return record
    except Exception as error:
        rollback_error: OSError | None = None
        try:
            if receipt_published or _lexists(receipt):
                receipt.unlink(missing_ok=True)
                _fsync_directory(receipt.parent)
            if history_published and history_stage is not None:
                os.replace(history_stage.rollback_temp, history)
                _fsync_directory(history.parent)
            if archive_published and _lexists(archive):
                archive.unlink()
                _fsync_directory(archive.parent)
        except OSError as caught:
            rollback_error = caught
        if rollback_error is not None:
            raise PreparationError(
                f"publication failed ({error}); rollback also failed: {rollback_error}"
            ) from error
        raise
    finally:
        temporary_paths = []
        if history_stage is not None:
            temporary_paths.append(history_stage.archive_temp)
            temporary_paths.append(history_stage.rollback_temp)
        if normalized_stage is not None:
            temporary_paths.append(normalized_stage.path)
        if receipt_stage is not None:
            temporary_paths.append(receipt_stage.path)
        for temporary in temporary_paths:
            temporary.unlink(missing_ok=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare an AthenaK history file for a Q011 restart continuation."
    )
    parser.add_argument("--restart", required=True, type=Path)
    parser.add_argument("--history", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _parser()
    arguments = parser.parse_args(argv)
    try:
        record = prepare_restart_history(
            restart=arguments.restart,
            history=arguments.history,
            archive=arguments.archive,
            receipt=arguments.receipt,
        )
    except (PreparationError, OSError) as error:
        parser.exit(1, f"{parser.prog}: error: {error}\n")
    sys.stdout.buffer.write(_canonical_json(record))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
