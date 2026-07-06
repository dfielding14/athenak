#!/usr/bin/env python3
"""Merge AthenaK rich tracked-particle shards into particle-major HDF5.

The source ``trk`` files are ownership shards: records for one tracked particle
can move between files as particle ownership moves between MPI ranks.  This tool
groups records by the rich ``output_tag`` and writes each trajectory contiguously
in a dense HDF5 layout.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import h5py
from mpi4py import MPI
import numpy as np


TRACK_HEADER_MARKER = b"# AthenaK tracked particle data at time="
COMPACT_FILE_MAGIC = b"AKTRK2F\0"
COMPACT_FRAME_MAGIC = b"AKTRK2R\0"
COMPACT_VERSION = 2
COMPACT_PROLOGUE = struct.Struct("<8sHHIqqIIiiiiiII")
COMPACT_FRAME = struct.Struct("<8sHHIqdII")
RICH_FIELDS = (
    "tag", "time", "x", "y", "z", "vx", "vy", "vz",
    "bx", "by", "bz", "k1", "k2", "k3", "db1", "db2", "db3", "jmag",
)
VALUE_FIELDS = RICH_FIELDS[2:]
N_SOURCE_FIELDS = len(RICH_FIELDS)
N_VALUE_FIELDS = len(VALUE_FIELDS)
FORMAT_NAME = "athenak_rich_trk_merged_v1"

TEMP_RECORD_DTYPE = np.dtype([
    ("tag", "<i8"),
    ("cycle", "<i8"),
    ("time", "<f8"),
    ("values", "<f4", (N_VALUE_FIELDS,)),
])

FRAME_META_DTYPE = np.dtype([
    ("cycle", "<i8"),
    ("time_min", "<f8"),
    ("time_max", "<f8"),
    ("record_count", "<i8"),
    ("header_count", "<i8"),
])

PARTICLE_DTYPE = np.dtype([
    ("output_tag", "<i8"),
    ("species", "<i4"),
    ("track_tag", "<i8"),
    ("row", "<i8"),
    ("count", "<i8"),
    ("first_time", "<f8"),
    ("last_time", "<f8"),
])

KEY_PATTERN = re.compile(r"([A-Za-z_]+)=\s*([^ \t\n]+)")
TIME_PATTERN = re.compile(r"# AthenaK tracked particle data at time=\s*([0-9.eE+-]+)")


@dataclass(frozen=True)
class FileEntry:
    path: str
    size: int


@dataclass
class Header:
    time: float
    cycle: int
    trk_format: str
    ntracked: int
    ntrack_per_species: int
    track_per_species: bool
    record_count: int
    nfields: int
    fields: Tuple[str, ...]
    layout: str


@dataclass
class RunMeta:
    trk_format: str
    ntracked: int
    ntrack_per_species: int
    track_per_species: bool
    nfields: int
    fields: Tuple[str, ...]
    layout: str

    @property
    def nspecies(self) -> int:
        if self.track_per_species:
            return self.ntracked // self.ntrack_per_species
        return 1


class MergeError(RuntimeError):
    pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge AthenaK rich trk shards into particle-major HDF5.")
    parser.add_argument("--run-dir", required=True, type=Path,
                        help="AthenaK run directory containing trk/.")
    parser.add_argument("--output", type=Path,
                        help="Output merged HDF5 path.")
    parser.add_argument("--tmp-dir", required=True, type=Path,
                        help="Temporary directory for MPI owner shards.")
    parser.add_argument("--layout", choices=("auto", "rank", "node", "shared"),
                        default="auto")
    parser.add_argument("--require-complete", dest="require_complete",
                        action="store_true", default=True,
                        help="Require every particle at every global frame.")
    parser.add_argument("--allow-missing", action="store_true",
                        help="Reserved for future sparse output support.")
    parser.add_argument("--keep-cycle", action="store_true", default=True,
                        help="Kept for CLI compatibility; cycles are always written.")
    parser.add_argument("--compression", choices=("none", "lzf", "gzip"),
                        default="none")
    parser.add_argument("--max-files-per-rank", type=int, default=0,
                        help="Testing limiter applied after file assignment.")
    parser.add_argument("--exchange-rows", type=int, default=250000,
                        help="Approximate records buffered before MPI exchange.")
    parser.add_argument("--assembly-rows", type=int, default=4,
                        help="Particle rows copied per HDF5 write during assembly.")
    parser.add_argument("--dry-run", action="store_true",
                        help="With --validate-only, validate fully but skip HDF5.")
    parser.add_argument("--validate-only", action="store_true",
                        help="Parse, redistribute, and sort/validate without HDF5.")
    parser.add_argument("--delete-temp", action="store_true",
                        help="Remove temp files after successful HDF5 validation.")
    parser.add_argument("--no-payload-finite-check", action="store_true",
                        help="Skip finite checks on rich value fields.")
    return parser.parse_args(argv)


def rank_print(comm: MPI.Comm, *items: object) -> None:
    if comm.Get_rank() == 0:
        print(*items, flush=True)


def fail_all(comm: MPI.Comm, message: str, code: int = 2) -> None:
    rank = comm.Get_rank()
    print(f"[rank {rank}] ERROR: {message}", file=sys.stderr, flush=True)
    comm.Abort(code)


def discover_files(run_dir: Path, layout: str) -> Tuple[str, List[FileEntry]]:
    trk_dir = run_dir / "trk"
    if not trk_dir.is_dir():
        raise MergeError(f"missing trk directory: {trk_dir}")

    patterns: List[Tuple[str, str]]
    if layout == "auto":
        patterns = [
            ("rank", "rank_*/*.trk"),
            ("node", "node_*/*.trk"),
            ("shared", "*.trk"),
        ]
    elif layout == "rank":
        patterns = [("rank", "rank_*/*.trk")]
    elif layout == "node":
        patterns = [("node", "node_*/*.trk")]
    else:
        patterns = [("shared", "*.trk")]

    for found_layout, pattern in patterns:
        paths = sorted(trk_dir.glob(pattern))
        if paths:
            entries = [FileEntry(str(path), path.stat().st_size) for path in paths]
            return found_layout, entries
    raise MergeError(f"no .trk files found under {trk_dir} for layout={layout}")


def assign_files(entries: Sequence[FileEntry], size: int) -> List[List[FileEntry]]:
    assigned: List[List[FileEntry]] = [[] for _ in range(size)]
    loads = [0] * size
    for entry in sorted(entries, key=lambda item: item.size, reverse=True):
        dest = min(range(size), key=lambda idx: loads[idx])
        assigned[dest].append(entry)
        loads[dest] += entry.size
    return assigned


def parse_header(header_lines: Sequence[str], path: Path) -> Header:
    header_text = " ".join(header_lines)
    time_match = TIME_PATTERN.search(header_text)
    if time_match is None:
        raise MergeError(f"{path}: unrecognized track header: {header_text!r}")
    values = {key: value for key, value in KEY_PATTERN.findall(header_text)}

    def need_int(key: str) -> int:
        if key not in values:
            raise MergeError(f"{path}: missing header key {key}")
        return int(values[key])

    cycle = need_int("cycle")
    if "trk_format" not in values:
        raise MergeError(f"{path}: missing header key trk_format")
    trk_format = values["trk_format"]
    ntracked = need_int("ntracked_prtcls")
    ntrack_per_species = need_int("ntrack_per_species")
    track_per_species = bool(need_int("track_per_species"))
    record_count = need_int("record_count")
    nfields = need_int("nfields")
    layout = values.get("layout", "unknown")
    fields = tuple(values.get("fields", "").split(","))
    if fields == ("",):
        fields = ()
    return Header(
        time=float(time_match.group(1)),
        cycle=cycle,
        trk_format=trk_format,
        ntracked=ntracked,
        ntrack_per_species=ntrack_per_species,
        track_per_species=track_per_species,
        record_count=record_count,
        nfields=nfields,
        fields=fields,
        layout=layout,
    )


def strict_header_check(header: Header, path: Path) -> None:
    if header.trk_format not in ("rich_v1", "rich_v2"):
        raise MergeError(
            f"{path}: expected trk_format=rich_v1/rich_v2, got {header.trk_format!r}")
    if header.nfields != N_SOURCE_FIELDS:
        raise MergeError(f"{path}: expected nfields=18, got {header.nfields}")
    if header.fields != RICH_FIELDS:
        raise MergeError(f"{path}: wrong rich field list: {header.fields}")
    if header.record_count < 0:
        raise MergeError(f"{path}: negative record_count={header.record_count}")
    if header.ntracked <= 0:
        raise MergeError(f"{path}: non-positive ntracked={header.ntracked}")
    if header.ntrack_per_species <= 0:
        raise MergeError(
            f"{path}: non-positive ntrack_per_species={header.ntrack_per_species}")
    if header.track_per_species and header.ntracked % header.ntrack_per_species != 0:
        raise MergeError(
            f"{path}: ntracked={header.ntracked} is not divisible by "
            f"ntrack_per_species={header.ntrack_per_species}")


def check_records(path: Path, header: Header, records: np.ndarray, check_finite: bool) -> None:
    if not header.record_count:
        return
    tags = records[:, 0]
    times = records[:, 1]
    rounded_tags = np.rint(tags)
    if not np.all(np.isfinite(tags)) or not np.all(np.isfinite(times)):
        raise MergeError(f"{path}: non-finite tag/time at cycle={header.cycle}")
    if not np.allclose(tags, rounded_tags, rtol=0.0, atol=1.0e-4):
        raise MergeError(f"{path}: non-integral output_tag at cycle={header.cycle}")
    if np.any(rounded_tags < 0) or np.any(rounded_tags >= header.ntracked):
        raise MergeError(
            f"{path}: output_tag outside [0,{header.ntracked}) at cycle={header.cycle}")
    time_min = float(np.min(times))
    time_max = float(np.max(times))
    if abs(time_max - time_min) > 1.0e-6 * max(1.0, abs(time_min)):
        raise MergeError(
            f"{path}: payload time is not constant at cycle={header.cycle}: "
            f"{time_min:g}-{time_max:g}")
    if abs(time_min - header.time) > 1.0e-4 * max(1.0, abs(header.time)):
        raise MergeError(
            f"{path}: payload time {time_min:g} differs too much from "
            f"header time {header.time:g} at cycle={header.cycle}")
    if check_finite and not np.all(np.isfinite(records[:, 2:])):
        raise MergeError(f"{path}: non-finite rich values at cycle={header.cycle}")


def iter_legacy_trk_frames(
    path: Path, blob: bytes, check_finite: bool = True
) -> Iterator[Tuple[Header, np.ndarray]]:
    offset = 0
    found = False
    while True:
        start = blob.find(TRACK_HEADER_MARKER, offset)
        if start < 0:
            break
        found = True
        cursor = start
        header_lines: List[str] = []
        while True:
            line_end = blob.find(b"\n", cursor)
            if line_end < 0:
                raise MergeError(f"{path}: unterminated track header")
            raw_line = blob[cursor:line_end]
            cursor = line_end + 1
            if raw_line.strip() == b"":
                break
            header_lines.append(raw_line.decode("ascii", errors="replace"))
        header = parse_header(header_lines, path)
        strict_header_check(header, path)

        payload_values = header.record_count * N_SOURCE_FIELDS
        payload_bytes = payload_values * np.dtype("<f4").itemsize
        payload_end = cursor + payload_bytes
        if payload_end > len(blob):
            raise MergeError(
                f"{path}: incomplete payload at cycle={header.cycle} "
                f"time={header.time}; need {payload_bytes} bytes")
        records = np.frombuffer(blob, dtype="<f4", count=payload_values,
                                offset=cursor).reshape(header.record_count,
                                                       N_SOURCE_FIELDS)
        check_records(path, header, records, check_finite)
        yield header, records
        offset = payload_end
    if not found:
        raise MergeError(f"{path}: no tracked-particle frames found")


def compact_layout_name(code: int) -> str:
    return {0: "shared", 1: "node", 2: "rank"}.get(code, "unknown")


def iter_compact_trk_frames(
    path: Path, blob: bytes, check_finite: bool = True
) -> Iterator[Tuple[Header, np.ndarray]]:
    if len(blob) < COMPACT_PROLOGUE.size:
        raise MergeError(f"{path}: truncated compact trk prologue")
    unpacked = COMPACT_PROLOGUE.unpack_from(blob, 0)
    (magic, version, prologue_bytes, nfields, ntracked, ntrack_per_species,
     track_per_species, layout_code, _rank, _node, _nranks, _nnodes,
     _ranks_per_node, fields_bytes, _reserved) = unpacked
    if magic != COMPACT_FILE_MAGIC or version != COMPACT_VERSION:
        raise MergeError(f"{path}: invalid compact trk prologue")
    if prologue_bytes != COMPACT_PROLOGUE.size:
        raise MergeError(f"{path}: unsupported compact prologue size {prologue_bytes}")
    fields_start = prologue_bytes
    fields_end = fields_start + fields_bytes
    if fields_end > len(blob):
        raise MergeError(f"{path}: truncated compact trk field list")
    fields = tuple(blob[fields_start:fields_end].decode("ascii").split(","))
    layout = compact_layout_name(layout_code)
    cursor = fields_end
    found = False
    while cursor < len(blob):
        if cursor + COMPACT_FRAME.size > len(blob):
            raise MergeError(f"{path}: truncated compact frame header")
        (frame_magic, frame_version, frame_bytes, record_count, cycle, time,
         payload_bytes, _frame_reserved) = COMPACT_FRAME.unpack_from(blob, cursor)
        if frame_magic != COMPACT_FRAME_MAGIC or frame_version != COMPACT_VERSION:
            raise MergeError(f"{path}: invalid compact frame magic/version at {cursor}")
        if frame_bytes != COMPACT_FRAME.size:
            raise MergeError(f"{path}: unsupported compact frame size {frame_bytes}")
        cursor += frame_bytes
        expected_payload = record_count * N_SOURCE_FIELDS * np.dtype("<f4").itemsize
        if payload_bytes != expected_payload:
            raise MergeError(
                f"{path}: compact payload_bytes={payload_bytes}, expected "
                f"{expected_payload} at cycle={cycle}")
        payload_end = cursor + payload_bytes
        if payload_end > len(blob):
            raise MergeError(f"{path}: truncated compact payload at cycle={cycle}")
        header = Header(
            time=float(time),
            cycle=int(cycle),
            trk_format="rich_v2",
            ntracked=int(ntracked),
            ntrack_per_species=int(ntrack_per_species),
            track_per_species=bool(track_per_species),
            record_count=int(record_count),
            nfields=int(nfields),
            fields=fields,
            layout=layout,
        )
        strict_header_check(header, path)
        records = np.frombuffer(blob, dtype="<f4", count=record_count*N_SOURCE_FIELDS,
                                offset=cursor).reshape(record_count, N_SOURCE_FIELDS)
        check_records(path, header, records, check_finite)
        yield header, records
        found = True
        cursor = payload_end
    if not found:
        raise MergeError(f"{path}: no compact tracked-particle frames found")


def iter_trk_frames(path: Path, check_finite: bool = True) -> Iterator[Tuple[Header, np.ndarray]]:
    blob = path.read_bytes()
    if blob.startswith(COMPACT_FILE_MAGIC):
        yield from iter_compact_trk_frames(path, blob, check_finite)
    else:
        yield from iter_legacy_trk_frames(path, blob, check_finite)


def first_header(path: Path) -> Header:
    for header, _records in iter_trk_frames(path, check_finite=False):
        return header
    raise MergeError(f"{path}: no header")


def run_meta_from_header(header: Header) -> RunMeta:
    return RunMeta(
        trk_format=header.trk_format,
        ntracked=header.ntracked,
        ntrack_per_species=header.ntrack_per_species,
        track_per_species=header.track_per_species,
        nfields=header.nfields,
        fields=header.fields,
        layout=header.layout,
    )


def meta_key(meta: RunMeta) -> Tuple[object, ...]:
    return (meta.trk_format, meta.ntracked, meta.ntrack_per_species,
            meta.track_per_species, meta.nfields, meta.fields)


def owner_range(rank: int, size: int, ntracked: int) -> Tuple[int, int]:
    start = (rank * ntracked) // size
    end = ((rank + 1) * ntracked) // size
    return start, end


def owners_for_tags(tags: np.ndarray, size: int, ntracked: int) -> np.ndarray:
    owners = (tags.astype(np.int64) * size) // ntracked
    np.minimum(owners, size - 1, out=owners)
    return owners.astype(np.int32, copy=False)


def records_to_temp(header: Header, records: np.ndarray) -> np.ndarray:
    out = np.empty(records.shape[0], dtype=TEMP_RECORD_DTYPE)
    out["tag"] = np.rint(records[:, 0]).astype("<i8")
    out["cycle"] = header.cycle
    out["time"] = records[:, 1].astype("<f8")
    out["values"] = records[:, 2:].astype("<f4", copy=False)
    return out


def alltoall_temp_records(comm: MPI.Comm, records: np.ndarray,
                          ntracked: int) -> np.ndarray:
    size = comm.Get_size()
    if size == 1:
        return np.ascontiguousarray(records)

    if records.size == 0:
        send_counts = np.zeros(size, dtype=np.int64)
        sendbuf = np.empty(0, dtype=np.uint8)
    else:
        owners = owners_for_tags(records["tag"], size, ntracked)
        parts: List[np.ndarray] = []
        send_counts = np.zeros(size, dtype=np.int64)
        for dest in range(size):
            part = np.ascontiguousarray(records[owners == dest])
            parts.append(part)
            send_counts[dest] = part.nbytes
        if parts:
            sendbuf = np.concatenate(parts).view(np.uint8)
        else:
            sendbuf = np.empty(0, dtype=np.uint8)

    recv_counts = np.empty(size, dtype=np.int64)
    comm.Alltoall([send_counts, MPI.LONG_LONG], [recv_counts, MPI.LONG_LONG])
    if np.any(send_counts > np.iinfo(np.int32).max):
        raise MergeError("one MPI send segment exceeds int32 byte count; "
                         "reduce --exchange-rows")
    if np.any(recv_counts > np.iinfo(np.int32).max):
        raise MergeError("one MPI receive segment exceeds int32 byte count; "
                         "reduce --exchange-rows")

    send_displs = np.zeros(size, dtype=np.int32)
    recv_displs = np.zeros(size, dtype=np.int32)
    if size > 1:
        send_displs[1:] = np.cumsum(send_counts[:-1], dtype=np.int64).astype(np.int32)
        recv_displs[1:] = np.cumsum(recv_counts[:-1], dtype=np.int64).astype(np.int32)
    recvbuf = np.empty(int(recv_counts.sum()), dtype=np.uint8)

    comm.Alltoallv(
        [sendbuf, (send_counts.astype(np.int32), send_displs), MPI.BYTE],
        [recvbuf, (recv_counts.astype(np.int32), recv_displs), MPI.BYTE],
    )
    if recvbuf.nbytes % TEMP_RECORD_DTYPE.itemsize != 0:
        raise MergeError("received byte count is not a whole temp record")
    return np.frombuffer(recvbuf, dtype=TEMP_RECORD_DTYPE).copy()


def update_frame_meta(frame_meta: Dict[int, List[float]], header: Header,
                      records: np.ndarray) -> None:
    if records.size:
        time_min = float(np.min(records[:, 1]))
        time_max = float(np.max(records[:, 1]))
    else:
        # Empty shard frames are normal.  They contribute to header accounting
        # but must not set the global payload-time table.
        time_min = math.nan
        time_max = math.nan
    entry = frame_meta.setdefault(
        header.cycle,
        [math.nan, math.nan, 0, 0],
    )
    if math.isfinite(time_min):
        entry[0] = time_min if not math.isfinite(entry[0]) else min(entry[0], time_min)
        entry[1] = time_max if not math.isfinite(entry[1]) else max(entry[1], time_max)
    entry[2] += header.record_count
    entry[3] += 1


def write_frame_meta(path: Path, frame_meta: Dict[int, List[float]]) -> None:
    arr = np.empty(len(frame_meta), dtype=FRAME_META_DTYPE)
    for index, (cycle, (time_min, time_max, record_count, header_count)) in enumerate(
            sorted(frame_meta.items())):
        arr[index] = (cycle, time_min, time_max, int(record_count), int(header_count))
    np.save(path, arr)


def combine_frame_meta(tmp_dir: Path, size: int, ntracked: int) -> Tuple[np.ndarray, np.ndarray]:
    combined: Dict[int, List[float]] = {}
    for rank in range(size):
        path = tmp_dir / f"frame_meta_rank_{rank:06d}.npy"
        if not path.exists():
            raise MergeError(f"missing frame metadata shard: {path}")
        arr = np.load(path)
        for row in arr:
            cycle = int(row["cycle"])
            entry = combined.setdefault(cycle, [math.nan, math.nan, 0, 0])
            row_time_min = float(row["time_min"])
            row_time_max = float(row["time_max"])
            if math.isfinite(row_time_min):
                entry[0] = (row_time_min if not math.isfinite(entry[0])
                            else min(entry[0], row_time_min))
                entry[1] = (row_time_max if not math.isfinite(entry[1])
                            else max(entry[1], row_time_max))
            entry[2] += int(row["record_count"])
            entry[3] += int(row["header_count"])

    cycles = np.array(sorted(combined), dtype="<i8")
    times = np.empty(cycles.size, dtype="<f8")
    bad_counts: List[Tuple[int, int]] = []
    bad_times: List[Tuple[int, float, float]] = []
    for idx, cycle in enumerate(cycles):
        time_min, time_max, record_count, _header_count = combined[int(cycle)]
        times[idx] = time_min
        if record_count != ntracked:
            bad_counts.append((int(cycle), int(record_count)))
        if not math.isfinite(time_min) or not math.isfinite(time_max):
            bad_times.append((int(cycle), time_min, time_max))
        if abs(time_max - time_min) > 1.0e-8 * max(1.0, abs(time_min)):
            bad_times.append((int(cycle), time_min, time_max))
    if bad_counts:
        preview = ", ".join(f"{cycle}:{count}" for cycle, count in bad_counts[:5])
        raise MergeError(
            f"{len(bad_counts)} global frames do not sum to ntracked={ntracked}; "
            f"first {preview}")
    if bad_times:
        preview = ", ".join(
            f"{cycle}:{tmin:g}-{tmax:g}" for cycle, tmin, tmax in bad_times[:5])
        raise MergeError(
            f"{len(bad_times)} global frames have inconsistent payload times; "
            f"first {preview}")
    return cycles, times


def parse_and_redistribute(comm: MPI.Comm, files: Sequence[FileEntry], tmp_dir: Path,
                           exchange_rows: int,
                           check_finite: bool) -> Tuple[Optional[RunMeta], Dict[str, int]]:
    rank = comm.Get_rank()
    unsorted_path = tmp_dir / f"owner_{rank:06d}.unsorted.bin"
    frame_meta_path = tmp_dir / f"frame_meta_rank_{rank:06d}.npy"
    unsorted_path.parent.mkdir(parents=True, exist_ok=True)
    if unsorted_path.exists():
        unsorted_path.unlink()

    batch: List[np.ndarray] = []
    batch_rows = 0
    frame_meta: Dict[int, List[float]] = {}
    local_meta: Optional[RunMeta] = None
    parsed_files = 0
    parsed_frames = 0
    empty_frames = 0
    parsed_records = 0
    received_records = 0

    def flush_batch() -> None:
        nonlocal batch, batch_rows, received_records
        if batch_rows == 0:
            outbound = np.empty(0, dtype=TEMP_RECORD_DTYPE)
        elif len(batch) == 1:
            outbound = np.ascontiguousarray(batch[0])
        else:
            outbound = np.concatenate(batch)
        if local_meta is None:
            received = alltoall_temp_records(comm, outbound, 1)
        else:
            received = alltoall_temp_records(comm, outbound, local_meta.ntracked)
        if received.size:
            with unsorted_path.open("ab") as handle:
                received.tofile(handle)
            received_records += int(received.size)
        batch = []
        batch_rows = 0

    # MPI_Alltoallv is collective.  Flush at synchronized file boundaries so
    # every rank enters the same number of exchanges, while keeping production
    # rank-layout memory bounded by roughly one source rank shard per MPI rank.
    max_files = comm.allreduce(len(files), op=MPI.MAX)
    for file_index in range(max_files):
        if file_index < len(files):
            entry = files[file_index]
            path = Path(entry.path)
            parsed_files += 1
            for header, records in iter_trk_frames(path, check_finite=check_finite):
                meta = run_meta_from_header(header)
                if local_meta is None:
                    local_meta = meta
                elif meta_key(local_meta) != meta_key(meta):
                    raise MergeError(
                        f"{path}: metadata differs from earlier frames: {meta}")
                update_frame_meta(frame_meta, header, records)
                parsed_frames += 1
                parsed_records += int(header.record_count)
                if header.record_count == 0:
                    empty_frames += 1
                if header.record_count:
                    temp = records_to_temp(header, records)
                    batch.append(temp)
                    batch_rows += int(temp.size)
            if batch_rows > exchange_rows and comm.Get_rank() == 0:
                print("WARNING: one-file exchange batch exceeded --exchange-rows; "
                      "rank-layout production is still expected to fit in memory",
                      flush=True)
        flush_batch()
    write_frame_meta(frame_meta_path, frame_meta)

    stats = {
        "parsed_files": parsed_files,
        "parsed_frames": parsed_frames,
        "empty_frames": empty_frames,
        "parsed_records": parsed_records,
        "received_records": received_records,
    }
    return local_meta, stats


def empty_owner_files(tmp_dir: Path, rank: int) -> Tuple[Path, Path]:
    values_path = tmp_dir / f"owner_{rank:06d}.values.f32"
    index_path = tmp_dir / f"owner_{rank:06d}.index.npy"
    if values_path.exists():
        values_path.unlink()
    if index_path.exists():
        index_path.unlink()
    return values_path, index_path


def process_owner_shard(rank: int, size: int, tmp_dir: Path, meta: RunMeta,
                        cycles: np.ndarray, times: np.ndarray,
                        require_complete: bool) -> Dict[str, int]:
    start_tag, end_tag = owner_range(rank, size, meta.ntracked)
    unsorted_path = tmp_dir / f"owner_{rank:06d}.unsorted.bin"
    values_path, index_path = empty_owner_files(tmp_dir, rank)
    ntimes = int(cycles.size)

    if unsorted_path.exists() and unsorted_path.stat().st_size:
        if unsorted_path.stat().st_size % TEMP_RECORD_DTYPE.itemsize != 0:
            raise MergeError(f"{unsorted_path}: size is not a whole temp record")
        records = np.fromfile(unsorted_path, dtype=TEMP_RECORD_DTYPE)
    else:
        records = np.empty(0, dtype=TEMP_RECORD_DTYPE)

    if records.size:
        if np.any(records["tag"] < start_tag) or np.any(records["tag"] >= end_tag):
            bad = records[(records["tag"] < start_tag) | (records["tag"] >= end_tag)]["tag"][:5]
            raise MergeError(
                f"owner rank {rank} received tags outside [{start_tag},{end_tag}): "
                f"{bad.tolist()}")
        order = np.lexsort((records["cycle"], records["tag"]))
        records = records[order]

    index_rows: List[Tuple[int, int, int, int, int, float, float]] = []
    cursor = 0
    row = start_tag
    with values_path.open("ab") as values_handle:
        for output_tag in range(start_tag, end_tag):
            begin = cursor
            while cursor < records.size and int(records["tag"][cursor]) == output_tag:
                cursor += 1
            group = records[begin:cursor]
            count = int(group.size)
            if require_complete and count != ntimes:
                raise MergeError(
                    f"output_tag={output_tag} has {count} records; expected {ntimes}")
            if count:
                dup_cycles = np.flatnonzero(group["cycle"][1:] == group["cycle"][:-1])
                if dup_cycles.size:
                    raise MergeError(
                        f"output_tag={output_tag} has duplicate cycle "
                        f"{int(group['cycle'][int(dup_cycles[0])])}")
                if np.any(group["cycle"][1:] <= group["cycle"][:-1]):
                    raise MergeError(
                        f"output_tag={output_tag} cycles are not strictly increasing")
                if np.any(group["time"][1:] <= group["time"][:-1]):
                    raise MergeError(
                        f"output_tag={output_tag} times are not strictly increasing")
                if require_complete and not np.array_equal(group["cycle"], cycles):
                    missing = np.setdiff1d(cycles, group["cycle"], assume_unique=True)
                    extra = np.setdiff1d(group["cycle"], cycles, assume_unique=True)
                    raise MergeError(
                        f"output_tag={output_tag} cycle grid mismatch; "
                        f"missing={missing[:3].tolist()} extra={extra[:3].tolist()}")
                group["values"].astype("<f4", copy=False).tofile(values_handle)
                first_time = float(group["time"][0])
                last_time = float(group["time"][-1])
            else:
                first_time = math.nan
                last_time = math.nan
            species = (output_tag // meta.ntrack_per_species
                       if meta.track_per_species else 0)
            track_tag = (output_tag % meta.ntrack_per_species
                         if meta.track_per_species else output_tag)
            index_rows.append((output_tag, species, track_tag, row, count,
                               first_time, last_time))
            row += 1

    if cursor != records.size:
        raise MergeError(
            f"owner rank {rank} did not consume all sorted records; cursor={cursor} "
            f"size={records.size}")

    index = np.array(index_rows, dtype=PARTICLE_DTYPE)
    np.save(index_path, index)
    return {
        "owner_start_tag": start_tag,
        "owner_end_tag": end_tag,
        "owner_records": int(records.size),
        "owner_particles": int(index.size),
    }


def h5_compression_args(name: str) -> Dict[str, object]:
    if name == "none":
        return {}
    return {"compression": name}


def source_commit(run_dir: Path) -> str:
    path = run_dir / "source.commit"
    if path.exists():
        return path.read_text().strip()
    return ""


def assemble_hdf5(output: Path, tmp_dir: Path, run_dir: Path, entries: Sequence[FileEntry],
                  source_layout: str, meta: RunMeta, cycles: np.ndarray,
                  times: np.ndarray, mpi_size: int, args: argparse.Namespace) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    work_output = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    if work_output.exists():
        work_output.unlink()

    all_index_parts = []
    for owner in range(mpi_size):
        index_path = tmp_dir / f"owner_{owner:06d}.index.npy"
        if not index_path.exists():
            raise MergeError(f"missing owner index: {index_path}")
        all_index_parts.append(np.load(index_path))
    particles = (np.concatenate(all_index_parts)
                 if all_index_parts else np.empty(0, dtype=PARTICLE_DTYPE))
    if particles.size != meta.ntracked:
        raise MergeError(
            f"particle index has {particles.size} rows; expected {meta.ntracked}")
    order = np.argsort(particles["output_tag"])
    particles = particles[order]
    if not np.array_equal(particles["output_tag"], np.arange(meta.ntracked)):
        raise MergeError("particle index does not cover contiguous output tags")

    compression_args = h5_compression_args(args.compression)
    value_chunks = (1, min(max(int(args.assembly_rows), 1) * 8192,
                           max(int(cycles.size), 1)), N_VALUE_FIELDS)
    if cycles.size:
        value_chunks = (1, min(32768, int(cycles.size)), N_VALUE_FIELDS)

    source_bytes = sum(entry.size for entry in entries)
    with h5py.File(work_output, "w") as handle:
        handle.attrs["format"] = FORMAT_NAME
        handle.attrs["source_trk_format"] = meta.trk_format
        handle.attrs["source_run_dir"] = str(run_dir)
        handle.attrs["source_commit"] = source_commit(run_dir)
        handle.attrs["source_layout"] = source_layout
        handle.attrs["source_nfiles"] = len(entries)
        handle.attrs["source_bytes"] = source_bytes
        handle.attrs["nfields_source"] = N_SOURCE_FIELDS
        handle.attrs["nvalue_fields"] = N_VALUE_FIELDS
        handle.attrs["ntracked_prtcls"] = meta.ntracked
        handle.attrs["ntrack_per_species"] = meta.ntrack_per_species
        handle.attrs["nspecies"] = meta.nspecies
        handle.attrs["track_per_species"] = int(meta.track_per_species)
        handle.attrs["merge_command"] = " ".join(sys.argv)
        handle.attrs["merge_mpi_size"] = mpi_size

        handle.create_dataset("particles", data=particles)
        handle.create_dataset("times", data=times)
        handle.create_dataset("cycles", data=cycles)
        values = handle.create_dataset(
            "values",
            shape=(meta.ntracked, int(cycles.size), N_VALUE_FIELDS),
            dtype="<f4",
            chunks=value_chunks,
            **compression_args,
        )
        values.attrs["fields"] = ",".join(VALUE_FIELDS)

        for owner in range(mpi_size):
            index = np.load(tmp_dir / f"owner_{owner:06d}.index.npy")
            if index.size == 0:
                continue
            values_path = tmp_dir / f"owner_{owner:06d}.values.f32"
            expected_floats = int(index["count"].sum()) * N_VALUE_FIELDS
            if not values_path.exists():
                raise MergeError(f"missing owner values: {values_path}")
            actual_floats = values_path.stat().st_size // np.dtype("<f4").itemsize
            if actual_floats != expected_floats:
                raise MergeError(
                    f"{values_path}: has {actual_floats} floats; "
                    f"expected {expected_floats}")
            if index.size and int(index["count"].min()) != int(cycles.size):
                raise MergeError(
                    f"{values_path}: sparse owner output is not supported in HDF5 assembly")
            mmap = np.memmap(values_path, dtype="<f4", mode="r",
                             shape=(index.size, int(cycles.size), N_VALUE_FIELDS))
            for begin in range(0, index.size, max(1, int(args.assembly_rows))):
                end = min(index.size, begin + max(1, int(args.assembly_rows)))
                row0 = int(index["row"][begin])
                row1 = int(index["row"][end - 1]) + 1
                values[row0:row1, :, :] = mmap[begin:end, :, :]
            del mmap
    os.replace(work_output, output)


def validate_hdf5_readback(output: Path, tmp_dir: Path, meta: RunMeta,
                           cycles: np.ndarray, mpi_size: int) -> None:
    ntimes = int(cycles.size)
    sample_tags = {0, max(0, meta.ntracked - 1)}
    if meta.track_per_species:
        for species in range(meta.nspecies):
            sample_tags.add(species * meta.ntrack_per_species)
            sample_tags.add(min(meta.ntracked - 1,
                                (species + 1) * meta.ntrack_per_species - 1))
    rng = np.random.default_rng(8675309)
    if meta.ntracked > 0:
        random_count = min(32, meta.ntracked)
        sample_tags.update(int(tag) for tag in rng.choice(
            meta.ntracked, size=random_count, replace=False))

    samples: List[Tuple[int, int]] = []
    for output_tag in sorted(sample_tags):
        owner = int((output_tag * mpi_size) // meta.ntracked)
        owner = min(owner, mpi_size - 1)
        start_tag, _end_tag = owner_range(owner, mpi_size, meta.ntracked)
        samples.append((owner, output_tag - start_tag))

    with h5py.File(output, "r") as handle:
        if handle.attrs.get("format", "") != FORMAT_NAME:
            raise MergeError(f"{output}: wrong format attr")
        values = handle["values"]
        particles = handle["particles"][:]
        if values.shape != (meta.ntracked, ntimes, N_VALUE_FIELDS):
            raise MergeError(f"{output}: wrong values shape {values.shape}")
        if particles.shape[0] != meta.ntracked:
            raise MergeError(f"{output}: wrong particle index length")
        for owner, local_row in samples:
            index = np.load(tmp_dir / f"owner_{owner:06d}.index.npy")
            values_path = tmp_dir / f"owner_{owner:06d}.values.f32"
            mmap = np.memmap(values_path, dtype="<f4", mode="r",
                             shape=(index.size, ntimes, N_VALUE_FIELDS))
            global_row = int(index["row"][local_row])
            if not np.array_equal(values[global_row, :, :], mmap[local_row, :, :]):
                raise MergeError(
                    f"{output}: read-back mismatch owner={owner} row={local_row}")
            del mmap


def dry_run_summary(comm: MPI.Comm, run_dir: Path, layout: str,
                    entries: Sequence[FileEntry]) -> RunMeta:
    header = first_header(Path(entries[0].path))
    meta = run_meta_from_header(header)
    source_bytes = sum(entry.size for entry in entries)
    # This estimates complete-frame science payload; short in-progress runs will
    # produce fewer frames.
    estimated_ntimes = None
    if meta.ntracked and source_bytes:
        estimated_ntimes = source_bytes / (
            max(meta.ntracked, 1) * N_SOURCE_FIELDS * np.dtype("<f4").itemsize)
    if comm.Get_rank() == 0:
        print(json.dumps({
            "run_dir": str(run_dir),
            "layout": layout,
            "source_nfiles": len(entries),
            "source_bytes": source_bytes,
            "first_file": entries[0].path,
            "trk_format": meta.trk_format,
            "ntracked_prtcls": meta.ntracked,
            "ntrack_per_species": meta.ntrack_per_species,
            "track_per_species": int(meta.track_per_species),
            "nspecies": meta.nspecies,
            "nfields": meta.nfields,
            "fields": list(meta.fields),
            "estimated_ntimes_from_source_bytes": estimated_ntimes,
            "estimated_dense_value_bytes_per_200k_frames":
                meta.ntracked * 200000 * N_VALUE_FIELDS * np.dtype("<f4").itemsize,
        }, indent=2), flush=True)
    return meta


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if args.allow_missing:
        fail_all(comm, "--allow-missing is reserved for future sparse HDF5 support")
    if args.exchange_rows <= 0:
        fail_all(comm, "--exchange-rows must be positive")
    if args.assembly_rows <= 0:
        fail_all(comm, "--assembly-rows must be positive")
    if args.output is None and not (args.validate_only or args.dry_run):
        fail_all(comm, "--output is required unless --validate-only or --dry-run is set")

    try:
        if rank == 0:
            layout, entries = discover_files(args.run_dir, args.layout)
            assigned = assign_files(entries, size)
            args.tmp_dir.mkdir(parents=True, exist_ok=True)
            (args.tmp_dir / "manifest.json").write_text(json.dumps({
                "run_dir": str(args.run_dir),
                "layout": layout,
                "source_nfiles": len(entries),
                "source_bytes": sum(entry.size for entry in entries),
                "mpi_size": size,
            }, indent=2))
        else:
            layout = ""
            entries = []
            assigned = []
        layout = comm.bcast(layout, root=0)
        entries = comm.bcast(entries, root=0)
        assigned = comm.bcast(assigned, root=0)

        local_files = assigned[rank]
        if args.max_files_per_rank > 0:
            local_files = local_files[:args.max_files_per_rank]
        local_bytes = sum(entry.size for entry in local_files)
        all_local_counts = comm.gather((len(local_files), local_bytes), root=0)
        if rank == 0:
            rank_print(comm, "file assignment:",
                       f"ranks={size}",
                       f"files={sum(c for c, _b in all_local_counts)}",
                       f"bytes={sum(b for _c, b in all_local_counts)}")

        if args.dry_run and not args.validate_only:
            dry_run_summary(comm, args.run_dir, layout, entries)
            return 0

        rank_print(comm, "parsing and redistributing rich trk records")
        local_meta, local_stats = parse_and_redistribute(
            comm=comm,
            files=local_files,
            tmp_dir=args.tmp_dir,
            exchange_rows=args.exchange_rows,
            check_finite=not args.no_payload_finite_check,
        )
        gathered_meta = comm.gather(local_meta, root=0)
        gathered_stats = comm.gather(local_stats, root=0)

        if rank == 0:
            metas = [item for item in gathered_meta if item is not None]
            if not metas:
                raise MergeError("no rank parsed any trk frames")
            meta = metas[0]
            for other in metas[1:]:
                if meta_key(meta) != meta_key(other):
                    raise MergeError(f"inconsistent run metadata: {meta} vs {other}")
            total_stats = {
                key: sum(int(stats[key]) for stats in gathered_stats)
                for key in gathered_stats[0]
            }
            print("parse stats:", json.dumps(total_stats, sort_keys=True), flush=True)
            cycles, times = combine_frame_meta(args.tmp_dir, size, meta.ntracked)
            np.save(args.tmp_dir / "global_cycles.npy", cycles)
            np.save(args.tmp_dir / "global_times.npy", times)
            print(json.dumps({
                "ntracked_prtcls": meta.ntracked,
                "ntrack_per_species": meta.ntrack_per_species,
                "track_per_species": int(meta.track_per_species),
                "nspecies": meta.nspecies,
                "ntimes": int(cycles.size),
                "first_cycle": int(cycles[0]) if cycles.size else None,
                "last_cycle": int(cycles[-1]) if cycles.size else None,
                "first_time": float(times[0]) if times.size else None,
                "last_time": float(times[-1]) if times.size else None,
                "estimated_values_bytes": int(meta.ntracked * cycles.size *
                                              N_VALUE_FIELDS *
                                              np.dtype("<f4").itemsize),
            }, indent=2), flush=True)
        else:
            meta = None
            cycles = None
            times = None
        meta = comm.bcast(meta, root=0)
        cycles = comm.bcast(cycles, root=0)
        times = comm.bcast(times, root=0)

        rank_print(comm, "sorting and validating owner shards")
        owner_stats = process_owner_shard(
            rank=rank,
            size=size,
            tmp_dir=args.tmp_dir,
            meta=meta,
            cycles=cycles,
            times=times,
            require_complete=args.require_complete,
        )
        gathered_owner_stats = comm.gather(owner_stats, root=0)
        comm.Barrier()
        if rank == 0:
            print("owner stats:", json.dumps({
                "owner_records": sum(item["owner_records"] for item in gathered_owner_stats),
                "owner_particles": sum(item["owner_particles"] for item in gathered_owner_stats),
            }, sort_keys=True), flush=True)

            if not args.validate_only and not args.dry_run:
                assert args.output is not None
                print(f"assembling HDF5: {args.output}", flush=True)
                assemble_hdf5(args.output, args.tmp_dir, args.run_dir, entries,
                              layout, meta, cycles, times, size, args)
                validate_hdf5_readback(args.output, args.tmp_dir, meta, cycles, size)
                print(f"HDF5 read-back validation passed: {args.output}", flush=True)
            else:
                print("validation-only path complete; HDF5 assembly skipped", flush=True)

            if args.delete_temp and not args.dry_run and not args.validate_only:
                shutil.rmtree(args.tmp_dir)
                print(f"deleted temp directory: {args.tmp_dir}", flush=True)
            elif args.delete_temp and args.validate_only:
                print("--delete-temp ignored for --validate-only; no HDF5 read-back "
                      "validation was performed", flush=True)
        comm.Barrier()
        return 0
    except Exception as exc:  # noqa: BLE001 - propagate through MPI_Abort.
        fail_all(comm, str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
