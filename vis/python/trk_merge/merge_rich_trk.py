#!/usr/bin/env python3
"""Merge AthenaK rich tracked-particle shards into particle-major HDF5."""

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

import h5py
from mpi4py import MPI
import numpy as np


HEADER_MARKER = b"# AthenaK tracked particle data at time="
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
FORMAT_NAME = "athenak_rich_trk_merged_v1"

N_SOURCE_FIELDS = len(RICH_FIELDS)
N_VALUE_FIELDS = len(VALUE_FIELDS)

RECORD_DTYPE = np.dtype([
    ("tag", "<i8"),
    ("cycle", "<i8"),
    ("time", "<f8"),
    ("values", "<f4", (N_VALUE_FIELDS,)),
])

FRAME_DTYPE = np.dtype([
    ("cycle", "<i8"),
    ("time_min", "<f8"),
    ("time_max", "<f8"),
    ("records", "<i8"),
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

KEY_RE = re.compile(r"([A-Za-z_]+)=\s*([^ \t\n]+)")
TIME_RE = re.compile(r"# AthenaK tracked particle data at time=\s*([0-9.eE+-]+)")


@dataclass(frozen=True)
class FileEntry:
    path: str
    size: int


@dataclass(frozen=True)
class Header:
    time: float
    cycle: int
    trk_format: str
    ntracked: int
    ntrack_per_species: int
    track_per_species: bool
    record_count: int
    layout: str


@dataclass(frozen=True)
class RunMeta:
    trk_format: str
    ntracked: int
    ntrack_per_species: int
    track_per_species: bool
    layout: str

    @property
    def nspecies(self) -> int:
        if self.track_per_species:
            return self.ntracked // self.ntrack_per_species
        return 1


class MergeError(RuntimeError):
    """Raised for malformed input that would produce an unusable merge."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge AthenaK rich trk shards into particle-major HDF5.")
    parser.add_argument("--run-dir", required=True, type=Path,
                        help="Run directory containing trk/.")
    parser.add_argument("--output", required=True, type=Path,
                        help="Output HDF5 path.")
    parser.add_argument("--tmp-dir", required=True, type=Path,
                        help="Temporary directory for MPI owner shards.")
    parser.add_argument("--layout", choices=("auto", "rank", "node", "shared"),
                        default="auto",
                        help="Source trk layout. Use rank for production shards.")
    parser.add_argument("--assembly-rows", type=int, default=4,
                        help="Particle rows per HDF5 write in final assembly.")
    parser.add_argument("--keep-temp", action="store_true",
                        help="Keep temporary owner shards after a successful merge.")
    return parser.parse_args()


def rank0(comm: MPI.Comm, *items: object) -> None:
    if comm.Get_rank() == 0:
        print(*items, flush=True)


def abort(comm: MPI.Comm, message: str) -> None:
    print(f"[rank {comm.Get_rank()}] ERROR: {message}", file=sys.stderr, flush=True)
    comm.Abort(2)


def discover_files(run_dir: Path, layout: str) -> tuple[str, list[FileEntry]]:
    trk_dir = run_dir / "trk"
    patterns = {
        "rank": "rank_*/*.trk",
        "node": "node_*/*.trk",
        "shared": "*.trk",
    }
    order = ("rank", "node", "shared") if layout == "auto" else (layout,)
    for name in order:
        paths = sorted(trk_dir.glob(patterns[name]))
        if paths:
            return name, [FileEntry(str(path), path.stat().st_size) for path in paths]
    raise MergeError(f"no trk files found under {trk_dir}")


def assign_files(entries: list[FileEntry], nranks: int) -> list[list[FileEntry]]:
    assigned: list[list[FileEntry]] = [[] for _ in range(nranks)]
    loads = [0] * nranks
    for entry in sorted(entries, key=lambda item: item.size, reverse=True):
        rank = min(range(nranks), key=loads.__getitem__)
        assigned[rank].append(entry)
        loads[rank] += entry.size
    return assigned


def parse_header(lines: list[str], path: Path) -> Header:
    text = " ".join(lines)
    time_match = TIME_RE.search(text)
    if time_match is None:
        raise MergeError(f"{path}: missing tracked-particle time header")

    values = {key: value for key, value in KEY_RE.findall(text)}
    fields = tuple(values.get("fields", "").split(","))
    trk_format = values.get("trk_format")
    if trk_format != "rich_v1":
        raise MergeError(f"{path}: expected trk_format=rich_v1")
    if int(values["nfields"]) != N_SOURCE_FIELDS or fields != RICH_FIELDS:
        raise MergeError(f"{path}: unexpected rich trk field list")

    return Header(
        time=float(time_match.group(1)),
        cycle=int(values["cycle"]),
        trk_format=trk_format,
        ntracked=int(values["ntracked_prtcls"]),
        ntrack_per_species=int(values["ntrack_per_species"]),
        track_per_species=bool(int(values["track_per_species"])),
        record_count=int(values["record_count"]),
        layout=values.get("layout", "unknown"),
    )


def check_frame(path: Path, header: Header, frame: np.ndarray) -> None:
    if not header.record_count:
        return
    tags = frame[:, 0]
    rounded = np.rint(tags)
    if not np.allclose(tags, rounded, rtol=0.0, atol=1.0e-4):
        raise MergeError(f"{path}: non-integral output_tag at cycle {header.cycle}")
    if np.any(rounded < 0) or np.any(rounded >= header.ntracked):
        raise MergeError(f"{path}: output_tag outside range at cycle {header.cycle}")
    times = frame[:, 1]
    if not np.all(np.isfinite(frame)):
        raise MergeError(f"{path}: non-finite rich trk payload at cycle {header.cycle}")
    if abs(float(np.max(times)) - float(np.min(times))) > 1.0e-6:
        raise MergeError(f"{path}: non-constant payload time at cycle {header.cycle}")


def iter_legacy_frames(path: Path, blob: bytes):
    offset = 0
    while True:
        start = blob.find(HEADER_MARKER, offset)
        if start < 0:
            return

        cursor = start
        lines: list[str] = []
        while True:
            line_end = blob.find(b"\n", cursor)
            if line_end < 0:
                raise MergeError(f"{path}: unterminated trk header")
            raw = blob[cursor:line_end]
            cursor = line_end + 1
            if raw.strip() == b"":
                break
            lines.append(raw.decode("ascii", errors="replace"))

        header = parse_header(lines, path)
        nvalues = header.record_count * N_SOURCE_FIELDS
        payload_end = cursor + nvalues * np.dtype("<f4").itemsize
        if payload_end > len(blob):
            raise MergeError(f"{path}: incomplete payload at cycle {header.cycle}")

        values = np.frombuffer(blob, dtype="<f4", count=nvalues, offset=cursor)
        frame = values.reshape(header.record_count, N_SOURCE_FIELDS)
        check_frame(path, header, frame)
        yield header, frame
        offset = payload_end


def compact_layout_name(code: int) -> str:
    return {0: "shared", 1: "node", 2: "rank"}.get(code, "unknown")


def iter_compact_frames(path: Path, blob: bytes):
    if len(blob) < COMPACT_PROLOGUE.size:
        raise MergeError(f"{path}: truncated compact trk prologue")
    (magic, version, prologue_bytes, nfields, ntracked, ntrack_per_species,
     track_per_species, layout_code, _rank, _node, _nranks, _nnodes,
     _ranks_per_node, fields_bytes, _reserved) = COMPACT_PROLOGUE.unpack_from(blob, 0)
    if magic != COMPACT_FILE_MAGIC or version != COMPACT_VERSION:
        raise MergeError(f"{path}: invalid compact trk prologue")
    if prologue_bytes != COMPACT_PROLOGUE.size:
        raise MergeError(f"{path}: unsupported compact prologue size {prologue_bytes}")
    fields_start = prologue_bytes
    fields_end = fields_start + fields_bytes
    if fields_end > len(blob):
        raise MergeError(f"{path}: truncated compact field list")
    fields = tuple(blob[fields_start:fields_end].decode("ascii").split(","))
    if nfields != N_SOURCE_FIELDS or fields != RICH_FIELDS:
        raise MergeError(f"{path}: unexpected compact rich trk field list")

    cursor = fields_end
    while cursor < len(blob):
        if cursor + COMPACT_FRAME.size > len(blob):
            raise MergeError(f"{path}: truncated compact frame header")
        (frame_magic, frame_version, frame_bytes, record_count, cycle, time,
         payload_bytes, _frame_reserved) = COMPACT_FRAME.unpack_from(blob, cursor)
        if frame_magic != COMPACT_FRAME_MAGIC or frame_version != COMPACT_VERSION:
            raise MergeError(f"{path}: invalid compact frame at offset {cursor}")
        if frame_bytes != COMPACT_FRAME.size:
            raise MergeError(f"{path}: unsupported compact frame size {frame_bytes}")
        cursor += frame_bytes
        expected = record_count * N_SOURCE_FIELDS * np.dtype("<f4").itemsize
        if payload_bytes != expected:
            raise MergeError(f"{path}: compact payload size mismatch at cycle {cycle}")
        payload_end = cursor + payload_bytes
        if payload_end > len(blob):
            raise MergeError(f"{path}: truncated compact payload at cycle {cycle}")
        header = Header(
            time=float(time),
            cycle=int(cycle),
            trk_format="rich_v2",
            ntracked=int(ntracked),
            ntrack_per_species=int(ntrack_per_species),
            track_per_species=bool(track_per_species),
            record_count=int(record_count),
            layout=compact_layout_name(layout_code),
        )
        frame = np.frombuffer(blob, dtype="<f4", count=record_count*N_SOURCE_FIELDS,
                              offset=cursor).reshape(record_count, N_SOURCE_FIELDS)
        check_frame(path, header, frame)
        yield header, frame
        cursor = payload_end


def iter_frames(path: Path):
    blob = path.read_bytes()
    if blob.startswith(COMPACT_FILE_MAGIC):
        yield from iter_compact_frames(path, blob)
    else:
        yield from iter_legacy_frames(path, blob)


def meta_from_header(header: Header) -> RunMeta:
    return RunMeta(
        trk_format=header.trk_format,
        ntracked=header.ntracked,
        ntrack_per_species=header.ntrack_per_species,
        track_per_species=header.track_per_species,
        layout=header.layout,
    )


def meta_key(meta: RunMeta) -> tuple[object, ...]:
    return (meta.trk_format, meta.ntracked, meta.ntrack_per_species, meta.track_per_species)


def records_to_owner_rows(header: Header, frame: np.ndarray) -> np.ndarray:
    rows = np.empty(frame.shape[0], dtype=RECORD_DTYPE)
    rows["tag"] = np.rint(frame[:, 0]).astype("<i8")
    rows["cycle"] = header.cycle
    rows["time"] = frame[:, 1].astype("<f8")
    rows["values"] = frame[:, 2:].astype("<f4", copy=False)
    return rows


def owner_range(rank: int, nranks: int, ntracked: int) -> tuple[int, int]:
    return (rank * ntracked) // nranks, ((rank + 1) * ntracked) // nranks


def owner_ranks(tags: np.ndarray, nranks: int, ntracked: int) -> np.ndarray:
    owners = (tags.astype(np.int64) * nranks) // ntracked
    np.minimum(owners, nranks - 1, out=owners)
    return owners.astype(np.int32, copy=False)


def alltoall_records(comm: MPI.Comm, rows: np.ndarray, ntracked: int) -> np.ndarray:
    nranks = comm.Get_size()
    if nranks == 1:
        return np.ascontiguousarray(rows)

    if rows.size:
        owners = owner_ranks(rows["tag"], nranks, ntracked)
        counts = np.bincount(owners, minlength=nranks).astype(np.int64)
        order = np.argsort(owners, kind="stable")
        sendbuf = np.ascontiguousarray(rows[order]).view(np.uint8)
        send_counts = counts * RECORD_DTYPE.itemsize
    else:
        sendbuf = np.empty(0, dtype=np.uint8)
        send_counts = np.zeros(nranks, dtype=np.int64)

    recv_counts = np.empty(nranks, dtype=np.int64)
    comm.Alltoall([send_counts, MPI.LONG_LONG], [recv_counts, MPI.LONG_LONG])

    if send_counts.max(initial=0) > np.iinfo(np.int32).max:
        raise MergeError("one send segment exceeds the MPI int count limit")
    if recv_counts.max(initial=0) > np.iinfo(np.int32).max:
        raise MergeError("one receive segment exceeds the MPI int count limit")

    send_displs = np.zeros(nranks, dtype=np.int32)
    recv_displs = np.zeros(nranks, dtype=np.int32)
    send_displs[1:] = np.cumsum(send_counts[:-1], dtype=np.int64).astype(np.int32)
    recv_displs[1:] = np.cumsum(recv_counts[:-1], dtype=np.int64).astype(np.int32)
    recvbuf = np.empty(int(recv_counts.sum()), dtype=np.uint8)

    comm.Alltoallv(
        [sendbuf, (send_counts.astype(np.int32), send_displs), MPI.BYTE],
        [recvbuf, (recv_counts.astype(np.int32), recv_displs), MPI.BYTE],
    )
    return np.frombuffer(recvbuf, dtype=RECORD_DTYPE).copy()


def add_frame_meta(frames: dict[int, list[float]], header: Header,
                   frame: np.ndarray) -> None:
    if header.record_count:
        time_min = float(np.min(frame[:, 1]))
        time_max = float(np.max(frame[:, 1]))
    else:
        time_min = math.nan
        time_max = math.nan

    item = frames.setdefault(header.cycle, [math.nan, math.nan, 0])
    if math.isfinite(time_min):
        item[0] = time_min if not math.isfinite(item[0]) else min(item[0], time_min)
        item[1] = time_max if not math.isfinite(item[1]) else max(item[1], time_max)
    item[2] += header.record_count


def write_frame_meta(path: Path, frames: dict[int, list[float]]) -> None:
    out = np.empty(len(frames), dtype=FRAME_DTYPE)
    for i, (cycle, (time_min, time_max, records)) in enumerate(sorted(frames.items())):
        out[i] = (cycle, time_min, time_max, int(records))
    np.save(path, out)


def parse_and_redistribute(comm: MPI.Comm, local_files: list[FileEntry],
                           tmp_dir: Path) -> tuple[RunMeta | None, dict[str, int]]:
    rank = comm.Get_rank()
    unsorted = tmp_dir / f"owner_{rank:06d}.unsorted.bin"
    frames_path = tmp_dir / f"frames_{rank:06d}.npy"
    unsorted.unlink(missing_ok=True)

    frames: dict[int, list[float]] = {}
    local_meta: RunMeta | None = None
    stats = {"files": 0, "frames": 0, "records": 0, "received": 0}

    max_files = comm.allreduce(len(local_files), op=MPI.MAX)
    for file_index in range(max_files):
        pending: list[np.ndarray] = []
        if file_index < len(local_files):
            path = Path(local_files[file_index].path)
            stats["files"] += 1
            for header, frame in iter_frames(path):
                meta = meta_from_header(header)
                if local_meta is None:
                    local_meta = meta
                elif meta_key(local_meta) != meta_key(meta):
                    raise MergeError(f"{path}: inconsistent trk metadata")
                add_frame_meta(frames, header, frame)
                if header.record_count:
                    pending.append(records_to_owner_rows(header, frame))
                stats["frames"] += 1
                stats["records"] += header.record_count

        if pending:
            outbound = np.concatenate(pending) if len(pending) > 1 else pending[0]
            assert local_meta is not None
            received = alltoall_records(comm, outbound, local_meta.ntracked)
        else:
            ntracked = local_meta.ntracked if local_meta is not None else 1
            received = alltoall_records(comm, np.empty(0, dtype=RECORD_DTYPE), ntracked)

        if received.size:
            with unsorted.open("ab") as handle:
                received.tofile(handle)
            stats["received"] += int(received.size)

    write_frame_meta(frames_path, frames)
    return local_meta, stats


def combine_frames(tmp_dir: Path, nranks: int, ntracked: int) -> tuple[np.ndarray, np.ndarray]:
    combined: dict[int, list[float]] = {}
    for rank in range(nranks):
        for row in np.load(tmp_dir / f"frames_{rank:06d}.npy"):
            cycle = int(row["cycle"])
            item = combined.setdefault(cycle, [math.nan, math.nan, 0])
            tmin = float(row["time_min"])
            tmax = float(row["time_max"])
            if math.isfinite(tmin):
                item[0] = tmin if not math.isfinite(item[0]) else min(item[0], tmin)
                item[1] = tmax if not math.isfinite(item[1]) else max(item[1], tmax)
            item[2] += int(row["records"])

    cycles = np.array(sorted(combined), dtype="<i8")
    times = np.empty(cycles.size, dtype="<f8")
    for i, cycle in enumerate(cycles):
        tmin, tmax, records = combined[int(cycle)]
        if records != ntracked:
            raise MergeError(f"cycle {cycle} has {records} records, expected {ntracked}")
        if not math.isfinite(tmin) or abs(tmax - tmin) > 1.0e-8 * max(1.0, abs(tmin)):
            raise MergeError(f"cycle {cycle} has inconsistent payload times")
        times[i] = tmin
    return cycles, times


def process_owner(rank: int, nranks: int, tmp_dir: Path, meta: RunMeta,
                  cycles: np.ndarray, times: np.ndarray) -> dict[str, int]:
    del times
    start_tag, end_tag = owner_range(rank, nranks, meta.ntracked)
    unsorted = tmp_dir / f"owner_{rank:06d}.unsorted.bin"
    values_path = tmp_dir / f"owner_{rank:06d}.values.f32"
    index_path = tmp_dir / f"owner_{rank:06d}.index.npy"
    values_path.unlink(missing_ok=True)
    index_path.unlink(missing_ok=True)

    if unsorted.exists() and unsorted.stat().st_size:
        rows = np.fromfile(unsorted, dtype=RECORD_DTYPE)
        rows = rows[np.lexsort((rows["cycle"], rows["tag"]))]
    else:
        rows = np.empty(0, dtype=RECORD_DTYPE)

    index_rows = []
    cursor = 0
    ntimes = int(cycles.size)
    with values_path.open("ab") as handle:
        for tag in range(start_tag, end_tag):
            begin = cursor
            while cursor < rows.size and int(rows["tag"][cursor]) == tag:
                cursor += 1
            group = rows[begin:cursor]
            if group.size != ntimes:
                raise MergeError(f"output_tag {tag} has {group.size} records, expected {ntimes}")
            if not np.array_equal(group["cycle"], cycles):
                raise MergeError(f"output_tag {tag} does not match the global cycle grid")

            group["values"].astype("<f4", copy=False).tofile(handle)
            species = tag // meta.ntrack_per_species if meta.track_per_species else 0
            track_tag = tag % meta.ntrack_per_species if meta.track_per_species else tag
            index_rows.append((tag, species, track_tag, tag, ntimes,
                               float(group["time"][0]), float(group["time"][-1])))

    if cursor != rows.size:
        raise MergeError(f"owner rank {rank} received records outside its tag range")

    np.save(index_path, np.array(index_rows, dtype=PARTICLE_DTYPE))
    return {"records": int(rows.size), "particles": len(index_rows)}


def source_commit(run_dir: Path) -> str:
    path = run_dir / "source.commit"
    return path.read_text().strip() if path.exists() else ""


def assemble_hdf5(output: Path, tmp_dir: Path, run_dir: Path,
                  entries: list[FileEntry], layout: str, meta: RunMeta,
                  cycles: np.ndarray, times: np.ndarray, nranks: int,
                  assembly_rows: int) -> None:
    if output.exists():
        raise MergeError(f"output already exists: {output}")

    indexes = [np.load(tmp_dir / f"owner_{rank:06d}.index.npy")
               for rank in range(nranks)]
    particles = np.concatenate(indexes) if indexes else np.empty(0, dtype=PARTICLE_DTYPE)
    particles = particles[np.argsort(particles["output_tag"])]
    if particles.size != meta.ntracked:
        raise MergeError(f"particle index has {particles.size} rows, expected {meta.ntracked}")
    if not np.array_equal(particles["output_tag"], np.arange(meta.ntracked)):
        raise MergeError("particle index does not cover every output_tag")

    output.parent.mkdir(parents=True, exist_ok=True)
    work = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    work.unlink(missing_ok=True)

    with h5py.File(work, "w") as handle:
        handle.attrs["format"] = FORMAT_NAME
        handle.attrs["source_trk_format"] = meta.trk_format
        handle.attrs["source_run_dir"] = str(run_dir)
        handle.attrs["source_commit"] = source_commit(run_dir)
        handle.attrs["source_layout"] = layout
        handle.attrs["source_nfiles"] = len(entries)
        handle.attrs["source_bytes"] = sum(entry.size for entry in entries)
        handle.attrs["nfields_source"] = N_SOURCE_FIELDS
        handle.attrs["nvalue_fields"] = N_VALUE_FIELDS
        handle.attrs["ntracked_prtcls"] = meta.ntracked
        handle.attrs["ntrack_per_species"] = meta.ntrack_per_species
        handle.attrs["nspecies"] = meta.nspecies
        handle.attrs["track_per_species"] = int(meta.track_per_species)
        handle.attrs["merge_mpi_size"] = nranks

        handle.create_dataset("particles", data=particles)
        handle.create_dataset("times", data=times)
        handle.create_dataset("cycles", data=cycles)
        values = handle.create_dataset(
            "values",
            shape=(meta.ntracked, int(cycles.size), N_VALUE_FIELDS),
            dtype="<f4",
            chunks=(1, min(32768, max(1, int(cycles.size))), N_VALUE_FIELDS),
        )
        values.attrs["fields"] = ",".join(VALUE_FIELDS)

        rows_per_write = max(1, int(assembly_rows))
        for rank, index in enumerate(indexes):
            if index.size == 0:
                continue
            path = tmp_dir / f"owner_{rank:06d}.values.f32"
            mmap = np.memmap(path, dtype="<f4", mode="r",
                             shape=(index.size, int(cycles.size), N_VALUE_FIELDS))
            for begin in range(0, index.size, rows_per_write):
                end = min(index.size, begin + rows_per_write)
                row0 = int(index["row"][begin])
                row1 = int(index["row"][end - 1]) + 1
                values[row0:row1, :, :] = mmap[begin:end, :, :]
            del mmap

    os.replace(work, output)


def main() -> int:
    args = parse_args()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.Get_size()

    try:
        if args.assembly_rows <= 0:
            raise MergeError("--assembly-rows must be positive")

        if rank == 0:
            layout, entries = discover_files(args.run_dir, args.layout)
            assigned = assign_files(entries, nranks)
            args.tmp_dir.mkdir(parents=True, exist_ok=True)
            (args.tmp_dir / "manifest.json").write_text(json.dumps({
                "run_dir": str(args.run_dir),
                "layout": layout,
                "source_nfiles": len(entries),
                "source_bytes": sum(entry.size for entry in entries),
                "mpi_size": nranks,
            }, indent=2))
        else:
            layout, entries, assigned = "", [], []

        layout = comm.bcast(layout, root=0)
        entries = comm.bcast(entries, root=0)
        assigned = comm.bcast(assigned, root=0)
        comm.Barrier()

        local_files = assigned[rank]
        rank0(comm, "file assignment:",
              f"ranks={nranks}",
              f"files={sum(len(item) for item in assigned)}",
              f"bytes={sum(entry.size for entry in entries)}")

        rank0(comm, "parsing and redistributing rich trk records")
        local_meta, local_stats = parse_and_redistribute(comm, local_files, args.tmp_dir)
        metas = comm.gather(local_meta, root=0)
        stats = comm.gather(local_stats, root=0)

        if rank == 0:
            present = [meta for meta in metas if meta is not None]
            if not present:
                raise MergeError("no trk frames were parsed")
            meta = present[0]
            if any(meta_key(item) != meta_key(meta) for item in present[1:]):
                raise MergeError("inconsistent run metadata across source files")
            cycles, times = combine_frames(args.tmp_dir, nranks, meta.ntracked)
            print("parse stats:", json.dumps({
                key: sum(item[key] for item in stats) for key in stats[0]
            }, sort_keys=True), flush=True)
            print(json.dumps({
                "ntracked_prtcls": meta.ntracked,
                "ntrack_per_species": meta.ntrack_per_species,
                "track_per_species": int(meta.track_per_species),
                "nspecies": meta.nspecies,
                "ntimes": int(cycles.size),
                "first_cycle": int(cycles[0]),
                "last_cycle": int(cycles[-1]),
                "first_time": float(times[0]),
                "last_time": float(times[-1]),
                "estimated_values_bytes": int(meta.ntracked * cycles.size *
                                              N_VALUE_FIELDS *
                                              np.dtype("<f4").itemsize),
            }, indent=2), flush=True)
        else:
            meta, cycles, times = None, None, None

        meta = comm.bcast(meta, root=0)
        cycles = comm.bcast(cycles, root=0)
        times = comm.bcast(times, root=0)

        rank0(comm, "sorting owner shards")
        owner_stats = process_owner(rank, nranks, args.tmp_dir, meta, cycles, times)
        all_owner_stats = comm.gather(owner_stats, root=0)
        comm.Barrier()

        if rank == 0:
            print("owner stats:", json.dumps({
                "owner_records": sum(item["records"] for item in all_owner_stats),
                "owner_particles": sum(item["particles"] for item in all_owner_stats),
            }, sort_keys=True), flush=True)
            print(f"assembling HDF5: {args.output}", flush=True)
            assemble_hdf5(args.output, args.tmp_dir, args.run_dir, entries, layout,
                          meta, cycles, times, nranks, args.assembly_rows)
            print(f"wrote {args.output}", flush=True)
            if not args.keep_temp:
                shutil.rmtree(args.tmp_dir)
                print(f"deleted temp directory: {args.tmp_dir}", flush=True)

        comm.Barrier()
        return 0
    except Exception as exc:  # noqa: BLE001
        abort(comm, str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
