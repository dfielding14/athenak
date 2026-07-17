#!/usr/bin/env python3
"""Read matching AthenaK MHD MeshBlocks and tracked-particle trajectories."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import struct
from typing import BinaryIO, Iterator, Sequence

import h5py
import numpy as np

import bin_convert


LEGACY_MARKER = b"# AthenaK tracked particle data at time="
COMPACT_FILE_MAGIC = b"AKTRK2F\0"
COMPACT_FRAME_MAGIC = b"AKTRK2R\0"
COMPACT_VERSION = 2
COMPACT_PROLOGUE = struct.Struct("<8sHHIqqIIiiiiiII")
COMPACT_FRAME = struct.Struct("<8sHHIqdII")
KEY_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=\s*([^ \t\n]+)")


@dataclass(frozen=True)
class TrackFrame:
    """One frame from an AthenaK tracked-particle shard."""

    cycle: int
    time: float
    fields: tuple[str, ...]
    values: np.ndarray


@dataclass(frozen=True)
class TrackVisits:
    """Particle identities observed inside a spatial box."""

    output_tags: np.ndarray
    frames_seen: int
    records_seen: int
    records_inside: int


def read_meshblock(
    filename: str | Path,
    meshblock_index: int = 0,
    quantities: Sequence[str] | None = None,
    dtype: np.dtype | type = np.float32,
) -> dict:
    """Read one MeshBlock with AthenaK's single-rank binary reader."""

    data = bin_convert.read_single_rank_binary_as_athdf(
        str(filename),
        meshblock_index=meshblock_index,
        quantities=None if quantities is None else list(quantities),
        dtype=dtype,
    )
    for name in data["VariableNames"]:
        data[name] = np.array(data[name], dtype=dtype, copy=True)
    data["SourceFile"] = str(filename)
    return data


def _meshblock_from_raw(
    filedata: dict,
    meshblock_index: int,
    quantities: Sequence[str],
    dtype: np.dtype | type,
    source: str,
) -> dict:
    """Convert one already-read raw MeshBlock without rereading its rank file."""

    geometry = filedata["mb_geometry"][meshblock_index]
    first = np.asarray(filedata["mb_data"][quantities[0]][meshblock_index])
    nz, ny, nx = first.shape
    shape = (nx, ny, nz)
    data: dict[str, object] = {}
    for axis, size in enumerate(shape, start=1):
        lower = geometry[2 * (axis - 1)]
        upper = geometry[2 * (axis - 1) + 1]
        faces = np.linspace(lower, upper, size + 1, dtype=dtype)
        data[f"x{axis}f"] = faces
        data[f"x{axis}v"] = 0.5 * (faces[:-1] + faces[1:])
    for name in quantities:
        data[name] = np.array(
            filedata["mb_data"][name][meshblock_index], dtype=dtype, copy=True
        )
    data["Time"] = filedata["time"]
    data["NumCycles"] = filedata["cycle"]
    data["MaxLevel"] = int(filedata["mb_logical"][meshblock_index, 3])
    data["MeshBlockIndex"] = int(meshblock_index)
    data["LogicalLocation"] = filedata["mb_logical"][meshblock_index].copy()
    data["VariableNames"] = tuple(quantities)
    data["Bounds"] = np.array(
        [[data[f"x{i}f"][0], data[f"x{i}f"][-1]] for i in range(1, 4)],
        dtype=dtype,
    )
    data["DomainBounds"] = np.array(
        [
            [filedata["x1min"], filedata["x1max"]],
            [filedata["x2min"], filedata["x2max"]],
            [filedata["x3min"], filedata["x3max"]],
        ],
        dtype=dtype,
    )
    data["PeriodicAxes"] = tuple(filedata.get("periodic_axes", ()))
    data["SourceFile"] = source
    return data


def read_rank_meshblocks(
    filename: str | Path,
    quantities: Sequence[str] | None = None,
    dtype: np.dtype | type = np.float32,
) -> list[dict]:
    """Read every MeshBlock in one rank file while touching the file once."""

    filedata = bin_convert.read_binary_as_athdf(str(filename), raw=True)
    names = tuple(filedata["var_names"] if quantities is None else quantities)
    missing = sorted(set(names) - set(filedata["var_names"]))
    if missing:
        raise KeyError(f"MHD quantities not present in {filename}: {missing}")
    return [
        _meshblock_from_raw(filedata, index, names, dtype, str(filename))
        for index in range(filedata["n_mbs"])
    ]


def _parse_fields(text: str, path: Path) -> tuple[str, ...]:
    values = {key: value for key, value in KEY_RE.findall(text)}
    fields = tuple(values.get("fields", "").split(","))
    if not fields or int(values.get("nfields", 0)) != len(fields):
        raise ValueError(f"{path}: invalid tracked-particle field list")
    return fields


def _iter_legacy_frames(handle: BinaryIO, path: Path) -> Iterator[TrackFrame]:
    while True:
        line = handle.readline()
        while line and not line.strip():
            line = handle.readline()
        if not line:
            return
        if not line.startswith(LEGACY_MARKER):
            raise ValueError(f"{path}: invalid legacy frame at byte {handle.tell()}")

        lines = [line.decode("ascii")]
        while True:
            line = handle.readline()
            if not line:
                raise ValueError(f"{path}: unterminated tracked-particle header")
            if not line.strip():
                break
            lines.append(line.decode("ascii"))

        text = " ".join(lines)
        values = {key: value for key, value in KEY_RE.findall(text)}
        fields = _parse_fields(text, path)
        count = int(values["record_count"])
        payload_bytes = count * len(fields) * np.dtype("<f4").itemsize
        payload = np.frombuffer(
            _read_exact(handle, payload_bytes, path, "legacy payload"),
            dtype="<f4",
        )
        yield TrackFrame(
            cycle=int(values["cycle"]),
            time=float(values["time"]),
            fields=fields,
            values=payload.reshape(count, len(fields)),
        )


def _read_exact(handle: BinaryIO, size: int, path: Path, what: str) -> bytes:
    data = handle.read(size)
    if len(data) != size:
        raise ValueError(f"{path}: truncated {what}")
    return data


def _iter_compact_frames(handle: BinaryIO, path: Path) -> Iterator[TrackFrame]:
    raw = _read_exact(handle, COMPACT_PROLOGUE.size, path, "compact prologue")
    values = COMPACT_PROLOGUE.unpack(raw)
    magic, version, prologue_bytes, nfields = values[:4]
    fields_bytes = values[-2]
    if magic != COMPACT_FILE_MAGIC or version != COMPACT_VERSION:
        raise ValueError(f"{path}: invalid compact track prologue")
    if prologue_bytes != COMPACT_PROLOGUE.size:
        raise ValueError(f"{path}: unsupported compact prologue size")
    fields = tuple(
        _read_exact(handle, fields_bytes, path, "compact field list")
        .decode("ascii")
        .split(",")
    )
    if len(fields) != nfields:
        raise ValueError(f"{path}: compact field count mismatch")

    while True:
        raw = handle.read(COMPACT_FRAME.size)
        if not raw:
            return
        if len(raw) != COMPACT_FRAME.size:
            raise ValueError(f"{path}: truncated compact frame")
        (magic, version, frame_bytes, count, cycle, time,
         payload_bytes, _reserved) = COMPACT_FRAME.unpack(raw)
        expected = count * nfields * np.dtype("<f4").itemsize
        if magic != COMPACT_FRAME_MAGIC or version != COMPACT_VERSION:
            raise ValueError(f"{path}: invalid compact frame at cycle {cycle}")
        if frame_bytes != COMPACT_FRAME.size or payload_bytes != expected:
            raise ValueError(f"{path}: compact frame size mismatch at cycle {cycle}")
        payload = np.frombuffer(
            _read_exact(handle, payload_bytes, path, "compact payload"),
            dtype="<f4",
        ).reshape(count, nfields)
        yield TrackFrame(int(cycle), float(time), fields, payload)


def iter_track_frames(filename: str | Path) -> Iterator[TrackFrame]:
    """Stream legacy or compact AthenaK track frames with bounded memory."""

    path = Path(filename)
    with path.open("rb") as handle:
        magic = handle.read(len(COMPACT_FILE_MAGIC))
        handle.seek(0)
        if magic == COMPACT_FILE_MAGIC:
            yield from _iter_compact_frames(handle, path)
        else:
            yield from _iter_legacy_frames(handle, path)


def points_in_bounds(points: np.ndarray, bounds: np.ndarray) -> np.ndarray:
    """Return a mask for points inside three inclusive coordinate intervals."""

    bounds = np.asarray(bounds)
    scale = max(1.0, float(np.max(np.abs(bounds))))
    tolerance = 8.0 * np.finfo(points.dtype).eps * scale
    return np.all(
        (points >= bounds[:, 0] - tolerance)
        & (points <= bounds[:, 1] + tolerance),
        axis=1,
    )


def find_track_visits(
    filenames: str | Path | Sequence[str | Path],
    bounds: np.ndarray,
    time_min: float | None = None,
    time_max: float | None = None,
) -> TrackVisits:
    """Find particle output tags recorded inside a MeshBlock-sized box."""

    paths = (
        [Path(filenames)]
        if isinstance(filenames, (str, Path))
        else [Path(filename) for filename in filenames]
    )
    tags: set[int] = set()
    frames_seen = records_seen = records_inside = 0
    for path in paths:
        for frame in iter_track_frames(path):
            if time_min is not None and frame.time < time_min:
                continue
            if time_max is not None and frame.time > time_max:
                break
            field = {name: index for index, name in enumerate(frame.fields)}
            tag_name = "output_tag" if "output_tag" in field else "tag"
            required = {tag_name, "time", "x", "y", "z"}
            missing = sorted(required - set(field))
            if missing:
                raise ValueError(f"{path}: required track fields missing: {missing}")
            points = frame.values[:, [field["x"], field["y"], field["z"]]]
            inside = points_in_bounds(points, bounds)
            raw_tags = frame.values[inside, field[tag_name]]
            rounded = np.rint(raw_tags)
            if not np.allclose(raw_tags, rounded, rtol=0.0, atol=1.0e-4):
                raise ValueError(f"{path}: non-integral particle output tag")
            tags.update(rounded.astype(np.int64).tolist())
            frames_seen += 1
            records_seen += frame.values.shape[0]
            records_inside += int(np.count_nonzero(inside))
    return TrackVisits(
        np.array(sorted(tags), dtype=np.int64),
        frames_seen,
        records_seen,
        records_inside,
    )


def _text(value: object) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def _field_names(dataset: h5py.Dataset) -> tuple[str, ...]:
    raw = dataset.attrs.get("fields")
    if raw is None:
        raise ValueError(f"{dataset.name} does not define its fields attribute")
    if isinstance(raw, (str, bytes)):
        return tuple(_text(raw).split(","))
    return tuple(_text(item) for item in raw)


def _time_slice(
    times: np.ndarray,
    time_min: float | None,
    time_max: float | None,
    stride: int,
) -> slice:
    if stride < 1:
        raise ValueError("time_stride must be at least one")
    start = 0 if time_min is None else int(np.searchsorted(times, time_min))
    stop = len(times) if time_max is None else int(
        np.searchsorted(times, time_max, side="right")
    )
    return slice(start, stop, stride)


def _selected_fields(
    available: tuple[str, ...], requested: Sequence[str] | None
) -> tuple[tuple[str, ...], np.ndarray]:
    names = available if requested is None else tuple(requested)
    missing = sorted(set(names) - set(available))
    if missing:
        raise KeyError(f"track fields not present: {missing}")
    return names, np.array([available.index(name) for name in names], dtype=int)


def read_merged_track_subset(
    filename: str | Path,
    output_tags: Sequence[int],
    bounds: np.ndarray | None = None,
    species: Sequence[int] | None = None,
    max_particles: int | None = None,
    seed: int = 1,
    time_min: float | None = None,
    time_max: float | None = None,
    time_stride: int = 1,
    fields: Sequence[str] | None = None,
    inside_only: bool = False,
) -> dict:
    """Read complete merged trajectories for particles selected by output tag."""

    requested_tags = np.asarray(output_tags, dtype=np.int64)
    with h5py.File(filename, "r") as handle:
        particles_all = handle["particles"][:]
        names = particles_all.dtype.names or ()
        if "output_tag" not in names:
            raise ValueError("merged particles dataset has no output_tag field")
        selected = np.isin(particles_all["output_tag"], requested_tags)
        if species is not None:
            if "species" not in names:
                raise ValueError("merged particles dataset has no species field")
            selected &= np.isin(particles_all["species"], np.asarray(species))
        source_rows = np.flatnonzero(selected)
        if not source_rows.size:
            raise ValueError("no merged particles match the requested selection")
        if max_particles is not None and source_rows.size > max_particles:
            if max_particles < 1:
                raise ValueError("max_particles must be at least one")
            source_rows = np.sort(
                np.random.default_rng(seed).choice(
                    source_rows, size=max_particles, replace=False
                )
            )

        values_ds = handle["values"]
        available = _field_names(values_ds)
        selected_names, selected_indices = _selected_fields(available, fields)
        coordinate_indices = [available.index(axis) for axis in ("x", "y", "z")]
        times_all = handle["times"][:]
        time_slice = _time_slice(times_all, time_min, time_max, time_stride)
        times = times_all[time_slice]
        cycles = handle["cycles"][time_slice]
        values = np.empty(
            (source_rows.size, times.size, len(selected_names)), dtype=values_ds.dtype
        )
        inside = np.ones((source_rows.size, times.size), dtype=bool)
        for output_row, source_row in enumerate(source_rows):
            complete = values_ds[
                source_row, time_slice.start:time_slice.stop, :
            ][::time_slice.step]
            values[output_row] = complete[:, selected_indices]
            if bounds is not None:
                inside[output_row] = points_in_bounds(
                    complete[:, coordinate_indices], bounds
                )
        if inside_only:
            values[~inside] = np.nan
        root_attrs = {name: handle.attrs[name] for name in handle.attrs}

    return {
        "particles": particles_all[source_rows],
        "source_rows": source_rows,
        "times": times,
        "cycles": cycles,
        "values": values,
        "fields": selected_names,
        "inside_meshblock": inside,
        "source_file": str(filename),
        "source_attrs": root_attrs,
    }


def read_merged_track_partition(
    filename: str | Path,
    partition: int,
    num_partitions: int,
    time_min: float | None = None,
    time_max: float | None = None,
    time_stride: int = 1,
    fields: Sequence[str] | None = None,
    particle_batch: int = 4,
) -> dict:
    """Read one contiguous particle-row partition from a merged track file."""

    with h5py.File(filename, "r") as handle:
        nrows = handle["values"].shape[0]
        row_start = partition * nrows // num_partitions
        row_stop = (partition + 1) * nrows // num_partitions
        times_all = handle["times"][:]
        time_slice = _time_slice(times_all, time_min, time_max, time_stride)
        values_ds = handle["values"]
        available = _field_names(values_ds)
        selected_names, selected_indices = _selected_fields(available, fields)
        if particle_batch < 1:
            raise ValueError("particle_batch must be at least one")
        values = np.empty(
            (row_stop - row_start, times_all[time_slice].size, len(selected_names)),
            dtype=values_ds.dtype,
        )
        for begin in range(row_start, row_stop, particle_batch):
            end = min(begin + particle_batch, row_stop)
            complete = values_ds[
                begin:end, time_slice.start:time_slice.stop, :
            ][:, ::time_slice.step, :]
            output = complete if fields is None else complete[..., selected_indices]
            values[begin - row_start:end - row_start] = output
        return {
            "particles": handle["particles"][row_start:row_stop],
            "source_rows": np.arange(row_start, row_stop, dtype=np.int64),
            "times": times_all[time_slice],
            "cycles": handle["cycles"][time_slice],
            "values": values,
            "fields": selected_names,
            "source_file": str(filename),
        }


def discover_rank_files(rank0_filename: str | Path) -> list[Path]:
    """Discover matching per-rank files from one rank-zero path."""

    rank0 = Path(rank0_filename)
    layout = re.fullmatch(r"(rank|node)_\d+", rank0.parent.name)
    if layout is None:
        return [rank0]
    files = sorted(
        rank0.parent.parent.glob(f"{layout.group(1)}_*/{rank0.name}")
    )
    if not files:
        raise FileNotFoundError(f"no rank files match {rank0}")
    return files


def write_meshblock_bundle(
    filename: str | Path, meshblock: dict, tracks: dict
) -> None:
    """Write one convenient HDF5 bundle containing a block and its tracks."""

    with h5py.File(filename, "x") as handle:
        handle.attrs["format"] = "athenak_meshblock_tracks_v1"
        mhd = handle.create_group("mhd")
        mhd.attrs["source_file"] = meshblock["SourceFile"]
        mhd.attrs["time"] = meshblock["Time"]
        mhd.attrs["cycle"] = meshblock["NumCycles"]
        mhd.attrs["meshblock_index"] = meshblock["MeshBlockIndex"]
        mhd.attrs["logical_location"] = meshblock["LogicalLocation"]
        mhd.attrs["periodic_axes"] = meshblock.get("PeriodicAxes", ())
        mhd.attrs["fields"] = ",".join(meshblock["VariableNames"])
        mhd.create_dataset("bounds", data=meshblock["Bounds"])
        mhd.create_dataset(
            "domain_bounds", data=meshblock.get("DomainBounds", meshblock["Bounds"])
        )
        for name in ("x1f", "x1v", "x2f", "x2v", "x3f", "x3v"):
            mhd.create_dataset(name, data=meshblock[name])
        for name in meshblock["VariableNames"]:
            mhd.create_dataset(name, data=meshblock[name])

        trk = handle.create_group("tracks")
        for name, value in tracks.get("source_attrs", {}).items():
            trk.attrs[name] = value
        trk.attrs["source_file"] = tracks["source_file"]
        trk.create_dataset("particles", data=tracks["particles"])
        trk.create_dataset("source_rows", data=tracks["source_rows"])
        trk.create_dataset("times", data=tracks["times"])
        trk.create_dataset("cycles", data=tracks["cycles"])
        values = trk.create_dataset("values", data=tracks["values"])
        values.attrs["fields"] = ",".join(tracks["fields"])
        trk.create_dataset("inside_meshblock", data=tracks["inside_meshblock"])
