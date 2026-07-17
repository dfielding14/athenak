"""Focused tests for the MHD and particle visualization readers."""

from pathlib import Path
import struct
import sys

import h5py
import numpy as np


PYTHON_VIS = Path(__file__).resolve().parents[1]
if str(PYTHON_VIS) not in sys.path:
    sys.path.insert(0, str(PYTHON_VIS))

import bin_convert  # noqa: E402
from cr_visualization import cr_data  # noqa: E402


TRACK_FIELDS = ("tag", "time", "x", "y", "z", "temperature")


def write_binary(path: Path, block_ids: tuple[int, ...] = (0, 1)) -> None:
    header = """<job>
problem = synthetic MHD with eta=3e-6
<mesh>
nx1 = 4
nx2 = 2
nx3 = 2
nghost = 0
x1min = 0
x1max = 1
x2min = 0
x2max = 1
x3min = 0
x3max = 1
<meshblock>
nx1 = 2
nx2 = 2
nx3 = 2
""".encode()
    with path.open("wb") as handle:
        handle.write(b"Athena binary output version=1.1\n")
        handle.write(b"  pheader_count=5\n")
        handle.write(b"  time=0.5\n")
        handle.write(b"  cycle=7\n")
        handle.write(b"  size of location=8\n")
        handle.write(b"  size of variable=4\n")
        handle.write(b"  number of variables=2\n")
        handle.write(b"  variables: dens bcc1\n")
        handle.write(f"  header offset={len(header)}\n".encode())
        handle.write(header)
        for block in block_ids:
            handle.write(struct.pack("<6i", 0, 1, 0, 1, 0, 1))
            handle.write(struct.pack("<4i", block, 0, 0, 0))
            handle.write(
                struct.pack("<6d", 0.5 * block, 0.5 * (block + 1), 0, 1, 0, 1)
            )
            values = np.arange(16, dtype="<f4").reshape(2, 2, 2, 2)
            (values + 100 * block).tofile(handle)


def track_values() -> np.ndarray:
    return np.array(
        [
            [0, 0, 0.10, 0.20, 0.30, 10],
            [1, 0, 0.60, 0.20, 0.30, 11],
            [2, 0, 0.90, 0.80, 0.70, 12],
        ],
        dtype="<f4",
    )


def write_legacy_track(path: Path) -> None:
    values = track_values()
    with path.open("wb") as handle:
        handle.write(
            b"# AthenaK tracked particle data at time=0 cycle=11 "
            b"record_count=3 nfields=6\n"
        )
        handle.write(b"# trk_format=custom_v1 layout=rank rank=0\n")
        handle.write(b"# fields=" + ",".join(TRACK_FIELDS).encode() + b"\n\n")
        values.tofile(handle)


def write_compact_track(path: Path) -> None:
    values = track_values()
    fields = ",".join(TRACK_FIELDS).encode()
    with path.open("wb") as handle:
        handle.write(
            cr_data.COMPACT_PROLOGUE.pack(
                cr_data.COMPACT_FILE_MAGIC,
                cr_data.COMPACT_VERSION,
                cr_data.COMPACT_PROLOGUE.size,
                len(TRACK_FIELDS),
                3,
                3,
                0,
                2,
                0,
                0,
                1,
                1,
                1,
                len(fields),
                0,
            )
        )
        handle.write(fields)
        payload_bytes = values.nbytes
        handle.write(
            cr_data.COMPACT_FRAME.pack(
                cr_data.COMPACT_FRAME_MAGIC,
                cr_data.COMPACT_VERSION,
                cr_data.COMPACT_FRAME.size,
                len(values),
                11,
                0.0,
                payload_bytes,
                0,
            )
        )
        values.tofile(handle)


def write_merged_tracks(path: Path) -> None:
    particle_dtype = np.dtype([("output_tag", "<i8"), ("species", "<i4")])
    particles = np.array([(0, 0), (1, 1), (2, 1)], dtype=particle_dtype)
    times = np.array([0.0, 1.0, 2.0])
    values = np.zeros((3, 3, 4), dtype="<f4")
    for particle in range(3):
        values[particle, :, 0] = [0.1 + 0.3 * particle, 0.6, 0.9]
        values[particle, :, 1] = 0.2 + 0.1 * particle
        values[particle, :, 2] = 0.3
        values[particle, :, 3] = 10 + particle
    with h5py.File(path, "w") as handle:
        handle.attrs["format"] = "test_merged_tracks"
        handle.create_dataset("particles", data=particles)
        handle.create_dataset("times", data=times)
        handle.create_dataset("cycles", data=np.array([11, 12, 13]))
        dataset = handle.create_dataset("values", data=values)
        dataset.attrs["fields"] = "x,y,z,temperature"


def test_binary_header_and_meshblock_selection(tmp_path: Path) -> None:
    binary = tmp_path / "test.bin"
    write_binary(binary)

    block = bin_convert.read_single_rank_binary_as_athdf(
        binary, meshblock_index=1, quantities=["dens"]
    )
    assert block["dens"].shape == (2, 2, 2)
    assert np.allclose(block["Bounds"], [[0.5, 1], [0, 1], [0, 1]])
    assert tuple(block["LogicalLocation"]) == (1, 0, 0, 0)
    assert block["dens"][0, 0, 0] == 100

    cropped = bin_convert.read_single_rank_binary_as_athdf(
        binary,
        meshblock_index=1,
        quantities=["dens"],
        x1_min=0.75,
    )
    assert cropped["dens"].shape == (2, 2, 1)
    assert np.allclose(cropped["x1f"], [0.75, 1.0])


def test_rank_reader_handles_multiple_meshblocks(tmp_path: Path) -> None:
    binary = tmp_path / "test.bin"
    write_binary(binary)
    blocks = cr_data.read_rank_meshblocks(binary, quantities=["dens"])
    assert len(blocks) == 2
    assert blocks[0]["dens"][0, 0, 0] == 0
    assert blocks[1]["dens"][0, 0, 0] == 100


def test_distributed_rank_parts_match_serial_assembly(tmp_path: Path) -> None:
    rank0 = tmp_path / "rank_00000000" / "snapshot.bin"
    rank1 = tmp_path / "rank_00000001" / "snapshot.bin"
    rank0.parent.mkdir()
    rank1.parent.mkdir()
    write_binary(rank0, block_ids=(0,))
    write_binary(rank1, block_ids=(1,))

    files = cr_data.discover_rank_files(rank0)
    blocks = [
        block
        for filename in files
        for block in cr_data.read_rank_meshblocks(filename, quantities=["dens"])
    ]
    serial = bin_convert.read_all_ranks_binary_as_athdf(
        rank0, quantities=["dens"]
    )
    assert len(files) == 2
    assert len(blocks) == 2
    assert sum(block["dens"].sum() for block in blocks) == serial["dens"].sum()


def test_legacy_and_compact_track_visits(tmp_path: Path) -> None:
    bounds = np.array([[0.5, 1.0], [0, 1], [0, 1]], dtype=np.float32)
    for writer in (write_legacy_track, write_compact_track):
        track = tmp_path / f"{writer.__name__}.trk"
        writer(track)
        visits = cr_data.find_track_visits(track, bounds)
        assert np.array_equal(visits.output_tags, [1, 2])
        assert visits.frames_seen == 1
        assert visits.records_seen == 3
        assert visits.records_inside == 2


def test_merged_subset_partition_and_bundle(tmp_path: Path) -> None:
    binary = tmp_path / "test.bin"
    merged = tmp_path / "tracks.h5"
    bundle = tmp_path / "bundle.h5"
    write_binary(binary)
    write_merged_tracks(merged)
    bounds = np.array([[0.5, 1.0], [0, 1], [0, 1]], dtype=np.float32)

    subset = cr_data.read_merged_track_subset(
        merged,
        [1, 2],
        bounds=bounds,
        species=[1],
        time_stride=2,
        fields=["x", "temperature"],
        inside_only=True,
    )
    assert subset["values"].shape == (2, 2, 2)
    assert np.array_equal(subset["particles"]["output_tag"], [1, 2])
    assert np.isnan(subset["values"][0, 0, 0])
    assert subset["values"][0, 1, 0] == np.float32(0.9)

    first = cr_data.read_merged_track_partition(merged, 0, 2)
    second = cr_data.read_merged_track_partition(merged, 1, 2)
    assert first["values"].shape[0] + second["values"].shape[0] == 3
    with h5py.File(merged, "r") as handle:
        assert np.array_equal(
            np.concatenate([first["values"], second["values"]]),
            handle["values"][:],
        )

    block = cr_data.read_meshblock(binary, 1, quantities=["dens"])
    cr_data.write_meshblock_bundle(bundle, block, subset)
    with h5py.File(bundle, "r") as handle:
        assert handle.attrs["format"] == "athenak_meshblock_tracks_v1"
        assert handle["mhd/dens"].shape == (2, 2, 2)
        assert handle["tracks"].attrs["format"] == "test_merged_tracks"
        assert handle["tracks/values"].shape == (2, 2, 2)
