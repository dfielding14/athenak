"""Pure-Python regression coverage for binary, PDF, and spherical-slice readers."""

from pathlib import Path
import shutil
import struct
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tst" / "fixtures" / "io" / "origin_main_886dd2a1"
sys.path.insert(0, str(ROOT / "vis" / "python"))

from bin_convert import (  # noqa: E402
    convert_file,
    read_all_ranks_binary,
    read_all_ranks_coarsened_binary,
    read_binary,
    read_coarsened_binary,
    write_athdf,
    write_xdmf_for,
)
from read_pdf import read_pdf  # noqa: E402
from read_sphslice import read_sphslice, read_sphslice_header  # noqa: E402


def _binary_fixture(kind, layout, dump, rank=None):
    if layout == "shared":
        name = f"io_legacy_shared.{kind}_shared.{dump}.{kind}"
        return FIXTURES / kind / layout / name
    name = f"io_legacy_per_rank.{kind}_per_rank.{dump}.{kind}"
    return FIXTURES / kind / layout / f"rank_{rank:08d}" / name


def _join_x1_blocks(filedata, variable):
    order = np.argsort(filedata["mb_logical"][:, 0])
    blocks = [filedata["mb_data"][variable][item] for item in order]
    return np.concatenate(blocks, axis=-1)


@pytest.mark.parametrize(
    ("kind", "reader", "all_ranks_reader"),
    (
        ("bin", read_binary, read_all_ranks_binary),
        ("cbin", read_coarsened_binary, read_all_ranks_coarsened_binary),
    ),
)
@pytest.mark.parametrize("dump", ("00000", "00001"))
def test_legacy_binary_fixtures_read_and_rank_assembly_matches_shared(
    kind, reader, all_ranks_reader, dump
):
    shared = reader(str(_binary_fixture(kind, "shared", dump)))
    shard_paths = [_binary_fixture(kind, "per_rank", dump, rank) for rank in (0, 1)]
    shards = [reader(str(path)) for path in shard_paths]
    assembled = reader(str(shard_paths[0]), assemble_shards=True)
    assembled_legacy_api = all_ranks_reader(str(shard_paths[0]))

    assert shared["time"] == assembled["time"]
    assert shared["cycle"] == assembled["cycle"]
    assert shared["var_names"] == assembled["var_names"]
    assert assembled["n_mbs"] == sum(shard["n_mbs"] for shard in shards) == 2
    np.testing.assert_array_equal(
        assembled_legacy_api["mb_logical"], assembled["mb_logical"]
    )

    for variable in shared["var_names"]:
        expected = shared["mb_data"][variable][0]
        np.testing.assert_allclose(_join_x1_blocks(assembled, variable), expected)
        np.testing.assert_allclose(
            assembled_legacy_api["mb_data"][variable],
            assembled["mb_data"][variable],
        )


@pytest.mark.parametrize("kind", ("bin", "cbin"))
def test_conversion_helpers_remain_public_and_convert_legacy_binary(tmp_path, kind):
    assert callable(write_athdf)
    assert callable(write_xdmf_for)
    assert callable(convert_file)

    source = _binary_fixture(kind, "shared", "00000")
    binary = tmp_path / source.name
    shutil.copyfile(source, binary)

    convert_file(str(binary))

    athdf = binary.with_suffix(".athdf")
    xdmf = Path(str(athdf) + ".xdmf")
    assert athdf.exists()
    assert xdmf.exists()
    assert athdf.name in xdmf.read_text()


def test_canonical_converter_cli_converts_shared_and_rank_shards(tmp_path):
    shared = tmp_path / _binary_fixture("bin", "shared", "00000").name
    shutil.copyfile(_binary_fixture("bin", "shared", "00000"), shared)
    rank_root = tmp_path / "ranked"
    for rank in (0, 1):
        target = rank_root / f"rank_{rank:08d}" / _binary_fixture(
            "bin", "per_rank", "00000", rank=rank
        ).name
        target.parent.mkdir(parents=True)
        shutil.copyfile(_binary_fixture("bin", "per_rank", "00000", rank=rank), target)
    rank0 = rank_root / "rank_00000000" / _binary_fixture(
        "bin", "per_rank", "00000", rank=0
    ).name

    subprocess.run(
        [sys.executable, str(ROOT / "vis" / "python" / "bin_convert.py"), str(shared)],
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "vis" / "python" / "bin_convert.py"),
            "--assemble-shards",
            str(rank0),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert shared.with_suffix(".athdf").exists()
    assert Path(str(shared.with_suffix(".athdf")) + ".xdmf").exists()
    assert rank0.with_suffix(".athdf").exists()


def test_repository_exposes_only_canonical_converter_imports():
    forbidden = "bin_convert_" + "new"
    assert not (ROOT / "vis" / "python" / (forbidden + ".py")).exists()
    for directory in ("vis", "inputs", "tst", "src"):
        for path in (ROOT / directory).rglob("*"):
            if path.is_file() and path.suffix in {".py", ".md", ".athinput"}:
                assert forbidden not in path.read_text(errors="ignore")


@pytest.mark.parametrize(
    ("dump", "time", "expected"),
    (
        ("00000", 0.0, [0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0]),
        (
            "00001",
            0.01,
            [0.0, 0.5, 0.0, 0.0, 0.125, 0.375, 0.0, 0.0, 0.0, 0.0],
        ),
    ),
)
def test_legacy_pdf_fixtures_are_read(dump, time, expected):
    source = FIXTURES / "pdf" / "legacy_1d" / f"io_legacy_shared.{dump}.pdf"
    result = read_pdf(str(source))

    assert result["header"]["format"] == "legacy_dense"
    assert result["header"]["dimensions"][0]["variable"] == "dens"
    assert result["time"] == time
    np.testing.assert_allclose(result["pdf"], expected)


def _write_pdf_header(path, fmt, distribution=None, v2=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"format = {fmt}",
        "ndim = 1",
        "variable_1 = dens",
        "nbin1 = 2",
        "bin1_min = 0.0",
        "bin1_max = 2.0",
        "scale1 = linear",
        "stride1 = 1",
        "total_bins = 4",
        "cycle = 5",
    ]
    if v2:
        lines.insert(1, "binary_magic = AKPDFV2")
    if distribution is not None:
        lines.insert(1, f"distribution = {distribution}")
    path.write_text("\n".join(lines) + "\n")


def _write_dense_pdf(path, values, time=0.25):
    path.write_bytes(np.asarray([time, *values], dtype=np.float64).tobytes())


def _write_sparse_pdf(path, entries, time=0.25):
    path.parent.mkdir(parents=True, exist_ok=True)
    indices = np.asarray([item[0] for item in entries], dtype=np.uint32)
    values = np.asarray([item[1] for item in entries], dtype=np.float64)
    payload = (
        struct.pack("=dI", time, len(entries))
        + indices.tobytes()
        + values.tobytes()
    )
    path.write_bytes(payload)


def _write_v2_sparse_pdf(path, entries, time=0.25, cycle=5, declared_count=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    count = len(entries) if declared_count is None else declared_count
    payload = struct.pack("=8sIIIIQdq", b"AKPDFV2\0", 2, 1, 1, 0, count, time, cycle)
    payload += b"".join(struct.pack("=Qd", index, value) for index, value in entries)
    path.write_bytes(payload)


def test_modern_pdf_sparse_empty_shard_matches_dense_shared(tmp_path):
    expected = np.array([1.0, 0.0, 0.0, 2.5])
    dense = tmp_path / "dense.00000.pdf"
    _write_pdf_header(tmp_path / "dense.header.pdf", "dense")
    _write_dense_pdf(dense, expected)

    sparse = tmp_path / "sparse"
    _write_pdf_header(sparse / "hist.header.pdf", "sparse_coo", "rank")
    rank0 = sparse / "rank_00000000" / "hist.00000.pdf"
    rank1 = sparse / "rank_00000001" / "hist.00000.pdf"
    _write_sparse_pdf(rank0, [(0, 1.0), (3, 2.5)])
    _write_sparse_pdf(rank1, [])

    np.testing.assert_allclose(read_pdf(str(rank0))["pdf"], read_pdf(str(dense))["pdf"])


def test_modern_pdf_rejects_truncated_sparse_payload(tmp_path):
    case = tmp_path / "truncated"
    _write_pdf_header(case / "hist.header.pdf", "sparse_coo", "rank", v2=True)
    shard = case / "rank_00000000" / "hist.00000.pdf"
    _write_v2_sparse_pdf(shard, [], declared_count=1)

    with pytest.raises(ValueError, match="sparse PDF shard .* truncated"):
        read_pdf(str(shard))


def test_modern_pdf_rejects_duplicate_sparse_indices(tmp_path):
    case = tmp_path / "duplicate"
    _write_pdf_header(case / "hist.header.pdf", "sparse_coo", "rank", v2=True)
    shard = case / "rank_00000000" / "hist.00000.pdf"
    _write_v2_sparse_pdf(shard, [(1, 2.0), (1, 3.0)])

    with pytest.raises(ValueError, match="duplicate COO indices"):
        read_pdf(str(shard))


def test_modern_pdf_rejects_malformed_header(tmp_path):
    data = tmp_path / "bad.00000.pdf"
    header = tmp_path / "bad.header.pdf"
    header.write_text("format = dense\nndim = 1\n")
    _write_dense_pdf(data, [0.0, 0.0, 0.0, 0.0])

    with pytest.raises(ValueError, match="missing dimension"):
        read_pdf(str(data))


def _write_sphslice(path, distribution, indices, values, npoints=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if npoints is None:
        npoints = 2 if distribution == "shared" else len(indices)
    lines = [
        "Athena spherical slice version=1.0",
        "time=0.25",
        "cycle=5",
        "radius=2.0",
        "ntheta=1",
        "nphi=2",
        "size_of_variable=4",
        "number_of_variables=1",
        f"npoints={npoints}",
        f"distribution={distribution}",
        "variables: dens",
        "header_offset=0",
    ]
    payload = b"\n".join(line.encode("ascii") for line in lines) + b"\n"
    if distribution == "shared":
        payload += np.asarray(values, dtype=np.float32).tobytes()
    else:
        payload += np.asarray(indices, dtype=np.int32).tobytes()
        payload += np.asarray(values, dtype=np.float32).tobytes()
    path.write_bytes(payload)


def test_sphslice_empty_shard_reassembles_to_shared_surface(tmp_path):
    shared = tmp_path / "surface.00000.sphslice"
    _write_sphslice(shared, "shared", [], [10.0, 20.0])

    shards = tmp_path / "shards"
    rank0 = shards / "rank_00000000" / "surface.00000.sphslice"
    rank1 = shards / "rank_00000001" / "surface.00000.sphslice"
    _write_sphslice(rank0, "rank", [0, 1], [10.0, 20.0])
    _write_sphslice(rank1, "rank", [], [])

    assert read_sphslice_header(str(rank1))["npoints"] == 0
    np.testing.assert_allclose(
        read_sphslice(str(rank0))["data"], read_sphslice(str(shared))["data"]
    )


def test_sphslice_rejects_truncated_shard_payload(tmp_path):
    shard = tmp_path / "truncated" / "rank_00000000" / "surface.00000.sphslice"
    _write_sphslice(shard, "rank", [0], [], npoints=1)

    with pytest.raises(ValueError, match="truncated or overlong spherical-slice shard"):
        read_sphslice(str(shard))


def test_sphslice_rejects_duplicate_ownership_across_shards(tmp_path):
    root = tmp_path / "duplicate"
    rank0 = root / "rank_00000000" / "surface.00000.sphslice"
    rank1 = root / "rank_00000001" / "surface.00000.sphslice"
    _write_sphslice(rank0, "rank", [0], [10.0])
    _write_sphslice(rank1, "rank", [0], [20.0])

    with pytest.raises(ValueError, match="duplicate ownership"):
        read_sphslice(str(rank0))


def test_sphslice_rejects_missing_angular_ownership(tmp_path):
    root = tmp_path / "missing"
    rank0 = root / "rank_00000000" / "surface.00000.sphslice"
    rank1 = root / "rank_00000001" / "surface.00000.sphslice"
    _write_sphslice(rank0, "rank", [0], [10.0])
    _write_sphslice(rank1, "rank", [], [])

    with pytest.raises(ValueError, match="missing 1 of 2 angular points"):
        read_sphslice(str(rank0))
