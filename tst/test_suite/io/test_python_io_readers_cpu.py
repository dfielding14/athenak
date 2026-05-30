"""Pure-Python regression coverage for binary, PDF, and spherical-slice readers."""

import builtins
import hashlib
import inspect
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import weakref

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tst" / "fixtures" / "io" / "origin_main_886dd2a1"
sys.path.insert(0, str(ROOT / "vis" / "python"))

import bin_convert  # noqa: E402
import read_pdf as read_pdf_module  # noqa: E402
import read_sphslice as read_sphslice_module  # noqa: E402
from bin_convert import (  # noqa: E402
    convert_file,
    read_all_ranks_binary,
    read_all_ranks_coarsened_binary,
    read_binary,
    read_coarsened_binary,
    write_athdf,
    write_xdmf_for,
)
from read_pdf import read_pdf, read_pdf_header  # noqa: E402
from read_sphslice import read_sphslice, read_sphslice_header  # noqa: E402


def _binary_fixture(kind, layout, dump, rank=None):
    if layout == "shared":
        name = f"io_legacy_shared.{kind}_shared.{dump}.{kind}"
        return FIXTURES / kind / layout / name
    name = f"io_legacy_per_rank.{kind}_per_rank.{dump}.{kind}"
    return FIXTURES / kind / layout / f"rank_{rank:08d}" / name


def test_frozen_fixture_manifest_and_checksums_match():
    metadata_files = {"MANIFEST.tsv", "README.md", "SHA256SUMS"}
    fixture_artifacts = {
        str(path.relative_to(FIXTURES))
        for path in FIXTURES.rglob("*")
        if path.is_file() and str(path.relative_to(FIXTURES)) not in metadata_files
    }
    manifest_lines = (FIXTURES / "MANIFEST.tsv").read_text().splitlines()
    assert manifest_lines[0] == "path\tbytes\tsha256"
    manifest = {}
    for line in manifest_lines[1:]:
        relative_path, size, digest = line.split("\t")
        assert relative_path not in manifest
        manifest[relative_path] = (int(size), digest)

    checksum_entries = {}
    for line in (FIXTURES / "SHA256SUMS").read_text().splitlines():
        digest, relative_path = line.split(maxsplit=1)
        assert relative_path not in checksum_entries
        checksum_entries[relative_path] = digest

    assert set(manifest) == fixture_artifacts
    assert set(checksum_entries) == set(manifest)
    for relative_path, (expected_size, expected_digest) in manifest.items():
        path = FIXTURES / relative_path
        assert path.stat().st_size == expected_size
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_digest
        assert checksum_entries[relative_path] == expected_digest


def _join_x1_blocks(filedata, variable):
    order = np.argsort(filedata["mb_logical"][:, 0])
    blocks = [filedata["mb_data"][variable][item] for item in order]
    return np.concatenate(blocks, axis=-1)


def _binary_layout(path):
    with path.open("rb") as fp:
        fp.readline()
        pheader_count = int(fp.readline().split(b"=")[-1])
        pheader = {}
        for _ in range(pheader_count - 1):
            key, value = [
                item.strip() for item in fp.readline().decode("utf-8").split("=")
            ]
            pheader[key] = value
        nvars = int(fp.readline().split(b"=")[-1])
        fp.readline()
        header_size = int(fp.readline().split(b"=")[-1])
        fp.seek(header_size, 1)
        return (
            fp.tell(),
            nvars,
            int(pheader["size of location"]),
            int(pheader["size of variable"]),
        )


def _copy_binary_shards(tmp_path, kind, sources):
    paths = []
    for rank, source in enumerate(sources):
        path = tmp_path / f"rank_{rank:08d}" / f"output.00000.{kind}"
        path.parent.mkdir(parents=True)
        shutil.copyfile(source, path)
        paths.append(path)
    return paths


def _write_empty_binary(path, source):
    payload_offset, _, _, _ = _binary_layout(source)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(source.read_bytes()[:payload_offset])


def _write_two_meshblock_binary(path, kind="bin"):
    rank0 = _binary_fixture(kind, "per_rank", "00000", rank=0)
    rank1 = _binary_fixture(kind, "per_rank", "00000", rank=1)
    payload_offset, _, _, _ = _binary_layout(rank1)
    path.write_bytes(rank0.read_bytes() + rank1.read_bytes()[payload_offset:])


def _write_duplicate_meshblock_binary(path):
    rank0 = _binary_fixture("bin", "per_rank", "00000", rank=0)
    payload_offset, _, _, _ = _binary_layout(rank0)
    path.write_bytes(rank0.read_bytes() + rank0.read_bytes()[payload_offset:])


def _shrink_first_meshblock_x1(path):
    payload_offset, nvars, locsizebytes, varsizebytes = _binary_layout(path)
    payload = bytearray(path.read_bytes())
    meshblock_index = list(struct.unpack_from("=6i", payload, payload_offset))
    assert meshblock_index[1] > meshblock_index[0]
    meshblock_index[1] -= 1
    struct.pack_into("=6i", payload, payload_offset, *meshblock_index)

    shape = (
        meshblock_index[1] - meshblock_index[0] + 1,
        meshblock_index[3] - meshblock_index[2] + 1,
        meshblock_index[5] - meshblock_index[4] + 1,
    )
    values_offset = payload_offset + 24 + 16 + 6 * locsizebytes
    values_size = nvars * shape[0] * shape[1] * shape[2] * varsizebytes
    path.write_bytes(payload[: values_offset + values_size])


def _shrink_second_meshblock_x1(path):
    payload_offset, nvars, locsizebytes, varsizebytes = _binary_layout(path)
    payload = bytearray(path.read_bytes())
    first_index = struct.unpack_from("=6i", payload, payload_offset)
    first_shape = (
        first_index[1] - first_index[0] + 1,
        first_index[3] - first_index[2] + 1,
        first_index[5] - first_index[4] + 1,
    )
    fixed_bytes = 24 + 16 + 6 * locsizebytes
    first_value_bytes = nvars * np.prod(first_shape) * varsizebytes
    second_offset = payload_offset + fixed_bytes + first_value_bytes
    meshblock_index = list(struct.unpack_from("=6i", payload, second_offset))
    assert meshblock_index[1] > meshblock_index[0]
    meshblock_index[1] -= 1
    struct.pack_into("=6i", payload, second_offset, *meshblock_index)

    shape = (
        meshblock_index[1] - meshblock_index[0] + 1,
        meshblock_index[3] - meshblock_index[2] + 1,
        meshblock_index[5] - meshblock_index[4] + 1,
    )
    values_offset = second_offset + fixed_bytes
    values_size = nvars * np.prod(shape) * varsizebytes
    path.write_bytes(payload[: values_offset + values_size])


def _oversize_first_meshblock_x1(path):
    payload_offset, _, _, _ = _binary_layout(path)
    payload = bytearray(path.read_bytes())
    meshblock_index = list(struct.unpack_from("=6i", payload, payload_offset))
    meshblock_index[1] = 10**9
    struct.pack_into("=6i", payload, payload_offset, *meshblock_index)
    path.write_bytes(payload)


def _oversize_parameter_header(path):
    payload = path.read_bytes()
    start = payload.index(b"  header offset=")
    end = payload.index(b"\n", start) + 1
    path.write_bytes(payload[:start] + b"  header offset=999999999\n" + payload[end:])


def _replace_first_meshblock_logical(path, *, i=None, level=None):
    payload_offset, _, _, _ = _binary_layout(path)
    payload = bytearray(path.read_bytes())
    logical_offset = payload_offset + 24
    logical = list(struct.unpack_from("=4i", payload, logical_offset))
    if i is not None:
        logical[0] = i
    if level is not None:
        logical[3] = level
    struct.pack_into("=4i", payload, logical_offset, *logical)
    path.write_bytes(payload)


def _replace_binary_line(path, prefix, replacement):
    payload = path.read_bytes()
    start = payload.index(prefix)
    end = payload.index(b"\n", start) + 1
    path.write_bytes(payload[:start] + replacement + b"\n" + payload[end:])


def _replace_binary_bytes_same_width(path, before, after):
    payload = path.read_bytes()
    assert payload.count(before) == 1
    assert len(after) <= len(before)
    path.write_bytes(payload.replace(before, after.ljust(len(before)), 1))


def _replace_first_meshblock_geometry(
    path, *, x1min=None, x1max=None, x2max=None, x3max=None
):
    payload_offset, _, locsizebytes, _ = _binary_layout(path)
    assert locsizebytes == 8
    payload = bytearray(path.read_bytes())
    geometry_offset = payload_offset + 24 + 16
    geometry = list(struct.unpack_from("=6d", payload, geometry_offset))
    for index, value in ((0, x1min), (1, x1max), (3, x2max), (5, x3max)):
        if value is not None:
            geometry[index] = value
    struct.pack_into("=6d", payload, geometry_offset, *geometry)
    path.write_bytes(payload)


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
    assert assembled["rank"] is None
    assert assembled["node"] is None
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


def test_binary_converter_public_api_and_signature_contracts():
    assert "read_rank_binary_as_athdf" in bin_convert.__all__
    assert not hasattr(bin_convert, "athinput")
    assert inspect.signature(bin_convert.read_rank_binary_as_athdf) == inspect.signature(
        bin_convert.read_binary_as_athdf
    )

    signature = inspect.signature(bin_convert.read_single_rank_binary_as_athdf)
    positional_names = [
        "filename",
        "raw",
        "data",
        "quantities",
        "dtype",
        "return_levels",
        "x1_min",
        "x1_max",
        "x2_min",
        "x2_max",
        "x3_min",
        "x3_max",
        "vol_func",
        "center_func_1",
        "center_func_2",
        "center_func_3",
    ]
    assert list(signature.parameters)[:-1] == positional_names
    selector = signature.parameters["meshblock_index_in_file"]
    assert selector.kind is inspect.Parameter.KEYWORD_ONLY
    assert selector.default == 0


def test_rank_binary_as_athdf_is_a_thin_canonical_delegate(monkeypatch):
    sentinel = object()
    calls = []

    def fake_read_binary_as_athdf(**kwargs):
        calls.append(kwargs)
        return sentinel

    monkeypatch.setattr(bin_convert, "read_binary_as_athdf", fake_read_binary_as_athdf)

    assert bin_convert.read_rank_binary_as_athdf("rank.bin", raw=True) is sentinel
    assert calls[0]["filename"] == "rank.bin"
    assert calls[0]["raw"] is True
    assert calls[0]["num_ghost"] == 0


def test_rank_binary_as_athdf_places_rank_one_at_its_logical_location():
    path = _binary_fixture("bin", "per_rank", "00000", rank=1)
    raw = read_binary(str(path))
    positioned = bin_convert.read_rank_binary_as_athdf(str(path))

    block_width = raw["nx1_out_mb"]
    logical_x1 = raw["mb_logical"][0, 0]
    start = logical_x1 * block_width
    np.testing.assert_allclose(positioned["dens"][..., :start], 0.0)
    np.testing.assert_allclose(
        positioned["dens"][..., start:start + block_width],
        raw["mb_data"]["dens"][0],
    )


def test_athdf_orchestration_helpers_match_frozen_shared_fixtures():
    shared_binary = bin_convert.read_binary_as_athdf(
        str(_binary_fixture("bin", "shared", "00000"))
    )
    ranked_binary = bin_convert.read_all_ranks_binary_as_athdf(
        str(_binary_fixture("bin", "per_rank", "00000", rank=0))
    )
    shared_coarsened = bin_convert.read_coarsened_binary_as_athdf(
        str(_binary_fixture("cbin", "shared", "00000"))
    )
    ranked_coarsened = bin_convert.read_all_ranks_coarsened_binary_as_athdf(
        str(_binary_fixture("cbin", "per_rank", "00000", rank=0))
    )

    for actual, expected in (
        (ranked_binary, shared_binary),
        (ranked_coarsened, shared_coarsened),
    ):
        assert actual.keys() == expected.keys()
        assert actual["dens"].size > 0
        for key in ("Time", "NumCycles", "MaxLevel"):
            assert actual[key] == expected[key]
        for key in actual.keys() - {"Time", "NumCycles", "MaxLevel"}:
            np.testing.assert_allclose(actual[key], expected[key])


def test_single_rank_binary_as_athdf_keyword_index_and_legacy_positionals(tmp_path):
    path = tmp_path / "two_meshblocks.bin"
    _write_two_meshblock_binary(path)
    raw = read_binary(str(path))

    legacy = bin_convert.read_single_rank_binary_as_athdf(
        str(path), False, None, ["dens"]
    )
    first = bin_convert.read_single_rank_binary_as_athdf(
        str(path), quantities=["dens"], meshblock_index_in_file=0
    )
    second = bin_convert.read_single_rank_binary_as_athdf(
        str(path), quantities=["dens"], meshblock_index_in_file=np.int64(1)
    )

    np.testing.assert_allclose(legacy["dens"], raw["mb_data"]["dens"][0])
    np.testing.assert_allclose(first["dens"], legacy["dens"])
    np.testing.assert_allclose(second["dens"], raw["mb_data"]["dens"][1])
    assert second["x1f"][0] == raw["mb_geometry"][1, 0]
    assert second["x1f"][-1] == raw["mb_geometry"][1, 1]


def test_single_rank_binary_as_athdf_preserves_destination_and_dtype():
    path = _binary_fixture("bin", "shared", "00000")
    raw = read_binary(str(path))
    destination = np.empty_like(raw["mb_data"]["dens"][0], dtype=np.float64)
    supplied = {"dens": destination}

    result = bin_convert.read_single_rank_binary_as_athdf(
        str(path), data=supplied, quantities=["dens"], dtype=np.float64
    )

    assert result["dens"] is destination
    assert result["dens"].dtype == np.float64
    np.testing.assert_allclose(result["dens"], raw["mb_data"]["dens"][0])


@pytest.mark.parametrize("index", (True, np.bool_(True), 1.5, "1", None))
def test_single_rank_binary_as_athdf_rejects_non_integral_indexes(tmp_path, index):
    path = tmp_path / "two_meshblocks.bin"
    _write_two_meshblock_binary(path)

    with pytest.raises(TypeError, match="meshblock_index_in_file must be an integer"):
        bin_convert.read_single_rank_binary_as_athdf(
            str(path), meshblock_index_in_file=index
        )


@pytest.mark.parametrize("index", (-1, 2))
def test_single_rank_binary_as_athdf_rejects_out_of_range_indexes(tmp_path, index):
    path = tmp_path / "two_meshblocks.bin"
    _write_two_meshblock_binary(path)

    with pytest.raises(IndexError, match="meshblock_index_in_file .* is out of range"):
        bin_convert.read_single_rank_binary_as_athdf(
            str(path), meshblock_index_in_file=index
        )


def test_single_rank_binary_as_athdf_rejects_nonzero_raw_index(tmp_path):
    path = tmp_path / "two_meshblocks.bin"
    _write_two_meshblock_binary(path)

    with pytest.raises(ValueError, match="meshblock_index_in_file must be 0"):
        bin_convert.read_single_rank_binary_as_athdf(
            str(path), raw=True, meshblock_index_in_file=1
        )


def test_single_rank_binary_as_athdf_rejects_unsupported_emitted_extents(monkeypatch):
    filedata = {
        "n_mbs": 1,
        "nx1_mb": 1,
        "nx2_mb": 1,
        "nx3_mb": 1,
        "nx1_out_mb": 5,
        "nx2_out_mb": 1,
        "nx3_out_mb": 1,
    }
    monkeypatch.setattr(bin_convert, "read_binary", lambda _filename: filedata)

    with pytest.raises(ValueError, match="does not support ghost-bearing or sliced"):
        bin_convert.read_single_rank_binary_as_athdf("malformed.bin")


@pytest.mark.parametrize(
    ("kind", "reader", "family"),
    (
        ("bin", read_binary, "binary"),
        ("cbin", read_coarsened_binary, "coarsened binary"),
    ),
)
def test_binary_readers_close_malformed_record_files(
    tmp_path, monkeypatch, kind, reader, family
):
    source = _binary_fixture(kind, "shared", "00000")
    malformed = tmp_path / source.name
    malformed.write_bytes(source.read_bytes()[:-1])

    actual_open = builtins.open
    opened = []

    def tracked_open(*args, **kwargs):
        fp = actual_open(*args, **kwargs)
        opened.append(fp)
        return fp

    monkeypatch.setattr(builtins, "open", tracked_open)

    with pytest.raises(ValueError, match=f"truncated {family} meshblock values"):
        reader(str(malformed))
    assert opened
    assert all(fp.closed for fp in opened)


@pytest.mark.parametrize(
    ("kind", "reader", "family"),
    (
        ("bin", read_binary, "binary"),
        ("cbin", read_coarsened_binary, "coarsened binary"),
    ),
)
def test_binary_shard_assembly_rejects_metadata_mismatch(tmp_path, kind, reader, family):
    paths = _copy_binary_shards(
        tmp_path,
        kind,
        [
            _binary_fixture(kind, "per_rank", "00000", rank=0),
            _binary_fixture(kind, "per_rank", "00001", rank=1),
        ],
    )

    with pytest.raises(ValueError, match=f"{family} shard metadata mismatch for 'time'"):
        reader(str(paths[0]), assemble_shards=True)


@pytest.mark.parametrize(
    ("kind", "reader", "family"),
    (
        ("bin", read_binary, "binary"),
        ("cbin", read_coarsened_binary, "coarsened binary"),
    ),
)
def test_binary_shard_assembly_rejects_output_shape_mismatch(
    tmp_path, kind, reader, family
):
    paths = _copy_binary_shards(
        tmp_path,
        kind,
        [
            _binary_fixture(kind, "per_rank", "00000", rank=0),
            _binary_fixture(kind, "per_rank", "00000", rank=1),
        ],
    )
    _shrink_first_meshblock_x1(paths[1])

    with pytest.raises(ValueError, match=f"{family} shard output-shape mismatch"):
        reader(str(paths[0]), assemble_shards=True)


@pytest.mark.parametrize(
    ("kind", "reader", "family"),
    (
        ("bin", read_binary, "binary"),
        ("cbin", read_coarsened_binary, "coarsened binary"),
    ),
)
def test_binary_readers_reject_nonuniform_meshblock_extents_in_one_file(
    tmp_path, kind, reader, family
):
    malformed = tmp_path / f"mixed_extents.{kind}"
    _write_two_meshblock_binary(malformed, kind)
    _shrink_second_meshblock_x1(malformed)

    with pytest.raises(
        ValueError, match=f"{family} file .* nonuniform MeshBlock extents"
    ):
        reader(str(malformed))


def test_binary_shard_assembly_rejects_duplicate_logical_meshblocks(tmp_path):
    rank0 = _binary_fixture("bin", "per_rank", "00000", rank=0)
    paths = _copy_binary_shards(tmp_path, "bin", [rank0, rank0])

    with pytest.raises(ValueError, match="duplicate logical MeshBlock"):
        read_binary(str(paths[0]), assemble_shards=True)


@pytest.mark.parametrize(
    ("kind", "reader", "family"),
    (
        ("bin", read_binary, "binary"),
        ("cbin", read_coarsened_binary, "coarsened binary"),
    ),
)
def test_binary_readers_reject_oversized_meshblock_payload(
    tmp_path, kind, reader, family
):
    source = _binary_fixture(kind, "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    _oversize_first_meshblock_x1(malformed)

    with pytest.raises(ValueError, match=f"{family} meshblock payload .* requires"):
        reader(str(malformed))


@pytest.mark.parametrize(
    ("kind", "reader", "family"),
    (
        ("bin", read_binary, "binary"),
        ("cbin", read_coarsened_binary, "coarsened binary"),
    ),
)
def test_binary_readers_reject_oversized_parameter_header(
    tmp_path, kind, reader, family
):
    source = _binary_fixture(kind, "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    _oversize_parameter_header(malformed)

    with pytest.raises(ValueError, match=f"{family} parameter header .* declares"):
        reader(str(malformed))


@pytest.mark.parametrize(
    ("kind", "reader", "family"),
    (
        ("bin", read_binary, "binary"),
        ("cbin", read_coarsened_binary, "coarsened binary"),
    ),
)
def test_binary_readers_reject_cumulative_metadata_above_practical_limit(
    tmp_path, monkeypatch, kind, reader, family
):
    source = _binary_fixture(kind, "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    monkeypatch.setattr(bin_convert, "_MAX_BINARY_HEADER_BYTES", 80)

    with pytest.raises(ValueError, match=f"{family} metadata records .* require"):
        reader(str(malformed))


@pytest.mark.parametrize(
    ("kind", "reader", "family"),
    (
        ("bin", read_binary, "binary"),
        ("cbin", read_coarsened_binary, "coarsened binary"),
    ),
)
def test_binary_readers_reject_nonpositive_variable_count(
    tmp_path, kind, reader, family
):
    source = _binary_fixture(kind, "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    _replace_binary_line(
        malformed, b"  number of variables=", b"  number of variables=0"
    )

    with pytest.raises(ValueError, match=f"{family} file .* invalid variable count"):
        reader(str(malformed))


@pytest.mark.parametrize(
    ("kind", "reader", "family"),
    (
        ("bin", read_binary, "binary"),
        ("cbin", read_coarsened_binary, "coarsened binary"),
    ),
)
def test_binary_readers_reject_oversized_meshblock_metadata(
    tmp_path, monkeypatch, kind, reader, family
):
    source = _binary_fixture(kind, "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    monkeypatch.setattr(bin_convert, "_MAX_MESHBLOCK_ALLOCATION_BYTES", 100)

    with pytest.raises(ValueError, match=f"{family} meshblock metadata .* requires"):
        reader(str(malformed))


def test_coarsened_binary_rejects_nonpositive_coarsening_factor(tmp_path):
    source = _binary_fixture("cbin", "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    _replace_binary_line(
        malformed, b"  coarsening factor=", b"  coarsening factor=0"
    )

    with pytest.raises(ValueError, match="invalid coarsening factor"):
        read_coarsened_binary(str(malformed))


@pytest.mark.parametrize("factor", (1, 3, 100))
def test_coarsened_binary_rejects_invalid_positive_coarsening_factor(tmp_path, factor):
    source = _binary_fixture("cbin", "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    _replace_binary_line(
        malformed,
        b"  coarsening factor=",
        f"  coarsening factor={factor}".encode("ascii"),
    )

    with pytest.raises(ValueError, match="invalid coarsening factor"):
        read_coarsened_binary(str(malformed))


@pytest.mark.parametrize(
    ("mutation", "expected"),
    (
        ({"level": -1}, "invalid logical level"),
        ({"level": 31}, "invalid logical level"),
        ({"i": -1}, "negative logical location"),
        ({"i": 10**9}, "out-of-range logical location"),
    ),
)
def test_binary_rejects_malformed_logical_locations(tmp_path, mutation, expected):
    source = _binary_fixture("bin", "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    _replace_first_meshblock_logical(malformed, **mutation)

    with pytest.raises(ValueError, match=expected):
        read_binary(str(malformed))


def test_binary_rejects_duplicate_variable_names(tmp_path):
    source = _binary_fixture("bin", "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    _replace_binary_line(
        malformed,
        b"  variables:",
        b"  variables:  velx  velx  vely  velz  eint  ",
    )

    with pytest.raises(ValueError, match="duplicate variable names"):
        read_binary(str(malformed))


def test_binary_rejects_nonfinite_time(tmp_path):
    source = _binary_fixture("bin", "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    _replace_binary_line(malformed, b"  time=", b"  time=nan")

    with pytest.raises(ValueError, match="non-finite metadata"):
        read_binary(str(malformed))


@pytest.mark.parametrize(
    ("before", "after", "expected"),
    (
        (b"nx1    = 16       ", b"nx1    = 0", "invalid root-grid dimensions"),
        (b"nghost = 2        ", b"nghost = -1", "negative ghost-zone count"),
        (b"x1max  = 0.5      ", b"x1max  = -0.5", "invalid coordinate bounds"),
    ),
)
def test_binary_rejects_malformed_grid_metadata(tmp_path, before, after, expected):
    source = _binary_fixture("bin", "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    _replace_binary_bytes_same_width(malformed, before, after)

    with pytest.raises(ValueError, match=expected):
        read_binary(str(malformed))


def test_binary_rejects_invalid_meshblock_geometry(tmp_path):
    source = _binary_fixture("bin", "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    first = read_binary(str(malformed))["mb_geometry"][0]
    _replace_first_meshblock_geometry(malformed, x1max=first[0])

    with pytest.raises(ValueError, match="invalid geometry"):
        read_binary(str(malformed))


def test_binary_rejects_meshblock_geometry_outside_root_domain(tmp_path):
    source = _binary_fixture("bin", "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    _replace_first_meshblock_geometry(malformed, x1max=10.0)

    with pytest.raises(ValueError, match="geometry outside its logical location"):
        read_binary(str(malformed))


def test_binary_rejects_shrunken_meshblock_geometry(tmp_path):
    source = _binary_fixture("bin", "shared", "00000")
    malformed = tmp_path / source.name
    shutil.copyfile(source, malformed)
    first = read_binary(str(malformed))["mb_geometry"][0]
    _replace_first_meshblock_geometry(malformed, x1min=first[0] + 0.01)

    with pytest.raises(ValueError, match="geometry outside its logical location"):
        read_binary(str(malformed))


def test_binary_rejects_large_offset_shifted_meshblock_geometry():
    with pytest.raises(ValueError, match="geometry outside its logical location"):
        bin_convert._validate_meshblock_metadata(
            (np.array((0, 0, 0, 0), dtype=np.int32),),
            (np.array((1000001.0, 1000005.0, 0.0, 1.0, 0.0, 1.0)),),
            "binary",
            "large-offset.bin",
            (1, 1, 1),
            (4, 1, 1),
            (1000000.0, 1000004.0, 0.0, 1.0, 0.0, 1.0),
            0,
        )


def test_binary_rejects_meshblock_geometry_from_another_logical_location(tmp_path):
    malformed = tmp_path / "wrong_logical_geometry.bin"
    _write_two_meshblock_binary(malformed)
    second = read_binary(str(malformed))["mb_geometry"][1]
    _replace_first_meshblock_geometry(malformed, x1max=second[1])

    with pytest.raises(ValueError, match="geometry outside its logical location"):
        read_binary(str(malformed))


def test_binary_rejects_duplicate_logical_meshblocks_in_direct_file(tmp_path):
    malformed = tmp_path / "duplicate.bin"
    _write_duplicate_meshblock_binary(malformed)

    with pytest.raises(ValueError, match="duplicate logical MeshBlocks"):
        read_binary(str(malformed))


@pytest.mark.parametrize("kwargs", ({}, {"fast_restrict": True}, {"subsample": True}))
def test_binary_as_athdf_restricts_fine_meshblock_and_sets_levels(tmp_path, kwargs):
    source = _binary_fixture("bin", "shared", "00000")
    refined = tmp_path / source.name
    shutil.copyfile(source, refined)
    _replace_first_meshblock_logical(refined, level=1)
    _replace_first_meshblock_geometry(refined, x1max=0.0, x2max=0.0, x3max=0.0)

    result = bin_convert.read_binary_as_athdf(
        str(refined), level=0, return_levels=True, **kwargs
    )

    assert result["dens"].shape == (4, 4, 16)
    assert np.count_nonzero(result["dens"]) > 0
    assert np.count_nonzero(result["Levels"] == 1) > 0


def test_binary_as_athdf_crops_selected_prolongation():
    path = _binary_fixture("bin", "shared", "00000")

    result = bin_convert.read_binary_as_athdf(
        str(path), level=1, x1_min=-0.5, x1_max=-0.4375
    )

    assert result["dens"].shape == (8, 8, 2)
    assert result["x1f"].shape == (3,)
    assert result["x1v"].shape == (2,)


def test_binary_as_athdf_lower_crop_retains_intersecting_cell():
    path = _binary_fixture("bin", "shared", "00000")

    result = bin_convert.read_binary_as_athdf(str(path), level=1, x1_min=-0.49)

    assert result["x1f"][0] == pytest.approx(-0.5)


def test_rank_binary_as_athdf_marks_uncovered_levels_deterministically():
    path = _binary_fixture("bin", "per_rank", "00000", rank=0)

    result = bin_convert.read_rank_binary_as_athdf(str(path), return_levels=True)

    assert np.count_nonzero(result["Levels"] == -1) > 0
    assert set(np.unique(result["Levels"])) == {-1, 0}


def test_athdf_like_ghost_placement_uses_interior_width_and_extended_coordinates():
    data = {
        "dens": np.zeros((1, 1, 6), dtype=np.float32),
        "Levels": np.full((1, 1, 6), -1, dtype=np.int32),
    }
    filedata = {
        "mb_logical": np.array(((0, 0, 0, 0), (1, 0, 0, 0)), dtype=np.int32),
        "mb_data": {
            "dens": [
                np.array([[[10.0, 11.0, 12.0, 13.0]]], dtype=np.float32),
                np.array([[[20.0, 21.0, 22.0, 23.0]]], dtype=np.float32),
            ]
        },
    }
    for block_num in range(2):
        bin_convert._copy_meshblock_to_athdf(
            data,
            filedata,
            ("dens",),
            block_num,
            0,
            (4, 1, 1),
            (6, 1, 1),
            (0, 6, 0, 1, 0, 1),
            128,
            True,
            False,
            False,
            None,
            None,
            1,
        )
    np.testing.assert_array_equal(
        data["dens"][0, 0], np.array([10.0, 11.0, 20.0, 21.0, 22.0, 23.0])
    )

    coordinates = {}
    bin_convert._populate_root_athdf_coordinates(
        coordinates,
        {"x1min": 0.0, "x1max": 4.0, "x2min": 0.0, "x2max": 1.0,
         "x3min": 0.0, "x3max": 1.0},
        (6, 1, 1),
        (4, 1, 1),
        0,
        1,
        np.float32,
        (lambda left, right: 0.5 * (left + right),) * 3,
    )
    assert coordinates["x1f"][0] == pytest.approx(-1.0)
    assert coordinates["x1f"][-1] == pytest.approx(5.0)


def test_athdf_like_rejects_detected_ghost_zones_when_count_is_omitted():
    filedata = {
        "nx1_mb": 2,
        "nx2_mb": 1,
        "nx3_mb": 1,
        "mb_logical": np.array(((0, 0, 0, 0),), dtype=np.int32),
    }

    with pytest.raises(ValueError, match="detected ghost zones"):
        bin_convert._validate_athdf_num_ghost(
            filedata, (4, 1, 1), (4, 1, 1), 0, 0, 0, (None,) * 6
        )


@pytest.mark.parametrize(
    "num_ghost, expected", ((0, "detected ghost zones"), (1, "does not match"))
)
def test_athdf_like_rejects_singleton_axis_ghost_width_mismatch(num_ghost, expected):
    filedata = {
        "nx1_mb": 1,
        "nx2_mb": 1,
        "nx3_mb": 1,
        "mb_logical": np.array(((0, 0, 0, 0),), dtype=np.int32),
    }

    with pytest.raises(ValueError, match=expected):
        bin_convert._validate_athdf_num_ghost(
            filedata, (5, 1, 1), (1, 1, 1), num_ghost, 0, 0, (None,) * 6
        )


def test_athdf_like_preflight_rejects_oversized_coordinates(monkeypatch):
    monkeypatch.setattr(bin_convert, "_MAX_ATHDF_ALLOCATION_BYTES", 16)

    with pytest.raises(ValueError, match="athdf-like coordinate arrays requires"):
        bin_convert._preflight_athdf_coordinates((3, 1, 1), np.float32)


def test_athdf_like_preflight_counts_coordinate_generation_temporary(monkeypatch):
    monkeypatch.setattr(bin_convert, "_MAX_ATHDF_ALLOCATION_BYTES", 60)

    with pytest.raises(
        ValueError, match="athdf-like coordinate generation peak requires"
    ):
        bin_convert._preflight_athdf_coordinates((3, 1, 1), np.float32)


def test_athdf_like_preflight_rejects_oversized_outputs(monkeypatch):
    monkeypatch.setattr(bin_convert, "_MAX_ATHDF_ALLOCATION_BYTES", 16)

    with pytest.raises(ValueError, match="athdf-like output arrays requires"):
        bin_convert._preflight_athdf_outputs(
            (3, 1, 1), ("dens", "eint"), np.float32, False, None
        )


def test_athdf_like_preflight_rejects_cumulative_coordinates_and_outputs(monkeypatch):
    monkeypatch.setattr(bin_convert, "_MAX_ATHDF_ALLOCATION_BYTES", 100)
    coordinate_bytes = bin_convert._preflight_athdf_coordinates((3, 1, 1), np.float32)

    with pytest.raises(ValueError, match="athdf-like output arrays requires"):
        bin_convert._preflight_athdf_outputs(
            (10, 1, 1), ("dens", "eint"), np.float32, False, None, coordinate_bytes
        )


def test_athdf_like_exact_restriction_preflights_temporary_arrays(monkeypatch):
    monkeypatch.setattr(bin_convert, "_MAX_ATHDF_ALLOCATION_BYTES", 100)
    data = {"dens": np.zeros((1, 1, 2), dtype=np.float32)}
    filedata = {
        "mb_logical": np.array(((0, 0, 0, 1),), dtype=np.int32),
        "mb_geometry": np.array(((0.0, 1.0, 0.0, 1.0, 0.0, 1.0),)),
        "mb_data": {"dens": [np.ones((1, 1, 2), dtype=np.float32)]},
    }

    with pytest.raises(ValueError, match="athdf-like exact restriction peak requires"):
        bin_convert._copy_meshblock_to_athdf(
            data,
            filedata,
            ("dens",),
            0,
            0,
            (2, 1, 1),
            (2, 1, 1),
            (0, 2, 0, 1, 0, 1),
            90,
            False,
            False,
            False,
            np.zeros((1, 1, 1), dtype=bool),
            lambda xm, xp, ym, yp, zm, zp: (xp - xm) * (yp - ym) * (zp - zm),
            0,
        )


def test_athdf_like_prolongation_preflights_index_arrays(monkeypatch):
    monkeypatch.setattr(bin_convert, "_MAX_ATHDF_ALLOCATION_BYTES", 40)
    data = {"dens": np.zeros((1, 1, 2), dtype=np.float32)}
    filedata = {
        "mb_logical": np.array(((0, 0, 0, 0),), dtype=np.int32),
        "mb_data": {"dens": [np.ones((1, 1, 1), dtype=np.float32)]},
    }

    with pytest.raises(ValueError, match="athdf-like prolongation index peak requires"):
        bin_convert._copy_meshblock_to_athdf(
            data,
            filedata,
            ("dens",),
            0,
            1,
            (1, 1, 1),
            (2, 1, 1),
            (0, 2, 0, 1, 0, 1),
            32,
            False,
            False,
            False,
            None,
            None,
            0,
        )


def test_athdf_like_prolongation_preflights_index_source_arrays(monkeypatch):
    monkeypatch.setattr(bin_convert, "_MAX_ATHDF_ALLOCATION_BYTES", 75)
    data = {"dens": np.zeros((1, 1, 2), dtype=np.float32)}
    filedata = {
        "mb_logical": np.array(((0, 0, 0, 0),), dtype=np.int32),
        "mb_data": {"dens": [np.ones((1, 1, 1), dtype=np.float32)]},
    }

    with pytest.raises(ValueError, match="athdf-like prolongation index peak requires"):
        bin_convert._copy_meshblock_to_athdf(
            data,
            filedata,
            ("dens",),
            0,
            1,
            (1, 1, 1),
            (2, 1, 1),
            (0, 2, 0, 1, 0, 1),
            32,
            False,
            False,
            False,
            None,
            None,
            0,
        )


def test_athdf_like_exact_restriction_preflights_repeat_source_arrays(monkeypatch):
    monkeypatch.setattr(bin_convert, "_MAX_ATHDF_ALLOCATION_BYTES", 192)
    data = {"dens": np.zeros((1, 1, 2), dtype=np.float32)}
    filedata = {
        "mb_logical": np.array(((0, 0, 0, 1),), dtype=np.int32),
        "mb_geometry": np.array(((0.0, 1.0, 0.0, 1.0, 0.0, 1.0),)),
        "mb_data": {"dens": [np.ones((1, 1, 2), dtype=np.float32)]},
    }

    with pytest.raises(ValueError, match="athdf-like exact restriction peak requires"):
        bin_convert._copy_meshblock_to_athdf(
            data,
            filedata,
            ("dens",),
            0,
            0,
            (2, 1, 1),
            (2, 1, 1),
            (0, 2, 0, 1, 0, 1),
            100,
            False,
            False,
            False,
            np.zeros((1, 1, 1), dtype=bool),
            lambda xm, xp, ym, yp, zm, zp: (xp - xm) * (yp - ym) * (zp - zm),
            0,
        )


@pytest.mark.parametrize(
    ("kind", "reader", "family"),
    (
        ("bin", read_binary, "binary"),
        ("cbin", read_coarsened_binary, "coarsened binary"),
    ),
)
def test_binary_shard_assembly_rejects_oversized_aggregate(
    tmp_path, monkeypatch, kind, reader, family
):
    paths = _copy_binary_shards(
        tmp_path,
        kind,
        [
            _binary_fixture(kind, "per_rank", "00000", rank=0),
            _binary_fixture(kind, "per_rank", "00000", rank=1),
        ],
    )
    first = reader(str(paths[0]))
    single_shard_bytes = sum(
        values.nbytes
        for variable in first["var_names"]
        for values in first["mb_data"][variable]
    )
    single_shard_metadata_bytes = sum(
        first[key].nbytes for key in ("mb_index", "mb_logical", "mb_geometry")
    )
    monkeypatch.setattr(
        bin_convert,
        "_MAX_MESHBLOCK_ALLOCATION_BYTES",
        max(single_shard_bytes, 2 * single_shard_metadata_bytes),
    )

    with pytest.raises(
        ValueError, match=f"{family} assembled shard (payload|metadata) requires"
    ):
        reader(str(paths[0]), assemble_shards=True)


@pytest.mark.parametrize(
    ("kind", "reader", "athdf_readers"),
    (
        (
            "bin",
            read_binary,
            (
                bin_convert.read_binary_as_athdf,
                bin_convert.read_rank_binary_as_athdf,
                bin_convert.read_single_rank_binary_as_athdf,
                bin_convert.read_all_ranks_binary_as_athdf,
            ),
        ),
        (
            "cbin",
            read_coarsened_binary,
            (
                bin_convert.read_coarsened_binary_as_athdf,
                bin_convert.read_all_ranks_coarsened_binary_as_athdf,
            ),
        ),
    ),
)
def test_athdf_like_readers_reject_valid_empty_shards(
    tmp_path, kind, reader, athdf_readers
):
    source = _binary_fixture(kind, "per_rank", "00000", rank=0)
    empty = tmp_path / "rank_00000000" / f"empty.00000.{kind}"
    _write_empty_binary(empty, source)

    assert reader(str(empty))["n_mbs"] == 0
    for athdf_reader in athdf_readers:
        with pytest.raises(
            ValueError, match="athdf-like data: binary output contains no"
        ):
            athdf_reader(str(empty))


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


def _write_pdf_header(
    path, fmt, distribution=None, v2=False, *, shard_id=None, sibling_count=None
):
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
    if shard_id is not None:
        lines.insert(1, f"{distribution} = {shard_id}")
    if sibling_count is not None:
        lines.insert(1, f"number_of_{distribution}s = {sibling_count}")
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


def _write_v2_sparse_pdf(
    path, entries, time=0.25, cycle=5, declared_count=None, rank=0
):
    path.parent.mkdir(parents=True, exist_ok=True)
    count = len(entries) if declared_count is None else declared_count
    payload = struct.pack(
        "=8sIIIIQdq", b"AKPDFV2\0", 2, 1, 1, rank, count, time, cycle
    )
    payload += b"".join(struct.pack("=Qd", index, value) for index, value in entries)
    path.write_bytes(payload)


def _write_v2_dense_pdf(path, values, time=0.25, cycle=5):
    values = np.asarray(values, dtype=np.float64)
    payload = struct.pack(
        "=8sIIIIQdq", b"AKPDFV2\0", 2, 0, 1, 0, values.size, time, cycle
    )
    path.write_bytes(payload + values.tobytes())


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


def test_modern_pdf_node_inventory_accepts_explicit_empty_shard(tmp_path):
    case = tmp_path / "node_inventory"
    node0_dir = case / "node_00000000"
    node1_dir = case / "node_00000001"
    for node_id, directory in enumerate((node0_dir, node1_dir)):
        _write_pdf_header(
            directory / "hist.header.pdf",
            "sparse_coo",
            "node",
            v2=True,
            shard_id=node_id,
            sibling_count=2,
        )
    node0 = node0_dir / "hist.00000.pdf"
    node1 = node1_dir / "hist.00000.pdf"
    _write_v2_sparse_pdf(node0, [(0, 1.0), (3, 2.5)])
    _write_v2_sparse_pdf(node1, [])

    result = read_pdf(str(node0))
    np.testing.assert_allclose(result["pdf"], [1.0, 0.0, 0.0, 2.5])
    assert "node" not in result["header"]
    assert result["header"]["number_of_nodes"] == 2


def test_modern_pdf_v2_dense_header_requires_v2_payload(tmp_path):
    data = tmp_path / "dense.00000.pdf"
    _write_pdf_header(tmp_path / "dense.header.pdf", "dense", v2=True)
    _write_dense_pdf(data, [1.0, 0.0, 0.0, 2.5])

    with pytest.raises(ValueError, match="missing its AKPDFV2 preamble"):
        read_pdf(str(data))


def test_modern_pdf_v2_sparse_header_requires_v2_payload(tmp_path):
    case = tmp_path / "sparse_missing_preamble"
    shard_dir = case / "rank_00000000"
    _write_pdf_header(
        shard_dir / "hist.header.pdf",
        "sparse_coo",
        "rank",
        v2=True,
        shard_id=0,
        sibling_count=1,
    )
    shard = shard_dir / "hist.00000.pdf"
    _write_sparse_pdf(shard, [(0, 1.0)])

    with pytest.raises(ValueError, match="missing its AKPDFV2 preamble"):
        read_pdf(str(shard))


def test_modern_pdf_v2_dense_cycle_must_match_header(tmp_path):
    data = tmp_path / "dense.00000.pdf"
    _write_pdf_header(tmp_path / "dense.header.pdf", "dense", v2=True)
    _write_v2_dense_pdf(data, [1.0, 0.0, 0.0, 2.5], cycle=99)

    with pytest.raises(ValueError, match="PDF V2 cycle mismatch"):
        read_pdf(str(data))


def test_modern_pdf_v2_sparse_cycle_must_match_header(tmp_path):
    case = tmp_path / "sparse_cycle"
    shard_dir = case / "rank_00000000"
    _write_pdf_header(
        shard_dir / "hist.header.pdf",
        "sparse_coo",
        "rank",
        v2=True,
        shard_id=0,
        sibling_count=1,
    )
    shard = shard_dir / "hist.00000.pdf"
    _write_v2_sparse_pdf(shard, [(0, 1.0)], cycle=99)

    with pytest.raises(ValueError, match="PDF V2 cycle mismatch"):
        read_pdf(str(shard))


def test_modern_pdf_sparse_sibling_headers_must_agree_on_v2_magic(tmp_path):
    case = tmp_path / "sparse_magic_mismatch"
    rank0_dir = case / "rank_00000000"
    rank1_dir = case / "rank_00000001"
    _write_pdf_header(
        rank0_dir / "hist.header.pdf",
        "sparse_coo",
        "rank",
        shard_id=0,
        sibling_count=2,
    )
    _write_pdf_header(
        rank1_dir / "hist.header.pdf",
        "sparse_coo",
        "rank",
        v2=True,
        shard_id=1,
        sibling_count=2,
    )
    rank0 = rank0_dir / "hist.00000.pdf"
    rank1 = rank1_dir / "hist.00000.pdf"
    _write_sparse_pdf(rank0, [(0, 1.0)])
    _write_sparse_pdf(rank1, [(3, 2.5)])

    with pytest.raises(ValueError, match="metadata mismatch for 'binary_magic'"):
        read_pdf(str(rank0))


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


def test_modern_pdf_rejects_alias_sibling_directory(tmp_path):
    case = tmp_path / "alias"
    _write_pdf_header(case / "hist.header.pdf", "sparse_coo", "rank")
    canonical = case / "rank_00000000" / "hist.00000.pdf"
    alias = case / "rank_0" / "hist.00000.pdf"
    _write_sparse_pdf(canonical, [(0, 1.0)])
    _write_sparse_pdf(alias, [(1, 1.0)])

    with pytest.raises(ValueError, match="invalid PDF shard directory"):
        read_pdf(str(canonical))


def test_modern_pdf_rejects_rank_payload_path_id_mismatch(tmp_path):
    case = tmp_path / "rank_mismatch"
    _write_pdf_header(
        case / "hist.header.pdf",
        "sparse_coo",
        "rank",
        v2=True,
        shard_id=0,
        sibling_count=1,
    )
    shard = case / "rank_00000000" / "hist.00000.pdf"
    _write_v2_sparse_pdf(shard, [(0, 1.0)], rank=1)

    with pytest.raises(ValueError, match="payload declares rank=1"):
        read_pdf(str(shard))


def test_modern_pdf_rejects_missing_declared_node_shard(tmp_path):
    case = tmp_path / "missing_node"
    shard_dir = case / "node_00000000"
    _write_pdf_header(
        shard_dir / "hist.header.pdf",
        "sparse_coo",
        "node",
        v2=True,
        shard_id=0,
        sibling_count=2,
    )
    shard = shard_dir / "hist.00000.pdf"
    _write_v2_sparse_pdf(shard, [(0, 1.0)])

    with pytest.raises(ValueError, match="node shard inventory is incomplete"):
        read_pdf(str(shard))


def test_modern_pdf_rejects_node_header_path_id_mismatch(tmp_path):
    case = tmp_path / "node_mismatch"
    shard_dir = case / "node_00000000"
    _write_pdf_header(
        shard_dir / "hist.header.pdf",
        "sparse_coo",
        "node",
        v2=True,
        shard_id=1,
        sibling_count=1,
    )
    shard = shard_dir / "hist.00000.pdf"
    _write_v2_sparse_pdf(shard, [(0, 1.0)])

    with pytest.raises(ValueError, match="declares node=1"):
        read_pdf(str(shard))


def test_modern_pdf_rejects_missing_observed_rank_id(tmp_path):
    case = tmp_path / "missing_rank_id"
    _write_pdf_header(case / "hist.header.pdf", "sparse_coo", "rank")
    rank0 = case / "rank_00000000" / "hist.00000.pdf"
    rank2 = case / "rank_00000002" / "hist.00000.pdf"
    _write_sparse_pdf(rank0, [(0, 1.0)])
    _write_sparse_pdf(rank2, [(1, 1.0)])

    with pytest.raises(ValueError, match="rank shard inventory is incomplete"):
        read_pdf(str(rank0))


def test_modern_pdf_rejects_oversized_declared_node_count(tmp_path):
    case = tmp_path / "oversized_node_count"
    shard_dir = case / "node_00000000"
    _write_pdf_header(
        shard_dir / "hist.header.pdf",
        "sparse_coo",
        "node",
        v2=True,
        shard_id=0,
        sibling_count=10**12,
    )
    shard = shard_dir / "hist.00000.pdf"
    _write_v2_sparse_pdf(shard, [(0, 1.0)])

    with pytest.raises(ValueError, match="expected 1000000000000 shards, found 1"):
        read_pdf(str(shard))


def test_modern_pdf_rejects_generated_edges_above_practical_limit(
    tmp_path, monkeypatch
):
    header = tmp_path / "edges.header.pdf"
    _write_pdf_header(header, "dense")
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 16)

    with pytest.raises(ValueError, match="PDF dimension 1 bin edges requires"):
        read_pdf_header(str(header))


def test_modern_pdf_preflights_generated_edges_before_numpy_generation(
    tmp_path, monkeypatch
):
    header = tmp_path / "edges.header.pdf"
    _write_pdf_header(header, "dense")
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 24)

    def fail_generated_edges(*args, **kwargs):
        raise AssertionError("PDF edges were generated before retained-peak preflight")

    monkeypatch.setattr(read_pdf_module, "_generated_edges", fail_generated_edges)
    with pytest.raises(ValueError, match="PDF generated bin-edge peak requires"):
        read_pdf_module._read_pdf_header(str(header), externally_retained_bytes=1)


def test_modern_pdf_rejects_dense_histogram_above_practical_limit(
    tmp_path, monkeypatch
):
    header = tmp_path / "dense.header.pdf"
    _write_pdf_header(header, "dense")
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 70)

    with pytest.raises(
        ValueError, match="PDF retained histogram and bin metadata requires"
    ):
        read_pdf_header(str(header))


def test_modern_pdf_rejects_cumulative_bin_metadata_above_practical_limit(
    tmp_path, monkeypatch
):
    header = tmp_path / "dense.header.pdf"
    _write_pdf_header(header, "dense")
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 55)

    with pytest.raises(ValueError, match="PDF retained bin metadata requires"):
        read_pdf_header(str(header))


def test_modern_pdf_preflights_explicit_edges_before_numpy_array(
    tmp_path, monkeypatch
):
    header = tmp_path / "explicit.header.pdf"
    _write_pdf_header(header, "dense")
    header.write_text(header.read_text() + "bin_edges_1 = 0.0 1.0 2.0\n")
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 100)

    def fail_array(*args, **kwargs):
        raise AssertionError("explicit edges were materialized before preflight")

    monkeypatch.setattr(read_pdf_module.np, "array", fail_array)
    with pytest.raises(ValueError, match="PDF bin_edges_1 parsing peak requires"):
        read_pdf_header(str(header))


def test_modern_pdf_preflights_edge_validation_before_numpy_diff(
    tmp_path, monkeypatch
):
    header = tmp_path / "dense.header.pdf"
    _write_pdf_header(header, "dense")
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 44)

    def fail_diff(*args, **kwargs):
        raise AssertionError("PDF edges were differenced before validation preflight")

    monkeypatch.setattr(read_pdf_module.np, "diff", fail_diff)
    with pytest.raises(ValueError, match="PDF bin-edge validation peak requires"):
        read_pdf_header(str(header))


def test_legacy_pdf_preflights_cumulative_parsing_before_second_numpy_array(
    tmp_path, monkeypatch
):
    header = tmp_path / "legacy.bins.pdf"
    header.write_text("# [1] = dens\n# [2] = velx\n0.0 1.0 2.0\n0.0 1.0 2.0\n")
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 486)
    array = read_pdf_module.np.array
    calls = 0

    def tracked_array(*args, **kwargs):
        nonlocal calls
        calls += 1
        return array(*args, **kwargs)

    monkeypatch.setattr(read_pdf_module.np, "array", tracked_array)
    with pytest.raises(ValueError, match="legacy PDF bin edges parsing peak requires"):
        read_pdf_header(str(header))
    assert calls == 1


def test_legacy_pdf_preflights_edge_validation_before_numpy_diff(
    tmp_path, monkeypatch
):
    header = tmp_path / "legacy.bins.pdf"
    header.write_text("# [1] = dens\n0.0 1.0 2.0\n")
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 44)

    def fail_diff(*args, **kwargs):
        raise AssertionError("legacy PDF edges were differenced before preflight")

    monkeypatch.setattr(
        read_pdf_module,
        "_parse_ascii_floats",
        lambda *args, **kwargs: np.array((0.0, 1.0, 2.0)),
    )
    monkeypatch.setattr(read_pdf_module.np, "diff", fail_diff)
    with pytest.raises(
        ValueError, match="legacy PDF bin-edge validation peak requires"
    ):
        read_pdf_header(str(header))


def test_modern_pdf_counts_future_explicit_dimensions_during_validation(
    tmp_path, monkeypatch
):
    header = tmp_path / "explicit_2d.header.pdf"
    header.write_text(
        "format = dense\n"
        "ndim = 2\n"
        "variable_1 = dens\nnbin1 = 1\nbin1_min = 0.0\nbin1_max = 1.0\n"
        "scale1 = linear\nstride1 = 3\nbin_edges_1 = 0.0 1.0\n"
        "variable_2 = velx\nnbin2 = 1\nbin2_min = 0.0\nbin2_max = 1.0\n"
        "scale2 = linear\nstride2 = 1\nbin_edges_2 = 0.0 1.0\n"
        "total_bins = 9\ncycle = 5\n"
    )
    array = read_pdf_module.np.array
    monkeypatch.setattr(
        read_pdf_module,
        "_parse_ascii_floats",
        lambda text, label, externally_retained_bytes=0: array(
            [float(token) for token in text.split()], dtype=np.float64
        ),
    )
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 42)

    with pytest.raises(ValueError, match="PDF bin-edge validation peak requires"):
        read_pdf_header(str(header))


def test_modern_pdf_explicit_edges_count_retained_reference_header(
    tmp_path, monkeypatch
):
    header = tmp_path / "explicit.header.pdf"
    _write_pdf_header(header, "dense")
    edge_text = "0.0 1.0 2.0"
    header.write_text(header.read_text() + f"bin_edges_1 = {edge_text}\n")
    peak = read_pdf_module._ascii_float_parse_peak_bytes(edge_text, 3)
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", peak)

    with pytest.raises(ValueError, match="PDF bin_edges_1 parsing peak requires"):
        read_pdf_module._read_pdf_header(str(header), externally_retained_bytes=1)


def test_legacy_pdf_edges_count_retained_reference_header(tmp_path, monkeypatch):
    header = tmp_path / "legacy.bins.pdf"
    edge_text = "0.0 1.0 2.0"
    header.write_text(f"# [1] = dens\n{edge_text}\n")
    peak = read_pdf_module._ascii_float_parse_peak_bytes(edge_text, 3)
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", peak)

    with pytest.raises(ValueError, match="legacy PDF bin edges parsing peak requires"):
        read_pdf_module._read_pdf_header(str(header), externally_retained_bytes=1)


def test_legacy_pdf_rejects_header_trailing_junk(tmp_path):
    header = tmp_path / "legacy.bins.pdf"
    header.write_text("# [1] = dens\n0.0 1.0 2.0 junk\n")

    with pytest.raises(ValueError, match="invalid numeric token"):
        read_pdf_header(str(header))


def test_legacy_pdf_rejects_non_ascii_numeric_tokens(tmp_path):
    header = tmp_path / "legacy.bins.pdf"
    header.write_text("# [1] = dens\n0.0 \U0001d7d9 2.0\n")

    with pytest.raises(ValueError, match="must contain ASCII numeric tokens"):
        read_pdf_header(str(header))


def test_legacy_pdf_rejects_payload_trailing_junk(tmp_path):
    header = tmp_path / "legacy.bins.pdf"
    data = tmp_path / "legacy.00000.pdf"
    header.write_text("# [1] = dens\n0.0 1.0 2.0\n")
    data.write_text("# time = 0.0\n0.0 1.0 2.0 3.0 junk\n")

    with pytest.raises(ValueError, match="invalid numeric token"):
        read_pdf(str(data))


def test_legacy_pdf_rejects_time_trailing_junk(tmp_path):
    header = tmp_path / "legacy.bins.pdf"
    data = tmp_path / "legacy.00000.pdf"
    header.write_text("# [1] = dens\n0.0 1.0 2.0\n")
    data.write_text("# time = 0.0 junk\n0.0 1.0 2.0 3.0\n")

    with pytest.raises(ValueError, match="missing its time header"):
        read_pdf(str(data))


def test_legacy_pdf_preflights_payload_rows_before_second_numpy_array(
    tmp_path, monkeypatch
):
    header_path = tmp_path / "legacy.bins.pdf"
    data = tmp_path / "legacy.00000.pdf"
    header_path.write_text("# [1] = dens\n# [2] = velx\n0 1 2\n0 1 2\n")
    data.write_text("# time = 0.0\n0 1 2 3\n0 1 2 3\n0 1 2 3\n0 1 2 3\n")
    header = read_pdf_header(str(header_path))
    retained = read_pdf_module._header_retained_bytes(header)
    first_peak = retained + read_pdf_module._ascii_float_parse_peak_bytes("0 1 2 3", 4)
    row_bytes = 4 * np.dtype(np.float64).itemsize
    monkeypatch.setattr(
        read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", first_peak + row_bytes - 1
    )
    array = read_pdf_module.np.array
    calls = 0

    def tracked_array(*args, **kwargs):
        nonlocal calls
        calls += 1
        return array(*args, **kwargs)

    monkeypatch.setattr(read_pdf_module.np, "array", tracked_array)
    with pytest.raises(ValueError, match="legacy PDF payload row parsing peak requires"):
        read_pdf_module._read_legacy_data(str(data), header)
    assert calls == 1


def test_legacy_pdf_preflights_vstack_peak(tmp_path, monkeypatch):
    header_path = tmp_path / "legacy.bins.pdf"
    data = tmp_path / "legacy.00000.pdf"
    second_edges = " ".join(str(value) for value in range(20))
    header_path.write_text(f"# [1] = dens\n# [2] = velx\n0 1 2\n{second_edges}\n")
    rows = ["0 1 2 3"] * 21
    data.write_text("# time = 0.0\n" + "\n".join(rows) + "\n")
    header = read_pdf_header(str(header_path))
    retained = read_pdf_module._header_retained_bytes(header)
    row_bytes = 4 * np.dtype(np.float64).itemsize
    parse_peak = (
        retained
        + (len(rows) - 1) * row_bytes
        + read_pdf_module._ascii_float_parse_peak_bytes(rows[0], 4)
    )
    stack_peak = retained + len(rows) * row_bytes + header["total_bins"] * 8
    assert parse_peak < stack_peak
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", stack_peak - 1)

    def fail_vstack(*args, **kwargs):
        raise AssertionError("legacy rows were stacked before peak preflight")

    monkeypatch.setattr(read_pdf_module.np, "vstack", fail_vstack)
    with pytest.raises(ValueError, match="legacy PDF stacked payload peak requires"):
        read_pdf_module._read_legacy_data(str(data), header)


@pytest.mark.parametrize("key", ("bin1_min", "bin1_max"))
def test_modern_pdf_rejects_nonfinite_bounds(tmp_path, key):
    header = tmp_path / "dense.header.pdf"
    _write_pdf_header(header, "dense")
    text = header.read_text().replace(f"{key} = 0.0", f"{key} = nan")
    text = text.replace(f"{key} = 2.0", f"{key} = nan")
    header.write_text(text)

    with pytest.raises(ValueError, match="non-finite bin bounds"):
        read_pdf_header(str(header))


def test_modern_pdf_rejects_nonfinite_symlog_threshold(tmp_path):
    header = tmp_path / "dense.header.pdf"
    _write_pdf_header(header, "dense")
    text = header.read_text().replace("scale1 = linear", "scale1 = symlog")
    header.write_text(text + "linthresh1 = nan\n")

    with pytest.raises(ValueError, match="invalid linthresh1"):
        read_pdf_header(str(header))


def test_modern_pdf_rejects_symlog_metadata_peak_above_limit(tmp_path, monkeypatch):
    header = tmp_path / "dense.header.pdf"
    _write_pdf_header(header, "dense")
    text = header.read_text().replace("scale1 = linear", "scale1 = symlog")
    header.write_text(text + "linthresh1 = 0.1\n")
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 100)

    with pytest.raises(ValueError, match="PDF generated bin-edge peak requires"):
        read_pdf_header(str(header))


def test_modern_pdf_rejects_nonfinite_v2_payload_time(tmp_path):
    data = tmp_path / "dense.00000.pdf"
    _write_pdf_header(tmp_path / "dense.header.pdf", "dense", v2=True)
    _write_v2_dense_pdf(data, [1.0, 0.0, 0.0, 2.5], time=np.nan)

    with pytest.raises(ValueError, match="payload time .* must be finite"):
        read_pdf(str(data))


def test_modern_pdf_rejects_invalid_dense_distribution(tmp_path):
    header = tmp_path / "dense.header.pdf"
    _write_pdf_header(header, "dense", "banana")

    with pytest.raises(ValueError, match="invalid distribution='banana'"):
        read_pdf_header(str(header))


def test_modern_pdf_rejects_dense_materialized_peak_above_limit(tmp_path, monkeypatch):
    data = tmp_path / "dense.00000.pdf"
    _write_pdf_header(tmp_path / "dense.header.pdf", "dense")
    _write_dense_pdf(data, [1.0, 0.0, 0.0, 2.5])
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 100)

    with pytest.raises(ValueError, match="PDF dense reconstruction peak requires"):
        read_pdf(str(data))


def test_modern_pdf_rejects_sparse_materialized_peak_above_limit(tmp_path, monkeypatch):
    case = tmp_path / "sparse_peak"
    _write_pdf_header(case / "hist.header.pdf", "sparse_coo", "rank")
    shard = case / "rank_00000000" / "hist.00000.pdf"
    _write_sparse_pdf(shard, [(0, 1.0)])
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 150)

    with pytest.raises(ValueError, match="PDF sparse reconstruction peak requires"):
        read_pdf(str(shard))


def test_modern_pdf_counts_retained_reference_state_while_parsing_shard_header(
    tmp_path, monkeypatch
):
    case = tmp_path / "retained_header_peak"
    shard_dir = case / "rank_00000000"
    _write_pdf_header(shard_dir / "hist.header.pdf", "sparse_coo", "rank")
    shard = shard_dir / "hist.00000.pdf"
    _write_sparse_pdf(shard, [(0, 1.0)])
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 140)

    with pytest.raises(
        ValueError, match="PDF retained histogram and bin metadata requires"
    ):
        read_pdf(str(shard))


def test_modern_pdf_rejects_sparse_duplicate_validation_peak_above_limit(
    tmp_path, monkeypatch
):
    case = tmp_path / "sparse_unique_peak"
    _write_pdf_header(case / "hist.header.pdf", "sparse_coo", "rank")
    shard = case / "rank_00000000" / "hist.00000.pdf"
    _write_sparse_pdf(shard, [(0, 1.0), (1, 2.0)])
    monkeypatch.setattr(read_pdf_module, "_MAX_DENSE_ALLOCATION_BYTES", 185)

    with pytest.raises(ValueError, match="PDF sparse duplicate-validation peak requires"):
        read_pdf(str(shard))


def test_modern_pdf_releases_sparse_arrays_before_reading_next_shard(
    tmp_path, monkeypatch
):
    case = tmp_path / "sparse_release"
    _write_pdf_header(case / "hist.header.pdf", "sparse_coo", "rank")
    rank0 = case / "rank_00000000" / "hist.00000.pdf"
    rank1 = case / "rank_00000001" / "hist.00000.pdf"
    _write_sparse_pdf(rank0, [(0, 1.0)])
    _write_sparse_pdf(rank1, [(1, 2.0)])
    read_sparse_file = read_pdf_module._read_sparse_file
    previous_arrays = []

    def tracked_read_sparse_file(*args, **kwargs):
        assert all(reference() is None for reference in previous_arrays)
        result = read_sparse_file(*args, **kwargs)
        previous_arrays[:] = [weakref.ref(result[2]), weakref.ref(result[3])]
        return result

    monkeypatch.setattr(read_pdf_module, "_read_sparse_file", tracked_read_sparse_file)
    np.testing.assert_allclose(read_pdf(str(rank0))["pdf"], [1.0, 2.0, 0.0, 0.0])
    assert all(reference() is None for reference in previous_arrays)


def test_modern_pdf_releases_local_header_arrays_before_parsing_next_shard(
    tmp_path, monkeypatch
):
    case = tmp_path / "header_release"
    rank0_dir = case / "rank_00000000"
    rank1_dir = case / "rank_00000001"
    _write_pdf_header(rank0_dir / "hist.header.pdf", "sparse_coo", "rank")
    _write_pdf_header(rank1_dir / "hist.header.pdf", "sparse_coo", "rank")
    rank0 = rank0_dir / "hist.00000.pdf"
    rank1 = rank1_dir / "hist.00000.pdf"
    _write_sparse_pdf(rank0, [(0, 1.0)])
    _write_sparse_pdf(rank1, [(1, 2.0)])
    read_header = read_pdf_module._read_pdf_header
    previous_arrays = []
    header_calls = 0

    def tracked_read_header(*args, **kwargs):
        nonlocal header_calls
        if header_calls >= 2:
            assert all(reference() is None for reference in previous_arrays)
        result = read_header(*args, **kwargs)
        if header_calls >= 1:
            previous_arrays[:] = [
                weakref.ref(array)
                for dimension in result["dimensions"]
                for array in (dimension["bin_edges"], dimension["bin_centers"])
            ]
        header_calls += 1
        return result

    monkeypatch.setattr(read_pdf_module, "_read_pdf_header", tracked_read_header)
    np.testing.assert_allclose(read_pdf(str(rank0))["pdf"], [1.0, 2.0, 0.0, 0.0])
    assert all(reference() is None for reference in previous_arrays)


def test_modern_pdf_v2_sparse_requires_inventory_metadata(tmp_path):
    case = tmp_path / "missing_inventory"
    _write_pdf_header(case / "hist.header.pdf", "sparse_coo", "rank", v2=True)
    shard = case / "rank_00000000" / "hist.00000.pdf"
    _write_v2_sparse_pdf(shard, [(0, 1.0)])

    with pytest.raises(ValueError, match="missing required inventory metadata"):
        read_pdf(str(shard))


def test_modern_pdf_rejects_header_above_practical_read_limit(tmp_path, monkeypatch):
    header = tmp_path / "dense.header.pdf"
    _write_pdf_header(header, "dense")
    monkeypatch.setattr(read_pdf_module, "_MAX_HEADER_READ_BYTES", 16)

    with pytest.raises(ValueError, match="PDF header .* practical file-read limit"):
        read_pdf_header(str(header))


def test_modern_pdf_rejects_payload_above_practical_read_limit(tmp_path, monkeypatch):
    data = tmp_path / "dense.00000.pdf"
    _write_pdf_header(tmp_path / "dense.header.pdf", "dense")
    _write_dense_pdf(data, [1.0, 0.0, 0.0, 2.5])
    monkeypatch.setattr(read_pdf_module, "_MAX_PAYLOAD_READ_BYTES", 16)

    with pytest.raises(ValueError, match="PDF payload .* practical file-read limit"):
        read_pdf(str(data))


def test_modern_pdf_rejects_malformed_header(tmp_path):
    data = tmp_path / "bad.00000.pdf"
    header = tmp_path / "bad.header.pdf"
    header.write_text("format = dense\nndim = 1\n")
    _write_dense_pdf(data, [0.0, 0.0, 0.0, 0.0])

    with pytest.raises(ValueError, match="missing dimension"):
        read_pdf(str(data))


def _write_sphslice(path, distribution, indices, values, npoints=None, dump=b""):
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
        f"header_offset={len(dump)}",
    ]
    payload = b"\n".join(line.encode("ascii") for line in lines) + b"\n" + dump
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


def test_sphslice_releases_sparse_arrays_before_reading_next_shard(
    tmp_path, monkeypatch
):
    root = tmp_path / "release"
    rank0 = root / "rank_00000000" / "surface.00000.sphslice"
    rank1 = root / "rank_00000001" / "surface.00000.sphslice"
    _write_sphslice(rank0, "rank", [0], [10.0])
    _write_sphslice(rank1, "rank", [1], [20.0])
    read_one = read_sphslice_module._read_one
    previous_arrays = []
    calls = 0

    def tracked_read_one(*args, **kwargs):
        nonlocal calls
        assert all(reference() is None for reference in previous_arrays)
        retained_bytes = args[1] if len(args) > 1 else kwargs.get(
            "externally_retained_bytes", 0
        )
        if calls:
            assert retained_bytes > 0
        result = read_one(*args, **kwargs)
        previous_arrays[:] = [weakref.ref(result[1]), weakref.ref(result[2])]
        calls += 1
        return result

    monkeypatch.setattr(read_sphslice_module, "_read_one", tracked_read_one)
    np.testing.assert_allclose(
        read_sphslice(str(rank0))["data"].reshape(-1), [10.0, 20.0]
    )
    assert all(reference() is None for reference in previous_arrays)


def test_sphslice_rejects_embedded_dump_materialized_peak(tmp_path, monkeypatch):
    shared = tmp_path / "surface.00000.sphslice"
    _write_sphslice(shared, "shared", [], [10.0, 20.0], dump=b"x" * 100)
    monkeypatch.setattr(read_sphslice_module, "_variable_token_peak_bytes", lambda *_: 0)
    monkeypatch.setattr(read_sphslice_module, "_MAX_DENSE_ALLOCATION_BYTES", 50)

    with pytest.raises(ValueError, match="dense reconstruction peak requires"):
        read_sphslice(str(shared))


def test_sphslice_rejects_sparse_duplicate_validation_peak(tmp_path, monkeypatch):
    shard = tmp_path / "rank_00000000" / "surface.00000.sphslice"
    _write_sphslice(shard, "rank", [0, 1], [10.0, 20.0])
    monkeypatch.setattr(read_sphslice_module, "_variable_token_peak_bytes", lambda *_: 0)
    monkeypatch.setattr(read_sphslice_module, "_MAX_DENSE_ALLOCATION_BYTES", 80)

    with pytest.raises(
        ValueError, match="spherical-slice sparse duplicate-validation peak requires"
    ):
        read_sphslice(str(shard))


def test_sphslice_missing_ownership_diagnostic_avoids_full_index_array(
    tmp_path, monkeypatch
):
    root = tmp_path / "missing_bounded"
    rank0 = root / "rank_00000000" / "surface.00000.sphslice"
    rank1 = root / "rank_00000001" / "surface.00000.sphslice"
    _write_sphslice(rank0, "rank", [0], [10.0])
    _write_sphslice(rank1, "rank", [], [])
    monkeypatch.setattr(
        read_sphslice_module.np,
        "flatnonzero",
        lambda *_args, **_kwargs: pytest.fail("must not materialize missing indexes"),
    )

    with pytest.raises(ValueError, match="missing 1 of 2 angular points"):
        read_sphslice(str(rank0))


def test_sphslice_rejects_coordinate_generation_peak(tmp_path, monkeypatch):
    shared = tmp_path / "surface.00000.sphslice"
    _write_sphslice(shared, "shared", [], [10.0, 20.0])
    monkeypatch.setattr(read_sphslice_module, "_variable_token_peak_bytes", lambda *_: 0)
    monkeypatch.setattr(read_sphslice_module, "_MAX_DENSE_ALLOCATION_BYTES", 50)

    with pytest.raises(ValueError, match="coordinate generation peak requires"):
        read_sphslice(str(shared))


def test_sphslice_coordinate_peak_counts_retained_variable_metadata(
    tmp_path, monkeypatch
):
    shared = tmp_path / "surface.00000.sphslice"
    _write_sphslice(shared, "shared", [], [10.0, 20.0])
    monkeypatch.setattr(read_sphslice_module, "_variable_token_peak_bytes", lambda *_: 10)
    retained_dense_bytes = 2 * np.dtype(bool).itemsize
    retained_dense_bytes += 2 * np.dtype(np.float32).itemsize
    retained_dense_bytes += 3 * np.dtype(np.float64).itemsize
    coordinate_bytes = 4 * 3 * np.dtype(np.float64).itemsize
    monkeypatch.setattr(
        read_sphslice_module,
        "_MAX_DENSE_ALLOCATION_BYTES",
        retained_dense_bytes + coordinate_bytes,
    )

    with pytest.raises(ValueError, match="coordinate generation peak requires"):
        read_sphslice(str(shared))


def test_sphslice_rejects_header_line_above_practical_limit(tmp_path, monkeypatch):
    shared = tmp_path / "surface.00000.sphslice"
    _write_sphslice(shared, "shared", [], [10.0, 20.0])
    payload = shared.read_bytes().replace(b"variables: dens", b"variables: " + b"x" * 100)
    shared.write_bytes(payload)
    monkeypatch.setattr(read_sphslice_module, "_MAX_HEADER_READ_BYTES", 64)

    with pytest.raises(ValueError, match="practical metadata limit"):
        read_sphslice_header(str(shared))


def test_sphslice_preflights_variable_tokens_before_split(tmp_path, monkeypatch):
    shared = tmp_path / "surface.00000.sphslice"
    _write_sphslice(shared, "shared", [], [10.0, 20.0])
    monkeypatch.setattr(read_sphslice_module, "_MAX_DENSE_ALLOCATION_BYTES", 1)

    with pytest.raises(ValueError, match="variable tokenization peak requires"):
        read_sphslice_header(str(shared))


def test_sphslice_rejects_duplicate_variable_metadata(tmp_path):
    shared = tmp_path / "surface.00000.sphslice"
    _write_sphslice(shared, "shared", [], [10.0, 20.0])
    payload = shared.read_bytes().replace(
        b"variables: dens\n", b"variables: dens\nvariables: dens\n"
    )
    shared.write_bytes(payload)

    with pytest.raises(ValueError, match="duplicate spherical-slice variables"):
        read_sphslice_header(str(shared))


def test_sphslice_counts_retained_reference_variable_metadata(
    tmp_path, monkeypatch
):
    shared = tmp_path / "surface.00000.sphslice"
    _write_sphslice(shared, "shared", [], [10.0, 20.0])
    with shared.open("rb") as handle:
        header = read_sphslice_module._read_header(handle)
    monkeypatch.setattr(
        read_sphslice_module,
        "_MAX_DENSE_ALLOCATION_BYTES",
        header["_retained_metadata_bytes"],
    )

    with shared.open("rb") as handle:
        with pytest.raises(ValueError, match="variable tokenization peak requires"):
            read_sphslice_module._read_header(handle, externally_retained_bytes=1)
