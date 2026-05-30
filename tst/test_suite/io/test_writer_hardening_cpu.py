"""Pure-reader regressions for CP-04 writer hardening contracts."""

from io import BytesIO
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tst" / "fixtures" / "io" / "origin_main_886dd2a1"
sys.path.insert(0, str(ROOT / "vis" / "python"))

import bin_convert  # noqa: E402
from bin_convert import read_binary, read_coarsened_binary  # noqa: E402
from read_sphslice import read_sphslice  # noqa: E402


def _binary_fixture(kind, rank):
    name = f"io_legacy_per_rank.{kind}_per_rank.00000.{kind}"
    return FIXTURES / kind / "per_rank" / f"rank_{rank:08d}" / name


def _binary_payload_offset(payload):
    fp = BytesIO(payload)
    fp.readline()
    pheader_count = int(fp.readline().split(b"=")[-1])
    for _ in range(pheader_count - 1):
        fp.readline()
    fp.readline()
    fp.readline()
    header_size = int(fp.readline().split(b"=")[-1])
    fp.seek(header_size, 1)
    return fp.tell()


def _make_node_binary(path, source, node, nnodes, empty=False):
    payload = source.read_bytes()
    fp = BytesIO(payload)
    first = fp.readline()
    count_line = fp.readline()
    pheader_count = int(count_line.split(b"=")[-1])
    pheader = [fp.readline() for _ in range(pheader_count - 1)]
    remainder = fp.read()
    additions = [
        b"  distribution=node\n",
        f"  node={node}\n".encode("ascii"),
        f"  number of nodes={nnodes}\n".encode("ascii"),
        f"  number of meshblocks={0 if empty else 1}\n".encode("ascii"),
    ]
    rewritten = (
        first
        + f"  size of preheader={pheader_count + len(additions)}\n".encode("ascii")
        + b"".join(pheader)
        + b"".join(additions)
        + remainder
    )
    if empty:
        rewritten = rewritten[: _binary_payload_offset(rewritten)]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(rewritten)


@pytest.mark.parametrize(
    ("kind", "reader"),
    (("bin", read_binary), ("cbin", read_coarsened_binary)),
)
def test_node_binary_inventory_accepts_explicit_empty_shards_and_rejects_gaps(
    tmp_path, kind, reader
):
    root = tmp_path / kind
    node0 = root / "node_00000000" / f"output.00000.{kind}"
    node1 = root / "node_00000001" / f"output.00000.{kind}"
    _make_node_binary(node0, _binary_fixture(kind, 0), 0, 2)
    _make_node_binary(node1, _binary_fixture(kind, 1), 1, 2, empty=True)

    assembled = reader(str(node0), assemble_shards=True)
    assert assembled["n_mbs"] == 1
    assert assembled["number_of_nodes"] == 2

    node1.unlink()
    with pytest.raises(ValueError, match="node shard inventory is incomplete"):
        reader(str(node0), assemble_shards=True)


@pytest.mark.parametrize(
    ("kind", "reader"),
    (("bin", read_binary), ("cbin", read_coarsened_binary)),
)
def test_node_binary_inventory_aggregates_meshblock_counts(tmp_path, kind, reader):
    root = tmp_path / kind
    node0 = root / "node_00000000" / f"output.00000.{kind}"
    node1 = root / "node_00000001" / f"output.00000.{kind}"
    _make_node_binary(node0, _binary_fixture(kind, 0), 0, 2)
    _make_node_binary(node1, _binary_fixture(kind, 1), 1, 2)

    assembled = reader(str(node0), assemble_shards=True)
    assert assembled["n_mbs"] == 2
    assert assembled["number_of_meshblocks"] == 2


@pytest.mark.parametrize(
    ("kind", "reader"),
    (("bin", read_binary), ("cbin", read_coarsened_binary)),
)
def test_node_binary_inventory_rejects_oversized_declared_count(tmp_path, kind, reader):
    path = tmp_path / kind / "node_00000000" / f"output.00000.{kind}"
    _make_node_binary(path, _binary_fixture(kind, 0), 0, 10**12)

    with pytest.raises(ValueError, match="node shard inventory is incomplete"):
        reader(str(path), assemble_shards=True)


@pytest.mark.parametrize(
    ("kind", "reader"),
    (("bin", read_binary), ("cbin", read_coarsened_binary)),
)
def test_node_binary_inventory_rejects_duplicate_integer_ids(tmp_path, kind, reader):
    root = tmp_path / kind
    canonical = root / "node_00000000" / f"output.00000.{kind}"
    duplicate = root / "node_0" / f"output.00000.{kind}"
    _make_node_binary(canonical, _binary_fixture(kind, 0), 0, 1)
    _make_node_binary(duplicate, _binary_fixture(kind, 1), 0, 1)

    with pytest.raises(ValueError, match="duplicate node IDs"):
        reader(str(canonical), assemble_shards=True)


def test_node_binary_rejects_header_path_id_mismatch(tmp_path):
    path = tmp_path / "node_00000001" / "output.00000.bin"
    _make_node_binary(path, _binary_fixture("bin", 0), 0, 1)

    with pytest.raises(ValueError, match="declares node=0"):
        read_binary(str(path))


def test_rank_binary_as_athdf_wrapper_accepts_node_shard(tmp_path):
    path = tmp_path / "node_00000000" / "output.00000.bin"
    _make_node_binary(path, _binary_fixture("bin", 0), 0, 1)

    converted = bin_convert.read_rank_binary_as_athdf(str(path))
    assert "dens" in converted
    assert converted["dens"].size > 0


def _write_sphslice(
    path,
    distribution,
    indices,
    values,
    *,
    shard_id=None,
    sibling_count=None,
    layout=None,
):
    path.parent.mkdir(parents=True, exist_ok=True)
    if layout is None:
        layout = "dense" if distribution == "shared" else "sparse_angles"
    lines = [
        "Athena spherical slice version=1.0",
        f"layout={layout}",
        f"distribution={distribution}",
        "time=0.25",
        "cycle=5",
        "radius=2.0",
        "ntheta=1",
        "nphi=2",
        "size_of_variable=4",
        "number_of_variables=1",
        f"npoints={2 if distribution == 'shared' else len(indices)}",
    ]
    if distribution == "rank" and shard_id is not None:
        lines += [f"rank={shard_id}", f"number of ranks={sibling_count}"]
    if distribution == "node" and shard_id is not None:
        lines += [f"node={shard_id}", f"number of nodes={sibling_count}"]
    lines += ["variables: dens", "header_offset=0"]
    payload = ("\n".join(lines) + "\n").encode("ascii")
    if distribution == "shared":
        payload += np.asarray(values, dtype=np.float32).tobytes()
    else:
        payload += np.asarray(indices, dtype=np.int32).tobytes()
        payload += np.asarray(values, dtype=np.float32).tobytes()
    path.write_bytes(payload)


def test_sphslice_node_inventory_accepts_explicit_empty_shard_and_rejects_gap(tmp_path):
    shared = tmp_path / "shared" / "surface.00000.sph.bin"
    node0 = tmp_path / "nodes" / "node_00000000" / "surface.00000.sph.bin"
    node1 = tmp_path / "nodes" / "node_00000001" / "surface.00000.sph.bin"
    _write_sphslice(shared, "shared", [], [10.0, 20.0])
    _write_sphslice(node0, "node", [0, 1], [10.0, 20.0], shard_id=0, sibling_count=2)
    _write_sphslice(node1, "node", [], [], shard_id=1, sibling_count=2)

    np.testing.assert_allclose(read_sphslice(str(node0))["data"],
                               read_sphslice(str(shared))["data"])
    node1.unlink()
    with pytest.raises(ValueError, match="node shard inventory is incomplete"):
        read_sphslice(str(node0))


@pytest.mark.parametrize(
    ("distribution", "layout", "directory"),
    (("shared", "sparse_angles", "shared"), ("rank", "dense", "rank_00000000")),
)
def test_sphslice_rejects_layout_distribution_mismatch(
    tmp_path, distribution, layout, directory
):
    path = tmp_path / directory / "surface.00000.sph.bin"
    kwargs = {}
    if distribution == "rank":
        kwargs = {"shard_id": 0, "sibling_count": 1}
    _write_sphslice(path, distribution, [0, 1], [10.0, 20.0], layout=layout, **kwargs)

    with pytest.raises(ValueError, match="declares layout"):
        read_sphslice(str(path))


def test_sphslice_rejects_rank_header_path_id_mismatch(tmp_path):
    path = tmp_path / "rank_00000001" / "surface.00000.sph.bin"
    _write_sphslice(path, "rank", [0, 1], [10.0, 20.0], shard_id=0, sibling_count=1)

    with pytest.raises(ValueError, match="declares rank=0"):
        read_sphslice(str(path))


def test_sphslice_rejects_malformed_sibling_directory(tmp_path):
    node0 = tmp_path / "node_00000000" / "surface.00000.sph.bin"
    malformed = tmp_path / "node_bad" / "surface.00000.sph.bin"
    _write_sphslice(node0, "node", [0, 1], [10.0, 20.0], shard_id=0, sibling_count=1)
    _write_sphslice(malformed, "node", [], [], shard_id=0, sibling_count=1)

    with pytest.raises(ValueError, match="invalid spherical-slice shard directory"):
        read_sphslice(str(node0))


def test_sphslice_rejects_duplicate_integer_sibling_ids(tmp_path):
    canonical = tmp_path / "node_00000000" / "surface.00000.sph.bin"
    duplicate = tmp_path / "node_0" / "surface.00000.sph.bin"
    _write_sphslice(canonical, "node", [0], [10.0], shard_id=0, sibling_count=1)
    _write_sphslice(duplicate, "node", [1], [20.0], shard_id=0, sibling_count=1)

    with pytest.raises(ValueError, match="duplicate IDs"):
        read_sphslice(str(canonical))


def test_sphslice_rejects_oversized_declared_sibling_count(tmp_path):
    path = tmp_path / "node_00000000" / "surface.00000.sph.bin"
    _write_sphslice(path, "node", [0, 1], [10.0, 20.0], shard_id=0,
                    sibling_count=10**12)

    with pytest.raises(ValueError, match="node shard inventory is incomplete"):
        read_sphslice(str(path))
