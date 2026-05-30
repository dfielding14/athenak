"""MPI writer regressions for CP-04 output hardening."""

import os
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
INPUT_FILE = ROOT / "tst" / "inputs" / "io_node_sharding.athinput"
sys.path.insert(0, str(ROOT / "vis" / "python"))

import bin_convert  # noqa: E402
from read_sphslice import read_sphslice, read_sphslice_header  # noqa: E402


def _run(
    tmp_path: Path, name: str, *overrides: str, check=True, input_file=INPUT_FILE,
    nranks=2
):
    run_dir = tmp_path / name
    run_dir.mkdir()
    env = os.environ.copy()
    env["ATHENAK_TEST_MAX_MPI_BYTES"] = "7"
    return run_dir, subprocess.run(
        [
            "mpirun",
            "-np",
            str(nranks),
            "./athena",
            "-i",
            str(input_file),
            "-d",
            str(run_dir),
            *overrides,
        ],
        check=check,
        capture_output=True,
        text=True,
        env=env,
        timeout=90,
    )


def _compare_binary(shared: Path, sharded: Path, *, coarsened=False):
    reader = bin_convert.read_coarsened_binary if coarsened else bin_convert.read_binary
    reference = reader(str(shared))
    reconstructed = reader(str(sharded), assemble_shards=True)
    assert reference["var_names"] == reconstructed["var_names"]
    assert reference["n_mbs"] == reconstructed["n_mbs"]
    np.testing.assert_array_equal(reference["mb_logical"], reconstructed["mb_logical"])
    for variable in reference["var_names"]:
        np.testing.assert_allclose(
            reference["mb_data"][variable], reconstructed["mb_data"][variable]
        )


def test_forced_tiny_chunks_preserve_shared_rank_and_node_writers(tmp_path):
    shared, _ = _run(tmp_path, "shared")
    rank, _ = _run(
        tmp_path,
        "rank",
        "output1/single_file_per_rank=true",
        "output2/single_file_per_rank=true",
        "output3/single_file_per_rank=true",
        "output5/single_file_per_rank=true",
    )
    node, _ = _run(
        tmp_path,
        "node",
        "output1/single_file_per_node=true",
        "output2/single_file_per_node=true",
        "output3/single_file_per_node=true",
        "output5/single_file_per_node=true",
    )

    shared_bin = shared / "bin" / "io_node.full.00000.bin"
    rank_bin = rank / "bin" / "rank_00000000" / "io_node.full.00000.bin"
    node_bin = node / "bin" / "node_00000000" / "io_node.full.00000.bin"
    _compare_binary(shared_bin, rank_bin)
    _compare_binary(shared_bin, node_bin)

    shared_slice = shared / "bin" / "io_node.slice.00000.bin"
    rank_slice = rank / "bin" / "rank_00000000" / "io_node.slice.00000.bin"
    _compare_binary(shared_slice, rank_slice)
    rank_slice_shards = [
        bin_convert.read_binary(
            str(rank / "bin" / f"rank_{rank_id:08d}" / "io_node.slice.00000.bin")
        )
        for rank_id in range(2)
    ]
    assert sorted(shard["n_mbs"] for shard in rank_slice_shards) == [0, 1]

    shared_cbin = shared / "cbin_coarse_2" / "io_node.coarse.00000.cbin"
    rank_cbin = rank / "cbin_coarse_2" / "rank_00000000" / "io_node.coarse.00000.cbin"
    node_cbin = node / "cbin_coarse_2" / "node_00000000" / "io_node.coarse.00000.cbin"
    _compare_binary(shared_cbin, rank_cbin, coarsened=True)
    _compare_binary(shared_cbin, node_cbin, coarsened=True)

    sph_name = "io_node.density.r_2.5000000000000000e-01.00000.sph.bin"
    shared_sph = shared / "bin" / sph_name
    rank_sph = rank / "bin" / "rank_00000000" / sph_name
    node_sph = node / "bin" / "node_00000000" / sph_name
    np.testing.assert_allclose(read_sphslice(str(rank_sph))["data"],
                               read_sphslice(str(shared_sph))["data"])
    np.testing.assert_allclose(read_sphslice(str(node_sph))["data"],
                               read_sphslice(str(shared_sph))["data"])

    assert bin_convert.read_binary(str(node_bin))["number_of_nodes"] == 1
    assert bin_convert.read_coarsened_binary(str(node_cbin))["number_of_nodes"] == 1
    assert read_sphslice_header(str(rank_sph))["number_of_ranks"] == 2
    assert read_sphslice_header(str(node_sph))["number_of_nodes"] == 1
    assert "dens" in bin_convert.read_rank_binary_as_athdf(str(node_bin))
    for run_dir in (shared, rank, node):
        assert not list(run_dir.rglob("*.tmp"))


def test_coarsened_binary_moments_round_trip_across_shared_rank_and_node(tmp_path):
    shared, _ = _run(tmp_path, "moments_shared", "output3/compute_moments=true")
    rank, _ = _run(
        tmp_path,
        "moments_rank",
        "output3/compute_moments=true",
        "output3/single_file_per_rank=true",
    )
    node, _ = _run(
        tmp_path,
        "moments_node",
        "output3/compute_moments=true",
        "output3/single_file_per_node=true",
    )

    shared_cbin = shared / "cbin_coarse_2" / "io_node.coarse.00000.cbin"
    rank_cbin = rank / "cbin_coarse_2" / "rank_00000000" / "io_node.coarse.00000.cbin"
    node_cbin = node / "cbin_coarse_2" / "node_00000000" / "io_node.coarse.00000.cbin"
    _compare_binary(shared_cbin, rank_cbin, coarsened=True)
    _compare_binary(shared_cbin, node_cbin, coarsened=True)

    shared_data = bin_convert.read_coarsened_binary(str(shared_cbin))
    assert shared_data["number_of_moments"] == 4
    assert shared_data["var_names"] == [
        "dens_1st",
        "dens_2nd",
        "dens_3rd",
        "dens_4th",
    ]
    expected_uov = np.asarray(
        [shared_data["mb_data"][variable] for variable in shared_data["var_names"]]
    )
    for path, assemble in ((shared_cbin, False), (rank_cbin, True), (node_cbin, True)):
        bin_convert.convert_file(str(path), assemble_shards=assemble)
        athdf = path.with_suffix(".athdf")
        assert athdf.exists()
        xdmf = Path(str(athdf) + ".xdmf")
        assert xdmf.exists()
        with h5py.File(athdf) as handle:
            variable_names = [
                name.decode() for name in handle.attrs["VariableNames"]
            ]
            assert variable_names == shared_data["var_names"]
            assert handle["uov"].shape == expected_uov.shape
            np.testing.assert_allclose(handle["uov"][:], expected_uov)
        xdmf_text = xdmf.read_text()
        for variable in shared_data["var_names"]:
            assert f'<Attribute Name="{variable}" Center="Cell">' in xdmf_text
        assert xdmf_text.count(":/uov") == (
            shared_data["n_mbs"] * len(shared_data["var_names"])
        )


def test_coarsened_binary_gid_filter_uses_global_meshblock_ids(tmp_path):
    gid_input = tmp_path / "cbin_gid.athinput"
    gid_input.write_text(
        INPUT_FILE.read_text().replace(
            "<output3>\nfile_type = cbin\n",
            "<output3>\nfile_type = cbin\ngid = 1\n",
        )
    )
    rank, _ = _run(
        tmp_path,
        "cbin_gid",
        "output3/single_file_per_rank=true",
        input_file=gid_input,
    )
    rank_one = (
        rank / "cbin_coarse_2" / "rank_00000001" / "io_node.coarse.00000.cbin"
    )
    rank_zero = (
        rank / "cbin_coarse_2" / "rank_00000000" / "io_node.coarse.00000.cbin"
    )
    data = bin_convert.read_coarsened_binary(str(rank_one))
    empty = bin_convert.read_coarsened_binary(str(rank_zero))
    assembled = bin_convert.read_coarsened_binary(str(rank_one), assemble_shards=True)
    assert empty["n_mbs"] == 0
    assert data["n_mbs"] == 1
    assert assembled["n_mbs"] == 1
    np.testing.assert_array_equal(data["mb_logical"][0], [1, 0, 0, 0])
    np.testing.assert_array_equal(assembled["mb_logical"][0], [1, 0, 0, 0])


@pytest.mark.parametrize("nranks", (1, 2))
@pytest.mark.parametrize(
    "distribution",
    (
        "output3/single_file_per_rank=true",
        "output3/single_file_per_node=true",
    ),
)
def test_sliced_sharded_coarsened_binary_remains_explicitly_unpromoted(
    tmp_path, nranks, distribution
):
    sliced_input = tmp_path / "sliced_sharded_cbin.athinput"
    sliced_input.write_text(
        INPUT_FILE.read_text().replace(
            "<output3>\nfile_type = cbin\n",
            "<output3>\nfile_type = cbin\nslice_x1 = 0.25\n",
        )
    )
    _, proc = _run(
        tmp_path,
        f"sliced_sharded_cbin_{nranks}",
        distribution,
        check=False,
        input_file=sliced_input,
        nranks=nranks,
    )
    assert proc.returncode != 0
    assert "Sliced coarsened-binary output is not supported." in (
        proc.stdout + proc.stderr
    )
