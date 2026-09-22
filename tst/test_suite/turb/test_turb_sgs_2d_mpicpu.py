"""MPI regression for fully 2D viscous turbulence and SGS output."""

import numpy as np
import pytest

from .test_turb_sgs_2d_cpu import (INPUT, check_favre_filter, latest,
                                  require_success, run_athena)

from bin_convert import (  # noqa: E402
    read_binary,
    read_coarsened_binary,
    read_all_ranks_coarsened_binary,
)


@pytest.mark.parametrize("factor, per_rank", [(2, False), (8, False), (8, True)])
def test_viscous_2d_sgs_output_under_mpi(tmp_path, factor, per_rank):
    """Four ranks preserve the 2D contract and write the viscous SGS state."""
    output_dir = tmp_path / "mpi_viscous_sgs"
    output_dir.mkdir()
    input_file = output_dir / "case.athinput"
    input_file.write_text(INPUT.read_text() +
                         f"\n<output4>\nsingle_file_per_rank = {str(per_rank).lower()}\n")
    require_success(run_athena(
        output_dir, f"output4/coarsen_factor={factor}", input_file=input_file,
        launcher=("mpirun", "-np", "4"),
    ))

    state = read_binary(str(latest(output_dir / "bin", "*.state.*.bin")))
    force = read_binary(str(latest(output_dir / "bin", "*.force.*.bin")))
    sgs_dir = output_dir / f"cbin_sgs_{factor}"
    reader = read_coarsened_binary
    if per_rank:
        sgs_dir /= "rank_00000000"
        reader = read_all_ranks_coarsened_binary
    sgs = reader(str(latest(sgs_dir, "*.sgs.*.cbin")))
    assert state["n_mbs"] == 4
    assert sgs["n_mbs"] == 4
    assert sgs["Nx3"] == sgs["nx3_mb"] == 1
    assert sgs["nx1_mb"] == sgs["nx2_mb"] == 8 // factor
    np.testing.assert_array_equal(state["mb_logical"], sgs["mb_logical"])
    assert sgs["var_names"] == [
        "dens",
        "velx",
        "vely",
        "tau_xx",
        "tau_xy",
        "tau_yy",
    ]
    assert np.all(np.asarray(force["mb_data"]["force3"]) == 0.0)
    assert np.all(np.asarray(state["mb_data"]["mom3"]) == 0.0)

    check_favre_filter(state, sgs, factor)


def test_sparse_2d_restart_with_changed_mpi_rank_count(tmp_path):
    """A four-to-two-rank restart preserves the driven state and SGS products."""
    input_file = tmp_path / "restart.athinput"
    input_file.write_text(INPUT.read_text().replace(
        "normalization = accel_rms\naccel_rms = 0.2",
        "normalization = edot\ndedt = 0.001",
    ) + "\n<hydro_srcterms>\nlinear_drag = true\ndrag_rate = 0.1\n"
        "\n<output5>\nfile_type = rst\ndcycle = 10\n")
    overrides = (
        "mesh/nghost=3", "mesh/nx1=64", "mesh/nx2=64",
        "meshblock/nx1=32", "meshblock/nx2=32",
        "time/integrator=rk3", "time/tlim=10", "hydro/reconstruct=wenoz",
        "turb_driving/mode_sampling=sparse_annulus",
        "turb_driving/sparse_mode_count=16", "turb_driving/nlow=5",
        "turb_driving/nhigh=7", "turb_driving/npeak=6",
        "turb_driving/dt_update=0.001", "output4/coarsen_factor=32",
    )

    def run(label, ranks, nlim, restart=None):
        output_dir = tmp_path / label
        require_success(run_athena(
            output_dir, *overrides, f"time/nlim={nlim}", input_file=input_file,
            launcher=("mpirun", "-np", str(ranks)), restart=restart,
        ))
        return output_dir

    reference = run("reference", 4, 20)
    split = run("split", 4, 10)
    resumed = run("resumed", 2, 20, latest(split / "rst", "*.rst"))
    for directory, pattern, reader in (
        ("bin", "*.state.*.bin", read_binary),
        ("bin", "*.force.*.bin", read_binary),
        ("cbin_sgs_32", "*.sgs.*.cbin", read_coarsened_binary),
    ):
        expected = reader(str(latest(reference / directory, pattern)))
        actual = reader(str(latest(resumed / directory, pattern)))
        assert actual["time"] == expected["time"]
        np.testing.assert_array_equal(actual["mb_logical"], expected["mb_logical"])
        for name in expected["var_names"]:
            np.testing.assert_array_equal(actual["mb_data"][name],
                                          expected["mb_data"][name])
