"""MPI interaction regression for CGL Landau fluid with turbulence driving."""

from pathlib import Path
import shutil

import numpy as np
import pytest

import test_suite.testutils as testutils
from test_suite.turb.test_turb_driving_cpu import (
    read_force_blocks,
    rms_acceleration,
)


INPUT_FILE = "../../../inputs/tests/cgl_lf_turb_driving_amr.athinput"


def _run(basename, threads, prolong_primitives):
    testutils.mpi_run(
        INPUT_FILE,
        [f"job/basename={basename}",
         f"mesh_refinement/prolong_primitives={str(prolong_primitives).lower()}"],
        threads=threads,
    )
    return (
        testutils.athena_read.hst(f"{basename}.mhd.hst"),
        testutils.athena_read.hst(f"{basename}.user.hst"),
    )


def _latest_force(basename):
    paths = sorted(Path("bin").glob(f"{basename}.force.*.bin"))
    assert paths, f"no rendered turbulence output found for {basename}"
    return paths[-1]


def _assert_clean(mhd, user, basename):
    assert user["ncell"][-1] > user["ncell"][0]
    assert np.max(user["max_ndiv"]) < 2.0e-12
    assert np.max(user["bad_state"]) == 0.0
    assert mhd["lf_nstage"][-1] > 0.0
    for column in (
        "lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd",
    ):
        assert mhd[column][-1] == 0.0
    assert np.isclose(
        rms_acceleration(read_force_blocks(_latest_force(basename))),
        0.02,
        rtol=2.0e-6,
    )


@pytest.mark.parametrize("prolong_primitives", [False, True])
def test_driven_cgl_lf_amr_is_reproducible_across_mpi_decomposition(prolong_primitives):
    """Driven strict LF evolution is clean and reduction-stable under MPI."""
    try:
        single_mhd, single_user = _run("cgl_turb_mpi_single", 1, prolong_primitives)
        multi_mhd, multi_user = _run("cgl_turb_mpi_four", 4, prolong_primitives)
        _assert_clean(single_mhd, single_user, "cgl_turb_mpi_single")
        _assert_clean(multi_mhd, multi_user, "cgl_turb_mpi_four")

        for column in ("mass", "tot-E", "aam-D", "lf_nstage"):
            assert np.isclose(
                single_mhd[column][-1],
                multi_mhd[column][-1],
                rtol=1.0e-12,
                atol=1.0e-12,
            )
    finally:
        for path in Path(".").glob("cgl_turb_mpi_*.hst"):
            path.unlink()
        shutil.rmtree("bin", ignore_errors=True)
        shutil.rmtree("rst", ignore_errors=True)
        testutils.cleanup()


@pytest.mark.parametrize("model", ["active", "passive"])
def test_forcing_work_ledger_closes_across_mpi_decomposition(tmp_path, model):
    """The first uniform-state kick records global total/kinetic energy gain."""
    input_file = (Path(__file__).resolve().parents[3] / "inputs" / "cgl_lf_paper" /
                  f"cgl_lf_paper_smoke_{model}_beta10.athinput")
    work_by_ranks = []
    for ranks in (1, 4):
        output_dir = tmp_path / f"ranks_{ranks}"
        output_dir.mkdir()
        assert testutils.run_command([
            "mpirun", "-np", str(ranks), "./athena", "-i", str(input_file),
            "-d", str(output_dir), "job/basename=forcing_work",
            "meshblock/nx1=4", "meshblock/nx2=4", "time/nlim=1",
        ])
        mhd = testutils.athena_read.hst(str(output_dir / "forcing_work.mhd.hst"))
        user = testutils.athena_read.hst(str(output_dir / "forcing_work.user.hst"))
        work = user["force_work"] - user["force_work"][0]
        energy = user["kinetic"] if model == "passive" else mhd["tot-E"]
        assert work[-1] > 0.0
        np.testing.assert_allclose(work, energy - energy[0], rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(work, 0.32 * 2.0 * (user["time"] - user["time"][0]),
                                   rtol=1e-10, atol=1e-12)
        work_by_ranks.append(work)
    np.testing.assert_allclose(*work_by_ranks, rtol=1e-12, atol=1e-12)
