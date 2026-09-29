"""MPI interaction regression for CGL Landau fluid with turbulence driving."""

from pathlib import Path
import shutil

import numpy as np

import test_suite.testutils as testutils
from test_suite.turb.test_turb_driving_cpu import (
    read_force_blocks,
    rms_acceleration,
)


INPUT_FILE = "../../../inputs/tests/cgl_lf_turb_driving_amr.athinput"


def _run(basename, threads):
    testutils.mpi_run(
        INPUT_FILE,
        [f"job/basename={basename}"],
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


def test_driven_cgl_lf_amr_is_reproducible_across_mpi_decomposition():
    """Driven strict LF evolution is clean and reduction-stable under MPI."""
    try:
        single_mhd, single_user = _run("cgl_turb_mpi_single", threads=1)
        multi_mhd, multi_user = _run("cgl_turb_mpi_four", threads=4)
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
