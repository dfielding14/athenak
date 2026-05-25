"""MPI reproducibility regression for CGL Landau-fluid AMR evolution."""

from pathlib import Path

import numpy as np

import test_suite.testutils as testutils


INPUT_FILE = "../../../inputs/tests/cgl_lf_amr_2d.athinput"


def _run_amr(basename, threads):
    testutils.mpi_run(
        INPUT_FILE,
        [f"job/basename={basename}"],
        threads=threads,
    )
    return (
        testutils.athena_read.hst(f"{basename}.mhd.hst"),
        testutils.athena_read.hst(f"{basename}.user.hst"),
    )


def _assert_admissible(mhd, user):
    assert user["ncell"][-1] > user["ncell"][0]
    assert np.max(user["max_ndiv"]) < 1.0e-12
    assert np.max(user["bad_state"]) == 0.0
    assert np.all(np.isfinite(user["abs_anis"]))
    assert mhd["lf_nstage"][-1] > 0.0
    for column in (
        "lf_dfloor",
        "lf_pfloor",
        "lf_nonfin",
        "lf_nonpos",
        "lf_hardbd",
    ):
        assert mhd[column][-1] == 0.0
    assert mhd["lf_qface"][-1] > 0.0
    for column in ("lf_qprcap", "lf_qpr10", "lf_qpecap", "lf_qpe10"):
        assert 0.0 <= mhd[column][-1] <= mhd["lf_qface"][-1]
    assert abs(mhd["lf_qprwrk"][-1]) + abs(mhd["lf_qpewrk"][-1]) > 0.0
    for column in ("lf_cpwrk", "lf_cawrk"):
        assert np.isfinite(mhd[column][-1])
    assert abs(mhd["lf_cpwrk"][-1]) + abs(mhd["lf_cawrk"][-1]) > 0.0
    scale = max(abs(mhd["tot-E"][0]), 1.0e-30)
    residual = abs(mhd["tot-E"][-1] - mhd["tot-E"][0]) / scale
    assert residual < 5.0e-3


def test_cgl_lf_amr_is_reproducible_across_mpi_decomposition():
    try:
        single_mhd, single_user = _run_amr("cgl_mpi_single", threads=1)
        multi_mhd, multi_user = _run_amr("cgl_mpi_four", threads=4)
        _assert_admissible(single_mhd, single_user)
        _assert_admissible(multi_mhd, multi_user)

        for column in (
            "mass",
            "tot-E",
            "aam-D",
            "lf_nstage",
            "lf_qface",
            "lf_qprcap",
            "lf_qpr10",
            "lf_qpecap",
            "lf_qpe10",
            "lf_qprwrk",
            "lf_qpewrk",
            "lf_cpwrk",
            "lf_cawrk",
        ):
            assert np.isclose(
                single_mhd[column][-1],
                multi_mhd[column][-1],
                rtol=1.0e-12,
                atol=1.0e-12,
            )
        for column in ("ncell", "bad_state", "abs_anis"):
            assert np.isclose(
                single_user[column][-1],
                multi_user[column][-1],
                rtol=1.0e-12,
                atol=1.0e-12,
            )
    finally:
        for path in Path(".").glob("cgl_mpi_*.hst"):
            path.unlink()
        testutils.cleanup()
