"""MPI reproducibility regression for CGL Landau-fluid AMR evolution."""

from pathlib import Path

import numpy as np

import test_suite.testutils as testutils


INPUT_FILE = "../../../inputs/tests/cgl_lf_amr_2d.athinput"
OBLIQUE_INPUT = "../../../inputs/unit_tests/cgl_pure_paper_oblique_wave.athinput"


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


def test_cgl_lf_quantitative_projection_is_global_across_mpi_ranks():
    flags = ["meshblock/nx1=32"]
    testutils.mpi_run(
        OBLIQUE_INPUT,
        ["job/basename=cgl_mpi_projection_single", *flags],
        threads=1,
    )
    testutils.mpi_run(
        OBLIQUE_INPUT,
        ["job/basename=cgl_mpi_projection_four", *flags],
        threads=4,
    )


def test_cgl_lf_post_sweep_timestep_refresh_agrees_across_mpi_ranks():
    input_file = "../../../inputs/unit_tests/cgl_lf_timestep_refresh.athinput"
    try:
        # Task1's aligned row gives dt_FE=h^2/(2*chi), chi=sqrt(8/pi)/(2*pi).
        # The cap is 20*0.9*dt_FE; p_after=1+(2/3)*heating_rate*dt_cycle.
        # A half-sweep ratio 10*sqrt(p_after) requires 7 or 13 RKL2 stages.
        for heating_rate, post_stages in ((0, 7), (3000, 13)):
            histories = []
            for nranks in (1, 4):
                Path("cgl_lf_timestep_refresh.mhd.hst").unlink(missing_ok=True)
                testutils.mpi_run(
                    input_file, [f"problem/heating_rate={heating_rate}"], threads=nranks
                )
                history = testutils.athena_read.hst("cgl_lf_timestep_refresh.mhd.hst")
                assert history["lf_nstage"][-1] == 64 * (7 + post_stages)
                histories.append(history)
            for name in ("time", "dt", "tot-E", "lf_nstage"):
                np.testing.assert_allclose(histories[0][name], histories[1][name],
                                           rtol=2.0e-12, atol=0.0)
    finally:
        Path("cgl_lf_timestep_refresh.mhd.hst").unlink(missing_ok=True)


def test_cgl_lf_oblique_decay_agrees_across_mpi_ranks():
    from test_suite.cgl.test_cgl_lf_oblique_decay_cpu import (
        run_oblique_decay, assert_oblique_agreement,
    )
    for axis in ("x", "y"):
        one_block = run_oblique_decay(axis, 64, nranks=1)
        four_blocks = run_oblique_decay(axis, 32, nranks=1)
        four_ranks = run_oblique_decay(axis, 32, nranks=4)
        assert_oblique_agreement(one_block, four_blocks)
        assert_oblique_agreement(one_block, four_ranks)


def test_cgl_lf_smr_decay_conserves_energy_across_mpi_ranks():
    from test_suite.cgl.test_cgl_lf_oblique_decay_cpu import run_smr_decay
    single = run_smr_decay(nranks=1)
    multi = run_smr_decay(nranks=4)
    for name in ("time", "dt", "tot-E", "lf_nstage"):
        np.testing.assert_allclose(single[name], multi[name], rtol=2.0e-13, atol=0)
