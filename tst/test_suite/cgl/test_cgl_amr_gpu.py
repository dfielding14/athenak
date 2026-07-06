"""GPU smoke coverage for CGL-aware AMR primitive transfers."""

from pathlib import Path
import shutil
import subprocess

import numpy as np

import test_suite.testutils as testutils


INPUT_ROOT = "../../../inputs/tests"


def _run(input_name, basename, *flags):
    testutils.run(
        f"{INPUT_ROOT}/{input_name}",
        [f"job/basename={basename}", *flags],
    )


def _cleanup(prefix="cgl_amr_gpu"):
    for path in Path(".").glob(f"{prefix}*.hst"):
        path.unlink()
    shutil.rmtree("rst", ignore_errors=True)
    testutils.cleanup()


def _user_history(basename):
    return testutils.athena_read.hst(f"{basename}.user.hst")


def _mhd_history(basename):
    return testutils.athena_read.hst(f"{basename}.mhd.hst")


def _restart(restart_file, basename, *flags):
    command = ["./athena", "-r", str(restart_file), f"job/basename={basename}", *flags]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(result.stdout + result.stderr)


def _assert_clean_user(history, max_ndiv=1.0e-10):
    assert np.max(history["bad_state"]) == 0.0
    assert np.all(np.isfinite(history["abs_anis"]))
    assert np.max(history["max_ndiv"]) < max_ndiv


def _assert_no_amr_repairs(history):
    for column in (
        "amr_cell",
        "amr_nan",
        "amr_rho",
        "amr_U",
        "amr_ppar",
        "amr_pperp",
        "amr_lowB",
        "amr_fh",
        "amr_mirr",
        "amr_dlt",
        "amr_int",
        "amr_slp",
    ):
        assert history[column][-1] == 0.0


def _assert_clean_lf(history):
    assert history["lf_nstage"][-1] > 0.0
    for column in (
        "lf_dfloor",
        "lf_pfloor",
        "lf_nonfin",
        "lf_nonpos",
        "lf_hardbd",
    ):
        assert history[column][-1] == 0.0


def test_cgl_amr_primitive_smooth_gpu():
    try:
        _run("cgl_amr_primitive_uniform.athinput", "cgl_amr_gpu_uniform")
        uniform = _user_history("cgl_amr_gpu_uniform")
        _assert_clean_user(uniform)
        _assert_no_amr_repairs(uniform)
        assert uniform["ncell"][-1] > uniform["ncell"][0]

        _run(
            "cgl_amr_primitive_strong_anisotropy.athinput",
            "cgl_amr_gpu_aniso",
        )
        aniso = _user_history("cgl_amr_gpu_aniso")
        _assert_clean_user(aniso)
        _assert_no_amr_repairs(aniso)
        assert np.max(aniso["ncell"]) > aniso["ncell"][0]
    finally:
        _cleanup()


def test_cgl_amr_current_churn_gpu():
    try:
        _run("cgl_amr_primitive_current_churn.athinput", "cgl_amr_gpu_churn")
        history = _user_history("cgl_amr_gpu_churn")
        _assert_clean_user(history)
        _assert_no_amr_repairs(history)
        assert np.max(history["ncell"]) > history["ncell"][0]
        assert history["ncell"][-1] < np.max(history["ncell"])
    finally:
        _cleanup()


def test_cgl_amr_restart_through_regrid_gpu():
    try:
        _run(
            "cgl_amr_primitive_current_churn_restart.athinput",
            "cgl_amr_gpu_restart_ref",
        )
        _run(
            "cgl_amr_primitive_current_churn_restart.athinput",
            "cgl_amr_gpu_restart_split",
            "time/nlim=3",
        )
        restarts = sorted(Path("rst").rglob("cgl_amr_gpu_restart_split*.rst"))
        assert restarts
        _restart(restarts[-1], "cgl_amr_gpu_restart_resume", "time/nlim=8")

        reference = _user_history("cgl_amr_gpu_restart_ref")
        resumed = _user_history("cgl_amr_gpu_restart_resume")
        _assert_clean_user(reference)
        _assert_clean_user(resumed)
        _assert_no_amr_repairs(reference)
        _assert_no_amr_repairs(resumed)
        assert np.max(reference["ncell"]) > reference["ncell"][0]
        assert reference["ncell"][-1] < np.max(reference["ncell"])
        assert np.isclose(resumed["time"][-1], reference["time"][-1], atol=1.0e-5)
        assert resumed["ncell"][-1] == reference["ncell"][-1]
    finally:
        _cleanup()


def test_cgl_amr_projection_stress_gpu():
    try:
        _run("cgl_amr_primitive_low_b.athinput", "cgl_amr_gpu_lowb")
        lowb = _user_history("cgl_amr_gpu_lowb")
        _assert_clean_user(lowb)
        assert lowb["amr_lowB"][-1] > 0.0

        _run("cgl_amr_primitive_firehose.athinput", "cgl_amr_gpu_firehose")
        firehose = _user_history("cgl_amr_gpu_firehose")
        _assert_clean_user(firehose)
        assert firehose["amr_fh"][-1] > 0.0

        _run("cgl_amr_primitive_mirror.athinput", "cgl_amr_gpu_mirror")
        mirror = _user_history("cgl_amr_gpu_mirror")
        _assert_clean_user(mirror)

        _run(
            "cgl_amr_primitive_mirror_slope_stress.athinput",
            "cgl_amr_gpu_mirror_slope",
        )
        mirror_slope = _user_history("cgl_amr_gpu_mirror_slope")
        _assert_clean_user(mirror_slope)
        assert mirror_slope["amr_mirr"][-1] > 0.0
        assert mirror_slope["amr_slp"][-1] > 0.0
    finally:
        _cleanup()


def test_cgl_lf_amr_primitive_churn_gpu():
    try:
        _run("cgl_lf_amr_primitive_churn.athinput", "cgl_amr_gpu_lf_churn")
        user = _user_history("cgl_amr_gpu_lf_churn")
        mhd = _mhd_history("cgl_amr_gpu_lf_churn")
        _assert_clean_user(user)
        _assert_clean_lf(mhd)
        assert np.max(user["ncell"]) > user["ncell"][0]
        assert user["ncell"][-1] < np.max(user["ncell"])
    finally:
        _cleanup()


def test_cgl_amr_3d_and_passive_gpu():
    try:
        _run(
            "cgl_amr_paper_oblique_wave_static.athinput",
            "cgl_amr_gpu_wave_static",
        )

        _run("cgl_amr_primitive_3d_current.athinput", "cgl_amr_gpu_3d")
        three_d = _user_history("cgl_amr_gpu_3d")
        _assert_clean_user(three_d, max_ndiv=1.0e-9)
        assert np.max(three_d["ncell"]) > three_d["ncell"][0]

        _run(
            "cgl_amr_passive_primitive_current.athinput",
            "cgl_amr_gpu_passive",
        )
        passive = _user_history("cgl_amr_gpu_passive")
        _assert_clean_user(passive)
        assert np.max(passive["ncell"]) > passive["ncell"][0]
    finally:
        _cleanup()
