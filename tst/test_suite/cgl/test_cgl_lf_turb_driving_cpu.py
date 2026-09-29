"""Interaction regression for CGL Landau fluid with modal turbulence driving."""

from pathlib import Path
import shutil

import numpy as np

import test_suite.testutils as testutils
from test_suite.turb.test_turb_driving_cpu import assert_equal_force_outputs
from test_suite.turb.test_turb_driving_cpu import read_force_blocks, rms_acceleration


INPUT_FILE = "../../../inputs/tests/cgl_lf_turb_driving_amr.athinput"


def _run(basename, *flags):
    testutils.run(INPUT_FILE, [f"job/basename={basename}", *flags])


def _latest_force(basename):
    outputs = sorted(Path("bin").glob(f"{basename}.force.*.bin"))
    assert outputs, f"no rendered turbulence output found for {basename}"
    return outputs[-1]


def _assert_clean_lf_history(history):
    assert history["lf_nstage"][-1] > 0.0
    for column in (
        "lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos", "lf_hardbd",
    ):
        assert history[column][-1] == 0.0


def test_cgl_lf_amr_with_turbulence_driving_restarts_reproducibly():
    """Forced LF evolution refines, stays admissible, and restores modal state."""
    try:
        _run("cgl_ci_turb_reference")
        reference_mhd = testutils.athena_read.hst("cgl_ci_turb_reference.mhd.hst")
        reference_user = testutils.athena_read.hst("cgl_ci_turb_reference.user.hst")
        reference_force = _latest_force("cgl_ci_turb_reference")
        assert reference_user["ncell"][-1] > reference_user["ncell"][0]
        assert np.max(reference_user["max_ndiv"]) < 1.0e-12
        assert np.max(reference_user["bad_state"]) == 0.0
        _assert_clean_lf_history(reference_mhd)
        assert np.isclose(
            rms_acceleration(read_force_blocks(reference_force)),
            0.02,
            rtol=2.0e-6,
        )

        _run("cgl_ci_turb_partial", "time/nlim=1")
        restart_paths = sorted(Path("rst").glob("cgl_ci_turb_partial.*.rst"))
        assert restart_paths, "forced partial run did not write a restart checkpoint"
        assert testutils.run_command([
            "./athena",
            "-r",
            str(restart_paths[-1]),
            "job/basename=cgl_ci_turb_resumed",
            "time/nlim=3",
        ])

        resumed_mhd = testutils.athena_read.hst("cgl_ci_turb_resumed.mhd.hst")
        resumed_user = testutils.athena_read.hst("cgl_ci_turb_resumed.user.hst")
        _assert_clean_lf_history(resumed_mhd)
        assert np.max(resumed_user["bad_state"]) == 0.0
        for column in ("mass", "tot-E", "aam-D"):
            assert np.isclose(
                reference_mhd[column][-1],
                resumed_mhd[column][-1],
                rtol=1.0e-12,
                atol=1.0e-12,
            )
        assert_equal_force_outputs(reference_force, _latest_force("cgl_ci_turb_resumed"))
    finally:
        for path in Path(".").glob("cgl_ci_turb_*.hst"):
            path.unlink()
        shutil.rmtree("bin", ignore_errors=True)
        shutil.rmtree("rst", ignore_errors=True)
        testutils.cleanup()
