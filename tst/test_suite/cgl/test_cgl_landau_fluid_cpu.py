"""CPU regressions for the CGL Landau-fluid closure and CGL FOFC path."""

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


def _cleanup():
    testutils.cleanup()
    for path in Path(".").glob("cgl_*.mhd.hst"):
        path.unlink()
    shutil.rmtree("bin", ignore_errors=True)


def test_cgl_lf_quantitative_decay_and_diagnostics():
    try:
        _run("cgl_lf_decay.athinput", "cgl_ci_decay")
        history = testutils.athena_read.hst("cgl_ci_decay.mhd.hst")
        assert history["lf_nstage"][-1] > 0.0
        assert history["lf_dfloor"][-1] == 0.0
        assert history["lf_pfloor"][-1] == 0.0
        assert history["lf_nonfin"][-1] == 0.0
        assert history["lf_nonpos"][-1] == 0.0
        assert history["lf_hardbd"][-1] == 0.0
    finally:
        _cleanup()


def test_cgl_lf_limiter_occupancy_remains_admissible():
    try:
        _run("cgl_lf_limiter.athinput", "cgl_ci_limiter")
        history = testutils.athena_read.hst("cgl_ci_limiter.mhd.hst")
        assert history["lf_mirror"][-1] > 0.0
        assert history["lf_nonfin"][-1] == 0.0
        assert history["lf_nonpos"][-1] == 0.0
        assert history["lf_hardbd"][-1] == 0.0
        assert np.all(np.diff(history["lf_nstage"]) >= 0.0)
    finally:
        _cleanup()


def test_cgl_fofc_live_flux_mutation():
    try:
        _run("cgl_fofc.athinput", "cgl_ci_fofc")
    finally:
        _cleanup()


def test_cgl_lf_invalid_integrator_is_rejected():
    command = [
        "./athena",
        "-i",
        f"{INPUT_ROOT}/cgl_lf_decay.athinput",
        "mhd/cgl_heat_flux_integrator=invalid",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "cgl_heat_flux_integrator" in result.stdout
