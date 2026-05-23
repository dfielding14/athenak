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
    for path in Path(".").glob("cgl_*.hst"):
        path.unlink()
    shutil.rmtree("bin", ignore_errors=True)


def _tab(basename):
    return testutils.athena_read.tab(f"tab/{basename}.mhd_w.00001.tab")


def _assert_clean_lf_history(history):
    assert history["lf_nstage"][-1] > 0.0
    assert history["lf_dfloor"][-1] == 0.0
    assert history["lf_pfloor"][-1] == 0.0
    assert history["lf_nonfin"][-1] == 0.0
    assert history["lf_nonpos"][-1] == 0.0
    assert history["lf_hardbd"][-1] == 0.0


def test_cgl_lf_quantitative_decay_and_diagnostics():
    try:
        _run("cgl_lf_decay.athinput", "cgl_ci_decay")
        history = testutils.athena_read.hst("cgl_ci_decay.mhd.hst")
        _assert_clean_lf_history(history)
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


def test_cgl_lf_explicit_reference_agrees_with_capped_sts():
    try:
        _run(
            "cgl_lf_decay.athinput",
            "cgl_ci_sts_reference",
            "time/sts_max_dt_ratio=1.0",
            "time/nlim=-1",
        )
        _run(
            "cgl_lf_decay.athinput",
            "cgl_ci_explicit_reference",
            "mhd/cgl_heat_flux_integrator=explicit",
            "time/sts_integrator=none",
            "time/nlim=-1",
        )
        sts = _tab("cgl_ci_sts_reference")
        explicit = _tab("cgl_ci_explicit_reference")
        fields = set(sts).intersection(explicit) - {"time", "cycle"}
        max_difference = max(np.max(np.abs(sts[field] - explicit[field]))
                             for field in fields)
        assert max_difference < 1.0e-3
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_sts_reference.mhd.hst")
        )
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_explicit_reference.mhd.hst")
        )
    finally:
        _cleanup()


def test_cgl_lf_explicit_reference_finite_collision_split():
    try:
        common = ("time/nlim=-1", "mhd/nu_coll=1.0")
        _run(
            "cgl_lf_decay.athinput",
            "cgl_ci_sts_collision",
            "time/sts_max_dt_ratio=1.0",
            *common,
        )
        _run(
            "cgl_lf_decay.athinput",
            "cgl_ci_explicit_collision",
            "mhd/cgl_heat_flux_integrator=explicit",
            "time/sts_integrator=none",
            *common,
        )
        sts = _tab("cgl_ci_sts_collision")
        explicit = _tab("cgl_ci_explicit_collision")
        assert np.max(np.abs(sts["eint"] - explicit["eint"])) < 1.0e-3
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_sts_collision.mhd.hst")
        )
        _assert_clean_lf_history(
            testutils.athena_read.hst("cgl_ci_explicit_collision.mhd.hst")
        )
    finally:
        _cleanup()


def test_cgl_fofc_live_flux_mutation():
    try:
        _run("cgl_fofc.athinput", "cgl_ci_fofc")
    finally:
        _cleanup()


def test_cgl_lf_amr_conserved_prolongation_stays_admissible():
    try:
        _run("cgl_lf_amr_2d.athinput", "cgl_ci_amr")
        user = testutils.athena_read.hst("cgl_ci_amr.user.hst")
        mhd = testutils.athena_read.hst("cgl_ci_amr.mhd.hst")
        assert user["ncell"][-1] > user["ncell"][0]
        assert np.max(user["max_ndiv"]) < 1.0e-12
        assert np.max(user["bad_state"]) == 0.0
        assert np.all(np.isfinite(user["abs_anis"]))
        _assert_clean_lf_history(mhd)
        energy_scale = max(abs(mhd["tot-E"][0]), 1.0e-30)
        energy_residual = abs(mhd["tot-E"][-1] - mhd["tot-E"][0]) / energy_scale
        assert energy_residual < 5.0e-3
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


def test_cgl_lf_explicit_reference_rejects_sts_configuration():
    command = [
        "./athena",
        "-i",
        f"{INPUT_ROOT}/cgl_lf_decay.athinput",
        "mhd/cgl_heat_flux_integrator=explicit",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "explicit reference integration requires" in result.stdout


def test_cgl_lf_amr_rejects_primitive_prolongation():
    command = [
        "./athena",
        "-i",
        f"{INPUT_ROOT}/cgl_lf_amr_2d.athinput",
        "mesh_refinement/prolong_primitives=true",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "use conserved prolongation for LF AMR runs" in result.stdout
