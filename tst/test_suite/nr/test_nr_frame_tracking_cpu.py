"""Regression coverage for serial non-relativistic frame tracking."""

from pathlib import Path
import shutil
from subprocess import PIPE, run

import athena_read
import numpy as np
import test_suite.testutils as testutils


HYDRO_INPUT = "inputs/frame_tracking_hydro.athinput"
MHD_INPUT = "inputs/frame_tracking_mhd.athinput"
AMBIGUOUS_INPUT = "inputs/frame_tracking_invalid_ambiguous.athinput"
LEGACY_INPUT = "inputs/frame_tracking_legacy.athinput"


def remove_outputs(basename: str) -> None:
    for path in Path(".").glob(f"{basename}.*.hst"):
        path.unlink()
    for path in Path("rst").glob(f"{basename}.*.rst") if Path("rst").exists() else []:
        path.unlink()


def test_hydro_structured_history() -> None:
    basename = "FrameTrackHydro"
    try:
        assert testutils.run(HYDRO_INPUT)
        data = athena_read.hst(f"{basename}.frame_tracker.hst")
        required = {
            "ft_vf_x1",
            "ft_dx_x1",
            "ft_pos_x1",
            "ft_err_x1",
            "ft_dv_x1",
            "ft_vf_x2",
            "ft_weight",
            "ft_misses",
            "ft_limit",
            "ft_skip",
        }
        assert required.issubset(data)
        assert np.all(np.isfinite(data["ft_vf_x1"]))
        assert np.any(np.abs(data["ft_vf_x1"]) > 0.0)
        assert np.max(data["ft_misses"]) == 0.0
    finally:
        remove_outputs(basename)
        shutil.rmtree("rst", ignore_errors=True)
        testutils.cleanup()


def test_mhd_structured_history() -> None:
    basename = "FrameTrackMHD"
    try:
        assert testutils.run(MHD_INPUT)
        data = athena_read.hst(f"{basename}.frame_tracker.hst")
        assert np.all(np.isfinite(data["ft_vf_x1"]))
        assert np.any(np.abs(data["ft_vf_x1"]) > 0.0)
        assert np.max(data["ft_misses"]) == 0.0
    finally:
        remove_outputs(basename)
        testutils.cleanup()


def test_relativistic_configuration_is_rejected() -> None:
    result = run(
        ["./athena", "-i", HYDRO_INPUT, "coord/special_rel=true"],
        check=False,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "supports non-relativistic hydro/MHD only" in output


def test_ambiguous_dual_fluid_configuration_is_rejected() -> None:
    result = run(
        ["./athena", "-i", AMBIGUOUS_INPUT],
        check=False,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "set tracked_fluid=hydro or tracked_fluid=mhd explicitly" in output


def test_legacy_frame_state_warns_and_remains_finite() -> None:
    basename = "FrameTrackLegacy"
    try:
        result = run(
            [
                "./athena",
                "-i",
                LEGACY_INPUT,
                f"job/basename={basename}",
            ],
            check=False,
            stdout=PIPE,
            stderr=PIPE,
            text=True,
        )
        assert result.returncode == 0
        assert "restoring legacy frame state without complete controller history" in (
            result.stdout + result.stderr
        )
        data = athena_read.hst(f"{basename}.frame_tracker.hst")
        assert np.all(np.isfinite(data["ft_vf_x1"]))
    finally:
        remove_outputs(basename)
        shutil.rmtree("rst", ignore_errors=True)
        testutils.cleanup()
