"""Restart-equivalence regression for non-relativistic frame tracking."""

from pathlib import Path
import shutil

import athena_read
import numpy as np
import test_suite.testutils as testutils


INPUT = "inputs/frame_tracking_hydro.athinput"
THREE_D_INPUT = "inputs/frame_tracking_3d.athinput"
EPS_TOL = 100.0 * np.finfo(float).eps


def clean_case(basename: str) -> None:
    for path in Path(".").glob(f"{basename}.*.hst"):
        path.unlink()
    for path in Path("rst").glob(f"{basename}.*.rst") if Path("rst").exists() else []:
        path.unlink()
    for path in Path("tab").glob(f"{basename}.*.tab") if Path("tab").exists() else []:
        path.unlink()


def test_restart_matches_continuous_controller_and_state() -> None:
    continuous = "FrameTrackContinuous"
    split = "FrameTrackSplit"
    try:
        assert testutils.run(INPUT, [f"job/basename={continuous}"])
        assert testutils.run(INPUT, [f"job/basename={split}", "time/nlim=3"])
        restart_files = sorted(Path("rst").glob(f"{split}.*.rst"))
        assert restart_files
        assert testutils.run(
            INPUT,
            ["-r", str(restart_files[-1]), f"job/basename={split}", "time/nlim=6"],
        )

        continuous_hist = athena_read.hst(f"{continuous}.frame_tracker.hst")
        split_hist = athena_read.hst(f"{split}.frame_tracker.hst")
        for field in (
            "ft_vf_x1",
            "ft_dx_x1",
            "ft_pos_x1",
            "ft_err_x1",
            "ft_dv_x1",
            "ft_weight",
            "ft_misses",
            "ft_recov",
            "ft_limit",
            "ft_skip",
        ):
            np.testing.assert_allclose(
                split_hist[field][-1],
                continuous_hist[field][-1],
                rtol=EPS_TOL,
                atol=EPS_TOL,
            )

        continuous_tabs = sorted(Path("tab").glob(f"{continuous}.hydro_u.*.tab"))
        split_tabs = sorted(Path("tab").glob(f"{split}.hydro_u.*.tab"))
        assert continuous_tabs and split_tabs
        continuous_state = np.loadtxt(continuous_tabs[-1])
        split_state = np.loadtxt(split_tabs[-1])
        np.testing.assert_allclose(
            split_state,
            continuous_state,
            rtol=EPS_TOL,
            atol=EPS_TOL,
        )
    finally:
        clean_case(continuous)
        clean_case(split)
        shutil.rmtree("rst", ignore_errors=True)
        testutils.cleanup()


def test_three_axis_restart_matches_continuous_controller_and_state() -> None:
    continuous = "FrameTrack3DContinuous"
    split = "FrameTrack3DSplit"
    try:
        assert testutils.run(THREE_D_INPUT, [f"job/basename={continuous}"])
        assert testutils.run(THREE_D_INPUT, [f"job/basename={split}", "time/nlim=3"])
        restart_files = sorted(Path("rst").glob(f"{split}.*.rst"))
        assert restart_files
        assert testutils.run(
            THREE_D_INPUT,
            ["-r", str(restart_files[-1]), f"job/basename={split}", "time/nlim=6"],
        )
        continuous_hist = athena_read.hst(f"{continuous}.frame_tracker.hst")
        split_hist = athena_read.hst(f"{split}.frame_tracker.hst")
        fields = ["ft_weight", "ft_misses", "ft_recov", "ft_limit", "ft_skip"]
        for axis in (1, 2, 3):
            fields.extend(
                [
                    f"ft_vf_x{axis}",
                    f"ft_dx_x{axis}",
                    f"ft_pos_x{axis}",
                    f"ft_err_x{axis}",
                    f"ft_dv_x{axis}",
                ]
            )
        for field in fields:
            np.testing.assert_allclose(
                split_hist[field][-1],
                continuous_hist[field][-1],
                rtol=EPS_TOL,
                atol=EPS_TOL,
            )
        continuous_tabs = sorted(Path("tab").glob(f"{continuous}.hydro_u.*.tab"))
        split_tabs = sorted(Path("tab").glob(f"{split}.hydro_u.*.tab"))
        assert continuous_tabs and split_tabs
        np.testing.assert_allclose(
            np.loadtxt(split_tabs[-1]),
            np.loadtxt(continuous_tabs[-1]),
            rtol=EPS_TOL,
            atol=EPS_TOL,
        )
    finally:
        clean_case(continuous)
        clean_case(split)
        shutil.rmtree("rst", ignore_errors=True)
        testutils.cleanup()
