"""Short integration coverage for shipped frame-aware problem generators."""

from pathlib import Path
from subprocess import PIPE, run
import tempfile

import athena_read
import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.parametrize(
    ("problem", "input_path", "basename", "overrides"),
    [
        (
            "cloud_crushing",
            "inputs/hydro/cloud_crushing_snr.athinput",
            "CloudCrushingSNR",
            [
                "mesh/nx1=16",
                "mesh/nx2=8",
                "mesh/nx3=8",
                "meshblock/nx1=16",
                "meshblock/nx2=8",
                "meshblock/nx3=8",
                "time/nlim=6",
                "time/tlim=0.004",
                "output1/dt=1.0e-20",
                "output2/dt=10",
                "frame_tracking/start_time=-1.0",
            ],
        ),
        (
            "TRML_frame_tracking",
            "inputs/hydro/TRML/TRML_frame_tracking.athinput",
            "TRMLFrameTracking",
            [
                "mesh/nx1=8",
                "mesh/nx2=8",
                "mesh/nx3=16",
                "meshblock/nx1=8",
                "meshblock/nx2=8",
                "meshblock/nx3=16",
                "time/nlim=6",
                "time/tlim=0.04",
                "output1/dt=1.0e-20",
                "output2/dt=10",
                "output3/dt=10",
            ],
        ),
    ],
)
def test_shipped_frame_aware_example_runs(
    problem: str, input_path: str, basename: str, overrides: list[str]
) -> None:
    with tempfile.TemporaryDirectory(prefix=f"athenak_{problem}_build_") as build_tmp:
        build_dir = Path(build_tmp) / "build"
        run_dir = Path(build_tmp) / "run"
        configure = run(
            ["cmake", "-S", str(REPO_ROOT), "-B", str(build_dir), f"-DPROBLEM={problem}"],
            check=False,
            stdout=PIPE,
            stderr=PIPE,
            text=True,
        )
        assert configure.returncode == 0, configure.stdout + configure.stderr
        compile_result = run(
            ["cmake", "--build", str(build_dir), "-j", "4"],
            check=False,
            stdout=PIPE,
            stderr=PIPE,
            text=True,
        )
        assert compile_result.returncode == 0, (
            compile_result.stdout + compile_result.stderr
        )
        run_dir.mkdir()
        result = run(
            [
                str(build_dir / "src" / "athena"),
                "-i",
                str(REPO_ROOT / input_path),
                "-d",
                str(run_dir),
                *overrides,
            ],
            check=False,
            stdout=PIPE,
            stderr=PIPE,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        history = athena_read.hst(str(run_dir / f"{basename}.frame_tracker.hst"))
        assert np.all(np.isfinite(history["ft_weight"]))
        assert np.all(np.isfinite(history["ft_skip"]))
        velocity_field = next(name for name in history if name.startswith("ft_vf_"))
        assert np.all(np.isfinite(history[velocity_field]))
        assert history["ft_weight"].size >= 2

        if problem == "TRML_frame_tracking":
            restart_input = Path(build_tmp) / "TRML_restart_check.athinput"
            restart_input.write_text(
                (REPO_ROOT / input_path).read_text().replace(
                    "variable    = hydro_w\n",
                    "variable    = hydro_w\n"
                    "data_format = %24.16e\n"
                    "slice_x1    = 0.0\n"
                    "slice_x2    = 0.0\n",
                    1,
                )
            )
            restart_overrides = [
                "mesh/nx1=8",
                "mesh/nx2=8",
                "mesh/nx3=16",
                "meshblock/nx1=8",
                "meshblock/nx2=8",
                "meshblock/nx3=16",
                "time/tlim=1.0",
                "output1/dt=1.0e-20",
                "output2/file_type=tab",
                "output2/variable=hydro_u",
                "output2/dt=1.0e-20",
                "output3/dt=1.0e-20",
                "frame_tracking/diagnostic_every=1000",
            ]
            continuous_dir = Path(build_tmp) / "continuous"
            split_dir = Path(build_tmp) / "split"
            continuous_dir.mkdir()
            split_dir.mkdir()
            for directory, basename, nlim in (
                (continuous_dir, "TRMLContinuous", "10"),
                (split_dir, "TRMLSplit", "5"),
            ):
                integration = run(
                    [
                        str(build_dir / "src" / "athena"),
                        "-i",
                        str(restart_input),
                        "-d",
                        str(directory),
                        f"job/basename={basename}",
                        f"time/nlim={nlim}",
                        *restart_overrides,
                    ],
                    check=False,
                    stdout=PIPE,
                    stderr=PIPE,
                    text=True,
                )
                assert integration.returncode == 0, (
                    integration.stdout + integration.stderr
                )
            restart_files = sorted((split_dir / "rst").glob("TRMLSplit.*.rst"))
            assert restart_files
            restarted = run(
                [
                    str(build_dir / "src" / "athena"),
                    "-r",
                    str(restart_files[-1]),
                    "-i",
                    str(restart_input),
                    "-d",
                    str(split_dir),
                    "job/basename=TRMLSplit",
                    "time/nlim=10",
                    *restart_overrides,
                ],
                check=False,
                stdout=PIPE,
                stderr=PIPE,
                text=True,
            )
            assert restarted.returncode == 0, restarted.stdout + restarted.stderr
            continuous_history = athena_read.hst(
                str(continuous_dir / "TRMLContinuous.frame_tracker.hst")
            )
            split_history = athena_read.hst(
                str(split_dir / "TRMLSplit.frame_tracker.hst")
            )
            tolerance = 100.0 * np.finfo(float).eps
            for field in (
                "ft_vf_x3",
                "ft_dx_x3",
                "ft_pos_x3",
                "ft_err_x3",
                "ft_dv_x3",
            ):
                np.testing.assert_allclose(
                    split_history[field][-1],
                    continuous_history[field][-1],
                    rtol=tolerance,
                    atol=tolerance,
                )
            continuous_tabs = sorted((continuous_dir / "tab").glob("*.tab"))
            split_tabs = sorted((split_dir / "tab").glob("*.tab"))
            assert continuous_tabs and split_tabs
            np.testing.assert_allclose(
                np.loadtxt(split_tabs[-1]),
                np.loadtxt(continuous_tabs[-1]),
                rtol=tolerance,
                atol=tolerance,
            )
