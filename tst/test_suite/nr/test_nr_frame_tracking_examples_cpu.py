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
