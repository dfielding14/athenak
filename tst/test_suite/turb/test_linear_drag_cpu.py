"""Regression tests for uniform linear Rayleigh drag."""

from pathlib import Path
import subprocess
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUT = REPO_ROOT / "tst" / "inputs" / "linear_drag.athinput"
ATHENA = Path.cwd() / "athena"
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))

from bin_convert import read_binary  # noqa: E402


def run_athena(output_dir, *overrides):
    """Run the focused input in an isolated output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [str(ATHENA), "-d", str(output_dir), "-i", str(INPUT), *overrides]
    return subprocess.run(command, capture_output=True, text=True, check=False)


def test_linear_drag_matches_uniform_flow_decay(tmp_path):
    """Uniform velocity follows dv/dt=-alpha*v without changing density."""
    run_dir = tmp_path / "run"
    result = run_athena(run_dir)
    assert result.returncode == 0, result.stdout + result.stderr

    outputs = sorted((run_dir / "bin").glob("*.prim.*.bin"))
    assert len(outputs) == 2
    initial = read_binary(str(outputs[0]))
    final = read_binary(str(outputs[-1]))
    decay = np.exp(-2.0 * final["time"])

    np.testing.assert_allclose(
        final["mb_data"]["dens"], initial["mb_data"]["dens"], rtol=0.0, atol=1.0e-14
    )
    np.testing.assert_allclose(
        final["mb_data"]["eint"], initial["mb_data"]["eint"], rtol=3.0e-6, atol=0.0
    )
    for name in ("velx", "vely", "velz"):
        np.testing.assert_allclose(
            np.asarray(final["mb_data"][name]),
            decay * np.asarray(initial["mb_data"][name]),
            rtol=2.0e-5,
            atol=1.0e-12,
        )


def test_negative_drag_rate_is_rejected(tmp_path):
    """Rayleigh drag cannot amplify velocity through a negative damping rate."""
    result = run_athena(tmp_path / "negative", "hydro_srcterms/drag_rate=-1.0")
    assert result.returncode != 0
    assert "drag_rate must not be negative" in result.stdout + result.stderr
