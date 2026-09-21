"""Regression tests for uniform linear Rayleigh drag."""

from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUT = REPO_ROOT / "tst" / "inputs" / "linear_drag.athinput"
ATHENA = Path.cwd() / "athena"
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))

from bin_convert import read_binary  # noqa: E402


def run_athena(output_dir, *overrides):
    """Run the focused input in an isolated output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    # Input blocks permit new parameters; command-line overrides require existing ones.
    input_text = INPUT.read_text()
    for override in overrides:
        block, assignment = override.split("/", 1)
        input_text += f"\n<{block}>\n{assignment}\n"
    input_path = output_dir / "drag.athinput"
    input_path.write_text(input_text)
    command = [str(ATHENA), "-d", str(output_dir), "-i", str(input_path)]
    return subprocess.run(command, capture_output=True, text=True, check=False, timeout=30)


@pytest.mark.parametrize("integrator", ["rk2", "rk3", "rk4"])
@pytest.mark.parametrize("eos", ["ideal", "isothermal"])
def test_linear_drag_matches_uniform_flow_decay(tmp_path, integrator, eos):
    """Uniform velocity follows dv/dt=-alpha*v without changing density."""
    run_dir = tmp_path / "run"
    overrides = [
        "mesh/nx2=8", "meshblock/nx2=4", "time/cfl_number=0.05",
        f"time/integrator={integrator}", f"hydro/eos={eos}",
    ]
    if eos == "isothermal":
        overrides.extend([
            "hydro/iso_sound_speed=1.0", "hydro/rsolver=llf",
            "problem/pgen_name=advection", "problem/flow_dir=2",
            "problem/iproblem=1", "problem/amplitude=0.0",
        ])
    result = run_athena(run_dir, *overrides)
    assert result.returncode == 0, result.stdout + result.stderr

    outputs = sorted((run_dir / "bin").glob("*.prim.*.bin"))
    assert len(outputs) == 2
    initial = read_binary(str(outputs[0]))
    final = read_binary(str(outputs[-1]))
    decay = np.exp(-2.0 * final["time"])

    np.testing.assert_allclose(
        final["mb_data"]["dens"], initial["mb_data"]["dens"], rtol=0.0, atol=1.0e-14
    )
    if eos == "ideal":
        np.testing.assert_allclose(
            final["mb_data"]["eint"], initial["mb_data"]["eint"],
            rtol=2.0e-5, atol=0.0,
        )
    for name in ("velx", "vely", "velz"):
        np.testing.assert_allclose(
            np.asarray(final["mb_data"][name]),
            decay * np.asarray(initial["mb_data"][name]),
            rtol=2.0e-5,
            atol=1.0e-12,
        )


@pytest.mark.parametrize("rate", ["-1.0", "nan", "inf", "-inf"])
def test_invalid_drag_rate_is_rejected(tmp_path, rate):
    """Invalid damping rates must fail before the first timestep."""
    result = run_athena(tmp_path / "invalid", f"hydro_srcterms/drag_rate={rate}")
    assert result.returncode != 0
    assert "drag_rate must be finite and nonnegative" in result.stdout + result.stderr


@pytest.mark.parametrize("cooling", [False, True])
def test_drag_limits_timestep_with_cooling(tmp_path, cooling):
    """A cooling reduction must retain the stronger drag timestep limit."""
    run_dir = tmp_path / "run"
    result = run_athena(
        run_dir, "time/nlim=1", "hydro_srcterms/drag_rate=100.0",
        f"hydro_srcterms/ism_cooling={str(cooling).lower()}",
        "hydro_srcterms/hrate=0.0", "units/length_cgs=1.0",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    final = read_binary(str(sorted((run_dir / "bin").glob("*.prim.*.bin"))[-1]))
    assert final["time"] == pytest.approx(0.4 / 100.0)
