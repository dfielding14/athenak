"""Regression tests for fully 2D turbulence and coarsened SGS output."""

from pathlib import Path
import subprocess
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUT = REPO_ROOT / "tst" / "inputs" / "turb_sgs_2d.athinput"
ATHENA = Path.cwd() / "athena"
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))

from bin_convert import read_binary, read_coarsened_binary  # noqa: E402


def run_athena(output_dir, *overrides):
    """Run the focused input in an isolated output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [str(ATHENA), "-d", str(output_dir), "-i", str(INPUT), *overrides]
    return subprocess.run(command, capture_output=True, text=True, check=False)


def require_success(result):
    """Report executable output when a regression run fails."""
    assert result.returncode == 0, result.stdout + result.stderr


def latest(path, pattern):
    """Return the lexically latest numbered output file."""
    outputs = sorted(path.glob(pattern))
    assert outputs
    return outputs[-1]


def square_mean(values, factor):
    """Average square, non-overlapping boxes in the active x-y plane."""
    nz, ny, nx = values.shape
    return values.reshape(nz, ny // factor, factor, nx // factor, factor).mean(
        axis=(2, 4)
    )


def test_2d_sgs_output_matches_direct_favre_filter(tmp_path):
    """The producer writes final Favre velocities and SGS stresses on a 2D mesh."""
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir))

    state = read_binary(str(latest(run_dir / "bin", "*.state.*.bin")))
    force = read_binary(str(latest(run_dir / "bin", "*.force.*.bin")))
    sgs = read_coarsened_binary(
        str(latest(run_dir / "cbin_sgs_2", "*.sgs.*.cbin"))
    )

    assert sgs["var_names"] == ["dens", "velx", "vely", "tau_xx", "tau_xy", "tau_yy"]
    assert sgs["Nx3"] == 1
    assert sgs["nx3_mb"] == 1
    np.testing.assert_array_equal(state["mb_logical"], sgs["mb_logical"])

    for block in range(state["n_mbs"]):
        rho = np.asarray(state["mb_data"]["dens"][block])
        mx = np.asarray(state["mb_data"]["mom1"][block])
        my = np.asarray(state["mb_data"]["mom2"][block])
        rho_bar = square_mean(rho, 2)
        mx_bar = square_mean(mx, 2)
        my_bar = square_mean(my, 2)
        expected = {
            "dens": rho_bar,
            "velx": mx_bar / rho_bar,
            "vely": my_bar / rho_bar,
            "tau_xx": square_mean(mx * mx / rho, 2) - mx_bar * mx_bar / rho_bar,
            "tau_xy": square_mean(mx * my / rho, 2) - mx_bar * my_bar / rho_bar,
            "tau_yy": square_mean(my * my / rho, 2) - my_bar * my_bar / rho_bar,
        }
        for name, values in expected.items():
            np.testing.assert_allclose(
                sgs["mb_data"][name][block], values, rtol=5.0e-6, atol=5.0e-8
            )

    assert np.all(np.asarray(force["mb_data"]["force3"]) == 0.0)
    assert np.all(np.asarray(state["mb_data"]["mom3"]) == 0.0)
    history = np.loadtxt(run_dir / "turb_sgs_2d_test.hydro.hst")
    assert np.all(history[:, 5] == 0.0)
    assert np.all(history[:, 8] == 0.0)


def test_2d_turbulence_rejects_kz_modes(tmp_path):
    """A nominally 2D mesh cannot silently enable out-of-plane forcing modes."""
    result = run_athena(tmp_path / "bad_kz", "turb_driving/max_kz=1")
    assert result.returncode != 0
    assert "two-dimensional meshes require min_kz = max_kz = 0" in (
        result.stdout + result.stderr
    )
