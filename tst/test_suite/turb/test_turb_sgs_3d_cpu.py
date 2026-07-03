"""Regression test for 3D coarsened SGS output."""

from pathlib import Path
import subprocess
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUT = REPO_ROOT / "tst" / "inputs" / "turb_sgs_3d.athinput"
ATHENA = Path.cwd() / "athena"
COARSEN_FACTOR = 2
SGS_NAMES = [
    "dens",
    "velx",
    "vely",
    "velz",
    "tau_xx",
    "tau_xy",
    "tau_xz",
    "tau_yy",
    "tau_yz",
    "tau_zz",
]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))

from bin_convert import read_binary, read_coarsened_binary  # noqa: E402


def run_athena(output_dir):
    """Run the focused input in an isolated output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    return subprocess.run(
        [str(ATHENA), "-d", str(output_dir), "-i", str(INPUT)],
        capture_output=True,
        text=True,
        check=False,
    )


def require_success(result):
    """Report executable output when a regression run fails."""
    assert result.returncode == 0, result.stdout + result.stderr


def latest(path, pattern):
    """Return the lexically latest numbered output file."""
    outputs = sorted(path.glob(pattern))
    assert outputs
    return outputs[-1]


def cube_mean(values, factor):
    """Average non-overlapping boxes in all active dimensions."""
    nz, ny, nx = values.shape
    return values.reshape(
        nz // factor, factor, ny // factor, factor, nx // factor, factor
    ).mean(axis=(1, 3, 5))


def expected_sgs_fields(state, block):
    """Compute the expected coarsened Favre state directly from hydro_u."""
    rho = np.asarray(state["mb_data"]["dens"][block])
    mx = np.asarray(state["mb_data"]["mom1"][block])
    my = np.asarray(state["mb_data"]["mom2"][block])
    mz = np.asarray(state["mb_data"]["mom3"][block])

    rho_bar = cube_mean(rho, COARSEN_FACTOR)
    mx_bar = cube_mean(mx, COARSEN_FACTOR)
    my_bar = cube_mean(my, COARSEN_FACTOR)
    mz_bar = cube_mean(mz, COARSEN_FACTOR)

    return {
        "dens": rho_bar,
        "velx": mx_bar / rho_bar,
        "vely": my_bar / rho_bar,
        "velz": mz_bar / rho_bar,
        "tau_xx": cube_mean(mx * mx / rho, COARSEN_FACTOR)
        - mx_bar * mx_bar / rho_bar,
        "tau_xy": cube_mean(mx * my / rho, COARSEN_FACTOR)
        - mx_bar * my_bar / rho_bar,
        "tau_xz": cube_mean(mx * mz / rho, COARSEN_FACTOR)
        - mx_bar * mz_bar / rho_bar,
        "tau_yy": cube_mean(my * my / rho, COARSEN_FACTOR)
        - my_bar * my_bar / rho_bar,
        "tau_yz": cube_mean(my * mz / rho, COARSEN_FACTOR)
        - my_bar * mz_bar / rho_bar,
        "tau_zz": cube_mean(mz * mz / rho, COARSEN_FACTOR)
        - mz_bar * mz_bar / rho_bar,
    }


def test_3d_sgs_output_matches_direct_favre_filter(tmp_path):
    """The producer writes final Favre velocities and SGS stresses on a 3D mesh."""
    run_dir = tmp_path / "run"
    require_success(run_athena(run_dir))

    state = read_binary(str(latest(run_dir / "bin", "*.state.*.bin")))
    sgs = read_coarsened_binary(
        str(latest(run_dir / "cbin_sgs_2", "*.sgs.*.cbin"))
    )

    assert sgs["var_names"] == SGS_NAMES
    assert state["n_mbs"] == 1
    assert sgs["n_mbs"] == state["n_mbs"]
    assert sgs["Nx1"] == state["Nx1"] // COARSEN_FACTOR
    assert sgs["Nx2"] == state["Nx2"] // COARSEN_FACTOR
    assert sgs["Nx3"] == state["Nx3"] // COARSEN_FACTOR
    assert sgs["nx1_mb"] == state["nx1_mb"] // COARSEN_FACTOR
    assert sgs["nx2_mb"] == state["nx2_mb"] // COARSEN_FACTOR
    assert sgs["nx3_mb"] == state["nx3_mb"] // COARSEN_FACTOR
    np.testing.assert_array_equal(state["mb_logical"], sgs["mb_logical"])

    for block in range(state["n_mbs"]):
        expected = expected_sgs_fields(state, block)
        for name in SGS_NAMES:
            np.testing.assert_allclose(
                sgs["mb_data"][name][block],
                expected[name],
                rtol=5.0e-6,
                atol=5.0e-8,
            )
