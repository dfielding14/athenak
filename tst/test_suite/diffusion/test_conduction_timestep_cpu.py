"""Focused CPU tests for the thermal-conduction stability bound."""

import re
import subprocess

import numpy as np

import test_suite.testutils as testutils


INPUT = "../../../inputs/tests/conduction_timestep.athinput"
GAMMA_MINUS_ONE = 2.0 / 3.0
CFL = 0.4


def _run(*flags):
    command = ["./athena", "-i", INPUT, *flags]
    return subprocess.run(command, capture_output=True, text=True, check=False)


def _initial_dt(result):
    assert result.returncode == 0, result.stdout + result.stderr
    match = re.search(r"cycle=0 .*dt=([0-9.eE+-]+)", result.stdout)
    assert match is not None, result.stdout
    return float(match.group(1))


def test_power_law_bound_uses_neighboring_face_coefficients():
    """A sharp coefficient jump uses the same face averages as the flux operator."""
    try:
        result = _run("job/basename=conduction_face_bound")

        nx1 = 32
        dx1 = 1.0 / nx1
        x1 = (np.arange(nx1) + 0.5) * dx1
        temperature = 1.0 + 1.5 * (
            1.0 - np.tanh((x1 - 0.02) / 0.001)
        )
        density = 1.0 / temperature
        kappa = np.clip(temperature**2, 0.01, 100.0)
        kappa_with_ghosts = np.pad(kappa, 1, mode="edge")
        kappa_left = 0.5 * (
            kappa_with_ghosts[:-2] + kappa_with_ghosts[1:-1]
        )
        kappa_right = 0.5 * (
            kappa_with_ghosts[1:-1] + kappa_with_ghosts[2:]
        )
        rate = GAMMA_MINUS_ONE * (kappa_left + kappa_right) / (
            density * dx1**2
        )
        expected = CFL / np.max(rate)

        np.testing.assert_allclose(_initial_dt(result), expected, rtol=2.0e-6)
    finally:
        testutils.cleanup()


def test_constant_bound_sums_anisotropic_directions():
    """An anisotropic grid uses the multidimensional row sum, not a scalar factor."""
    try:
        result = _run(
            "job/basename=conduction_anisotropic_bound",
            "mesh/nx2=8",
            "meshblock/nx2=8",
            "hydro/conductivity_model=constant",
            "hydro/conductivity=10.0",
            "problem/thot=1.0",
        )

        dx1 = 1.0 / 32.0
        dx2 = 1.0 / 8.0
        rate = GAMMA_MINUS_ONE * 2.0 * 10.0 * (
            1.0 / dx1**2 + 1.0 / dx2**2
        )
        expected = CFL / rate

        np.testing.assert_allclose(_initial_dt(result), expected, rtol=2.0e-6)
    finally:
        testutils.cleanup()
