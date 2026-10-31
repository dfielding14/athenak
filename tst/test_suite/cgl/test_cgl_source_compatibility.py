"""CGL rejects scalar-temperature sources before their kernels can read pressures."""

import os
from pathlib import Path
import subprocess

import numpy as np
import pytest

from test_suite.cgl import passive_acceptance


def run(directory, settings, mode="active"):
    directory.mkdir()
    source = passive_acceptance.input_text(
        passive=mode != "active", active_cgl=mode == "active",
        adiabatic=mode == "adiabatic", nx=16, nlim=2, amp=0, thermal=0,
    )
    deck = directory / "input.athinput"
    deck.write_text(source + settings)
    executable = str(Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve())
    result = subprocess.run([executable, "-i", str(deck), "-d", str(directory)],
                            capture_output=True, text=True, timeout=60, check=False)
    (directory / "run.log").write_text(result.stdout + result.stderr)
    return result


@pytest.mark.parametrize("mode", ["active", "isothermal", "adiabatic"])
@pytest.mark.parametrize("settings, message", [
    ("\n<mhd>\nviscosity = 0.1\ntdep_viscosity = true\nviscosity_exponent = 1\n",
     "CGL does not support temperature-dependent <mhd>/tdep_viscosity"),
    ("\n<mhd>\nviscosity = 0\nviscosity_floor = 0.1\ntdep_viscosity = true\n",
     "CGL does not support temperature-dependent <mhd>/tdep_viscosity"),
    ("\n<mhd_srcterms>\nism_cooling = true\nhrate = 0\n",
     "CGL does not support <mhd_srcterms>/ism_cooling or rel_cooling"),
    ("\n<mhd_srcterms>\nrel_cooling = true\ncrate_rel = 0\n",
     "CGL does not support <mhd_srcterms>/ism_cooling or rel_cooling"),
])
def test_cgl_rejects_scalar_temperature_sources(tmp_path, mode, settings, message):
    result = run(tmp_path / "run", settings, mode)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr


@pytest.mark.parametrize("settings", [
    "\n<mhd>\nviscosity = 0.1\ntdep_viscosity = false\n",
    "\n<mhd>\nviscosity = 0\ntdep_viscosity = true\n",
    "\n<mhd>\ntdep_viscosity = true\n",
    "\n<mhd_srcterms>\nism_cooling = false\nrel_cooling = false\n",
])
def test_active_cgl_accepts_constant_or_inactive_sources(tmp_path, settings):
    reference = run(tmp_path / "reference", "")
    result = run(tmp_path / "source", settings)
    assert reference.returncode == result.returncode == 0, result.stdout + result.stderr
    states = []
    for path in (tmp_path / "reference", tmp_path / "source"):
        checkpoint = sorted((path / "rst").glob("*.rst"))[-1]
        state = passive_acceptance.decode_restart(checkpoint, False, active_cgl=True)
        states.append(state)
    np.testing.assert_array_equal(states[0]["active"], states[1]["active"])
    for first, second in zip(states[0]["faces"], states[1]["faces"]):
        np.testing.assert_array_equal(first, second)
