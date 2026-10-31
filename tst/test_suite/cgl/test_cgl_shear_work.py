"""Homogeneous CGL shear work follows the double-adiabatic invariants."""

import os
from pathlib import Path
import re
import subprocess

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]


def run(directory, dimension, ppar, pperp, dt):
    directory.mkdir()
    source = (ROOT / "tst/inputs/cgl_lf_sts_sbox.athinput").read_text()
    source += f"""
<mesh>
nx1 = 8
nx2 = 8
nx3 = {8 if dimension == 3 else 1}
<meshblock>
nx1 = 8
nx2 = 8
nx3 = {8 if dimension == 3 else 1}
<time>
nlim = 1
tlim = {dt}
<shearing_box>
qshear = 1
<problem>
amp = 0
bx0 = 1
by0 = {1 if dimension == 3 else 0}
bz0 = {0 if dimension == 3 else 1}
ppar0 = {ppar}
pperp0 = {pperp}
"""
    deck = directory / "input.athinput"
    deck.write_text(source)
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()
    result = subprocess.run([str(executable), "-i", str(deck)], cwd=directory,
                            capture_output=True, text=True, timeout=60)
    (directory / "run.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    lines = (directory / "cgl_lf_sts_sbox.mhd.hst").read_text().splitlines()
    labels = re.findall(r"\[\d+\]=(\S+)", lines[1])
    values = dict(zip(labels, np.array(lines[-1].split(), dtype=float)))
    # The box has volume 1/8, and this state remains spatially uniform.
    energy = values["tot-E"] / values["mass"]
    magnetic = sum(values[f"{axis}-ME"] for axis in (1, 2, 3)) / values["mass"]
    for field in ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos"):
        assert values[field] == 0, (field, values[field])
    assert values["time"] == dt
    return energy, energy - magnetic


@pytest.mark.parametrize("dimension", [3, 2])
@pytest.mark.parametrize("ppar,pperp", [(1, 2), (2, 1), (1, 1)])
def test_cgl_homogeneous_shear_work_converges(tmp_path, dimension, ppar, pperp):
    errors = []
    for dt in (0.0002, 0.0001):
        energy, thermal = run(tmp_path / str(dt), dimension, ppar, pperp, dt)
        # q*Omega=1 winds the toroidal field from 1 to 1-t in either orientation.
        b2 = 1 + (1 - dt)**2
        b_ratio = np.sqrt(b2 / 2)
        exact_thermal = pperp*b_ratio + 0.5*ppar/(b_ratio*b_ratio)
        exact_energy = exact_thermal + 0.5*b2
        errors.append(max(abs(thermal - exact_thermal), abs(energy - exact_energy)))
    assert errors[-1] < 2e-8, errors
    # Orbital remapping in 3D is split from the RK source and gives local O(dt^2).
    assert errors[0] / errors[1] > 3.5, errors
