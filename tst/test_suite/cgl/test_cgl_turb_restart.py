"""The shipped pure-CGL turbulence run restores user history and forcing."""

import os
from pathlib import Path
import subprocess

import numpy as np
import pytest

from test_suite.cgl import passive_acceptance


ROOT = Path(__file__).resolve().parents[3]
BASENAME = "cgl_pure_beta25_accel05_gamma10001_256"


def run(directory, *, mode="cgl", restart=None, nlim=4):
    directory.mkdir()
    command = [str(Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()),
               "-d", str(directory), f"time/nlim={nlim}"]
    if restart is None:
        source = (ROOT / "inputs/cgl_energy_compare" / f"{BASENAME}.athinput").read_text()
        source += f"""
<mesh>
nx1 = 8
nx2 = 8
nx3 = 8
<meshblock>
nx1 = 8
nx2 = 8
nx3 = 8
<mhd>
eos = {mode}
<output1>
dcycle = 1
<output2>
single_file_per_rank = false
dcycle = 1
"""
        deck = directory / "input.athinput"
        deck.write_text(source)
        command += ["-i", str(deck)]
    else:
        command += ["-r", str(restart)]
    result = subprocess.run(command, cwd=directory, capture_output=True,
                            text=True, check=False, timeout=30)
    output = result.stdout + result.stderr
    (directory / "run.log").write_text(output)
    assert result.returncode == 0, output
    checkpoint = sorted((directory / "rst").glob("*.rst"))[-1]
    state = passive_acceptance.decode_restart(checkpoint, False,
                                              active_cgl=mode == "cgl",
                                              ideal_mhd=mode == "ideal")
    history = passive_acceptance.history(directory / f"{BASENAME}.user.hst")
    return checkpoint, state, history


@pytest.mark.parametrize("mode", ["cgl", "ideal"])
def test_turb_restart_preserves_history_and_forcing(tmp_path, mode):
    _, full, full_history = run(tmp_path / "full", mode=mode)
    checkpoint, _, _ = run(tmp_path / "split", mode=mode, nlim=2)
    _, continued, continued_history = run(
        tmp_path / "continued", mode=mode, restart=checkpoint)
    assert full["cycle"] == continued["cycle"] == 4
    assert full["time"] == continued["time"]
    np.testing.assert_array_equal(full["payload"], continued["payload"])
    assert set(full_history) == set(continued_history)
    for key in full_history:
        np.testing.assert_equal(full_history[key][-1], continued_history[key][-1])
