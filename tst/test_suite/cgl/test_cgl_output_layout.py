"""CGL history capacity at the configured column limit."""

import os
from pathlib import Path
import subprocess

import numpy as np
import pytest

from test_suite.cgl import passive_acceptance


ROOT = Path(__file__).resolve().parents[3]
TEMPLATE = ROOT / "inputs/cgl_lf_paper/cgl_lf_paper_smoke_active_beta10.athinput"


def run(directory, source):
    directory.mkdir()
    deck = directory / "input.athinput"
    deck.write_text(source)
    executable = Path(os.environ.get("ATHENAK_TEST_EXE", "athena")).resolve()
    result = subprocess.run([str(executable), "-i", str(deck)], cwd=directory,
                            capture_output=True, text=True, timeout=90, check=False)
    (directory / "run.log").write_text(result.stdout + result.stderr)
    return result


@pytest.mark.parametrize("mode,scalar_limit", [("active", 34), ("isothermal", 31),
                                              ("adiabatic", 30)])
def test_cgl_history_checks_complete_capacity(tmp_path, mode, scalar_limit):
    source = TEMPLATE.read_text().split("<output1>", 1)[0] + """
<mesh>
nx3 = 8
x3max = 1
<meshblock>
nx3 = 8
<time>
nlim = 0
<turb_driving>
dedt = 0
<output1>
file_type = hst
dcycle = 1
data_format = %24.16e
"""
    if mode != "active":
        source += f"""
<mhd>
passive = true
passive_eos = {mode}
passive_max_pressure_ratio = 100
iso_sound_speed = 1
<problem>
passive_delta = true
"""
    for extra in (0, 1):
        directory = tmp_path / str(extra)
        result = run(directory, source + f"\n<mhd>\nnscalars = {scalar_limit + extra}\n")
        log = result.stdout + result.stderr
        if extra:
            assert result.returncode != 0, log
            assert "MHD history requires 65 columns, exceeding NHISTORY_VARIABLES=64" in log
        else:
            assert result.returncode == 0, log
            history = passive_acceptance.history(next(directory.glob("*.mhd.hst")))
            assert len(history) == 66  # time and dt precede the 64 history variables.
            assert np.all(np.isfinite(np.array(list(history.values()))))

