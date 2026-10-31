"""CGL history capacity and complete conserved-field output."""

import os
from pathlib import Path
import subprocess

import numpy as np
import pytest

from test_suite.cgl import passive_acceptance
from test_suite.cgl.test_cgl_anisotropy_outputs import READER


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


@pytest.mark.parametrize("mode,scalar_limit", [
    ("active", 34), ("isothermal", 31), ("adiabatic", 30),
])
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
            message = "MHD history requires 65 columns, exceeding NHISTORY_VARIABLES=64"
            assert message in log
        else:
            assert result.returncode == 0, log
            history = passive_acceptance.history(next(directory.glob("*.mhd.hst")))
            assert len(history) == 66  # time and dt precede the 64 history variables.
            assert np.all(np.isfinite(np.array(list(history.values()))))


def test_active_cgl_conserved_bundles_retain_both_pressures(tmp_path):
    fields = []
    for ppar, pperp in ((1, 2), (3, 1)):
        source = TEMPLATE.read_text().split("<output1>", 1)[0] + f"""
<time>
nlim = 0
<mhd>
nscalars = 1
<problem>
p_parallel0 = {ppar}
p_perp0 = {pperp}
"""
        for number, variable in enumerate(("mhd_u", "mhd_u_bcc", "mhd_u_e"), 1):
            source += f"""
<output{number}>
file_type = bin
variable = {variable}
id = {variable}
dcycle = 1
"""
        directory = tmp_path / str(ppar)
        result = run(directory, source)
        assert result.returncode == 0, result.stdout + result.stderr
        outputs = {name: READER.read_binary(str(sorted(
            (directory / "bin").glob(f"*.{name}.*.bin"))[0]))
            for name in ("mhd_u", "mhd_u_bcc", "mhd_u_e")}
        assert outputs["mhd_u_e"]["var_names"] == ["ener"]
        assert outputs["mhd_u"]["var_names"] == [
            "dens", "mom1", "mom2", "mom3", "ener", "r_00", "cgl_A"]
        assert outputs["mhd_u_bcc"]["var_names"] == [
            "dens", "mom1", "mom2", "mom3", "ener", "r_00",
            "bcc1", "bcc2", "bcc3", "cgl_A"]
        raw = {name: np.asarray(value, dtype=float) for name, value in
               outputs["mhd_u_bcc"]["mb_data"].items()}
        for name in outputs["mhd_u"]["var_names"]:
            np.testing.assert_array_equal(outputs["mhd_u"]["mb_data"][name], raw[name])
        bmag = np.sqrt(sum(raw[name]**2 for name in ("bcc1", "bcc2", "bcc3")))
        energy = raw["ener"] - 0.5*bmag**2
        ratio = np.exp(raw["cgl_A"]/raw["dens"] -
                       2*np.log(raw["dens"]) + 3*np.log(bmag))
        decoded_parallel = energy/(ratio + 0.5)
        np.testing.assert_allclose(decoded_parallel, ppar, rtol=2e-6)
        np.testing.assert_allclose(decoded_parallel*ratio, pperp, rtol=2e-6)
        fields.append(raw)
    np.testing.assert_array_equal(fields[0]["ener"], fields[1]["ener"])
    assert not np.array_equal(fields[0]["cgl_A"], fields[1]["cgl_A"])
