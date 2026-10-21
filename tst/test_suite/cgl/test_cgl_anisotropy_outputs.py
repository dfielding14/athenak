"""Check signed anisotropy fields and their 1D/2D PDFs with a supplied executable."""
import importlib.util
import os
from pathlib import Path
import subprocess

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("bin_convert", ROOT / "vis/python/bin_convert.py")
READER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(READER)


def deck(delta_p=0.023, bfield=0.3535533905932738, bfloor=1.e-10):
    text = f"""<job>
basename = anisotropy
<mesh>
nx1 = 16
nx2 = 8
nx3 = 8
x1min = 0
x1max = 1
x2min = 0
x2max = 1
x3min = 0
x3max = 1
ix1_bc = periodic
ox1_bc = periodic
ix2_bc = periodic
ox2_bc = periodic
ix3_bc = periodic
ox3_bc = periodic
nghost = 2
<meshblock>
nx1 = 8
nx2 = 8
nx3 = 8
<time>
evolution = dynamic
integrator = rk2
cfl_number = 0.3
tlim = 0.002
nlim = 2
ndiag = 1
<mhd>
eos = cgl
passive = false
reconstruct = plm
rsolver = hlle
gamma = 1.6666666666666667
nu_coll = 0
mirror_limiter = false
firehose_limiter = false
bfloor = {bfloor}
<problem>
pgen_name = cgl_lf_paper
paper_mode = turbulence
passive_delta = false
rho0 = 1
b0 = {bfield}
beta0 = 16
p_parallel0 = 1
p_perp0 = {1 + delta_p}
"""
    for number, variable in enumerate(("mhd_w_bcc", "mhd_delta_p", "mhd_beta_delta"), 1):
        text += f"""<output{number}>
file_type = bin
variable = {variable}
id = {variable}
dt = 0.001
ghost_zones = false
"""
    for number, variable in enumerate(("mhd_delta_p", "mhd_beta_delta"), 4):
        lo, hi = (-0.1, 0.1) if number == 4 else (-2., 2.)
        text += f"""<output{number}>
file_type = pdf
variable = {variable}
id = {variable}
dt = 0.001
data_format = %25.17e
bin_min = {lo}
bin_max = {hi}
nbin = 32
logscale = false
"""
    # Both derived axes must survive allocation and reset correctly on later writes.
    text += """<output6>
file_type = pdf
variable = mhd_delta_p
variable_2 = mhd_beta_delta
id = joint
dt = 0.001
data_format = %25.17e
bin_min = -0.1
bin_max = 0.1
nbin = 32
logscale = false
bin2_min = -2
bin2_max = 2
nbin2 = 32
logscale2 = false
"""
    return text


def run(tmp_path, text):
    binary = Path(os.environ.get("ATHENAK_TEST_EXE", "./athena")).resolve()
    path = tmp_path / "input.athinput"
    path.write_text(text)
    result = subprocess.run([str(binary), "-i", str(path)], cwd=tmp_path,
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=90)
    (tmp_path / "run.log").write_text(result.stdout)
    return result


@pytest.mark.parametrize("delta_p,bfield,bfloor", [
    (0.023, np.sqrt(1/8), 1.e-10),
    (-0.017, np.sqrt(1/8), 1.e-10),
    (0., np.sqrt(1/8), 1.e-10),
    (0.023, 0.05, 0.1),
])
def test_anisotropy_fields_and_pdfs(tmp_path, delta_p, bfield, bfloor):
    result = run(tmp_path, deck(delta_p, bfield, bfloor))
    assert result.returncode == 0, result.stdout
    fields = sorted((tmp_path / "bin").glob("*.mhd_w_bcc.*.bin"))
    assert len(fields) >= 2
    for path in fields:
        raw = {key: np.asarray(value) for key, value in
               READER.read_binary(str(path))["mb_data"].items()}
        dp = raw["p_perp"] - raw["eint"]
        b2 = sum(raw[key]**2 for key in ("bcc1", "bcc2", "bcc3"))
        normalized = 2*dp/np.maximum(b2, bfloor**2)
        for variable, label, expected in (("mhd_delta_p", "delta_p", dp),
                                          ("mhd_beta_delta", "beta_delta", normalized)):
            output = Path(str(path).replace("mhd_w_bcc", variable))
            actual = np.asarray(READER.read_binary(str(output))["mb_data"][label])
            np.testing.assert_allclose(actual, expected, rtol=2.e-5, atol=2.e-5)
            assert np.all(np.isfinite(actual))
    for directory in tmp_path.glob("pdf_*"):
        bins = np.loadtxt(next(directory.glob("*.bins.pdf")))
        for path in directory.glob("*.pdf"):
            if path.name.endswith(".bins.pdf"):
                continue
            histogram = np.loadtxt(path)
            assert np.isclose(histogram.sum(), 1.)
            assert np.count_nonzero(histogram) == 1
            # Check the populated bin against the saved scalar field, including sign.
            index = path.name.split(".")[-2]
            dp_field = READER.read_binary(str(tmp_path / "bin" /
                f"anisotropy.mhd_delta_p.{index}.bin"))["mb_data"]["delta_p"][0].flat[0]
            bd_field = READER.read_binary(str(tmp_path / "bin" /
                f"anisotropy.mhd_beta_delta.{index}.bin"))["mb_data"]["beta_delta"][0].flat[0]
            if bins.ndim == 2:
                x = np.searchsorted(bins[0], dp_field, side="right")
                y = np.searchsorted(bins[1], bd_field, side="right")
                assert histogram[y, x] == pytest.approx(1.)
            else:
                value = dp_field if directory.name == "pdf_mhd_delta_p" else bd_field
                assert histogram[np.searchsorted(bins, value, side="right")] == pytest.approx(1.)


@pytest.mark.parametrize("second_axis", [False, True])
def test_anisotropy_requires_cgl_on_either_pdf_axis(tmp_path, second_axis):
    text = deck().split("<output1>")[0]
    text = text.replace("eos = cgl", "eos = ideal")
    text = text.replace("pgen_name = cgl_lf_paper", "pgen_name = turb")
    text += """<output1>
file_type = pdf
variable = mhd_delta_p
dt = 1
bin_min = -1
bin_max = 1
nbin = 16
logscale = false
"""
    if second_axis:
        text = text.replace("variable = mhd_delta_p", "variable = mhd_w_d")
        text += """variable_2 = mhd_beta_delta
bin2_min = -2
bin2_max = 2
nbin2 = 16
logscale2 = false
"""
    result = run(tmp_path, text)
    assert result.returncode != 0
    assert "requires <mhd>/eos = cgl" in result.stdout
