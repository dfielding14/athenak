"""The benchmark must remove all SGS block fragments without altering the physics."""

from pathlib import Path
import runpy

import pytest


@pytest.mark.parametrize("variable", ["hydro_sgs_2d", "hydro_sgs_3d", "mhd_sgs"])
def test_benchmark_inputs_preserve_overrides_and_remove_only_sgs(variable):
    script = Path(__file__).resolve().parents[3] / "scripts" / "benchmark_sgs.py"
    prepare = runpy.run_path(str(script))["prepare_inputs"]
    source = """<time>
nlim = 100
<output1>
file_type = hst
dt = 0.1
<output2> # SGS
file_type = cbin
variable = SGS_VARIABLE # SGS fields
dt = 0.1
<output3>
file_type = cbin
variable = hydro_u
<par_end>
ignored restart payload
""".replace("SGS_VARIABLE", variable)
    inputs = prepare(source, ["time/nlim=3", "output2/coarsen_factor=512"])
    assert "nlim = 3" in inputs["on"] and "nlim = 3" in inputs["off"]
    assert "coarsen_factor = 512" in inputs["on"]
    assert "<output2>" not in inputs["off"]
    assert "<output1>\nfile_type = hst\ndt = 0.1" in inputs["off"]
    assert "<output3>\nfile_type = cbin\nvariable = hydro_u" in inputs["off"]
    assert "ignored restart payload" not in inputs["on"]


def test_benchmark_removes_all_sgs_products_together():
    script = Path(__file__).resolve().parents[3] / "scripts" / "benchmark_sgs.py"
    prepare = runpy.run_path(str(script))["prepare_inputs"]
    source = "<time>\nnlim = 3\n" + "".join(
        f"<output{i}>\nfile_type = cbin\nvariable = {variable}\n"
        for i, variable in enumerate(("hydro_sgs_2d", "hydro_sgs_3d", "mhd_sgs"), 1)
    )
    assert prepare(source, [])["off"] == "<time>\nnlim = 3\n"
