"""The benchmark must remove all SGS block fragments without altering the physics."""

from pathlib import Path
import runpy


def test_benchmark_inputs_preserve_overrides_and_remove_only_sgs():
    script = Path(__file__).resolve().parents[3] / "scripts" / "benchmark_sgs.py"
    prepare = runpy.run_path(str(script))["prepare_inputs"]
    source = """<time>
nlim = 100
<output1>
file_type = hst
dt = 0.1
<output2> # SGS
file_type = cbin
variable = hydro_sgs_2d # final SGS labels
dt = 0.1
<par_end>
ignored restart payload
"""
    inputs = prepare(source, ["time/nlim=3", "output2/coarsen_factor=512"])
    assert "nlim = 3" in inputs["on"] and "nlim = 3" in inputs["off"]
    assert "coarsen_factor = 512" in inputs["on"]
    assert "<output2>" not in inputs["off"]
    assert "<output1>\nfile_type = hst\ndt = 0.1" in inputs["off"]
    assert "ignored restart payload" not in inputs["on"]
