"""GPU-selectable smoke regression for derived PDF axes and scalar weighting."""

from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "vis" / "python"))

from read_pdf import read_pdf  # noqa: E402


def test_gpu_pdf_derived_axes_and_scalar_weight_round_trip(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    subprocess.run(
        ["./athena", "-i", "inputs/io_pdf_extended.athinput", "-d", str(run_dir)],
        check=True,
        capture_output=True,
        text=True,
    )

    pdf = read_pdf(
        str(
            run_dir
            / "pdf_nd4_coord_r_hydro_w_s_0_vel_sph_r"
            / "io_pdf_extended.00000.pdf"
        )
    )
    assert pdf["header"]["ndim"] == 4
    assert pdf["header"]["weight"] == "variable"
    assert pdf["header"]["weight_variable"] == "hydro_u_s_0"
    assert [entry["variable"] for entry in pdf["header"]["dimensions"]] == [
        "coord_x",
        "coord_r",
        "hydro_w_s_0",
        "vel_sph_r",
    ]
    assert pdf["pdf"].shape == (4, 4, 4, 4)
    assert np.isfinite(pdf["pdf"]).all()
    assert pdf["pdf"].sum() > 0.0
