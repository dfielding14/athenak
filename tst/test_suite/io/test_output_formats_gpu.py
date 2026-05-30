"""GPU-targeted derived-variable PDF regression.

This test is named for GPU-suite selection. It may also be run against a CPU
binary to check the producer/readback assertions, but a GPU build is required
to qualify device-memory safety.
"""

from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "vis" / "python"))

from read_pdf import read_pdf  # noqa: E402


def test_scalar_and_spherical_derived_fields_are_device_safe(tmp_path: Path):
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
    assert pdf["pdf"].shape == (4, 4, 4, 4)
    assert np.isfinite(pdf["pdf"]).all()
    assert pdf["pdf"].sum() > 0.0
