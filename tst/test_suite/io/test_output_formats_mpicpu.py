"""MPI shared-versus-sharded writer-reader comparison for new output formats."""

from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "vis" / "python"))

from read_pdf import read_pdf  # noqa: E402
from read_sphslice import read_sphslice  # noqa: E402


def _run(tmp_path: Path, name: str, *overrides: str):
    run_dir = tmp_path / name
    run_dir.mkdir()
    subprocess.run(
        [
            "mpirun",
            "-np",
            "2",
            "./athena",
            "-i",
            "inputs/io_formats.athinput",
            "-d",
            str(run_dir),
            "mesh/nx1=16",
            "meshblock/nx1=8",
            *overrides,
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return run_dir


def test_rank_sharded_pdf_and_sphslice_match_shared_output(tmp_path):
    shared = _run(tmp_path, "shared")
    sharded = _run(
        tmp_path,
        "sharded",
        "output1/single_file_per_rank=true",
        "output2/single_file_per_rank=true",
    )
    shared_pdf = read_pdf(
        str(shared / "pdf_nd3_coord_abscostheta_vel_sph_r" / "io_formats.00000.pdf")
    )
    sharded_pdf = read_pdf(
        str(
            sharded
            / "pdf_nd3_coord_abscostheta_vel_sph_r"
            / "rank_00000000"
            / "io_formats.00000.pdf"
        )
    )
    shared_surface = read_sphslice(
        str(shared / "bin" / "io_formats.density.r_0.25.00000.sph.bin")
    )
    sharded_surface = read_sphslice(
        str(
            sharded
            / "bin"
            / "rank_00000000"
            / "io_formats.density.r_0.25.00000.sph.bin"
        )
    )
    np.testing.assert_allclose(sharded_pdf["pdf"], shared_pdf["pdf"])
    np.testing.assert_allclose(sharded_surface["data"], shared_surface["data"])
