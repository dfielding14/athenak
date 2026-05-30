"""MPI shared-versus-sharded writer-reader comparison for new output formats."""

from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "vis" / "python"))

from read_pdf import read_pdf, read_pdf_header  # noqa: E402
from read_sphslice import read_sphslice  # noqa: E402


def _run(
    tmp_path: Path, name: str, *overrides: str, input_file="inputs/io_formats.athinput"
):
    run_dir = tmp_path / name
    run_dir.mkdir()
    subprocess.run(
        [
            "mpirun",
            "-np",
            "2",
            "./athena",
            "-i",
            str(input_file),
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


def test_rank_and_node_sharded_pdf_and_sphslice_match_shared_output(tmp_path):
    shared = _run(tmp_path, "shared")
    rank = _run(
        tmp_path,
        "rank",
        "output1/single_file_per_rank=true",
        "output2/single_file_per_rank=true",
    )
    node_input = tmp_path / "io_formats_node.athinput"
    node_input.write_text(
        Path("inputs/io_formats.athinput").read_text()
        .replace("<output1>\n", "<output1>\nsingle_file_per_node = true\n")
        .replace("<output2>\n", "<output2>\nsingle_file_per_node = true\n")
    )
    node = _run(tmp_path, "node", input_file=node_input)
    shared_pdf = read_pdf(
        str(shared / "pdf_nd3_coord_abscostheta_vel_sph_r" / "io_formats.00000.pdf")
    )
    rank_pdf_dir = rank / "pdf_nd3_coord_abscostheta_vel_sph_r"
    rank_pdf = read_pdf(
        str(
            rank_pdf_dir
            / "rank_00000000"
            / "io_formats.00000.pdf"
        )
    )
    node_pdf_dir = node / "pdf_nd3_coord_abscostheta_vel_sph_r"
    node_pdf = read_pdf(
        str(
            node_pdf_dir
            / "node_00000000"
            / "io_formats.00000.pdf"
        )
    )
    shared_surface = read_sphslice(
        str(shared / "bin" /
            "io_formats.density.r_2.5000000000000000e-01.00000.sph.bin")
    )
    rank_surface = read_sphslice(
        str(
            rank
            / "bin"
            / "rank_00000000"
            / "io_formats.density.r_2.5000000000000000e-01.00000.sph.bin"
        )
    )
    node_surface = read_sphslice(
        str(
            node
            / "bin"
            / "node_00000000"
            / "io_formats.density.r_2.5000000000000000e-01.00000.sph.bin"
        )
    )
    np.testing.assert_allclose(rank_pdf["pdf"], shared_pdf["pdf"])
    np.testing.assert_allclose(node_pdf["pdf"], shared_pdf["pdf"])
    np.testing.assert_array_equal(rank_surface["data"], shared_surface["data"])
    np.testing.assert_array_equal(node_surface["data"], shared_surface["data"])
    assert read_pdf_header(
        str(rank_pdf_dir / "rank_00000000" / "io_formats.header.pdf")
    )["rank"] == 0
    assert rank_pdf["header"]["number_of_ranks"] == 2
    assert "rank" not in rank_pdf["header"]
    assert read_pdf_header(
        str(node_pdf_dir / "node_00000000" / "io_formats.header.pdf")
    )["node"] == 0
    assert node_pdf["header"]["number_of_nodes"] == 1
    assert "node" not in node_pdf["header"]
    assert not list(rank_pdf_dir.rglob("*.tmp"))
    assert not list(node_pdf_dir.rglob("*.tmp"))
