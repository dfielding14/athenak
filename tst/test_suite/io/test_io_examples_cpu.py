"""Executable smoke tests for promoted IO examples."""

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tst" / "fixtures" / "io" / "origin_main_886dd2a1"
sys.path.insert(0, str(ROOT / "vis" / "python"))

from read_pdf import read_pdf  # noqa: E402
from read_sphslice import read_sphslice  # noqa: E402


def _run(tmp_path: Path, input_file: Path, *overrides: str):
    run_dir = tmp_path / input_file.stem
    run_dir.mkdir()
    proc = subprocess.run(
        ["./athena", "-i", str(input_file), "-d", str(run_dir), *overrides],
        check=True,
        capture_output=True,
        text=True,
    )
    return run_dir, proc.stdout


def test_output_formats_example_generates_readable_pdf_and_slice(tmp_path):
    run_dir, _ = _run(tmp_path, ROOT / "inputs" / "io" / "output_formats.athinput")
    pdf = read_pdf(
        str(
            run_dir
            / "pdf_nd3_coord_abscostheta_vel_sph_r"
            / "io_formats.00000.pdf"
        )
    )
    surface = read_sphslice(
        str(run_dir / "bin" /
            "io_formats.density.r_2.5000000000000000e-01.00000.sph.bin")
    )
    assert pdf["pdf"].shape == (26, 14, 22)
    assert surface["data"].shape == (32, 64, 1)
    assert (
        run_dir / "sph" / "io_formats.r=0.25.legacy_surface.00000.vtk"
    ).exists()


def test_scalar_weight_pdf_example_generates_readable_four_dimensional_output(tmp_path):
    run_dir, _ = _run(
        tmp_path, ROOT / "inputs" / "io" / "output_pdf_scalar_weight.athinput"
    )
    pdf = read_pdf(
        str(
            run_dir
            / "pdf_nd4_coord_r_hydro_w_s_0_vel_sph_r"
            / "io_pdf_scalar.00000.pdf"
        )
    )
    assert pdf["header"]["weight"] == "variable"
    assert pdf["header"]["weight_variable"] == "hydro_u_s_0"
    assert pdf["pdf"].shape == (10, 10, 6, 10)


def test_runtime_policy_example_demonstrates_timing_and_restart_only(tmp_path):
    run_dir, stdout = _run(
        tmp_path,
        ROOT / "inputs" / "io" / "runtime_policy.athinput",
        "time/output_timing=true",
        "time/final_output_policy=restart_only",
    )
    assert "[output-io] event=initial " in stdout
    assert "[output-io] event=final " in stdout
    assert list((run_dir / "rst").glob("*.rst"))


def test_shipped_readback_example_reads_legacy_cbin_fixture():
    fixture = FIXTURES / "cbin" / "shared" / "io_legacy_shared.cbin_shared.00000.cbin"
    proc = subprocess.run(
        [
            sys.executable,
            str(ROOT / "vis" / "python" / "examples" / "read_io_outputs.py"),
            "cbin",
            str(fixture),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "coarsened meshblocks=1 variables=['dens']" in proc.stdout
