"""Serial writer-reader regression for legacy PDFs, N-D PDFs, and spherical slices."""

from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "tst" / "fixtures" / "io" / "origin_main_886dd2a1"
sys.path.insert(0, str(ROOT / "vis" / "python"))

from bin_convert import read_binary  # noqa: E402
from read_pdf import read_pdf  # noqa: E402
from read_sphslice import read_sphslice  # noqa: E402


def _run(tmp_path: Path, input_file: str):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    subprocess.run(
        ["./athena", "-i", input_file, "-d", str(run_dir)],
        check=True,
        capture_output=True,
        text=True,
    )
    return run_dir


def test_legacy_pdf_writer_remains_byte_compatible(tmp_path):
    run_dir = _run(
        tmp_path,
        str(FIXTURES / "producer" / "origin_main_legacy_shared.athinput"),
    )
    generated = run_dir / "pdf_pdf_legacy"
    frozen = FIXTURES / "pdf" / "legacy_1d"
    for filename in (
        "io_legacy_shared.bins.pdf",
        "io_legacy_shared.00000.pdf",
        "io_legacy_shared.00001.pdf",
    ):
        assert (generated / filename).read_bytes() == (frozen / filename).read_bytes()


def test_legacy_two_dimensional_pdf_writer_remains_byte_compatible(tmp_path):
    run_dir = _run(
        tmp_path,
        str(FIXTURES / "producer" / "origin_main_legacy_pdf_2d.athinput"),
    )
    generated = run_dir / "pdf_pdf_legacy_2d_hydro_w_vx"
    frozen = FIXTURES / "pdf" / "legacy_2d"
    for filename in (
        "io_legacy_pdf_2d.bins.pdf",
        "io_legacy_pdf_2d.00000.pdf",
        "io_legacy_pdf_2d.00001.pdf",
    ):
        assert (generated / filename).read_bytes() == (frozen / filename).read_bytes()

    pdf = read_pdf(str(generated / "io_legacy_pdf_2d.00000.pdf"))
    assert pdf["header"]["format"] == "legacy_dense"
    assert pdf["pdf"].shape == (10, 6)


def test_modern_pdf_and_sphslice_round_trip(tmp_path):
    run_dir = _run(tmp_path, "inputs/io_formats.athinput")
    pdf = read_pdf(
        str(
            run_dir
            / "pdf_nd3_coord_abscostheta_vel_sph_r"
            / "io_formats.00000.pdf"
        )
    )
    surface = read_sphslice(
        str(run_dir / "bin" / "io_formats.density.r_0.25.00000.sph.bin")
    )
    legacy = read_pdf(str(run_dir / "pdf_legacy" / "io_formats.00000.pdf"))

    assert pdf["header"]["format"] == "dense"
    assert [axis["scale"] for axis in pdf["header"]["dimensions"]] == [
        "log",
        "linear",
        "symlog",
    ]
    assert pdf["pdf"].shape == (6, 6, 6)
    assert np.isfinite(pdf["pdf"]).all()
    assert surface["data"].shape == (4, 8, 1)
    assert surface["variables"] == ["dens"]
    assert surface["radius"] == 0.25
    assert legacy["header"]["format"] == "legacy_dense"
    assert (
        run_dir / "sph" / "io_formats.r=0.25.legacy_surface.00000.vtk"
    ).exists()


def test_four_dimensional_scalar_weighted_and_volume_pdf_outputs(tmp_path):
    run_dir = _run(tmp_path, "inputs/io_pdf_extended.athinput")
    scalar_weighted = read_pdf(
        str(
            run_dir
            / "pdf_nd4_coord_r_hydro_w_s_0_vel_sph_r"
            / "io_pdf_extended.00000.pdf"
        )
    )
    volume_weighted = read_pdf(
        str(run_dir / "pdf_volume" / "io_pdf_extended.00000.pdf")
    )

    assert scalar_weighted["header"]["ndim"] == 4
    assert scalar_weighted["header"]["weight"] == "variable"
    assert scalar_weighted["header"]["weight_variable"] == "hydro_u_s_0"
    assert [entry["variable"] for entry in scalar_weighted["header"]["dimensions"]] == [
        "coord_x",
        "coord_r",
        "hydro_w_s_0",
        "vel_sph_r",
    ]
    assert scalar_weighted["pdf"].shape == (4, 4, 4, 4)
    assert scalar_weighted["pdf"].sum() > 0.0
    assert volume_weighted["header"]["weight"] == "volume"
    assert volume_weighted["pdf"].shape == (6,)


@pytest.mark.parametrize(
    ("overrides", "expected"),
    (
        (("output1/scale2=bogus",), "invalid scale2"),
        (
            ("output1/bin2_min=-1.0",),
            "requires positive bounds for logarithmic dimension 2",
        ),
        (
            ("output1/linthresh4=-0.1",),
            "requires positive linthresh for symlog dimension 4",
        ),
        (
            ("output1/linthresh4=nan",),
            "requires positive linthresh for symlog dimension 4",
        ),
        (
            ("output1/bin4_max=inf",),
            "requires finite bounds for dimension 4",
        ),
    ),
)
def test_pdf_invalid_axis_configuration_is_rejected(tmp_path, overrides, expected):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "./athena",
            "-i",
            "inputs/io_pdf_extended.athinput",
            "-d",
            str(run_dir),
            *overrides,
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert expected in proc.stdout


def test_sphslice_rejects_derived_field_until_sampling_is_ghost_safe(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "./athena",
            "-i",
            "inputs/io_formats.athinput",
            "-d",
            str(run_dir),
            "output2/variable=coord_r",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "requires derived-field interpolation" in proc.stdout
    assert "ghost-zone-safe sampling" in proc.stdout


@pytest.mark.parametrize("factor", (0, 1, 3, 16))
def test_cbin_rejects_invalid_coarsen_factor_before_writer_construction(
    tmp_path, factor
):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "./athena",
            "-i",
            "inputs/io_node_sharding.athinput",
            "-d",
            str(run_dir),
            f"output3/coarsen_factor={factor}",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "requires coarsen_factor to be a power of two between 2" in proc.stdout


@pytest.mark.parametrize(
    ("extra", "factor", "expected"),
    (
        (
            "ghost_zones = true\n",
            8,
            "Coarsened-binary output extents must be divisible by coarsen_factor.",
        ),
        (
            "slice_x1 = 0.25\n",
            2,
            "Sliced coarsened-binary output is not supported.",
        ),
    ),
)
def test_cbin_rejects_incompatible_emitted_extents_during_construction(
    tmp_path, extra, factor, expected
):
    input_file = tmp_path / "invalid_cbin_extent.athinput"
    text = Path("inputs/io_node_sharding.athinput").read_text()
    text = text.replace(
        "<output3>\nfile_type = cbin\n",
        "<output3>\nfile_type = cbin\n" + extra,
        1,
    )
    input_file.write_text(text)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "./athena",
            "-i",
            str(input_file),
            "-d",
            str(run_dir),
            f"output3/coarsen_factor={factor}",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert expected in proc.stdout


def test_binary_passive_scalar_labels_do_not_wrap_after_99(tmp_path):
    input_file = tmp_path / "scalar_labels.athinput"
    text = Path("inputs/io_pdf_extended.athinput").read_text()
    text = text.replace(
        "<output1>\nfile_type = pdf\n",
        "<output1>\nfile_type = pdf\nvariable = hydro_w_s\n",
        1,
    )
    input_file.write_text(text)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    subprocess.run(
        [
            "./athena",
            "-i",
            str(input_file),
            "-d",
            str(run_dir),
            "hydro/nscalars=101",
            "output1/file_type=bin",
            "time/nlim=0",
            "time/tlim=0.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    output = read_binary(str(run_dir / "bin" / "io_pdf_extended.nd4.00000.bin"))
    assert len(output["var_names"]) == 101
    assert len(set(output["var_names"])) == 101
    assert output["var_names"][-1] == "s_100"


def test_two_fluid_pdf_rejects_unqualified_generic_flux_diagnostic(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "./athena",
            "-i",
            str(ROOT / "tst" / "inputs" / "io_twofluid_diagnostics.athinput"),
            "-d",
            str(run_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "is ambiguous for <ion-neutral> two-fluid runs" in proc.stdout


def test_shipped_readback_example_reads_generated_pdf_and_slice(tmp_path):
    run_dir = _run(tmp_path, "inputs/io_formats.athinput")
    example = ROOT / "vis" / "python" / "examples" / "read_io_outputs.py"
    pdf_path = (
        run_dir
        / "pdf_nd3_coord_abscostheta_vel_sph_r"
        / "io_formats.00000.pdf"
    )
    surface_path = run_dir / "bin" / "io_formats.density.r_0.25.00000.sph.bin"
    pdf = subprocess.run(
        [sys.executable, str(example), "pdf", str(pdf_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    surface = subprocess.run(
        [sys.executable, str(example), "sphslice", str(surface_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "pdf shape=(6, 6, 6)" in pdf.stdout
    assert "sphslice shape=(4, 8, 1)" in surface.stdout
