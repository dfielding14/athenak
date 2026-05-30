"""Construction-time output namespace and path-component regressions."""

from pathlib import Path
import subprocess

import pytest


BASE_INPUT = "inputs/io_node_sharding.athinput"


def _run(tmp_path: Path, text: str):
    input_file = tmp_path / "namespace.athinput"
    input_file.write_text(text)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    return run_dir, subprocess.run(
        ["./athena", "-i", str(input_file), "-d", str(run_dir)],
        capture_output=True,
        text=True,
        timeout=90,
    )


def test_negative_file_number_is_rejected_before_directory_creation(tmp_path):
    text = Path(BASE_INPUT).read_text().replace(
        "<output1>\n", "<output1>\nfile_number = -1\n", 1
    )
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert (
        "requires file_number to be a non-negative integer"
        in proc.stdout + proc.stderr
    )
    assert not (run_dir / "bin").exists()


def test_exhausted_file_number_reaches_writer_rejection(tmp_path):
    text = Path(BASE_INPUT).read_text().replace(
        "<output1>\n", "<output1>\nfile_number = 2147483647\n", 1
    )
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert "binary output file number 2147483647 is outside the publishable range" in (
        proc.stdout + proc.stderr
    )
    assert (run_dir / "bin").is_dir()


def test_conflicting_feature_output_directory_is_rejected(tmp_path):
    input_file = tmp_path / "namespace.athinput"
    input_file.write_text(Path(BASE_INPUT).read_text())
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "bin").write_text("not a directory")

    proc = subprocess.run(
        ["./athena", "-i", str(input_file), "-d", str(run_dir)],
        capture_output=True,
        text=True,
        timeout=90,
    )

    assert proc.returncode != 0
    assert "exists but is not a directory" in proc.stdout + proc.stderr


def test_duplicate_binary_family_is_rejected_before_directory_creation(tmp_path):
    text = Path(BASE_INPUT).read_text().replace(
        "<output2>\nfile_type = bin\nid = slice\n",
        "<output2>\nfile_type = bin\nid = full\n",
    )
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert "resolve to the same public target family" in proc.stdout + proc.stderr
    assert not (run_dir / "bin").exists()


def test_duplicate_binary_family_is_rejected_independent_of_cadence(tmp_path):
    text = Path(BASE_INPUT).read_text() + """
<output7>
file_type = bin
id = full
variable = hydro_w_d
dt = 2.0
"""
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert "resolve to the same public target family" in proc.stdout + proc.stderr
    assert not (run_dir / "bin").exists()


def test_lexically_equivalent_cross_directory_families_are_rejected(tmp_path):
    text = Path(BASE_INPUT).read_text().replace(
        "basename = io_node", "basename = ../cart/io_node", 1
    ) + """
<output7>
file_type = cart
id = full
variable = hydro_w_d
dt = 1.0
"""
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert "resolve to the same public target family" in proc.stdout + proc.stderr
    assert not (run_dir / "bin").exists()
    assert not (run_dir / "cart").exists()


def test_explicit_negative_vtk_gid_matches_omitted_gid(tmp_path):
    text = Path(BASE_INPUT).read_text() + """
<output7>
file_type = vtk
id = density
variable = hydro_w_d
dt = 1.0

<output8>
file_type = vtk
id = density
gid = -1
variable = hydro_w_d
dt = 1.0
"""
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert "resolve to the same public target family" in proc.stdout + proc.stderr
    assert not (run_dir / "vtk").exists()


def test_duplicate_pdf_family_is_rejected_before_directory_creation(tmp_path):
    text = Path(BASE_INPUT).read_text() + """
<output7>
file_type = pdf
id = node_pdf
variable_1 = coord_y
nbin1 = 8
bin1_min = -0.5
bin1_max = 0.5
scale1 = linear
variable_2 = hydro_w_d
nbin2 = 8
bin2_min = 0.0
bin2_max = 2.0
scale2 = linear
weight = mass
dt = 1.0
"""
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert "pdf_node_pdf_hydro_w_d/io_node.header.pdf" in proc.stdout + proc.stderr
    assert not (run_dir / "pdf_node_pdf_hydro_w_d").exists()


@pytest.mark.parametrize(
    "changed_line",
    (
        "variable_1 = coord_y",
        "nbin1 = 8",
        "scale1 = log",
        "weight = mass",
    ),
)
def test_duplicate_pdf_family_is_rejected_despite_histogram_configuration(
    tmp_path, changed_line
):
    block = """
<output7>
file_type = pdf
id = node_pdf
variable_1 = coord_x
nbin1 = 4
bin1_min = -0.5
bin1_max = 0.5
scale1 = linear
variable_2 = hydro_w_d
nbin2 = 4
bin2_min = 0.0
bin2_max = 2.0
scale2 = linear
weight = volume
dt = 2.0
"""
    original = {
        "variable_1 = coord_y": "variable_1 = coord_x",
        "nbin1 = 8": "nbin1 = 4",
        "scale1 = log": "scale1 = linear",
        "weight = mass": "weight = volume",
    }[changed_line]
    text = Path(BASE_INPUT).read_text() + block.replace(original, changed_line, 1)
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert "resolve to the same public target family" in proc.stdout + proc.stderr
    assert not (run_dir / "pdf_node_pdf_hydro_w_d").exists()


def test_duplicate_cbin_family_is_rejected_before_directory_creation(tmp_path):
    text = Path(BASE_INPUT).read_text() + """
<output7>
file_type = cbin
id = coarse
variable = hydro_w_d
coarsen_factor = 2
compute_moments = true
dt = 2.0
"""
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert "resolve to the same public target family" in proc.stdout + proc.stderr
    assert not (run_dir / "cbin_coarse_2").exists()


def test_duplicate_sphslice_family_is_rejected_before_directory_creation(tmp_path):
    text = Path(BASE_INPUT).read_text() + """
<output7>
file_type = sphslice
id = density
variable = hydro_w_d
slice_r = 0.25
ntheta = 8
nphi = 16
dt = 2.0
"""
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert "resolve to the same public target family" in proc.stdout + proc.stderr
    assert not (run_dir / "bin").exists()


def test_adjacent_representable_sphslice_radii_remain_disjoint(tmp_path):
    text = Path(BASE_INPUT).read_text() + """
<output7>
file_type = sphslice
id = density
variable = hydro_w_d
slice_r = 0.25000000000000006
ntheta = 4
nphi = 8
dt = 1.0
"""
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    outputs = sorted((run_dir / "bin").glob("*.sph.bin"))
    assert len(outputs) == 2
    assert outputs[0].name != outputs[1].name


def test_explicit_negative_pvtk_gid_matches_omitted_gid(tmp_path):
    text = Path(BASE_INPUT).read_text() + """
<output7>
file_type = pvtk
id = density
variable = hydro_w_d
dt = 1.0

<output8>
file_type = pvtk
id = density
gid = -1
variable = hydro_w_d
dt = 2.0
"""
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert "resolve to the same public target family" in proc.stdout + proc.stderr
    assert not (run_dir / "pvtk").exists()


def test_unsafe_generated_id_is_rejected_before_directory_creation(tmp_path):
    text = Path(BASE_INPUT).read_text().replace("id = full", "id = ../escape", 1)
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode != 0
    assert "contains an unsafe path component" in proc.stdout + proc.stderr
    assert not (run_dir / "bin").exists()
    assert not (tmp_path / "escape").exists()


def test_distinct_cbin_factors_remain_disjoint(tmp_path):
    text = Path(BASE_INPUT).read_text().replace(
        "<output4>\nfile_type = pdf\n",
        """<output4>
file_type = cbin
id = coarse
variable = hydro_w_d
coarsen_factor = 4
compute_moments = false
single_file_per_rank = false
single_file_per_node = false
dt = 1.0

<output14>
file_type = pdf
""",
    )
    run_dir, proc = _run(tmp_path, text)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (run_dir / "cbin_coarse_2").is_dir()
    assert (run_dir / "cbin_coarse_4").is_dir()
