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


TOTAL_ENERGY_DIAGNOSTICS = (
    "edot_sph",
    "edot_sph_out",
    "edot_sph_in",
    "edot_sph_th",
    "edot_vert",
    "edot_vert_out",
    "edot_vert_in",
)
GENERIC_FLUID_DIAGNOSTICS = (
    "vel_sph_r",
    "vel_sph_theta",
    "vel_sph_phi",
    "vel_cyl_R",
    "vel_cyl_phi",
    "mdot_sph",
    "mdot_sph_out",
    "mdot_sph_in",
    "mdot_vert",
    "mdot_vert_out",
    "mdot_vert_in",
    *TOTAL_ENERGY_DIAGNOSTICS,
    "edot_sph_kin",
    "edot_sph_mag",
)


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


def _isothermal_linear_wave_input(module):
    text = (ROOT / "tst" / "inputs" / f"lwave_{module}.athinput").read_text()
    text = text.replace("eos         = ideal", "eos         = isothermal")
    text = text.replace("basename  = LinWave", f"basename  = isothermal_{module}")
    return text


def _analytic_hydro_input():
    text = Path("inputs/io_formats.athinput").read_text().split("<output1>\n", 1)[0]
    text = text.replace("nx1 = 8", "nx1 = 5")
    text = text.replace("nx2 = 8", "nx2 = 5")
    text = text.replace("nx3 = 8", "nx3 = 5")
    replacements = {
        "dl = 1.0": "dl = 2.0",
        "pl = 1.0": "pl = 5.0",
        "ul = 0.0": "ul = 2.0",
        "vl = 0.0": "vl = -3.0",
        "wl = 0.0": "wl = 4.0",
        "dr = 0.125": "dr = 2.0",
        "pr = 0.1": "pr = 5.0",
        "ur = 0.0": "ur = 2.0",
        "vr = 0.0": "vr = -3.0",
        "wr = 0.0": "wr = 4.0",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    variables = (
        "coord_x", "coord_y", "coord_z", "coord_r", "coord_theta", "coord_phi",
        "coord_costheta", "coord_abscostheta", "coord_cyl_R", "coord_cyl_phi",
        "coord_cyl_z", "vel_sph_r", "vel_sph_theta", "vel_sph_phi", "vel_cyl_R",
        "vel_cyl_phi", "mdot_sph", "mdot_sph_out", "mdot_sph_in", "mdot_vert",
        "mdot_vert_out", "mdot_vert_in", "edot_sph", "edot_sph_out",
        "edot_sph_in", "edot_sph_kin", "edot_sph_th", "edot_vert",
        "edot_vert_out", "edot_vert_in",
    )
    for number, variable in enumerate(variables, start=1):
        text += f"""
<output{number}>
file_type = bin
id = {variable}
variable = {variable}
dt = 1.0
"""
    return text, variables


def _read_analytic_field(run_dir, variable):
    data = read_binary(str(run_dir / "bin" / f"diagnostics.{variable}.00000.bin"))
    return data["mb_data"][variable][0]


def test_uniform_hydro_diagnostics_match_analytic_oracles(tmp_path):
    input_file = tmp_path / "diagnostics.athinput"
    text, variables = _analytic_hydro_input()
    input_file.write_text(text.replace("basename = io_formats", "basename = diagnostics"))
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    subprocess.run(
        [
            "./athena",
            "-i",
            str(input_file),
            "-d",
            str(run_dir),
            "time/nlim=0",
            "time/tlim=0.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    centers = np.linspace(-0.4, 0.4, 5)
    z, y, x = np.meshgrid(centers, centers, centers, indexing="ij")
    cylindrical_radius = np.sqrt(x*x + y*y)
    radius = np.sqrt(cylindrical_radius*cylindrical_radius + z*z)
    costheta = np.divide(z, radius, out=np.ones_like(radius), where=radius > 0.0)
    theta = np.arccos(np.clip(costheta, -1.0, 1.0))
    phi = np.mod(np.arctan2(y, x), 2.0*np.pi)
    vx, vy, vz, rho = 2.0, -3.0, 4.0, 2.0
    radial_velocity = np.divide(
        vx*x + vy*y + vz*z, radius, out=np.zeros_like(radius), where=radius > 0.0
    )
    cylindrical_velocity = np.divide(
        vx*x + vy*y, cylindrical_radius, out=np.zeros_like(radius),
        where=cylindrical_radius > 0.0,
    )
    phi_velocity = np.divide(
        -vx*y + vy*x, cylindrical_radius, out=np.zeros_like(radius),
        where=cylindrical_radius > 0.0,
    )
    theta_velocity = np.divide(
        z*(vx*x + vy*y), radius*cylindrical_radius,
        out=np.zeros_like(radius), where=(radius > 0.0) & (cylindrical_radius > 0.0),
    ) - np.divide(
        vz*cylindrical_radius, radius, out=np.zeros_like(radius), where=radius > 0.0
    )
    sign_z = np.sign(z)
    radial_mass_flux = rho*radial_velocity
    vertical_mass_flux = rho*vz*sign_z
    kinetic_radial = 29.0*radial_velocity
    thermal_radial = 17.5*radial_velocity
    total_radial = 46.5*radial_velocity
    total_vertical = 46.5*vz*sign_z
    expected = {
        "coord_x": x, "coord_y": y, "coord_z": z, "coord_r": radius,
        "coord_theta": theta, "coord_phi": phi, "coord_costheta": costheta,
        "coord_abscostheta": np.abs(costheta), "coord_cyl_R": cylindrical_radius,
        "coord_cyl_phi": phi, "coord_cyl_z": z, "vel_sph_r": radial_velocity,
        "vel_sph_theta": theta_velocity, "vel_sph_phi": phi_velocity,
        "vel_cyl_R": cylindrical_velocity, "vel_cyl_phi": phi_velocity,
        "mdot_sph": radial_mass_flux, "mdot_sph_out": np.maximum(radial_mass_flux, 0.0),
        "mdot_sph_in": np.minimum(radial_mass_flux, 0.0), "mdot_vert": vertical_mass_flux,
        "mdot_vert_out": np.maximum(vertical_mass_flux, 0.0),
        "mdot_vert_in": np.minimum(vertical_mass_flux, 0.0), "edot_sph": total_radial,
        "edot_sph_out": np.where(radial_velocity > 0.0, total_radial, 0.0),
        "edot_sph_in": np.where(radial_velocity < 0.0, total_radial, 0.0),
        "edot_sph_kin": kinetic_radial, "edot_sph_th": thermal_radial,
        "edot_vert": total_vertical,
        "edot_vert_out": np.where(vz*sign_z > 0.0, total_vertical, 0.0),
        "edot_vert_in": np.where(vz*sign_z < 0.0, total_vertical, 0.0),
    }
    assert set(expected) == set(variables)
    for variable in variables:
        np.testing.assert_allclose(_read_analytic_field(run_dir, variable),
                                   expected[variable], rtol=1.0e-6, atol=1.0e-6)


def test_uniform_mhd_energy_diagnostics_match_analytic_oracles(tmp_path):
    text, _ = _analytic_hydro_input()
    text = text.split("<output1>\n", 1)[0]
    text = text.replace("basename = io_formats", "basename = mhd_diagnostics")
    text = text.replace("<hydro>", "<mhd>")
    text += """
bxl = 1.0
byl = 2.0
bzl = -1.0
bxr = 1.0
byr = 2.0
bzr = -1.0

<output1>
file_type = bin
id = edot_sph_mag
variable = edot_sph_mag
dt = 1.0

<output2>
file_type = bin
id = edot_sph
variable = edot_sph
dt = 1.0

<output3>
file_type = bin
id = edot_vert
variable = edot_vert
dt = 1.0
"""
    input_file = tmp_path / "mhd_diagnostics.athinput"
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
            "time/nlim=0",
            "time/tlim=0.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    centers = np.linspace(-0.4, 0.4, 5)
    z, y, x = np.meshgrid(centers, centers, centers, indexing="ij")
    radius = np.sqrt(x*x + y*y + z*z)
    radial_velocity = np.divide(
        2.0*x - 3.0*y + 4.0*z, radius, out=np.zeros_like(radius),
        where=radius > 0.0,
    )
    radial_magnetic = np.divide(
        x + 2.0*y - z, radius, out=np.zeros_like(radius), where=radius > 0.0
    )
    expected = {
        "edot_sph_mag": 6.0*radial_velocity + 8.0*radial_magnetic,
        "edot_sph": 52.5*radial_velocity + 8.0*radial_magnetic,
        "edot_vert": 202.0*np.sign(z),
    }
    for variable, values in expected.items():
        data = read_binary(
            str(run_dir / "bin" / f"mhd_diagnostics.{variable}.00000.bin")
        )
        np.testing.assert_allclose(data["mb_data"][variable][0], values,
                                   rtol=1.0e-6, atol=1.0e-6)


def test_uniform_pdf_matches_numpy_histogram_oracle(tmp_path):
    text, _ = _analytic_hydro_input()
    text = text.split("<output1>\n", 1)[0]
    text = text.replace("basename = io_formats", "basename = pdf_oracle")
    text += """
<output1>
file_type = pdf
id = coord_x
variable_1 = coord_x
nbin1 = 4
bin1_min = -0.5
bin1_max = 0.5
scale1 = linear
weight = volume
dt = 1.0
"""
    input_file = tmp_path / "pdf_oracle.athinput"
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
            "time/nlim=0",
            "time/tlim=0.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    result = read_pdf(str(run_dir / "pdf_coord_x" / "pdf_oracle.00000.pdf"))
    edges = np.linspace(-0.5, 0.5, 5)
    samples = np.repeat(np.linspace(-0.4, 0.4, 5), 25)
    expected = np.zeros(6)
    expected[1:-1] = np.histogram(samples, bins=edges)[0] * 0.2**3
    np.testing.assert_allclose(result["header"]["dimensions"][0]["bin_edges"], edges)
    np.testing.assert_allclose(result["pdf"], expected)


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
        str(run_dir / "bin" /
            "io_formats.density.r_2.5000000000000000e-01.00000.sph.bin")
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
        (
            ("output1/linthresh4=4.9406564584124654e-324",),
            "requires finite transformed bounds and a positive finite bin step "
            "for dimension 4",
        ),
        (
            (
                "output2/bin1_min=0.0",
                "output2/bin1_max=4.9406564584124654e-324",
                "output2/nbin1=4",
            ),
            "requires finite transformed bounds and a positive finite bin step "
            "for dimension 1",
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
    assert expected in (proc.stdout + proc.stderr)


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
    assert "requires derived-field interpolation" in (proc.stdout + proc.stderr)
    assert "ghost-zone-safe sampling" in (proc.stdout + proc.stderr)


@pytest.mark.parametrize(
    ("block", "expected"),
    (
        ("output1",
         "PDF persistent allocation requires"),
        ("output2",
         "sphslice allocation requires"),
    ),
)
def test_writers_reject_reduced_allocation_caps_before_large_allocations(
    tmp_path, block, expected
):
    input_file = tmp_path / "reduced_cap.athinput"
    input_file.write_text(
        Path("inputs/io_formats.athinput").read_text().replace(
            f"<{block}>\n", f"<{block}>\nmax_writer_allocation_bytes = 1\n", 1
        )
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "./athena",
            "-i",
            str(input_file),
            "-d",
            str(run_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert expected in (proc.stdout + proc.stderr)


def test_pdf_rejects_load_cap_before_derived_field_allocation(tmp_path):
    input_file = tmp_path / "derived_load_cap.athinput"
    input_file.write_text(
        Path("inputs/io_pdf_extended.athinput").read_text().replace(
            "<output1>\n", "<output1>\nmax_writer_allocation_bytes = 7000\n", 1
        )
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        ["./athena", "-i", str(input_file), "-d", str(run_dir)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "PDF load allocation requires" in (proc.stdout + proc.stderr)


def test_pdf_rejects_explicit_ghost_zone_sampling(tmp_path):
    input_file = tmp_path / "pdf_ghosts.athinput"
    input_file.write_text(
        Path("inputs/io_formats.athinput").read_text().replace(
            "<output3>\n", "<output3>\nghost_zones = true\n", 1
        )
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "./athena",
            "-i",
            str(input_file),
            "-d",
            str(run_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "PDFs sample active zones only" in (proc.stdout + proc.stderr)


@pytest.mark.parametrize(
    ("input_file", "overrides", "expected"),
    (
        (
            "inputs/io_node_sharding.athinput",
            ("output1/variable=coord_r", "output1/ghost_zones=true"),
            "cannot set ghost_zones=true for derived diagnostic output",
        ),
        (
            "inputs/io_formats.athinput",
            ("output4/variable=coord_r",),
            "Spherical-surface derived-field interpolation is not supported",
        ),
    ),
)
def test_unsafe_derived_sampling_paths_are_rejected(
    tmp_path, input_file, overrides, expected
):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    configured_input = input_file
    if "output1/ghost_zones=true" in overrides:
        configured_input = str(tmp_path / "derived_ghosts.athinput")
        Path(configured_input).write_text(
            Path(input_file).read_text().replace(
                "<output1>\n", "<output1>\nghost_zones = true\n", 1
            )
        )
        overrides = tuple(
            item for item in overrides if item != "output1/ghost_zones=true"
        )
    proc = subprocess.run(
        ["./athena", "-i", configured_input, "-d", str(run_dir), *overrides],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert expected in (proc.stdout + proc.stderr)


def test_cartesian_grid_rejects_derived_interpolation(tmp_path):
    input_file = tmp_path / "cart_derived.athinput"
    input_file.write_text(
        Path("inputs/io_formats.athinput").read_text()
        + """
<output5>
file_type = cart
id = derived
variable = coord_r
dt = 1.0
"""
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        ["./athena", "-i", str(input_file), "-d", str(run_dir)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "Cartesian-grid derived-field interpolation is not supported" in (
        proc.stdout + proc.stderr
    )


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
    assert "requires coarsen_factor to be a power of two between 2" in (
        proc.stdout + proc.stderr
    )


@pytest.mark.parametrize(
    ("extra", "factor", "expected"),
    (
        (
            "ghost_zones = true\n",
            8,
            "Coarsened-binary output does not support ghost_zones=true.",
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
    assert expected in (proc.stdout + proc.stderr)
    assert not (run_dir / f"cbin_coarse_{factor}").exists()


@pytest.mark.parametrize(
    ("overrides", "expected"),
    (
        (
            ("mesh/nx2=1", "meshblock/nx2=1", "mesh/nx3=1", "meshblock/nx3=1"),
            "Coarsened-binary output supports three-dimensional meshes only.",
        ),
        (
            ("mesh/nx3=1", "meshblock/nx3=1"),
            "Coarsened-binary output supports three-dimensional meshes only.",
        ),
    ),
)
def test_cbin_rejects_unsupported_mesh_contract_before_publication(
    tmp_path, overrides, expected
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
            *overrides,
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert expected in (proc.stdout + proc.stderr)
    assert not (run_dir / "cbin_coarse_2").exists()


@pytest.mark.parametrize(
    "refinement",
    (
        "static",
        (
            "adaptive\nnum_levels = 2\nmax_nmb_per_rank = 64\n"
            "\n<amr_criterion0>\nmethod = slope\nvariable = hydro_w_d\n"
            "value_max = 0.1"
        ),
    ),
)
def test_cbin_rejects_refinement_before_publication(tmp_path, refinement):
    input_file = tmp_path / "refined_cbin.athinput"
    input_file.write_text(
        Path("inputs/io_node_sharding.athinput").read_text()
        + f"\n<mesh_refinement>\nrefinement = {refinement}\n"
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        ["./athena", "-i", str(input_file), "-d", str(run_dir)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "Coarsened-binary output supports uniform meshes only" in (
        proc.stdout + proc.stderr
    )
    assert not (run_dir / "cbin_coarse_2").exists()


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


@pytest.mark.parametrize("variable", GENERIC_FLUID_DIAGNOSTICS)
def test_two_fluid_pdf_rejects_unqualified_generic_flux_diagnostic(tmp_path, variable):
    input_file = tmp_path / f"twofluid_{variable}.athinput"
    input_file.write_text(
        (ROOT / "tst" / "inputs" / "io_twofluid_diagnostics.athinput").read_text()
        .replace("variable_1 = mdot_sph", f"variable_1 = {variable}")
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        [
            "./athena",
            "-i",
            str(input_file),
            "-d",
            str(run_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "is ambiguous for <ion-neutral> two-fluid runs" in (proc.stdout + proc.stderr)


@pytest.mark.parametrize("variable", TOTAL_ENERGY_DIAGNOSTICS)
def test_isothermal_hydro_rejects_total_energy_diagnostic_before_publication(
    tmp_path, variable
):
    input_file = tmp_path / f"isothermal_{variable}.athinput"
    text = _isothermal_linear_wave_input("hydro")
    text += f"""
<output1>
file_type = bin
id = {variable}
variable = {variable}
dt = 1.0
"""
    input_file.write_text(text)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        ["./athena", "-i", str(input_file), "-d", str(run_dir)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert (
        f"Variable '{variable}' requires an ideal-gas total-energy fluid module"
        in (proc.stdout + proc.stderr)
    )
    assert not (run_dir / "bin").exists()


@pytest.mark.parametrize(("module", "variable"), (("hydro", "edot_sph_kin"),
                                                  ("mhd", "edot_sph_mag")))
def test_isothermal_module_accepts_energy_diagnostic_without_total_energy(
    tmp_path, module, variable
):
    input_file = tmp_path / f"isothermal_{module}_{variable}.athinput"
    text = _isothermal_linear_wave_input(module)
    text += f"""
<output1>
file_type = bin
id = {variable}
variable = {variable}
dt = 1.0
"""
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
            "time/nlim=0",
            "time/tlim=0.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    output = read_binary(
        str(run_dir / "bin" / f"isothermal_{module}.{variable}.00000.bin")
    )
    assert np.isfinite(output["mb_data"][variable]).all()


@pytest.mark.parametrize("module", ("hydro", "mhd"))
def test_passive_scalar_indices_are_distinct_and_numerically_correct(tmp_path, module):
    input_file = tmp_path / f"{module}_scalar_indices.athinput"
    text = Path("inputs/io_formats.athinput").read_text().split("<output1>\n", 1)[0]
    text = text.replace("basename = io_formats", f"basename = {module}_scalar_indices")
    text = text.replace("<hydro>", f"<{module}>")
    text = text.replace("gamma = 1.4", "gamma = 1.4\nnscalars = 2")
    text = text.replace("wr = 0.0", "wr = 0.0\nyl = 0.75\nyr = 0.25")
    if module == "mhd":
        text = text.replace(
            "yr = 0.25",
            "yr = 0.25\nbxl = 1.0\nbyl = 2.0\nbzl = -1.0\n"
            "bxr = 1.0\nbyr = 2.0\nbzr = -1.0",
        )
    variables = tuple(
        f"{module}_{kind}_s_{index}" for kind in ("u", "w") for index in (0, 1)
    )
    for number, variable in enumerate(variables, start=1):
        text += f"""
<output{number}>
file_type = bin
id = {variable}
variable = {variable}
dt = 1.0
"""
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
            "time/nlim=0",
            "time/tlim=0.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    x = np.linspace(-0.4375, 0.4375, 8)
    primitive = np.broadcast_to(np.where(x < 0.0, 0.75, 0.25), (8, 8, 8))
    conserved = primitive * np.broadcast_to(
        np.where(x < 0.0, 1.0, 0.125), (8, 8, 8)
    )
    expected = {"w": primitive, "u": conserved}
    for variable in variables:
        kind = variable.split("_")[1]
        index = variable.rsplit("_", 1)[1]
        output = read_binary(
            str(run_dir / "bin" / f"{module}_scalar_indices.{variable}.00000.bin")
        )
        label = f"{'r' if kind == 'u' else 's'}_{index}"
        assert output["var_names"] == [label]
        np.testing.assert_allclose(output["mb_data"][label][0],
                                   expected[kind])


def test_two_fluid_pdf_rejects_unqualified_mass_weight(tmp_path):
    input_file = tmp_path / "twofluid_mass.athinput"
    input_file.write_text(
        (ROOT / "tst" / "inputs" / "io_twofluid_diagnostics.athinput").read_text()
        .replace("variable_1 = mdot_sph", "variable_1 = coord_x")
        .replace("weight = volume", "weight = mass")
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        ["./athena", "-i", str(input_file), "-d", str(run_dir)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "Mass-weighted PDF output block" in (proc.stdout + proc.stderr)
    assert "ambiguous for <ion-neutral> two-fluid runs" in (proc.stdout + proc.stderr)


def test_shipped_readback_example_reads_generated_pdf_and_slice(tmp_path):
    run_dir = _run(tmp_path, "inputs/io_formats.athinput")
    example = ROOT / "vis" / "python" / "examples" / "read_io_outputs.py"
    pdf_path = (
        run_dir
        / "pdf_nd3_coord_abscostheta_vel_sph_r"
        / "io_formats.00000.pdf"
    )
    surface_path = (
        run_dir / "bin" /
        "io_formats.density.r_2.5000000000000000e-01.00000.sph.bin"
    )
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
