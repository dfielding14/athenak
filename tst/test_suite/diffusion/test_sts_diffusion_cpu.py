"""CPU regression coverage for RKL2 super time stepping."""

from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

import test_suite.testutils as testutils


INPUT_ROOT = "../../../inputs/tests"


def _run(input_name, basename, *flags):
    testutils.run(
        f"{INPUT_ROOT}/{input_name}",
        [f"job/basename={basename}", *flags],
    )


def _scalar_errors(basename):
    return np.loadtxt(f"{basename}-scalar-errors.dat")


def _tab(basename, variable):
    return testutils.athena_read.tab(f"tab/{basename}.{variable}.00001.tab")


def _cleanup():
    testutils.cleanup()
    shutil.rmtree("bin", ignore_errors=True)


def test_sts_independent_scalar_diffusivities_converge():
    """Each scalar coefficient participates in a second-order analytic check."""
    try:
        errors = {}
        for resolution in (32, 64, 128):
            basename = f"sts_scalar_{resolution}"
            _run(
                "sts_scalar_modes.athinput",
                basename,
                f"mesh/nx1={resolution}",
                f"meshblock/nx1={min(64, resolution)}",
            )
            errors[resolution] = _scalar_errors(basename)

        for column in (2, 3):
            assert errors[128][column] < 3.0e-6
            assert errors[64][column] / errors[32][column] < 0.27
            assert errors[128][column] / errors[64][column] < 0.27

        _run(
            "sts_scalar_modes.athinput",
            "explicit_scalar_128",
            "mesh/nx1=128",
            "meshblock/nx1=64",
            "hydro/scalar_diffusivity_integrator=explicit",
            "time/sts_integrator=none",
        )
        explicit = _scalar_errors("explicit_scalar_128")
        assert np.max(np.abs(errors[128][2:4] - explicit[2:4])) < 5.0e-9
    finally:
        _cleanup()


def test_sts_bounded_power_law_conduction_and_viscosity():
    """Power-law coefficients run with STS and agree with explicit small steps."""
    try:
        common = (
            "mesh/nx1=64",
            "meshblock/nx1=64",
            "time/tlim=0.01",
            "output1/dt=0.01",
        )
        _run(
            "sts_thermal_front.athinput",
            "sts_front_cap",
            *common,
            "time/sts_max_dt_ratio=1.0",
        )
        _run(
            "sts_thermal_front.athinput",
            "explicit_front",
            *common,
            "hydro/conductivity_integrator=explicit",
            "time/sts_integrator=none",
        )
        front_sts = _tab("sts_front_cap", "hydro_w")
        front_explicit = _tab("explicit_front", "hydro_w")
        assert np.max(np.abs(front_sts["eint"] - front_explicit["eint"])) < 1.0e-3

        _run("sts_viscous_shear.athinput", "sts_shear", *common)
        _run(
            "sts_viscous_shear.athinput",
            "explicit_shear",
            *common,
            "hydro/viscosity_integrator=explicit",
            "time/sts_integrator=none",
        )
        shear_sts = _tab("sts_shear", "hydro_w")
        shear_explicit = _tab("explicit_shear", "hydro_w")
        assert np.max(np.abs(shear_sts["vely"] - shear_explicit["vely"])) < 1.0e-5
    finally:
        _cleanup()


def test_sts_resistive_field_update_and_advected_blob_smoke():
    """Exercise magnetic CT STS updates and a two-dimensional scalar visualization."""
    try:
        _run("sts_resistivity.athinput", "sts_resistive")
        resistive = _tab("sts_resistive", "mhd_bcc")
        assert np.all(np.isfinite(resistive["bcc2"]))

        _run(
            "sts_resistivity.athinput",
            "sts_resistive_cap",
            "time/sts_max_dt_ratio=1.0",
        )
        _run(
            "sts_resistivity.athinput",
            "explicit_resistive",
            "mhd/ohmic_resistivity_integrator=explicit",
            "time/sts_integrator=none",
        )
        capped = _tab("sts_resistive_cap", "mhd_bcc")
        explicit = _tab("explicit_resistive", "mhd_bcc")
        assert np.max(np.abs(capped["bcc2"] - explicit["bcc2"])) < 1.0e-2

        _run(
            "sts_scalar_blob.athinput",
            "sts_blob",
            "mesh/nx1=32",
            "mesh/nx2=32",
            "meshblock/nx1=32",
            "meshblock/nx2=32",
            "time/tlim=0.02",
            "output1/dt=0.02",
            "output2/dt=0.02",
        )
        start = testutils.athena_read.tab("tab/sts_blob.hydro_w.00000.tab")
        final = _tab("sts_blob", "hydro_w")
        assert Path("bin/sts_blob.hydro_w.00001.bin").exists()
        assert np.max(final["s_00"]) < np.max(start["s_00"])
    finally:
        _cleanup()


def test_mhd_sts_independent_scalar_diffusivities():
    """Exercise MHD scalar STS tasks with separate diffusion coefficients."""
    try:
        _run("sts_mhd_scalar_advection.athinput", "sts_mhd_scalars")
        _run(
            "sts_mhd_scalar_advection.athinput",
            "explicit_mhd_scalars",
            "mhd/scalar_diffusivity_integrator=explicit",
            "time/sts_integrator=none",
        )
        sts = _tab("sts_mhd_scalars", "mhd_w")
        explicit = _tab("explicit_mhd_scalars", "mhd_w")
        assert np.max(np.abs(sts["s_00"] - explicit["s_00"])) < 1.0e-7
        assert np.max(np.abs(sts["s_01"] - explicit["s_01"])) < 1.0e-7
        amplitude_0 = np.max(sts["s_00"]) - np.min(sts["s_00"])
        amplitude_1 = np.max(sts["s_01"]) - np.min(sts["s_01"])
        assert amplitude_1 < amplitude_0
    finally:
        _cleanup()


@pytest.mark.parametrize("saturated", ("false", "true"))
def test_spitzer_requires_units_and_applies_heat_flux(tmp_path, saturated):
    source = Path(f"{INPUT_ROOT}/sts_thermal_front.athinput").read_text()
    source = source.replace("conductivity_model = power_law",
                            f"conductivity_model = spitzer\nsat_hflux = {saturated}")
    source = source.replace("conductivity_integrator = sts",
                            "conductivity_integrator = explicit")
    source = source.replace("sts_integrator = rkl2", "sts_integrator = none")
    source = source.replace("nx1 = 256", "nx1 = 32").replace("nx1 = 128", "nx1 = 32")
    source = source.replace("nlim = -1", "nlim = 1")
    input_path = tmp_path / "spitzer.athinput"
    input_path.write_text(source)
    command = ["./athena", "-i", str(input_path), "-d", str(tmp_path)]
    missing = subprocess.run(command, capture_output=True, text=True, check=False)
    assert missing.returncode != 0
    assert "Spitzer conduction requires a <units> block" in missing.stdout

    input_path.write_text(source +
                          "\n<units>\nlength_cgs = 1\ntime_cgs = 1e-6\nmass_cgs = 1e-8\n")
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    paths = sorted((tmp_path / "tab").glob("*.tab"))
    assert len(paths) == 2
    initial = testutils.athena_read.tab(str(paths[0]))
    final = testutils.athena_read.tab(str(paths[-1]))
    assert np.all(np.isfinite(final["eint"]))
    assert np.min(final["eint"]) > 0.0
    assert np.max(np.abs(final["eint"] - initial["eint"])) > 1e-6
