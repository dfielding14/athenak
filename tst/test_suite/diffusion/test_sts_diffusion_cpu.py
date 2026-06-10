"""CPU regression coverage for RKL2 super time stepping."""

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

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


def _command(input_name, basename, *flags):
    command = [
        "./athena",
        "-i",
        f"{INPUT_ROOT}/{input_name}",
        f"job/basename={basename}",
        *flags,
    ]
    return subprocess.run(command, capture_output=True, text=True, check=False)


def _diagnostic_dt(result, cycle=1):
    matches = re.findall(
        rf"cycle={cycle}\s+time=[0-9.eE+-]+\s+dt=([0-9.eE+-]+)",
        result.stdout,
    )
    assert matches, result.stdout + result.stderr
    return float(matches[0])


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


def test_constant_viscosity_anisotropic_timestep_bound():
    """The constant-viscosity limit includes every active mesh direction."""
    try:
        cases = (
            (8, 1, 1, 2.0, 1.0, 1.0),
            (8, 4, 1, 2.0, 2.0, 1.0),
            (8, 4, 4, 2.0, 2.0, 4.0),
        )
        cfl = 0.4
        nu = 1.0
        for nx1, nx2, nx3, lx1, lx2, lx3 in cases:
            result = _command(
                "sts_viscous_shear.athinput",
                f"constant_viscosity_{nx1}_{nx2}_{nx3}",
                f"mesh/nx1={nx1}",
                f"mesh/nx2={nx2}",
                f"mesh/nx3={nx3}",
                "meshblock/nx1=8",
                f"meshblock/nx2={nx2}",
                f"meshblock/nx3={nx3}",
                f"mesh/x1max={lx1}",
                f"mesh/x2max={lx2}",
                f"mesh/x3max={lx3}",
                "time/cfl_number=0.4",
                "time/tlim=1.0",
                "time/nlim=1",
                "time/ndiag=1",
                "time/sts_integrator=none",
                "hydro/viscosity=1.0",
                "hydro/tdep_viscosity=false",
                "hydro/viscosity_integrator=explicit",
            )
            assert result.returncode == 0, result.stdout + result.stderr

            inverse_dx2 = [(nx1/lx1) ** 2]
            if nx2 > 1:
                inverse_dx2.append((nx2/lx2) ** 2)
            if nx3 > 1:
                inverse_dx2.append((nx3/lx3) ** 2)
            expected = cfl/(2.0*nu*(sum(inverse_dx2) + max(inverse_dx2)/3.0))
            assert np.isclose(_diagnostic_dt(result), expected, rtol=2.0e-6)
    finally:
        _cleanup()


def test_conduction_model_activation_and_saturated_timestep():
    """Model-only Spitzer conduction is active and retains its parabolic timestep."""
    try:
        _run(
            "conduction_configuration.athinput",
            "spitzer_model",
        )
        _run(
            "conduction_configuration.athinput",
            "no_conduction",
            "hydro/conductivity_model=constant",
        )
        _run(
            "conduction_configuration.athinput",
            "capped_spitzer",
            "hydro/cond_ceiling=0.01",
        )
        spitzer = _tab("spitzer_model", "hydro_w")
        no_conduction = _tab("no_conduction", "hydro_w")
        capped = _tab("capped_spitzer", "hydro_w")
        spitzer_change = np.max(np.abs(spitzer["eint"] - no_conduction["eint"]))
        capped_change = np.max(np.abs(capped["eint"] - no_conduction["eint"]))
        assert spitzer_change > 1.0e-4
        assert 0.0 < capped_change < spitzer_change

        saturated = _command(
            "conduction_configuration.athinput",
            "saturated_spitzer",
            "hydro/sat_hflux=true",
        )
        assert saturated.returncode == 0, saturated.stdout + saturated.stderr
        match = re.search(r"time=.* cycle=(\d+)", saturated.stdout)
        assert match is not None
        assert int(match.group(1)) > 1
    finally:
        _cleanup()


def test_conduction_rejects_invalid_configurations():
    """Reject unsupported saturation and invalid Spitzer unit or ceiling choices."""
    try:
        no_units = _command(
            "sts_thermal_front.athinput",
            "spitzer_no_units",
            "hydro/conductivity_model=spitzer",
            "hydro/conductivity_integrator=explicit",
            "time/sts_integrator=none",
        )
        assert no_units.returncode != 0
        assert "requires a <units> block" in no_units.stdout + no_units.stderr

        for model in ("constant", "power_law"):
            unsupported = _command(
                "conduction_configuration.athinput",
                f"saturated_{model}",
                f"hydro/conductivity_model={model}",
                "hydro/sat_hflux=true",
            )
            assert unsupported.returncode != 0
            assert "sat_hflux is implemented only" in unsupported.stdout + unsupported.stderr

        negative_ceiling = _command(
            "conduction_configuration.athinput",
            "negative_conduction_ceiling",
            "hydro/cond_ceiling=-1.0",
        )
        assert negative_ceiling.returncode != 0
        assert "Invalid conductivity" in negative_ceiling.stdout + negative_ceiling.stderr

        sts_saturation = _command(
            "conduction_configuration.athinput",
            "sts_saturated_spitzer",
            "hydro/sat_hflux=true",
            "hydro/conductivity_integrator=sts",
            "time/sts_integrator=rkl2",
        )
        assert sts_saturation.returncode != 0
        assert "sat_hflux is not compatible with STS" in (
            sts_saturation.stdout + sts_saturation.stderr
        )
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
