"""CPU regression coverage for fourth-derivative hyperviscosity and RKL2 STS."""

import shutil
import subprocess

import numpy as np

import test_suite.testutils as testutils


INPUT_ROOT = "../../../inputs/tests"


def _run(input_name, basename, *flags):
    testutils.run(
        f"{INPUT_ROOT}/{input_name}",
        [f"job/basename={basename}", *flags],
    )


def _errors(basename):
    return np.loadtxt(f"{basename}-hypervisc-errors.dat")


def _tab(basename, variable):
    return testutils.athena_read.tab(f"tab/{basename}.{variable}.00001.tab")


def _failure(input_name, *flags):
    command = ["./athena", "-i", f"{INPUT_ROOT}/{input_name}", *flags]
    return subprocess.run(command, capture_output=True, text=True, check=False)


def _cleanup():
    testutils.cleanup()
    shutil.rmtree("bin", ignore_errors=True)


def test_hydro_hyperviscosity_convergence_and_energy_transfer():
    """The analytic shear mode converges and conservatively heats an ideal gas."""
    try:
        errors = {}
        for resolution in (32, 64, 128):
            basename = f"hypervisc_hydro_{resolution}"
            _run(
                "sts_hyperviscous_shear.athinput",
                basename,
                f"mesh/nx1={resolution}",
                f"meshblock/nx1={min(64, resolution)}",
                "output2/dt=0.05",
            )
            errors[resolution] = _errors(basename)

        assert errors[128][2] < 2.0e-5
        assert errors[64][2] / errors[32][2] < 0.27
        assert errors[128][2] / errors[64][2] < 0.27

        history = testutils.athena_read.hst("hypervisc_hydro_128.hydro.hst")
        assert abs(history["tot-E"][-1] - history["tot-E"][0]) < 2.0e-12
        assert history["2-KE"][-1] < history["2-KE"][0]
    finally:
        _cleanup()


def test_hyperviscosity_explicit_sts_selectivity_and_speedup():
    """Capped STS matches explicit while uncapped STS avoids the dx^4 cycle limit."""
    try:
        common = (
            "mesh/nx1=64",
            "meshblock/nx1=64",
            "time/tlim=0.005",
            "output1/dt=0.005",
            "output2/dt=0.005",
        )
        _run(
            "sts_hyperviscous_shear.athinput",
            "hypervisc_capped",
            *common,
            "problem/mode=4.0",
            "time/sts_max_dt_ratio=1.0",
        )
        _run(
            "sts_hyperviscous_shear.athinput",
            "hypervisc_explicit",
            *common,
            "problem/mode=4.0",
            "hydro/hyperviscosity_integrator=explicit",
            "time/sts_integrator=none",
        )
        capped = _tab("hypervisc_capped", "hydro_w")
        explicit = _tab("hypervisc_explicit", "hydro_w")
        assert np.max(np.abs(capped["vely"] - explicit["vely"])) < 5.0e-8

        _run(
            "sts_hyperviscous_shear.athinput",
            "hypervisc_low",
            *common,
            "problem/mode=1.0",
        )
        _run(
            "sts_hyperviscous_shear.athinput",
            "hypervisc_high",
            *common,
            "problem/mode=4.0",
        )
        low = _errors("hypervisc_low")
        high = _errors("hypervisc_high")
        assert high[4] < low[4]
        assert high[1] < _errors("hypervisc_explicit")[1]

        expected_low = 0.1*np.exp(-1.0e-4*(2.0*np.pi)**4*0.005)
        expected_high = 0.1*np.exp(-1.0e-4*(8.0*np.pi)**4*0.005)
        assert abs(low[4] - expected_low) / expected_low < 1.0e-4
        assert abs(high[4] - expected_high) / expected_high < 1.0e-2
    finally:
        _cleanup()


def test_mhd_mixed_and_isothermal_hyperviscosity_paths():
    """MHD, mixed viscosity, and isothermal momentum-only paths execute correctly."""
    try:
        _run("sts_mhd_hyperviscous_shear.athinput", "hypervisc_mhd")
        mhd = _errors("hypervisc_mhd")
        assert mhd[2] < 1.0e-5

        _run("sts_viscosity_plus_hyperviscosity.athinput", "hypervisc_mixed")
        mixed = _errors("hypervisc_mixed")
        assert mixed[2] < 2.0e-5

        _run(
            "sts_hyperviscous_shear.athinput",
            "hypervisc_isothermal",
            "hydro/eos=isothermal",
            "hydro/iso_sound_speed=1.0",
            "hydro/rsolver=hlle",
        )
        isothermal = _errors("hypervisc_isothermal")
        assert np.all(np.isfinite(isothermal))
    finally:
        _cleanup()


def test_hyperviscosity_rejects_unsupported_configurations():
    """Reject coefficients and discretizations outside the supported contract."""
    try:
        negative = _failure(
            "sts_hyperviscous_shear.athinput",
            "hydro/hyperviscosity=-1.0e-4",
        )
        assert negative.returncode != 0
        assert "hyperviscosity must be non-negative" in negative.stdout + negative.stderr

        multilevel = _failure(
            "unsupported_hyperviscosity_smr.athinput",
        )
        assert multilevel.returncode != 0
        assert "not supported with SMR or AMR" in multilevel.stdout + multilevel.stderr

        relativistic = _failure(
            "sts_hyperviscous_shear.athinput",
            "coord/special_rel=true",
            "hydro/rsolver=hlle",
        )
        assert relativistic.returncode != 0
        assert "supported only for Newtonian" in relativistic.stdout + relativistic.stderr
    finally:
        _cleanup()
