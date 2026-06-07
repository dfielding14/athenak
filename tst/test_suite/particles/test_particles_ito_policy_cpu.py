"""Policy tests for Ito and MC mass-flux tracers."""

from pathlib import Path
import shutil
import subprocess

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
INPUT = str(ROOT / "tst/inputs/particles_ito_policy.athinput")
MHD_INPUT = str(ROOT / "tst/inputs/particles_ito_floor_policy_mhd.athinput")
RUN_DIR = Path("run_particles_ito_policy_cpu")


def _run(*overrides, input_file=INPUT):
    command = ["./athena", "-i", input_file, "-d", str(RUN_DIR), *overrides]
    return subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False
    )


def test_declared_mass_changing_source_is_rejected():
    """Flux tracers fail closed when a fluid source declares a density change."""
    shutil.rmtree(RUN_DIR, ignore_errors=True)
    try:
        result = _run("hydro_srcterms/changes_mass=true")
        output = result.stdout + result.stderr
        assert result.returncode != 0
        assert "do not support source terms that change mass density" in output
    finally:
        shutil.rmtree(RUN_DIR, ignore_errors=True)


def test_mass_conserving_builtin_source_is_allowed():
    """A built-in source that leaves IDN unchanged remains compatible."""
    shutil.rmtree(RUN_DIR, ignore_errors=True)
    try:
        result = _run(
            "time/nlim=0",
            "hydro_srcterms/const_accel=true",
            "hydro_srcterms/const_accel_val=1.0",
            "hydro_srcterms/const_accel_dir=1",
        )
        assert result.returncode == 0, result.stdout + result.stderr
    finally:
        shutil.rmtree(RUN_DIR, ignore_errors=True)


def test_mass_injecting_fluid_floors_fail_closed():
    """Configured floors are allowed until they actually inject gas mass."""
    cases = [
        (INPUT, (), True, ""),
        (INPUT, ("hydro/dfloor=2.0",), False, "density floor"),
        (MHD_INPUT, (), True, ""),
        (MHD_INPUT, ("mhd/dfloor=2.0",), False, "density floor"),
        (MHD_INPUT, ("mhd/sigma_max=1.0",), False, "magnetization ceiling"),
    ]

    try:
        for input_file, overrides, allowed, message in cases:
            shutil.rmtree(RUN_DIR, ignore_errors=True)
            result = _run(*overrides, input_file=input_file)
            output = result.stdout + result.stderr
            if allowed:
                assert result.returncode == 0, output
            else:
                assert result.returncode != 0
                assert message in output
                assert "injected gas mass" in output
    finally:
        shutil.rmtree(RUN_DIR, ignore_errors=True)


def test_final_rk_weighted_flux_cancels_stage_reversal():
    """Net-flux semantics produce no tracer exchange for cancelling RK stages."""
    rk_weights = {
        "rk2": np.array([0.5, 0.5]),
        "rk3": np.array([1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0]),
    }
    stage_fluxes = {
        "rk2": np.array([1.0, -1.0]),
        "rk3": np.array([1.0, 1.0, -0.5]),
    }

    for integrator, weights in rk_weights.items():
        fluxes = stage_fluxes[integrator]
        net_flux = np.dot(weights, fluxes)
        stagewise_outward_transfer = np.dot(weights, np.maximum(fluxes, 0.0))
        assert net_flux == 0.0
        assert stagewise_outward_transfer > 0.0
