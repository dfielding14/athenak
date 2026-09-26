"""Check explicit diffusion bounds against uniform modes and a coefficient jump."""

import os
from pathlib import Path
import re
import subprocess

import pytest


INPUT = Path(__file__).resolve().parents[2] / "inputs/diffusion_timestep.athinput"
CFL = 0.4
GM1 = 2.0 / 3.0
# These units make T_cgs = T_code and the low-temperature conductivity sqrt(T_code).
UNITS = f"\n<units>\nmu = {1.3806488e-16 / 1.660538921e-24:.17g}\nmass_cgs = 2500\n"


def run_case(tmp_path, fluid, *flags, extra=""):
    deck = tmp_path / "test.athinput"
    deck.write_text(INPUT.read_text().replace("<hydro>", f"<{fluid}>") + extra)
    binary = Path(os.environ.get("ATHENA", "./athena")).resolve()
    return subprocess.run(
        [str(binary), "-i", str(deck), *flags], cwd=tmp_path,
        capture_output=True, text=True, check=False, timeout=30,
    )


def timesteps(result):
    assert result.returncode == 0, result.stdout + result.stderr
    values = re.findall(r"cycle=\d+ .*dt=([0-9.eE+-]+)", result.stdout)
    assert values, result.stdout
    return [float(value) for value in values]


@pytest.mark.parametrize("fluid", ["hydro", "mhd"])
@pytest.mark.parametrize("shape", [(32, 1, 1), (32, 8, 1), (32, 8, 4)])
@pytest.mark.parametrize("coefficient", ["conductivity", "viscosity"])
def test_uniform_diffusion_bound(tmp_path, fluid, shape, coefficient):
    flags = [f"{fluid}/{coefficient}=1"]
    for axis in (2, 3):
        flags += [f"mesh/nx{axis}={shape[axis-1]}",
                  f"meshblock/nx{axis}={shape[axis-1]}"]
    inverse_spacing = [n*n for n in shape if n > 1]
    if coefficient == "conductivity":
        rate = 2 * GM1 * sum(inverse_spacing)
    else:
        rate = 2 * (sum(inverse_spacing) + max(inverse_spacing) / 3)
    assert timesteps(run_case(tmp_path, fluid, *flags))[0] == pytest.approx(
        CFL / rate, rel=2e-6
    )


@pytest.mark.parametrize("fluid", ["hydro", "mhd"])
@pytest.mark.parametrize("saturated", ["false", "true"])
@pytest.mark.parametrize("axis", [1, 2, 3])
def test_temperature_jump_uses_face_coefficients(tmp_path, fluid, saturated, axis):
    # Hot side: rho=100, T=100, kappa=10. Cold side: rho=T=kappa=1.
    # The cold cell next to the interface limits dt with face coefficients 5.5 and 1.
    shape = (32, 8, 4)
    result = run_case(
        tmp_path, fluid, f"{fluid}/tdep_conductivity=true",
        f"{fluid}/sat_hflux={saturated}", "problem/dl=100", "problem/pl=10000",
        f"problem/shock_dir={axis}", "mesh/nx2=8", "meshblock/nx2=8",
        "mesh/nx3=4", "meshblock/nx3=4", extra=UNITS,
    )
    row_sum = 2 * sum(n*n for n in shape) + 4.5 * shape[axis-1]**2
    assert timesteps(result)[0] == pytest.approx(CFL / (GM1 * row_sum), rel=2e-6)


@pytest.mark.parametrize("fluid", ["hydro", "mhd"])
def test_viscosity_bound_refreshes_after_amr(tmp_path, fluid):
    extra = f"""
<mesh_refinement>
refinement = adaptive
num_levels = 2
refinement_interval = 1
max_nmb_per_rank = 16
<amr_criterion0>
method = min_max
variable = {fluid}_w_d
value_max = 0.5
"""
    result = run_case(tmp_path, fluid, f"{fluid}/viscosity=1", "time/nlim=2", extra=extra)
    dt = timesteps(result)
    assert len(dt) >= 2
    assert dt[1] == pytest.approx(dt[0] / 4, rel=2e-6)


def test_temperature_dependent_conduction_requires_units(tmp_path):
    result = run_case(tmp_path, "hydro", "hydro/tdep_conductivity=true")
    assert result.returncode != 0
    assert "Temperature-dependent conduction requires a <units> block" in result.stdout
