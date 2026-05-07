#!/usr/bin/env python3
"""Analytic and optional runtime checks for the TRML cooling refactor.

The analytic checks mirror src/utils/trml_cooling.hpp.  Runtime smoke tests run only
when ATHENA_TRML_EXE and ATHENA_TRML_NOSHEAR_EXE point at built executables.
"""

from __future__ import annotations

import math
import os
import re
import subprocess
import struct
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AMU_CODE = 1.67262192369e-24


def number_density(params: dict[str, float], dens: float) -> float:
    return dens / max(params["mass_per_particle"], 1.0e-300)


def cooling_lambda(params: dict[str, float], temp: float) -> float:
    if temp <= 0.0 or params["T_peak"] <= 0.0 or params["lambda_peak"] <= 0.0:
        return 0.0

    x = math.log(temp / params["T_peak"])
    width = max(params["lambda_smooth_width"], 1.0e-12)
    if x < 0.0:
        ln_shape = params["lambda_slope_lo"] * (x + width * (1.0 - math.exp(x / width)))
    else:
        ln_shape = params["lambda_slope_hi"] * (
            x + width * (math.exp(-x / width) - 1.0)
        )
    return params["lambda_peak"] * math.exp(ln_shape)


def evaluate(params: dict[str, float], dens: float, temp: float) -> tuple[float, float, float]:
    if dens <= 0.0 or temp <= 0.0:
        return 0.0, 0.0, 0.0
    if temp < params["cool_T_min"] or temp > params["cool_T_max"]:
        return 0.0, 0.0, 0.0
    if params.get("use_dens_ceiling", 0.0) and dens > params["dens_ceiling"]:
        return 0.0, 0.0, 0.0

    n = number_density(params, dens)
    edot_cool = n * n * cooling_lambda(params, temp)
    heat_shape = (temp / params["heat_gamma_T_ref"]) ** params["heat_gamma_T_slope"]
    edot_heat = n * params["heat_gamma"] * heat_shape
    return edot_cool, edot_heat, edot_cool - edot_heat


def log_slope(params: dict[str, float], temp: float, h: float = 1.0e-4) -> float:
    lam_p = cooling_lambda(params, temp * math.exp(h))
    lam_m = cooling_lambda(params, temp * math.exp(-h))
    return (math.log(lam_p) - math.log(lam_m)) / (2.0 * h)


def cooling_time_at_peak(params: dict[str, float], rho_0: float, pgas_0: float,
                         gamma: float) -> float:
    dens_peak = pgas_0 / params["T_peak"]
    eint_peak = pgas_0 / (gamma - 1.0)
    n_peak = number_density(params, dens_peak)
    return eint_peak / (n_peak * n_peak * params["lambda_peak"])


def baseline_params() -> dict[str, float]:
    return {
        "T_peak": 0.1,
        "cool_T_min": 0.01,
        "cool_T_max": 1.0,
        "lambda_peak": 2.5e-4,
        "lambda_slope_lo": 2.0,
        "lambda_slope_hi": -3.0,
        "lambda_smooth_width": 0.05,
        "heat_gamma": 0.0,
        "heat_gamma_T_ref": 0.01,
        "heat_gamma_T_slope": 0.0,
        "mass_per_particle": 1.0,
        "use_dens_ceiling": 0.0,
        "dens_ceiling": 1.0e99,
    }


def test_curve_shape_and_gates() -> None:
    params = baseline_params()

    assert math.isclose(cooling_lambda(params, params["T_peak"]),
                        params["lambda_peak"], rel_tol=1.0e-14)
    assert abs(log_slope(params, params["T_peak"])) < 2.0e-3
    assert math.isclose(log_slope(params, params["T_peak"] * 1.0e-4),
                        params["lambda_slope_lo"], rel_tol=1.0e-3)
    assert math.isclose(log_slope(params, params["T_peak"] * 1.0e4),
                        params["lambda_slope_hi"], rel_tol=1.0e-3)

    params["heat_gamma"] = 3.0e-5
    assert evaluate(params, 1.0, params["cool_T_min"] * 0.5) == (0.0, 0.0, 0.0)
    assert evaluate(params, 1.0, params["cool_T_max"] * 2.0) == (0.0, 0.0, 0.0)


def test_flat_lambda_curve() -> None:
    params = baseline_params()
    params["lambda_slope_lo"] = 0.0
    params["lambda_slope_hi"] = 0.0

    for temp in (params["cool_T_min"], 0.5 * params["T_peak"], params["T_peak"],
                 2.0 * params["T_peak"], params["cool_T_max"]):
        assert math.isclose(cooling_lambda(params, temp), params["lambda_peak"],
                            rel_tol=1.0e-14)
        assert math.isclose(log_slope(params, temp), 0.0, abs_tol=1.0e-10)


def test_balance_cold_and_xi_derivation() -> None:
    params = baseline_params()
    rho_0 = 1.0
    pgas_0 = 1.0
    contrast = 100.0
    gamma = 5.0 / 3.0
    T_cold = pgas_0 / rho_0 / contrast
    rho_cold = rho_0 * contrast

    n_cold = number_density(params, rho_cold)
    heat_shape_cold = (T_cold / params["heat_gamma_T_ref"]) ** params["heat_gamma_T_slope"]
    params["heat_gamma"] = n_cold * cooling_lambda(params, T_cold) / heat_shape_cold
    _, _, edot_net = evaluate(params, rho_cold, T_cold)
    assert abs(edot_net) < 1.0e-14

    temp_dependent = baseline_params()
    temp_dependent["heat_gamma"] = 2.0
    temp_dependent["heat_gamma_T_ref"] = T_cold
    temp_dependent["heat_gamma_T_slope"] = 1.5
    _, heat_ref, _ = evaluate(temp_dependent, 1.0, T_cold)
    _, heat_high, _ = evaluate(temp_dependent, 1.0, 2.0*T_cold)
    assert math.isclose(heat_high/heat_ref, 2.0**1.5, rel_tol=1.0e-14)

    temp_dependent["heat_gamma_T_slope"] = -0.75
    n_cold = number_density(temp_dependent, rho_cold)
    heat_shape_cold = (
        T_cold / temp_dependent["heat_gamma_T_ref"]
    ) ** temp_dependent["heat_gamma_T_slope"]
    temp_dependent["heat_gamma"] = (
        n_cold * cooling_lambda(temp_dependent, T_cold) / heat_shape_cold
    )
    _, _, edot_net = evaluate(temp_dependent, rho_cold, T_cold)
    assert abs(edot_net) < 1.0e-14

    t_shear = 1.0 / 0.645
    xi = 100.0
    t_cool_0 = t_shear / xi
    dens_peak = pgas_0 / params["T_peak"]
    eint_peak = pgas_0 / (gamma - 1.0)
    n_peak = number_density(params, dens_peak)
    params["lambda_peak"] = eint_peak / (t_cool_0 * n_peak * n_peak)
    assert math.isclose(cooling_time_at_peak(params, rho_0, pgas_0, gamma),
                        t_cool_0, rel_tol=1.0e-14)


def run_command(cmd: list[str], cwd: Path) -> None:
    completed = subprocess.run(cmd, cwd=cwd, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, check=False)
    if completed.returncode != 0:
        print(completed.stdout)
        raise SystemExit(completed.returncode)


def minimal_input(problem_block: str) -> str:
    return f"""
<job>
basename = TRML_cooling_smoke

<mesh>
nghost = 4
nx1 = 8
x1min = -0.5
x1max = 0.5
ix1_bc = periodic
ox1_bc = periodic
nx2 = 8
x2min = -0.5
x2max = 0.5
ix2_bc = periodic
ox2_bc = periodic
nx3 = 8
x3min = -0.5
x3max = 0.5
ix3_bc = outflow
ox3_bc = user

<meshblock>
nx1 = 8
nx2 = 8
nx3 = 8

<time>
evolution = dynamic
integrator = rk3
cfl_number = 0.3
tlim = 1.0e-3
ndiag = 1000000

<hydro>
eos = ideal
reconstruct = plm
rsolver = hllc
gamma = 1.666666667
nscalars = 1
fofc = true
dfloor = 1.0e-4
pfloor = 1.0e-5

<units>
length_cgs = 1.0
mass_cgs = 1.0
time_cgs = 1.0
mu = 0.6

<problem>
rho_0 = 1.0
pgas_0 = 1.0
density_contrast = 100.0
z_interface = 0.0
smoothing_thickness = 0.01
T_floor = 0.1
t_cool_start = 0.0
dt_cutoff = 1.0e-4
cfl_cool = 0.1
dt_min_terminate = -1.0
user_srcs = true
user_dt = true
user_hist = true
ndiag = -1
analysis_on = false
use_frame_tracking = false
user_work_in_loop = false
{problem_block}

<output1>
file_type = log
dcycle = 1000000
"""


def minimal_particle_input(problem_block: str) -> str:
    return f"""
<job>
basename = TRML_particle_cooling_smoke

<mesh>
nghost = 4
nx1 = 8
x1min = -0.5
x1max = 0.5
ix1_bc = periodic
ox1_bc = periodic
nx2 = 8
x2min = -0.5
x2max = 0.5
ix2_bc = periodic
ox2_bc = periodic
nx3 = 8
x3min = -0.5
x3max = 0.5
ix3_bc = outflow
ox3_bc = user

<meshblock>
nx1 = 8
nx2 = 8
nx3 = 8

<time>
evolution = dynamic
integrator = rk3
cfl_number = 0.3
tlim = 1.0e-3
ndiag = 1000000

<hydro>
eos = ideal
reconstruct = plm
rsolver = hllc
gamma = 1.666666667
nscalars = 1
fofc = true
dfloor = 1.0e-4
pfloor = 1.0e-5

<particles>
particle_type = lagrangian_mc
pusher = lagrangian_mc
target_count = 64
uniform_by_volume = false
pos_init_seed = 280496
random_seed = 12345
min_radius = 0.0
inject_at_inflow = false
assign_tag = rank_order

<units>
length_cgs = 1.0
mass_cgs = 1.0
time_cgs = 1.0
mu = 0.6

<problem>
rho_0 = 1.0
pgas_0 = 1.0
density_contrast = 100.0
velocity = 0.645
z_interface = 0.0
smoothing_thickness = 0.01
T_floor = 0.1
t_cool_start = 0.0
dt_cutoff = 1.0e-4
cfl_cool = 0.1
dt_min_terminate = -1.0
user_srcs = true
user_dt = true
user_hist = true
ndiag = -1
analysis_on = false
use_frame_tracking = false
user_work_in_loop = false
{problem_block}

<output1>
file_type = log
dcycle = 1000000

<output2>
file_type = trk
id = trk
variable = prtcl_all
dt = 1.0e-3
nparticles = 16
"""


def read_last_track_record(
    run_dir: Path,
) -> tuple[float, tuple[str, ...], list[dict[str, float]]]:
    paths = sorted(run_dir.rglob("*.trk"))
    if not paths:
        raise AssertionError(f"no tracked-particle output found in {run_dir}")
    data = paths[-1].read_bytes()
    marker = b"# AthenaK tracked particle data"
    header_start = data.rfind(marker)
    if header_start < 0:
        raise AssertionError(f"no tracked-particle header found in {paths[-1]}")
    header_end = data.find(b"\n", header_start)
    header = data[header_start:header_end].decode("ascii")
    match = re.search(
        r"time=\s*([-+0-9.eE]+).*ntracked_prtcls=\s*(\d+).*nvalues=\s*(\d+).*"
        r"variables=\s*(.*)$",
        header,
    )
    if match is None:
        raise AssertionError(f"unrecognized tracked-particle header: {header}")
    time = float(match.group(1))
    ntrack = int(match.group(2))
    nvalues = int(match.group(3))
    variables = tuple(match.group(4).split())
    payload = data[header_end + 1:header_end + 1 + ntrack * nvalues * 4]
    values = struct.unpack(f"<{ntrack * nvalues}f", payload)
    rows = []
    for row_index in range(ntrack):
        offset = row_index * nvalues
        rows.append({name: values[offset + col] for col, name in enumerate(variables)})
    return time, variables, rows


def check_particle_cooling_output(run_dir: Path) -> None:
    time, variables, rows = read_last_track_record(run_dir)
    required = {"active", "rho", "temp", "edot_cool", "edot_heat", "edot_net"}
    missing = required.difference(variables)
    if missing:
        raise AssertionError(f"tracked-particle output missing columns: {sorted(missing)}")
    if time <= 0.0:
        raise AssertionError("tracked-particle smoke needs a post-start output")

    params = baseline_params()
    params.update({
        "T_peak": 0.0464,
        "cool_T_min": 0.0,
        "cool_T_max": 0.85,
        "lambda_slope_lo": 2.0,
        "lambda_slope_hi": -3.0,
        "lambda_smooth_width": 0.05,
        "heat_gamma": 0.0,
        "heat_gamma_T_ref": 0.01,
        "heat_gamma_T_slope": 0.0,
        "mass_per_particle": 0.6 * AMU_CODE,
    })
    t_cool_0 = (1.0 / 0.645) / 100.0
    dens_peak = 1.0 / params["T_peak"]
    eint_peak = 1.0 / (5.0 / 3.0 - 1.0)
    n_peak = number_density(params, dens_peak)
    params["lambda_peak"] = eint_peak / (t_cool_0 * n_peak * n_peak)

    checked = 0
    for row in rows:
        if int(round(row["active"])) != 1 or row["rho"] <= 0.0 or row["temp"] <= 0.0:
            continue
        expected = evaluate(params, row["rho"], row["temp"])
        got = (row["edot_cool"], row["edot_heat"], row["edot_net"])
        for got_value, expected_value in zip(got, expected):
            scale = max(abs(expected_value), 1.0e-30)
            if abs(got_value - expected_value) / scale > 5.0e-5:
                raise AssertionError(
                    f"particle cooling mismatch: got={got}, expected={expected}"
                )
        checked += 1
    if checked == 0:
        raise AssertionError("no active tracked particles were checked")


def run_optional_smoke_tests() -> None:
    trml_exe = os.environ.get("ATHENA_TRML_EXE")
    noshear_exe = os.environ.get("ATHENA_TRML_NOSHEAR_EXE")
    if not trml_exe or not noshear_exe:
        print(
            "Skipping runtime smoke tests; set ATHENA_TRML_EXE and "
            "ATHENA_TRML_NOSHEAR_EXE."
        )
        return

    with tempfile.TemporaryDirectory(prefix="trml_cooling_smoke_") as tmp:
        tmp_path = Path(tmp)
        for run_name in ("direct", "xi", "balance", "legacy", "particles"):
            (tmp_path / run_name).mkdir()

        direct_input = tmp_path / "athinput.direct_lambda"
        direct_input.write_text(minimal_input("""
velocity = 0.645
lambda_slope_lo = 2.0
lambda_slope_hi = -3.0
lambda_smooth_width = 0.05
T_peak_over_T_cold = 4.64
cool_T_min = 0.0
cool_T_max = 0.85
lambda_peak = 1.0e-50
heating_mode = off
"""), encoding="utf-8")
        run_command([trml_exe, "-i", str(direct_input), "-d",
                     str(tmp_path / "direct")], ROOT)

        xi_input = tmp_path / "athinput.xi_derived"
        xi_input.write_text(minimal_input("""
velocity = 0.645
lambda_slope_lo = 2.0
lambda_slope_hi = -3.0
lambda_smooth_width = 0.05
T_peak_over_T_cold = 4.64
cool_T_min = 0.0
cool_T_max = 0.85
xi = 100.0
heating_mode = off
"""), encoding="utf-8")
        run_command([trml_exe, "-i", str(xi_input), "-d",
                     str(tmp_path / "xi")], ROOT)

        balance_input = tmp_path / "athinput.balance_cold"
        balance_input.write_text(minimal_input("""
velocity = 0.0
lambda_slope_lo = 2.0
lambda_slope_hi = -3.0
lambda_smooth_width = 0.05
T_peak_over_T_cold = 4.64
cool_T_min = 0.01
cool_T_max = 0.85
heat_gamma_T_ref = 0.01
heat_gamma_T_slope = -0.75
t_cool_0 = 0.05
heating_mode = balance_cold
"""), encoding="utf-8")
        run_command([noshear_exe, "-i", str(balance_input), "-d",
                     str(tmp_path / "balance")], ROOT)

        legacy_input = tmp_path / "athinput.legacy_alias"
        legacy_input.write_text(minimal_input("""
velocity = 0.645
beta_lo = -2.0
beta_hi = 3.0
T_peak_over_T_cold = 4.64
T_cutoff_over_T_hot = 0.85
cooling_below_T_cold = false
heating_on = false
xi = 100.0
"""), encoding="utf-8")
        run_command([trml_exe, "-i", str(legacy_input), "-d",
                     str(tmp_path / "legacy")], ROOT)

        particle_input = tmp_path / "athinput.particle_cooling"
        particle_input.write_text(minimal_particle_input("""
lambda_slope_lo = 2.0
lambda_slope_hi = -3.0
lambda_smooth_width = 0.05
T_peak_over_T_cold = 4.64
cool_T_min = 0.0
cool_T_max = 0.85
xi = 100.0
heating_mode = off
"""), encoding="utf-8")
        particle_run_dir = tmp_path / "particles"
        run_command([trml_exe, "-i", str(particle_input), "-d",
                     str(particle_run_dir)], ROOT)
        check_particle_cooling_output(particle_run_dir)


def main() -> None:
    test_curve_shape_and_gates()
    test_flat_lambda_curve()
    test_balance_cold_and_xi_derivation()
    run_optional_smoke_tests()
    print("TRML cooling refactor checks passed.")


if __name__ == "__main__":
    main()
