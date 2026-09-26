"""Stage-1 current-limited diffusion: analytic limits and nonlinear evolution."""

import math
import os
from pathlib import Path
from functools import lru_cache

import numpy as np
import pytest

from test_suite.diffusion.test_sts_rkl2_cpu import (
    cells, flags_for, read_binary, read_state, run_case,
)


def run_model(directory, flags=(), mode="explicit", **kwargs):
    settings = [f"mhd/resistivity_integrator={mode}",
                "time/sts_integrator=" + ("rkl2" if mode == "sts" else "none")]
    return run_case(directory, "mhd", settings + list(flags),
                    input_file="current_limited.athinput", **kwargs)


def energy_budget(data):
    _, volume, state = cells(data)
    total = np.dot(volume, state["ener"])
    magnetic = np.dot(volume, sum(state[f"bcc{d}"]**2 for d in (1, 2, 3))) / 2
    kinetic = np.dot(volume, sum(state[f"mom{d}"]**2 for d in (1, 2, 3))
                     / state["dens"]) / 2
    return total, magnetic, kinetic


def exact_forcefree_amplitude(time, k2, q0=3.0, eta0=0.1, eta_max=1.0):
    """Integrate dq/dt=-k^2 eta(q)q analytically, then invert the low-q branch."""
    qs = 1 - math.sqrt(eta0 / eta_max)
    equilibrium = qs * (1 - math.sqrt(eta0 / eta_max))
    transition = math.log((q0 - equilibrium) / (qs - equilibrium)) / (k2 * eta_max)
    if time <= transition:
        return (equilibrium + (q0 - equilibrium) * math.exp(-k2 * eta_max * time)) / q0
    target = math.log(qs) - qs - k2 * eta0 * (time - transition)
    low, high = 0.0, qs
    for _ in range(80):
        mid = (low + high) / 2
        if math.log(mid) - mid > target:
            high = mid
        else:
            low = mid
    return (low + high) / (2 * q0)


@pytest.mark.parametrize("mode", ["explicit", "sts"])
def test_constant_matches_saved_baseline(tmp_path, mode):
    baseline = os.environ.get("ATHENA_BASELINE")
    if not baseline:
        pytest.skip("set ATHENA_BASELINE to the unmodified STS executable")
    assert Path(baseline).is_file()
    fluid, flags = flags_for("ohmic", mode)
    flags += ["time/nlim=8", "time/tlim=10"]
    run_case(tmp_path / "baseline", fluid, flags, binary=baseline)
    run_case(tmp_path / "current", fluid, flags)
    a, b = (read_state(tmp_path / name) for name in ("baseline", "current"))
    assert (a["time"], a["cycle"]) == (b["time"], b["cycle"])
    for field in a["var_names"]:
        np.testing.assert_array_equal(a["mb_data"][field], b["mb_data"][field])


@pytest.mark.parametrize("mode", ["explicit", "sts"])
def test_degenerate_limit_is_constant(tmp_path, mode):
    for model in ("constant", "current_limited"):
        run_model(tmp_path / model, [f"mhd/resistivity_model={model}",
                  "mhd/eta_max=0.1", "time/nlim=8", "time/tlim=10"], mode=mode)
    a, b = (read_state(tmp_path / name) for name in ("constant", "current_limited"))
    assert (a["time"], a["cycle"]) == (b["time"], b["cycle"])
    for field in a["var_names"]:
        np.testing.assert_array_equal(a["mb_data"][field], b["mb_data"][field])


@pytest.mark.parametrize("mode", ["explicit", "sts"])
def test_weak_gaussian_matches_background_diffusion(tmp_path, mode):
    flags = ["problem/test=gaussian", "problem/b0=1e-6", "problem/width=0.1",
             "mesh/nx1=128", "time/tlim=0.02"]
    run_model(tmp_path, flags, mode=mode)
    data = read_state(tmp_path)
    xyz, _, state = cells(data)
    # The pgen initializes B_y = b0 exp[-x^2/(2 width^2)].
    width2 = 0.1**2 + 2 * 0.1 * data["time"]
    expected = 1e-6 * 0.1 / math.sqrt(width2) * np.exp(-xyz[:, 0]**2 / (2 * width2))
    assert np.mean(np.abs(state["bcc2"] - expected)) / 1e-6 < 8e-4


@pytest.mark.parametrize("guide", [0, 1])
def test_harris_initial_pressure_balance(tmp_path, guide):
    run_model(tmp_path, ["problem/test=harris", f"problem/guide_field={guide}",
              "mesh/nx1=64", "mesh/nx2=64", "meshblock/nx2=16", "time/nlim=0"])
    _, _, state = cells(read_state(tmp_path))
    magnetic = sum(state[f"bcc{d}"]**2 for d in (1, 2, 3)) / 2
    kinetic = sum(state[f"mom{d}"]**2 for d in (1, 2, 3)) / (2 * state["dens"])
    internal = state["ener"] - magnetic - kinetic
    assert internal.min() > 0
    np.testing.assert_allclose((2/3)*internal + magnetic, 1.5 + 0.5*guide**2,
                               rtol=2e-14, atol=2e-14)


@pytest.mark.parametrize("mode", ["explicit", "sts"])
@pytest.mark.parametrize("diagonal", [False, True])
def test_nonlinear_forcefree_spatial_convergence(tmp_path, mode, diagonal):
    errors, kinetic = [], []
    for nx in (16, 32) if diagonal else (32, 64):
        directory = tmp_path / str(nx)
        k2 = (2 * math.pi)**2 * (2 if diagonal else 1)
        flags = [f"mesh/nx1={nx}", f"mhd/d_i={3 / math.sqrt(k2)}",
                 f"time/tlim={0.12 / (2 if diagonal else 1)}"]
        if diagonal:
            flags += [f"mesh/nx2={nx}", "meshblock/nx2=16", "problem/wave_n2=1"]
        run_model(directory, flags, mode=mode)
        initial, final = read_state(directory, first=True), read_state(directory)
        _, vol, start = cells(initial)
        _, _, end = cells(final)
        # Project on the actual staggered initial field; its O(dx^2) amplitude
        # correction must converge with the discretization, not be fit away.
        amplitude = sum(np.dot(vol, end[f"bcc{d}"] * start[f"bcc{d}"])
                        for d in (1, 2, 3)) / sum(
            np.dot(vol, start[f"bcc{d}"]**2) for d in (1, 2, 3))
        expected = exact_forcefree_amplitude(final["time"], k2)
        errors.append(abs(amplitude - expected))
        e0, m0, k0 = energy_budget(initial)
        e1, m1, k1 = energy_budget(final)
        np.testing.assert_allclose(e1, e0, rtol=2e-12, atol=2e-12)
        assert m1 < 0.2 * m0
        assert e1 - m1 - k1 > e0 - m0 - k0
        kinetic.append(k1)
    assert errors[1] < 0.5 * errors[0], errors
    assert errors[1] < 0.01, errors
    # A rotated staggered field is only discretely force-free to truncation error.
    assert kinetic[1] < 0.6 * kinetic[0] + 1e-24, kinetic


def test_sts_temporal_convergence(tmp_path):
    flags = ["time/tlim=0.08", "time/cfl_number=0.02"]
    run_model(tmp_path / "reference", flags)
    _, _, reference = cells(read_state(tmp_path / "reference"))
    errors, cycles = [], []
    for ratio in (8, 4, 2):
        directory = tmp_path / str(ratio)
        run_model(directory, ["time/tlim=0.08", f"time/sts_max_dt_ratio={ratio}"],
                  mode="sts")
        data = read_state(directory)
        _, _, state = cells(data)
        errors.append(max(np.max(abs(state[f"bcc{d}"] - reference[f"bcc{d}"]))
                          for d in (1, 2, 3)))
        cycles.append(data["cycle"])
    assert cycles[0] < cycles[1] < cycles[2], cycles
    assert errors[1] < 0.6 * errors[0], errors
    assert errors[2] < 0.6 * errors[1], errors


@pytest.mark.parametrize("three_d", [False, True])
def test_nonlinear_operator_at_explicit_bound(tmp_path, three_d):
    flags = ["mesh/nx1=16", "mesh/nx2=16", "meshblock/nx2=16",
             "problem/wave_n2=1", "time/cfl_number=1", "time/nlim=32",
             "time/tlim=10", "mhd/d_i=0.1"]
    if three_d:
        flags += ["mesh/nx3=8", "meshblock/nx3=8", "problem/wave_n3=1"]
    run_model(tmp_path, flags)
    initial, final = read_state(tmp_path, first=True), read_state(tmp_path)
    _, _, state = cells(final)
    assert all(np.isfinite(field).all() for field in state.values())
    assert state["dens"].min() > 0.9
    e0, m0, _ = energy_budget(initial)
    e1, m1, _ = energy_budget(final)
    assert m1 < m0
    np.testing.assert_allclose(e1, e0, rtol=2e-12, atol=2e-12)


@lru_cache(maxsize=2)
def sheet_reference(time, nx=2048):
    """Independent conservative RK2 solution of dB/dt=d[eta(dB/dx)dB/dx]/dx."""
    dx = 1 / nx
    x = np.linspace(-0.5 + dx / 2, 0.5 - dx / 2, nx)
    field = 0.02 * np.tanh(x / 0.03)
    eta0, eta_max, di, brec = 0.01, 0.1, 0.09, 0.02
    qs, etas = 1 - math.sqrt(eta0 / eta_max), math.sqrt(eta0 * eta_max)

    def rhs(b):
        gradient = np.diff(np.pad(b, 1, mode="edge")) / dx
        q = np.abs(gradient) * di / brec
        eta = np.empty_like(q)
        low = q <= qs
        eta[low] = eta0 / (1 - q[low])
        eta[~low] = eta_max - qs / q[~low] * (eta_max - etas)
        return np.diff(eta * gradient) / dx

    steps = math.ceil(time / (0.4 * dx**2 / eta_max))
    dt = time / steps
    for _ in range(steps):
        first = field + dt * rhs(field)
        field = 0.5 * (field + first + dt * rhs(first))
    return x, field


@pytest.mark.parametrize("mode", ["explicit", "sts"])
def test_nonlinear_sheet_spreading(tmp_path, mode):
    flags = ["problem/test=sheet", "problem/b0=0.02", "problem/width=0.03",
             "mhd/b_rec=0.02", "mhd/d_i=0.09", "mhd/eta_max=0.1",
             "mhd/ohmic_resistivity=0.01", "time/tlim=0.1", "mesh/nx1=256",
             "mesh/ix1_bc=outflow", "mesh/ox1_bc=outflow", "output1/dt=0.005"]
    run_model(tmp_path, flags, mode=mode)
    data = read_state(tmp_path)
    xyz, _, state = cells(data)
    rx, rb = sheet_reference(data["time"])
    reference = np.interp(xyz[:, 0], rx, rb)
    assert np.mean(abs(state["bcc2"] - reference)) / 0.02 < 2e-3
    # High beta makes the static-density scalar equation an accurate limit.
    assert np.max(abs(state["dens"] - 1)) < 1e-3
    files = sorted((tmp_path / "bin").glob("*.state.*.bin"))
    widths, times = [], []
    for sample in (read_binary(files[0]), read_binary(files[1]), data):
        x, _, fields = cells(sample)
        rho = 0.5 * (fields["dens"][1:] + fields["dens"][:-1])
        q = np.max(abs(np.diff(fields["bcc2"]) / np.diff(x[:, 0])) / np.sqrt(rho))
        q *= 0.09 / 0.02
        widths.append(0.09 / q)
        times.append(sample["time"])
    assert 0.09 / widths[-1] < 1
    early_growth = (widths[1] - widths[0]) / (times[1] - times[0])
    late_growth = (widths[2] - widths[1]) / (times[2] - times[1])
    assert 0 < late_growth < 0.5 * early_growth


def test_sheet_reference_is_resolved():
    x, field = sheet_reference(0.1, 2048)
    xlo, fieldlo = sheet_reference(0.1, 1024)
    assert np.mean(abs(np.interp(xlo, x, field) - fieldlo)) / 0.02 < 2e-5


@pytest.mark.parametrize("setting, message", [
    ("mhd/ohmic_resistivity=0", "finite 0 < ohmic_resistivity <= eta_max"),
    ("mhd/eta_max=0.01", "finite 0 < ohmic_resistivity <= eta_max"),
    ("mhd/eta_max=nan", "finite 0 < ohmic_resistivity <= eta_max"),
    ("mhd/d_i=-1", "finite positive d_i and b_rec"),
    ("mhd/b_rec=inf", "finite positive d_i and b_rec"),
    ("mhd/b_rec_method=invalid", "b_rec_method must be constant or jump"),
    ("mhd/eta_ad=1", "cannot be combined with eta_ad"),
])
def test_invalid_current_limited_inputs(tmp_path, setting, message):
    result = run_model(tmp_path, [setting], check=False)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
