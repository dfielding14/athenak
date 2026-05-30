#!/usr/bin/env python3
"""Exercise the serial prescribed-field CR production path against independent oracles."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from pathlib import Path


FIELDS = ("cycle time tag x y z vx vy vz wx wy wz cE1 cE2 cE3 "
          "b1 b2 b3 work alpha gamma kinetic dx dy dz db legacy_ipm").split()
PHYSICS_FIELDS = tuple(field for field in FIELDS if field != "legacy_ipm")
WITNESS_ROWS = []


def _norm(values):
    return math.sqrt(sum(value * value for value in values))


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _vector(row, prefix):
    return tuple(row[f"{prefix}{axis}"] for axis in "xyz")


def _near(actual, expected, tolerance, label):
    error = _norm(_sub(actual, expected))
    scale = max(1.0, _norm(expected))
    if error > tolerance * scale:
        raise RuntimeError(f"{label}: error={error:.6e}, tolerance={tolerance:.6e}")
    return error


def _scalar_near(actual, expected, tolerance, label):
    error = abs(actual - expected)
    scale = max(1.0, abs(expected))
    if error > tolerance * scale:
        raise RuntimeError(f"{label}: error={error:.6e}, tolerance={tolerance:.6e}")
    return error


def _read_witness(path):
    with path.open() as stream:
        rows = list(csv.DictReader(
            (line for line in stream if not line.startswith("#")),
            fieldnames=FIELDS, delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"expected one witness row, found {len(rows)}")
    return {key: float(value) for key, value in rows[0].items()}


def _run(binary, input_path, root, name, overrides, c_model=1.0):
    case_dir = root / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    command = [str(binary), "-i", str(input_path), *overrides]
    result = subprocess.run(command, cwd=case_dir, text=True, capture_output=True)
    (case_dir / "athena.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"{name} failed:\n{result.stdout}\n{result.stderr}")
    row = _read_witness(case_dir / "cr_relativistic_pusher_runtime.tsv")
    gamma = math.sqrt(
        1.0 + sum(row[field] ** 2 for field in ("wx", "wy", "wz")) / c_model**2)
    row["_gamma_reconstruction_error"] = _scalar_near(
        row["gamma"], gamma, 3.0e-13, f"{name} independently reconstructed gamma")
    row["_velocity_shadow_error"] = _near(
        (row["vx"], row["vy"], row["vz"]),
        tuple(row[field] / gamma for field in ("wx", "wy", "wz")),
        3.0e-13, f"{name} velocity shadow")
    WITNESS_ROWS.append(row)
    return row


def _run_expect_failure(binary, input_path, root, name, overrides, expected):
    case_dir = root / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    command = [str(binary), "-i", str(input_path), *overrides]
    result = subprocess.run(command, cwd=case_dir, text=True, capture_output=True)
    output = result.stdout + result.stderr
    (case_dir / "athena.log").write_text(output)
    if result.returncode == 0 or expected.lower() not in output.lower():
        raise RuntimeError(
            f"{name} did not fail closed with {expected!r}:\n{output}")


def _rk4(w0, x0, c_e, b, alpha, c_model, total_time, steps=120000):
    """Independent high-accuracy reference for dw/dt and dx/dt."""
    dt = total_time / steps

    def rhs(state):
        w = state[:3]
        gamma = math.sqrt(1.0 + sum(value * value for value in w) / c_model**2)
        v = tuple(value / gamma for value in w)
        vx_b = (v[1] * b[2] - v[2] * b[1],
                v[2] * b[0] - v[0] * b[2],
                v[0] * b[1] - v[1] * b[0])
        return tuple(alpha * (c_e[n] + vx_b[n]) for n in range(3)) + v

    state = tuple(w0) + tuple(x0)
    for _ in range(steps):
        k1 = rhs(state)
        k2 = rhs(tuple(value + 0.5 * dt * slope
                       for value, slope in zip(state, k1)))
        k3 = rhs(tuple(value + 0.5 * dt * slope
                       for value, slope in zip(state, k2)))
        k4 = rhs(tuple(value + dt * slope for value, slope in zip(state, k3)))
        state = tuple(value + dt * (a + 2.0 * b2 + 2.0 * c + d) / 6.0
                      for value, a, b2, c, d in zip(state, k1, k2, k3, k4))
    return state[:3], state[3:]


def _hc_magnetic_step(w, b, alpha, c_model, dt):
    """Independent one-step HC transcription for the negative-sigma production probe."""
    beta = tuple(0.5 * dt * alpha * value for value in b)
    gamma_minus = math.sqrt(1.0 + sum(value * value for value in w) / c_model**2)
    beta2 = sum(value * value for value in beta)
    sigma = gamma_minus**2 - beta2
    beta_dot_u = sum(left * right for left, right in zip(beta, w)) / c_model
    positive = beta2 + beta_dot_u**2
    radical = math.sqrt(sigma**2 + 4.0 * positive)
    gamma_rotation2 = (
        0.5 * (sigma + radical) if sigma >= 0.0
        else 2.0 * positive / (radical - sigma))
    t = tuple(value / math.sqrt(gamma_rotation2) for value in beta)
    t2 = sum(value * value for value in t)
    s = tuple(2.0 * value / (1.0 + t2) for value in t)

    def cross(a, right):
        return (a[1] * right[2] - a[2] * right[1],
                a[2] * right[0] - a[0] * right[2],
                a[0] * right[1] - a[1] * right[0])

    w_cross_t = cross(w, t)
    w_prime = tuple(left + right for left, right in zip(w, w_cross_t))
    w_new = tuple(left + right for left, right in zip(w, cross(w_prime, s)))
    gamma_new = math.sqrt(1.0 + sum(value * value for value in w_new) / c_model**2)
    velocity_old = tuple(value / gamma_minus for value in w)
    velocity_new = tuple(value / gamma_new for value in w_new)
    displacement = tuple(
        0.5 * dt * (left + right) for left, right in zip(velocity_old, velocity_new))
    return w_new, displacement, sigma


def _w(row):
    return row["wx"], row["wy"], row["wz"]


def _x(row):
    return row["x"], row["y"], row["z"]


def _dx(row):
    return row["dx"], row["dy"], row["dz"]


def _flags(**values):
    return [f"{key.replace('__', '/')}={value}" for key, value in values.items()]


def _fit_slope(step_sizes, errors):
    pairs = [(math.log(h), math.log(error)) for h, error in zip(step_sizes, errors)
             if error > 0.0]
    mean_x = sum(pair[0] for pair in pairs) / len(pairs)
    mean_y = sum(pair[1] for pair in pairs) / len(pairs)
    return sum((x - mean_x) * (y - mean_y) for x, y in pairs) / sum(
        (x - mean_x) ** 2 for x, _ in pairs)


def _validate_criteria(path, metrics):
    criteria = json.loads(path.read_text())
    registered = {item["metric"] for item in criteria["criteria"]}
    emitted = set(metrics)
    if emitted != registered:
        raise RuntimeError(
            f"criteria completeness mismatch: unbound={sorted(emitted - registered)}, "
            f"missing={sorted(registered - emitted)}")
    for item in criteria["criteria"]:
        actual = metrics[item["metric"]]
        limit = item["limit"]
        operator = item["operator"]
        accepted = ((operator == "<=" and actual <= limit) or
                    (operator == ">=" and actual >= limit) or
                    (operator == "==" and actual == limit))
        if not accepted:
            raise RuntimeError(
                f"{item['id']} failed: {item['metric']}={actual} "
                f"{operator} {limit} is false")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--metrics", type=Path)
    parser.add_argument("--criteria", required=True, type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    input_path = args.input.resolve()
    metrics = {}
    WITNESS_ROWS.clear()

    gamma0 = 1.0 / math.sqrt(1.0 - (0.2**2 + 0.1**2 + 0.05**2))
    w0 = (gamma0 * 0.2, -gamma0 * 0.1, gamma0 * 0.05)
    zero = _run(binary, input_path, args.work_dir, "zero", _flags(
        time__nlim=4, problem__w0x=w0[0], problem__w0y=w0[1], problem__w0z=w0[2]))
    metrics["zero_field_momentum_error"] = _near(
        _w(zero), w0, 2.0e-13, "zero-field momentum")
    metrics["zero_field_trajectory_error"] = _near(
        _dx(zero), tuple(zero["time"] * v for v in (0.2, -0.1, 0.05)),
        2.0e-13, "zero-field displacement")

    electric = _run(binary, input_path, args.work_dir, "electric", _flags(
        time__nlim=4, problem__cE0x=0.3))
    expected_w = (0.3 * electric["time"], 0.0, 0.0)
    gamma = math.sqrt(1.0 + expected_w[0] ** 2)
    expected_x = ((gamma - 1.0) / 0.3, 0.0, 0.0)
    metrics["pure_e_momentum_error"] = _near(
        _w(electric), expected_w, 2.0e-13, "pure-E momentum")
    metrics["pure_e_trajectory_error"] = _near(
        _dx(electric), expected_x, 4.0e-4, "pure-E trajectory")
    metrics["pure_e_work_closure_error"] = _scalar_near(
        electric["work"], electric["kinetic"], 2.0e-13, "pure-E direct work")

    electric_c2 = _run(binary, input_path, args.work_dir, "electric_cmodel_2", _flags(
        time__nlim=4, particles__c_model=2.0, problem__cE0x=0.3), c_model=2.0)
    expected_w_c2 = (0.3 * electric_c2["time"], 0.0, 0.0)
    expected_gamma_c2 = math.sqrt(1.0 + expected_w_c2[0] ** 2 / 4.0)
    expected_x_c2 = (4.0 * (expected_gamma_c2 - 1.0) / 0.3, 0.0, 0.0)
    metrics["nonunit_cmodel_momentum_error"] = _near(
        _w(electric_c2), expected_w_c2, 2.0e-13, "non-unit-c_model pure-E momentum")
    metrics["nonunit_cmodel_trajectory_error"] = _near(
        _dx(electric_c2), expected_x_c2, 4.0e-4, "non-unit-c_model pure-E trajectory")
    metrics["nonunit_cmodel_work_closure_error"] = _scalar_near(
        electric_c2["work"], electric_c2["kinetic"], 2.0e-13,
        "non-unit-c_model pure-E direct work")

    decel = _run(binary, input_path, args.work_dir, "deceleration", _flags(
        time__nlim=1, problem__w0x=0.2, problem__cE0x=-0.1))
    if not decel["work"] < 0.0:
        raise RuntimeError("controlled deceleration did not record negative direct work")
    metrics["deceleration_work"] = decel["work"]

    gamma_b = 1.0 / math.sqrt(1.0 - 0.2**2)
    w_b = gamma_b * 0.2
    magnetic = _run(binary, input_path, args.work_dir, "magnetic", _flags(
        time__nlim=8, problem__w0x=w_b, problem__B0z=0.4))
    metrics["magnetic_energy_error"] = _scalar_near(
        magnetic["kinetic"], gamma_b - 1.0, 3.0e-13, "B-only kinetic energy")
    reference_w, reference_x = _rk4(
        (w_b, 0.0, 0.0), (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0), (0.0, 0.0, 0.4), 1.0, 1.0, magnetic["time"])
    metrics["magnetic_rk4_momentum_error"] = _near(
        _w(magnetic), reference_w, 8.0e-5, "B-only RK4 momentum")
    metrics["magnetic_rk4_trajectory_error"] = _near(
        _dx(magnetic), reference_x, 8.0e-5, "B-only RK4 displacement")
    metrics["magnetic_phase_error"] = abs(
        math.atan2(magnetic["wy"], magnetic["wx"]) + 0.4 * magnetic["time"] / gamma_b)
    if metrics["magnetic_phase_error"] > 8.0e-5:
        raise RuntimeError("B-only phase oracle failed")

    magnetic_c2 = _run(binary, input_path, args.work_dir, "magnetic_cmodel_2", _flags(
        time__nlim=8, particles__c_model=2.0, problem__w0x=0.4,
        problem__B0z=0.4), c_model=2.0)
    reference_w, reference_x = _rk4(
        (0.4, 0.0, 0.0), (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0), (0.0, 0.0, 0.4), 1.0, 2.0, magnetic_c2["time"])
    metrics["nonunit_cmodel_magnetic_rk4_momentum_error"] = _near(
        _w(magnetic_c2), reference_w, 8.0e-5, "non-unit-c_model B-only RK4 momentum")
    metrics["nonunit_cmodel_magnetic_rk4_trajectory_error"] = _near(
        _dx(magnetic_c2), reference_x, 8.0e-5,
        "non-unit-c_model B-only RK4 displacement")

    equilibrium = _run(binary, input_path, args.work_dir, "equilibrium", _flags(
        time__nlim=4, problem__w0x=w_b, problem__B0z=0.4, problem__cE0y=0.08))
    metrics["equilibrium_momentum_error"] = _near(
        _w(equilibrium), (w_b, 0.0, 0.0), 2.0e-13,
        "crossed-field equilibrium momentum")
    metrics["equilibrium_trajectory_error"] = _near(
        _dx(equilibrium), (0.2 * equilibrium["time"], 0.0, 0.0), 2.0e-13,
        "crossed-field equilibrium displacement")

    crossed = _run(binary, input_path, args.work_dir, "crossed", _flags(
        time__nlim=8, problem__B0z=0.4, problem__cE0y=0.05))
    reference_w, reference_x = _rk4(
        (0.0, 0.0, 0.0), (0.0, 0.0, 0.0),
        (0.0, 0.05, 0.0), (0.0, 0.0, 0.4), 1.0, 1.0, crossed["time"])
    metrics["crossed_rk4_momentum_error"] = _near(
        _w(crossed), reference_w, 8.0e-5, "crossed-field RK4 momentum")
    metrics["crossed_rk4_trajectory_error"] = _near(
        _dx(crossed), reference_x, 8.0e-5, "crossed-field RK4 displacement")
    metrics["crossed_work_closure_error"] = _scalar_near(
        crossed["work"], crossed["kinetic"], 2.0e-13, "crossed-field direct work")

    negative = _run(binary, input_path, args.work_dir, "negative_alpha", _flags(
        time__nlim=4, particles__alpha_s=-1.0, problem__cE0x=0.3))
    metrics["negative_alpha_mirror_error"] = _near(
        _w(negative), tuple(-value for value in _w(electric)), 2.0e-13,
        "negative-alpha pure-E momentum")
    _scalar_near(negative["alpha"], -1.0, 0.0, "stored negative alpha")

    reverse_forward = _run(binary, input_path, args.work_dir, "reverse_forward", _flags(
        time__nlim=6, problem__w0x=0.31, problem__w0y=-0.17, problem__w0z=0.08,
        problem__B0x=0.11, problem__B0y=-0.07, problem__B0z=0.43))
    reverse_back = _run(binary, input_path, args.work_dir, "reverse_back", _flags(
        time__nlim=6, particles__alpha_s=-1.0,
        problem__w0x=reverse_forward["wx"], problem__w0y=reverse_forward["wy"],
        problem__w0z=reverse_forward["wz"], problem__B0x=0.11,
        problem__B0y=-0.07, problem__B0z=0.43))
    metrics["magnetic_reversibility_momentum_error"] = _near(
        _w(reverse_back), (0.31, -0.17, 0.08), 2.0e-12,
        "production magnetic reversibility")

    high = _run(binary, input_path, args.work_dir, "high_gamma", _flags(
        time__nlim=8, problem__w0x=1000.0, problem__B0z=0.4))
    if not math.isfinite(high["gamma"]):
        raise RuntimeError("high-gamma production run produced nonfinite gamma")
    metrics["high_gamma_relative_energy_drift"] = abs(
        high["kinetic"] - (math.sqrt(1.0 + 1000.0**2) - 1.0)) / high["kinetic"]
    if metrics["high_gamma_relative_energy_drift"] > 2.0e-13:
        raise RuntimeError("high-gamma B-only energy drift exceeded tolerance")
    metrics["high_gamma_phase_error"] = abs(
        math.atan2(high["wy"], high["wx"]) +
        0.4 * high["time"] / math.sqrt(1.0 + 1000.0**2))
    if metrics["high_gamma_phase_error"] > 2.0e-10:
        raise RuntimeError("high-gamma B-only phase oracle failed")

    conjugate = _run(binary, input_path, args.work_dir, "conjugate_root", _flags(
        time__tlim=0.25, time__nlim=1, time__cfl_number=100.0,
        particles__cfl_part=100.0, problem__w0x=0.2, problem__B0z=12.0))
    expected_w, expected_dx, sigma = _hc_magnetic_step(
        (0.2, 0.0, 0.0), (0.0, 0.0, 12.0), 1.0, 1.0, conjugate["time"])
    if not sigma < 0.0:
        raise RuntimeError("conjugate-root production oracle did not select sigma < 0")
    metrics["conjugate_root_sigma"] = sigma
    metrics["conjugate_root_momentum_error"] = _near(
        _w(conjugate), expected_w, 2.0e-13, "conjugate-root production momentum")
    metrics["conjugate_root_trajectory_error"] = _near(
        _dx(conjugate), expected_dx, 2.0e-13, "conjugate-root production displacement")

    ipm_a = _run(binary, input_path, args.work_dir, "legacy_ipm_a", _flags(
        time__tlim=0.005, time__nlim=100, problem__w0x=0.2, problem__B0z=0.4,
        particles__min_mass=1.0))
    ipm_b = _run(binary, input_path, args.work_dir, "legacy_ipm_b", _flags(
        time__tlim=0.005, time__nlim=100, problem__w0x=0.2, problem__B0z=0.4,
        particles__min_mass=99.0))
    if ipm_a["legacy_ipm"] == ipm_b["legacy_ipm"]:
        raise RuntimeError("legacy-IPM poison control did not vary IPM")
    metrics["legacy_ipm_independence_max_error"] = max(
        abs(ipm_a[field] - ipm_b[field]) for field in PHYSICS_FIELDS)
    if metrics["legacy_ipm_independence_max_error"] > 2.0e-13:
        raise RuntimeError("relativistic production path depends on legacy IPM")

    live_a = _run(binary, input_path, args.work_dir, "live_grid_a", _flags(
        time__tlim=0.005, time__nlim=100, problem__w0x=0.2, problem__B0z=0.4))
    live_b = _run(binary, input_path, args.work_dir, "live_grid_b", _flags(
        time__tlim=0.005, time__nlim=100, problem__w0x=0.2, problem__B0z=0.4,
        problem__runtime_poison_live_mhd="true"))
    metrics["live_grid_isolation_max_error"] = max(
        abs(live_a[field] - live_b[field]) for field in PHYSICS_FIELDS)
    if metrics["live_grid_isolation_max_error"] > 2.0e-13:
        raise RuntimeError(
            "Phase-4a prescribed push depends on poisoned live-grid arrays")

    expected_failures = [
        ("failure_initial_shadow_cap",
         _flags(time__nlim=0, problem__w0x=1.0e20),
         "could not initialize a finite representable authoritative w state"),
        ("failure_electric_range",
         _flags(time__nlim=1, particles__alpha_s=1.0e308, problem__cE0x=1.0e308),
         "status 103"),
        ("failure_magnetic_range",
         _flags(time__nlim=1, particles__alpha_s=1.0e308, problem__B0z=1.0),
         "status 103"),
        ("failure_output_shadow_cap",
         _flags(time__nlim=1, problem__cE0x=1.0e12),
         "status 105"),
    ]
    for name, overrides, expected in expected_failures:
        _run_expect_failure(binary, input_path, args.work_dir, name, overrides, expected)
    metrics["production_failure_probe_count"] = len(expected_failures)

    step_sizes = []
    convergence_errors = []
    convergence_momentum_errors = []
    convergence_trajectory_errors = []
    reference_w, reference_x = _rk4(
        (0.23, -0.11, 0.04), (0.0, 0.0, 0.0),
        (0.0, 0.05, 0.0), (0.0, 0.0, 0.4), 1.0, 1.0, 0.2)
    for n, cfl in enumerate((0.08, 0.04, 0.02, 0.01)):
        row = _run(binary, input_path, args.work_dir, f"convergence_{n}", _flags(
            time__tlim=0.2, time__nlim=100, time__cfl_number=cfl,
            particles__cfl_part=cfl,
            problem__w0x=0.23, problem__w0y=-0.11, problem__w0z=0.04,
            problem__cE0y=0.05, problem__B0z=0.4))
        momentum_error = _norm(_sub(_w(row), reference_w))
        trajectory_error = _norm(_sub(_dx(row), reference_x))
        error = momentum_error + trajectory_error
        convergence_errors.append(error)
        convergence_momentum_errors.append(momentum_error)
        convergence_trajectory_errors.append(trajectory_error)
        step_sizes.append(cfl)
    metrics["production_convergence_slope"] = _fit_slope(step_sizes, convergence_errors)
    metrics["production_convergence_finest_error"] = convergence_errors[-1]
    metrics["production_convergence_momentum_slope"] = _fit_slope(
        step_sizes, convergence_momentum_errors)
    metrics["production_convergence_trajectory_slope"] = _fit_slope(
        step_sizes, convergence_trajectory_errors)
    metrics["production_convergence_finest_momentum_error"] = (
        convergence_momentum_errors[-1])
    metrics["production_convergence_finest_trajectory_error"] = (
        convergence_trajectory_errors[-1])
    if metrics["production_convergence_slope"] < 1.8:
        raise RuntimeError(
            "production-path convergence slope is below second-order target")

    metrics["production_velocity_shadow_max_error"] = max(
        row["_velocity_shadow_error"] for row in WITNESS_ROWS)
    metrics["production_gamma_reconstruction_max_error"] = max(
        row["_gamma_reconstruction_error"] for row in WITNESS_ROWS)
    _validate_criteria(args.criteria, metrics)
    if args.metrics:
        args.metrics.parent.mkdir(parents=True, exist_ok=True)
        args.metrics.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    print("CR relativistic prescribed production-path runtime checks passed")


if __name__ == "__main__":
    main()
