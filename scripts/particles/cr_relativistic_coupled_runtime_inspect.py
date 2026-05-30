#!/usr/bin/env python3
"""Qualify serial experimental mhd_ideal frozen_tn relativistic CR coupling."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from pathlib import Path


FIELDS = (
    "cycle time gid tag x y z vx vy vz wx wy wz cE1 cE2 cE3 "
    "b1 b2 b3 work alpha gamma kinetic dx dy dz db legacy_ipm "
    "final_live_u1 final_live_u2 final_live_u3 final_live_b1 final_live_b2 "
    "final_live_b3 final_live_cE1 final_live_cE2 final_live_cE3 "
    "final_live_cE_dot_b"
).split()
COMPARISON_FIELDS = (
    "cycle time x y z vx vy vz wx wy wz cE1 cE2 cE3 "
    "b1 b2 b3 work alpha gamma kinetic dx dy dz db"
).split()
WITNESS_ROWS: list[dict[str, float]] = []


def _add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def _sub(left, right):
    return tuple(a - b for a, b in zip(left, right))


def _scale(factor, vector):
    return tuple(factor * value for value in vector)


def _dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def _cross(left, right):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _norm(values):
    return math.sqrt(sum(value * value for value in values))


def _gamma(w, c_model=1.0):
    return math.sqrt(1.0 + _dot(w, w) / c_model**2)


def _velocity(w, c_model=1.0):
    gamma = _gamma(w, c_model)
    return _scale(1.0 / gamma, w)


def _kinetic(w, c_model=1.0):
    return c_model**2 * (_gamma(w, c_model) - 1.0)


def _ideal_ce(u, b):
    return _scale(-1.0, _cross(u, b))


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


def _hc_step(w, c_e, b, alpha, c_model, dt):
    """Independent Higuera-Cary transcription with drift-kick-drift displacement."""
    velocity_old = _velocity(w, c_model)
    half_scale = 0.5 * dt * alpha
    epsilon = _scale(half_scale, c_e)
    beta = _scale(half_scale, b)
    w_minus = _add(w, epsilon)
    gamma_minus = _gamma(w_minus, c_model)
    beta2 = _dot(beta, beta)
    sigma = gamma_minus**2 - beta2
    beta_dot_u = _dot(beta, w_minus) / c_model
    positive = beta2 + beta_dot_u**2
    radical = math.hypot(sigma, 2.0 * math.sqrt(positive))
    gamma_rotation2 = (
        0.5 * (sigma + radical) if sigma >= 0.0
        else 2.0 * positive / (radical - sigma)
    )
    t = _scale(1.0 / math.sqrt(gamma_rotation2), beta)
    s = _scale(2.0 / (1.0 + _dot(t, t)), t)
    w_prime = _add(w_minus, _cross(w_minus, t))
    w_plus = _add(w_minus, _cross(w_prime, s))
    w_new = _add(w_plus, epsilon)
    velocity_new = _velocity(w_new, c_model)
    work = alpha * dt * _dot(c_e, _add(w, w_new)) / (
        _gamma(w, c_model) + _gamma(w_new, c_model)
    )
    displacement = _scale(0.5 * dt, _add(velocity_old, velocity_new))
    return w_new, displacement, work


def _flags(**values):
    return [f"{key.replace('__', '/')}={value}" for key, value in values.items()]


def _read_witness(path):
    lines = path.read_text().splitlines()
    if len(lines) != 3 or lines[0] != (
            "# coupled-runtime\tmhd_ideal\tfrozen_tn\tcE=-u_cross_B"):
        raise RuntimeError(f"{path}: malformed coupled-runtime witness")
    reader = csv.DictReader([lines[1][2:], lines[2]], delimiter="\t")
    rows = list(reader)
    if len(rows) != 1:
        raise RuntimeError(f"{path}: expected one witness row, found {len(rows)}")
    return {key: float(value) for key, value in rows[0].items()}


def _read_prepush_witness(path):
    lines = path.read_text().splitlines()
    if len(lines) != 2:
        raise RuntimeError(f"{path}: malformed prepush witness")
    reader = csv.DictReader([lines[0][2:], lines[1]], delimiter="\t")
    rows = list(reader)
    if len(rows) != 1:
        raise RuntimeError(f"{path}: expected one prepush row, found {len(rows)}")
    return {key: float(value) for key, value in rows[0].items()}


def _run(binary, input_path, root, name, overrides):
    case_dir = root / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    command = [str(binary), "-i", str(input_path), *overrides]
    result = subprocess.run(command, cwd=case_dir, text=True, capture_output=True)
    (case_dir / "athena.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"{name} failed:\n{result.stdout}\n{result.stderr}")
    row = _read_witness(case_dir / "cr_relativistic_coupled_runtime.tsv")
    prepush_path = case_dir / "cr_relativistic_coupled_prepush.tsv"
    row["_prepush"] = (
        _read_prepush_witness(prepush_path) if prepush_path.exists() else None)
    gamma = _gamma((row["wx"], row["wy"], row["wz"]))
    row["_gamma_reconstruction_error"] = _scalar_near(
        row["gamma"], gamma, 3.0e-13, f"{name} independently reconstructed gamma")
    row["_velocity_shadow_error"] = _near(
        (row["vx"], row["vy"], row["vz"]),
        tuple(row[field] / gamma for field in ("wx", "wy", "wz")),
        3.0e-13, f"{name} velocity shadow")
    WITNESS_ROWS.append(row)
    return row


def _run_expected_failure(binary, input_path, root, name, overrides, expected):
    case_dir = root / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    command = [str(binary), "-i", str(input_path), *overrides]
    result = subprocess.run(command, cwd=case_dir, text=True, capture_output=True)
    output = result.stdout + result.stderr
    (case_dir / "athena.log").write_text(output)
    if result.returncode == 0 or expected not in output:
        raise RuntimeError(
            f"{name}: expected failure containing {expected!r}, got:\n{output}")


def _w(row):
    return row["wx"], row["wy"], row["wz"]


def _dx(row):
    return row["dx"], row["dy"], row["dz"]


def _ce(row):
    return row["cE1"], row["cE2"], row["cE3"]


def _b(row):
    return row["b1"], row["b2"], row["b3"]


def _fit_slope(step_sizes, errors):
    pairs = [(math.log(h), math.log(error)) for h, error in zip(step_sizes, errors)
             if error > 0.0]
    mean_x = sum(pair[0] for pair in pairs) / len(pairs)
    mean_y = sum(pair[1] for pair in pairs) / len(pairs)
    return sum((x - mean_x) * (y - mean_y) for x, y in pairs) / sum(
        (x - mean_x) ** 2 for x, _ in pairs)


def _rk4_manufactured(w0, x0, u0, acceleration, b, alpha, c_model, tlim,
                      nsteps=200000):
    state = tuple(w0) + tuple(x0)
    dt = tlim / nsteps

    def rhs(time, values):
        w = values[:3]
        velocity = _velocity(w, c_model)
        fluid = (u0[0], u0[1] + acceleration*time, u0[2])
        force = _scale(alpha, _add(_ideal_ce(fluid, b), _cross(velocity, b)))
        return force + velocity

    time = 0.0
    for _ in range(nsteps):
        k1 = rhs(time, state)
        k2 = rhs(time + 0.5*dt, _add(state, _scale(0.5*dt, k1)))
        k3 = rhs(time + 0.5*dt, _add(state, _scale(0.5*dt, k2)))
        k4 = rhs(time + dt, _add(state, _scale(dt, k3)))
        weighted = _add(_add(k1, _scale(2.0, k2)),
                        _add(_scale(2.0, k3), k4))
        state = _add(state, _scale(dt/6.0, weighted))
        time += dt
    return state


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
    parser.add_argument("--prescribed-input", required=True, type=Path)
    parser.add_argument("--criteria", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--metrics", type=Path)
    parser.add_argument("--manufactured-oracle-acceleration-scale", type=float,
                        default=1.0)
    parser.add_argument("--manufactured-one-step-oracle-time-level",
                        choices=("tn", "right_end"), default="tn")
    args = parser.parse_args()
    binary = args.binary.resolve()
    input_path = args.input.resolve()
    prescribed_input = args.prescribed_input.resolve()
    metrics = {}
    WITNESS_ROWS.clear()

    one_step = _flags(time__tlim=0.05, time__nlim=1, time__cfl_number=100.0,
                      particles__cfl_part=100.0)
    live_u = (0.0, -0.2, 0.0)
    live_b = (0.0, 0.0, 0.4)
    expected_ce = _ideal_ce(live_u, live_b)
    uniform = _run(binary, input_path, args.work_dir, "uniform_live", one_step)
    expected_w, expected_dx, _ = _hc_step(
        (0.0, 0.0, 0.0), expected_ce, live_b, 1.0, 1.0, uniform["time"])
    metrics["uniform_live_field_momentum_error"] = _near(
        _w(uniform), expected_w, 2.0e-13, "uniform live-field momentum")
    metrics["uniform_live_field_trajectory_error"] = _near(
        _dx(uniform), expected_dx, 2.0e-13, "uniform live-field displacement")
    metrics["ideal_ce_sign_error"] = _near(
        _ce(uniform), expected_ce, 2.0e-13, "pointwise ideal cE sign")
    metrics["ideal_ce_dot_b_error"] = abs(_dot(_ce(uniform), _b(uniform)))

    poison = _run(binary, input_path, args.work_dir, "poison_control", one_step + _flags(
        problem__poison_B1=-31.0, problem__poison_B2=37.0, problem__poison_B3=-41.0,
        problem__poison_cE1=43.0, problem__poison_cE2=-47.0,
        problem__poison_cE3=53.0))
    metrics["stale_field_poison_overwrite_max_error"] = max(
        abs(uniform[field] - poison[field]) for field in COMPARISON_FIELDS)
    prepush = uniform["_prepush"]
    prepush_tuple = tuple(
        prepush[field] for field in ("b1", "b2", "b3", "cE1", "cE2", "cE3"))
    postpush_tuple = _b(uniform) + _ce(uniform)
    metrics["runtime_prepush_poison_error"] = _near(
        prepush_tuple, (7.0, -11.0, 13.0, -17.0, 19.0, -23.0), 0.0,
        "runtime prepush poison tuple")
    poison_prepush = poison["_prepush"]
    poison_prepush_tuple = tuple(
        poison_prepush[field]
        for field in ("b1", "b2", "b3", "cE1", "cE2", "cE3"))
    metrics["alternate_runtime_prepush_poison_error"] = _near(
        poison_prepush_tuple, (-31.0, 37.0, -41.0, 43.0, -47.0, 53.0), 0.0,
        "alternate runtime prepush poison tuple")
    metrics["stale_field_poison_sentinel_separation"] = _norm(
        _sub(prepush_tuple, postpush_tuple))

    cancel_u = (0.2, -0.1, 0.05)
    cancel_b = (0.11, -0.07, 0.4)
    cancel_gamma = 1.0 / math.sqrt(1.0 - _dot(cancel_u, cancel_u))
    cancel_w = _scale(cancel_gamma, cancel_u)
    cancel = _run(binary, input_path, args.work_dir, "fluid_velocity_cancellation",
                  one_step + _flags(
                      problem__w0x=cancel_w[0], problem__w0y=cancel_w[1],
                      problem__w0z=cancel_w[2], problem__live_u1=cancel_u[0],
                      problem__live_u2=cancel_u[1], problem__live_u3=cancel_u[2],
                      problem__live_B1=cancel_b[0], problem__live_B2=cancel_b[1],
                      problem__live_B3=cancel_b[2], problem__B0x=cancel_b[0],
                      problem__B0y=cancel_b[1], problem__B0z=cancel_b[2]))
    metrics["fluid_velocity_cancellation_momentum_error"] = _near(
        _w(cancel), cancel_w, 2.0e-13, "fluid-velocity force cancellation momentum")
    metrics["fluid_velocity_cancellation_trajectory_error"] = _near(
        _dx(cancel), _scale(cancel["time"], cancel_u), 2.0e-13,
        "fluid-velocity force cancellation displacement")

    parity_w = (0.23, -0.11, 0.04)
    parity_b = (0.0, 0.0, 0.4)
    parity_ce = (0.0, 0.05, 0.0)
    coupled_parity = _run(
        binary, input_path, args.work_dir, "coupled_parity", one_step + _flags(
            problem__w0x=parity_w[0], problem__w0y=parity_w[1],
            problem__w0z=parity_w[2], problem__live_u1=0.125,
            problem__live_u2=0.0, problem__live_u3=0.0,
            problem__live_B1=parity_b[0], problem__live_B2=parity_b[1],
            problem__live_B3=parity_b[2]))
    prescribed_parity = _run(
        binary, prescribed_input, args.work_dir, "prescribed_parity", one_step + _flags(
            problem__w0x=parity_w[0], problem__w0y=parity_w[1],
            problem__w0z=parity_w[2], problem__cE0x=parity_ce[0],
            problem__cE0y=parity_ce[1], problem__cE0z=parity_ce[2],
            problem__B0x=parity_b[0], problem__B0y=parity_b[1],
            problem__B0z=parity_b[2]))
    metrics["prescribed_coupled_parity_max_error"] = max(
        abs(coupled_parity[field] - prescribed_parity[field])
        for field in COMPARISON_FIELDS)

    metrics["acceleration_work_closure_error"] = _scalar_near(
        uniform["work"], uniform["kinetic"], 2.0e-13,
        "accelerating coupled direct-work closure")
    decel_w0 = (0.2, 0.0, 0.0)
    decel = _run(binary, input_path, args.work_dir, "deceleration", one_step + _flags(
        problem__w0x=decel_w0[0], problem__live_u1=0.0,
        problem__live_u2=0.25, problem__live_u3=0.0))
    metrics["deceleration_work_closure_error"] = _scalar_near(
        decel["work"], decel["kinetic"] - _kinetic(decel_w0), 2.0e-13,
        "decelerating coupled direct-work closure")
    metrics["deceleration_work"] = decel["work"]

    midpoint_w0 = (0.2, 0.0, 0.0)
    midpoint = _run(binary, input_path, args.work_dir, "midpoint_sampling", _flags(
        time__tlim=0.1, time__nlim=1, time__cfl_number=100.0,
        particles__cfl_part=100.0, problem__w0x=midpoint_w0[0],
        problem__live_velocity_profile="affine", problem__live_u1=0.0,
        problem__live_u2=0.0, problem__live_u3=0.0,
        problem__live_u2_x1=1.0))
    midpoint_x = 0.5 * midpoint["time"] * _velocity(midpoint_w0)[0]
    midpoint_ce = _ideal_ce((0.0, midpoint_x, 0.0), live_b)
    midpoint_expected_w, _, _ = _hc_step(
        midpoint_w0, midpoint_ce, live_b, 1.0, 1.0, midpoint["time"])
    metrics["midpoint_live_sampling_momentum_error"] = _near(
        _w(midpoint), midpoint_expected_w, 2.0e-13, "spatial-midpoint coupled momentum")
    metrics["midpoint_live_sampling_field_error"] = _near(
        _ce(midpoint), midpoint_ce, 2.0e-13, "spatial-midpoint coupled field")
    metrics["midpoint_old_position_separation"] = _norm(
        _sub(midpoint_ce, _ideal_ce((0.0, 0.0, 0.0), live_b)))

    temporal_common = _flags(
        time__tlim=0.15, time__nlim=1000, particles__cfl_part=100.0,
        problem__particle_x1=0.17, problem__particle_x2=-0.13,
        problem__particle_x3=0.11, problem__w0x=0.21, problem__w0y=-0.09,
        problem__w0z=0.04, problem__live_velocity_profile="periodic",
        problem__live_u1=0.04, problem__live_u2=-0.03, problem__live_u3=0.02,
        problem__live_velocity_amplitude=0.08, problem__live_velocity_wave_number=1.0,
        problem__live_B1=0.08, problem__live_B2=-0.05, problem__live_B3=0.35,
        problem__B0x=0.08, problem__B0y=-0.05, problem__B0z=0.35)
    temporal_ref = _run(binary, input_path, args.work_dir, "temporal_reference",
                        temporal_common + _flags(time__cfl_number=0.00375))
    step_sizes = []
    temporal_errors = []
    temporal_states = []
    reference_state = _w(temporal_ref) + (temporal_ref["x"], temporal_ref["y"],
                                          temporal_ref["z"])
    for n, cfl in enumerate((0.06, 0.03, 0.015, 0.0075)):
        row = _run(binary, input_path, args.work_dir, f"temporal_{n}",
                   temporal_common + _flags(time__cfl_number=cfl))
        state = _w(row) + (row["x"], row["y"], row["z"])
        temporal_errors.append(_norm(_sub(state, reference_state)))
        temporal_states.append(state)
        step_sizes.append(cfl)
    self_convergence_errors = [
        _norm(_sub(left, right))
        for left, right in zip(temporal_states, temporal_states[1:] + [reference_state])
    ]
    metrics["frozen_tn_temporal_convergence_slope"] = _fit_slope(
        step_sizes, self_convergence_errors)
    metrics["frozen_tn_temporal_finest_error"] = temporal_errors[-1]

    manufactured_u0 = (0.04, -0.03, 0.02)
    manufactured_acceleration = 0.2
    manufactured_b = (0.08, -0.05, 0.35)
    manufactured_w0 = (0.21, -0.09, 0.04)
    manufactured_x0 = (0.17, -0.13, 0.11)
    manufactured_common = _flags(
        time__tlim=0.15, time__nlim=1000, particles__cfl_part=100.0,
        problem__pgen_name="cr_relativistic_coupled_runtime_test",
        problem__user_srcs="true", problem__relativistic_coupled_manufactured_test="true",
        problem__manufactured_acceleration=manufactured_acceleration,
        problem__particle_x1=manufactured_x0[0], problem__particle_x2=manufactured_x0[1],
        problem__particle_x3=manufactured_x0[2], problem__w0x=manufactured_w0[0],
        problem__w0y=manufactured_w0[1], problem__w0z=manufactured_w0[2],
        problem__live_u1=manufactured_u0[0], problem__live_u2=manufactured_u0[1],
        problem__live_u3=manufactured_u0[2], problem__live_B1=manufactured_b[0],
        problem__live_B2=manufactured_b[1], problem__live_B3=manufactured_b[2],
        problem__B0x=manufactured_b[0], problem__B0y=manufactured_b[1],
        problem__B0z=manufactured_b[2])
    manufactured_one_step = _run(
        binary, input_path, args.work_dir, "manufactured_one_step_orientation",
        manufactured_common + _flags(time__tlim=0.05, time__nlim=1,
                                     time__cfl_number=100.0))
    manufactured_tn_ce = _ideal_ce(manufactured_u0, manufactured_b)
    manufactured_right_end_u = (
        manufactured_u0[0],
        manufactured_u0[1] + manufactured_acceleration*manufactured_one_step["time"],
        manufactured_u0[2])
    manufactured_right_end_ce = _ideal_ce(manufactured_right_end_u, manufactured_b)
    manufactured_oracle_ce = (
        manufactured_tn_ce
        if args.manufactured_one_step_oracle_time_level == "tn"
        else manufactured_right_end_ce)
    manufactured_one_step_expected_w, _, _ = _hc_step(
        manufactured_w0, manufactured_oracle_ce, manufactured_b, 1.0, 1.0,
        manufactured_one_step["time"])
    metrics["manufactured_one_step_tn_sampled_ce_error"] = _near(
        _ce(manufactured_one_step), manufactured_oracle_ce, 2.0e-13,
        "manufactured one-step frozen-t^n sampled cE")
    metrics["manufactured_one_step_tn_momentum_error"] = _near(
        _w(manufactured_one_step), manufactured_one_step_expected_w, 2.0e-13,
        "manufactured one-step frozen-t^n momentum")
    metrics["manufactured_one_step_live_right_end_ce_error"] = _near(
        tuple(manufactured_one_step[field] for field in
              ("final_live_cE1", "final_live_cE2", "final_live_cE3")),
        manufactured_right_end_ce, 2.0e-13,
        "manufactured one-step evolved live right-end cE")
    metrics["manufactured_one_step_tn_vs_right_end_separation"] = _norm(
        _sub(manufactured_tn_ce, manufactured_right_end_ce))
    manufactured_reference = _rk4_manufactured(
        manufactured_w0, manufactured_x0, manufactured_u0,
        manufactured_acceleration*args.manufactured_oracle_acceleration_scale,
        manufactured_b, 1.0, 1.0, 0.15)
    manufactured_steps = []
    manufactured_errors = []
    manufactured_rows = []
    for n, cfl in enumerate((0.06, 0.03, 0.015, 0.0075)):
        row = _run(binary, input_path, args.work_dir, f"manufactured_{n}",
                   manufactured_common + _flags(time__cfl_number=cfl))
        state = _w(row) + (row["x"], row["y"], row["z"])
        manufactured_steps.append(cfl)
        manufactured_errors.append(_norm(_sub(state, manufactured_reference)))
        manufactured_rows.append(row)
    metrics["manufactured_temporal_convergence_slope"] = _fit_slope(
        manufactured_steps, manufactured_errors)
    metrics["manufactured_temporal_finest_error"] = manufactured_errors[-1]
    initial_ce = _ideal_ce(manufactured_u0, manufactured_b)
    metrics["manufactured_live_state_evolution_sentinel"] = _norm(_sub(
        tuple(manufactured_rows[-1][field] for field in
              ("final_live_cE1", "final_live_cE2", "final_live_cE3")),
        initial_ce))
    stale_t0_reference = _rk4_manufactured(
        manufactured_w0, manufactured_x0, manufactured_u0, 0.0, manufactured_b,
        1.0, 1.0, 0.15)
    metrics["manufactured_stale_t0_discriminator"] = _norm(_sub(
        _w(manufactured_rows[-1]), stale_t0_reference[:3]))
    metrics["manufactured_multicycle_work_closure_error"] = _scalar_near(
        manufactured_rows[-1]["work"],
        manufactured_rows[-1]["kinetic"] - _kinetic(manufactured_w0), 2.0e-12,
        "manufactured multi-cycle direct-work closure")

    _run_expected_failure(
        binary, input_path, args.work_dir, "failure_invalid_gid",
        one_step + _flags(problem__runtime_invalid_gid="true"), "status 207")
    _run_expected_failure(
        binary, input_path, args.work_dir, "failure_superluminal_fluid",
        one_step + _flags(problem__live_u1=1.0), "status 201")
    _run_expected_failure(
        binary, input_path, args.work_dir, "failure_stencil_guard",
        one_step + _flags(problem__runtime_invalid_stencil_position="true"),
        "status 208")
    _run_expected_failure(
        binary, input_path, args.work_dir, "failure_extreme_stencil_guard",
        one_step + _flags(problem__runtime_invalid_stencil_position="true",
                          problem__runtime_invalid_stencil_x1=1e300),
        "status 208")
    _run_expected_failure(
        binary, input_path, args.work_dir, "failure_negative_extreme_stencil_guard",
        one_step + _flags(problem__runtime_invalid_stencil_position="true",
                          problem__runtime_invalid_stencil_x1=-1e300),
        "status 208")
    metrics["production_failure_probe_count"] = 5
    metrics["final_live_ce_dot_b_error"] = max(
        abs(row["final_live_cE_dot_b"]) for row in WITNESS_ROWS)

    metrics["production_velocity_shadow_max_error"] = max(
        row["_velocity_shadow_error"] for row in WITNESS_ROWS)
    metrics["production_gamma_reconstruction_max_error"] = max(
        row["_gamma_reconstruction_error"] for row in WITNESS_ROWS)
    _validate_criteria(args.criteria, metrics)
    if args.metrics:
        args.metrics.parent.mkdir(parents=True, exist_ok=True)
        args.metrics.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    print("CR relativistic coupled frozen_tn runtime checks passed")


if __name__ == "__main__":
    main()
