#!/usr/bin/env python3
"""Qualify acceleration-aware relativistic CR subcycling and timestep refresh."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import cr_relativistic_coupled_runtime_inspect as coupled


SCHEDULE_RE = re.compile(
    r"Particle relativistic subcycling: steps=(?P<steps>\d+) "
    r"active_constraint=(?P<constraint>[a-z_]+) "
    r"cell_steps=(?P<cell>\d+) meshblock_steps=(?P<meshblock>\d+) "
    r"gyro_steps=(?P<gyro>\d+) electric_steps=(?P<electric>\d+)")
CYCLE_RE = re.compile(
    r"elapsed=\S+ cycle=(?P<cycle>\d+) time=(?P<time>\S+) dt=(?P<dt>\S+)")


def _flags(**values):
    return coupled._flags(**values)


def _run(binary, input_path, root, name, overrides):
    row = coupled._run(binary, input_path, root, name, overrides)
    text = (root / name / "athena.log").read_text()
    schedules = [
        {key: int(value) if key != "constraint" else value
         for key, value in match.groupdict().items()}
        for match in SCHEDULE_RE.finditer(text)
    ]
    cycles = [
        {"cycle": int(match.group("cycle")), "time": float(match.group("time")),
         "dt": float(match.group("dt"))}
        for match in CYCLE_RE.finditer(text)
    ]
    return row, schedules, cycles


def _expected_failure(binary, input_path, root, name, overrides, expected):
    case_dir = root / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    command = [str(binary), "-i", str(input_path), *overrides]
    result = subprocess.run(
        command, cwd=case_dir, text=True, capture_output=True, timeout=10)
    output = result.stdout + result.stderr
    (case_dir / "athena.log").write_text(output)
    if result.returncode == 0 or expected not in output:
        raise RuntimeError(
            f"{name}: expected failure containing {expected!r}, got:\n{output}")


def _particle_input_variant(input_path, root, name, line):
    path = root / f"{name}.athinput"
    text = input_path.read_text()
    marker = "subcycle_electric_kick_max     = 0.1\n"
    if marker not in text:
        raise RuntimeError("could not locate Phase-5 particle-input insertion point")
    path.write_text(text.replace(marker, marker + line + "\n"))
    return path


def _only_schedule(schedules, name):
    if len(schedules) != 1:
        raise RuntimeError(f"{name}: expected one schedule row, found {schedules}")
    return schedules[0]


def _repeat_hc(w0, x0, b, alpha, dt, nsub, fluid):
    w = tuple(w0)
    x = tuple(x0)
    dt_sub = dt / nsub
    first_ce = None
    for _ in range(nsub):
        midpoint = coupled._add(x, coupled._scale(0.5 * dt_sub,
                                                  coupled._velocity(w)))
        ce = coupled._ideal_ce(fluid(midpoint), b)
        if first_ce is None:
            first_ce = ce
        w, dx, _ = coupled._hc_step(w, ce, b, alpha, 1.0, dt_sub)
        x = coupled._add(x, dx)
    return w, x, first_ce


def _repeat_hc_stale(w0, x0, b, alpha, dt, nsub, ce):
    w = tuple(w0)
    x = tuple(x0)
    for _ in range(nsub):
        w, dx, _ = coupled._hc_step(w, ce, b, alpha, 1.0, dt / nsub)
        x = coupled._add(x, dx)
    return w, x


def _scalar_error(actual, expected):
    return abs(actual - expected)


def _assert_source_order(root):
    source = {}
    for name in (
            "src/particles/particles_pushers.cpp",
            "src/particles/particles_tasks.cpp",
            "src/particles/particles.cpp",
            "src/mesh/mesh_refinement.cpp",
            "src/mesh/mesh.cpp"):
        source[name] = (root / name).read_text()

    pushers = source["src/particles/particles_pushers.cpp"]
    if "if (static_cast<Real>(n) < ratio) {++n;}" not in pushers:
        raise RuntimeError("ceiling helper no longer preserves exact integer ratios")
    if "envelope.w_min - kick" not in pushers:
        raise RuntimeError("gyro schedule no longer uses the conservative momentum floor")
    if pushers.count("envelope.alpha_max, envelope.e_max") < 2:
        raise RuntimeError("relativistic electric bound lost checked multiplication")
    if "relativistic::CheckedProduct(alpha_e, abs_dt, kick)" not in pushers:
        raise RuntimeError("relativistic electric kick lost checked multiplication")
    if pushers.count("envelope.alpha_max, envelope.b_max") < 2:
        raise RuntimeError("relativistic gyro bound lost checked multiplication")
    push_tail = source["src/particles/particles_pushers.cpp"].split(
        'CheckMotionBounds("particle push");', 1)[0].rsplit(
            'MarkRelativisticTimestepBoundDirty();', 1)[1]
    if "RefreshRelativisticTimestepBound" in push_tail:
        raise RuntimeError("push recomputes relativistic dtnew before final exchange")
    push_loop = pushers.split("for (int sub=0; sub<nsub; ++sub)", 1)[1]
    if "RunRelativisticMHDIdealStep(" not in push_loop:
        raise RuntimeError(
            "relativistic coupled fields are not re-gathered in each substep")
    relativistic_schedule = source["src/particles/particles_pushers.cpp"].split(
        "int Particles::ComputeRelativisticSubcycleSteps", 1)[1].split(
            "// LogSubcycle()", 1)[0]
    if "if (nsub > subcycle_max_steps)" not in relativistic_schedule:
        raise RuntimeError("relativistic cap preflight is missing before the first drift")
    tasks = source["src/particles/particles_tasks.cpp"]
    newdt = tasks.split("TaskStatus Particles::NewTimeStep", 1)[1].split(
        "// TaskList Particles::NewGID", 1)[0]
    if "RefreshRelativisticTimestepBound();" not in newdt:
        raise RuntimeError("post-integrator relativistic refresh hook is missing")
    if newdt.index("MarkRelativisticTimestepBoundDirty();") > newdt.index(
            "RefreshRelativisticTimestepBound();"):
        raise RuntimeError("post-integrator task refreshes before marking dirty")
    remap = source["src/particles/particles.cpp"]
    if remap.index("MarkRelativisticTimestepBoundDirty();",
                   remap.index("void Particles::RemapAfterAMR")) < remap.index(
                       'CheckConsistency("AMR remap");',
                       remap.index("void Particles::RemapAfterAMR")):
        raise RuntimeError("AMR remap dirty marker occurs before completed remap")
    refinement = source["src/mesh/mesh_refinement.cpp"]
    refresh = refinement.index("RefreshRelativisticTimestepBound();")
    if refresh < refinement.index("RepairAMRFC(") or refresh < refinement.index(
            "ConToPrim("):
        raise RuntimeError(
            "AMR refresh occurs before repair and primitive reconstruction")
    mesh = source["src/mesh/mesh.cpp"]
    if mesh.index("RefreshRelativisticTimestepBound();") > mesh.index(
            "pmb_pack->ppart->dtnew"):
        raise RuntimeError(
            "mesh consumer reads relativistic dtnew before defensive refresh")
    coupled_pgen = (root / "src/pgen/unit_tests/"
                    "cr_relativistic_coupled_runtime_test.cpp").read_text().split(
                        "InitializeLiveMHDFields(pin, pmy_mesh_);", 1)[1].split(
                            "PoisonStoredParticleFields(pin, pmy_mesh_);", 1)[0]
    if "MarkRelativisticTimestepBoundDirty();" not in coupled_pgen:
        raise RuntimeError(
            "coupled initialization does not mark the timestep bound dirty")
    if "RefreshRelativisticTimestepBound();" in coupled_pgen:
        raise RuntimeError("coupled initialization refreshes before driver boundary fill")
    return 12


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--criteria", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--metrics", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    input_path = args.input.resolve()
    root = args.work_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    coupled.WITNESS_ROWS.clear()
    metrics = {}

    quiet = _flags(time__nlim=1, time__cfl_number=100.0,
                   particles__subcycle="true", particles__log_performance="true")
    zero_field = _flags(problem__live_u1=0.0, problem__live_u2=0.0,
                        problem__live_u3=0.0, problem__live_B1=0.0,
                        problem__live_B2=0.0, problem__live_B3=0.0,
                        problem__B0x=0.0, problem__B0y=0.0, problem__B0z=0.0)

    ceil_cases = (
        ("ceil_below", 0.249999, 2),
        ("ceil_exact", 0.25, 2),
        ("ceil_above", 0.250001, 3),
    )
    for name, dt, expected in ceil_cases:
        _, schedules, _ = _run(
            binary, input_path, root, name,
            quiet + zero_field + _flags(
                time__tlim=dt, particles__subcycle_max_steps=16))
        schedule = _only_schedule(schedules, name)
        metrics[f"{name}_cell_steps"] = schedule["cell"]
        if schedule["cell"] != expected:
            raise RuntimeError(f"{name}: cell schedule {schedule} != {expected}")

    electric = quiet + _flags(
        time__tlim=0.1, particles__subcycle_max_steps=32,
        particles__subcycle_gyro_fraction=0.5,
        particles__subcycle_electric_kick_max=0.02,
        problem__live_u1=0.0, problem__live_u2=-0.7, problem__live_u3=0.0,
        problem__live_B1=0.0, problem__live_B2=0.0, problem__live_B3=1.0,
        problem__B0x=0.0, problem__B0y=0.0, problem__B0z=1.0)
    electric_row, schedules, _ = _run(binary, input_path, root, "electric", electric)
    electric_schedule = _only_schedule(schedules, "electric")
    metrics["electric_cell_steps"] = electric_schedule["cell"]
    metrics["electric_meshblock_steps"] = electric_schedule["meshblock"]
    metrics["electric_gyro_steps"] = electric_schedule["gyro"]
    metrics["electric_kick_steps"] = electric_schedule["electric"]
    metrics["electric_total_steps"] = electric_schedule["steps"]
    metrics["near_rest_total_steps"] = electric_schedule["steps"]
    _, poison_schedules, _ = _run(
        binary, input_path, root, "electric_poison", electric + _flags(
            problem__poison_B1=31.0, problem__poison_B2=-37.0,
            problem__poison_B3=41.0, problem__poison_cE1=-43.0,
            problem__poison_cE2=47.0, problem__poison_cE3=-53.0))
    poison_schedule = _only_schedule(poison_schedules, "electric_poison")
    metrics["electric_poison_schedule_difference"] = sum(
        abs(electric_schedule[key] - poison_schedule[key])
        for key in ("steps", "cell", "meshblock", "gyro", "electric"))

    strict_row, _, _ = _run(
        binary, input_path, root, "electric_strict",
        electric + _flags(particles__subcycle_electric_kick_max=0.01))
    reference_row, _, _ = _run(
        binary, input_path, root, "electric_reference",
        electric + _flags(particles__subcycle_max_steps=256,
                          particles__subcycle_electric_kick_max=0.001))
    metrics["electric_strict_reference_momentum_error"] = coupled._norm(
        coupled._sub(coupled._w(strict_row), coupled._w(reference_row)))

    momentum_variant = _particle_input_variant(
        input_path, root, "deferred_momentum",
        "subcycle_momentum_fraction       = 0.1")
    kinetic_variant = _particle_input_variant(
        input_path, root, "deferred_kinetic",
        "subcycle_kinetic_energy_fraction = 0.1")
    _expected_failure(
        binary, momentum_variant, root, "failure_deferred_momentum", [],
        "separately reviewed nonsingular near-rest scale")
    _expected_failure(
        binary, kinetic_variant, root, "failure_deferred_kinetic", [],
        "separately reviewed nonsingular near-rest scale")
    metrics["deferred_fraction_failure_probe_count"] = 2

    gyro = quiet + _flags(
        time__tlim=0.2, particles__subcycle_max_steps=64,
        particles__subcycle_gyro_fraction=0.01,
        particles__subcycle_electric_kick_max=1.0,
        problem__live_u1=0.0, problem__live_u2=0.0, problem__live_u3=0.0,
        problem__live_B1=0.0, problem__live_B2=0.0, problem__live_B3=1.0,
        problem__B0x=0.0, problem__B0y=0.0, problem__B0z=1.0)
    _, schedules, _ = _run(binary, input_path, root, "gyro_low", gyro)
    gyro_low = _only_schedule(schedules, "gyro_low")
    for key in ("cell", "meshblock", "gyro", "electric", "steps"):
        metrics[f"gyro_low_{'total' if key == 'steps' else key}_steps"] = gyro_low[key]
    _, schedules, _ = _run(
        binary, input_path, root, "gyro_high",
        gyro + _flags(problem__w0x=math.sqrt(9999.0)))
    gyro_high = _only_schedule(schedules, "gyro_high")
    for key in ("cell", "meshblock", "gyro", "electric", "steps"):
        metrics[f"gyro_high_{'total' if key == 'steps' else key}_steps"] = gyro_high[key]

    _, schedules, _ = _run(
        binary, input_path, root, "gamma_floor",
        quiet + _flags(
            time__tlim=1.0, particles__subcycle_max_steps=128,
            particles__subcycle_gyro_fraction=0.1,
            particles__subcycle_electric_kick_max=10.0,
            problem__w0x=10.0, problem__live_u1=0.0, problem__live_u2=-0.9,
            problem__live_u3=0.0, problem__live_B1=0.0, problem__live_B2=0.0,
            problem__live_B3=10.0, problem__B0x=0.0, problem__B0y=0.0,
            problem__B0z=10.0))
    metrics["gamma_floor_gyro_steps"] = _only_schedule(
        schedules, "gamma_floor")["gyro"]

    _, schedules, _ = _run(
        binary, input_path, root, "cell_bound",
        quiet + zero_field + _flags(
            time__tlim=0.3, particles__subcycle_max_steps=16,
            particles__subcycle_cell_fraction=0.25))
    metrics["cell_bound_cell_steps"] = _only_schedule(schedules, "cell_bound")["cell"]

    block_w = 100.0
    block_velocity = block_w / math.sqrt(1.0 + block_w**2)
    block_row, schedules, _ = _run(
        binary, input_path, root, "meshblock_bound",
        quiet + zero_field + _flags(
            time__tlim=0.3, particles__subcycle_max_steps=16,
            particles__subcycle_meshblock_fraction=0.02,
            problem__particle_x1=0.9, problem__w0x=block_w))
    meshblock_schedule = _only_schedule(schedules, "meshblock_bound")
    metrics["meshblock_bound_cell_steps"] = meshblock_schedule["cell"]
    metrics["meshblock_bound_meshblock_steps"] = meshblock_schedule["meshblock"]
    expected_x = ((0.9 + 0.3*block_velocity + 1.0) % 2.0) - 1.0
    metrics["meshblock_bound_wrapped_position_error"] = _scalar_error(
        block_row["x"], expected_x)
    metrics["meshblock_bound_displacement_error"] = _scalar_error(
        block_row["dx"], 0.3*block_velocity)

    cap_case = electric + _flags(
        particles__subcycle_max_steps=3,
        problem__runtime_particle_dtnew_override=0.1)
    _expected_failure(
        binary, input_path, root, "failure_cap", cap_case,
        "cap clipping is unsupported")
    _expected_failure(
        binary, input_path, root, "failure_cap_clipping_parser",
        electric + _flags(particles__subcycle_strict="false"),
        "cap clipping is unsupported")
    metrics["cap_failure_probe_count"] = 2

    regather_w0 = (0.8, 0.0, 0.0)
    regather_x0 = (-0.4, 0.0, 0.0)
    regather_b = (0.0, 0.0, 1.0)
    regather_row, schedules, _ = _run(
        binary, input_path, root, "midpoint_regather",
        quiet + _flags(
            time__tlim=0.4, particles__subcycle_max_steps=16,
            problem__particle_x1=regather_x0[0], problem__w0x=regather_w0[0],
            problem__live_velocity_profile="affine", problem__live_u1=0.0,
            problem__live_u2=0.0, problem__live_u3=0.0,
            problem__live_u2_x1=0.4, problem__live_B1=0.0,
            problem__live_B2=0.0, problem__live_B3=1.0,
            problem__B0x=0.0, problem__B0y=0.0, problem__B0z=1.0))
    regather_steps = _only_schedule(schedules, "midpoint_regather")["steps"]
    oracle_w, _, first_ce = _repeat_hc(
        regather_w0, regather_x0, regather_b, 1.0, regather_row["time"],
        regather_steps, lambda x: (0.0, 0.4*x[0], 0.0))
    stale_w, _ = _repeat_hc_stale(
        regather_w0, regather_x0, regather_b, 1.0, regather_row["time"],
        regather_steps, first_ce)
    metrics["midpoint_regather_reference_momentum_error"] = coupled._norm(
        coupled._sub(coupled._w(regather_row), oracle_w))
    metrics["midpoint_regather_stale_field_separation"] = coupled._norm(
        coupled._sub(oracle_w, stale_w))

    dtnew_common = _flags(
        time__cfl_number=100.0, particles__subcycle="true",
        particles__subcycle_max_steps=8,
        particles__subcycle_electric_kick_max=0.01,
        problem__live_u1=0.0, problem__live_u2=-0.4, problem__live_u3=0.0,
        problem__live_B1=0.0, problem__live_B2=0.0, problem__live_B3=1.0,
        problem__B0x=0.0, problem__B0y=0.0, problem__B0z=1.0)
    initial_row, _, _ = _run(
        binary, input_path, root, "dtnew_initial",
        dtnew_common + _flags(time__nlim=0, time__tlim=0.5))
    metrics["initial_outer_bound_error"] = _scalar_error(
        initial_row["particle_dtnew"], 0.2)
    evolved_row, _, cycles = _run(
        binary, input_path, root, "dtnew_evolved",
        dtnew_common + _flags(
            time__nlim=2, time__tlim=0.5, particles__log_performance="true",
            problem__user_srcs="true",
            problem__relativistic_coupled_manufactured_test="true",
            problem__manufactured_acceleration=-1.0))
    if len(cycles) < 2:
        raise RuntimeError(f"dtnew_evolved: missing cycle diagnostics: {cycles}")
    metrics["next_step_consumption_error"] = _scalar_error(
        evolved_row["last_push_dt"], 0.01/0.6*8.0)
    final_e = 0.4 + evolved_row["time"]
    metrics["evolved_outer_bound_error"] = _scalar_error(
        evolved_row["particle_dtnew"], 0.01/final_e*8.0)

    _expected_failure(
        binary, input_path, root, "failure_overflow",
        quiet + _flags(
            time__tlim=0.1, particles__subcycle_max_steps=8,
            particles__alpha_s=1e308, problem__live_u2=-0.5,
            problem__live_B3=1e10, problem__B0z=1e10,
            problem__runtime_particle_dtnew_override=0.1),
        "relativistic electric-kick subcycle request overflowed")
    metrics["overflow_failure_probe_count"] = 1
    metrics["lifecycle_source_order_assertion_count"] = _assert_source_order(
        Path(__file__).resolve().parents[2])

    coupled._validate_criteria(args.criteria, metrics)
    if args.metrics:
        args.metrics.parent.mkdir(parents=True, exist_ok=True)
        args.metrics.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    print("CR relativistic Phase-5 subcycling and timestep checks passed")


if __name__ == "__main__":
    main()
