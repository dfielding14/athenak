#!/usr/bin/env python3
"""Build and run deterministic lagrangian_mc particle-tracking tests."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

from oracle import (
    MeshSpec,
    analytic_state,
    initial_cell_center_particles,
    outer_x1_particles,
    predicted_x1_mc_push,
    uniform_state,
)
from read_trk import TrackRecord, read_trk


ROOT = Path(__file__).resolve().parents[2]
TEST_DIR = Path(__file__).resolve().parent


class TestFailure(AssertionError):
    pass


def assert_close(
    actual: float, expected: float, label: str, tol: float = 2.0e-5,
    rel: float = 1.0e-6,
) -> None:
    allowed = max(tol, rel * max(1.0, abs(expected)))
    if abs(actual - expected) > allowed:
        raise TestFailure(f"{label}: got {actual:.9g}, expected {expected:.9g}")


def assert_int(actual: float, expected: int, label: str) -> None:
    got = int(round(actual))
    if got != expected:
        raise TestFailure(f"{label}: got {got}, expected {expected}")


def cmake_build(build_dir: Path, jobs: int) -> Path:
    if not (build_dir / "CMakeCache.txt").exists():
        subprocess.run(
            [
                "cmake",
                "-S",
                str(ROOT),
                "-B",
                str(build_dir),
                "-D",
                "PROBLEM=particle_tracking_test",
                "-D",
                "CMAKE_BUILD_TYPE=Release",
            ],
            check=True,
        )
    subprocess.run(
        ["cmake", "--build", str(build_dir), "-j", str(jobs)],
        check=True,
    )
    exe = build_dir / "src" / "athena"
    if not exe.exists():
        exe = build_dir / "athena"
    if not exe.exists():
        raise FileNotFoundError(f"could not find built athena executable under {build_dir}")
    return exe


def run_athena(exe: Path, input_file: Path, run_dir: Path) -> None:
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    log = run_dir / "athena.log"
    cmd = [str(exe), "-i", str(input_file), "-d", str(run_dir)]
    with log.open("w") as stream:
        proc = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise TestFailure(f"Athena failed for {input_file.name}; see {log}")


def load_records(run_dir: Path, basename: str) -> list[TrackRecord]:
    trk = run_dir / "trk" / f"{basename}.trk"
    if not trk.exists():
        raise TestFailure(f"missing track file {trk}")
    records = read_trk(trk)
    if not records:
        raise TestFailure(f"no records parsed from {trk}")
    expected = ("tag", "x", "y", "z", "vx", "vy", "vz", "rho", "press",
                "temp", "eint", "scalar0", "gid", "level", "active",
                "edot_cool", "edot_heat", "edot_net", "dedt_rad_mass",
                "dTdt_rad", "tcool", "entropy", "ln_entropy", "divv",
                "dTdt_ad", "T_mix_scalar", "T_minus_T_mix_scalar",
                "T_label_mix", "T_minus_T_label_mix", "gradT_mag",
                "grad_scalar_mag", "strain_mag")
    if records[-1].variables != expected:
        raise TestFailure(
            f"unexpected variables in {trk}: {' '.join(records[-1].variables)}"
        )
    return records


def validate_static(records: list[TrackRecord]) -> None:
    mesh = MeshSpec(nx1=4, nx2=4, nx3=4)
    record = records[0]
    expected_positions = initial_cell_center_particles(mesh)
    if record.ntrack != len(expected_positions):
        raise TestFailure(f"static ntrack={record.ntrack}, expected {len(expected_positions)}")
    for tag, (x, y, z) in enumerate(expected_positions):
        row = record.rows[tag]
        assert_int(row["active"], 1, f"static tag {tag} active")
        assert_int(row["gid"], 0, f"static tag {tag} gid")
        assert_int(row["level"], 0, f"static tag {tag} level")
        assert_close(row["x"], x, f"static tag {tag} x")
        assert_close(row["y"], y, f"static tag {tag} y")
        assert_close(row["z"], z, f"static tag {tag} z")
        assert_close(row["vx"], 0.0, f"static tag {tag} vx")
        assert_close(row["vy"], 0.0, f"static tag {tag} vy")
        assert_close(row["vz"], 0.0, f"static tag {tag} vz")
        state = analytic_state(mesh, x, y, z)
        for name, expected in state.items():
            assert_close(row[name], expected, f"static tag {tag} {name}", tol=7.0e-5)


def validate_mc_push(records: list[TrackRecord]) -> None:
    mesh = MeshSpec(nx1=4, nx2=4, nx3=1)
    final = records[-1]
    initial_positions = initial_cell_center_particles(mesh)
    state = uniform_state()
    for tag, (x, y, z) in enumerate(initial_positions):
        row = final.rows[tag]
        assert_int(row["active"], 1, f"mc_push tag {tag} active")
        expected_x = predicted_x1_mc_push(
            tag, x, mesh, probability=0.5, random_seed=12345, cycle_used_by_pusher=0
        )
        assert_close(row["x"], expected_x, f"mc_push tag {tag} x")
        assert_close(row["y"], y, f"mc_push tag {tag} y")
        assert_close(row["z"], z, f"mc_push tag {tag} z")
        assert_close(row["vx"], 1.0, f"mc_push tag {tag} vx")
        for name, expected in state.items():
            assert_close(row[name], expected, f"mc_push tag {tag} {name}")


def validate_exit_x1(records: list[TrackRecord]) -> None:
    initial = records[0]
    final = records[-1]
    if initial.ntrack != 4:
        raise TestFailure(f"exit_x1 ntrack={initial.ntrack}, expected 4")
    for tag in range(4):
        assert_int(initial.rows[tag]["active"], 1, f"exit_x1 initial tag {tag} active")
        assert_int(final.rows[tag]["active"], 0, f"exit_x1 final tag {tag} active")


def validate_inflow_injection(records: list[TrackRecord]) -> None:
    mesh = MeshSpec(nx1=4, nx2=4, nx3=1)
    initial = records[0]
    final = records[-1]
    if any(int(round(row["active"])) for row in initial.rows):
        raise TestFailure("inflow initial record should contain no active particles")
    expected_positions = [mesh.center(0, j, 0) for j in range(4)]
    state = uniform_state()
    for tag, (x, y, z) in enumerate(expected_positions):
        row = final.rows[tag]
        assert_int(row["active"], 1, f"inflow tag {tag} active")
        assert_int(row["gid"], 0, f"inflow tag {tag} gid")
        assert_int(row["level"], 0, f"inflow tag {tag} level")
        assert_close(row["x"], x, f"inflow tag {tag} x")
        assert_close(row["y"], y, f"inflow tag {tag} y")
        assert_close(row["z"], z, f"inflow tag {tag} z")
        for name, expected in state.items():
            assert_close(row[name], expected, f"inflow tag {tag} {name}")
    for tag in range(len(expected_positions), final.ntrack):
        assert_int(final.rows[tag]["active"], 0, f"inflow unused tag {tag} active")


def validate_scheduled_injection(records: list[TrackRecord]) -> None:
    if records[-1].ntrack != 4:
        raise TestFailure(f"scheduled ntrack={records[-1].ntrack}, expected 4")

    # The slots are due at t=0.1, but the pgen suppresses top inflow until t=0.3.
    # They must remain empty while only non-selected faces inject particles.
    pre_top_records = [record for record in records if record.time < 0.3 - 1.0e-8]
    if len(pre_top_records) < 2:
        raise TestFailure("scheduled test did not produce a pre-top-inflow record")
    for record in pre_top_records:
        for slot, row in enumerate(record.rows):
            assert_int(row["active"], 0, f"scheduled pre-top slot {slot} active")
            assert_int(row["tag"], -1, f"scheduled pre-top slot {slot} tag")

    filled = next(
        (record for record in records if all(int(round(row["active"])) == 1
                                             for row in record.rows)),
        None,
    )
    if filled is None:
        raise TestFailure("scheduled slots never all became active")
    if filled.time < 0.3 - 1.0e-8:
        raise TestFailure(f"scheduled slots filled too early at t={filled.time:g}")

    state = uniform_state()
    for slot, row in enumerate(filled.rows):
        assert_int(row["active"], 1, f"scheduled filled slot {slot} active")
        if int(round(row["tag"])) < 0:
            raise TestFailure(f"scheduled filled slot {slot} has no tag")
        assert_int(row["gid"], 0, f"scheduled filled slot {slot} gid")
        assert_int(row["level"], 0, f"scheduled filled slot {slot} level")
        assert_close(row["z"], 0.875, f"scheduled filled slot {slot} z")
        assert_close(row["vz"], -1.0, f"scheduled filled slot {slot} vz")
        for name, expected in state.items():
            assert_close(row[name], expected, f"scheduled filled slot {slot} {name}")


def validate_amr_refine_derefine(records: list[TrackRecord]) -> None:
    mesh = MeshSpec(nx1=4, nx2=4, nx3=1)
    initial = records[0]
    final = records[-1]
    expected_positions = initial_cell_center_particles(mesh)
    state = uniform_state()
    if final.ntrack != len(expected_positions):
        raise TestFailure(f"amr ntrack={final.ntrack}, expected {len(expected_positions)}")
    for tag, (x, y, z) in enumerate(expected_positions):
        assert_int(initial.rows[tag]["active"], 1, f"amr initial tag {tag} active")
        row = final.rows[tag]
        assert_int(row["active"], 1, f"amr final tag {tag} active")
        assert_int(row["level"], 0, f"amr final tag {tag} level")
        assert_close(row["x"], x, f"amr final tag {tag} x")
        assert_close(row["y"], y, f"amr final tag {tag} y")
        assert_close(row["z"], z, f"amr final tag {tag} z")
        for name, expected in state.items():
            assert_close(row[name], expected, f"amr final tag {tag} {name}")


CASES = {
    "static": ("athinput.particle_static", "particle_static", validate_static),
    "mc_push": ("athinput.particle_mc_push", "particle_mc_push", validate_mc_push),
    "exit_x1": ("athinput.particle_exit_x1", "particle_exit_x1", validate_exit_x1),
    "inflow_injection": (
        "athinput.particle_inflow_injection",
        "particle_inflow_injection",
        validate_inflow_injection,
    ),
    "scheduled_injection": (
        "athinput.particle_scheduled_injection",
        "particle_scheduled_injection",
        validate_scheduled_injection,
    ),
    "amr_refine_derefine": (
        "athinput.particle_amr_refine_derefine",
        "particle_amr_refine_derefine",
        validate_amr_refine_derefine,
    ),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--athena", type=Path, help="Existing Athena executable to use")
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=TEST_DIR / "build_particle_tracking",
        help="CMake build directory used when --athena is not provided",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=TEST_DIR / "runs_particle_tracking",
        help="Scratch directory for test runs",
    )
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument(
        "--case",
        choices=["all", *CASES.keys()],
        default="all",
        help="Run one case or the full implemented suite",
    )
    args = parser.parse_args()

    try:
        exe = args.athena if args.athena is not None else cmake_build(args.build_dir, args.jobs)
        selected = list(CASES) if args.case == "all" else [args.case]
        for case in selected:
            input_name, basename, validator = CASES[case]
            run_dir = args.work_dir / case
            run_athena(exe, TEST_DIR / input_name, run_dir)
            records = load_records(run_dir, basename)
            validator(records)
            print(f"PASS {case}")
    except Exception as exc:
        print(f"FAIL {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
