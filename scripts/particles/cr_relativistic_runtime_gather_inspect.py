#!/usr/bin/env python3
"""Validate diagnostic-only relativistic CR runtime gather TSV files."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


FLOAT_FIELDS = (
    "time", "x", "y", "z", "u1", "u2", "u3", "b1", "b2", "b3",
    "cE1", "cE2", "cE3", "cE_dot_b",
)
INT_FIELDS = ("cycle", "rank", "gid", "species", "tag")


RuntimeRows = tuple[list[dict[str, float | int]], list[int], list[int]]


def _read_rows(path: Path) -> RuntimeRows:
    lines = path.read_text().splitlines()
    if len(lines) < 4 or lines[0] != (
            "# diagnostic-only\tpost-finalize\ttrilinear_cell_centered"
            "\tcE=-u_cross_B"):
        raise ValueError(f"{path}: missing diagnostic-only metadata header")
    if not lines[1].startswith("# rank_counts\t"):
        raise ValueError(f"{path}: missing rank-count metadata")
    rank_counts = []
    for expected_rank, token in enumerate(lines[1].split("\t")[1:]):
        raw_rank, raw_count = token.split("=", maxsplit=1)
        if int(raw_rank) != expected_rank:
            raise ValueError(f"{path}: rank-count metadata is not contiguous")
        rank_counts.append(int(raw_count))
    if not lines[2].startswith("# rank_participants\tworld_size="):
        raise ValueError(f"{path}: missing rank-participant witness")
    participant_tokens = lines[2].split("\t")
    world_size = int(participant_tokens[1].split("=", maxsplit=1)[1])
    participants = [int(token) for token in participant_tokens[2:]]
    if participants != list(range(world_size)):
        raise ValueError(f"{path}: rank-participant witness is not contiguous")
    if len(rank_counts) != world_size:
        raise ValueError(f"{path}: rank counts and participant witness disagree")
    if not lines[3].startswith("# "):
        raise ValueError(f"{path}: missing column header")
    reader = csv.DictReader([lines[3][2:]] + lines[4:], delimiter="\t")
    rows: list[dict[str, float | int]] = []
    for raw in reader:
        row: dict[str, float | int] = {}
        for field in INT_FIELDS:
            row[field] = int(raw[field])
        for field in FLOAT_FIELDS:
            row[field] = float(raw[field])
        rows.append(row)
    return rows, rank_counts, participants


def _cross_after_gather(u: tuple[float, float, float],
                        b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        -u[1]*b[2] + u[2]*b[1],
        -u[2]*b[0] + u[0]*b[2],
        -u[0]*b[1] + u[1]*b[0],
    )


def _expected_u(args: argparse.Namespace, row: dict[str, float | int]
                ) -> tuple[float, float, float]:
    x, y, z = (float(row["x"]), float(row["y"]), float(row["z"]))
    u = [args.u1, args.u2, args.u3]
    if args.profile == "affine":
        u[0] += args.u1_x1*x + args.u1_x2*y + args.u1_x3*z
        u[1] += args.u2_x1*x + args.u2_x2*y + args.u2_x3*z
        u[2] += args.u3_x1*x + args.u3_x2*y + args.u3_x3*z
    elif args.profile == "periodic":
        k = 2.0*math.pi*args.u_wave
        u[0] += args.u_amplitude*math.sin(k*x)*math.cos(k*y)
        u[1] += args.u_amplitude*math.sin(k*y)*math.cos(k*z)
        u[2] += args.u_amplitude*math.sin(k*z)*math.cos(k*x)
    return tuple(u)


def _expected_b(args: argparse.Namespace, row: dict[str, float | int]
                ) -> tuple[float, float, float]:
    x, y, z = (float(row["x"]), float(row["y"]), float(row["z"]))
    b = [args.b1, args.b2, args.b3]
    if args.b_profile == "sinusoidal_divb_free":
        k = 2.0*math.pi*args.b_wave
        b[0] += args.b_amplitude*k*math.sin(k*x)*(math.cos(k*y) - math.cos(k*z))
        b[1] += args.b_amplitude*k*math.sin(k*y)*(math.cos(k*z) - math.cos(k*x))
        b[2] += args.b_amplitude*k*math.sin(k*z)*(math.cos(k*x) - math.cos(k*y))
    return tuple(b)


def _near(actual: float, expected: float, tolerance: float) -> bool:
    return abs(actual - expected) <= tolerance*max(1.0, abs(expected))


def _require_contiguous_tags(args: argparse.Namespace,
                             rows: list[dict[str, float | int]]) -> None:
    if args.expected_contiguous_tags is None:
        return
    expected = {
        (args.expected_species, tag)
        for tag in range(args.expected_tag_start,
                         args.expected_tag_start + args.expected_contiguous_tags)
    }
    actual = {(int(row["species"]), int(row["tag"])) for row in rows}
    if actual != expected:
        raise ValueError("particle identity set is not the expected contiguous tag set")


def _require_boundary_tag_probes(args: argparse.Namespace,
                                 rows: list[dict[str, float | int]]) -> None:
    if not args.require_boundary_tag_probes:
        return
    mid = tuple(0.5*(lo + hi) for lo, hi in zip(args.mesh_min, args.mesh_max))
    epsilon = tuple(args.probe_epsilon*(hi - lo)
                    for lo, hi in zip(args.mesh_min, args.mesh_max))
    anchor = (
        mid[0] + 0.17*(args.mesh_max[0] - args.mesh_min[0]),
        mid[1] - 0.13*(args.mesh_max[1] - args.mesh_min[1]),
        mid[2] + 0.11*(args.mesh_max[2] - args.mesh_min[2]),
    )
    expected_by_pattern = (
        (args.mesh_min[0] + epsilon[0], anchor[1], anchor[2]),
        (args.mesh_max[0] - epsilon[0], anchor[1], anchor[2]),
        (mid[0] - epsilon[0], anchor[1], anchor[2]),
        (mid[0] + epsilon[0], anchor[1], anchor[2]),
        (anchor[0], args.mesh_min[1] + epsilon[1], anchor[2]),
        (anchor[0], args.mesh_max[1] - epsilon[1], anchor[2]),
        (anchor[0], mid[1] - epsilon[1], anchor[2]),
        (anchor[0], mid[1] + epsilon[1], anchor[2]),
        (anchor[0], anchor[1], args.mesh_min[2] + epsilon[2]),
        (anchor[0], anchor[1], args.mesh_max[2] - epsilon[2]),
        (anchor[0], anchor[1], mid[2] - epsilon[2]),
        (anchor[0], anchor[1], mid[2] + epsilon[2]),
        tuple(lo + eps for lo, eps in zip(args.mesh_min, epsilon)),
        tuple(hi - eps for hi, eps in zip(args.mesh_max, epsilon)),
    )
    observed_patterns = set()
    for row in rows:
        pattern = int(row["tag"]) % len(expected_by_pattern)
        observed_patterns.add(pattern)
        actual = (float(row["x"]), float(row["y"]), float(row["z"]))
        expected = expected_by_pattern[pattern]
        if not all(abs(a - e) <= args.probe_tolerance
                   for a, e in zip(actual, expected)):
            raise ValueError(
                f"tag {row['tag']}: boundary probe {actual} expected {expected}")
    if observed_patterns != set(range(len(expected_by_pattern))):
        raise ValueError("boundary probe layout does not cover every required pattern")


def _require_refined_region_tag_probes(args: argparse.Namespace,
                                       rows: list[dict[str, float | int]]) -> None:
    if not args.require_refined_region_tag_probes:
        return
    mid = tuple(0.5*(lo + hi) for lo, hi in zip(args.mesh_min, args.mesh_max))
    epsilon = tuple(args.probe_epsilon*(hi - lo)
                    for lo, hi in zip(args.mesh_min, args.mesh_max))
    anchor = (
        mid[0] + 0.17*(args.mesh_max[0] - args.mesh_min[0]),
        mid[1] - 0.13*(args.mesh_max[1] - args.mesh_min[1]),
        mid[2] + 0.11*(args.mesh_max[2] - args.mesh_min[2]),
    )
    h = args.refined_half_width
    expected_by_pattern = (
        (mid[0] - h - epsilon[0], anchor[1], anchor[2]),
        (mid[0] - h + epsilon[0], anchor[1], anchor[2]),
        (mid[0] + h - epsilon[0], anchor[1], anchor[2]),
        (mid[0] + h + epsilon[0], anchor[1], anchor[2]),
        (anchor[0], mid[1] - h - epsilon[1], anchor[2]),
        (anchor[0], mid[1] - h + epsilon[1], anchor[2]),
        (anchor[0], mid[1] + h - epsilon[1], anchor[2]),
        (anchor[0], mid[1] + h + epsilon[1], anchor[2]),
        (anchor[0], anchor[1], mid[2] - h - epsilon[2]),
        (anchor[0], anchor[1], mid[2] - h + epsilon[2]),
        (anchor[0], anchor[1], mid[2] + h - epsilon[2]),
        (anchor[0], anchor[1], mid[2] + h + epsilon[2]),
    )
    observed_patterns = set()
    for row in rows:
        pattern = int(row["tag"]) % len(expected_by_pattern)
        observed_patterns.add(pattern)
        actual = (float(row["x"]), float(row["y"]), float(row["z"]))
        expected = expected_by_pattern[pattern]
        if not all(abs(a - e) <= args.probe_tolerance
                   for a, e in zip(actual, expected)):
            raise ValueError(
                f"tag {row['tag']}: refined-region probe {actual} expected {expected}")
    if observed_patterns != set(range(len(expected_by_pattern))):
        raise ValueError("refined-region probe layout lacks a required pattern")


def _require_interface_planes(args: argparse.Namespace,
                              rows: list[dict[str, float | int]]) -> None:
    coordinate_index = {"x": 0, "y": 1, "z": 2}
    for spec in args.require_interface_plane:
        axis, raw_plane = spec.split("=", maxsplit=1)
        if axis not in coordinate_index:
            raise ValueError(f"invalid interface-plane coordinate '{axis}'")
        plane = float(raw_plane)
        values = [float(row[axis]) for row in rows]
        lower = sum(0.0 < plane - value <= args.interface_window for value in values)
        upper = sum(0.0 < value - plane <= args.interface_window for value in values)
        if lower < args.interface_min_rows or upper < args.interface_min_rows:
            raise ValueError(
                f"interface plane {spec} lacks two-sided adjacent coverage: "
                f"lower={lower} upper={upper}")


def _validate_structure(rows: list[dict[str, float | int]],
                        rank_counts: list[int], participants: list[int]) -> None:
    if participants != list(range(len(rank_counts))):
        raise ValueError("rank-participant witness does not match recorded communicator")
    keys = [(int(row["species"]), int(row["tag"])) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate (species, tag) rows")
    sort_keys = [
        (int(row["species"]), int(row["tag"]), int(row["rank"]), int(row["gid"]))
        for row in rows
    ]
    if sort_keys != sorted(sort_keys):
        raise ValueError("rows are not deterministically sorted")
    observed_counts = [0]*len(rank_counts)
    for row in rows:
        rank = int(row["rank"])
        if rank < 0 or rank >= len(observed_counts):
            raise ValueError("row contains rank outside recorded communicator")
        observed_counts[rank] += 1
    if observed_counts != rank_counts:
        raise ValueError(
            f"rank-count metadata {rank_counts} does not match rows {observed_counts}")


def _hash_particle_state(value: int) -> int:
    mask = (1 << 64) - 1
    value = (value + 0x9E3779B97F4A7C15) & mask
    value = ((value ^ (value >> 30))*0xBF58476D1CE4E5B9) & mask
    value = ((value ^ (value >> 27))*0x94D049BB133111EB) & mask
    return value ^ (value >> 31)


def _tag_random_unit(seed: int, species: int, tag: int, stream: int) -> float:
    mask = (1 << 64) - 1
    value = seed ^ (species << 32)
    value ^= tag & 0xFFFFFFFF
    value ^= ((stream + 1)*0xD2B74407B1CE6E93) & mask
    return (_hash_particle_state(value) >> 11)/9007199254740992.0


def _require_tag_random_motion(args: argparse.Namespace,
                               rows: list[dict[str, float | int]]) -> None:
    if not args.require_tag_random_motion:
        return
    moved = 0
    for row in rows:
        initial = tuple(
            lo + _tag_random_unit(args.particle_seed, int(row["species"]),
                                  int(row["tag"]), stream)*(hi - lo)
            for stream, (lo, hi) in enumerate(zip(args.mesh_min, args.mesh_max))
        )
        actual = (float(row["x"]), float(row["y"]), float(row["z"]))
        distances = []
        for position, start, lo, hi in zip(actual, initial, args.mesh_min,
                                           args.mesh_max):
            length = hi - lo
            raw = abs(position - start)
            distances.append(min(raw, length - raw))
        if math.sqrt(sum(value*value for value in distances)) >= args.min_motion:
            moved += 1
    if moved < args.min_moved_particles:
        raise ValueError(
            f"expected at least {args.min_moved_particles} moved particles, "
            f"found {moved}")


def _validate(args: argparse.Namespace,
              rows: list[dict[str, float | int]], rank_counts: list[int],
              participants: list[int]) -> None:
    _validate_structure(rows, rank_counts, participants)
    if args.expected_count is not None and len(rows) != args.expected_count:
        raise ValueError(f"expected {args.expected_count} rows, found {len(rows)}")
    if args.expected_cycle is not None:
        cycles = {int(row["cycle"]) for row in rows}
        if cycles != {args.expected_cycle}:
            raise ValueError(f"expected cycle {args.expected_cycle}, found {cycles}")
    if args.expected_time is not None:
        times = {float(row["time"]) for row in rows}
        if not times or not all(_near(time, args.expected_time,
                                      args.time_tolerance) for time in times):
            raise ValueError(f"expected time {args.expected_time}, found {times}")
    if args.nranks is not None:
        if len(rank_counts) != args.nranks:
            raise ValueError(
                f"expected {args.nranks} rank counts, found {len(rank_counts)}")
        if args.require_empty_particle_rank and not any(count == 0
                                                        for count in rank_counts):
            raise ValueError("expected at least one particle-empty MPI rank")
    _require_contiguous_tags(args, rows)
    _require_boundary_tag_probes(args, rows)
    _require_refined_region_tag_probes(args, rows)
    _require_interface_planes(args, rows)
    _require_tag_random_motion(args, rows)

    for row in rows:
        if not all(math.isfinite(float(row[field])) for field in FLOAT_FIELDS):
            raise ValueError(f"tag {row['tag']}: nonfinite diagnostic")
        expected_u = _expected_u(args, row)
        expected_b = _expected_b(args, row)
        actual_u = (float(row["u1"]), float(row["u2"]), float(row["u3"]))
        actual_b = (float(row["b1"]), float(row["b2"]), float(row["b3"]))
        expected_ce = _cross_after_gather(actual_u, actual_b)
        actual_ce = (float(row["cE1"]), float(row["cE2"]), float(row["cE3"]))
        comparisons = [("cE construction", actual_ce, expected_ce, args.tolerance)]
        if not args.skip_profile_check:
            comparisons = [
                ("u", actual_u, expected_u, args.profile_tolerance),
                ("B", actual_b, expected_b, args.profile_tolerance),
            ] + comparisons
        for label, actual, expected, tolerance in comparisons:
            if not all(_near(a, e, tolerance) for a, e in zip(actual, expected)):
                raise ValueError(
                    f"tag {row['tag']}: {label}={actual} expected {expected}")
        scale = max(1.0, math.sqrt(sum(value*value for value in actual_ce)) *
                    math.sqrt(sum(value*value for value in actual_b)))
        if abs(float(row["cE_dot_b"])) > args.tolerance*scale:
            raise ValueError(f"tag {row['tag']}: cE dot B invariant failed")


def _compare(left: RuntimeRows, right: RuntimeRows, tolerance: float) -> None:
    def keyed(rows: list[dict[str, float | int]]
              ) -> dict[tuple[int, int], dict[str, float | int]]:
        return {(int(row["species"]), int(row["tag"])): row for row in rows}

    left_rows, left_counts, left_participants = left
    right_rows, right_counts, right_participants = right
    _validate_structure(left_rows, left_counts, left_participants)
    _validate_structure(right_rows, right_counts, right_participants)
    left_by_key = keyed(left_rows)
    right_by_key = keyed(right_rows)
    if left_by_key.keys() != right_by_key.keys():
        raise ValueError("decomposition comparison changed particle identity set")
    physical_fields = ("cycle", "time", "x", "y", "z", "u1", "u2", "u3",
                       "b1", "b2", "b3", "cE1", "cE2", "cE3", "cE_dot_b")
    for key, left_row in left_by_key.items():
        right_row = right_by_key[key]
        for field in physical_fields:
            if not _near(float(left_row[field]), float(right_row[field]), tolerance):
                raise ValueError(
                    f"decomposition comparison changed {field} for {key}: "
                    f"{left_row[field]} != {right_row[field]}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("tsv", type=Path)
    parser.add_argument("--compare", type=Path)
    parser.add_argument("--profile", choices=("uniform", "affine", "periodic"),
                        default="uniform")
    parser.add_argument("--expected-count", type=int)
    parser.add_argument("--expected-cycle", type=int)
    parser.add_argument("--expected-time", type=float)
    parser.add_argument("--time-tolerance", type=float, default=2.0e-12)
    parser.add_argument("--expected-contiguous-tags", type=int)
    parser.add_argument("--expected-tag-start", type=int, default=0)
    parser.add_argument("--expected-species", type=int, default=0)
    parser.add_argument("--nranks", type=int)
    parser.add_argument("--require-empty-particle-rank", action="store_true")
    parser.add_argument("--require-boundary-tag-probes", action="store_true")
    parser.add_argument("--require-refined-region-tag-probes", action="store_true")
    parser.add_argument("--refined-half-width", type=float, default=0.25)
    parser.add_argument("--mesh-min", type=float, nargs=3, default=(-0.5, -0.5, -0.5))
    parser.add_argument("--mesh-max", type=float, nargs=3, default=(0.5, 0.5, 0.5))
    parser.add_argument("--probe-epsilon", type=float, default=1.0e-10)
    parser.add_argument("--probe-tolerance", type=float, default=2.0e-12)
    parser.add_argument("--require-interface-plane", action="append", default=[])
    parser.add_argument("--interface-window", type=float, default=0.04)
    parser.add_argument("--interface-min-rows", type=int, default=1)
    parser.add_argument("--tolerance", type=float, default=2.0e-12)
    parser.add_argument("--profile-tolerance", type=float, default=2.0e-12)
    parser.add_argument("--skip-profile-check", action="store_true")
    parser.add_argument("--require-tag-random-motion", action="store_true")
    parser.add_argument("--particle-seed", type=int, default=0)
    parser.add_argument("--min-motion", type=float, default=1.0e-6)
    parser.add_argument("--min-moved-particles", type=int, default=1)
    parser.add_argument("--u1", type=float, default=0.0)
    parser.add_argument("--u2", type=float, default=0.0)
    parser.add_argument("--u3", type=float, default=0.0)
    parser.add_argument("--u-amplitude", type=float, default=0.0)
    parser.add_argument("--u-wave", type=float, default=1.0)
    parser.add_argument("--b1", type=float, default=0.0)
    parser.add_argument("--b2", type=float, default=0.0)
    parser.add_argument("--b3", type=float, default=1.0)
    parser.add_argument("--b-profile", choices=("uniform", "sinusoidal_divb_free"),
                        default="uniform")
    parser.add_argument("--b-amplitude", type=float, default=0.0)
    parser.add_argument("--b-wave", type=float, default=1.0)
    for component in ("u1", "u2", "u3"):
        for coordinate in ("x1", "x2", "x3"):
            parser.add_argument(f"--{component}-{coordinate}", type=float, default=0.0)
    args = parser.parse_args()
    rows, rank_counts, participants = _read_rows(args.tsv)
    _validate(args, rows, rank_counts, participants)
    if args.compare is not None:
        _compare((rows, rank_counts, participants), _read_rows(args.compare),
                 args.tolerance)
    print(f"validated {len(rows)} CR runtime gather rows from {args.tsv}")


if __name__ == "__main__":
    main()
