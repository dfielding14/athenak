#!/usr/bin/env python3
"""Extract diffusion operator and rocprof timing tables."""

import argparse
import csv
import re
import statistics
from pathlib import Path


MATRIX_BEGIN = re.compile(
    r"DIFFUSION_MATRIX_BEGIN version=(?P<version>\S+) case=(?P<case>\S+)"
)
RATE_PATTERN = re.compile(
    r"cycle=(?P<cycle>\d+).*zone-cycles/s=(?P<rate>[0-9.eE+-]+)"
)


def parse_matrix(log_path, first_cycle, last_cycle, nodes):
    rows = []
    current = None
    rates = []
    for line in log_path.read_text(errors="replace").splitlines():
        begin = MATRIX_BEGIN.search(line)
        if begin:
            current = begin.groupdict()
            rates = []
            continue
        if current is None:
            continue
        match = RATE_PATTERN.search(line)
        if match and first_cycle <= int(match.group("cycle")) <= last_cycle:
            rates.append(float(match.group("rate")) / nodes)
        if line.startswith("DIFFUSION_MATRIX_END"):
            expected = last_cycle - first_cycle + 1
            if len(rates) != expected:
                raise ValueError(
                    f"{current['case']} has {len(rates)} measurements; "
                    f"expected {expected}"
                )
            rows.append(
                {
                    **current,
                    "window_aggregate": len(rates)
                    / sum(1.0 / value for value in rates),
                    "mean": statistics.mean(rates),
                    "minimum": min(rates),
                    "maximum": max(rates),
                }
            )
            current = None
            rates = []
    return rows


def parse_profiles(profile_root):
    rows = []
    for directory in sorted(path for path in profile_root.iterdir() if path.is_dir()):
        version, method = directory.name.split("_", maxsplit=1)
        stats_path = directory / f"{directory.name}_kernel_stats.csv"
        with stats_path.open(newline="") as stream:
            for record in csv.DictReader(stream):
                rows.append(
                    {
                        "version": version,
                        "method": method,
                        "name": record["Name"],
                        "calls": int(record["Calls"]),
                        "total_ms": int(record["TotalDurationNs"]) / 1.0e6,
                        "average_ms": float(record["AverageNs"]) / 1.0e6,
                        "percentage": float(record["Percentage"]),
                    }
                )
    return rows


def write_csv(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def print_matrix(rows):
    rates = {
        (row["version"], row["case"]): row["window_aggregate"] for row in rows
    }
    for version in ("baseline", "optimized"):
        hydro = rates[(version, "hydro")]
        for case in sorted(case for item_version, case in rates if item_version == version):
            rate = rates[(version, case)]
            relative_cost = hydro / rate
            print(
                f"{version:9s} {case:20s} {rate:.6e} zcs/s/node "
                f"cost={relative_cost:.4f}"
            )


def print_profile_changes(rows):
    selected = {
        "hydro_sts_update",
        "visc1",
        "visc2",
        "visc3",
        "conduct1",
        "conduct2",
        "conduct3",
        "hyd_c2p",
        "SendBuff",
        "RecvBuff",
        "Kokkos::Initialization Complete",
    }
    timings = {
        (row["version"], row["method"], row["name"]): row["total_ms"]
        for row in rows
    }
    for method in ("explicit", "sts"):
        print(f"\n{method}:")
        for name in sorted(selected):
            baseline = timings.get(("baseline", method, name), 0.0)
            optimized = timings.get(("optimized", method, name), 0.0)
            if baseline == 0.0 and optimized == 0.0:
                continue
            change = 100.0 * (optimized / baseline - 1.0) if baseline else 0.0
            print(
                f"{name:32s} {baseline:9.3f} -> {optimized:9.3f} ms "
                f"({change:+7.2f}%)"
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("log", type=Path)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--nodes", type=int, default=2)
    parser.add_argument("--first-cycle", type=int, default=6)
    parser.add_argument("--last-cycle", type=int, default=10)
    parser.add_argument("--matrix-csv", type=Path)
    parser.add_argument("--profile-csv", type=Path, required=True)
    args = parser.parse_args()

    matrix = parse_matrix(
        args.log, args.first_cycle, args.last_cycle, args.nodes
    )
    profiles = parse_profiles(args.profiles)
    if args.matrix_csv:
        if not matrix:
            raise ValueError("matrix CSV requested, but the log has no matrix cases")
        write_csv(matrix, args.matrix_csv)
    write_csv(profiles, args.profile_csv)
    if matrix:
        print_matrix(matrix)
    print_profile_changes(profiles)


if __name__ == "__main__":
    main()
