#!/usr/bin/env python3
"""Run reproducible CPU FrameTracker overhead benchmarks."""

from __future__ import annotations

import argparse
import csv
import re
import statistics
import subprocess
import tempfile
from pathlib import Path


TIME_PATTERN = re.compile(r"cpu time used\s*=\s*([0-9.eE+-]+)")
CASES = (("disabled", False), ("x1", True), ("all", True))
FIELDS = (
    "revision",
    "backend",
    "ranks",
    "axes",
    "mesh_layout",
    "block_layout",
    "cycles",
    "update_count",
    "repeat",
    "elapsed_seconds",
    "disabled_baseline_seconds",
    "incremental_seconds_per_update",
    "pass_fail_criterion",
    "result",
)


def run_checked(command: list[str], cwd: Path) -> str:
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(command)}\n"
            f"{result.stdout}\n{result.stderr}"
        )
    return result.stdout + result.stderr


def revision(tree: Path) -> str:
    return run_checked(["git", "rev-parse", "--short=8", "HEAD"], tree).strip()


def require_clean_revision(tree: Path) -> None:
    for diff_command in (["git", "diff", "--quiet", "--exit-code", "HEAD", "--"],
                         ["git", "diff", "--cached", "--quiet", "--exit-code"]):
        result = subprocess.run(diff_command, cwd=tree, check=False)
        if result.returncode != 0:
            raise RuntimeError(
                f"Tracked modifications exist in {tree}; commit them before benchmarking."
            )


def build_binary(tree: Path, backend: str, build: bool) -> Path:
    build_dir = tree / f"build_frame_tracking_benchmark_{backend}"
    if build:
        configure = ["cmake", "-S", str(tree), "-B", str(build_dir),
                     "-DCMAKE_BUILD_TYPE=Release"]
        if backend == "mpicpu":
            configure.append("-DAthena_ENABLE_MPI=ON")
        run_checked(configure, tree)
        run_checked(["cmake", "--build", str(build_dir), "-j", "4"], tree)
    binary = build_dir / "src" / "athena"
    if not binary.exists():
        raise FileNotFoundError(f"Benchmark binary not found: {binary}")
    return binary


def run_once(binary: Path, input_path: Path, backend: str, ranks: int, axes: str,
             enabled: bool, cycles: int, mesh: int, block: int, run_dir: Path) -> float:
    command = [str(binary), "-i", str(input_path), "-d", str(run_dir),
               f"time/nlim={cycles}", f"mesh/nx1={mesh}", f"mesh/nx2={mesh}",
               f"mesh/nx3={mesh}", f"meshblock/nx1={block}",
               f"meshblock/nx2={block}", f"meshblock/nx3={block}",
               f"frame_tracking/enabled={'true' if enabled else 'false'}",
               f"frame_tracking/axes={axes if enabled else 'all'}"]
    if backend == "mpicpu":
        command = ["mpiexec", "-n", str(ranks)] + command
    output = run_checked(command, run_dir)
    match = TIME_PATTERN.search(output)
    if match is None:
        raise RuntimeError(f"No CPU timing found in output for {backend}/{axes}.")
    return float(match.group(1))


def measure_revision(tree: Path, input_path: Path, backend: str, ranks: int, warmups: int,
                     repeats: int, cycles: int, mesh: int, block: int,
                     build: bool) -> tuple[str, dict[str, list[float]]]:
    rev = revision(tree)
    binary = build_binary(tree, backend, build)
    measurements: dict[str, list[float]] = {case[0]: [] for case in CASES}
    with tempfile.TemporaryDirectory(prefix=f"frame_benchmark_{rev}_{backend}_") as root:
        root_path = Path(root)
        for name, enabled in CASES:
            for repeat in range(warmups + repeats):
                run_dir = root_path / f"{name}_{repeat}"
                run_dir.mkdir()
                elapsed = run_once(binary, input_path, backend, ranks, name,
                                   enabled, cycles, mesh, block, run_dir)
                if repeat >= warmups:
                    measurements[name].append(elapsed)
    return rev, measurements


def median_incremental(values: dict[str, list[float]], axes: str, updates: int) -> float:
    return (
        statistics.median(values[axes]) - statistics.median(values["disabled"])
    ) / updates


def add_rows(rows: list[dict[str, object]], rev: str, backend: str, ranks: int,
             values: dict[str, list[float]], cycles: int, mesh: int, block: int) -> None:
    updates = cycles
    disabled_median = statistics.median(values["disabled"])
    for axes, enabled in CASES:
        for repeat, elapsed in enumerate(values[axes], start=1):
            disabled = values["disabled"][repeat - 1]
            incremental = "" if not enabled else (elapsed - disabled) / updates
            rows.append({
                "revision": rev, "backend": backend, "ranks": ranks, "axes": axes,
                "mesh_layout": f"{mesh}x{mesh}x{mesh}",
                "block_layout": f"{block}x{block}x{block}",
                "cycles": cycles, "update_count": updates if enabled else 0,
                "repeat": repeat, "elapsed_seconds": elapsed,
                "disabled_baseline_seconds": disabled,
                "incremental_seconds_per_update": incremental,
                "pass_fail_criterion": "raw measurement", "result": "measured",
            })
        median_elapsed = statistics.median(values[axes])
        rows.append({
            "revision": rev, "backend": backend, "ranks": ranks, "axes": axes,
            "mesh_layout": f"{mesh}x{mesh}x{mesh}",
            "block_layout": f"{block}x{block}x{block}",
            "cycles": cycles, "update_count": updates if enabled else 0,
            "repeat": "median", "elapsed_seconds": median_elapsed,
            "disabled_baseline_seconds": disabled_median,
            "incremental_seconds_per_update": (
                "" if not enabled else (median_elapsed - disabled_median) / updates
            ),
            "pass_fail_criterion": "summary row", "result": "measured",
        })


def apply_criteria(rows: list[dict[str, object]], baseline_rev: str,
                   candidate_rev: str,
                   results: dict[tuple[str, str], dict[str, list[float]]],
                   cycles: int) -> None:
    for backend in ("cpu",):
        baseline = results[(baseline_rev, backend)]
        candidate = results[(candidate_rev, backend)]
        baseline_x1 = median_incremental(baseline, "x1", cycles)
        baseline_all = median_incremental(baseline, "all", cycles)
        candidate_x1 = median_incremental(candidate, "x1", cycles)
        candidate_all = median_incremental(candidate, "all", cycles)
        for row in rows:
            if row["revision"] != candidate_rev or row["backend"] != backend:
                continue
            if row["repeat"] != "median":
                continue
            if row["axes"] == "x1":
                row["pass_fail_criterion"] = (
                    "incremental cost <= baseline x1 cost * 1.05"
                )
                row["result"] = (
                    "pass" if baseline_x1 > 0.0 and candidate_x1 <= 1.05 * baseline_x1
                    else "fail"
                )
            if row["axes"] == "all":
                row["pass_fail_criterion"] = (
                    "incremental cost <= baseline all cost * 0.75"
                )
                row["result"] = (
                    "pass" if baseline_all > 0.0 and candidate_all <= 0.75 * baseline_all
                    else "fail"
                )
    for row in rows:
        if (row["revision"] == candidate_rev and row["backend"] == "mpicpu" and
                row["repeat"] == "median" and row["axes"] in ("x1", "all")):
            row["pass_fail_criterion"] = (
                "recorded comparison; agreement tested separately"
            )
            row["result"] = "measured"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_plot(path: Path, baseline_rev: str, candidate_rev: str,
               results: dict[tuple[str, str], dict[str, list[float]]],
               cycles: int) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Install matplotlib to write the benchmark plot.") from exc
    labels = ("x1", "all")
    x = range(len(labels))
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5), sharey=True)
    for panel, backend in zip(axes, ("cpu", "mpicpu")):
        baseline = [
            1.0e6 * median_incremental(results[(baseline_rev, backend)], label, cycles)
            for label in labels
        ]
        candidate = [
            1.0e6 * median_incremental(results[(candidate_rev, backend)], label, cycles)
            for label in labels
        ]
        panel.bar([value - 0.18 for value in x], baseline,
                  width=0.36, label=baseline_rev)
        panel.bar([value + 0.18 for value in x], candidate,
                  width=0.36, label=candidate_rev)
        panel.set_xticks(list(x), labels)
        panel.set_title(backend)
        panel.set_xlabel("Active axes")
        panel.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("Incremental seconds/update (microseconds)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-tree", required=True, type=Path)
    parser.add_argument("--candidate-tree", required=True, type=Path)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--cycles", type=int, default=200)
    parser.add_argument("--mesh", type=int, default=64)
    parser.add_argument("--block", type=int, default=16)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    candidate_tree = args.candidate_tree.resolve()
    baseline_tree = args.baseline_tree.resolve()
    if not args.allow_dirty:
        require_clean_revision(baseline_tree)
        require_clean_revision(candidate_tree)
    input_path = (
        args.input or candidate_tree / "inputs/hydro/frame_tracking_benchmark.athinput"
    )
    output = (
        args.output or candidate_tree / "docs/source/_static/frame_tracking_benchmark.csv"
    )
    plot = (
        args.plot or candidate_tree / "docs/source/_static/frame_tracking_benchmark.png"
    )

    results: dict[tuple[str, str], dict[str, list[float]]] = {}
    rows: list[dict[str, object]] = []
    baseline_rev = revision(baseline_tree)
    candidate_rev = revision(candidate_tree)
    revisions = ((baseline_tree, baseline_rev), (candidate_tree, candidate_rev))
    for tree, rev in revisions:
        for backend, ranks in (("cpu", 1), ("mpicpu", args.ranks)):
            _, values = measure_revision(
                tree, input_path, backend, ranks, args.warmups, args.repeats,
                args.cycles, args.mesh, args.block, not args.skip_build
            )
            results[(rev, backend)] = values
            add_rows(rows, rev, backend, ranks, values,
                     args.cycles, args.mesh, args.block)

    apply_criteria(rows, baseline_rev, candidate_rev, results, args.cycles)
    write_csv(output, rows)
    write_plot(plot, baseline_rev, candidate_rev, results, args.cycles)
    print(f"Wrote {output}")
    print(f"Wrote {plot}")


if __name__ == "__main__":
    main()
