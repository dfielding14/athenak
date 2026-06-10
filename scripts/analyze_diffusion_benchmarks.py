#!/usr/bin/env python3
"""Summarize labeled AthenaK diffusion benchmarks and plot baseline versus optimized."""

import argparse
import csv
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


BEGIN_PATTERN = re.compile(
    r"DIFFUSION_OPT_BEGIN version=(?P<version>\S+) case=(?P<case>\S+)"
)
RATE_PATTERN = re.compile(
    r"cycle=(?P<cycle>\d+).*zone-cycles/s=(?P<rate>[0-9.eE+-]+)"
)


def parse_logs(paths, first_cycle, last_cycle, nodes):
    records = []
    for path in paths:
        current = None
        rates = []
        for line in path.read_text(errors="replace").splitlines():
            begin = BEGIN_PATTERN.search(line)
            if begin:
                current = begin.groupdict()
                current["source"] = str(path)
                rates = []
                continue
            if current is None:
                continue
            match = RATE_PATTERN.search(line)
            if match:
                cycle = int(match.group("cycle"))
                if first_cycle <= cycle <= last_cycle:
                    rates.append(float(match.group("rate")) / nodes)
            if line.startswith("DIFFUSION_OPT_END"):
                expected = last_cycle - first_cycle + 1
                if len(rates) != expected:
                    raise ValueError(
                        f"{path}: {current['case']} has {len(rates)} measurements; "
                        f"expected {expected}"
                    )
                harmonic = len(rates) / sum(1.0 / value for value in rates)
                record = {
                    **current,
                    "method": current["case"].split("_", maxsplit=1)[0],
                    "window_aggregate": harmonic,
                    "mean": statistics.mean(rates),
                    "median": statistics.median(rates),
                    "minimum": min(rates),
                    "maximum": max(rates),
                    "stdev": statistics.stdev(rates),
                    "cv_percent": 100.0 * statistics.stdev(rates)
                    / statistics.mean(rates),
                    "samples": len(rates),
                }
                records.append(record)
                current = None
                rates = []
    return records


def write_csv(records, path):
    fields = [
        "version",
        "case",
        "method",
        "window_aggregate",
        "mean",
        "median",
        "minimum",
        "maximum",
        "stdev",
        "cv_percent",
        "samples",
        "source",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def aggregate(records):
    grouped = defaultdict(list)
    for record in records:
        grouped[(record["method"], record["version"])].append(record)

    summary = {}
    for key, values in grouped.items():
        window = [record["window_aggregate"] for record in values]
        cycle_min = min(record["minimum"] for record in values)
        cycle_max = max(record["maximum"] for record in values)
        summary[key] = {
            "value": statistics.median(window),
            "minimum": cycle_min,
            "maximum": cycle_max,
            "runs": len(values),
        }
    return summary


def plot_summary(summary, path):
    methods = ["explicit", "sts"]
    versions = ["baseline", "optimized"]
    colors = {"baseline": "#4C78A8", "optimized": "#E45756"}
    markers = {"baseline": "o", "optimized": "s"}
    offsets = {"baseline": -0.12, "optimized": 0.12}

    fig, ax = plt.subplots(figsize=(3.25, 3.0))
    for version in versions:
        xs = []
        ys = []
        lower = []
        upper = []
        for index, method in enumerate(methods):
            key = (method, version)
            if key not in summary:
                continue
            value = summary[key]["value"]
            xs.append(index + offsets[version])
            ys.append(value / 1.0e9)
            lower.append((value - summary[key]["minimum"]) / 1.0e9)
            upper.append((summary[key]["maximum"] - value) / 1.0e9)
        if xs:
            ax.errorbar(
                xs,
                ys,
                yerr=np.array([lower, upper]),
                fmt=markers[version],
                color=colors[version],
                capsize=3,
                markersize=5,
                linewidth=1.2,
                label=version.capitalize(),
            )

    for index, method in enumerate(methods):
        baseline = summary.get((method, "baseline"))
        optimized = summary.get((method, "optimized"))
        if baseline and optimized:
            speedup = optimized["value"] / baseline["value"]
            speedup_format = ".3f" if speedup < 1.01 else ".2f"
            top = max(baseline["maximum"], optimized["maximum"]) / 1.0e9
            ax.annotate(
                f"{speedup:{speedup_format}}$\\times$",
                (index, top),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
            )

    ax.set_xticks(range(len(methods)), ["Explicit", "STS"])
    ax.set_ylabel(r"$10^9 \times$ zone cycles s$^{-1}$ node$^{-1}$")
    ax.margins(y=0.12)
    ax.grid(axis="y", color="0.88", linewidth=0.6)
    ax.legend(frameon=False, fontsize=8)
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout(pad=0.35)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("logs", nargs="+", type=Path)
    parser.add_argument("--nodes", type=int, default=2)
    parser.add_argument("--first-cycle", type=int, default=6)
    parser.add_argument("--last-cycle", type=int, default=20)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()

    records = parse_logs(
        args.logs, args.first_cycle, args.last_cycle, args.nodes
    )
    write_csv(records, args.csv)
    summary = aggregate(records)
    if args.plot:
        plot_summary(summary, args.plot)

    for key in sorted(summary):
        item = summary[key]
        print(
            f"{key[0]:8s} {key[1]:9s} "
            f"{item['value']:.6e} zcs/s/node "
            f"[{item['minimum']:.6e}, {item['maximum']:.6e}] "
            f"runs={item['runs']}"
        )


if __name__ == "__main__":
    main()
