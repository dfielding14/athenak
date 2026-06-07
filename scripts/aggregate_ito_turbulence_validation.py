#!/usr/bin/env python3
"""Aggregate paired Ito turbulence summaries across particle seeds."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


def aggregate_summaries(paths: list[Path]) -> dict[str, Any]:
    if not paths:
        raise ValueError("at least one summary is required")

    summaries = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    comparison_keys = sorted(
        key
        for key, value in summaries[0]["comparison"].items()
        if isinstance(value, (int, float))
    )
    for summary in summaries[1:]:
        keys = {
            key
            for key, value in summary["comparison"].items()
            if isinstance(value, (int, float))
        }
        if keys != set(comparison_keys):
            raise ValueError("comparison metrics differ between summaries")

    metrics = {}
    for key in comparison_keys:
        values = np.asarray(
            [summary["comparison"][key] for summary in summaries], dtype=np.float64
        )
        metrics[key] = {
            "values": values.tolist(),
            "mean": float(np.mean(values)),
            "standard_deviation": (
                float(np.std(values, ddof=1)) if values.size > 1 else math.nan
            ),
            "standard_error": (
                float(np.std(values, ddof=1) / math.sqrt(values.size))
                if values.size > 1
                else math.nan
            ),
        }

    runtime_ratios = np.asarray(
        [
            summary["corrected"]["runtime"]["elapsed_seconds"]
            / summary["old"]["runtime"]["elapsed_seconds"]
            for summary in summaries
        ],
        dtype=np.float64,
    )
    return {
        "summary_files": [str(path.resolve()) for path in paths],
        "seed_count": len(paths),
        "metrics": metrics,
        "runtime_wall_ratio": {
            "values": runtime_ratios.tolist(),
            "mean": float(np.mean(runtime_ratios)),
            "standard_error": (
                float(np.std(runtime_ratios, ddof=1) / math.sqrt(runtime_ratios.size))
                if runtime_ratios.size > 1
                else math.nan
            ),
        },
    }


def write_csv(path: Path, aggregate: dict[str, Any]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["metric", "mean", "standard_deviation", "standard_error"])
        for key, values in aggregate["metrics"].items():
            writer.writerow(
                [
                    key,
                    values["mean"],
                    values["standard_deviation"],
                    values["standard_error"],
                ]
            )
        runtime = aggregate["runtime_wall_ratio"]
        writer.writerow(
            ["runtime_wall_ratio", runtime["mean"], "", runtime["standard_error"]]
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summaries", type=Path, nargs="+")
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    args = parser.parse_args()

    aggregate = aggregate_summaries(args.summaries)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(aggregate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_csv(args.output_csv, aggregate)
    print(json.dumps(aggregate, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
