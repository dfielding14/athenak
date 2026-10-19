#!/usr/bin/env python3
"""Audit rebuilt-PDF bin ranges, resolution, occupancy, and refinement costs."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


DEFAULT_PRODUCTION_ROOT = Path(
    "/lustre/orion/ast207/proj-shared/gotham/analysis/products/pdf_rebuild/production"
)
DEFAULT_ANALYSIS_ROOT = Path("/lustre/orion/ast207/proj-shared/gotham/analysis")


def _json_ready(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, Mapping):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_ready(item) for item in value]
    return value


def _write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(_json_ready(value), stream, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(path)


def _percentiles(values: Iterable[float]) -> dict[str, float]:
    array = np.asarray(list(values), dtype=np.float64)
    if array.size == 0:
        return {}
    return {
        "min": float(np.min(array)),
        "p05": float(np.percentile(array, 5)),
        "median": float(np.median(array)),
        "p95": float(np.percentile(array, 95)),
        "max": float(np.max(array)),
    }


def _weighted_quantile_indices(values: np.ndarray, lower: float, upper: float) -> list[int]:
    total = float(np.sum(values, dtype=np.float64))
    if total <= 0.0:
        return [-1, -1]
    cumulative = np.cumsum(values, dtype=np.float64) / total
    first = int(np.searchsorted(cumulative, lower, side="left"))
    last = int(np.searchsorted(cumulative, upper, side="left"))
    return [min(first, values.size - 1), min(last, values.size - 1)]


def _axis_resolution(dimension: Mapping[str, Any]) -> dict[str, float]:
    edges = np.asarray(dimension["bin_edges"], dtype=np.float64)
    scale = str(dimension["scale"])
    if scale == "log":
        dex = float(np.log10(edges[-1] / edges[0]) / (edges.size - 1))
        return {
            "dex_per_bin": dex,
            "multiplicative_step": float(10.0**dex),
        }
    if scale == "symlog":
        return {
            "linear_core_width": float(2.0 * dimension["linthresh"]),
            "transformed_bins": float(edges.size - 1),
        }
    return {"linear_width": float((edges[-1] - edges[0]) / (edges.size - 1))}


def _axis_signature(dimension: Mapping[str, Any]) -> str:
    return "|".join(
        (
            str(dimension["variable"]),
            str(dimension["nbin"]),
            f"{float(dimension['bin_edges'][0]):.17g}",
            f"{float(dimension['bin_edges'][-1]):.17g}",
            str(dimension["scale"]),
            f"{float(dimension.get('linthresh', 1.0)):.17g}",
        )
    )


def _axis_record(
    *,
    product: str,
    weight: str,
    snapshot: str,
    axis_index: int,
    dimension: Mapping[str, Any],
    marginal: np.ndarray,
) -> dict[str, Any]:
    total = float(np.sum(marginal, dtype=np.float64))
    interior = np.asarray(marginal[1:-1], dtype=np.float64)
    interior_total = float(np.sum(interior, dtype=np.float64))
    positive = interior[interior > 0.0]
    if positive.size and interior_total > 0.0:
        probability = positive / interior_total
        effective_bins = float(np.exp(-np.sum(probability * np.log(probability))))
    else:
        effective_bins = 0.0
    significant_floor = interior_total * 1.0e-8
    nbin = int(dimension["nbin"])
    return {
        "snapshot": snapshot,
        "product": product,
        "weight": weight,
        "axis_index": axis_index,
        "axis_signature": _axis_signature(dimension),
        "variable": str(dimension["variable"]),
        "nbin": nbin,
        "minimum": float(dimension["bin_edges"][0]),
        "maximum": float(dimension["bin_edges"][-1]),
        "scale": str(dimension["scale"]),
        "linthresh": float(dimension.get("linthresh", 1.0)),
        "resolution": _axis_resolution(dimension),
        "absolute_total": total,
        "absolute_interior": interior_total,
        "underflow_absolute_fraction": float(marginal[0] / total) if total else 0.0,
        "overflow_absolute_fraction": float(marginal[-1] / total) if total else 0.0,
        "excluded_absolute_fraction": (
            float((marginal[0] + marginal[-1]) / total) if total else 0.0
        ),
        "first_interior_absolute_fraction": (
            float(interior[0] / interior_total) if interior_total else 0.0
        ),
        "last_interior_absolute_fraction": (
            float(interior[-1] / interior_total) if interior_total else 0.0
        ),
        "populated_bins": int(np.count_nonzero(interior)),
        "significant_bins_1e-8": int(np.count_nonzero(interior > significant_floor)),
        "effective_bins": effective_bins,
        "effective_bin_fraction": effective_bins / nbin,
        "central_99p8_bin_indices": _weighted_quantile_indices(interior, 0.001, 0.999),
    }


def _audit_snapshot(arguments: tuple[str, str]) -> dict[str, Any]:
    snapshot_path_text, analysis_root_text = arguments
    snapshot_path = Path(snapshot_path_text)
    analysis_root = Path(analysis_root_text)
    if str(analysis_root) not in sys.path:
        sys.path.insert(0, str(analysis_root))
    from gotham_analysis.readers.pdf import read_pdf_header

    manifest = json.loads((snapshot_path / "rebuild_manifest.json").read_text())
    output_number = str(manifest["output_number"])
    snapshot = str(snapshot_path)
    records = []
    definitions = {}
    for entry in manifest["products"]:
        product = str(entry["id"])
        product_dir = snapshot_path / product
        header = read_pdf_header(product_dir / "gotham.header.pdf")
        weight = str(
            header.get("weight_variable")
            if header.get("weight") == "variable"
            else header.get("weight")
        )
        payload_path = product_dir / f"gotham.{output_number}.pdf"
        raw = np.memmap(payload_path, mode="r", dtype=np.float64)
        data = raw[1:].reshape(header["shape"], order="C")
        absolute = np.abs(data)
        for axis_index, dimension in enumerate(header["dimensions"]):
            sum_axes = tuple(index for index in range(data.ndim) if index != axis_index)
            marginal = (
                np.sum(absolute, axis=sum_axes, dtype=np.float64)
                if sum_axes
                else np.asarray(absolute)
            )
            records.append(
                _axis_record(
                    product=product,
                    weight=weight,
                    snapshot=snapshot,
                    axis_index=axis_index,
                    dimension=dimension,
                    marginal=marginal,
                )
            )
            definitions[_axis_signature(dimension)] = {
                "variable": str(dimension["variable"]),
                "nbin": int(dimension["nbin"]),
                "minimum": float(dimension["bin_edges"][0]),
                "maximum": float(dimension["bin_edges"][-1]),
                "scale": str(dimension["scale"]),
                "linthresh": float(dimension.get("linthresh", 1.0)),
                "resolution": _axis_resolution(dimension),
            }
        del data
        del raw
    return {
        "snapshot": snapshot,
        "manifest": {
            "simulation": snapshot_path.parts[-3],
            "phase": snapshot_path.parts[-2],
            "sequence": snapshot_path.parts[-1],
            "cells_processed": int(manifest["cells_processed"]),
            "elapsed_seconds": float(manifest["elapsed_seconds"]),
            "histogram_bytes_per_rank": int(manifest["histogram_bytes_per_rank"]),
        },
        "records": records,
        "axis_definitions": definitions,
    }


def _aggregate(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, int, str], list[Mapping[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[
            (str(record["product"]), int(record["axis_index"]), str(record["axis_signature"]))
        ].append(record)
    metrics = (
        "underflow_absolute_fraction",
        "overflow_absolute_fraction",
        "excluded_absolute_fraction",
        "first_interior_absolute_fraction",
        "last_interior_absolute_fraction",
        "populated_bins",
        "significant_bins_1e-8",
        "effective_bins",
        "effective_bin_fraction",
    )
    aggregates = []
    for (product, axis_index, signature), members in grouped.items():
        representative = members[0]
        aggregates.append(
            {
                "product": product,
                "weight": representative["weight"],
                "axis_index": axis_index,
                "axis_signature": signature,
                "variable": representative["variable"],
                "nbin": representative["nbin"],
                "minimum": representative["minimum"],
                "maximum": representative["maximum"],
                "scale": representative["scale"],
                "resolution": representative["resolution"],
                "snapshot_count": len(members),
                "metrics": {
                    metric: _percentiles(float(member[metric]) for member in members)
                    for metric in metrics
                },
            }
        )
    return sorted(aggregates, key=lambda item: (item["product"], item["axis_index"]))


def _top_rows(
    aggregates: Sequence[Mapping[str, Any]], metric: str, count: int = 20
) -> list[Mapping[str, Any]]:
    return sorted(
        aggregates,
        key=lambda item: float(item["metrics"][metric]["median"]),
        reverse=True,
    )[:count]


def _markdown(report: Mapping[str, Any]) -> str:
    aggregates = report["aggregates"]
    lines = [
        "# Rebuilt PDF binning audit",
        "",
        f"Scanned {report['snapshot_count']} snapshots and "
        f"{report['axis_record_count']} product-axis instances.",
        "",
        "## Unique axis definitions",
        "",
        "| variable | bins | range | scale | resolution |",
        "|---|---:|---:|---|---|",
    ]
    for item in sorted(
        report["axis_definitions"].values(),
        key=lambda value: (value["variable"], value["nbin"], value["minimum"]),
    ):
        resolution = ", ".join(
            f"{key}={value:.5g}" for key, value in item["resolution"].items()
        )
        lines.append(
            f"| `{item['variable']}` | {item['nbin']} | "
            f"{item['minimum']:.5g} to {item['maximum']:.5g} | "
            f"{item['scale']} | {resolution} |"
        )

    for heading, metric in (
        ("Largest median underflow", "underflow_absolute_fraction"),
        ("Largest median overflow", "overflow_absolute_fraction"),
        ("Largest median total excluded fraction", "excluded_absolute_fraction"),
        ("Lowest median effective-bin fraction", "effective_bin_fraction"),
    ):
        rows = (
            sorted(
                aggregates,
                key=lambda item: float(item["metrics"][metric]["median"]),
            )[:20]
            if metric == "effective_bin_fraction"
            else _top_rows(aggregates, metric)
        )
        lines.extend(
            [
                "",
                f"## {heading}",
                "",
                "| product | axis | bins | median | p95 |",
                "|---|---|---:|---:|---:|",
            ]
        )
        for row in rows:
            values = row["metrics"][metric]
            lines.append(
                f"| `{row['product']}` | `{row['variable']}` | {row['nbin']} | "
                f"{values['median']:.5g} | {values['p95']:.5g} |"
            )
    lines.append("")
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production-root", type=Path, default=DEFAULT_PRODUCTION_ROOT)
    parser.add_argument("--analysis-root", type=Path, default=DEFAULT_ANALYSIS_ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--snapshot-limit", type=int)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    production_root = args.production_root.resolve()
    analysis_root = args.analysis_root.resolve()
    output_dir = args.output_dir.resolve()
    snapshots = sorted(path.parent for path in production_root.glob("*/*/*/rebuild_manifest.json"))
    if args.snapshot_limit is not None:
        snapshots = snapshots[: args.snapshot_limit]
    if not snapshots:
        raise SystemExit(f"No rebuilt snapshots found beneath {production_root}")
    if args.workers <= 0:
        raise SystemExit("--workers must be positive")

    results = []
    tasks = [(str(snapshot), str(analysis_root)) for snapshot in snapshots]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        for result in executor.map(_audit_snapshot, tasks):
            results.append(result)
            print(f"DONE {result['snapshot']}", flush=True)

    records = [record for result in results for record in result["records"]]
    definitions = {}
    for result in results:
        definitions.update(result["axis_definitions"])
    report = {
        "schema_version": 1,
        "kind": "rebuilt_pdf_binning_audit",
        "production_root": str(production_root),
        "analysis_root": str(analysis_root),
        "snapshot_count": len(results),
        "axis_record_count": len(records),
        "snapshot_manifests": [result["manifest"] for result in results],
        "axis_definitions": definitions,
        "aggregates": _aggregate(records),
        "records": records,
    }
    _write_json(output_dir / "binning_audit.json", report)
    (output_dir / "binning_audit.md").write_text(_markdown(report), encoding="utf-8")
    print(f"Wrote audit to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
