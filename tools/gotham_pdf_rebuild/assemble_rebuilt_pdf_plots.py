#!/usr/bin/env python3
"""Audit staged rebuilt-PDF plots and publish one shallow production tree."""

from __future__ import annotations

import argparse
import json
import os
import re
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Tuple

from repath_rebuilt_pdf_plot_metadata import repath_value


OUTPUT_LABEL_RE = re.compile(r"(?:^|[._-])output\d+(?:$|[._-])")


def _write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(path)


def _link_unique(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(source, destination)
    except FileExistsError as exc:
        raise FileExistsError(f"Campaign filename collision: {destination}") from exc


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("staging_root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--snapshot-count", type=int, required=True)
    parser.add_argument("--expected-products", type=int, default=64)
    parser.add_argument("--expected-plots", type=int, default=239)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    staging_root = args.staging_root.resolve()
    output_dir = args.output_dir.resolve()
    if not staging_root.is_dir():
        raise SystemExit(f"Staging root does not exist: {staging_root}")
    if any(output_dir.iterdir()) if output_dir.exists() else False:
        raise SystemExit(f"Final output directory is not empty: {output_dir}")

    snapshots = []
    png_sources: Dict[str, Path] = {}
    metadata_sources: Dict[str, Tuple[Path, Path]] = {}
    series_counts: Counter[str] = Counter()
    for index in range(args.snapshot_count):
        snapshot_root = staging_root / str(index)
        metadata_path = snapshot_root / "metadata.json"
        if not metadata_path.is_file():
            raise FileNotFoundError(f"Missing staged snapshot metadata: {metadata_path}")
        metadata = json.loads(metadata_path.read_text())
        if metadata["failed_product_count"] != 0:
            raise ValueError(f"Snapshot {index} has plot failures: {metadata['failures']}")
        if metadata["product_count"] != args.expected_products:
            raise ValueError(
                f"Snapshot {index} has {metadata['product_count']} products; "
                f"expected {args.expected_products}"
            )
        if metadata["plot_count"] != args.expected_plots:
            raise ValueError(
                f"Snapshot {index} has {metadata['plot_count']} plots; "
                f"expected {args.expected_plots}"
            )

        pngs = sorted((snapshot_root / "png").glob("*.png"))
        if len(pngs) != args.expected_plots:
            raise ValueError(
                f"Snapshot {index} has {len(pngs)} PNGs; expected {args.expected_plots}"
            )
        for path in pngs:
            if OUTPUT_LABEL_RE.search(path.name):
                raise ValueError(f"PNG has forbidden outputXX label: {path}")
            if path.name in png_sources:
                raise ValueError(
                    f"Duplicate campaign PNG name {path.name}: "
                    f"{png_sources[path.name]} and {path}"
                )
            png_sources[path.name] = path

        for path in sorted((snapshot_root / "metadata").glob("*.json")):
            if path.name in metadata_sources:
                raise ValueError(
                    f"Duplicate campaign metadata name {path.name}: "
                    f"{metadata_sources[path.name]} and {path}"
                )
            metadata_sources[path.name] = (path, snapshot_root)

        identity = metadata["snapshot_identity"]
        series_counts[str(identity["simulation"])] += 1
        snapshots.append(
            {
                "index": index,
                "input_dir": metadata["input_dir"],
                "simulation": identity["simulation"],
                "phase": identity["phase"],
                "time_token": identity["time_token"],
                "time_title": identity["time_title"],
                "product_count": metadata["product_count"],
                "plot_count": metadata["plot_count"],
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    for name, source in png_sources.items():
        _link_unique(source, output_dir / "png" / name)
    for name, (source, snapshot_root) in metadata_sources.items():
        value = json.loads(source.read_text())
        rewritten = repath_value(value, snapshot_root, output_dir)
        _write_json_atomic(output_dir / "metadata" / name, rewritten)

    manifest = {
        "schema_version": 1,
        "kind": "rebuilt_pdf_plot_campaign",
        "staging_root": str(staging_root),
        "output_dir": str(output_dir),
        "snapshot_count": len(snapshots),
        "product_count": len(snapshots) * args.expected_products,
        "png_count": len(png_sources),
        "metadata_count": len(metadata_sources),
        "series_counts": dict(sorted(series_counts.items())),
        "forbidden_output_label_png_count": 0,
        "snapshots": snapshots,
    }
    _write_json_atomic(output_dir / "plot_campaign_manifest.json", manifest)
    print(
        f"Published {len(png_sources)} descriptive PNGs from {len(snapshots)} "
        f"snapshots to {output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
