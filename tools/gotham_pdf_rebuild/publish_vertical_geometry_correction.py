#!/usr/bin/env python3
"""Publish corrected cylindrical-radius/absolute-height plot axes."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping

from repath_rebuilt_pdf_plot_metadata import repath_value


PRODUCTS = {
    "science_vertical_geometry_edot_out",
    "science_vertical_geometry_mdot_in_abs",
    "science_vertical_geometry_mdot_out",
    "science_vertical_geometry_ram_out",
}
OLD_PREFIX = "absolute_z_cylindrical_radius_"
NEW_PREFIX = "cylindrical_radius_absolute_z_"


def _write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(path)


def _publish_png(source: Path, destination: Path) -> None:
    if destination.exists() and os.path.samefile(source, destination):
        return
    temporary = destination.with_name(destination.name + ".partial")
    temporary.unlink(missing_ok=True)
    try:
        os.link(source, temporary)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("staging_root", type=Path)
    parser.add_argument("--final-root", type=Path, required=True)
    parser.add_argument("--snapshot-count", type=int, required=True)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    staging_root = args.staging_root.resolve()
    final_root = args.final_root.resolve()
    published_pngs = 0
    updated_products = 0
    updated_batches = 0

    for index in range(args.snapshot_count):
        snapshot_root = staging_root / str(index)
        batch_path = snapshot_root / "metadata.json"
        if not batch_path.is_file():
            raise FileNotFoundError(f"Missing correction batch: {batch_path}")
        staged_batch = json.loads(batch_path.read_text())
        if staged_batch["failed_product_count"] != 0:
            raise ValueError(f"Correction snapshot {index} has failures")
        if staged_batch["product_count"] != 4 or staged_batch["plot_count"] != 4:
            raise ValueError(
                f"Correction snapshot {index} has unexpected counts: "
                f"{staged_batch['product_count']} products, "
                f"{staged_batch['plot_count']} plots"
            )
        if {product["product"] for product in staged_batch["products"]} != PRODUCTS:
            raise ValueError(f"Correction snapshot {index} selected wrong products")

        staged_products: Dict[str, Dict[str, Any]] = {}
        for staged_product_path in sorted((snapshot_root / "metadata").glob("product.*.json")):
            staged = json.loads(staged_product_path.read_text())
            product = staged["product"]
            if product not in PRODUCTS:
                raise ValueError(f"Unexpected correction product: {product}")
            rewritten = repath_value(staged, staging_root, final_root)
            plot = rewritten["plots"][0]
            if plot["x_axis"] != "cylindrical_radius" or plot["y_axis"] != "absolute_z":
                raise ValueError(f"Incorrect corrected axes for {product}: {plot}")
            new_path = Path(plot["paths"][0])
            if not new_path.name.startswith(NEW_PREFIX):
                raise ValueError(f"Incorrect corrected filename: {new_path}")

            final_product_path = final_root / "metadata" / staged_product_path.name
            old = json.loads(final_product_path.read_text())
            for old_plot in old["plots"]:
                for old_path_value in old_plot["paths"]:
                    old_name = Path(old_path_value).name
                    if old_name.startswith(OLD_PREFIX):
                        (final_root / "png" / old_name).unlink(missing_ok=True)

            source_png = next((snapshot_root / "png").glob(new_path.name))
            _publish_png(source_png, new_path)
            _write_json_atomic(final_product_path, rewritten)
            staged_products[product] = rewritten
            published_pngs += 1
            updated_products += 1

        staged_batch_metadata = next((snapshot_root / "metadata").glob("batch.*.json"))
        final_batch_path = final_root / "metadata" / staged_batch_metadata.name
        final_batch = json.loads(final_batch_path.read_text())
        final_batch["products"] = [
            staged_products.get(product["product"], product)
            for product in final_batch["products"]
        ]
        _write_json_atomic(final_batch_path, final_batch)
        updated_batches += 1

    manifest_path = final_root / "plot_campaign_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["vertical_geometry_axis_correction"] = {
        "products": sorted(PRODUCTS),
        "snapshot_count": args.snapshot_count,
        "png_count": published_pngs,
        "x_axis": "cylindrical_radius",
        "y_axis": "absolute_z",
    }
    _write_json_atomic(manifest_path, manifest)
    print(
        f"Published {published_pngs} corrected vertical-geometry PNGs, "
        f"{updated_products} product metadata files, and {updated_batches} batches"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
