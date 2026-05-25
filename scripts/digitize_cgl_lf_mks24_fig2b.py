#!/usr/bin/env python3
"""Extract MKS24 Figure 2(b) unstable-fraction curves from its vector PDF."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import zlib


EXPECTED_FIGURE_SHA256 = (
    "828287244aa0ec72aad6e802268fbcb8d77af65f11b03fcf9634eb9cf741ab0e"
)
PLOT_X_MIN = 59.0
PLOT_X_MAX = 407.0
PLOT_Y_MIN = 30.0
PLOT_Y_MAX = 351.0
X_DATA_MAX = 10.0
Y_DATA_MAX = 0.8
Y_UNCERTAINTY = 0.0025
GROUP_RE = re.compile(rb"q\n6 0 0 6 1098 1968 cm\n(.*?)\nQ", re.DOTALL)
POINT_RE = re.compile(
    r"^([-+]?[0-9]*\.?[0-9]+) ([-+]?[0-9]*\.?[0-9]+) ([ml])$",
    re.MULTILINE,
)
COLOR_RE = re.compile(r"^([0-9.]+ [0-9.]+ [0-9.]+) RG$", re.MULTILINE)
CURVES = {
    ("0 0.4471 0.7412", "solid"): (
        "fig2b_active_alfvenic_beta10",
        "paper_standard_active_alfvenic_beta10",
    ),
    ("0.851 0.3255 0.098", "solid"): (
        "fig2b_active_alfvenic_beta100",
        "paper_standard_active_alfvenic_beta100",
    ),
    ("0.9294 0.6941 0.1255", "solid"): (
        "fig2b_active_random_beta10",
        "paper_standard_active_random_beta10",
    ),
    ("0.4941 0.1843 0.5569", "solid"): (
        "fig2b_active_random_beta100",
        "paper_standard_active_random_beta100",
    ),
    ("0 0.4471 0.7412", "dashed"): (
        "fig2b_passive_alfvenic_beta10",
        "paper_standard_passive_alfvenic_beta10",
    ),
    ("0.851 0.3255 0.098", "dashed"): (
        "fig2b_passive_alfvenic_beta100",
        "paper_standard_passive_alfvenic_beta100",
    ),
    ("0.9294 0.6941 0.1255", "dashed"): (
        "fig2b_passive_random_beta10",
        "paper_standard_passive_random_beta10",
    ),
    ("0.4941 0.1843 0.5569", "dashed"): (
        "fig2b_passive_random_beta100",
        "paper_standard_passive_random_beta100",
    ),
}


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pdf_content(pdf: bytes) -> bytes:
    """Return the decompressed vector-content stream containing the axes."""

    for compressed in re.findall(
        rb"stream\r?\n(.*?)\r?\nendstream", pdf, re.DOTALL
    ):
        try:
            content = zlib.decompress(compressed)
        except zlib.error:
            continue
        if b"59 351 m" in content and b"0 0.4471 0.7412 RG" in content:
            return content
    raise ValueError("did not locate the Figure 2(b) vector-content stream")


def extract_curves(
    content: bytes,
) -> tuple[dict[tuple[str, str], list[tuple[float, float]]], dict[str, int]]:
    """Extract visible vector vertices and convert axis coordinates to data."""

    extracted = {key: [] for key in CURVES}
    clipped_vertices = 0
    for raw_group in GROUP_RE.findall(content):
        group = raw_group.decode("ascii", errors="ignore")
        color = COLOR_RE.search(group)
        if color is None:
            continue
        style = "dashed" if "[10 6] 0 d" in group else "solid"
        key = (color.group(1), style)
        if key not in CURVES:
            continue
        points = POINT_RE.findall(group)
        # Legend samples contain two vertices; plotted series contain at least three.
        if len(points) < 3:
            continue
        for raw_x, raw_y, _ in points:
            x = (float(raw_x) - PLOT_X_MIN) * X_DATA_MAX / (
                PLOT_X_MAX - PLOT_X_MIN
            )
            y = (PLOT_Y_MAX - float(raw_y)) * Y_DATA_MAX / (
                PLOT_Y_MAX - PLOT_Y_MIN
            )
            if x < -1.0e-8 or x > X_DATA_MAX + 1.0e-8:
                continue
            if y < -1.0e-8 or y > Y_DATA_MAX + 1.0e-8:
                clipped_vertices += 1
                continue
            extracted[key].append(
                (max(0.0, min(X_DATA_MAX, x)), max(0.0, min(Y_DATA_MAX, y)))
            )
    for key, points in extracted.items():
        points.sort()
        if len(points) < 2:
            raise ValueError(f"no usable plotted vector curve found for {key}")
        if any(x1 >= x2 for (x1, _), (x2, _) in zip(points, points[1:])):
            raise ValueError(f"non-increasing coordinates extracted for {key}")
    return extracted, {"clipped_vertices_omitted": clipped_vertices}


def write_curves(source_figure: Path, output_dir: Path) -> Path:
    """Write CSVs, copied figure, metadata, and analyzer-compatible manifest."""

    digest = sha256_file(source_figure)
    if digest != EXPECTED_FIGURE_SHA256:
        raise ValueError(
            "Figure 2(b) SHA-256 does not match pinned arXiv:2405.02418v2 source"
        )
    curves, extraction = extract_curves(pdf_content(source_figure.read_bytes()))
    output_dir.mkdir(parents=True, exist_ok=True)
    archived_figure = output_dir / "fig2b.pdf"
    shutil.copy2(source_figure, archived_figure)
    entries: list[dict[str, str]] = []
    point_counts: dict[str, int] = {}
    for key, (curve_id, case) in CURVES.items():
        csv_path = output_dir / f"{curve_id}.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("x", "y", "y_uncertainty"))
            for x, y in curves[key]:
                writer.writerow(
                    (f"{x:.10f}", f"{y:.10f}", f"{Y_UNCERTAINTY:.10f}")
                )
        point_counts[curve_id] = len(curves[key])
        entries.append(
            {
                "id": curve_id,
                "case": case,
                "product": "history.unstable_fraction",
                "data_file": csv_path.name,
                "data_sha256": sha256_file(csv_path),
                "interpolation": "linear",
            }
        )
    manifest = {
        "schema_version": 1,
        "provenance": {
            "method": "digitized",
            "source_description": (
                "MKS24 arXiv:2405.02418v2 source/fig2b.pdf; vector polyline "
                "vertices extracted in plotted axis coordinates"
            ),
            "source_figure": archived_figure.name,
            "source_figure_sha256": digest,
            "digitization_tool": Path(__file__).name,
            "uncertainty_description": (
                "Absolute unstable-fraction uncertainty 0.0025 at each visible "
                "vector vertex, conservatively covering plotted stroke width and "
                "coordinate mapping; vertices clipped above y=0.8 are omitted."
            ),
        },
        "curves": entries,
    }
    manifest_path = output_dir / "curves.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    metadata = {
        "source_figure_sha256": digest,
        "axis_mapping": {
            "plot_x": [PLOT_X_MIN, PLOT_X_MAX],
            "plot_y": [PLOT_Y_MIN, PLOT_Y_MAX],
            "data_x": [0.0, X_DATA_MAX],
            "data_y": [0.0, Y_DATA_MAX],
        },
        "y_uncertainty": Y_UNCERTAINTY,
        "point_counts": point_counts,
        **extraction,
        "curve_manifest": manifest_path.name,
        "curve_manifest_sha256": sha256_file(manifest_path),
    }
    (output_dir / "digitization_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest_path


def main() -> int:
    """Command-line entry point."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_figure", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    manifest = write_curves(args.source_figure, args.output_dir)
    print(f"Wrote digitized Figure 2(b) curve manifest to {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
