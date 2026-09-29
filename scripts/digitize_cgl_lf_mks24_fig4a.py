#!/usr/bin/env python3
"""Extract normalized-density PDF curves from MKS24 Figure 4(a)."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import zlib


EXPECTED_FIGURE_SHA256 = (
    "18a5c7119f2cb316e0283b354d25479d006aa3f27df33031035b16810937c024"
)
PLOT = (60.0, 362.0, 26.0, 305.0)
X_LINEAR_TICKS = ((69.14786530, 0.8), (212.45257568, 1.0))
Y_LOG_TICKS = ((36.76585007, 100.0), (163.52516174, 1.0))
PDF_RELATIVE_Y_UNCERTAINTY = 0.05
GROUP_RE = re.compile(rb"q\n6 0 0 6 1248 2118 cm\n(.*?)\nQ", re.DOTALL)
POINT_RE = re.compile(
    r"^([-+]?[0-9]*\.?[0-9]+) ([-+]?[0-9]*\.?[0-9]+) ([ml])$",
    re.MULTILINE,
)
COLOR_RE = re.compile(r"^([0-9.]+ [0-9.]+ [0-9.]+) RG$", re.MULTILINE)
CURVES = {
    ("0 0.4471 0.7412", "solid"): (
        "fig4a_density_active_alfvenic_beta10",
        "paper_standard_active_alfvenic_beta10",
    ),
    ("0.851 0.3255 0.098", "solid"): (
        "fig4a_density_active_alfvenic_beta100",
        "paper_standard_active_alfvenic_beta100",
    ),
    ("0.9294 0.6941 0.1255", "solid"): (
        "fig4a_density_active_random_beta10",
        "paper_standard_active_random_beta10",
    ),
    ("0.4941 0.1843 0.5569", "solid"): (
        "fig4a_density_active_random_beta100",
        "paper_standard_active_random_beta100",
    ),
    ("0 0.4471 0.7412", "dashed"): (
        "fig4a_density_passive_alfvenic_beta10",
        "paper_standard_passive_alfvenic_beta10",
    ),
    ("0.851 0.3255 0.098", "dashed"): (
        "fig4a_density_passive_alfvenic_beta100",
        "paper_standard_passive_alfvenic_beta100",
    ),
    ("0.9294 0.6941 0.1255", "dashed"): (
        "fig4a_density_passive_random_beta10",
        "paper_standard_passive_random_beta10",
    ),
    ("0.4941 0.1843 0.5569", "dashed"): (
        "fig4a_density_passive_random_beta100",
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
    """Return the vector-content stream containing Figure 4(a)."""

    for compressed in re.findall(
        rb"stream\r?\n(.*?)\r?\nendstream", pdf, re.DOTALL
    ):
        try:
            content = zlib.decompress(compressed)
        except zlib.error:
            continue
        if b"60 305 m" in content and b"0 0.4471 0.7412 RG" in content:
            return content
    raise ValueError("did not locate the Figure 4(a) vector-content stream")


def map_coordinate(raw_x: float, raw_y: float) -> tuple[float, float]:
    """Map plotted axes to density fluctuation and probability density."""

    (x0, value0), (x1, value1) = X_LINEAR_TICKS
    density_ratio = value0 + (raw_x - x0) * (value1 - value0) / (x1 - x0)
    (y0, yvalue0), (y1, yvalue1) = Y_LOG_TICKS
    probability = yvalue0 * 10.0 ** (
        (raw_y - y0) * math.log10(yvalue1 / yvalue0) / (y1 - y0)
    )
    return density_ratio - 1.0, probability


def extract_curves(
    content: bytes,
) -> tuple[dict[tuple[str, str], list[tuple[float, float]]], int]:
    """Extract visible vector paths and omit boundary-clipped vertices."""

    extracted = {key: [] for key in CURVES}
    clipped_vertices = 0
    xmin, xmax, ymin, ymax = PLOT
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
        if len(points) < 3:
            continue
        for raw_x, raw_y, _ in points:
            x_coordinate = float(raw_x)
            y_coordinate = float(raw_y)
            if (
                x_coordinate <= xmin or x_coordinate >= xmax
                or y_coordinate <= ymin or y_coordinate >= ymax
            ):
                clipped_vertices += 1
                continue
            extracted[key].append(map_coordinate(x_coordinate, y_coordinate))
    for key, points in extracted.items():
        points.sort()
        if len(points) < 2:
            raise ValueError(f"no usable plotted vector curve found for {key}")
        if any(x1 >= x2 for (x1, _), (x2, _) in zip(points, points[1:])):
            raise ValueError(f"non-increasing coordinates extracted for {key}")
        if any(not math.isfinite(x) or not math.isfinite(y) for x, y in points):
            raise ValueError(f"nonfinite coordinates extracted for {key}")
    return extracted, clipped_vertices


def write_curves(source_figure: Path, output_dir: Path) -> Path:
    """Write admitted density CSVs, copied figure, metadata, and manifest."""

    digest = sha256_file(source_figure)
    if digest != EXPECTED_FIGURE_SHA256:
        raise ValueError(
            "Figure 4(a) SHA-256 does not match pinned arXiv:2405.02418v2 source"
        )
    curves, clipped_vertices = extract_curves(pdf_content(source_figure.read_bytes()))
    output_dir.mkdir(parents=True, exist_ok=True)
    archived_figure = output_dir / "fig4a.pdf"
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
                    (
                        f"{x:.17g}",
                        f"{y:.17g}",
                        f"{PDF_RELATIVE_Y_UNCERTAINTY * y:.17g}",
                    )
                )
        point_counts[curve_id] = len(curves[key])
        entries.append({
            "id": curve_id,
            "case": case,
            "product": "pdf.density_fluctuation",
            "data_file": csv_path.name,
            "data_sha256": sha256_file(csv_path),
            "interpolation": "linear",
        })
    manifest = {
        "schema_version": 1,
        "provenance": {
            "method": "digitized",
            "source_description": (
                "MKS24 arXiv:2405.02418v2 source/fig4a.pdf; vector "
                "polyline vertices extracted through tick-calibrated axes "
                "and mapped from rho/<rho> to rho/<rho>-1"
            ),
            "source_figure": archived_figure.name,
            "source_figure_sha256": digest,
            "digitization_tool": Path(__file__).name,
            "uncertainty_description": (
                "Five-percent relative probability-density uncertainty "
                "conservatively covers plotted stroke width and "
                "tick-calibrated coordinate mapping; vertices clipped at "
                "the plotted boundaries are omitted. The panel's explicit "
                "rho/<rho> coordinate maps to the analyzer's "
                "density_fluctuation product by subtracting one."
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
            "x_linear_ticks_rho_over_mean_rho": [
                list(tick) for tick in X_LINEAR_TICKS
            ],
            "x_product_transform": "density_fluctuation = rho/<rho> - 1",
            "y_log_ticks": [list(tick) for tick in Y_LOG_TICKS],
        },
        "relative_y_uncertainty": {
            "pdf.density_fluctuation": PDF_RELATIVE_Y_UNCERTAINTY,
        },
        "point_counts": point_counts,
        "clipped_vertices_omitted": clipped_vertices,
        "excluded_curves": {
            "fig4b_density_spectrum": (
                "The companion panel labels E_rho without an explicit "
                "rho/<rho>-1 spectrum normalization; it is not emitted "
                "pending a qualified spectral ordinate transform."
            ),
        },
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
    print(f"Wrote digitized Figure 4(a) curve manifest to {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
