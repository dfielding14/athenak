#!/usr/bin/env python3
"""Extract selected-shell MKS24 Figure 8 alignment PDFs from its raster panels."""

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
    "0e3e292406bc27ce24b2ef449bd95cb8addf240f02b5ac60f60286459a299d8b"
)
PLOT_Y_TOP = 22.0
PLOT_Y_BOTTOM = 244.0
PDF_BINS = 64
ALIGNMENT_SHELLS = (8, 16, 32, 64, 128)
PDF_DENSITY_UNCERTAINTY = 0.05
NORMALIZATION_TOLERANCE = 0.03
COLORBAR_TICKS = (
    (225.18644714, 0.4),
    (187.55932617, 0.6),
    (149.93220520, 0.8),
    (112.30508423, 1.0),
    (74.67796326, 1.2),
    (37.05084229, 1.4),
)
IMAGE_RE = re.compile(
    rb"<<(.*?)>>\s*stream\r?\n(.*?)\r?\nendstream", re.DOTALL
)
PANELS = {
    "active": {
        "image_size": (1616, 1332),
        "case": "paper_standard_active_alfvenic_beta10",
        "curve_prefix": "fig8_alignment_pdf_active_alfvenic_beta10",
        "plot_x": (90.96793365, 269.064148),
        "x_log_tick": (215.02635193, 100.0),
        "x_log_decade_width": (
            (271.80395508 - 215.02635193) / math.log10(2.0)
        ),
    },
    "passive": {
        "image_size": (1610, 1332),
        "case": "paper_standard_passive_alfvenic_beta10",
        "curve_prefix": "fig8_alignment_pdf_passive_alfvenic_beta10",
        "plot_x": (359.96804810, 268.063873),
        "x_log_tick": (483.56527710, 100.0),
        "x_log_decade_width": (
            (540.13183594 - 483.56527710) / math.log10(2.0)
        ),
    },
}
COLORBAR_IMAGE_SIZE = (96, 1332)


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pdf_images(pdf: bytes) -> dict[tuple[int, int], bytes]:
    """Extract DeviceRGB Flate image streams indexed by their dimensions."""

    images: dict[tuple[int, int], bytes] = {}
    for header, compressed in IMAGE_RE.findall(pdf):
        if b"/Subtype /Image" not in header:
            continue
        if (
            b"/Filter /FlateDecode" not in header
            or b"/ColorSpace /DeviceRGB" not in header
        ):
            raise ValueError("Figure 8 raster image encoding is not recognized")
        width_match = re.search(rb"/Width\s+(\d+)", header)
        height_match = re.search(rb"/Height\s+(\d+)", header)
        if width_match is None or height_match is None:
            raise ValueError("Figure 8 raster image dimensions are absent")
        size = (int(width_match.group(1)), int(height_match.group(1)))
        image = zlib.decompress(compressed)
        if len(image) != size[0] * size[1] * 3:
            raise ValueError(f"Figure 8 image byte count is invalid for {size}")
        if size in images:
            raise ValueError(f"Figure 8 contains duplicate image dimensions {size}")
        images[size] = image
    expected_sizes = {COLORBAR_IMAGE_SIZE}
    expected_sizes.update(
        panel["image_size"] for panel in PANELS.values()
    )
    if expected_sizes.difference(images):
        raise ValueError("Figure 8 does not contain all pinned raster panels")
    return images


def pixel(image: bytes, width: int, row: int, column: int) -> bytes:
    """Return one three-byte RGB raster value."""

    offset = (row * width + column) * 3
    return image[offset:offset + 3]


def linear_calibration() -> tuple[float, float, float]:
    """Fit the plotted colorbar tick values as a function of PDF ordinate."""

    mean_y = math.fsum(point[0] for point in COLORBAR_TICKS) / len(COLORBAR_TICKS)
    mean_value = (
        math.fsum(point[1] for point in COLORBAR_TICKS) / len(COLORBAR_TICKS)
    )
    numerator = math.fsum(
        (point[0] - mean_y) * (point[1] - mean_value)
        for point in COLORBAR_TICKS
    )
    denominator = math.fsum(
        (point[0] - mean_y) ** 2 for point in COLORBAR_TICKS
    )
    slope = numerator / denominator
    intercept = mean_value - slope * mean_y
    residual = max(
        abs(intercept + slope * point[0] - point[1])
        for point in COLORBAR_TICKS
    )
    return intercept, slope, residual


def colorbar_values(
    colorbar: bytes, width: int, height: int
) -> tuple[dict[bytes, float], dict[str, float]]:
    """Map each exact colorbar RGB value onto its tick-calibrated density."""

    intercept, slope, residual = linear_calibration()
    rows_by_color: dict[bytes, list[int]] = {}
    for row in range(height):
        rgb = pixel(colorbar, width, row, 0)
        if any(pixel(colorbar, width, row, column) != rgb
               for column in range(1, width)):
            raise ValueError("Figure 8 colorbar is not horizontally uniform")
        rows_by_color.setdefault(rgb, []).append(row)
    value_by_color: dict[bytes, float] = {}
    for rgb, rows in rows_by_color.items():
        values = []
        for row in rows:
            pdf_y = PLOT_Y_TOP + (row + 0.5) * (
                PLOT_Y_BOTTOM - PLOT_Y_TOP
            ) / height
            values.append(intercept + slope * pdf_y)
        value_by_color[rgb] = math.fsum(values) / len(values)
    return value_by_color, {
        "intercept": intercept,
        "slope": slope,
        "maximum_tick_fit_residual": residual,
        "distinct_rgb_values": len(value_by_color),
    }


def verify_panel_colors(image: bytes, values: dict[bytes, float]) -> int:
    """Require every heatmap RGB sample to occur in the pinned colorbar."""

    distinct_colors = {image[offset:offset + 3]
                       for offset in range(0, len(image), 3)}
    if distinct_colors.difference(values):
        raise ValueError("Figure 8 panel contains a color absent from its colorbar")
    return len(distinct_colors)


def bilinear_density(
    image: bytes, width: int, height: int, column: float, row: float,
    values: dict[bytes, float]
) -> float:
    """Sample calibrated density from a heatmap at fractional pixel position."""

    if not (0.0 <= column <= width - 1 and 0.0 <= row <= height - 1):
        raise ValueError("Figure 8 requested sample lies outside a raster panel")
    column0 = int(column)
    row0 = int(row)
    column1 = min(column0 + 1, width - 1)
    row1 = min(row0 + 1, height - 1)
    column_fraction = column - column0
    row_fraction = row - row0

    def value(sample_row: int, sample_column: int) -> float:
        return values[pixel(image, width, sample_row, sample_column)]

    lower = (
        (1.0 - column_fraction) * value(row0, column0)
        + column_fraction * value(row0, column1)
    )
    upper = (
        (1.0 - column_fraction) * value(row1, column0)
        + column_fraction * value(row1, column1)
    )
    return (1.0 - row_fraction) * lower + row_fraction * upper


def sample_panel(
    panel: dict[str, object], image: bytes, values: dict[bytes, float]
) -> tuple[dict[int, list[tuple[float, float]]], dict[int, float]]:
    """Extract and normalize the selected-shell PDF slices in one panel."""

    width, height = panel["image_size"]
    plot_x, plot_width = panel["plot_x"]
    x_log_tick, x_log_tick_value = panel["x_log_tick"]
    decade_width = panel["x_log_decade_width"]
    curves: dict[int, list[tuple[float, float]]] = {}
    raw_integrals: dict[int, float] = {}
    for shell in ALIGNMENT_SHELLS:
        k_perp = shell * math.pi
        pdf_x = x_log_tick + math.log10(k_perp / x_log_tick_value) * decade_width
        column = (pdf_x - plot_x) * width / plot_width - 0.5
        sampled = []
        for index in range(PDF_BINS):
            cosine = (index + 0.5) / PDF_BINS
            row = (1.0 - cosine) * height - 0.5
            sampled.append(
                (cosine, bilinear_density(image, width, height, column, row, values))
            )
        integral = math.fsum(point[1] for point in sampled) / PDF_BINS
        if abs(integral - 1.0) > NORMALIZATION_TOLERANCE:
            raise ValueError(
                f"Figure 8 shell {shell} decoded integral {integral} "
                "does not satisfy the published per-shell normalization"
            )
        raw_integrals[shell] = integral
        curves[shell] = [(x, y / integral) for x, y in sampled]
    return curves, raw_integrals


def write_curves(source_figure: Path, output_dir: Path) -> Path:
    """Write selected-shell CSVs, copied figure, metadata, and manifest."""

    digest = sha256_file(source_figure)
    if digest != EXPECTED_FIGURE_SHA256:
        raise ValueError(
            "Figure 8 SHA-256 does not match pinned arXiv:2405.02418v2 source"
        )
    images = pdf_images(source_figure.read_bytes())
    colorbar = images[COLORBAR_IMAGE_SIZE]
    values, calibration = colorbar_values(colorbar, *COLORBAR_IMAGE_SIZE)
    output_dir.mkdir(parents=True, exist_ok=True)
    archived_figure = output_dir / "fig8.pdf"
    shutil.copy2(source_figure, archived_figure)
    entries: list[dict[str, str]] = []
    raw_integrals: dict[str, dict[str, float]] = {}
    distinct_panel_colors: dict[str, int] = {}
    for panel_name, panel in PANELS.items():
        size = panel["image_size"]
        image = images[size]
        distinct_panel_colors[panel_name] = verify_panel_colors(image, values)
        curves, integrals = sample_panel(panel, image, values)
        raw_integrals[panel_name] = {
            str(shell): integral for shell, integral in integrals.items()
        }
        for shell in ALIGNMENT_SHELLS:
            curve_id = f"{panel['curve_prefix']}_shell{shell}"
            csv_path = output_dir / f"{curve_id}.csv"
            with csv_path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(("x", "y", "y_uncertainty"))
                for x, y in curves[shell]:
                    writer.writerow((
                        f"{x:.17g}",
                        f"{y:.17g}",
                        f"{PDF_DENSITY_UNCERTAINTY:.17g}",
                    ))
            entries.append({
                "id": curve_id,
                "case": panel["case"],
                "product": f"alignment.{shell}",
                "data_file": csv_path.name,
                "data_sha256": sha256_file(csv_path),
                "interpolation": "linear",
            })
    manifest = {
        "schema_version": 1,
        "provenance": {
            "method": "digitized",
            "source_description": (
                "MKS24 arXiv:2405.02418v2 source/fig8.pdf; embedded RGB "
                "alignment heatmaps decoded through their linear labeled "
                "colorbar and sampled at selected k_perp shells"
            ),
            "source_figure": archived_figure.name,
            "source_figure_sha256": digest,
            "digitization_tool": Path(__file__).name,
            "uncertainty_description": (
                "Absolute probability-density uncertainty 0.05 covers RGB "
                "color quantization, tick-calibrated color mapping, raster "
                "slice interpolation, and the at-most 0.03 pre-correction "
                "unit-integral discrepancy. Each extracted slice is "
                "renormalized to unit integral as required by the published "
                "per-k_perp-bin normalization."
            ),
        },
        "curves": entries,
    }
    manifest_path = output_dir / "curves.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    maximum_normalization_error = max(
        abs(integral - 1.0)
        for panel_integrals in raw_integrals.values()
        for integral in panel_integrals.values()
    )
    metadata = {
        "source_figure_sha256": digest,
        "colorbar": {
            "pdf_y_value_ticks": [list(tick) for tick in COLORBAR_TICKS],
            **calibration,
        },
        "panel_geometry": {
            name: {
                "image_size": list(panel["image_size"]),
                "plot_x": list(panel["plot_x"]),
                "x_log_tick": list(panel["x_log_tick"]),
                "x_log_decade_width": panel["x_log_decade_width"],
                "distinct_rgb_values": distinct_panel_colors[name],
            }
            for name, panel in PANELS.items()
        },
        "plot_y": [PLOT_Y_TOP, PLOT_Y_BOTTOM],
        "pdf_bins": PDF_BINS,
        "recommended_alignment_shells": list(ALIGNMENT_SHELLS),
        "pdf_density_uncertainty": PDF_DENSITY_UNCERTAINTY,
        "normalization_tolerance": NORMALIZATION_TOLERANCE,
        "unnormalized_slice_integrals": raw_integrals,
        "maximum_unnormalized_integral_error": maximum_normalization_error,
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
    print(f"Wrote digitized Figure 8 curve manifest to {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
