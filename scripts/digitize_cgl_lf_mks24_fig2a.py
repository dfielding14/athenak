#!/usr/bin/env python3
"""Extract raw donor-coordinate MKS24 Figure 2(a) surfaces from its raster PDF."""

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
    "277545ff1a86b5ca46e97fa034ce189f72f149ecbb515470c30996494c1de129"
)
IMAGE_RE = re.compile(
    rb"<<(.*?)>>\s*stream\r?\n(.*?)\r?\nendstream", re.DOTALL
)
TILE_SIZES = (
    (450, 400), (450, 400), (450, 400), (56, 400),
    (450, 400), (450, 400), (450, 400), (56, 400),
    (450, 400), (450, 400), (450, 400), (56, 400),
    (450, 50), (450, 50), (450, 50), (56, 50),
)
MOSAIC_WIDTH = 1406
MOSAIC_HEIGHT = 1250
COLORBAR_Y = 210
COLORBAR_X = (151, 1263)
COLORBAR_TICKS = (
    (148.0, 0.0),
    (438.0, 10.0),
    (728.0, 20.0),
    (1018.0, 30.0),
)
X_TICKS = (
    (253.0, -0.2),
    (478.0, -0.1),
    (703.0, 0.0),
    (927.0, 0.1),
    (1152.0, 0.2),
)
SAMPLE_STRIDE = 8
SURFACE_Z_UNCERTAINTY = 0.5
PANELS = {
    "parallel": {
        "intended_case": "paper_standard_active_alfvenic_beta10",
        "intended_product": "pressure_density_joint.parallel",
        "bounds": (144, 1261, 253, 652),
        "y_ticks": ((290.0, 1.0), (452.0, 0.0), (615.0, -1.0)),
    },
    "perpendicular": {
        "intended_case": "paper_standard_active_alfvenic_beta10",
        "intended_product": "pressure_density_joint.perpendicular",
        "bounds": (144, 1261, 659, 1057),
        "y_ticks": ((696.0, 1.0), (859.0, 0.0), (1021.0, -1.0)),
    },
}


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pixel(image: bytes, width: int, row: int, column: int) -> bytes:
    """Return one three-byte RGB raster value."""

    offset = (row * width + column) * 3
    return image[offset:offset + 3]


def pdf_images(pdf: bytes) -> list[bytes]:
    """Return the pinned tiled ICC/sRGB raster image streams in object order."""

    images: list[bytes] = []
    sizes: list[tuple[int, int]] = []
    for header, compressed in IMAGE_RE.findall(pdf):
        if b"/Subtype /Image" not in header:
            continue
        if (
            b"/Filter /FlateDecode" not in header
            or b"/ColorSpace [/ICCBased" not in header
        ):
            raise ValueError("Figure 2(a) raster image encoding is not recognized")
        width_match = re.search(rb"/Width\s+(\d+)", header)
        height_match = re.search(rb"/Height\s+(\d+)", header)
        if width_match is None or height_match is None:
            raise ValueError("Figure 2(a) raster image dimensions are absent")
        size = (int(width_match.group(1)), int(height_match.group(1)))
        image = zlib.decompress(compressed)
        if len(image) != size[0] * size[1] * 3:
            raise ValueError(f"Figure 2(a) image byte count is invalid for {size}")
        sizes.append(size)
        images.append(image)
    if tuple(sizes) != TILE_SIZES:
        raise ValueError("Figure 2(a) does not contain the pinned tiled raster layout")
    return images


def tiled_mosaic(images: list[bytes]) -> bytes:
    """Reassemble the source-resolution plotted mosaic from its PDF tiles."""

    rows = (images[12:16], images[8:12], images[4:8], images[0:4])
    output = bytearray()
    for row_tiles in rows:
        tile_height = len(row_tiles[0]) // (TILE_SIZES[0][0] * 3)
        widths = [TILE_SIZES[index][0] for index in (
            (12, 13, 14, 15) if tile_height == 50 else (0, 1, 2, 3)
        )]
        for row in range(tile_height):
            for image, width in zip(row_tiles, widths):
                start = row * width * 3
                output.extend(image[start:start + width * 3])
    if len(output) != MOSAIC_WIDTH * MOSAIC_HEIGHT * 3:
        raise ValueError("Figure 2(a) reconstructed mosaic has invalid dimensions")
    return bytes(output)


def linear_calibration(
    ticks: tuple[tuple[float, float], ...],
) -> tuple[float, float, float]:
    """Fit plotted values as a linear function of a raster coordinate."""

    mean_coordinate = sum(point[0] for point in ticks) / len(ticks)
    mean_value = sum(point[1] for point in ticks) / len(ticks)
    numerator = sum(
        (point[0] - mean_coordinate) * (point[1] - mean_value)
        for point in ticks
    )
    denominator = sum((point[0] - mean_coordinate) ** 2 for point in ticks)
    slope = numerator / denominator
    intercept = mean_value - slope * mean_coordinate
    residual = max(
        abs(intercept + slope * coordinate - value)
        for coordinate, value in ticks
    )
    return intercept, slope, residual


def colorbar_values(mosaic: bytes) -> tuple[dict[bytes, float], dict[str, float]]:
    """Decode plotted colorbar RGB values into joint-PDF density values."""

    intercept, slope, residual = linear_calibration(COLORBAR_TICKS)
    values_by_color: dict[bytes, list[float]] = {}
    for column in range(COLORBAR_X[0], COLORBAR_X[1] + 1):
        rgb = pixel(mosaic, MOSAIC_WIDTH, COLORBAR_Y, column)
        values_by_color.setdefault(rgb, []).append(intercept + slope * column)
    values = {
        rgb: sum(samples) / len(samples)
        for rgb, samples in values_by_color.items()
    }
    if len(values) != 256:
        raise ValueError("Figure 2(a) colorbar does not retain its pinned palette")
    return values, {
        "intercept": intercept,
        "slope": slope,
        "maximum_tick_fit_residual": residual,
        "distinct_rgb_values": len(values),
        "minimum_decoded_density": min(values.values()),
        "maximum_decoded_density": max(values.values()),
    }


def sample_surface(
    mosaic: bytes, panel: dict[str, object], values: dict[bytes, float]
) -> tuple[list[tuple[float, float, float]], dict[str, object]]:
    """Sample one plotted surface while dropping overlaid guide-line pixels."""

    x_intercept, x_slope, x_residual = linear_calibration(X_TICKS)
    y_intercept, y_slope, y_residual = linear_calibration(panel["y_ticks"])
    x_min, x_max, y_min, y_max = panel["bounds"]
    non_colorbar = {
        pixel(mosaic, MOSAIC_WIDTH, row, column)
        for row in range(y_min, y_max + 1)
        for column in range(x_min, x_max + 1)
        if pixel(mosaic, MOSAIC_WIDTH, row, column) not in values
    }
    recognized_overlays = {b"\x00\x00\x00", b"\x26\x26\x26"}
    if non_colorbar.difference(recognized_overlays):
        raise ValueError(
            "Figure 2(a) panel has non-colorbar pixels other than black overlays"
        )
    samples: list[tuple[float, float, float]] = []
    omitted = 0
    for row in range(y_min + SAMPLE_STRIDE // 2, y_max + 1, SAMPLE_STRIDE):
        for column in range(
            x_min + SAMPLE_STRIDE // 2, x_max + 1, SAMPLE_STRIDE
        ):
            rgb = pixel(mosaic, MOSAIC_WIDTH, row, column)
            if rgb not in values:
                omitted += 1
                continue
            samples.append((
                x_intercept + x_slope * column,
                y_intercept + y_slope * row,
                values[rgb],
            ))
    if len(samples) < 100:
        raise ValueError("Figure 2(a) surface does not contain enough valid samples")
    return samples, {
        "x_tick_fit_residual": x_residual,
        "y_tick_fit_residual": y_residual,
        "recognized_non_colorbar_overlay_rgb": ["#000000", "#262626"],
        "sample_count": len(samples),
        "omitted_non_colorbar_samples": omitted,
    }


def write_surfaces(source_figure: Path, output_dir: Path) -> Path:
    """Write donor-coordinate CSV surfaces and an explicitly excluded record."""

    digest = sha256_file(source_figure)
    if digest != EXPECTED_FIGURE_SHA256:
        raise ValueError(
            "Figure 2(a) SHA-256 does not match pinned arXiv:2405.02418v2 source"
        )
    mosaic = tiled_mosaic(pdf_images(source_figure.read_bytes()))
    values, colorbar = colorbar_values(mosaic)
    output_dir.mkdir(parents=True, exist_ok=True)
    archived_figure = output_dir / "fig2a.pdf"
    shutil.copy2(source_figure, archived_figure)
    entries: list[dict[str, object]] = []
    panel_metadata: dict[str, object] = {}
    for panel_name, panel in PANELS.items():
        samples, extraction = sample_surface(mosaic, panel, values)
        csv_path = output_dir / f"fig2a_{panel_name}_paper_pressure_units.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("paper_x", "paper_y", "z", "z_uncertainty"))
            for x, y, z in samples:
                writer.writerow((
                    f"{x:.17g}",
                    f"{y:.17g}",
                    f"{z:.17g}",
                    f"{SURFACE_Z_UNCERTAINTY:.17g}",
                ))
        entries.append({
            "id": f"fig2a_{panel_name}_paper_pressure_units",
            "intended_case": panel["intended_case"],
            "intended_product": panel["intended_product"],
            "data_file": csv_path.name,
            "data_sha256": sha256_file(csv_path),
            "sample_count": extraction["sample_count"],
        })
        panel_metadata[panel_name] = {
            "bounds": list(panel["bounds"]),
            "y_ticks": [list(point) for point in panel["y_ticks"]],
            **extraction,
        }
    record = {
        "admission_status": "excluded_pending_paper_to_athenak_pressure_transform",
        "source_description": (
            "MKS24 arXiv:2405.02418v2 source/fig2a.pdf; source-resolution "
            "tiled raster decoded through its labeled linear colorbar and axes"
        ),
        "source_figure": archived_figure.name,
        "source_figure_sha256": digest,
        "digitization_tool": Path(__file__).name,
        "coordinate_units": "MKS24 plotted donor pressure units",
        "exclusion_reason": (
            "Both plotted axes carry donor pressure units. The current "
            "AthenaK normalized deck has no qualified donor-pressure conversion; "
            "if a pressure scale s is later established, conversion requires "
            "x_A=s*x_paper, y_A=s*y_paper, z_A=z_paper/s^2, and "
            "sigma_z_A=sigma_z_paper/s^2."
        ),
        "uncertainty_description": (
            "Absolute joint-PDF density uncertainty 0.5 conservatively covers "
            "color quantization and raster sampling. Black dotted-guide pixels "
            "and axis-overlay pixels absent from the colorbar palette are omitted."
        ),
        "surfaces": entries,
    }
    record_path = output_dir / "excluded_surfaces.json"
    record_path.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    metadata = {
        "source_figure_sha256": digest,
        "mosaic_size": [MOSAIC_WIDTH, MOSAIC_HEIGHT],
        "tile_sizes_in_pdf_order": [list(size) for size in TILE_SIZES],
        "colorbar": {
            "raster_y": COLORBAR_Y,
            "raster_x": list(COLORBAR_X),
            "x_value_ticks": [list(point) for point in COLORBAR_TICKS],
            **colorbar,
        },
        "x_ticks": [list(point) for point in X_TICKS],
        "sample_stride_pixels": SAMPLE_STRIDE,
        "surface_z_uncertainty": SURFACE_Z_UNCERTAINTY,
        "panels": panel_metadata,
        "excluded_surface_record": record_path.name,
        "excluded_surface_record_sha256": sha256_file(record_path),
    }
    (output_dir / "digitization_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return record_path


def main() -> int:
    """Command-line entry point."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_figure", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    record = write_surfaces(args.source_figure, args.output_dir)
    print(f"Wrote excluded Figure 2(a) surface record to {record}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
