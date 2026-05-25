#!/usr/bin/env python3
"""Extract the dimensionless MKS24 Figure 7 pressure-transfer curves."""

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
    "5411da17391d19fceb26f617f91b4fb1153c7638ccca8afee6c8e0d013651afa"
)
PLOT = (91.0, 634.0, 187.0, 344.0)
X_LOG_TICKS = ((154.19473267, 10.0), (562.1607666, 100.0))
Y_LINEAR_TICKS = ((225.85353088, 0.0), (265.5, -0.5))
TRANSFER_ABSOLUTE_Y_UNCERTAINTY = 0.025
GROUP_RE = re.compile(rb"q\n6 0 0 6 348 1968 cm\n(.*?)\nQ", re.DOTALL)
POINT_RE = re.compile(
    r"^([-+]?[0-9]*\.?[0-9]+) ([-+]?[0-9]*\.?[0-9]+) ([ml])$",
    re.MULTILINE,
)
CURVES = {
    "solid": (
        "fig7_transfer_active_random_beta10",
        "paper_standard_active_random_beta10",
    ),
    "dashdot": (
        "fig7_transfer_active_alfvenic_beta10",
        "paper_standard_active_alfvenic_beta10",
    ),
    "dashed": (
        "fig7_transfer_passive_alfvenic_beta10",
        "paper_standard_passive_alfvenic_beta10",
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
    """Return the vector-content stream containing both Figure 7 panels."""

    for compressed in re.findall(
        rb"stream\r?\n(.*?)\r?\nendstream", pdf, re.DOTALL
    ):
        try:
            content = zlib.decompress(compressed)
        except zlib.error:
            continue
        if b"91 344 m" in content and b"[8 3 2 3] 0 d" in content:
            return content
    raise ValueError("did not locate the Figure 7 vector-content stream")


def map_coordinate(raw_x: float, raw_y: float) -> tuple[float, float]:
    """Map lower-panel vector coordinates onto the physical axes."""

    (x0, value0), (x1, value1) = X_LOG_TICKS
    x = value0 * 10.0 ** (
        (raw_x - x0) * math.log10(value1 / value0) / (x1 - x0)
    )
    (y0, yvalue0), (y1, yvalue1) = Y_LINEAR_TICKS
    y = yvalue0 + (raw_y - y0) * (yvalue1 - yvalue0) / (y1 - y0)
    return x, y


def extract_curves(content: bytes) -> tuple[dict[str, list[tuple[float, float]]], int]:
    """Extract lower-panel transfer paths and omit boundary-clipped vertices."""

    extracted = {style: [] for style in CURVES}
    clipped = 0
    xmin, xmax, ymin, ymax = PLOT
    for raw_group in GROUP_RE.findall(content):
        group = raw_group.decode("ascii", errors="ignore")
        if "0 G\n0 g" not in group:
            continue
        if "[10 6] 0 d" in group:
            style = "dashed"
        elif "[8 3 2 3] 0 d" in group:
            style = "dashdot"
        else:
            style = "solid"
        points = POINT_RE.findall(group)
        if style not in CURVES or len(points) < 3:
            continue
        if not any(float(y) > ymin for _, y, _ in points):
            continue
        for raw_x, raw_y, _ in points:
            x_coord = float(raw_x)
            y_coord = float(raw_y)
            if (
                x_coord <= xmin or x_coord >= xmax
                or y_coord <= ymin or y_coord >= ymax
            ):
                clipped += 1
                continue
            extracted[style].append(map_coordinate(x_coord, y_coord))
    for style, points in extracted.items():
        points.sort()
        if len(points) < 2:
            raise ValueError(f"no usable plotted vector curve found for {style}")
        if any(x1 >= x2 for (x1, _), (x2, _) in zip(points, points[1:])):
            raise ValueError(f"non-increasing coordinates extracted for {style}")
        if any(not math.isfinite(x) or not math.isfinite(y) for x, y in points):
            raise ValueError(f"nonfinite coordinates extracted for {style}")
    return extracted, clipped


def write_curves(source_figure: Path, output_dir: Path) -> Path:
    """Write admitted transfer CSVs, copied figure, metadata, and manifest."""

    digest = sha256_file(source_figure)
    if digest != EXPECTED_FIGURE_SHA256:
        raise ValueError(
            "Figure 7 SHA-256 does not match pinned arXiv:2405.02418v2 source"
        )
    curves, clipped = extract_curves(pdf_content(source_figure.read_bytes()))
    output_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_figure, output_dir / "fig7.pdf")
    entries: list[dict[str, str]] = []
    point_counts: dict[str, int] = {}
    for style, (curve_id, case) in CURVES.items():
        csv_path = output_dir / f"{curve_id}.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("x", "y", "y_uncertainty"))
            for x, y in curves[style]:
                writer.writerow(
                    (
                        f"{x:.17g}",
                        f"{y:.17g}",
                        f"{TRANSFER_ABSOLUTE_Y_UNCERTAINTY:.17g}",
                    )
                )
        point_counts[curve_id] = len(curves[style])
        entries.append({
            "id": curve_id,
            "case": case,
            "product": "pressure_transfer.transfer_normalized_by_total",
            "data_file": csv_path.name,
            "data_sha256": sha256_file(csv_path),
            "interpolation": "linear",
        })
    manifest = {
        "schema_version": 1,
        "provenance": {
            "method": "digitized",
            "source_description": (
                "MKS24 arXiv:2405.02418v2 source/fig7.pdf lower panel; "
                "vector polyline vertices extracted through "
                "tick-calibrated plotted axes"
            ),
            "source_figure": "fig7.pdf",
            "source_figure_sha256": digest,
            "digitization_tool": Path(__file__).name,
            "uncertainty_description": (
                "Absolute T_Delta_p/T_total uncertainty 0.025, "
                "conservatively covering plotted stroke width and "
                "tick-calibrated mapping; vertices clipped at plotting "
                "boundaries are omitted. The dimensionful upper-panel "
                "spectra are excluded pending a paper-to-AthenaK transform."
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
            "lower_transfer_panel": {
                "x_log_ticks": [list(tick) for tick in X_LOG_TICKS],
                "y_linear_ticks": [list(tick) for tick in Y_LINEAR_TICKS],
            },
        },
        "absolute_y_uncertainty": {
            "pressure_transfer.transfer_normalized_by_total": (
                TRANSFER_ABSOLUTE_Y_UNCERTAINTY
            ),
        },
        "point_counts": point_counts,
        "clipped_vertices_omitted": clipped,
        "excluded_curves": {
            "upper_panel_gradient_spectra": (
                "The plotted gradient spectra are dimensionful; MKS24's "
                "reported p0 code-unit scale is not yet reconciled with "
                "AthenaK's v_A=1 normalization."
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
    print(f"Wrote digitized Figure 7 curve manifest to {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
