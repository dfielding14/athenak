#!/usr/bin/env python3
"""Extract the dimensionless MKS24 Figure 11 alignment-peak curves."""

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
    "efab7ae64d5837ee621249acfedf79cd77bb221e37504f8c41a79bf5ba30f5f8"
)
PLOT_X_MIN = 55.0
PLOT_X_MAX = 493.0
X_LOG_TICK = 256.63327026
X_LOG_TICK_VALUE = 100.0
X_LOG_DECADE_WIDTH = (348.93823242 - X_LOG_TICK) / math.log10(2.0)
Y_ZERO = 246.0
Y_POINT_TWO = 210.92308044
ALIGNMENT_Y_UNCERTAINTY = 0.005
GROUP_RE = re.compile(
    rb"q\n6 0 0 6 804 2310\.00024414 cm\n(.*?)\nQ", re.DOTALL
)
POINT_RE = re.compile(
    r"^([-+]?[0-9]*\.?[0-9]+) ([-+]?[0-9]*\.?[0-9]+) ([ml])$",
    re.MULTILINE,
)
COLOR_RE = re.compile(r"^([0-9.]+ [0-9.]+ [0-9.]+) RG$", re.MULTILINE)
CURVES = {
    "0 0.4471 0.7412": (
        "fig11_alignment_active_alfvenic_beta10_nperp96",
        "paper_scale_separation_active_alfvenic_beta10_nperp96",
    ),
    "0.851 0.3255 0.098": (
        "fig11_alignment_active_alfvenic_beta10_nperp192",
        "paper_standard_active_alfvenic_beta10",
    ),
    "0.4941 0.1843 0.5569": (
        "fig11_alignment_active_alfvenic_beta10_nperp384",
        "paper_scale_separation_active_alfvenic_beta10_nperp384",
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
    """Return the decompressed vector-content stream containing both panels."""

    for compressed in re.findall(
        rb"stream\r?\n(.*?)\r?\nendstream", pdf, re.DOTALL
    ):
        try:
            content = zlib.decompress(compressed)
        except zlib.error:
            continue
        if b"55 246 m" in content and b"0.4941 0.1843 0.5569 RG" in content:
            return content
    raise ValueError("did not locate the Figure 11 vector-content stream")


def map_x(raw_x: float) -> float:
    """Map log-axis PDF coordinates to physical perpendicular wavenumber."""

    return X_LOG_TICK_VALUE * 10.0 ** (
        (raw_x - X_LOG_TICK) / X_LOG_DECADE_WIDTH
    )


def map_y(raw_y: float) -> float:
    """Map the lower-panel vertical coordinate to the alignment peak cosine."""

    return (Y_ZERO - raw_y) * 0.2 / (Y_ZERO - Y_POINT_TWO)


def extract_curves(
    content: bytes,
) -> tuple[dict[str, list[tuple[float, float]]], dict[str, object]]:
    """Extract lower-panel vector paths and convert plotted axes to data."""

    extracted = {color: [] for color in CURVES}
    x_clipped_vertices = 0
    alignment_shells: set[int] = set()
    for raw_group in GROUP_RE.findall(content):
        group = raw_group.decode("ascii", errors="ignore")
        color = COLOR_RE.search(group)
        if color is None or color.group(1) not in CURVES:
            continue
        points = POINT_RE.findall(group)
        if len(points) < 3:
            continue
        mean_y = sum(float(point[1]) for point in points) / len(points)
        if mean_y <= 132.0:
            continue
        mapped: list[tuple[float, float]] = []
        for raw_x, raw_y, _ in points:
            x_coordinate = float(raw_x)
            if x_coordinate <= PLOT_X_MIN or x_coordinate >= PLOT_X_MAX:
                x_clipped_vertices += 1
                continue
            x = map_x(x_coordinate)
            shell = round(x / math.pi)
            if abs(x / math.pi - shell) > 1.0e-3:
                raise ValueError("alignment vertex is not on a k_perp shell")
            alignment_shells.add(shell)
            mapped.append((shell * math.pi, map_y(float(raw_y))))
        extracted[color.group(1)] = mapped
    for color, points in extracted.items():
        if len(points) < 2:
            raise ValueError(f"no usable plotted vector curve found for {color}")
        if any(x1 >= x2 for (x1, _), (x2, _) in zip(points, points[1:])):
            raise ValueError(f"non-increasing coordinates extracted for {color}")
        if any(not math.isfinite(x) or not math.isfinite(y) for x, y in points):
            raise ValueError(f"nonfinite coordinates extracted for {color}")
    return extracted, {
        "x_clipped_vertices_omitted": x_clipped_vertices,
        "recommended_alignment_shells": sorted(alignment_shells),
    }


def write_curves(source_figure: Path, output_dir: Path) -> Path:
    """Write admitted alignment CSVs, copied figure, metadata, and manifest."""

    digest = sha256_file(source_figure)
    if digest != EXPECTED_FIGURE_SHA256:
        raise ValueError(
            "Figure 11 SHA-256 does not match pinned arXiv:2405.02418v2 source"
        )
    curves, extraction = extract_curves(pdf_content(source_figure.read_bytes()))
    output_dir.mkdir(parents=True, exist_ok=True)
    archived_figure = output_dir / "fig11.pdf"
    shutil.copy2(source_figure, archived_figure)
    entries: list[dict[str, str]] = []
    point_counts: dict[str, int] = {}
    for color, (curve_id, case) in CURVES.items():
        csv_path = output_dir / f"{curve_id}.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("x", "y", "y_uncertainty"))
            for x, y in curves[color]:
                writer.writerow(
                    (f"{x:.17g}", f"{y:.17g}", f"{ALIGNMENT_Y_UNCERTAINTY:.17g}")
                )
        point_counts[curve_id] = len(curves[color])
        entries.append({
            "id": curve_id,
            "case": case,
            "product": "alignment_peak.cos_theta",
            "data_file": csv_path.name,
            "data_sha256": sha256_file(csv_path),
            "interpolation": "linear",
        })
    manifest = {
        "schema_version": 1,
        "provenance": {
            "method": "digitized",
            "source_description": (
                "MKS24 arXiv:2405.02418v2 source/fig11.pdf lower panel; "
                "vector polyline vertices extracted through tick-calibrated "
                "plotted axes"
            ),
            "source_figure": archived_figure.name,
            "source_figure_sha256": digest,
            "digitization_tool": Path(__file__).name,
            "uncertainty_description": (
                "Absolute |cos(theta)| uncertainty 0.005 for the lower-panel "
                "alignment curves, conservatively covering plotted stroke "
                "width and tick-calibrated coordinate mapping; vertices "
                "clipped at the horizontal plot limits are omitted. The "
                "dimensionful upper-panel kinetic spectra are excluded "
                "pending a paper-to-AthenaK code-unit transform."
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
            "x_log_tick": [X_LOG_TICK, X_LOG_TICK_VALUE],
            "x_log_decade_width": X_LOG_DECADE_WIDTH,
            "lower_y_ticks": [[Y_ZERO, 0.0], [Y_POINT_TWO, 0.2]],
        },
        "alignment_y_uncertainty": ALIGNMENT_Y_UNCERTAINTY,
        "point_counts": point_counts,
        **extraction,
        "excluded_curves": {
            "upper_panel_spectra.velocity": (
                "The plotted kinetic energy spectrum is dimensionful; "
                "MKS24's code-unit scale is not yet reconciled with AthenaK's "
                "v_A=1 normalization."
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
    print(f"Wrote digitized Figure 11 curve manifest to {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
