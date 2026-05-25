#!/usr/bin/env python3
"""Extract normalized eddy-anisotropy curves from MKS24 Figure 5(b)."""

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
    "8886c425749f97af6f50b7960db071179758e6ea68b79ca547dfe1b163c2e0fe"
)
PLOT = (57.0, 346.0, 24.0, 272.0)
X_LOG_TICKS = ((155.61149597, 0.05), (243.22183228, 0.1))
Y_LOG_TICKS = ((215.74073792, 0.1), (165.60038757, 0.15))
LENGTH_RELATIVE_Y_UNCERTAINTY = 0.05
GROUP_RE = re.compile(rb"q\n6 0 0 6 1302 2214 cm\n(.*?)\nQ", re.DOTALL)
POINT_RE = re.compile(
    r"^([-+]?[0-9]*\.?[0-9]+) ([-+]?[0-9]*\.?[0-9]+) ([ml])$",
    re.MULTILINE,
)
COLOR_RE = re.compile(r"^([0-9.]+ [0-9.]+ [0-9.]+) RG$", re.MULTILINE)
CURVES = {
    ("0 0.4471 0.7412", "solid"): (
        "fig5b_eddy_velocity_perp_active_random_beta10",
        "paper_standard_active_random_beta10",
        "eddy_anisotropy.velocity_perp",
    ),
    ("0.851 0.3255 0.098", "solid"): (
        "fig5b_eddy_magnetic_perp_active_random_beta10",
        "paper_standard_active_random_beta10",
        "eddy_anisotropy.magnetic_perp",
    ),
    ("0 0.4471 0.7412", "dash_dotted"): (
        "fig5b_eddy_velocity_perp_active_alfvenic_beta10",
        "paper_standard_active_alfvenic_beta10",
        "eddy_anisotropy.velocity_perp",
    ),
    ("0.851 0.3255 0.098", "dash_dotted"): (
        "fig5b_eddy_magnetic_perp_active_alfvenic_beta10",
        "paper_standard_active_alfvenic_beta10",
        "eddy_anisotropy.magnetic_perp",
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
    """Return the decompressed vector-content stream containing Figure 5(b)."""

    for compressed in re.findall(
        rb"stream\r?\n(.*?)\r?\nendstream", pdf, re.DOTALL
    ):
        try:
            content = zlib.decompress(compressed)
        except zlib.error:
            continue
        if b"57 272 m" in content and b"[8 3 2 3] 0 d" in content:
            return content
    raise ValueError("did not locate the Figure 5(b) vector-content stream")


def map_coordinate(raw_x: float, raw_y: float) -> tuple[float, float]:
    """Map plotted axes to normalized perpendicular and parallel lengths."""

    (x0, value0), (x1, value1) = X_LOG_TICKS
    perpendicular = value0 * 10.0 ** (
        (raw_x - x0) * math.log10(value1 / value0) / (x1 - x0)
    )
    (y0, yvalue0), (y1, yvalue1) = Y_LOG_TICKS
    parallel = yvalue0 * 10.0 ** (
        (raw_y - y0) * math.log10(yvalue1 / yvalue0) / (y1 - y0)
    )
    return perpendicular, parallel


def extract_curves(
    content: bytes,
) -> tuple[dict[tuple[str, str], list[tuple[float, float]]], int]:
    """Extract visible curve paths and omit plot-boundary-clipped vertices."""

    extracted = {key: [] for key in CURVES}
    clipped_vertices = 0
    xmin, xmax, ymin, ymax = PLOT
    for raw_group in GROUP_RE.findall(content):
        group = raw_group.decode("ascii", errors="ignore")
        color = COLOR_RE.search(group)
        if color is None:
            continue
        style = "dash_dotted" if "[8 3 2 3] 0 d" in group else "solid"
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
        if any(
            not math.isfinite(x) or not math.isfinite(y) or x <= 0.0 or y <= 0.0
            for x, y in points
        ):
            raise ValueError(f"invalid normalized length extracted for {key}")
    return extracted, clipped_vertices


def write_curves(source_figure: Path, output_dir: Path) -> Path:
    """Write admitted eddy curves, copied figure, metadata, and manifest."""

    digest = sha256_file(source_figure)
    if digest != EXPECTED_FIGURE_SHA256:
        raise ValueError(
            "Figure 5(b) SHA-256 does not match pinned arXiv:2405.02418v2 source"
        )
    curves, clipped_vertices = extract_curves(pdf_content(source_figure.read_bytes()))
    output_dir.mkdir(parents=True, exist_ok=True)
    archived_figure = output_dir / "fig5b.pdf"
    shutil.copy2(source_figure, archived_figure)
    entries: list[dict[str, str]] = []
    point_counts: dict[str, int] = {}
    for key, (curve_id, case, product) in CURVES.items():
        csv_path = output_dir / f"{curve_id}.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("x", "y", "y_uncertainty"))
            for x, y in curves[key]:
                writer.writerow(
                    (
                        f"{x:.17g}",
                        f"{y:.17g}",
                        f"{LENGTH_RELATIVE_Y_UNCERTAINTY * y:.17g}",
                    )
                )
        point_counts[curve_id] = len(curves[key])
        entries.append({
            "id": curve_id,
            "case": case,
            "product": product,
            "data_file": csv_path.name,
            "data_sha256": sha256_file(csv_path),
            "interpolation": "loglog",
        })
    manifest = {
        "schema_version": 1,
        "provenance": {
            "method": "digitized",
            "source_description": (
                "MKS24 arXiv:2405.02418v2 source/fig5b.pdf; vector "
                "polyline vertices extracted through tick-calibrated log "
                "axes for normalized eddy lengths"
            ),
            "source_figure": archived_figure.name,
            "source_figure_sha256": digest,
            "digitization_tool": Path(__file__).name,
            "uncertainty_description": (
                "Five-percent relative normalized-parallel-length uncertainty "
                "conservatively covers plotted stroke width, log-axis "
                "calibration, and vector extraction; vertices clipped at "
                "the plotted boundaries are omitted. Blue and orange legend "
                "curves map to u_perp and B_perp; solid and dash-dotted "
                "curves map to random and Alfvenic forcing."
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
            "plot": list(PLOT),
            "x_log_ticks_ell_perp_over_l0": [
                list(tick) for tick in X_LOG_TICKS
            ],
            "y_log_ticks_ell_parallel_over_l0": [
                list(tick) for tick in Y_LOG_TICKS
            ],
            "analysis_normalization": (
                "l0 identified with L_perp for the shared [1,1,2] domain; "
                "the analyzer emits lengths divided by sqrt(Lx Ly)"
            ),
        },
        "relative_y_uncertainty": LENGTH_RELATIVE_Y_UNCERTAINTY,
        "point_counts": point_counts,
        "clipped_vertices_omitted": clipped_vertices,
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
    print(f"Wrote digitized Figure 5(b) curve manifest to {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
