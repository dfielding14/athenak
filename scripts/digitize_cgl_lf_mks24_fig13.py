#!/usr/bin/env python3
"""Extract dimensionless, directly comparable MKS24 Figure 13 curves."""

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


EXPECTED_FIGURE_SHA256 = {
    "fig13a.pdf": (
        "51cb2fef25ce29c649cb7f172471917276687f51e5dae2be26d54e6f7f7e1c92"
    ),
    "fig13b.pdf": (
        "92be89eb69ed03efd754c53c44da66a9a5afa34c561428bbcedd1cc287fc3013"
    ),
    "fig13c.pdf": (
        "1e67f92b91381ab634add89a4a2e9cdaddd181338b7aef3f6669fad7c2881afc"
    ),
    "fig13d.pdf": (
        "bf1c063fa5a54c2dbbe9dce417490b5bf821fcd2e1cb1926d78c923b90387a24"
    ),
}
PANEL_CONFIG = {
    "fig13b.pdf": {
        "plot": (60.0, 279.0, 19.0, 215.0),
        "group": re.compile(
            rb"q\n6 0 0 6 1524 2387\.99975586 cm\n(.*?)\nQ", re.DOTALL
        ),
        "marker": b"60 215 m",
    },
    "fig13d.pdf": {
        "plot": (82.0, 571.0, 20.0, 169.0),
        "group": re.compile(
            rb"q\n6 0 0 6 555 2493 cm\n(.*?)\nQ", re.DOTALL
        ),
        "marker": b"82 169 m",
    },
}
POINT_RE = re.compile(
    r"^([-+]?[0-9]*\.?[0-9]+) ([-+]?[0-9]*\.?[0-9]+) ([ml])$",
    re.MULTILINE,
)
COLOR_RE = re.compile(r"^([0-9.]+ [0-9.]+ [0-9.]+) RG$", re.MULTILINE)
PDF_RELATIVE_Y_UNCERTAINTY = 0.05
TRANSFER_ABSOLUTE_Y_UNCERTAINTY = 0.025
CURVES = {
    ("fig13b.pdf", "0 0.4471 0.7412", "solid"): (
        "fig13b_beta_delta_nulim_20",
        "paper_nulim_beta100_20",
        "pdf.beta_delta",
        "linear",
    ),
    ("fig13b.pdf", "0.851 0.3255 0.098", "solid"): (
        "fig13b_beta_delta_nulim_200",
        "paper_nulim_beta100_200",
        "pdf.beta_delta",
        "linear",
    ),
    ("fig13b.pdf", "0.4941 0.1843 0.5569", "solid"): (
        "fig13b_beta_delta_nulim_hardwall",
        "paper_nulim_beta100_hardwall",
        "pdf.beta_delta",
        "linear",
    ),
    ("fig13b.pdf", "0.4941 0.1843 0.5569", "dashed"): (
        "fig13b_beta_delta_passive_hardwall",
        "paper_standard_passive_alfvenic_beta100",
        "pdf.beta_delta",
        "linear",
    ),
    ("fig13d.pdf", "0 0.4471 0.7412", "solid"): (
        "fig13d_transfer_nulim_20",
        "paper_nulim_beta100_20",
        "pressure_transfer.transfer_normalized_by_total",
        "linear",
    ),
    ("fig13d.pdf", "0.851 0.3255 0.098", "solid"): (
        "fig13d_transfer_nulim_200",
        "paper_nulim_beta100_200",
        "pressure_transfer.transfer_normalized_by_total",
        "linear",
    ),
    ("fig13d.pdf", "0.4941 0.1843 0.5569", "solid"): (
        "fig13d_transfer_nulim_hardwall",
        "paper_nulim_beta100_hardwall",
        "pressure_transfer.transfer_normalized_by_total",
        "linear",
    ),
    ("fig13d.pdf", "0.4941 0.1843 0.5569", "dashed"): (
        "fig13d_transfer_passive_hardwall",
        "paper_standard_passive_alfvenic_beta100",
        "pressure_transfer.transfer_normalized_by_total",
        "linear",
    ),
}


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pdf_content(pdf: bytes, marker: bytes) -> bytes:
    """Return the compressed vector-content stream containing a panel axis."""

    stream_header = re.compile(
        rb"/Filter /FlateDecode /Length [0-9]+ >>\s*stream\r?\n"
    )
    for match in stream_header.finditer(pdf):
        end = pdf.find(b"endstream", match.end())
        if end < 0:
            continue
        try:
            content = zlib.decompress(pdf[match.end():end].rstrip(b"\r\n"))
        except zlib.error:
            continue
        if marker in content and b"0 0.4471 0.7412 RG" in content:
            return content
    raise ValueError("did not locate the Figure 13 vector-content stream")


def map_coordinate(panel: str, raw_x: float, raw_y: float) -> tuple[float, float]:
    """Map plotted vector coordinates onto the panel's physical axes."""

    if panel == "fig13b.pdf":
        x = -3.0 + (raw_x - 60.0) / (279.0 - 60.0) * 5.0
        y = 10.0 ** (
            -(raw_y - 36.87377167) * 2.0 / (155.62458801 - 36.87377167)
        )
        return x, y
    if panel == "fig13d.pdf":
        x = 10.0 * 10.0 ** (
            (raw_x - 138.91017151) / (506.30499268 - 138.91017151)
        )
        y = 0.5 - (raw_y - 20.0) * 0.5 / (53.11111069 - 20.0)
        return x, y
    raise ValueError(f"no admitted coordinate mapping for {panel}")


def extract_curves(
    contents: dict[str, bytes],
) -> tuple[dict[tuple[str, str, str], list[tuple[float, float]]], dict[str, int]]:
    """Extract supported visible curves and convert axis coordinates to data."""

    extracted = {key: [] for key in CURVES}
    clipped = {panel: 0 for panel in PANEL_CONFIG}
    for panel, config in PANEL_CONFIG.items():
        xmin, xmax, ymin, ymax = config["plot"]
        for raw_group in config["group"].findall(contents[panel]):
            group = raw_group.decode("ascii", errors="ignore")
            color = COLOR_RE.search(group)
            if color is None:
                continue
            style = "dashed" if "[10 6] 0 d" in group else "solid"
            key = (panel, color.group(1), style)
            if key not in CURVES:
                continue
            points = POINT_RE.findall(group)
            if len(points) < 3:
                continue
            for raw_x, raw_y, _ in points:
                x_coord = float(raw_x)
                y_coord = float(raw_y)
                if (
                    x_coord <= xmin or x_coord >= xmax
                    or y_coord <= ymin or y_coord >= ymax
                ):
                    clipped[panel] += 1
                    continue
                extracted[key].append(map_coordinate(panel, x_coord, y_coord))
    for key, points in extracted.items():
        points.sort()
        if len(points) < 2:
            raise ValueError(f"no usable plotted vector curve found for {key}")
        if any(x1 >= x2 for (x1, _), (x2, _) in zip(points, points[1:])):
            raise ValueError(f"non-increasing coordinates extracted for {key}")
        if any(not math.isfinite(x) or not math.isfinite(y) for x, y in points):
            raise ValueError(f"nonfinite coordinates extracted for {key}")
    return extracted, clipped


def write_curves(source_dir: Path, output_dir: Path) -> Path:
    """Write supported panel CSVs, copied figures, metadata, and a manifest."""

    digests: dict[str, str] = {}
    for name, expected in EXPECTED_FIGURE_SHA256.items():
        source = source_dir / name
        digest = sha256_file(source)
        if digest != expected:
            raise ValueError(
                f"{name} SHA-256 does not match pinned arXiv:2405.02418v2 source"
            )
        digests[name] = digest
    contents = {
        panel: pdf_content((source_dir / panel).read_bytes(), config["marker"])
        for panel, config in PANEL_CONFIG.items()
    }
    curves, clipped = extract_curves(contents)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in EXPECTED_FIGURE_SHA256:
        shutil.copy2(source_dir / name, output_dir / name)
    entries: list[dict[str, str]] = []
    point_counts: dict[str, int] = {}
    for key, (curve_id, case, product, interpolation) in CURVES.items():
        csv_path = output_dir / f"{curve_id}.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("x", "y", "y_uncertainty"))
            for x, y in curves[key]:
                uncertainty = (
                    PDF_RELATIVE_Y_UNCERTAINTY * y
                    if key[0] == "fig13b.pdf"
                    else TRANSFER_ABSOLUTE_Y_UNCERTAINTY
                )
                writer.writerow(
                    (f"{x:.17g}", f"{y:.17g}", f"{uncertainty:.17g}")
                )
        point_counts[curve_id] = len(curves[key])
        entries.append(
            {
                "id": curve_id,
                "case": case,
                "product": product,
                "data_file": csv_path.name,
                "data_sha256": sha256_file(csv_path),
                "interpolation": interpolation,
            }
        )
    manifest = {
        "schema_version": 1,
        "provenance": {
            "method": "digitized",
            "source_description": (
                "MKS24 arXiv:2405.02418v2 source/fig13b.pdf and "
                "source/fig13d.pdf; "
                "vector polyline vertices extracted "
                "through tick-calibrated plotted axes"
            ),
            "source_figures": [
                {
                    "source_figure": name,
                    "source_figure_sha256": digests[name],
                }
                for name in ("fig13b.pdf", "fig13d.pdf")
            ],
            "digitization_tool": Path(__file__).name,
            "uncertainty_description": (
                "Five-percent relative y uncertainty for dimensionless "
                "beta-Delta PDF curves and absolute 0.025 uncertainty for "
                "dimensionless T_Delta_p/T_total transfer curves, "
                "conservatively covering plotted stroke width and "
                "tick-calibrated mapping; vertices clipped at plotting "
                "boundaries are omitted. Dimensionful panels are excluded "
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
        "source_figure_sha256": digests,
        "axis_mapping": {
            "fig13b.pdf": {
                "x_linear_ticks": [[60.0, -3.0], [279.0, 2.0]],
                "y_log_ticks": [[36.87377167, 1.0], [155.62458801, 1.0e-2]],
            },
            "fig13d.pdf": {
                "x_log_ticks": [[138.91017151, 10.0], [506.30499268, 100.0]],
                "y_linear_ticks": [[20.0, 0.5], [53.11111069, 0.0]],
            },
        },
        "relative_y_uncertainty": {
            "pdf.beta_delta": PDF_RELATIVE_Y_UNCERTAINTY,
        },
        "absolute_y_uncertainty": {
            "pressure_transfer.transfer_normalized_by_total": (
                TRANSFER_ABSOLUTE_Y_UNCERTAINTY
            ),
        },
        "point_counts": point_counts,
        "clipped_vertices_omitted": clipped,
        "excluded_panels": {
            "fig13a.pdf": (
                "The plotted kinetic spectrum is dimensionful; MKS24's "
                "reported p0 code-unit scale is not yet reconciled with "
                "AthenaK's v_A=1 normalization."
            ),
            "fig13c.pdf": (
                "The plotted strain PDF uses "
                "bhat bhat : grad(u) / <B^2>; its scale requires the "
                "unresolved paper-to-AthenaK code-unit transform."
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
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    manifest = write_curves(args.source_dir, args.output_dir)
    print(f"Wrote digitized Figure 13 curve manifest to {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
