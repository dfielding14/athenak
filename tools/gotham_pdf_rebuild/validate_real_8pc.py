#!/usr/bin/env python3.11
"""Validate a full res_8pc/phase2 PDF rebuild against archived sequence 00062."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable

import numpy as np


ARCHIVED_SEQUENCE = "00062"
EXPECTED_SOURCE_SEQUENCE = "00028"
EXPECTED_SOURCE_CYCLE: int | None = 202829
EXPECTED_SHARDS = 1024
EXPECTED_DOMAIN_VOLUME = 400.0**3
EXPECTED_OUTPUT_NUMBER = "00062"
EXPECTED_OUTPUT_TIME = 10.520002963096713

GLOBAL_REL_TOL = 2.0e-9
CONTROL_TOTAL_TOL = 5.0e-6
HEADER_EDGE_RTOL = 5.0e-13
HEADER_EDGE_ATOL = 5.0e-15

ARCHIVED_STREAMS = {
    "output29": "pdf_coord_abscostheta_coord_r",
    "output31": "pdf_coord_costheta_coord_r_temperature_vel_sph_r",
}

EXPECTED_CONTROLS = {
    "output29": {
        "variables": ("coord_abscostheta", "coord_r"),
        "shape": (18, 130),
        "weight": "edot_sph",
    },
    "output31": {
        "variables": ("coord_costheta", "coord_r", "temperature", "vel_sph_r"),
        "shape": (18, 34, 130, 130),
        "weight": "edot_sph",
    },
}

COMPARISON_LIMITS = {
    "output29_radial_marginal": {"e_l1": 2.0e-5, "e_total": CONTROL_TOTAL_TOL},
    "output31_full_4d": {"e_l1": 2.0e-3, "e_total": CONTROL_TOTAL_TOL},
    "output31_geometry_costheta_radius": {
        "e_l1": 5.0e-5,
        "e_total": CONTROL_TOTAL_TOL,
    },
    "output31_temperature_1d": {"e_l1": 5.0e-4, "e_total": CONTROL_TOTAL_TOL},
    "output31_radial_velocity_1d": {
        "e_l1": 5.0e-4,
        "e_total": CONTROL_TOTAL_TOL,
    },
    "output31_coarsened_temperature_velocity": {
        "e_l1": 5.0e-4,
        "e_total": CONTROL_TOTAL_TOL,
    },
}


class ValidationError(RuntimeError):
    """Raised when the requested validation cannot be completed reliably."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate a complete GOTHAM rebuild against internal invariants and, "
            "unless disabled, archived sparse output29/output31 controls."
        )
    )
    parser.add_argument(
        "rebuilt_output_dir",
        type=Path,
        help="Reducer output directory containing rebuild_manifest.json and output29/output31",
    )
    parser.add_argument(
        "archived_phase2_root",
        type=Path,
        help="Archived res_8pc/phase2 root containing pdf/",
    )
    parser.add_argument("--archived-sequence", default=ARCHIVED_SEQUENCE)
    parser.add_argument("--expected-source-sequence", default=EXPECTED_SOURCE_SEQUENCE)
    parser.add_argument("--expected-source-cycle", type=int, default=EXPECTED_SOURCE_CYCLE)
    parser.add_argument("--skip-source-cycle-check", action="store_true")
    parser.add_argument("--expected-shards", type=int, default=EXPECTED_SHARDS)
    parser.add_argument(
        "--expected-domain-volume", type=float, default=EXPECTED_DOMAIN_VOLUME
    )
    parser.add_argument("--expected-output-number", default=EXPECTED_OUTPUT_NUMBER)
    parser.add_argument(
        "--expected-output-time", type=float, default=EXPECTED_OUTPUT_TIME
    )
    parser.add_argument(
        "--no-controls",
        action="store_true",
        help="Validate only manifest, payload, geometry, and global product invariants",
    )
    parser.add_argument(
        "--report-path",
        type=Path,
        help="Atomically publish the complete JSON report to this path",
    )
    return parser.parse_args()


def configure_expectations(args: argparse.Namespace) -> None:
    global ARCHIVED_SEQUENCE
    global EXPECTED_SOURCE_SEQUENCE
    global EXPECTED_SOURCE_CYCLE
    global EXPECTED_SHARDS
    global EXPECTED_DOMAIN_VOLUME
    global EXPECTED_OUTPUT_NUMBER
    global EXPECTED_OUTPUT_TIME

    ARCHIVED_SEQUENCE = args.archived_sequence
    EXPECTED_SOURCE_SEQUENCE = args.expected_source_sequence
    EXPECTED_SOURCE_CYCLE = (
        None if args.skip_source_cycle_check else args.expected_source_cycle
    )
    EXPECTED_SHARDS = args.expected_shards
    EXPECTED_DOMAIN_VOLUME = args.expected_domain_volume
    EXPECTED_OUTPUT_NUMBER = args.expected_output_number
    EXPECTED_OUTPUT_TIME = args.expected_output_time


def new_report(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "status": "running",
        "archived_sequence": ARCHIVED_SEQUENCE,
        "expected_source_sequence": EXPECTED_SOURCE_SEQUENCE,
        "inputs": {
            "rebuilt_output_dir": str(args.rebuilt_output_dir.resolve()),
            "archived_phase2_root": str(args.archived_phase2_root.resolve()),
        },
        "controls_required": not args.no_controls,
        "thresholds": {
            "comparisons": COMPARISON_LIMITS,
            "global_relative": GLOBAL_REL_TOL,
            "header_edge_rtol": HEADER_EDGE_RTOL,
            "header_edge_atol": HEADER_EDGE_ATOL,
            "domain_volume": EXPECTED_DOMAIN_VOLUME,
            "expected_shards": EXPECTED_SHARDS,
            "expected_output_number": EXPECTED_OUTPUT_NUMBER,
            "expected_output_time": EXPECTED_OUTPUT_TIME,
        },
        "checks": [],
        "failures": [],
    }


def add_check(
    report: dict[str, Any],
    name: str,
    passed: bool,
    *,
    metrics: dict[str, Any] | None = None,
    limits: dict[str, Any] | None = None,
    details: dict[str, Any] | None = None,
) -> None:
    check: dict[str, Any] = {"name": name, "status": "pass" if passed else "fail"}
    if metrics:
        check["metrics"] = metrics
    if limits:
        check["limits"] = limits
    if details:
        check["details"] = details
    report["checks"].append(check)
    if not passed:
        report["failures"].append(name)


def add_skip(report: dict[str, Any], name: str, reason: str) -> None:
    report["checks"].append({"name": name, "status": "skipped", "reason": reason})


def finalize_report(report: dict[str, Any]) -> None:
    statuses = [check["status"] for check in report["checks"]]
    report["summary"] = {
        "passed": statuses.count("pass"),
        "failed": statuses.count("fail"),
        "skipped": statuses.count("skipped"),
    }
    report["status"] = "pass" if not report["failures"] and "fatal_error" not in report else "fail"


def load_pdf_reader() -> tuple[Callable[..., dict[str, Any]], Callable[..., dict[str, Any]], Path]:
    analysis_root = Path(
        os.environ.get(
            "GOTHAM_ANALYSIS_ROOT",
            "/lustre/orion/ast207/proj-shared/gotham/analysis",
        )
    ).resolve()
    reader_path = analysis_root / "gotham_analysis" / "readers" / "pdf.py"
    if not reader_path.is_file():
        raise ValidationError(
            f"production PDF reader not found beneath GOTHAM_ANALYSIS_ROOT: {reader_path}"
        )
    sys.path.insert(0, str(analysis_root))
    try:
        from gotham_analysis.readers.pdf import read_pdf, read_pdf_header
    except Exception as exc:
        raise ValidationError(f"unable to import production PDF reader: {exc}") from exc
    return read_pdf, read_pdf_header, analysis_root


def canonical_weight(header: dict[str, Any]) -> str:
    weight = header.get("weight")
    if weight == "variable":
        weight = header.get("weight_variable")
    if not isinstance(weight, str) or not weight:
        raise ValidationError("PDF header does not identify its weight")
    return weight


def stream_dense_payload(path: Path, total_bins: int) -> dict[str, float]:
    expected_bytes = (total_bins + 1) * np.dtype(np.float64).itemsize
    if not path.is_file():
        raise ValidationError(f"missing rebuilt payload: {path}")
    if path.stat().st_size != expected_bytes:
        raise ValidationError(
            f"unexpected dense payload size for {path}: "
            f"expected {expected_bytes}, found {path.stat().st_size}"
        )

    raw = np.memmap(path, dtype=np.float64, mode="r")
    embedded_time = float(raw[0])
    if not math.isfinite(embedded_time):
        raise ValidationError(f"non-finite embedded time in {path}")

    total = 0.0
    absolute_total = 0.0
    chunk_values = 1_048_576
    for first in range(1, raw.size, chunk_values):
        chunk = np.asarray(raw[first : first + chunk_values])
        if not np.all(np.isfinite(chunk)):
            raise ValidationError(f"non-finite PDF value in {path}")
        total += float(np.sum(chunk, dtype=np.float64))
        absolute_total += float(np.sum(np.abs(chunk), dtype=np.float64))
    del raw
    if not math.isfinite(total) or not math.isfinite(absolute_total):
        raise ValidationError(f"non-finite aggregate PDF weight in {path}")
    return {"time": embedded_time, "total": total, "abs_sum": absolute_total}


def read_manifest_and_products(
    rebuilt_root: Path,
    read_pdf_header: Callable[..., dict[str, Any]],
    report: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    manifest_path = rebuilt_root / "rebuild_manifest.json"
    if not manifest_path.is_file():
        raise ValidationError(f"missing rebuild manifest: {manifest_path}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(f"unable to read rebuild manifest {manifest_path}: {exc}") from exc

    output_number = manifest.get("output_number")
    source_time = manifest.get("source_time")
    output_time = manifest.get("output_time", source_time)
    entries = manifest.get("products")
    if not isinstance(output_number, str) or not output_number:
        raise ValidationError("rebuild manifest has no valid output_number")
    if not isinstance(source_time, (int, float)) or not math.isfinite(source_time):
        raise ValidationError("rebuild manifest has no finite source_time")
    if not isinstance(output_time, (int, float)) or not math.isfinite(output_time):
        raise ValidationError("rebuild manifest has no finite output_time")
    if not isinstance(entries, list) or not entries:
        raise ValidationError("rebuild manifest has no products")

    products: dict[str, dict[str, Any]] = {}
    worst_manifest_error = 0.0
    worst_manifest_product = ""
    max_payload_time_error = 0.0
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValidationError("rebuild manifest contains a non-object product entry")
        product_id = entry.get("id")
        if not isinstance(product_id, str) or not product_id:
            raise ValidationError("rebuild manifest contains a product without an id")
        if product_id in products:
            raise ValidationError(f"duplicate product id in rebuild manifest: {product_id}")

        product_dir = rebuilt_root / product_id
        header_path = product_dir / "gotham.header.pdf"
        payload_path = product_dir / f"gotham.{output_number}.pdf"
        header = read_pdf_header(header_path)
        if header.get("format") != "dense":
            raise ValidationError(f"rebuilt product is not dense: {header_path}")
        if header.get("distribution") != "global_rebuild":
            raise ValidationError(f"rebuilt product is not a global rebuild: {header_path}")
        if entry.get("total_bins") != header["total_bins"]:
            raise ValidationError(f"manifest/header total_bins mismatch for {product_id}")

        weight = canonical_weight(header)
        if entry.get("weight") != weight:
            raise ValidationError(f"manifest/header weight mismatch for {product_id}")
        payload = stream_dense_payload(payload_path, int(header["total_bins"]))
        manifest_sum = entry.get("sum")
        if not isinstance(manifest_sum, (int, float)) or not math.isfinite(manifest_sum):
            raise ValidationError(f"manifest has a non-finite sum for {product_id}")
        scale = max(payload["abs_sum"], np.finfo(np.float64).tiny)
        manifest_error = abs(payload["total"] - float(manifest_sum)) / scale
        if manifest_error > worst_manifest_error:
            worst_manifest_error = manifest_error
            worst_manifest_product = product_id
        max_payload_time_error = max(
            max_payload_time_error, abs(payload["time"] - float(output_time))
        )

        products[product_id] = {
            "id": product_id,
            "weight": weight,
            "header": header,
            "header_path": header_path,
            "payload_path": payload_path,
            "time": payload["time"],
            "total": payload["total"],
            "abs_sum": payload["abs_sum"],
            "manifest_sum": float(manifest_sum),
            "manifest_sum_relative_error": manifest_error,
        }

    add_check(
        report,
        "manifest_product_sums_match_payloads",
        worst_manifest_error <= GLOBAL_REL_TOL,
        metrics={
            "worst_relative_error": worst_manifest_error,
            "worst_product": worst_manifest_product,
        },
        limits={"relative_error_max": GLOBAL_REL_TOL},
    )
    add_check(
        report,
        "rebuilt_payload_times_match_manifest",
        max_payload_time_error == 0.0,
        metrics={"maximum_absolute_time_error": max_payload_time_error},
        limits={"absolute_time_error_max": 0.0},
    )

    report["manifest"] = {
        key: manifest.get(key)
        for key in (
            "source_input_dir",
            "source_sequence",
            "output_number",
            "source_time",
            "output_time",
            "source_cycle",
            "gamma",
            "shards_available",
            "shards_processed",
            "meshblocks_processed",
            "cells_processed",
            "payload_bytes_read",
            "product_set",
            "geometry_to_logical_key_validated",
            "domain_min",
            "domain_max",
            "root_meshblocks",
            "expected_domain_volume",
            "summed_leaf_volume",
            "elapsed_seconds",
            "histogram_bytes_per_rank",
            "raw_buffer_bytes_per_rank",
            "planned_peak_buffer_bytes_per_rank",
        )
    }
    report["rebuilt_products"] = {
        product_id: {
            "weight": product["weight"],
            "shape": list(product["header"]["shape"]),
            "variables": list(product["header"]["variables"]),
            "time": product["time"],
            "total": product["total"],
            "abs_sum": product["abs_sum"],
            "manifest_sum": product["manifest_sum"],
            "manifest_sum_relative_error": product["manifest_sum_relative_error"],
        }
        for product_id, product in products.items()
    }
    return manifest, products


def add_manifest_checks(
    report: dict[str, Any], manifest: dict[str, Any], archived_root: Path
) -> None:
    add_check(
        report,
        "source_sequence_matches_required_identity",
        manifest.get("source_sequence") == EXPECTED_SOURCE_SEQUENCE,
        metrics={"source_sequence": manifest.get("source_sequence")},
        limits={"required": EXPECTED_SOURCE_SEQUENCE},
    )
    if EXPECTED_SOURCE_CYCLE is None:
        add_skip(report, "source_cycle_is_matching_cycle", "no expected cycle declared")
    else:
        add_check(
            report,
            "source_cycle_is_matching_cycle",
            manifest.get("source_cycle") == EXPECTED_SOURCE_CYCLE,
            metrics={"source_cycle": manifest.get("source_cycle")},
            limits={"required": EXPECTED_SOURCE_CYCLE},
        )
    for key in ("shards_available", "shards_processed"):
        add_check(
            report,
            f"manifest_{key}_matches_required_full_snapshot",
            manifest.get(key) == EXPECTED_SHARDS,
            metrics={key: manifest.get(key)},
            limits={"required": EXPECTED_SHARDS},
        )

    source_input_dir = manifest.get("source_input_dir")
    expected_input_dir = (archived_root / "bin").resolve()
    source_matches = isinstance(source_input_dir, str) and Path(source_input_dir).resolve() == expected_input_dir
    add_check(
        report,
        "manifest_source_input_matches_archived_phase2_bin",
        source_matches,
        metrics={"source_input_dir": source_input_dir},
        limits={"required": str(expected_input_dir)},
    )
    add_check(
        report,
        "manifest_output_number_matches_required_identity",
        manifest.get("output_number") == EXPECTED_OUTPUT_NUMBER,
        metrics={"output_number": manifest.get("output_number")},
        limits={"required": EXPECTED_OUTPUT_NUMBER},
    )
    add_check(
        report,
        "manifest_output_time_matches_required_identity",
        manifest.get("output_time") == EXPECTED_OUTPUT_TIME,
        metrics={"output_time": manifest.get("output_time")},
        limits={"required": EXPECTED_OUTPUT_TIME},
    )
    expected_volume = manifest.get("expected_domain_volume")
    summed_volume = manifest.get("summed_leaf_volume")
    volume_values_valid = (
        isinstance(expected_volume, (int, float))
        and isinstance(summed_volume, (int, float))
        and math.isfinite(expected_volume)
        and math.isfinite(summed_volume)
        and expected_volume > 0.0
    )
    relative_volume_error = (
        abs(float(summed_volume) - float(expected_volume)) / float(expected_volume)
        if volume_values_valid
        else math.inf
    )
    add_check(
        report,
        "manifest_geometry_to_logical_key_validation_passed",
        manifest.get("geometry_to_logical_key_validated") is True,
        metrics={
            "geometry_to_logical_key_validated": manifest.get(
                "geometry_to_logical_key_validated"
            )
        },
        limits={"required": True},
    )
    add_check(
        report,
        "manifest_leaf_volume_closure",
        volume_values_valid
        and float(expected_volume) == EXPECTED_DOMAIN_VOLUME
        and relative_volume_error <= GLOBAL_REL_TOL,
        metrics={
            "expected_domain_volume": expected_volume,
            "summed_leaf_volume": summed_volume,
            "relative_error": relative_volume_error,
        },
        limits={
            "required_domain_volume": EXPECTED_DOMAIN_VOLUME,
            "relative_error_max": GLOBAL_REL_TOL,
        },
    )


def find_archived_payload(archived_root: Path, stream: str) -> Path:
    stream_dir = archived_root / "pdf" / stream
    if not stream_dir.is_dir():
        raise ValidationError(f"missing archived PDF stream directory: {stream_dir}")
    preferred = stream_dir / "node_00000000" / f"gotham.{ARCHIVED_SEQUENCE}.pdf"
    if preferred.is_file() and preferred.with_name("gotham.header.pdf").is_file():
        return preferred
    for candidate in sorted(stream_dir.glob(f"node_*/gotham.{ARCHIVED_SEQUENCE}.pdf")):
        if candidate.with_name("gotham.header.pdf").is_file():
            return candidate
    raise ValidationError(
        f"no archived sequence {ARCHIVED_SEQUENCE} payload with header beneath {stream_dir}"
    )


def validate_control_header(product_id: str, header: dict[str, Any], label: str) -> None:
    expected = EXPECTED_CONTROLS[product_id]
    if tuple(header["variables"]) != expected["variables"]:
        raise ValidationError(
            f"{label} {product_id} variables are {header['variables']}, "
            f"expected {expected['variables']}"
        )
    if tuple(header["shape"]) != expected["shape"]:
        raise ValidationError(
            f"{label} {product_id} shape is {header['shape']}, expected {expected['shape']}"
        )
    if canonical_weight(header) != expected["weight"]:
        raise ValidationError(
            f"{label} {product_id} weight is {canonical_weight(header)}, "
            f"expected {expected['weight']}"
        )


def add_header_alignment_check(
    report: dict[str, Any],
    product_id: str,
    rebuilt_header: dict[str, Any],
    archived_header: dict[str, Any],
) -> None:
    if rebuilt_header["variables"] != archived_header["variables"]:
        raise ValidationError(f"rebuilt/archive variable ordering differs for {product_id}")
    if rebuilt_header["shape"] != archived_header["shape"]:
        raise ValidationError(f"rebuilt/archive shape differs for {product_id}")

    maximum_absolute = 0.0
    maximum_relative = 0.0
    passed = True
    for rebuilt_dim, archived_dim in zip(
        rebuilt_header["dimensions"], archived_header["dimensions"]
    ):
        if rebuilt_dim["scale"] != archived_dim["scale"]:
            raise ValidationError(f"rebuilt/archive scale differs for {product_id}")
        rebuilt_edges = np.asarray(rebuilt_dim["bin_edges"], dtype=np.float64)
        archived_edges = np.asarray(archived_dim["bin_edges"], dtype=np.float64)
        differences = np.abs(rebuilt_edges - archived_edges)
        maximum_absolute = max(maximum_absolute, float(np.max(differences)))
        nonzero = np.abs(archived_edges) > 0.0
        if np.any(nonzero):
            maximum_relative = max(
                maximum_relative,
                float(np.max(differences[nonzero] / np.abs(archived_edges[nonzero]))),
            )
        passed = passed and bool(
            np.allclose(
                rebuilt_edges,
                archived_edges,
                rtol=HEADER_EDGE_RTOL,
                atol=HEADER_EDGE_ATOL,
            )
        )
    add_check(
        report,
        f"{product_id}_rebuilt_archived_headers_align",
        passed,
        metrics={
            "maximum_edge_absolute_error": maximum_absolute,
            "maximum_edge_relative_error": maximum_relative,
        },
        limits={"rtol": HEADER_EDGE_RTOL, "atol": HEADER_EDGE_ATOL},
    )


def load_pdf_array(
    path: Path,
    read_pdf: Callable[..., dict[str, Any]],
    *,
    label: str,
) -> dict[str, Any]:
    try:
        loaded = read_pdf(path)
    except Exception as exc:
        raise ValidationError(f"unable to read {label} PDF {path}: {exc}") from exc
    values = np.asarray(loaded["pdf"], dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValidationError(f"{label} PDF contains non-finite values: {path}")
    loaded["pdf"] = values
    return loaded


def comparison_metrics(candidate: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    if candidate.shape != reference.shape:
        raise ValidationError(
            f"comparison shape mismatch: candidate {candidate.shape}, reference {reference.shape}"
        )
    if not np.all(np.isfinite(candidate)) or not np.all(np.isfinite(reference)):
        raise ValidationError("comparison input contains non-finite values")
    reference_abs_sum = float(np.sum(np.abs(reference), dtype=np.float64))
    if not math.isfinite(reference_abs_sum) or reference_abs_sum <= 0.0:
        raise ValidationError("comparison reference has no finite nonzero absolute weight")
    candidate_total = float(np.sum(candidate, dtype=np.float64))
    reference_total = float(np.sum(reference, dtype=np.float64))
    metrics = {
        "e_l1": float(np.sum(np.abs(candidate - reference), dtype=np.float64))
        / reference_abs_sum,
        "e_total": abs(candidate_total - reference_total) / reference_abs_sum,
        "candidate_total": candidate_total,
        "reference_total": reference_total,
        "reference_abs_sum": reference_abs_sum,
    }
    if not all(math.isfinite(value) for value in metrics.values()):
        raise ValidationError("comparison produced a non-finite aggregate metric")
    return metrics


def add_comparison(
    report: dict[str, Any],
    name: str,
    candidate: np.ndarray,
    reference: np.ndarray,
) -> None:
    limits = COMPARISON_LIMITS[name]
    metrics = comparison_metrics(candidate, reference)
    add_check(
        report,
        name,
        metrics["e_l1"] <= limits["e_l1"] and metrics["e_total"] <= limits["e_total"],
        metrics=metrics,
        limits=limits,
    )


def coarsen_interior_pairs(values: np.ndarray, axis: int) -> np.ndarray:
    moved = np.moveaxis(values, axis, -1)
    interior = moved.shape[-1] - 2
    if interior <= 0 or interior % 2:
        raise ValidationError(
            f"axis {axis} cannot be adjacent-pair coarsened with shape {values.shape}"
        )
    result = np.empty((*moved.shape[:-1], interior // 2 + 2), dtype=np.float64)
    result[..., 0] = moved[..., 0]
    result[..., -1] = moved[..., -1]
    result[..., 1:-1] = np.sum(
        moved[..., 1:-1].reshape(*moved.shape[:-1], interior // 2, 2),
        axis=-1,
        dtype=np.float64,
    )
    return np.moveaxis(result, -1, axis)


def six_significant_digit_time_tolerance(value: float) -> float:
    if value == 0.0:
        return 0.51e-6
    exponent = math.floor(math.log10(abs(value)))
    return 0.51 * 10.0 ** (exponent - 5)


def add_control_time_checks(
    report: dict[str, Any],
    manifest: dict[str, Any],
    rebuilt: dict[str, dict[str, Any]],
    archived: dict[str, dict[str, Any]],
) -> None:
    archived_times = [float(archived[product_id]["time"]) for product_id in EXPECTED_CONTROLS]
    archive_delta = max(archived_times) - min(archived_times)
    add_check(
        report,
        "archived_control_times_agree",
        archive_delta == 0.0,
        metrics={"maximum_time_difference": archive_delta, "times": archived_times},
        limits={"maximum_time_difference": 0.0},
    )

    source_time = float(manifest["source_time"])
    tolerance = six_significant_digit_time_tolerance(source_time)
    maximum_delta = max(abs(time - source_time) for time in archived_times)
    add_check(
        report,
        f"rebuilt_source_time_matches_archived_sequence_{ARCHIVED_SEQUENCE}",
        maximum_delta <= tolerance,
        metrics={
            "rebuilt_source_time": source_time,
            "archived_time": archived_times[0],
            "absolute_time_difference": maximum_delta,
        },
        limits={"absolute_time_difference_max": tolerance},
    )

    output_time = float(manifest.get("output_time", source_time))
    rebuilt_delta = max(
        abs(float(rebuilt[product_id]["time"]) - output_time)
        for product_id in EXPECTED_CONTROLS
    )
    add_check(
        report,
        "required_rebuilt_control_times_match_manifest",
        rebuilt_delta == 0.0,
        metrics={
            "manifest_output_time": output_time,
            "maximum_absolute_time_error": rebuilt_delta,
        },
        limits={"maximum_absolute_time_error": 0.0},
    )


def add_control_comparisons(
    report: dict[str, Any],
    rebuilt: dict[str, dict[str, Any]],
    archived: dict[str, dict[str, Any]],
) -> None:
    rebuilt29 = rebuilt["output29"]["pdf"]
    archived29 = archived["output29"]["pdf"]
    add_comparison(
        report,
        "output29_radial_marginal",
        np.sum(rebuilt29, axis=0, dtype=np.float64),
        np.sum(archived29, axis=0, dtype=np.float64),
    )

    rebuilt31 = rebuilt["output31"]["pdf"]
    archived31 = archived["output31"]["pdf"]
    add_comparison(report, "output31_full_4d", rebuilt31, archived31)
    add_comparison(
        report,
        "output31_geometry_costheta_radius",
        np.sum(rebuilt31, axis=(2, 3), dtype=np.float64),
        np.sum(archived31, axis=(2, 3), dtype=np.float64),
    )
    add_comparison(
        report,
        "output31_temperature_1d",
        np.sum(rebuilt31, axis=(0, 1, 3), dtype=np.float64),
        np.sum(archived31, axis=(0, 1, 3), dtype=np.float64),
    )
    add_comparison(
        report,
        "output31_radial_velocity_1d",
        np.sum(rebuilt31, axis=(0, 1, 2), dtype=np.float64),
        np.sum(archived31, axis=(0, 1, 2), dtype=np.float64),
    )
    rebuilt_coarse = coarsen_interior_pairs(
        coarsen_interior_pairs(rebuilt31, axis=3), axis=2
    )
    archived_coarse = coarsen_interior_pairs(
        coarsen_interior_pairs(archived31, axis=3), axis=2
    )
    add_comparison(
        report,
        "output31_coarsened_temperature_velocity",
        rebuilt_coarse,
        archived_coarse,
    )


def normalized_identity_error(
    products: dict[str, dict[str, Any]], total_id: str, component_ids: tuple[str, ...]
) -> tuple[float, float, float]:
    total = products[total_id]["total"]
    components = sum(products[product_id]["total"] for product_id in component_ids)
    scale = max(
        products[total_id]["abs_sum"],
        sum(products[product_id]["abs_sum"] for product_id in component_ids),
        np.finfo(np.float64).tiny,
    )
    return abs(total - components) / scale, total, components


def add_global_product_checks(
    report: dict[str, Any],
    products: dict[str, dict[str, Any]],
    read_pdf: Callable[..., dict[str, Any]],
) -> None:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for product in products.values():
        grouped[product["weight"]].append(product)

    for weight, members in sorted(grouped.items()):
        name = f"same_weight_totals_{weight}"
        if len(members) < 2:
            add_skip(report, name, "fewer than two rebuilt products use this weight")
            continue
        totals = [member["total"] for member in members]
        scale = max(
            max(member["abs_sum"] for member in members),
            np.finfo(np.float64).tiny,
        )
        relative_spread = (max(totals) - min(totals)) / scale
        add_check(
            report,
            name,
            relative_spread <= GLOBAL_REL_TOL,
            metrics={
                "relative_spread": relative_spread,
                "minimum_total": min(totals),
                "maximum_total": max(totals),
                "products": [member["id"] for member in members],
            },
            limits={"relative_spread_max": GLOBAL_REL_TOL},
        )

    volume_products = grouped.get("volume", [])
    if volume_products:
        volume = volume_products[0]["total"]
        relative_error = abs(volume - EXPECTED_DOMAIN_VOLUME) / EXPECTED_DOMAIN_VOLUME
        add_check(
            report,
            "domain_volume_is_400_cubed",
            relative_error <= GLOBAL_REL_TOL,
            metrics={"volume": volume, "relative_error": relative_error},
            limits={
                "expected_volume": EXPECTED_DOMAIN_VOLUME,
                "relative_error_max": GLOBAL_REL_TOL,
            },
        )
    else:
        add_skip(report, "domain_volume_is_400_cubed", "no volume-weighted product exists")

    identities = (
        (
            "mdot_equals_out_plus_in",
            "output28",
            ("science_r_theta_mdot_out", "science_r_theta_mdot_in"),
        ),
        (
            "edot_equals_out_plus_in",
            "output29",
            ("science_r_theta_edot_out", "science_r_theta_edot_in"),
        ),
        (
            "edot_equals_kinetic_plus_thermal",
            "output29",
            ("science_r_theta_edot_kin", "science_r_theta_edot_th"),
        ),
    )
    for name, total_id, component_ids in identities:
        required = (total_id, *component_ids)
        missing = [product_id for product_id in required if product_id not in products]
        if missing:
            add_skip(report, name, f"missing products: {', '.join(missing)}")
            continue
        relative_error, total, components = normalized_identity_error(
            products, total_id, component_ids
        )
        add_check(
            report,
            name,
            relative_error <= GLOBAL_REL_TOL,
            metrics={
                "relative_error": relative_error,
                "total": total,
                "component_sum": components,
            },
            limits={"relative_error_max": GLOBAL_REL_TOL},
            details={"total_product": total_id, "component_products": list(component_ids)},
        )

    geometry_product = next(
        (
            product
            for product in products.values()
            if product["weight"] == "volume"
            and "coord_r" in product["header"]["variables"]
            and "coord_abscostheta" in product["header"]["variables"]
        ),
        None,
    )
    if geometry_product is None:
        add_skip(
            report,
            "volume_weighted_coordinate_invariants",
            "no volume product contains coord_r and coord_abscostheta",
        )
    else:
        loaded = load_pdf_array(
            geometry_product["payload_path"],
            read_pdf,
            label=geometry_product["id"],
        )
        values = loaded["pdf"]
        variables = geometry_product["header"]["variables"]
        radius_axis = variables.index("coord_r")
        angle_axis = variables.index("coord_abscostheta")
        axes_except_radius = tuple(axis for axis in range(values.ndim) if axis != radius_axis)
        axes_except_angle = tuple(axis for axis in range(values.ndim) if axis != angle_axis)
        radius = np.sum(values, axis=axes_except_radius, dtype=np.float64)
        angle = np.sum(values, axis=axes_except_angle, dtype=np.float64)
        total = float(np.sum(values, dtype=np.float64))
        if total <= 0.0:
            raise ValidationError(
                f"volume coordinate invariant product has non-positive total: "
                f"{geometry_product['id']}"
            )
        radius_underflow_fraction = abs(float(radius[0])) / total
        radius_populated = int(np.count_nonzero(radius[1:-1]))
        angle_underflow_fraction = abs(float(angle[0])) / total
        angle_overflow_fraction = abs(float(angle[-1])) / total
        angle_populated = int(np.count_nonzero(angle[1:-1]))
        angle_interior_bins = int(
            geometry_product["header"]["dimensions"][angle_axis]["nbin"]
        )
        passed = (
            radius_underflow_fraction < 1.0e-6
            and radius_populated >= 16
            and angle_underflow_fraction <= 1.0e-12
            and angle_overflow_fraction <= 1.0e-10
            and angle_populated == angle_interior_bins
        )
        add_check(
            report,
            "volume_weighted_coordinate_invariants",
            passed,
            metrics={
                "product": geometry_product["id"],
                "radius_underflow_fraction": radius_underflow_fraction,
                "radius_populated_interior_bins": radius_populated,
                "abs_costheta_underflow_fraction": angle_underflow_fraction,
                "abs_costheta_overflow_fraction": angle_overflow_fraction,
                "abs_costheta_populated_interior_bins": angle_populated,
                "abs_costheta_declared_interior_bins": angle_interior_bins,
            },
            limits={
                "radius_underflow_fraction_strict_max": 1.0e-6,
                "radius_populated_interior_bins_min": 16,
                "abs_costheta_underflow_fraction_max": 1.0e-12,
                "abs_costheta_overflow_fraction_max": 1.0e-10,
                "abs_costheta_populated_interior_bins_required": angle_interior_bins,
            },
        )


def add_required_control_geometry_checks(
    report: dict[str, Any], rebuilt: dict[str, dict[str, Any]]
) -> None:
    output29 = rebuilt["output29"]["pdf"]
    absolute29 = np.abs(output29)
    angle_weight = np.sum(absolute29, axis=1, dtype=np.float64)
    radius_weight = np.sum(absolute29, axis=0, dtype=np.float64)
    total_abs = float(np.sum(absolute29, dtype=np.float64))
    if total_abs <= 0.0:
        raise ValidationError("rebuilt output29 has no nonzero absolute weight")
    angle_underflow = float(angle_weight[0]) / total_abs
    angle_overflow = float(angle_weight[-1]) / total_abs
    angle_populated = int(np.count_nonzero(angle_weight[1:-1]))
    radius_populated = int(np.count_nonzero(radius_weight[1:-1]))
    add_check(
        report,
        "output29_rebuilt_coordinate_population",
        angle_underflow <= 1.0e-12
        and angle_overflow <= 1.0e-10
        and angle_populated == 16
        and radius_populated >= 16,
        metrics={
            "abs_costheta_underflow_abs_weight_fraction": angle_underflow,
            "abs_costheta_overflow_abs_weight_fraction": angle_overflow,
            "abs_costheta_populated_interior_bins": angle_populated,
            "radius_populated_interior_bins": radius_populated,
        },
        limits={
            "abs_costheta_underflow_abs_weight_fraction_max": 1.0e-12,
            "abs_costheta_overflow_abs_weight_fraction_max": 1.0e-10,
            "abs_costheta_populated_interior_bins_required": 16,
            "radius_populated_interior_bins_min": 16,
        },
    )

    output31 = rebuilt["output31"]["pdf"]
    radius_weight31 = np.sum(np.abs(output31), axis=(0, 2, 3), dtype=np.float64)
    radius_populated31 = int(np.count_nonzero(radius_weight31[1:-1]))
    add_check(
        report,
        "output31_rebuilt_radius_population",
        radius_populated31 >= 16,
        metrics={"radius_populated_interior_bins": radius_populated31},
        limits={"radius_populated_interior_bins_min": 16},
    )


def validate(args: argparse.Namespace, report: dict[str, Any]) -> None:
    rebuilt_root = args.rebuilt_output_dir.resolve()
    archived_root = args.archived_phase2_root.resolve()
    if not rebuilt_root.is_dir():
        raise ValidationError(f"rebuilt output directory does not exist: {rebuilt_root}")
    if not archived_root.is_dir():
        raise ValidationError(f"archived phase2 root does not exist: {archived_root}")

    read_pdf, read_pdf_header, analysis_root = load_pdf_reader()
    report["analysis_root"] = str(analysis_root)
    manifest, products = read_manifest_and_products(
        rebuilt_root, read_pdf_header, report
    )
    add_manifest_checks(report, manifest, archived_root)
    if args.no_controls:
        add_skip(
            report,
            "archived_output29_output31_controls",
            "snapshot is explicitly declared to have no archived PDF control",
        )
        add_global_product_checks(report, products, read_pdf)
        return

    missing = [product_id for product_id in EXPECTED_CONTROLS if product_id not in products]
    if missing:
        raise ValidationError(f"rebuilt output is missing required products: {', '.join(missing)}")

    rebuilt_controls: dict[str, dict[str, Any]] = {}
    archived_controls: dict[str, dict[str, Any]] = {}
    archived_sources: dict[str, Any] = {}
    for product_id, stream in ARCHIVED_STREAMS.items():
        validate_control_header(product_id, products[product_id]["header"], "rebuilt")
        rebuilt_controls[product_id] = load_pdf_array(
            products[product_id]["payload_path"],
            read_pdf,
            label=f"rebuilt {product_id}",
        )

        archived_path = find_archived_payload(archived_root, stream)
        archived_controls[product_id] = load_pdf_array(
            archived_path,
            read_pdf,
            label=f"archived {product_id}",
        )
        validate_control_header(
            product_id, archived_controls[product_id]["header"], "archived"
        )
        source_count = len(archived_controls[product_id]["source_files"])
        add_check(
            report,
            f"{product_id}_archived_sparse_shard_count",
            source_count == EXPECTED_SHARDS,
            metrics={"source_files": source_count},
            limits={"required": EXPECTED_SHARDS},
        )
        add_header_alignment_check(
            report,
            product_id,
            rebuilt_controls[product_id]["header"],
            archived_controls[product_id]["header"],
        )
        archived_sources[product_id] = {
            "representative_file": archived_controls[product_id]["representative_file"],
            "header_file": archived_controls[product_id]["header_file"],
            "source_files": source_count,
            "time": archived_controls[product_id]["time"],
        }
    report["archived_controls"] = archived_sources

    add_control_time_checks(report, manifest, rebuilt_controls, archived_controls)
    add_control_comparisons(report, rebuilt_controls, archived_controls)
    add_required_control_geometry_checks(report, rebuilt_controls)
    add_global_product_checks(report, products, read_pdf)


def main() -> int:
    args = parse_args()
    configure_expectations(args)
    report = new_report(args)
    try:
        validate(args, report)
    except Exception as exc:
        report["fatal_error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
        }
    finalize_report(report)
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.report_path is not None:
        report_path = args.report_path.resolve()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        partial = report_path.with_name(report_path.name + ".partial")
        partial.write_text(encoded, encoding="utf-8")
        os.chmod(partial, 0o444)
        os.replace(partial, report_path)
    print(encoded, end="")
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
