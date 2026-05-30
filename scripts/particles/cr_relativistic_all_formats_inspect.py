#!/usr/bin/env python3
"""Export deterministic Phase-7 relativistic CR all-format evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Optional
from typing import Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
GENERATOR = "scripts/particles/cr_relativistic_all_formats_inspect.py"
QUALIFICATION = "cr_relativistic_diagnostics_runtime_phase7"
sys.path.insert(0, str(SCRIPT_DIR))
import cr_relativistic_diagnostics_runtime_inspect as diagnostics  # noqa: E402
import cr_tracer_inspect as tracer_inspect  # noqa: E402


def _require_file(path: Path, label: str) -> None:
    diagnostics._require(path.is_file(), f"{label} is not a file: {path}")


def _require_matching_path(report: dict, key: str, expected: Path) -> None:
    value = report.get(key)
    diagnostics._require(
        isinstance(value, str) and Path(value).resolve() == expected,
        f"metrics {key}={value!r}, expected {str(expected)!r}",
    )


def _require_matching_sha256(report: dict, key: str, expected: Path) -> str:
    digest = hashlib.sha256(expected.read_bytes()).hexdigest()
    diagnostics._require(
        report.get(key) == digest,
        f"metrics {key}={report.get(key)!r}, expected {digest!r}",
    )
    return digest


def _validate_replay_provenance(metrics_path: Path, binary: Path,
                                input_path: Path, criteria_path: Path,
                                runtime_root: Path) -> dict[str, str]:
    try:
        report = json.loads(metrics_path.read_text())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RuntimeError(
            f"could not read metrics {metrics_path}: {error}") from error
    diagnostics._require(
        isinstance(report, dict), f"{metrics_path}: expected object")
    diagnostics._require(
        report.get("schema_version") == 1,
        f"{metrics_path}: unsupported metrics schema",
    )
    diagnostics._require(
        report.get("qualification") == QUALIFICATION,
        f"{metrics_path}: unexpected qualification "
        f"{report.get('qualification')!r}",
    )
    diagnostics._require(
        report.get("status") == "PASS",
        f"{metrics_path}: qualification status is not PASS",
    )
    _require_matching_path(report, "binary", binary)
    _require_matching_path(report, "input", input_path)
    _require_matching_path(report, "criteria", criteria_path)
    sha256 = {
        "binary": _require_matching_sha256(report, "binary_sha256", binary),
        "input": _require_matching_sha256(report, "input_sha256", input_path),
        "criteria": _require_matching_sha256(
            report, "criteria_sha256", criteria_path),
    }
    runtime_cases = report.get("runtime_cases")
    diagnostics._require(
        isinstance(runtime_cases, list)
        and all(isinstance(path, str) for path in runtime_cases),
        f"{metrics_path}: runtime_cases must be a list of paths",
    )
    diagnostics._require(
        runtime_root in {Path(path).resolve() for path in runtime_cases},
        f"{metrics_path}: runtime root {runtime_root} is not a qualified case",
    )
    metrics = report.get("metrics")
    diagnostics._require(
        isinstance(metrics, dict),
        f"{metrics_path}: metrics must be an object",
    )
    diagnostics._validate_criteria(criteria_path, metrics)
    return sha256


def _retained_metadata(runtime_root: Path,
                       checkpoints: Sequence[diagnostics.Checkpoint]) -> dict:
    diagnostics._retained_output_dictionary(runtime_root, checkpoints)

    ppd = tracer_inspect.summarize_ppd(runtime_root)
    diagnostics._require(
        ppd is not None, f"{runtime_root}: missing ppd output")
    histograms = tracer_inspect.validate_histograms(runtime_root, [1], 8)
    tracer_inspect.validate_moments(runtime_root, [1])
    pmom = tracer_inspect.read_pmom_record(
        diagnostics._latest_output(runtime_root, "pmom", "pmom"), 1)
    joint = tracer_inspect.validate_joint_spectra(runtime_root, [1], 8, 8)
    _, sample = diagnostics._sample(runtime_root)
    spectrum = diagnostics._spectrum(runtime_root)

    return {
        "df": histograms["df_record"]["metadata"],
        "dparh": histograms["dparh"]["metadata"],
        "drh": histograms["drh"]["metadata"],
        "dxh": histograms["dxh_record"]["metadata"],
        "pmom": pmom["metadata"],
        "ppd": ppd["metadata"],
        "psamp": sample["metadata"],
        "pspec": spectrum["metadata"],
        "pspec2": joint["metadata"],
    }


def _restart_inventory(checkpoints: Sequence[diagnostics.Checkpoint]) -> dict:
    artifacts = set()
    formats = Counter()
    shards = []
    for checkpoint in checkpoints:
        shard = checkpoint.shard_relative
        decoded = tracer_inspect.read_prst_file(checkpoint.root / shard)
        diagnostics._require(
            decoded["format_version"] == "AKPRST-v2.0",
            f"{shard}: expected typed-v2 restart shard",
        )
        formats[decoded["format_version"]] += 1
        shards.append(shard.as_posix())
        for relative in (
                checkpoint.manifest_relative, checkpoint.shard_relative,
                checkpoint.mesh_relative, checkpoint.witness_relative):
            artifacts.add(relative.as_posix())
    inspected_shards = sorted(shards)
    return {
        "restart_artifacts": sorted(artifacts),
        "restart_formats_all_checkpoints": dict(sorted(formats.items())),
        "restart_inspected_shard_count": len(inspected_shards),
        "restart_inspected_shards": inspected_shards,
    }


def _inspect(runtime_root: Path, binary: Path, input_path: Path,
             criteria_path: Path, metrics_path: Path) -> dict:
    diagnostics._require(
        runtime_root.is_dir(),
        f"runtime root is not a directory: {runtime_root}",
    )
    _require_file(binary, "binary")
    _require_file(input_path, "input")
    _require_file(criteria_path, "criteria")
    _require_file(metrics_path, "metrics")
    sha256 = _validate_replay_provenance(
        metrics_path, binary, input_path, criteria_path, runtime_root)
    checkpoints = diagnostics._checkpoints(runtime_root)
    metadata = _retained_metadata(runtime_root, checkpoints)
    inventory = _restart_inventory(checkpoints)
    return {
        **metadata,
        "provenance": {
            "binary": str(binary),
            "binary_sha256": sha256["binary"],
            "criteria": str(criteria_path),
            "criteria_sha256": sha256["criteria"],
            "generator": GENERATOR,
            "input": str(input_path),
            "input_sha256": sha256["input"],
            "metrics": str(metrics_path),
            "runtime_root": str(runtime_root),
        },
        **inventory,
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Export deterministic Phase-7 relativistic CR all-format "
            "evidence."))
    parser.add_argument("--runtime-root", required=True, type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--criteria", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        report = _inspect(
            args.runtime_root.resolve(),
            args.binary.resolve(),
            args.input.resolve(),
            args.criteria.resolve(),
            args.metrics.resolve(),
        )
    except Exception as error:
        print(f"CR relativistic Phase-7 all-format inspection FAIL: {error}",
              file=sys.stderr)
        return 1
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
