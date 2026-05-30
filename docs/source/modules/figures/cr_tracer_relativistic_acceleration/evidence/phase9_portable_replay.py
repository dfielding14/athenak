#!/usr/bin/env python3
"""Replay strict all-format validation from an extracted Phase-9 package."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _verify_manifest(root: Path) -> None:
    manifest = root / "portable_manifest.sha256"
    for line in manifest.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = root / relative
        actual = _sha256(path)
        if actual != expected:
            raise RuntimeError(
                f"portable manifest mismatch for {relative}: {actual} != {expected}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    _verify_manifest(root)

    binary = (root / "qualified_binary" / "athena").resolve()
    runtime_root = (root / "runtime_root" / "acceleration_uninterrupted").resolve()
    input_path = (
        root / "inputs" / "cr_relativistic_diagnostics_runtime.athinput").resolve()
    criteria = (
        root / "criteria"
        / "cr_relativistic_phase7_preregistered_criteria.json").resolve()
    metrics = json.loads(
        (root / "metrics" / "diagnostics_runtime_metrics_template.json").read_text())
    metrics["binary"] = str(binary)
    metrics["input"] = str(input_path)
    metrics["criteria"] = str(criteria)
    metrics["runtime_cases"] = [str(runtime_root)]
    relocated_metrics = root / "metrics" / "diagnostics_runtime_metrics_relocated.json"
    relocated_metrics.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")

    output = args.output or (root / "portable_all_formats_report.json")
    command = [
        sys.executable,
        str(root / "scripts" / "particles" / "cr_relativistic_all_formats_inspect.py"),
        "--runtime-root", str(runtime_root),
        "--binary", str(binary),
        "--input", str(input_path),
        "--criteria", str(criteria),
        "--metrics", str(relocated_metrics),
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    output.write_text(result.stdout.replace(str(root), "$PACKAGE_ROOT"))
    print("CR relativistic Phase-9 portable all-format replay PASS")
    print("  package root: $PACKAGE_ROOT")
    print("  report: portable_all_formats_report.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
