#!/usr/bin/env python3
"""Require the Phase-4b coupled analyzer to reject deliberately broken criteria."""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path


def _run(command, label):
    print(f"\n=== mutation: {label} ===")
    result = subprocess.run(command, text=True, capture_output=True)
    print(result.stdout, end="")
    print(result.stderr, end="")
    print(f"exit_code={result.returncode}")
    if result.returncode == 0:
        raise RuntimeError(f"{label}: analyzer unexpectedly accepted mutated criteria")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analyzer", required=True, type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--prescribed-input", required=True, type=Path)
    parser.add_argument("--criteria", required=True, type=Path)
    args = parser.parse_args()

    base_command = [
        "python3", str(args.analyzer), "--binary", str(args.binary),
        "--input", str(args.input), "--prescribed-input", str(args.prescribed_input),
    ]
    criteria = json.loads(args.criteria.read_text())
    with tempfile.TemporaryDirectory(prefix="cr-rel-phase4b-oracle-") as root:
        root = Path(root)

        tight = json.loads(json.dumps(criteria))
        for item in tight["criteria"]:
            if item["id"] == "P4B-14":
                item["limit"] = -1.0
                break
        tight_path = root / "tight.json"
        tight_path.write_text(json.dumps(tight, indent=2) + "\n")
        _run(base_command + [
            "--criteria", str(tight_path), "--work-dir", str(root / "tight")], "tight")

        missing = json.loads(json.dumps(criteria))
        missing["criteria"] = [
            item for item in missing["criteria"] if item["id"] != "P4B-R08"]
        missing_path = root / "missing.json"
        missing_path.write_text(json.dumps(missing, indent=2) + "\n")
        _run(base_command + [
            "--criteria", str(missing_path), "--work-dir", str(root / "missing")],
            "missing")

        _run(base_command + [
            "--criteria", str(args.criteria), "--work-dir", str(root / "source-drift"),
            "--manufactured-oracle-acceleration-scale", "0.0"], "source-drift")

        _run(base_command + [
            "--criteria", str(args.criteria), "--work-dir", str(root / "right-end"),
            "--manufactured-one-step-oracle-time-level", "right_end"], "right-end")

    print("\nmutation negative controls rejected as expected")


if __name__ == "__main__":
    main()
