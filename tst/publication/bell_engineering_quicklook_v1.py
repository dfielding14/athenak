#!/usr/bin/env python3
"""Raw-output quicklooks for compact Bell engineering runs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Any

import numpy as np

from tst.publication import analyze_q011_section54_outputs as binary
from tst.publication import analyze_q023_paper_bell_linear_joverc as linear
from tst.publication import q019_registered_raw_reduction_v1 as nonlinear_raw


LINEAR_MEMBER_ID = "d1-fine-reference_x1-eps0p4"
NONLINEAR_CASE_ID = "q019-q023-carrier-s2-resolution-coarse-s0"
_OUTPUT_RE = re.compile(r"[.]([0-9]+)[.]bin$")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact(path: Path, root: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(root)),
        "sha256": _sha256(path),
        "byte_count": path.stat().st_size,
    }


def analyze_linear(run_root: Path) -> dict[str, Any]:
    member = next(
        item
        for item in linear.expected_deck_members()
        if item["member_id"] == LINEAR_MEMBER_ID
    )
    basename = "q023_joverc_" + LINEAR_MEMBER_ID.replace("-", "_")
    paths = sorted((run_root / "bin").glob(f"{basename}.mhd_w_bcc.*.bin"))
    if len(paths) < 8:
        raise RuntimeError(f"linear quicklook needs at least 8 snapshots, found {len(paths)}")
    datasets = [binary.read_athenak_binary(path) for path in paths]
    trace = linear.physics_trace_from_raw_datasets(datasets, member=member)
    physics = linear._analyze_physics_trace(trace, float(member["epsilon"]))
    return {
        "record_type": "q023_linear_bell_engineering_quicklook_v1",
        "case_id": LINEAR_MEMBER_ID,
        "run_root": str(run_root),
        "snapshot_count": len(paths),
        "physics": physics,
        "engineering_ready": bool(physics["passed"]),
        "claim_scope": "single_case_engineering_evidence_not_registered_qualification",
        "raw_artifacts": [_artifact(path, run_root) for path in paths],
    }


def _output_index(path: Path) -> int:
    match = _OUTPUT_RE.search(path.name)
    if match is None:
        raise RuntimeError(f"cannot parse output index: {path}")
    return int(match.group(1))


def analyze_nonlinear(run_root: Path) -> dict[str, Any]:
    products_by_index: dict[int, dict[str, Path]] = {}
    for product in nonlinear_raw.REQUIRED_BINARY_PRODUCTS:
        for path in sorted((run_root / "bin").glob(f"{NONLINEAR_CASE_ID}.{product}.*.bin")):
            products_by_index.setdefault(_output_index(path), {})[product] = path
    complete = [
        (index, products)
        for index, products in sorted(products_by_index.items())
        if set(products) == set(nonlinear_raw.REQUIRED_BINARY_PRODUCTS)
    ]
    if len(complete) < 2:
        raise RuntimeError(
            f"nonlinear quicklook needs at least 2 complete ten-product snapshots, found {len(complete)}"
        )

    snapshots = []
    raw_paths: list[Path] = []
    for index, paths in complete:
        raw_paths.extend(paths.values())
        snapshot = nonlinear_raw.compose_snapshot(
            NONLINEAR_CASE_ID,
            {name: path.read_bytes() for name, path in paths.items()},
        )
        fields = snapshot["fields"]
        b1 = np.asarray(fields["bcc1"], dtype=np.float64)
        b2 = np.asarray(fields["bcc2"], dtype=np.float64)
        b3 = np.asarray(fields["bcc3"], dtype=np.float64)
        density = np.asarray(fields["dens"], dtype=np.float64)
        snapshots.append(
            {
                "output_index": index,
                "cycle": int(snapshot["cycle"]),
                "time": float(snapshot["time"]),
                "mean_b1": float(np.mean(b1)),
                "bperp_rms": float(np.sqrt(np.mean(b2 * b2 + b3 * b3))),
                "btotal_rms": float(np.sqrt(np.mean(b1 * b1 + b2 * b2 + b3 * b3))),
                "density_std_over_mean": float(np.std(density) / np.mean(density)),
                "mean_deposited_jx": float(np.mean(fields["prtcl_jx"])),
            }
        )

    b0 = abs(snapshots[0]["mean_b1"])
    if not np.isfinite(b0) or b0 <= 0.0:
        raise RuntimeError("initial mean B1 is not positive and finite")
    for snapshot in snapshots:
        snapshot["bperp_rms_over_B0"] = snapshot["bperp_rms"] / b0
        snapshot["btotal_rms_over_B0"] = snapshot["btotal_rms"] / b0
    maximum = max(item["bperp_rms_over_B0"] for item in snapshots)
    return {
        "record_type": "q019_nonlinear_bell_engineering_quicklook_v1",
        "case_id": NONLINEAR_CASE_ID,
        "run_root": str(run_root),
        "snapshot_count": len(snapshots),
        "initial_B0": b0,
        "maximum_bperp_rms_over_B0": maximum,
        "final_bperp_rms_over_B0": snapshots[-1]["bperp_rms_over_B0"],
        "nonlinear_amplitude_reached": bool(maximum >= 1.0),
        "snapshots": snapshots,
        "claim_scope": (
            "coarse_single_case_engineering_evidence_not_saturation_or_registered_qualification"
        ),
        "raw_artifacts": [_artifact(path, run_root) for path in sorted(raw_paths)],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("linear", "nonlinear"))
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-pass", action="store_true")
    args = parser.parse_args()
    run_root = args.run_root.resolve(strict=True)
    report = analyze_linear(run_root) if args.mode == "linear" else analyze_nonlinear(run_root)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload, encoding="utf-8")
    passed = (
        report["engineering_ready"]
        if args.mode == "linear"
        else report["nonlinear_amplitude_reached"]
    )
    return 0 if (not args.require_pass or passed) else 1


if __name__ == "__main__":
    raise SystemExit(main())
