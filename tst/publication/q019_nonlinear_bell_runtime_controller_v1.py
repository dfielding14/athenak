#!/usr/bin/env python3
"""Build sealed Q019 runtime-controller regression and excluded-pilot overlays."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Mapping

from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as q019


REPO_ROOT = Path(__file__).resolve().parents[2]
CONTROLLER_BLOCK = "q019_nonlinear_bell_runtime_controller_v1"
CHECKED_IN_ROOT = REPO_ROOT / "inputs/publication/q019_nonlinear_bell_runtime_controller_v1"
CHECKED_IN_MANIFEST = CHECKED_IN_ROOT / "deck_manifest.json"


def _float_token(value: float) -> str:
    return format(value, ".17g")


def controller_identity_payload(parameters: Mapping[str, str]) -> str:
    entries = sorted(
        (name, value)
        for name, value in parameters.items()
        if name != "controller_identity_fingerprint" and not name.startswith("runtime_")
    )
    return "".join(
        f"{CONTROLLER_BLOCK}/{name}={value}\n" for name, value in entries
    )


def controller_identity_fingerprint(parameters: Mapping[str, str]) -> str:
    return hashlib.sha256(controller_identity_payload(parameters).encode("utf-8")).hexdigest()


def _controller_parameters(
    *,
    source_case_id: str,
    authority: str,
    monitor_dt: float,
    box_edge_monitor_enabled: bool,
    resolution_monitor_enabled: bool,
    diagnostic_failure_stop_armed: bool,
    resolution_stop_armed: bool,
    resolution_stop_b_over_b0: float,
    box_edge_stop_armed: bool,
    box_edge_stop_ppm: int,
    pilot_cycle_limit: int,
) -> dict[str, str]:
    if authority not in {"runtime_regression_only", "excluded_pilot_only"}:
        raise ValueError(f"unsupported controller authority: {authority}")
    parameters = {
        "schema": "1",
        "authority": authority,
        "contract_id": (
            "q019-runtime-controller-v1-regression"
            if authority == "runtime_regression_only"
            else "q019-runtime-controller-v1-excluded-pilot"
        ),
        "cadence_status": (
            "accelerated_runtime_regression"
            if authority == "runtime_regression_only"
            else "excluded_pilot_candidate_not_frozen"
        ),
        "stop_disposition": "accepted_guard_stop_or_diagnostic_failure",
        "source_case_id": source_case_id,
        "monitor_dt": _float_token(monitor_dt),
        "box_edge_monitor_enabled": str(box_edge_monitor_enabled).lower(),
        "resolution_monitor_enabled": str(resolution_monitor_enabled).lower(),
        "diagnostic_failure_stop_armed": str(
            diagnostic_failure_stop_armed
        ).lower(),
        "resolution_stop_armed": str(resolution_stop_armed).lower(),
        "resolution_stop_B_over_B0": _float_token(resolution_stop_b_over_b0),
        "box_edge_stop_armed": str(box_edge_stop_armed).lower(),
        "box_edge_stop_ppm": str(box_edge_stop_ppm),
        "pilot_cycle_limit": str(pilot_cycle_limit),
    }
    parameters["controller_identity_fingerprint"] = controller_identity_fingerprint(
        parameters
    )
    return parameters


def expected_overlays() -> tuple[dict[str, object], ...]:
    definitions = (
        {
            "artifact_id": "q019-controller-regression-box-stop",
            "source_case_id": "q019-fr-runtime-initializer-ppc24-s0",
            "authority": "runtime_regression_only",
            "monitor_dt": 1.0e-6,
            "box_edge_monitor_enabled": True,
            "resolution_monitor_enabled": True,
            "diagnostic_failure_stop_armed": True,
            "resolution_stop_armed": False,
            "resolution_stop_b_over_b0": -1.0,
            "box_edge_stop_armed": True,
            "box_edge_stop_ppm": 1,
            "pilot_cycle_limit": 2,
            "expected_stop_reason": 1902,
        },
        {
            "artifact_id": "q019-controller-regression-cycle-stop",
            "source_case_id": "q019-fr-runtime-initializer-ppc24-s0",
            "authority": "runtime_regression_only",
            "monitor_dt": 1.0e-6,
            "box_edge_monitor_enabled": False,
            "resolution_monitor_enabled": False,
            "diagnostic_failure_stop_armed": False,
            "resolution_stop_armed": False,
            "resolution_stop_b_over_b0": -1.0,
            "box_edge_stop_armed": False,
            "box_edge_stop_ppm": -1,
            "pilot_cycle_limit": 1,
            "expected_stop_reason": 1903,
        },
        {
            "artifact_id": "q019-controller-pilot-2d-baseline",
            "source_case_id": "q019-fr-grid-k8-rho1em05-s0",
            "authority": "excluded_pilot_only",
            "monitor_dt": 1.0e-6,
            "box_edge_monitor_enabled": False,
            "resolution_monitor_enabled": False,
            "diagnostic_failure_stop_armed": False,
            "resolution_stop_armed": False,
            "resolution_stop_b_over_b0": -1.0,
            "box_edge_stop_armed": False,
            "box_edge_stop_ppm": -1,
            "pilot_cycle_limit": 20,
            "expected_stop_reason": 1903,
        },
        {
            "artifact_id": "q019-controller-pilot-2d-instrumented",
            "source_case_id": "q019-fr-grid-k8-rho1em05-s0",
            "authority": "excluded_pilot_only",
            "monitor_dt": 1.0e-6,
            "box_edge_monitor_enabled": True,
            "resolution_monitor_enabled": True,
            "diagnostic_failure_stop_armed": True,
            "resolution_stop_armed": False,
            "resolution_stop_b_over_b0": -1.0,
            "box_edge_stop_armed": False,
            "box_edge_stop_ppm": -1,
            "pilot_cycle_limit": 20,
            "expected_stop_reason": 1903,
        },
        {
            "artifact_id": "q019-controller-pilot-3d-baseline",
            "source_case_id": "q019-fr-3d-onset-small-s0",
            "authority": "excluded_pilot_only",
            "monitor_dt": 1.0e-5,
            "box_edge_monitor_enabled": False,
            "resolution_monitor_enabled": False,
            "diagnostic_failure_stop_armed": False,
            "resolution_stop_armed": False,
            "resolution_stop_b_over_b0": -1.0,
            "box_edge_stop_armed": False,
            "box_edge_stop_ppm": -1,
            "pilot_cycle_limit": 20,
            "expected_stop_reason": 1903,
        },
        {
            "artifact_id": "q019-controller-pilot-3d-instrumented",
            "source_case_id": "q019-fr-3d-onset-small-s0",
            "authority": "excluded_pilot_only",
            "monitor_dt": 1.0e-5,
            "box_edge_monitor_enabled": True,
            "resolution_monitor_enabled": True,
            "diagnostic_failure_stop_armed": True,
            "resolution_stop_armed": False,
            "resolution_stop_b_over_b0": -1.0,
            "box_edge_stop_armed": False,
            "box_edge_stop_ppm": -1,
            "pilot_cycle_limit": 20,
            "expected_stop_reason": 1903,
        },
    )
    overlays: list[dict[str, object]] = []
    for definition in definitions:
        parameters = _controller_parameters(
            source_case_id=str(definition["source_case_id"]),
            authority=str(definition["authority"]),
            monitor_dt=float(definition["monitor_dt"]),
            box_edge_monitor_enabled=bool(definition["box_edge_monitor_enabled"]),
            resolution_monitor_enabled=bool(definition["resolution_monitor_enabled"]),
            diagnostic_failure_stop_armed=bool(
                definition["diagnostic_failure_stop_armed"]
            ),
            resolution_stop_armed=bool(definition["resolution_stop_armed"]),
            resolution_stop_b_over_b0=float(
                definition["resolution_stop_b_over_b0"]
            ),
            box_edge_stop_armed=bool(definition["box_edge_stop_armed"]),
            box_edge_stop_ppm=int(definition["box_edge_stop_ppm"]),
            pilot_cycle_limit=int(definition["pilot_cycle_limit"]),
        )
        overlays.append({**definition, "controller_parameters": parameters})
    return tuple(overlays)


def render_overlay(overlay: Mapping[str, object]) -> str:
    source_case_id = str(overlay["source_case_id"])
    base = (q019.CHECKED_IN_DECK_ROOT / f"{source_case_id}.athinput").read_text(
        encoding="utf-8"
    )
    parameters = dict(overlay["controller_parameters"])
    lines = ["", f"<{CONTROLLER_BLOCK}>"]
    lines.extend(f"{name} = {value}" for name, value in parameters.items())
    return base.rstrip() + "\n" + "\n".join(lines) + "\n"


def build_manifest() -> dict[str, object]:
    artifacts = []
    for overlay in expected_overlays():
        source_case_id = str(overlay["source_case_id"])
        source_path = q019.CHECKED_IN_DECK_ROOT / f"{source_case_id}.athinput"
        rendered = render_overlay(overlay)
        artifacts.append(
            {
                **overlay,
                "filename": f"{overlay['artifact_id']}.athinput",
                "source_deck": str(source_path.relative_to(REPO_ROOT)),
                "source_deck_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                "rendered_sha256": hashlib.sha256(rendered.encode("utf-8")).hexdigest(),
                "saturation_evidence_eligible": False,
                "production_promotion_authorized": False,
            }
        )
    return {
        "schema_version": 1,
        "record_type": "q019_nonlinear_bell_runtime_controller_packet",
        "controller_block": CONTROLLER_BLOCK,
        "authority": "runtime_regression_and_excluded_pilot_only",
        "saturation_evidence_eligible": False,
        "production_promotion_authorized": False,
        "artifacts": artifacts,
    }


def materialize(root: Path = CHECKED_IN_ROOT) -> None:
    root.mkdir(parents=True, exist_ok=True)
    manifest = build_manifest()
    expected_names = set()
    for overlay in expected_overlays():
        filename = f"{overlay['artifact_id']}.athinput"
        expected_names.add(filename)
        (root / filename).write_text(render_overlay(overlay), encoding="utf-8")
    for path in root.glob("*.athinput"):
        if path.name not in expected_names:
            path.unlink()
    (root / "deck_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=CHECKED_IN_ROOT)
    args = parser.parse_args()
    materialize(args.output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
