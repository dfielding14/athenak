#!/usr/bin/env python3
"""Build the excluded six-run Q019 Q023-carrier resource calibration packet."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Mapping

from tst.publication import q019_nonlinear_bell_runtime_controller_v1 as controller
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as q019
from tst.publication import q019_q023_carrier_nonlinear_bell_redesign_v1 as carrier


REPO_ROOT = Path(__file__).resolve().parents[2]
CHECKED_IN_ROOT = (
    REPO_ROOT
    / "inputs/publication/"
    "q019_q023_carrier_resource_calibration_runtime_controller_v1"
)
CHECKED_IN_MANIFEST = CHECKED_IN_ROOT / "deck_manifest.json"
SCHEMA_VERSION = 1
RECORD_TYPE = "q019_q023_carrier_resource_calibration_packet_v1"
SUCCESSOR_ID = "q019_q023_carrier_resource_calibration_runtime_controller_v1"
EXPECTED_ARTIFACT_COUNT = 6
TASKS_PER_NODE = 8
CALIBRATION_NODES = {2: 1, 3: 2}
CALIBRATION_AUTHORITY = {
    "launch_authorized": False,
    "scheduler_submission_authorized": False,
    "policy_mutation_authorized": False,
    "production_resource_freeze_authorized": False,
    "q019_qualification_authorized": False,
    "nonlinear_saturation_claim_authorized": False,
    "scientific_claim_authorized": False,
    "publication_authorized": False,
}


class CalibrationError(ValueError):
    """Reject a drifted or authority-bearing carrier calibration packet."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise CalibrationError(message)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _definition(
    *,
    artifact_id: str,
    source_case_id: str,
    dimension: int,
    cycle_limit: int,
    instrumentation: str,
    calibration_role: str,
) -> dict[str, object]:
    _require(dimension in {2, 3}, "carrier calibration dimension is invalid")
    _require(cycle_limit > 0, "carrier calibration cycle limit is invalid")
    _require(
        instrumentation in {"baseline", "instrumented"},
        "carrier calibration instrumentation role is invalid",
    )
    enabled = instrumentation == "instrumented"
    monitor_dt = 1.0e-6 if dimension == 2 else 1.0e-5
    nodes = CALIBRATION_NODES[dimension]
    parameters = controller._controller_parameters(
        source_case_id=source_case_id,
        authority="excluded_pilot_only",
        monitor_dt=monitor_dt,
        box_edge_monitor_enabled=enabled,
        resolution_monitor_enabled=enabled,
        diagnostic_failure_stop_armed=enabled,
        resolution_stop_armed=False,
        resolution_stop_b_over_b0=-1.0,
        box_edge_stop_armed=False,
        box_edge_stop_ppm=-1,
        pilot_cycle_limit=cycle_limit,
    )
    return {
        "artifact_id": artifact_id,
        "source_case_id": source_case_id,
        "dimension": dimension,
        "cycle_limit": cycle_limit,
        "instrumentation": instrumentation,
        "calibration_role": calibration_role,
        "authority": "excluded_pilot_only",
        "monitor_dt": monitor_dt,
        "expected_stop_reason": 1903,
        "nodes": nodes,
        "tasks": nodes * TASKS_PER_NODE,
        "controller_parameters": parameters,
    }


def expected_overlays() -> tuple[dict[str, object], ...]:
    definitions = (
        (
            "q019-carrier-calibration-2d-onset-c0032-baseline",
            "q019-q023-carrier-s1-onset-s0",
            2,
            32,
            "baseline",
            "startup_and_peak_memory",
        ),
        (
            "q019-carrier-calibration-2d-onset-c0256-baseline",
            "q019-q023-carrier-s1-onset-s0",
            2,
            256,
            "baseline",
            "steady_cycle_cost",
        ),
        (
            "q019-carrier-calibration-2d-onset-c0256-instrumented",
            "q019-q023-carrier-s1-onset-s0",
            2,
            256,
            "instrumented",
            "steady_cycle_cost_and_diagnostic_overhead",
        ),
        (
            "q019-carrier-calibration-3d-fiducial-c0016-baseline",
            "q019-q023-carrier-s3-3d-fiducial-s0",
            3,
            16,
            "baseline",
            "startup_and_peak_memory",
        ),
        (
            "q019-carrier-calibration-3d-fiducial-c0064-baseline",
            "q019-q023-carrier-s3-3d-fiducial-s0",
            3,
            64,
            "baseline",
            "steady_cycle_cost",
        ),
        (
            "q019-carrier-calibration-3d-fiducial-c0064-instrumented",
            "q019-q023-carrier-s3-3d-fiducial-s0",
            3,
            64,
            "instrumented",
            "steady_cycle_cost_and_diagnostic_overhead",
        ),
    )
    return tuple(
        _definition(
            artifact_id=artifact_id,
            source_case_id=source_case_id,
            dimension=dimension,
            cycle_limit=cycle_limit,
            instrumentation=instrumentation,
            calibration_role=calibration_role,
        )
        for (
            artifact_id,
            source_case_id,
            dimension,
            cycle_limit,
            instrumentation,
            calibration_role,
        ) in definitions
    )


def render_overlay(overlay: Mapping[str, object]) -> str:
    source_case_id = str(overlay["source_case_id"])
    source_path = carrier.CHECKED_IN_DECK_ROOT / f"{source_case_id}.athinput"
    _require(source_path.is_file(), f"carrier source deck is absent: {source_case_id}")
    base = source_path.read_text(encoding="utf-8")
    parameters = dict(overlay["controller_parameters"])
    _require(
        parameters.get("source_case_id") == source_case_id
        and parameters.get("controller_identity_fingerprint")
        == controller.controller_identity_fingerprint(parameters),
        "carrier calibration controller identity fingerprint drifted",
    )
    case = next(
        (
            item
            for item in carrier.expected_cases()
            if item["case_id"] == source_case_id
        ),
        None,
    )
    _require(case is not None, f"carrier source case is absent: {source_case_id}")
    q019.validate_rendered_deck(case, base)
    lines = ["", f"<{controller.CONTROLLER_BLOCK}>"]
    lines.extend(f"{name} = {value}" for name, value in parameters.items())
    rendered = base.rstrip() + "\n" + "\n".join(lines) + "\n"
    blocks = q019.parse_athinput_text(rendered)
    _require(
        blocks.get(controller.CONTROLLER_BLOCK) == parameters,
        "carrier calibration controller block drifted while rendering",
    )
    return rendered


def build_manifest() -> dict[str, object]:
    carrier_manifest = carrier.validate_checked_in_decks()
    cases = {str(case["case_id"]): case for case in carrier.expected_cases()}
    artifacts = []
    for overlay in expected_overlays():
        source_case_id = str(overlay["source_case_id"])
        case = cases[source_case_id]
        source_path = carrier.CHECKED_IN_DECK_ROOT / f"{source_case_id}.athinput"
        source_payload = source_path.read_bytes()
        rendered = render_overlay(overlay).encode("utf-8")
        meshblocks = 1
        for global_nx, block_nx in zip(case["nx"], case["meshblock_nx"]):
            meshblocks *= int(global_nx) // int(block_nx)
        _require(
            meshblocks <= int(overlay["tasks"]),
            f"{overlay['artifact_id']}: calibration tasks do not cover meshblocks",
        )
        artifacts.append(
            {
                **overlay,
                "filename": f"{overlay['artifact_id']}.athinput",
                "source_deck": str(source_path.relative_to(REPO_ROOT)),
                "source_deck_sha256": _sha256(source_payload),
                "source_matrix_identity_fingerprint": case[
                    "matrix_identity_fingerprint"
                ],
                "root_cells": int(case["nx"][0])
                * int(case["nx"][1])
                * int(case["nx"][2]),
                "meshblocks": meshblocks,
                "macro_particles": int(case["nx"][0])
                * int(case["nx"][1])
                * int(case["nx"][2])
                * int(case["ppc"]),
                "rendered_sha256": _sha256(rendered),
                "required_measurements": [
                    "elapsed_seconds",
                    "completed_cycles",
                    "peak_resident_memory_bytes",
                    "artifact_payload_bytes",
                    "raw_payload_bytes",
                    "output_slot_count",
                ],
                "saturation_evidence_eligible": False,
                "production_resource_freeze_authorized": False,
            }
        )
    _require(
        len(artifacts) == EXPECTED_ARTIFACT_COUNT,
        "carrier calibration artifact count drifted",
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "status": "excluded_resource_calibration_packet_execution_prohibited",
        "controller_block": controller.CONTROLLER_BLOCK,
        "carrier_manifest": {
            "path": str(carrier.CHECKED_IN_MANIFEST.relative_to(REPO_ROOT)),
            "sha256": _sha256(carrier.CHECKED_IN_MANIFEST.read_bytes()),
            "record_type": carrier_manifest["record_type"],
        },
        "calibration_design": {
            "attempt_count": EXPECTED_ARTIFACT_COUNT,
            "tasks_per_node": TASKS_PER_NODE,
            "nodes_per_dimension": {
                str(dimension): nodes
                for dimension, nodes in CALIBRATION_NODES.items()
            },
            "short_runs_measure_startup_and_peak_memory": True,
            "long_baseline_instrumented_pairs_measure_diagnostic_overhead": True,
            "measured_results_required_before_production_resource_freeze": True,
        },
        "authorization": dict(CALIBRATION_AUTHORITY),
        "artifacts": artifacts,
    }


def materialize(root: Path = CHECKED_IN_ROOT) -> None:
    root.mkdir(parents=True, exist_ok=True)
    manifest = build_manifest()
    expected_names = {"deck_manifest.json"}
    for overlay in expected_overlays():
        filename = f"{overlay['artifact_id']}.athinput"
        expected_names.add(filename)
        (root / filename).write_text(render_overlay(overlay), encoding="utf-8")
    for path in root.iterdir():
        if path.is_file() and path.name not in expected_names:
            path.unlink()
    (root / "deck_manifest.json").write_bytes(_json_bytes(manifest))


def validate_checked_in_packet() -> dict[str, object]:
    expected = build_manifest()
    _require(CHECKED_IN_MANIFEST.is_file(), "carrier calibration manifest is absent")
    actual = json.loads(CHECKED_IN_MANIFEST.read_text(encoding="utf-8"))
    _require(actual == expected, "carrier calibration manifest drifted")
    expected_names = {
        "deck_manifest.json",
        *(f"{overlay['artifact_id']}.athinput" for overlay in expected_overlays()),
    }
    actual_names = {
        path.name for path in CHECKED_IN_ROOT.iterdir() if path.is_file()
    }
    _require(
        actual_names == expected_names,
        "carrier calibration deck inventory drifted",
    )
    for overlay in expected_overlays():
        path = CHECKED_IN_ROOT / f"{overlay['artifact_id']}.athinput"
        _require(
            path.read_text(encoding="utf-8") == render_overlay(overlay),
            f"{overlay['artifact_id']}: carrier calibration deck drifted",
        )
    return actual


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=CHECKED_IN_ROOT)
    parser.add_argument("--validate", action="store_true")
    arguments = parser.parse_args()
    if arguments.validate:
        _require(
            arguments.output_root == CHECKED_IN_ROOT,
            "checked-in validation does not accept another output root",
        )
        validate_checked_in_packet()
    else:
        materialize(arguments.output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
