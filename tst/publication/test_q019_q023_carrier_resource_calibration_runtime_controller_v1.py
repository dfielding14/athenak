#!/usr/bin/env python3
"""Tests for the Q019 Q023-carrier six-run resource calibration packet."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from tst.publication import q019_nonlinear_bell_runtime_controller_v1 as controller
from tst.publication import q019_q023_carrier_nonlinear_bell_redesign_v1 as carrier
from tst.publication import (
    q019_q023_carrier_resource_calibration_runtime_controller_v1 as calibration,
)


EXPECTED = (
    (
        "q019-carrier-calibration-2d-onset-c0032-baseline",
        "q019-q023-carrier-s1-onset-s0",
        32,
        "baseline",
    ),
    (
        "q019-carrier-calibration-2d-onset-c0256-baseline",
        "q019-q023-carrier-s1-onset-s0",
        256,
        "baseline",
    ),
    (
        "q019-carrier-calibration-2d-onset-c0256-instrumented",
        "q019-q023-carrier-s1-onset-s0",
        256,
        "instrumented",
    ),
    (
        "q019-carrier-calibration-3d-fiducial-c0016-baseline",
        "q019-q023-carrier-s3-3d-fiducial-s0",
        16,
        "baseline",
    ),
    (
        "q019-carrier-calibration-3d-fiducial-c0064-baseline",
        "q019-q023-carrier-s3-3d-fiducial-s0",
        64,
        "baseline",
    ),
    (
        "q019-carrier-calibration-3d-fiducial-c0064-instrumented",
        "q019-q023-carrier-s3-3d-fiducial-s0",
        64,
        "instrumented",
    ),
)


def test_exact_six_run_calibration_inventory_and_no_authority() -> None:
    overlays = calibration.expected_overlays()
    assert tuple(
        (
            item["artifact_id"],
            item["source_case_id"],
            item["cycle_limit"],
            item["instrumentation"],
        )
        for item in overlays
    ) == EXPECTED
    assert all(
        item["nodes"] == calibration.CALIBRATION_NODES[item["dimension"]]
        and item["tasks"] == item["nodes"] * calibration.TASKS_PER_NODE
        and item["expected_stop_reason"] == 1903
        and item["authority"] == "excluded_pilot_only"
        for item in overlays
    )
    manifest = calibration.build_manifest()
    assert all(value is False for value in manifest["authorization"].values())
    assert manifest["status"].endswith("execution_prohibited")
    assert all(
        item["saturation_evidence_eligible"] is False
        and item["production_resource_freeze_authorized"] is False
        for item in manifest["artifacts"]
    )


def test_matched_long_pairs_change_only_diagnostic_flags() -> None:
    overlays = {
        str(item["artifact_id"]): item
        for item in calibration.expected_overlays()
    }
    for prefix in (
        "q019-carrier-calibration-2d-onset-c0256",
        "q019-carrier-calibration-3d-fiducial-c0064",
    ):
        baseline = dict(overlays[f"{prefix}-baseline"]["controller_parameters"])
        instrumented = dict(
            overlays[f"{prefix}-instrumented"]["controller_parameters"]
        )
        baseline.pop("controller_identity_fingerprint")
        instrumented.pop("controller_identity_fingerprint")
        differing = {
            name for name in baseline if baseline[name] != instrumented[name]
        }
        assert differing == {
            "box_edge_monitor_enabled",
            "diagnostic_failure_stop_armed",
            "resolution_monitor_enabled",
        }


def test_rendered_decks_preserve_exact_carrier_base_and_bind_controller() -> None:
    for overlay in calibration.expected_overlays():
        rendered = calibration.render_overlay(overlay)
        base, block = rendered.split(
            f"\n<{controller.CONTROLLER_BLOCK}>\n", maxsplit=1
        )
        source = (
            carrier.CHECKED_IN_DECK_ROOT
            / f"{overlay['source_case_id']}.athinput"
        ).read_text(encoding="utf-8")
        assert base.rstrip() == source.rstrip()
        assert (
            "controller_identity_fingerprint = "
            + overlay["controller_parameters"]["controller_identity_fingerprint"]
        ) in block


def test_checked_in_packet_is_exact(tmp_path: Path) -> None:
    calibration.materialize(tmp_path)
    expected_names = {
        path.name for path in calibration.CHECKED_IN_ROOT.iterdir() if path.is_file()
    }
    assert {path.name for path in tmp_path.iterdir() if path.is_file()} == expected_names
    for name in expected_names:
        assert (tmp_path / name).read_bytes() == (
            calibration.CHECKED_IN_ROOT / name
        ).read_bytes()
    assert calibration.validate_checked_in_packet() == json.loads(
        calibration.CHECKED_IN_MANIFEST.read_text(encoding="utf-8")
    )


def test_controller_or_source_identity_tampering_fails_closed() -> None:
    overlay = copy.deepcopy(calibration.expected_overlays()[0])
    overlay["controller_parameters"]["pilot_cycle_limit"] = "33"
    with pytest.raises(
        calibration.CalibrationError,
        match="controller identity fingerprint drifted",
    ):
        calibration.render_overlay(overlay)

    missing = copy.deepcopy(calibration.expected_overlays()[0])
    missing["source_case_id"] = "q019-q023-carrier-does-not-exist"
    with pytest.raises(calibration.CalibrationError, match="source deck is absent"):
        calibration.render_overlay(missing)
