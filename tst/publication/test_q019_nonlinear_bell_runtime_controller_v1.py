#!/usr/bin/env python3
"""Tests for the sealed Q019 nonlinear Bell runtime-controller packet."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess

import pytest

from tst.publication import q019_nonlinear_bell_runtime_controller_v1 as controller
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as q019


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE = REPO_ROOT / "src/pgen/tests/q019_physics_first_nonlinear_bell_successor_v2.cpp"


def test_checked_packet_is_exact_and_non_authorizing(tmp_path: Path) -> None:
    controller.materialize(tmp_path)
    expected_files = {
        path.name for path in controller.CHECKED_IN_ROOT.iterdir() if path.is_file()
    }
    assert {path.name for path in tmp_path.iterdir() if path.is_file()} == expected_files
    for name in expected_files:
        assert (tmp_path / name).read_bytes() == (
            controller.CHECKED_IN_ROOT / name
        ).read_bytes()

    manifest = json.loads(controller.CHECKED_IN_MANIFEST.read_text())
    assert manifest == controller.build_manifest()
    assert manifest["saturation_evidence_eligible"] is False
    assert manifest["production_promotion_authorized"] is False
    assert len(manifest["artifacts"]) == 25
    assert all(
        not artifact["saturation_evidence_eligible"]
        and not artifact["production_promotion_authorized"]
        for artifact in manifest["artifacts"]
    )


def test_overlay_identity_and_base_matrix_are_independently_bound() -> None:
    for overlay in controller.expected_overlays():
        parameters = dict(overlay["controller_parameters"])
        assert parameters["controller_identity_fingerprint"] == (
            controller.controller_identity_fingerprint(parameters)
        )
        rendered = controller.render_overlay(overlay)
        base, block = rendered.split(
            f"\n<{controller.CONTROLLER_BLOCK}>\n", maxsplit=1
        )
        source_case_id = str(overlay["source_case_id"])
        checked_base = (
            q019.CHECKED_IN_DECK_ROOT / f"{source_case_id}.athinput"
        ).read_text()
        assert base.rstrip() == checked_base.rstrip()
        assert "controller_identity_fingerprint =" in block

        mutated = dict(parameters)
        mutated["pilot_cycle_limit"] = str(int(mutated["pilot_cycle_limit"]) + 1)
        assert mutated["controller_identity_fingerprint"] != (
            controller.controller_identity_fingerprint(mutated)
        )


def test_matched_pilot_pairs_change_only_diagnostic_overlay() -> None:
    overlays = {
        str(overlay["artifact_id"]): overlay
        for overlay in controller.expected_overlays()
    }
    for dimension in ("2d", "3d"):
        baseline = overlays[f"q019-controller-pilot-{dimension}-baseline"]
        instrumented = overlays[f"q019-controller-pilot-{dimension}-instrumented"]
        assert baseline["source_case_id"] == instrumented["source_case_id"]
        left = dict(baseline["controller_parameters"])
        right = dict(instrumented["controller_parameters"])
        left.pop("controller_identity_fingerprint")
        right.pop("controller_identity_fingerprint")
        differing = {name for name in left if left[name] != right[name]}
        assert differing == {
            "box_edge_monitor_enabled",
            "diagnostic_failure_stop_armed",
            "resolution_monitor_enabled",
        }


def test_physical_pilot_overlays_are_exact_monitor_only_stage_inventory() -> None:
    physical = [
        overlay
        for overlay in controller.expected_overlays()
        if str(overlay["artifact_id"]).startswith("q019-physical-pilot-")
    ]
    assert len(physical) == 19
    assert [overlay["physical_pilot_stage"] for overlay in physical] == [
        *([1] * 7),
        *([2] * 8),
        *([3] * 4),
    ]
    assert all(
        overlay["monitor_dt"] == q019.BOX_EDGE_MONITOR_DT
        and overlay["box_edge_monitor_enabled"] is True
        and overlay["resolution_monitor_enabled"] is True
        and overlay["diagnostic_failure_stop_armed"] is True
        and overlay["resolution_stop_armed"] is False
        and overlay["box_edge_stop_armed"] is False
        and overlay["pilot_cycle_limit"] == -1
        and overlay["expected_stop_reason"] is None
        and overlay["saturation_evidence_eligible"] is False
        and overlay["production_promotion_authorized"] is False
        for overlay in controller.build_manifest()["artifacts"]
        if str(overlay["artifact_id"]).startswith("q019-physical-pilot-")
    )


def test_source_contract_is_fail_closed_and_restart_persistent() -> None:
    source = SOURCE.read_text()
    assert "Q019RuntimeControllerIdentityPayload" in source
    assert "controller_identity_fingerprint" in source
    assert "Q019StoreRuntimeControllerState();" in source
    assert "RequestUserStop(reason, failure)" in source
    assert "kQ019ResolutionStopReason = 1901" in source
    assert "kQ019BoxEdgeStopReason = 1902" in source
    assert "kQ019PilotCycleStopReason = 1903" in source
    assert "kQ019DiagnosticFailureReason = 1991" in source
    assert "Q019_SATURATION_EVIDENCE_ELIGIBLE=false" in source


def _binary() -> Path:
    executable_dir = os.environ.get("ATHENA_Q019_EXE_DIR")
    if not executable_dir:
        pytest.skip("ATHENA_Q019_EXE_DIR is required")
    binary = Path(executable_dir) / "athena"
    if not binary.is_file():
        raise RuntimeError(f"Q019 runtime executable is missing: {binary}")
    return binary


def _run(
    tmp_path: Path,
    artifact_id: str,
    *,
    restart: Path | None = None,
    deck: Path | None = None,
) -> tuple[subprocess.CompletedProcess[str], Path]:
    run_dir = tmp_path / f"{artifact_id}-{len(list(tmp_path.iterdir()))}"
    run_dir.mkdir()
    command = [str(_binary())]
    if restart is None:
        selected = deck or controller.CHECKED_IN_ROOT / f"{artifact_id}.athinput"
        command.extend(["-i", str(selected)])
    else:
        command.extend(["-r", str(restart)])
    result = subprocess.run(
        command,
        cwd=run_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=180,
        check=False,
    )
    return result, run_dir


def test_runtime_box_stop_and_stopped_restart_are_idempotent(tmp_path: Path) -> None:
    artifact = "q019-controller-regression-box-stop"
    first, run_dir = _run(tmp_path, artifact)
    assert first.returncode == 0, first.stdout
    assert "user_stop_reason_code=1902 user_stop_failure=false" in first.stdout
    assert "Q019_RUNTIME_CONTROLLER_TRIGGER_REASON=1902" in first.stdout
    assert "Q019_RUNTIME_RESOLUTION_SAMPLES=1" in first.stdout
    restarts = sorted((run_dir / "rst").glob("*.rst"))
    assert len(restarts) >= 2
    assert all(Path(f"{path}.complete").is_file() for path in restarts)

    replay, replay_dir = _run(tmp_path, artifact, restart=restarts[-1])
    assert replay.returncode == 0, replay.stdout
    assert "user_stop_reason_code=1902 user_stop_failure=false" in replay.stdout
    assert "Q019_RUNTIME_CONTROLLER_TRIGGER_REASON=1902" in replay.stdout
    assert len(list((replay_dir / "rst").glob("*.rst"))) == 1


def test_runtime_cycle_cap_and_controller_hash_rejection(tmp_path: Path) -> None:
    artifact = "q019-controller-regression-cycle-stop"
    result, run_dir = _run(tmp_path, artifact)
    assert result.returncode == 0, result.stdout
    assert "user_stop_reason_code=1903 user_stop_failure=false" in result.stdout
    assert "Q019_RUNTIME_RESOLUTION_SAMPLES=0" in result.stdout
    assert "Q019_BOX_EDGE_MONITOR_COMPLETED_SLOTS=0" in result.stdout
    restarts = sorted((run_dir / "rst").glob("*.rst"))
    assert len(restarts) >= 2
    replay, _ = _run(tmp_path, artifact, restart=restarts[-1])
    assert replay.returncode == 0, replay.stdout
    assert "user_stop_reason_code=1903 user_stop_failure=false" in replay.stdout
    assert "Q019_RUNTIME_RESOLUTION_SAMPLES=0" in replay.stdout

    original = controller.CHECKED_IN_ROOT / f"{artifact}.athinput"
    hostile = tmp_path / "hostile.athinput"
    hostile.write_text(
        original.read_text().replace("pilot_cycle_limit = 1", "pilot_cycle_limit = 2", 1)
    )
    rejected, _ = _run(tmp_path, "hostile", deck=hostile)
    assert rejected.returncode != 0
    assert "runtime controller identity fingerprint drifted" in rejected.stdout


@pytest.mark.parametrize("parameter", ("unbound_parameter", "runtime_unbound"))
def test_runtime_controller_rejects_unknown_parameters(
    tmp_path: Path, parameter: str
) -> None:
    artifact = "q019-controller-regression-cycle-stop"
    original = controller.CHECKED_IN_ROOT / f"{artifact}.athinput"
    hostile = tmp_path / f"{parameter}.athinput"
    hostile.write_text(original.read_text() + f"{parameter} = 0\n")
    rejected, _ = _run(tmp_path, parameter, deck=hostile)
    assert rejected.returncode != 0
    assert "runtime controller contains an unknown parameter" in rejected.stdout


def test_runtime_controller_rejects_partial_mutable_inventory(
    tmp_path: Path,
) -> None:
    artifact = "q019-controller-regression-cycle-stop"
    original = controller.CHECKED_IN_ROOT / f"{artifact}.athinput"
    hostile = tmp_path / "partial-runtime-state.athinput"
    hostile.write_text(original.read_text() + "runtime_resolution_samples = 0\n")
    rejected, _ = _run(tmp_path, "partial-runtime-state", deck=hostile)
    assert rejected.returncode != 0
    assert "runtime controller mutable inventory drifted" in rejected.stdout


def test_runtime_controller_rejects_inconsistent_complete_trigger_state(
    tmp_path: Path,
) -> None:
    artifact = "q019-controller-regression-cycle-stop"
    original = controller.CHECKED_IN_ROOT / f"{artifact}.athinput"
    runtime_state = {
        "runtime_resolution_samples": "0",
        "runtime_resolution_last_cycle": "0",
        "runtime_resolution_last_time": "0",
        "runtime_resolution_last_B_over_B0": "-1",
        "runtime_resolution_max_B_over_B0": "-1",
        "runtime_controller_triggered": "false",
        "runtime_controller_trigger_failure": "false",
        "runtime_controller_trigger_reason": "1903",
        "runtime_controller_trigger_cycle": "0",
        "runtime_controller_trigger_time": "0",
        "runtime_controller_trigger_metric": "-1",
    }
    hostile = tmp_path / "inconsistent-runtime-state.athinput"
    hostile.write_text(
        original.read_text()
        + "".join(f"{name} = {value}\n" for name, value in runtime_state.items())
    )
    rejected, _ = _run(tmp_path, "inconsistent-runtime-state", deck=hostile)
    assert rejected.returncode != 0
    assert "runtime controller restart state drifted" in rejected.stdout
