#!/usr/bin/env python3
"""Focused tests for Q019 excluded-pilot launch and policy preparation."""

from __future__ import annotations

import copy
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

import pytest

from tst.publication import q019_excluded_pilot_launch_policy_preparation_v1 as prep


def test_pending_materialization_is_exact_paired_and_non_authorizing() -> None:
    manifest, files = prep.build_materialization()
    prep.validate_materialization(manifest, files)
    assert manifest["pilot_attempt_count"] == 4
    assert len(manifest["policy_slices"]) == 4
    assert len(files) == 10
    assert [item["pair_role"] for item in manifest["attempt_records"]] == [
        "baseline",
        "instrumented",
        "baseline",
        "instrumented",
    ]
    assert [item["pair_dimension"] for item in manifest["attempt_records"]] == [
        2,
        2,
        3,
        3,
    ]
    assert len({item["authorization_id"] for item in manifest["attempt_records"]}) == 4
    assert len({item["test_id"] for item in manifest["attempt_records"]}) == 2
    assert manifest["paired_pilot_contract"]["saturation_evidence_eligible"] is False
    assert all(value is False for value in prep.AUTHORIZATION_BOUNDARY.values())


def test_launch_contracts_bind_resource_classes_and_source_case_ids() -> None:
    manifest, files = prep.build_materialization()
    expected = {
        2: (4, 32, 1800, 64 * 1024**3),
        3: (16, 128, 3600, 256 * 1024**3),
    }
    for record in manifest["attempt_records"]:
        launch = json.loads(files[record["launch_candidate"]["path"]])
        action = launch["launch_contract"]["actions"][0]
        nodes, tasks, walltime, storage = expected[record["pair_dimension"]]
        assert action["resources"]["nodes"] == nodes
        assert action["resources"]["tasks"] == tasks
        assert launch["resource_ceiling"]["maximum_walltime_seconds"] == walltime
        assert launch["resource_ceiling"]["maximum_storage_bytes"] == storage
        assert launch["test_id"] == record["test_id"]
        assert launch["artifact_id"] == record["artifact_id"]
        assert launch["stdout_wrapper_evidence"]["required_controller_trigger_reason"] == 1903
        assert launch["stdout_wrapper_evidence"]["required_controller_trigger_cycle"] == 20
        assert launch["stdout_wrapper_evidence"][
            "required_saturation_evidence_eligible"
        ] == "false"


def test_budget_is_exact_and_one_attempt_without_retry() -> None:
    manifest, files = prep.build_materialization()
    budget = json.loads(files["batch_budget_accounting_input.json"])
    assert budget["maximum_batch_node_hours"] == 36.0
    assert budget["maximum_batch_storage_bytes"] == 640 * 1024**3
    assert budget["maximum_attempts_per_pilot"] == 1
    assert budget["maximum_retries_per_pilot"] == 0
    altered = copy.deepcopy(manifest)
    altered["pilot_attempt_count"] = 3
    with pytest.raises(prep.PreparationError, match="manifest drifted"):
        prep.validate_materialization(altered, files)


def test_overlay_byte_or_identity_drift_fails_closed() -> None:
    original = prep.controller.render_overlay
    with patch.object(
        prep.controller,
        "render_overlay",
        side_effect=lambda item: original(item) + "# drift\n",
    ):
        with pytest.raises(
            prep.PreparationError, match="manifest drifted|deck drifted"
        ):
            prep.build_materialization()


def test_duplicate_policy_authorization_is_rejected_by_exact_validation() -> None:
    manifest, files = prep.build_materialization()
    altered = copy.deepcopy(manifest)
    altered["attempt_records"][1]["authorization_id"] = altered["attempt_records"][0][
        "authorization_id"
    ]
    with pytest.raises(prep.PreparationError, match="manifest drifted"):
        prep.validate_materialization(altered, files)


def test_final_binding_structure_requires_both_predecessor_matrices() -> None:
    version = "a" * 64
    candidate = prep.AUTHORIZED_ORION_ROOT / "clean_candidates/q019-final"
    controller = prep.AUTHORIZED_ORION_ROOT / "control_plane" / version
    project = prep.CANONICAL_PROJECT_HOME_ROOT / "control_plane" / version
    final = {
        key: "b" * 64
        for key in prep.FINAL_BINDING_KEYS
        if key.endswith("_sha256")
    }
    final.update(
        {
            "record_type": prep.FINAL_BINDING_RECORD_TYPE,
            "schema_version": 1,
            "source_commit": "c" * 40,
            "source_archive_path": str(candidate / "source.tar"),
            "clean_candidate_manifest_path": str(
                candidate / "clean_candidate_manifest.json"
            ),
            "executable_path": str(candidate / "athena"),
            "installed_control_plane_version": version,
            "orion_installed_control_plane_root": str(controller),
            "project_home_installed_control_plane_root": str(project),
            "environment_profile_path": str(
                controller / "frontier_pic_environment.sh"
            ),
            "job_script_path": str(controller / "frontier_job.sh"),
            "analysis_script_paths": [
                str(controller / "reconcile_q019_registered_execution.py")
            ],
            "analysis_script_sha256": ["b" * 64],
            "reconcile_q019_registered_execution_path": str(
                controller / "reconcile_q019_registered_execution.py"
            ),
            "q043_registered_matrix_path": str(
                prep.AUTHORIZED_ORION_ROOT / "analysis/q043/matrix.json"
            ),
            "q043_registered_matrix_record_type": prep.q043.MATRIX_RECORD_TYPE,
            "q023_registered_matrix_path": str(
                prep.AUTHORIZED_ORION_ROOT / "analysis/q023/matrix.json"
            ),
            "q023_registered_matrix_record_type": prep.q023.MATRIX_RECORD_TYPE,
        }
    )
    prep.validate_final_bindings(final)
    del final["q023_registered_matrix_path"]
    with pytest.raises(prep.PreparationError, match="keys drifted"):
        prep.validate_final_bindings(final)


def test_timeout_artifact_is_dimension_specific_and_fresh() -> None:
    final = dict(prep.PENDING_FINAL_BINDINGS)
    final["environment_profile_sha256"] = "a" * 64
    now = datetime(2026, 6, 8, 12, 0, tzinfo=timezone.utc)
    with patch.object(prep, "validate_final_bindings", return_value=final):
        value = prep.materialize_q019_timeout_margin(
            artifact_id="q019-controller-pilot-3d-instrumented",
            final_bindings=final,
            measured_utc="2026-06-08T11:59:00Z",
            expires_utc="2026-06-08T12:59:00Z",
        )
        assert value["athena_walltime_seconds"] == 3300
        assert value["scheduler_walltime_seconds"] == 3600
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "timeout.json"
            path.write_text(json.dumps(value), encoding="utf-8")
            path.chmod(0o444)
            assert (
                prep.validate_q019_timeout_margin_artifact(
                    path,
                    artifact_id="q019-controller-pilot-3d-instrumented",
                    final_bindings=final,
                    now=now,
                )
                == str(path)
            )
            with pytest.raises(prep.PreparationError, match="stale"):
                prep.validate_q019_timeout_margin_artifact(
                    path,
                    artifact_id="q019-controller-pilot-3d-instrumented",
                    final_bindings=final,
                    now=now + timedelta(hours=2),
                )


def test_json_output_is_exclusive_and_read_only() -> None:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "bindings.json"
        value = {"record_type": "test", "schema_version": 1}
        prep._write_json_exclusive(path, value)
        assert json.loads(path.read_text(encoding="utf-8")) == value
        assert path.stat().st_mode & 0o777 == 0o444
        with pytest.raises(FileExistsError):
            prep._write_json_exclusive(path, value)
