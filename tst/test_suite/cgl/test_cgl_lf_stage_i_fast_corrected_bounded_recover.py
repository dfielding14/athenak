"""Tests for duplicate-safe R14/R15 zero-field recovery."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
RECOVERY = (
    REPOSITORY
    / "scripts/frontier/cgl_lf_stage_i_fast_corrected_bounded_recover.py"
)


def load_recovery():
    name = "cgl_lf_stage_i_fast_corrected_bounded_recover_for_tests"
    spec = importlib.util.spec_from_file_location(name, RECOVERY)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def recovery():
    return load_recovery()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def test_recovery_contract_is_narrow_and_physics_preserving(recovery):
    assert recovery.FAILURE_SIGNATURE == (
        "cannot inject non-zero dedt with a zero forcing field"
    )
    assert recovery.FAILURE_REASON == "turbulence_driver_zero_forcing_field_guard"
    assert recovery.RECOVERY_POLICY.endswith("zero_forcing_guard_v1")
    assert "FAILED" in recovery.FAILED_STATES


def test_failed_attempt_requires_exact_signature(
    recovery, tmp_path, monkeypatch
):
    segment = tmp_path / "fast_s006_t7_to_t10"
    manifest = {
        "case_id": "R14",
        "sequence": 6,
        "job_id": "123",
        "slurm_log": str(tmp_path / "%x.%j.log"),
    }
    write_json(segment / "manifest/fast_run.json", manifest)
    (segment / "manifest/run_exit_code").write_text("143\n", encoding="utf-8")
    log = tmp_path / "cglc_R14_s006.123.log"
    log.write_text("different failure\n", encoding="utf-8")
    monkeypatch.setattr(
        recovery.bounded, "validate_bounded_segment", lambda _path: manifest
    )
    monkeypatch.setattr(
        recovery.bounded.corrected.fast, "job_state", lambda _job: "FAILED"
    )

    with pytest.raises(recovery.BoundedRecoveryError, match="zero-field signature"):
        recovery.failed_attempt_record(segment)


def test_failed_attempt_record_binds_manifest_log_and_exit(
    recovery, tmp_path, monkeypatch
):
    segment = tmp_path / "fast_s006_t7_to_t10"
    manifest = {
        "case_id": "R14",
        "sequence": 6,
        "job_id": "123",
        "slurm_log": str(tmp_path / "%x.%j.log"),
    }
    write_json(segment / "manifest/fast_run.json", manifest)
    (segment / "manifest/run_exit_code").write_text("143\n", encoding="utf-8")
    log = tmp_path / "cglc_R14_s006.123.log"
    log.write_text(recovery.FAILURE_SIGNATURE + "\n", encoding="utf-8")
    monkeypatch.setattr(
        recovery.bounded, "validate_bounded_segment", lambda _path: manifest
    )
    monkeypatch.setattr(
        recovery.bounded.corrected.fast, "job_state", lambda _job: "FAILED"
    )

    record = recovery.failed_attempt_record(segment)

    assert record["sequence"] == 6
    assert record["job_id"] == "123"
    assert record["exit_code"] == 143
    assert record["failure_signature"] == recovery.FAILURE_SIGNATURE
    assert len(str(record["manifest_sha256"])) == 64
    assert len(str(record["slurm_log_sha256"])) == 64


def test_prepared_recovery_requires_failed_attempt_evidence(
    recovery, tmp_path, monkeypatch
):
    segment = tmp_path / "fast_s007_t7_to_t10"
    manifest = {
        "recovery_policy": recovery.RECOVERY_POLICY,
        "recovery_reason": recovery.FAILURE_REASON,
        "recovery_changes_physics": False,
        "supersedes_failed_segments": [],
    }
    monkeypatch.setattr(
        recovery.bounded, "validate_bounded_segment", lambda _path: manifest
    )
    with pytest.raises(recovery.BoundedRecoveryError, match="failed-attempt evidence"):
        recovery.validate_prepared_recovery(segment)


def test_active_case_blocks_recovery(recovery, monkeypatch):
    class Completed:
        stdout = "cglc_R14_s007\n"

    monkeypatch.setattr(
        recovery.subprocess,
        "check_output",
        lambda *_args, **_kwargs: "user\n",
    )
    monkeypatch.setattr(
        recovery.subprocess, "run", lambda *_args, **_kwargs: Completed()
    )
    with pytest.raises(recovery.BoundedRecoveryError, match="already has an active"):
        recovery.require_inactive_case("R14")
