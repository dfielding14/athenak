from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import struct
import tarfile

import numpy as np
import pytest

from tst.publication import (
    q019_hardened_installed_control_plane_registered_admission_v1 as admission,
)


def _marker(payload: bytes) -> bytes:
    return (
        "ATHENAK_RESTART_COMPLETE_V1\n"
        f"size={len(payload)}\n"
        f"fnv1a64={admission._fnv1a64(payload):016x}\n"
    ).encode("ascii")


def _restart(cycle: int = 7, time: float = 0.5) -> bytes:
    mesh_header = struct.pack(
        admission._RESTART_MESH_HEADER_FORMAT,
        2,
        1,
        *([0.0] * 9),
        *([0] * 19),
        *([0] * 19),
        time,
        0.01,
        cycle,
        1,
    )
    return b"<time>\ntlim = 1\n<par_end>\n" + mesh_header + b"schema-7-payload"


def _publication() -> tuple[str, dict[str, bytes]]:
    stem = "rst/q019-fr-runtime-initializer-ppc24-s0.00001.rst"
    restart = _restart()
    manifest = (
        json.dumps(
            {
                "schema": "ATHENAK_RESTART_MANIFEST_V1",
                "members": [
                    {
                        "path": stem,
                        "size": len(restart),
                        "fnv1a64": f"{admission._fnv1a64(restart):016x}",
                    }
                ],
            },
            indent=2,
        )
        + "\n"
    ).encode("utf-8")
    return stem, {
        stem: restart,
        stem + ".complete": _marker(restart),
        stem + ".manifest": manifest,
        stem + ".manifest.complete": _marker(manifest),
    }


def _source_archive(*, omit: str | None = None) -> bytes:
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w") as archive:
        for path in sorted(admission.REQUIRED_SOURCE_PATHS):
            if path == omit:
                continue
            payload = f"fixture:{path}\n".encode("ascii")
            member = tarfile.TarInfo(path)
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
    return output.getvalue()


def test_restart_publication_binds_marker_manifest_cycle_and_time() -> None:
    stem, payloads = _publication()
    assert admission._validate_restart_publication(stem=stem, payloads=payloads) == (
        7,
        0.5,
    )
    payloads[stem + ".complete"] = payloads[stem + ".complete"].replace(
        b"size=", b"size=1"
    )
    with pytest.raises(admission.RegisteredAdmissionError, match="marker"):
        admission._validate_restart_publication(stem=stem, payloads=payloads)


def test_restart_manifest_substitution_fails_closed() -> None:
    stem, payloads = _publication()
    manifest = json.loads(payloads[stem + ".manifest"])
    manifest["members"][0]["path"] = "rst/substituted.rst"
    replacement = (json.dumps(manifest, indent=2) + "\n").encode("utf-8")
    payloads[stem + ".manifest"] = replacement
    payloads[stem + ".manifest.complete"] = _marker(replacement)
    with pytest.raises(admission.RegisteredAdmissionError, match="manifest"):
        admission._validate_restart_publication(stem=stem, payloads=payloads)


def test_source_archive_requires_complete_runtime_and_analysis_closure() -> None:
    bindings = admission._source_archive_bindings(_source_archive())
    assert set(bindings) == admission.REQUIRED_SOURCE_PATHS
    assert all(
        len(item["sha256"]) == 64 and item["byte_count"] > 0
        for item in bindings.values()
    )
    missing = next(iter(admission.REQUIRED_SOURCE_PATHS))
    with pytest.raises(admission.RegisteredAdmissionError, match="complete"):
        admission._source_archive_bindings(_source_archive(omit=missing))


def test_candidate_analysis_closure_binds_manifest_deck_and_executing_sources() -> None:
    case = admission._case("q019-fr-runtime-initializer-ppc24-s0")
    required = {
        *admission.REQUIRED_SOURCE_PATHS,
        *admission.EXECUTING_QUALIFICATION_SOURCE_PATHS,
        *admission.design.ANALYSIS_BINDING_PATHS,
        str(case["path"]),
    }
    payloads = {
        relative: (admission.REPO_ROOT / relative).read_bytes()
        for relative in required
    }
    closure = admission._candidate_analysis_closure(payloads, case=case)
    assert closure["registered_deck"]["sha256"] == case["sha256"]
    assert [item["path"] for item in closure["analysis_bindings"]] == list(
        admission.design.ANALYSIS_BINDING_PATHS
    )
    assert set(closure["executing_qualification_source_bindings"]) == (
        admission.EXECUTING_QUALIFICATION_SOURCE_PATHS
    )
    assert closure["execution_deck"]["kind"] == "base_matrix_case"

    substituted = dict(payloads)
    analyzer = "tst/publication/analyze_q019_physics_first_nonlinear_bell_successor_v2.py"
    substituted[analyzer] += b"\n# substituted\n"
    with pytest.raises(admission.RegisteredAdmissionError, match="analysis binding"):
        admission._candidate_analysis_closure(substituted, case=case)


def test_candidate_analysis_closure_accepts_only_exact_controller_overlay() -> None:
    case = admission._case("q019-fr-grid-k8-rho1em05-s0")
    required = {
        *admission.REQUIRED_SOURCE_PATHS,
        *admission.EXECUTING_QUALIFICATION_SOURCE_PATHS,
        *admission.design.ANALYSIS_BINDING_PATHS,
        str(case["path"]),
    }
    payloads = {
        relative: (admission.REPO_ROOT / relative).read_bytes()
        for relative in required
    }
    overlay = next(
        item
        for item in admission.controller.build_manifest()["artifacts"]
        if item["artifact_id"] == "q019-controller-pilot-2d-instrumented"
    )
    closure = admission._candidate_analysis_closure(
        payloads,
        case=case,
        execution_deck_sha256=str(overlay["rendered_sha256"]),
    )
    assert closure["execution_deck"]["artifact_id"] == overlay["artifact_id"]
    assert closure["execution_deck"]["authority"] == "excluded_pilot_only"
    assert closure["execution_deck"]["saturation_evidence_eligible"] is False

    with pytest.raises(
        admission.RegisteredAdmissionError,
        match="not one exact runtime-controller overlay",
    ):
        admission._candidate_analysis_closure(
            payloads,
            case=case,
            execution_deck_sha256="0" * 64,
        )

    substituted = dict(payloads)
    dependency = "tst/publication/q023_registered_execution_linear_qualification_successor_v1.py"
    substituted[dependency] += b"\n# substituted\n"
    with pytest.raises(
        admission.RegisteredAdmissionError,
        match="executing qualification source differs",
    ):
        admission._candidate_analysis_closure(substituted, case=case)


def test_reduction_binding_hashes_exact_array_bytes() -> None:
    reduction = {
        "record_type": "q019_registered_raw_reduction_v1",
        "case_id": "q019-fr-runtime-initializer-ppc24-s0",
        "campaign_id": "Q019",
        "matched_checkpoint_count": 1,
        "chronology": [{"cycle": 0, "time": 0.0}],
        "reference_budget": {"total_momentum": [0.0, 0.0, 0.0], "total_energy": 1.0},
        "execution_profile": {
            "kind": "base_matrix_case",
            "source_case_id": "q019-fr-runtime-initializer-ppc24-s0",
            "artifact_id": None,
            "authority": "matrix_case",
            "saturation_evidence_eligible": False,
        },
        "runtime_controller_states": [None],
        "snapshots": [
            {
                "cycle": 0,
                "time": 0.0,
                "x1_faces": np.asarray([0.0, 1.0]),
                "x2_faces": np.asarray([0.0, 1.0]),
                "x3_faces": np.asarray([0.0, 1.0]),
                "fields": {"dens": np.ones((1, 1, 1), dtype=np.float64)},
            }
        ],
        "particle_states": [{"cycle": 0, "authority": {"claim_authorized": False}}],
    }
    binding = admission._reduction_binding(reduction)
    expected = hashlib.sha256(np.ones((1, 1, 1)).tobytes()).hexdigest()
    assert binding["snapshots"][0]["fields"]["dens"]["sha256"] == expected
    changed = dict(reduction)
    changed["snapshots"] = [dict(reduction["snapshots"][0])]
    changed["snapshots"][0]["fields"] = {"dens": np.zeros((1, 1, 1))}
    assert (
        admission._reduction_binding(changed)["snapshots"][0]["fields"]["dens"][
            "sha256"
        ]
        != expected
    )


def test_completion_status_is_bound_to_trusted_wrapper_evidence() -> None:
    record = admission._completion_record(
        {
            "command_evidence": {
                "trusted_wrapper_evidence": {
                    "termination_reason": "Terminating on time limit",
                    "problem_final_evidence_status": "completed_saturation_eligible",
                    "problem_saturation_evidence_eligible": "true",
                }
            },
            "slurm_terminal_state": "COMPLETED",
        }
    )
    assert record["trusted_execution_binding_present"]
    assert record["problem_saturation_evidence_eligible"]
    assert record["stop_reason_code"] == "Terminating on time limit"
    assert record["runtime_controller_trigger_cycle"] is None

    controller_record = admission._completion_record(
        {
            "command_evidence": {
                "trusted_wrapper_evidence": {
                    "termination_reason": "Terminating on user request",
                    "problem_final_evidence_status": (
                        "completed_not_acceptance_eligible"
                    ),
                    "problem_saturation_evidence_eligible": "false",
                }
            },
            "slurm_terminal_state": "COMPLETED",
        },
        execution_profile={
            "kind": "runtime_controller_overlay",
            "expected_stop_reason": 1903,
        },
        runtime_controller_states=[
            {
                "runtime_controller_triggered": False,
                "runtime_controller_trigger_failure": False,
                "runtime_controller_trigger_reason": 0,
                "runtime_controller_trigger_cycle": 0,
                "runtime_controller_trigger_time": 0.0,
                "runtime_controller_trigger_metric": -1.0,
            },
            {
                "runtime_controller_triggered": True,
                "runtime_controller_trigger_failure": False,
                "runtime_controller_trigger_reason": 1903,
                "runtime_controller_trigger_cycle": 20,
                "runtime_controller_trigger_time": 0.01,
                "runtime_controller_trigger_metric": 20.0,
            },
        ],
    )
    assert controller_record["stop_reason_code"] == "1903"
    assert controller_record["runtime_controller_trigger_cycle"] == 20
    assert controller_record["runtime_controller_trigger_time"] == 0.01

    monitor_only = admission._completion_record(
        {
            "command_evidence": {
                "trusted_wrapper_evidence": {
                    "termination_reason": "Terminating on time limit",
                    "problem_final_evidence_status": (
                        "completed_not_acceptance_eligible"
                    ),
                    "problem_saturation_evidence_eligible": "false",
                }
            },
            "slurm_terminal_state": "COMPLETED",
        },
        execution_profile={
            "kind": "runtime_controller_overlay",
            "expected_stop_reason": None,
        },
        runtime_controller_states=[
            {
                "runtime_controller_triggered": False,
                "runtime_controller_trigger_failure": False,
                "runtime_controller_trigger_reason": 0,
                "runtime_controller_trigger_cycle": 0,
                "runtime_controller_trigger_time": 0.0,
                "runtime_controller_trigger_metric": -1.0,
            }
        ],
    )
    assert monitor_only["stop_reason_code"] == "Terminating on time limit"
    assert monitor_only["runtime_controller_trigger_cycle"] is None

    with pytest.raises(
        admission.RegisteredAdmissionError,
        match="completion status or state is invalid",
    ):
        admission._completion_record(
            {
                "command_evidence": {
                    "trusted_wrapper_evidence": {
                        "termination_reason": "Terminating on time limit",
                        "problem_final_evidence_status": "completed_saturation_eligible",
                        "problem_saturation_evidence_eligible": "true",
                    }
                },
                "slurm_terminal_state": "COMPLETED",
            },
            execution_profile={
                "kind": "runtime_controller_overlay",
                "expected_stop_reason": None,
            },
            runtime_controller_states=[
                {
                    "runtime_controller_triggered": False,
                    "runtime_controller_trigger_failure": False,
                    "runtime_controller_trigger_reason": 0,
                    "runtime_controller_trigger_cycle": 0,
                    "runtime_controller_trigger_time": 0.0,
                    "runtime_controller_trigger_metric": -1.0,
                }
            ],
        )

    with pytest.raises(
        admission.RegisteredAdmissionError,
        match="monitor-only controller completion state drifted",
    ):
        admission._completion_record(
            {
                "command_evidence": {
                    "trusted_wrapper_evidence": {
                        "termination_reason": "Terminating on cycle limit",
                        "problem_final_evidence_status": (
                            "completed_not_acceptance_eligible"
                        ),
                        "problem_saturation_evidence_eligible": "false",
                    }
                },
                "slurm_terminal_state": "COMPLETED",
            },
            execution_profile={
                "kind": "runtime_controller_overlay",
                "expected_stop_reason": None,
            },
            runtime_controller_states=[
                {
                    "runtime_controller_triggered": False,
                    "runtime_controller_trigger_failure": False,
                    "runtime_controller_trigger_reason": 0,
                    "runtime_controller_trigger_cycle": 0,
                    "runtime_controller_trigger_time": 0.0,
                    "runtime_controller_trigger_metric": -1.0,
                }
            ],
        )

    with pytest.raises(
        admission.RegisteredAdmissionError,
        match="completion reason drifted",
    ):
        admission._completion_record(
            {
                "command_evidence": {
                    "trusted_wrapper_evidence": {
                        "termination_reason": "Terminating on user request",
                        "problem_final_evidence_status": (
                            "completed_not_acceptance_eligible"
                        ),
                        "problem_saturation_evidence_eligible": "false",
                    }
                },
                "slurm_terminal_state": "COMPLETED",
            },
            execution_profile={
                "kind": "runtime_controller_overlay",
                "expected_stop_reason": 1903,
            },
            runtime_controller_states=[
                {
                    "runtime_controller_triggered": True,
                    "runtime_controller_trigger_failure": False,
                    "runtime_controller_trigger_reason": 1902,
                    "runtime_controller_trigger_cycle": 20,
                    "runtime_controller_trigger_time": 0.01,
                    "runtime_controller_trigger_metric": 0.5,
                }
            ],
        )
    with pytest.raises(admission.RegisteredAdmissionError, match="status"):
        admission._completion_record(
            {
                "command_evidence": {
                    "trusted_wrapper_evidence": {
                        "termination_reason": "Terminating on time limit",
                        "problem_final_evidence_status": "completed",
                        "problem_saturation_evidence_eligible": True,
                    }
                }
            }
        )


def test_caller_supplied_admission_cannot_self_attest() -> None:
    with pytest.raises(admission.RegisteredAdmissionError, match="path anchor"):
        admission.validate_admission(
            {
                "record_type": admission.RECORD_TYPE,
                "case_id": "q019-fr-runtime-initializer-ppc24-s0",
                "artifact_root": "/tmp/self-attested",
                "q043_dependency": {"path": "/tmp/q043.json"},
                "q043_artifact_root": "/tmp/q043",
                "q023_dependency": {"path": "/tmp/q023.json"},
                "raw_science_admission_eligible": True,
            }
        )


def test_bound_analysis_inputs_reject_array_or_completion_substitution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    record = {"record_type": admission.RECORD_TYPE}
    bundle = {
        "snapshots": [{"fields": {"dens": np.ones((1, 1, 1), dtype=np.float64)}}],
        "particle_states": [{"cycle": 0, "time": 0.0}],
        "completion_record": {"record_type": "q019_runtime_completion_status_v1"},
    }
    monkeypatch.setattr(
        admission, "validate_analysis_bundle", lambda value: (record, bundle)
    )
    assert (
        admission.validate_bound_analysis_inputs(
            record,
            snapshots=bundle["snapshots"],
            particle_states=bundle["particle_states"],
            completion_record=bundle["completion_record"],
        )
        == record
    )
    substituted = [{"fields": {"dens": np.zeros((1, 1, 1), dtype=np.float64)}}]
    with pytest.raises(admission.RegisteredAdmissionError, match="differ"):
        admission.validate_bound_analysis_inputs(
            record,
            snapshots=substituted,
            particle_states=bundle["particle_states"],
            completion_record=bundle["completion_record"],
        )


def test_derive_case_bundle_cross_binds_q023_to_selected_q043(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    orion = tmp_path / "orion"
    project_home = tmp_path / "project-home"
    q043_root = orion / "analysis/q043"
    artifact_root = (
        orion
        / "runs"
        / admission.REGISTERED_CAMPAIGN
        / "12345678-1234-4123-8123-123456789abc"
    )
    for path in (q043_root, artifact_root, project_home):
        path.mkdir(parents=True, exist_ok=True)
    q043_record = {
        "record_type": "q043_registered_execution_raw_oracle_matrix_qualification",
        "case_bindings_sha256": "a" * 64,
    }
    q043_path = q043_root / "matrix.json"
    q043_payload = (json.dumps(q043_record) + "\n").encode("utf-8")
    q043_path.write_bytes(q043_payload)
    q043_path.chmod(0o444)
    q043_sha256 = hashlib.sha256(q043_payload).hexdigest()
    q023_path = orion / "analysis/q023.json"
    q023_path.parent.mkdir(parents=True, exist_ok=True)
    q023_path.write_text("{}\n", encoding="utf-8")
    q023_path.chmod(0o444)
    q023_record = {
        "q043_dependency": {
            "binding_kind": "registered_matrix_qualification",
            "registered_matrix_path": "matrix.json",
            "registered_matrix_sha256": q043_sha256,
            "registered_matrix_record_type": q043_record["record_type"],
            "registered_matrix_case_bindings_sha256": q043_record[
                "case_bindings_sha256"
            ],
        }
    }
    receipt = {
        "member_id": "q019-fr-runtime-initializer-ppc24-s0",
        "campaign_id": "Q019-PHYSICS-FIRST-NONLINEAR-BELL",
        "slurm_terminal_state": "COMPLETED",
        "slurm_exit_code": "0:0",
        "reservation_id": "reservation",
        "submission_id": artifact_root.name,
        "reconciliation_event_sha256": "b" * 64,
        "reconciliation_mirror_ack_sha256": "c" * 64,
        "control_plane_version": "d" * 64,
        "registered_science_authorization_id": "q019-fixture",
        "slurm_job_id": "12345",
        "pre_submit_manifest_sha256": "e" * 64,
        "clean_candidate_manifest_sha256": "1" * 64,
        "source_commit": "2" * 40,
        "source_bundle_sha256": "3" * 64,
        "source_archive_sha256": "4" * 64,
        "executable_sha256": "5" * 64,
        "environment_sha256": "6" * 64,
        "deck_sha256": "7" * 64,
        "raw_inventory_sha256": "8" * 64,
        "producer": {
            "entrypoint": admission.PRODUCER_ENTRYPOINT,
            "entrypoint_sha256": "9" * 64,
            "launch_trampoline_sha256": "a" * 64,
            "control_plane_version": "d" * 64,
        },
        "command_evidence": {
            "trusted_wrapper_evidence": {
                "termination_reason": "Terminating on cycle limit",
                "problem_final_evidence_status": "completed_not_acceptance_eligible",
                "problem_saturation_evidence_eligible": "true",
            }
        },
    }
    receipt_payload = b"receipt"
    terminal_payload = b"terminal"
    paired = {
        "execution_receipt": {
            "sha256": hashlib.sha256(receipt_payload).hexdigest()
        },
        "terminal_receipt": {
            "sha256": hashlib.sha256(terminal_payload).hexdigest()
        },
        "_receipt_payload": receipt_payload,
        "_terminal_payload": terminal_payload,
    }
    reduction = {
        "snapshots": [{"cycle": 0, "time": 0.0}],
        "particle_states": [{"cycle": 0, "time": 0.0}],
        "execution_profile": {
            "kind": "base_matrix_case",
            "source_case_id": receipt["member_id"],
            "artifact_id": None,
            "authority": "matrix_case",
            "saturation_evidence_eligible": False,
        },
    }
    monkeypatch.setattr(
        admission.q043,
        "validate_downstream_q023_q019_prerequisite",
        lambda value: q043_record,
    )
    monkeypatch.setattr(
        admission.q023,
        "validate_downstream_q019_prerequisite",
        lambda *args, **kwargs: q023_record,
    )
    monkeypatch.setattr(admission, "_receipt_pair", lambda **kwargs: (receipt, paired))
    monkeypatch.setattr(
        admission,
        "_installed_rederivation",
        lambda *args, **kwargs: {
            "control_plane_version": receipt["control_plane_version"],
            "entrypoint": admission.PRODUCER_ENTRYPOINT,
            "entrypoint_sha256": receipt["producer"]["entrypoint_sha256"],
            "launch_trampoline_sha256": receipt["producer"][
                "launch_trampoline_sha256"
            ],
            "reconciliation_event_sha256": receipt[
                "reconciliation_event_sha256"
            ],
            "reconciliation_mirror_ack_sha256": receipt[
                "reconciliation_mirror_ack_sha256"
            ],
            "exact_byte_rederivation_passed": True,
        },
    )
    monkeypatch.setattr(
        admission,
        "_manifest_and_candidate",
        lambda *args, **kwargs: {
            "source_bindings": {},
            "pre_submit_manifest": {
                "sha256": receipt["pre_submit_manifest_sha256"]
            },
            "clean_candidate_manifest": {
                "sha256": receipt["clean_candidate_manifest_sha256"]
            },
            "source_archive": {"sha256": receipt["source_archive_sha256"]},
            "executable_snapshot": {"sha256": receipt["executable_sha256"]},
            "analysis_closure": {
                "execution_deck": {
                    "kind": "base_matrix_case",
                    "source_case_id": receipt["member_id"],
                    "artifact_id": None,
                    "authority": "matrix_case",
                    "saturation_evidence_eligible": False,
                }
            },
        },
    )
    monkeypatch.setattr(
        admission,
        "_raw_bundle",
        lambda *args, **kwargs: (
            reduction,
            {
                "raw_inventory_sha256": receipt["raw_inventory_sha256"],
                "retained_raw_bindings": [],
                "reduction_binding": {"record_type": "fixture"},
            },
        ),
    )
    artifact_root.chmod(0o555)
    try:
        record, bundle = admission.derive_case_bundle(
            case_id=receipt["member_id"],
            artifact_root=artifact_root,
            q043_qualification_path=q043_path,
            q043_artifact_root=q043_root,
            q023_qualification_path=q023_path,
            authorized_orion_root=orion,
            authorized_project_home_root=project_home,
        )
        assert record["q043_dependency"]["sha256"] == q043_sha256
        assert record["execution_lineage"]["source_commit"] == receipt["source_commit"]
        assert record["execution_lineage"]["installed_producer"][
            "launch_trampoline_sha256"
        ] == receipt["producer"]["launch_trampoline_sha256"]
        assert record["problem_reported_saturation_evidence_eligible"] is True
        assert record["saturation_evidence_eligible"] is False
        assert bundle["snapshots"] == reduction["snapshots"]
        q023_record["q043_dependency"]["registered_matrix_sha256"] = "2" * 64
        with pytest.raises(admission.RegisteredAdmissionError, match="selected Q043"):
            admission.derive_case_bundle(
                case_id=receipt["member_id"],
                artifact_root=artifact_root,
                q043_qualification_path=q043_path,
                q043_artifact_root=q043_root,
                q023_qualification_path=q023_path,
                authorized_orion_root=orion,
                authorized_project_home_root=project_home,
            )
    finally:
        artifact_root.chmod(0o755)
