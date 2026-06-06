#!/usr/bin/env python3
"""Filesystem-backed adversarial tests for the Q-011 admission successor."""

from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest
import uuid

try:
    from tst.publication import (
        q011_section54_production_science_admission_orchestration_successor_v1 as successor,
    )
except ModuleNotFoundError:
    import q011_section54_production_science_admission_orchestration_successor_v1 as successor


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_production_science_admission_orchestration_successor_v1_2026-06-06.json"
)


def _bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _write(path: Path, payload: bytes, *, executable: bool = False) -> dict[str, str]:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    if executable:
        path.chmod(0o755)
    return {"path": str(path), "sha256": _sha(payload)}


def _fnv1a64(payload: bytes) -> int:
    value = 0xCBF29CE484222325
    for byte in payload:
        value ^= byte
        value = (value * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return value


def _marker(payload: bytes) -> bytes:
    return (
        "ATHENAK_RESTART_COMPLETE_V1\n"
        f"size={len(payload)}\n"
        f"fnv1a64={_fnv1a64(payload):016x}\n"
    ).encode("ascii")


def _closure(root: Path, paths: list[str]) -> dict[str, object]:
    members = []
    for relative in sorted(paths):
        path = root / relative
        if not path.exists():
            _write(path, f"fixture:{relative}\n".encode("ascii"))
        members.append({"path": relative, "sha256": _sha(path.read_bytes())})
    return {"members": members, "sha256": successor.canonical_sha256(members)}


def _identity(variant: str, seed: int) -> dict[str, object]:
    return {
        "campaign_id": successor.CAMPAIGN_ID,
        "attempt_id": f"q011-ps-v1-{variant}-{seed}",
        "variant": variant,
        "seed": seed,
        "physical_mode": successor.PHYSICAL_MODE,
    }


class Fixture:
    def __init__(self, root: Path):
        self.root = root
        self.source = root / "source"
        self.evidence = root / "evidence"
        deck_source = REPO_ROOT / successor.model.BASE_DECK_PATH
        deck_payload = deck_source.read_bytes()
        deck = _write(self.source / successor.SUCCESSOR_DECK, deck_payload)
        environment_payload = b"#!/bin/bash\n# fixture environment\n"
        reviewed_environment_path = (
            "tst/publication/frontier_control_plane/frontier_pic_environment.sh"
        )
        _write(self.source / reviewed_environment_path, environment_payload)
        source_closure = _closure(self.source, list(successor.REQUIRED_SOURCE_PATHS))
        reducer_closure = _closure(self.source, list(successor.REQUIRED_REDUCER_PATHS))
        executable = _write(self.evidence / "candidate" / "athena", b"\x7fELFfixture", executable=True)
        environment = _write(
            self.evidence / "controller" / "frontier_pic_environment.sh",
            environment_payload,
            executable=True,
        )
        environment["control_plane_version"] = "4" * 64
        environment["reviewed_source"] = {
            "path": reviewed_environment_path,
            "sha256": _sha((self.source / reviewed_environment_path).read_bytes()),
        }
        clean_manifest_value = {
            "schema_version": 4,
            "source": {
                "git_commit": "1" * 40,
                "source_bundle_sha256": "2" * 64,
            },
            "build": {"executable_path": executable["path"]},
        }
        clean_manifest = _write(
            self.evidence / "candidate" / "clean_candidate_manifest.json",
            _bytes(clean_manifest_value),
        )
        self.candidate = {
            "git_commit": "1" * 40,
            "source_bundle_sha256": "2" * 64,
            "clean_candidate_manifest": clean_manifest,
            "executable": executable,
            "environment_profile": environment,
            "deck": {"path": successor.SUCCESSOR_DECK, "sha256": deck["sha256"]},
            "source_closure": source_closure,
            "reducer_closure": reducer_closure,
        }
        self.pressure = _write(
            self.evidence / "pressure_selection_receipt.json",
            _bytes(
                {
                    "schema_version": 3,
                    "record_type": "q011_section54_pressure_selection_receipt",
                    "selection_method": "human_review_only",
                    "selected_case": {
                        "case_id": "ps_p0_1p00",
                        "problem_ps_p0": 1.0,
                    }
                }
            ),
        )

    def _execution(self, identity: dict[str, object]) -> dict[str, object]:
        attempt = str(identity["attempt_id"])
        raw = self.root / "runs" / attempt / "raw"
        raw.mkdir(parents=True)
        artifact_dir = self.root / "artifacts" / attempt
        attempt_digest = hashlib.sha256(attempt.encode("ascii")).hexdigest()
        argv = [
            "-i",
            successor.SUCCESSOR_DECK_LAUNCH_PATH,
            "-d",
            str(raw),
            f"job/basename={attempt}",
            "problem/ps_p0=1.0",
            *(f"{name}={identity['seed']}" for name in successor.SEED_OVERRIDE_NAMES),
            *successor.model.variant_binding(identity["variant"]).model_launch_overrides,
        ]
        contract = {
            "attempt_id": attempt,
            "variant": identity["variant"],
            "qualifying_seed": identity["seed"],
            "selected_problem_ps_p0": 1.0,
            "launch_authorized": False,
            "scheduler_submission_authorized": False,
            "executable": self.candidate["executable"],
            "environment_profile": self.candidate["environment_profile"],
            "paper_deck": self.candidate["deck"],
            "authorized_orion_attempt_root": str(raw.parent),
            "argv": argv,
        }
        contract_binding = _write(
            self.evidence / attempt / "attempt_contract.json", _bytes(contract)
        )
        receipt = {
            "record_type": "q011_section54_reconciled_registered_execution_receipt",
            "schema_version": 1,
            "receipt_role": "immutable_reconciled_registered_execution",
            "registration_scope": "registered_science",
            "reconciled": True,
            "reservation_id": str(uuid.uuid5(uuid.NAMESPACE_URL, f"reservation:{attempt}")),
            "submission_id": str(uuid.uuid5(uuid.NAMESPACE_URL, f"submission:{attempt}")),
            "reconciliation_event_sha256": attempt_digest,
            "attempt_id": attempt,
            "source_commit": self.candidate["git_commit"],
            "executable_sha256": self.candidate["executable"]["sha256"],
            "deck_sha256": self.candidate["deck"]["sha256"],
            "environment_sha256": self.candidate["environment_profile"]["sha256"],
            "control_plane_version": self.candidate["environment_profile"][
                "control_plane_version"
            ],
            "argv": argv,
            "slurm_terminal_state": "COMPLETED",
            "slurm_job_id": str(int(attempt_digest[:12], 16) + 1),
            "raw_output_root": str(raw),
            "artifact_dir": str(artifact_dir),
            "planner_retention": {
                "attempt_id": attempt,
                "authorized_orion_raw_root": str(raw),
                "argv": argv,
            },
            "pre_submit_manifest_sha256": hashlib.sha256(
                f"pre-submit:{attempt}".encode("ascii")
            ).hexdigest(),
        }
        receipt_binding = _write(
            self.evidence / attempt / "registered_execution_receipt.json",
            _bytes(receipt),
        )
        return {
            "raw_output_root": str(raw),
            "artifact_dir": str(artifact_dir),
            "attempt_contract": contract_binding,
            "registered_execution_receipt": receipt_binding,
            "selected_pressure_receipt": self.pressure,
        }

    def _raw(
        self, identity: dict[str, object], raw_root: Path
    ) -> list[dict[str, object]]:
        artifacts = []
        for index, slot in enumerate(successor.REQUIRED_NOMINAL_SLOTS):
            observed = 500.05 if slot == 500.0 else slot
            cycle = 10 * index
            suffix = f"{int(slot):05d}"
            for kind in successor.SNAPSHOT_PRODUCT_KINDS:
                if kind == "prtcl_all":
                    relative = f"pvtk/q011.prtcl_all.{suffix}.part.vtk"
                else:
                    relative = f"bin/q011.{kind}.{suffix}.bin"
                payload = f"{identity['attempt_id']}:{kind}:{slot}\n".encode("ascii")
                _write(raw_root / relative, payload)
                artifacts.append(
                    self._artifact(identity, relative, kind, payload, slot, cycle, observed)
                )
            restart_relative = f"rst/rank_00000000/q011.{suffix}.rst"
            restart_payload = (
                f"<time>\ncycle={cycle}\n<par_end>\n".encode("ascii") + b"binary"
            )
            restart_marker = _marker(restart_payload)
            manifest_relative = f"rst/q011.{suffix}.rst.manifest"
            manifest_payload = _bytes(
                {
                    "schema": "ATHENAK_RESTART_MANIFEST_V1",
                    "members": [
                        {
                            "path": restart_relative,
                            "size": len(restart_payload),
                            "fnv1a64": f"{_fnv1a64(restart_payload):016x}",
                        }
                    ],
                }
            )
            manifest_marker = _marker(manifest_payload)
            for relative, payload in (
                (restart_relative, restart_payload),
                (f"{restart_relative}.complete", restart_marker),
                (manifest_relative, manifest_payload),
                (f"{manifest_relative}.complete", manifest_marker),
            ):
                _write(raw_root / relative, payload)
                artifacts.append(
                    self._artifact(identity, relative, "rst", payload, slot, cycle, observed)
                )
        hst = f"{identity['attempt_id']}:hst\n".encode("ascii")
        stdout = b"time=1200 cycle=120\ntlim=1200 nlim=-1\nTerminating on time limit\n"
        _write(raw_root / "hst/q011.hst", hst)
        _write(raw_root / "stdout.txt", stdout)
        artifacts.append(self._artifact(identity, "hst/q011.hst", "hst", hst, None, None, None))
        artifacts.append(self._artifact(identity, "stdout.txt", "stdout", stdout, None, None, None))
        return artifacts

    @staticmethod
    def _artifact(
        identity: dict[str, object],
        path: str,
        kind: str,
        payload: bytes,
        slot: float | None,
        cycle: int | None,
        observed: float | None,
    ) -> dict[str, object]:
        return {
            "path": path,
            "kind": kind,
            "sha256": _sha(payload),
            "byte_count": len(payload),
            "attempt_id": identity["attempt_id"],
            "variant": identity["variant"],
            "seed": identity["seed"],
            "nominal_slot_time": slot,
            "cycle": cycle,
            "observed_committed_time": observed,
        }

    def attempt(
        self,
        variant: str = successor.GRID_VARIANTS[0],
        seed: int = successor.QUALIFYING_SEEDS[0],
    ) -> dict[str, object]:
        identity = _identity(variant, seed)
        execution = self._execution(identity)
        raw = self._raw(identity, Path(execution["raw_output_root"]))
        return successor.build_attempt_admission(
            attempt_identity=identity,
            candidate_binding=self.candidate,
            execution_binding=execution,
            raw_artifacts=raw,
            source_root=self.source,
        )

    def inputs(
        self, variant: str, seed: int
    ) -> tuple[dict[str, object], dict[str, object], list[dict[str, object]]]:
        identity = _identity(variant, seed)
        execution = self._execution(identity)
        raw = self._raw(identity, Path(execution["raw_output_root"]))
        return identity, execution, raw


def _rewrite_binding(binding: dict[str, str], value: object) -> None:
    payload = _bytes(value) if not isinstance(value, bytes) else value
    Path(binding["path"]).write_bytes(payload)
    binding["sha256"] = _sha(payload)


def _rewrite_raw_artifact(
    execution: dict[str, object], raw: list[dict[str, object]], path: str, payload: bytes
) -> None:
    Path(execution["raw_output_root"], path).write_bytes(payload)
    artifact = next(item for item in raw if item["path"] == path)
    artifact["sha256"] = _sha(payload)
    artifact["byte_count"] = len(payload)


class Q011ProductionScienceAdmissionSuccessorV1Tests(unittest.TestCase):
    def test_readiness_binds_new_files_and_refuses_all_authority(self) -> None:
        record = json.loads(READINESS.read_text(encoding="utf-8"))
        for binding in record["source_bindings"].values():
            path = REPO_ROOT / binding["path"]
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), binding["sha256"])
        self.assertTrue(
            all(
                value is False
                for key, value in record["authorization"].items()
                if key.endswith("_authorized")
            )
        )

    def test_attempt_verifies_bytes_lineage_slot_publication_and_terminal_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            record = fixture.attempt()
            self.assertEqual(len(record["raw_artifacts"]), 145)
            self.assertEqual(len(record["snapshot_bindings"]), 13)
            self.assertEqual(
                record["run_bindings"]["terminal_completion"]["termination_reason"],
                "Terminating on time limit",
            )
            self.assertEqual(
                successor.validate_attempt_admission(record, source_root=fixture.source),
                record,
            )

    def test_policy_uses_explicit_deposited_j_over_c_and_complete_hooks(self) -> None:
        policy = successor.analysis_policy()
        self.assertEqual(policy["current_profiles"]["deposited_representation"], "J_CR_over_c")
        self.assertIn("conserved_energy_closure", policy)
        self.assertIn("morphology", policy)
        self.assertEqual(
            policy["shock_front_preregistration_supersession"]["successor_detector"],
            "unique_strongest_negative_density_gradient",
        )

    def test_candidate_execution_pressure_and_launch_provenance_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            bad_candidate = copy.deepcopy(fixture.candidate)
            bad_candidate["executable"]["path"] = str(fixture.root / "missing-athena")
            manifest = json.loads(
                Path(bad_candidate["clean_candidate_manifest"]["path"]).read_text()
            )
            manifest["build"]["executable_path"] = bad_candidate["executable"]["path"]
            _rewrite_binding(bad_candidate["clean_candidate_manifest"], manifest)
            with self.assertRaisesRegex(
                successor.ProductionScienceAdmissionError, "unavailable"
            ):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=bad_candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            bad_candidate = copy.deepcopy(fixture.candidate)
            bad_candidate["executable"]["sha256"] = "9" * 64
            with self.assertRaisesRegex(
                successor.ProductionScienceAdmissionError, "SHA-256 drifted"
            ):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=bad_candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            contract = json.loads(Path(execution["attempt_contract"]["path"]).read_text())
            contract["argv"][-1:] = ["problem/forged=true"]
            _rewrite_binding(execution["attempt_contract"], contract)
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "attempt contract"):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

    def test_deck_environment_source_closure_and_clean_manifest_bytes_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            (fixture.source / successor.SUCCESSOR_DECK).write_bytes(b"forged deck\n")
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "SHA-256 drifted"):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            Path(fixture.candidate["environment_profile"]["path"]).write_bytes(b"forged environment\n")
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "SHA-256 drifted"):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            source_member = next(iter(successor.REQUIRED_SOURCE_PATHS))
            (fixture.source / source_member).write_bytes(b"forged source\n")
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "SHA-256 drifted"):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            manifest = json.loads(
                Path(fixture.candidate["clean_candidate_manifest"]["path"]).read_text()
            )
            manifest["source"]["git_commit"] = "9" * 40
            _rewrite_binding(fixture.candidate["clean_candidate_manifest"], manifest)
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "cross-link drifted"):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

    def test_registered_receipt_and_selected_pressure_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[1], successor.QUALIFYING_SEEDS[0]
            )
            receipt = json.loads(
                Path(execution["registered_execution_receipt"]["path"]).read_text()
            )
            receipt["deck_sha256"] = "9" * 64
            _rewrite_binding(execution["registered_execution_receipt"], receipt)
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "registered execution"):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[1], successor.QUALIFYING_SEEDS[0]
            )
            pressure = {
                "schema_version": 3,
                "record_type": "q011_section54_pressure_selection_receipt",
                "selection_method": "human_review_only",
                "selected_case": {"case_id": "ps_p0_0p10", "problem_ps_p0": 0.1},
            }
            _rewrite_binding(execution["selected_pressure_receipt"], pressure)
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "ps_p0_1p00"):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

    def test_relabeling_one_run_as_another_cell_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            relabeled = _identity(successor.GRID_VARIANTS[2], successor.QUALIFYING_SEEDS[1])
            for artifact in raw:
                artifact["attempt_id"] = relabeled["attempt_id"]
                artifact["variant"] = relabeled["variant"]
                artifact["seed"] = relabeled["seed"]
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "attempt contract"):
                successor.build_attempt_admission(
                    attempt_identity=relabeled,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

    def test_every_interior_float32_slot_assignment_is_enforced(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            for artifact in raw:
                if artifact["nominal_slot_time"] == 100.0:
                    artifact["observed_committed_time"] = 1.0
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "float32 slot assignment"):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

    def test_restart_header_marker_manifest_and_stdout_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            target = next(
                item
                for item in raw
                if item["path"].endswith("00500.rst") and not item["path"].endswith(".complete")
            )
            target["cycle"] = 999
            for item in raw:
                if item["nominal_slot_time"] == 500.0:
                    item["cycle"] = 999
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "restart header cycle"):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            manifest_path = "rst/q011.00500.rst.manifest"
            manifest = json.loads(Path(execution["raw_output_root"], manifest_path).read_text())
            manifest["members"][0]["path"] = "rst/rank_00000000/unrelated.00500.rst"
            payload = _bytes(manifest)
            _rewrite_raw_artifact(execution, raw, manifest_path, payload)
            _rewrite_raw_artifact(
                execution, raw, f"{manifest_path}.complete", _marker(payload)
            )
            with self.assertRaisesRegex(
                successor.ProductionScienceAdmissionError,
                "does not bind every payload",
            ):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            marker_path = "rst/rank_00000000/q011.00500.rst.complete"
            _rewrite_raw_artifact(execution, raw, marker_path, b"forged marker\n")
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "completion marker"):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            identity, execution, raw = fixture.inputs(
                successor.GRID_VARIANTS[0], successor.QUALIFYING_SEEDS[0]
            )
            _rewrite_raw_artifact(execution, raw, "stdout.txt", b"time=1200 cycle=120\n")
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "terminal"):
                successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )

    def test_core_and_telemetry_only_expansion_emit_complete_validated_graph(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            core = [
                fixture.attempt(variant, seed)
                for variant in successor.GRID_VARIANTS
                for seed in successor.CORE_QUALIFYING_SEEDS
            ]
            record = successor.build_campaign_orchestration(
                core, source_root=fixture.source
            )
            graph = record["work_graph"]
            self.assertEqual(len(record["attempt_admission_bindings"]), 9)
            self.assertEqual(len(graph["downstream_spectrum_jobs"]), 18)
            self.assertEqual(len(graph["current_profile_jobs"]), 18)
            self.assertEqual(len(graph["paired_grid_amr_residual_jobs"]), 3)
            self.assertEqual(len(graph["conserved_energy_closure_jobs"]), 9)
            self.assertEqual(
                successor.validate_campaign_orchestration(
                    record, core, source_root=fixture.source
                ),
                record,
            )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            seeds = successor.QUALIFYING_SEEDS[:4]
            attempts = [
                fixture.attempt(variant, seed)
                for variant in successor.GRID_VARIANTS
                for seed in seeds
            ]
            evidence = {
                "record_type": "q011_resource_telemetry_expansion_decision_v1",
                "selected_phase_id": "phase_4_paired_triads",
                "previous_phase_id": "phase_3_paired_triads",
                "decision_basis_fields": list(successor.RESOURCE_TELEMETRY_FIELDS),
                "qualifying_output_inspected": False,
                "resource_model_sha256": "8" * 64,
            }
            expanded = successor.build_campaign_orchestration(
                attempts,
                source_root=fixture.source,
                phase_id="phase_4_paired_triads",
                resource_expansion_evidence=evidence,
            )
            self.assertEqual(expanded["campaign_matrix"]["expected_attempt_count"], 12)
            bad = copy.deepcopy(evidence)
            bad["qualifying_output_inspected"] = True
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "telemetry-only"):
                successor.build_campaign_orchestration(
                    attempts,
                    source_root=fixture.source,
                    phase_id="phase_4_paired_triads",
                    resource_expansion_evidence=bad,
                )

    def test_incomplete_phase_reused_execution_and_persisted_graph_tampering_fail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            attempts = [
                fixture.attempt(variant, seed)
                for variant in successor.GRID_VARIANTS
                for seed in successor.CORE_QUALIFYING_SEEDS
            ]
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "complete 9-attempt"):
                successor.build_campaign_orchestration(
                    attempts[:-1], source_root=fixture.source
                )
            record = successor.build_campaign_orchestration(
                attempts, source_root=fixture.source
            )
            record["work_graph"]["current_profile_jobs"].pop()
            with self.assertRaisesRegex(successor.ProductionScienceAdmissionError, "persisted work graph"):
                successor.validate_campaign_orchestration(
                    record, attempts, source_root=fixture.source
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            attempts = []
            first_registered_identity = None
            cells = [
                (variant, seed)
                for variant in successor.GRID_VARIANTS
                for seed in successor.CORE_QUALIFYING_SEEDS
            ]
            for index, (variant, seed) in enumerate(cells):
                identity, execution, raw = fixture.inputs(variant, seed)
                if index == len(cells) - 1:
                    receipt = json.loads(
                        Path(execution["registered_execution_receipt"]["path"]).read_text()
                    )
                    receipt.update(first_registered_identity)
                    _rewrite_binding(execution["registered_execution_receipt"], receipt)
                attempt = successor.build_attempt_admission(
                    attempt_identity=identity,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=fixture.source,
                )
                if first_registered_identity is None:
                    first_registered_identity = dict(
                        attempt["execution_binding"]["registered_execution_identity"]
                    )
                attempts.append(attempt)
            with self.assertRaisesRegex(
                successor.ProductionScienceAdmissionError,
                "reused execution evidence",
            ):
                successor.build_campaign_orchestration(
                    attempts, source_root=fixture.source
                )


if __name__ == "__main__":
    unittest.main()
