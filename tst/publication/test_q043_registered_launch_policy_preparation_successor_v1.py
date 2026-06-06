#!/usr/bin/env python3
"""Focused adversarial tests for Q043 registered-launch/policy preparation."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import stat
import tempfile
import unittest

from tst.publication import (
    q043_registered_launch_policy_preparation_successor_v1 as preparation,
)


def _json(payload: bytes) -> dict[str, object]:
    value = json.loads(payload)
    assert isinstance(value, dict)
    return value


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _final_bindings() -> dict[str, object]:
    version = "a" * 64
    candidate = preparation.AUTHORIZED_ORION_ROOT / "clean_candidates/final-q043"
    controller = preparation.AUTHORIZED_ORION_ROOT / "control_plane" / version
    project_controller = preparation.CANONICAL_PROJECT_HOME_ROOT / "control_plane" / version
    return {
        "record_type": preparation.FINAL_BINDING_RECORD_TYPE,
        "schema_version": 1,
        "source_commit": "b" * 40,
        "source_bundle_sha256": "c" * 64,
        "source_archive_path": str(candidate / "source.tar"),
        "source_archive_sha256": "d" * 64,
        "clean_candidate_manifest_path": str(candidate / "clean_candidate_manifest.json"),
        "clean_candidate_manifest_sha256": "e" * 64,
        "executable_path": str(candidate / "athena"),
        "executable_sha256": "f" * 64,
        "installed_control_plane_version": version,
        "orion_installed_control_plane_root": str(controller),
        "project_home_installed_control_plane_root": str(project_controller),
        "environment_profile_path": str(controller / "frontier_pic_environment.sh"),
        "environment_profile_sha256": "1" * 64,
        "job_script_path": str(controller / "frontier_job.sh"),
        "job_script_sha256": "2" * 64,
        "analysis_script_paths": [
            str(candidate / "q043_registered_execution_raw_oracle_qualification_successor_v1.py"),
            str(candidate / "q043_bell_current_volume_aware_deposited_current_oracle.py"),
        ],
        "analysis_script_sha256": ["4" * 64, "5" * 64],
        "reconcile_q043_registered_execution_path": str(
            controller / "reconcile_q043_registered_execution.py"
        ),
        "reconcile_q043_registered_execution_sha256": "3" * 64,
    }


class Q043RegisteredLaunchPolicyPreparationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest, cls.files = preparation.build_materialization()
        cls.launches = [
            _json(cls.files[record["launch_candidate"]["path"]])
            for record in cls.manifest["case_records"]
        ]
        cls.policies = [
            _json(cls.files[record["policy_candidate"]["path"]])
            for record in cls.manifest["case_records"]
        ]

    def test_exact_132_case_candidates_are_deterministic_and_schema_valid(self) -> None:
        self.assertEqual(self.manifest["case_count"], 132)
        self.assertEqual(self.manifest["launch_candidate_count"], 132)
        self.assertEqual(self.manifest["policy_candidate_count"], 132)
        self.assertEqual(len(self.files), 267)
        second_manifest, second_files = preparation.build_materialization()
        self.assertEqual(second_manifest, self.manifest)
        self.assertEqual(second_files, self.files)
        preparation.validate_materialization(self.manifest, self.files)
        compatibility = self.manifest["launch_contract_schema_compatibility"]
        self.assertEqual(
            compatibility["checked_in_control_plane_schema_source_validation"],
            "passed_for_all_132_candidates",
        )
        self.assertTrue(
            compatibility["final_installed_generation_schema_validation_required"]
        )
        self.assertFalse(compatibility["final_installed_generation_assumed"])
        for launch in self.launches:
            contract = launch["launch_contract_candidate"]
            self.assertEqual(preparation.validate_launch_contract(contract), contract)
            self.assertEqual(
                launch["launch_contract_sha256"],
                preparation.launch_contract_sha256(contract),
            )

    def test_case_identity_rank_decomposition_command_and_deck_are_exact(self) -> None:
        expected = preparation.validate_case_matrix(
            preparation._checked_in_deck_manifest()["cases"]
        )
        for index, (case, launch) in enumerate(zip(expected, self.launches), 1):
            self.assertEqual(launch["identity"]["case_index"], index)
            self.assertEqual(launch["identity"]["case_id"], case["case_id"])
            self.assertEqual(launch["identity"]["maximum_registered_attempts"], 1)
            self.assertEqual(launch["identity"]["maximum_retries"], 0)
            mpi = launch["mpi_and_decomposition"]
            self.assertEqual(mpi["mpi_ranks"], case["mpi_ranks"])
            self.assertEqual(mpi["meshblock_grid"], case["meshblock_grid"])
            self.assertEqual(len(mpi["rank_meshblock_ids"]), case["mpi_ranks"])
            command = launch["command_template"]
            self.assertIn(f"--ntasks={case['mpi_ranks']}", command)
            self.assertIn(f"--ntasks-per-node={case['mpi_ranks']}", command)
            deck = preparation.REPO_ROOT / launch["checked_in_deck"]["path"]
            self.assertEqual(_sha256(deck), launch["checked_in_deck"]["sha256"])

    def test_output_topology_is_exact_single_file_per_rank_contract(self) -> None:
        expected_total = 0
        for launch in self.launches:
            ranks = launch["case_contract"]["mpi_ranks"]
            topology = launch["output_topology"]
            self.assertEqual(topology["single_file_per_rank"], ranks > 1)
            self.assertEqual(topology["single_rank_file_topology"], ranks == 1)
            self.assertTrue(topology["one_meshblock_per_rank"])
            self.assertEqual(topology["expected_raw_artifact_count"], 8 * ranks)
            self.assertEqual(
                len(topology["expected_relative_raw_paths"]),
                topology["expected_raw_artifact_count"],
            )
            self.assertEqual(
                len(set(topology["expected_relative_raw_paths"])),
                topology["expected_raw_artifact_count"],
            )
            expected_total += topology["maximum_case_storage_bytes"]
        budget = _json(self.files["batch_budget_accounting_input.json"])
        self.assertEqual(budget["maximum_batch_storage_bytes"], expected_total)
        self.assertLessEqual(
            budget["maximum_batch_storage_bytes"],
            budget["maximum_storage_cap_bytes"],
        )

    def test_wrapper_and_future_trusted_paired_receipt_interface_are_bound(self) -> None:
        for launch in self.launches:
            case_id = launch["identity"]["case_id"]
            ranks = launch["case_contract"]["mpi_ranks"]
            wrapper = launch["stdout_wrapper_evidence"]
            self.assertEqual(
                wrapper["required_exact_rank_line"],
                (
                    f"Q043_REGISTERED_EXECUTION case_id={case_id} "
                    f"mpi_world_size={ranks} rank_ids="
                    + ",".join(str(rank) for rank in range(ranks))
                ),
            )
            self.assertEqual(
                wrapper["required_exact_exit_line"],
                "Q043_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0",
            )
            reconcile = launch["trusted_reconciliation_interface"]
            self.assertEqual(
                reconcile["producer_name"], "reconcile_q043_registered_execution.py"
            )
            self.assertTrue(reconcile["paired_receipts_required"])
            self.assertTrue(reconcile["paired_receipts_must_be_byte_identical"])
            self.assertIn(
                str(preparation.CANONICAL_PROJECT_HOME_ROOT),
                reconcile["canonical_project_home_receipt_path_template"],
            )
            self.assertIn(case_id, reconcile["orion_receipt_path_template"])

    def test_budget_inputs_are_exact_and_fail_closed(self) -> None:
        budget = _json(self.files["batch_budget_accounting_input.json"])
        self.assertEqual(budget["maximum_registered_attempts"], 132)
        self.assertEqual(budget["maximum_batch_node_hours"], 22.0)
        self.assertEqual(budget["maximum_retries_per_case"], 0)
        self.assertEqual(budget["maximum_live_q043_submissions"], 1)
        self.assertFalse(budget["ledger_mutation_authorized"])
        preparation.validate_budget_accounting(budget, self.launches)
        drifted = copy.deepcopy(budget)
        drifted["maximum_batch_node_hours"] = 10001.0
        with self.assertRaisesRegex(preparation.PreparationError, "node-hour ceiling"):
            preparation.validate_budget_accounting(drifted, self.launches)
        drifted = copy.deepcopy(budget)
        drifted["maximum_batch_storage_bytes"] = preparation.MAXIMUM_BATCH_STORAGE_BYTES + 1
        with self.assertRaisesRegex(preparation.PreparationError, "storage ceiling"):
            preparation.validate_budget_accounting(drifted, self.launches)

    def test_explicit_queue_clean_candidate_and_installed_controller_blockers(self) -> None:
        required = set(preparation.REQUIRED_BLOCKERS)
        for token in {
            "empty_user_queue_proved_at_fresh_reservation_boundary",
            "final_clean_candidate_manifest_source_archive_and_executable_bound",
            "final_paired_installed_control_plane_generation_bound_and_independently_verified",
        }:
            self.assertIn(token, required)
        self.assertFalse(self.manifest["execution_boundary"]["empty_user_queue_assumed"])
        self.assertFalse(self.manifest["execution_boundary"]["final_clean_candidate_assumed"])
        self.assertFalse(
            self.manifest["execution_boundary"]["final_installed_control_plane_assumed"]
        )
        for launch, policy in zip(self.launches, self.policies):
            self.assertEqual(set(launch["required_blockers"]), required)
            self.assertEqual(set(policy["required_blockers"]), required)

    def test_missing_duplicate_relabelled_and_wrong_rank_matrix_fail_closed(self) -> None:
        cases = [dict(case) for case in preparation.oracle.expected_cases()]
        for label, mutation in {
            "missing": cases[:-1],
            "duplicate": [*cases[:-1], copy.deepcopy(cases[-2])],
        }.items():
            with self.subTest(label), self.assertRaisesRegex(
                preparation.PreparationError, "exact ordered 132-case"
            ):
                preparation.validate_case_matrix(mutation)
        relabelled = copy.deepcopy(cases)
        relabelled[0]["case_id"] = "q043-current-oracle-relabelled"
        with self.assertRaisesRegex(preparation.PreparationError, "exact ordered 132-case"):
            preparation.validate_case_matrix(relabelled)
        wrong_rank = copy.deepcopy(cases)
        wrong_rank[-1]["mpi_ranks"] = 4
        with self.assertRaisesRegex(preparation.PreparationError, "exact ordered 132-case"):
            preparation.validate_case_matrix(wrong_rank)
        wrong_decomposition = copy.deepcopy(cases)
        wrong_decomposition[-1]["meshblock_grid"] = [2, 2, 1]
        with self.assertRaisesRegex(preparation.PreparationError, "exact ordered 132-case"):
            preparation.validate_case_matrix(wrong_decomposition)

    def test_wrong_deck_and_wrong_output_topology_fail_closed(self) -> None:
        files = dict(self.files)
        path = self.manifest["case_records"][0]["launch_candidate"]["path"]
        launch = _json(files[path])
        launch["checked_in_deck"]["sha256"] = "0" * 64
        files[path] = preparation._json_bytes(launch)
        with self.assertRaisesRegex(preparation.PreparationError, "file drifted"):
            preparation.validate_materialization(self.manifest, files)
        files = dict(self.files)
        launch = _json(files[path])
        launch["output_topology"]["single_file_per_rank"] = True
        files[path] = preparation._json_bytes(launch)
        with self.assertRaisesRegex(preparation.PreparationError, "file drifted"):
            preparation.validate_materialization(self.manifest, files)

    def test_authority_drift_fails_closed(self) -> None:
        files = dict(self.files)
        path = self.manifest["case_records"][0]["policy_candidate"]["path"]
        policy = _json(files[path])
        policy["authorization"]["launch_authorized"] = True
        files[path] = preparation._json_bytes(policy)
        with self.assertRaisesRegex(preparation.PreparationError, "file drifted"):
            preparation.validate_materialization(self.manifest, files)
        for value in [self.manifest, *self.launches, *self.policies]:
            text = json.dumps(value, sort_keys=True)
            self.assertNotIn('"launch_authorized": true', text)
            self.assertNotIn('"publication_authorized": true', text)

    def test_final_digest_placeholder_rules_fail_closed(self) -> None:
        final = _final_bindings()
        preparation.validate_final_bindings(final)
        manifest, files = preparation.build_materialization(final)
        self.assertEqual(
            manifest["binding_stage"], "final_digests_bound_review_candidate_still_blocked"
        )
        self.assertNotIn(
            preparation.FINAL_PLACEHOLDER_PREFIX,
            json.dumps([manifest, *[payload.decode() for payload in files.values()]]),
        )
        preparation.validate_materialization(manifest, files, final)
        drifted = copy.deepcopy(final)
        drifted["executable_sha256"] = "PENDING_FINAL_EXECUTABLE_SHA256"
        with self.assertRaisesRegex(preparation.PreparationError, "retain"):
            preparation.validate_final_bindings(drifted)
        drifted = copy.deepcopy(final)
        drifted["project_home_installed_control_plane_root"] = (
            "/ccs/proj/ast207/proj-shared/PIC/control_plane/" + "a" * 64
        )
        with self.assertRaisesRegex(preparation.PreparationError, "paired installed"):
            preparation.validate_final_bindings(drifted)

    def test_policy_candidates_are_incomplete_non_authorizing_fragments(self) -> None:
        for launch, policy in zip(self.launches, self.policies):
            self.assertEqual(
                launch["binding_stage"], "source_local_preparation_pending_final_digests"
            )
            self.assertTrue(launch["unresolved_final_bindings"])
            self.assertEqual(
                policy["status"], "review_required_not_authorized_not_live_policy"
            )
            self.assertEqual(
                policy["policy_slice_candidate"]["status"],
                "review_required_not_authorized",
            )
            self.assertEqual(
                policy["policy_slice_candidate"]["analysis_script_sha256"],
                ["PENDING_FINAL_ANALYSIS_SCRIPT_SHA256"],
            )
            self.assertIn(
                "paired_reconciliation_receipt_sha256",
                policy["fresh_runtime_bindings_required"],
            )
            self.assertFalse(policy["authorization"]["policy_mutation_authorized"])
        fragment = _json(self.files["registered_policy_slice_review_fragment.json"])
        self.assertEqual(fragment["case_count"], 132)
        self.assertEqual(len(fragment["registered_science_slice_candidates"]), 132)
        self.assertFalse(fragment["execution_boundary"]["complete_storage_policy"])
        self.assertFalse(
            fragment["execution_boundary"]["live_policy_mutation_authorized"]
        )

    def test_materialized_bundle_is_read_only_and_non_authorizing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory).resolve() / "q043-review"
            result = preparation.materialize_review_bundle(output)
            self.assertEqual(result["case_count"], 132)
            self.assertEqual(result["maximum_batch_node_hours"], 22.0)
            self.assertTrue(result["recursively_read_only"])
            self.assertFalse(result["launch_authorized"])
            self.assertFalse(result["policy_mutation_authorized"])
            for path in [output, *output.rglob("*")]:
                self.assertFalse(path.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
            for path in [output, *output.rglob("*")]:
                path.chmod(path.stat().st_mode | stat.S_IWUSR)

    def test_readiness_record_binds_exact_source_and_non_authorizing_contract(self) -> None:
        path = (
            preparation.READINESS_ROOT
            / "q043_registered_launch_policy_preparation_successor_v1_2026-06-06.json"
        )
        readiness = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(readiness["successor_id"], preparation.SUCCESSOR_ID)
        self.assertEqual(readiness["exact_matrix"]["case_count"], 132)
        self.assertEqual(readiness["resource_ceiling"]["maximum_batch_node_hours"], 22.0)
        self.assertFalse(readiness["authorization"]["launch_authorized"])
        self.assertFalse(readiness["authorization"]["policy_mutation_authorized"])
        for binding in readiness["source_bindings"].values():
            source = preparation.REPO_ROOT / binding["path"]
            self.assertEqual(_sha256(source), binding["sha256"])


if __name__ == "__main__":
    unittest.main()
