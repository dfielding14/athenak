#!/usr/bin/env python3
"""Focused tests for Q023 registered launch and policy preparation."""

from __future__ import annotations

import copy
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tst.publication import q023_registered_launch_policy_preparation_successor_v1 as prep


class Q023RegisteredLaunchPreparationTests(unittest.TestCase):
    def test_pending_materialization_is_exact_and_non_authorizing(self) -> None:
        manifest, files = prep.build_materialization()
        prep.validate_materialization(manifest, files)
        self.assertEqual(manifest["case_count"], 55)
        self.assertEqual(len(manifest["policy_slices"]), 55)
        self.assertEqual(len(files), 113)
        self.assertEqual(
            sum(
                record["mpi_ranks"] for record in manifest["case_records"]
            ),
            150,
        )
        self.assertEqual(
            manifest["exact_matrix_contract"]["rank_distribution"],
            {"2": 45, "4": 5, "8": 5},
        )
        self.assertEqual(
            manifest["exact_matrix_contract"]["raw_artifacts_total"], 24_475
        )
        self.assertTrue(manifest["unresolved_final_bindings"])
        self.assertTrue(
            all(value is False for value in prep.AUTHORIZATION_BOUNDARY.values())
        )

    def test_launch_contracts_bind_exact_rank_count_and_no_override(self) -> None:
        manifest, files = prep.build_materialization()
        for record in manifest["case_records"]:
            launch = __import__("json").loads(
                files[record["launch_candidate"]["path"]]
            )
            action = launch["launch_contract"]["actions"][0]
            self.assertEqual(action["resources"]["tasks"], record["mpi_ranks"])
            self.assertEqual(
                launch["output_topology"]["expected_raw_artifact_count"], 445
            )
            self.assertEqual(
                len(set(launch["output_topology"]["expected_relative_raw_paths"])),
                445,
            )
            self.assertFalse(launch["output_topology"]["single_file_per_rank"])
            self.assertEqual(
                launch["stdout_wrapper_evidence"]["required_terminal_time"], 1.75
            )
            literals = [
                item["literal"]
                for item in action["arguments"]
                if set(item) == {"literal"}
            ]
            self.assertFalse(any(value.startswith("time/") for value in literals))

    def test_budget_and_policy_slices_fail_closed(self) -> None:
        manifest, files = prep.build_materialization()
        budget = __import__("json").loads(files["batch_budget_accounting_input.json"])
        self.assertEqual(budget["maximum_batch_node_hours"], 55.0)
        self.assertEqual(
            budget["maximum_batch_storage_bytes"], 220 * 1024**3
        )
        self.assertEqual(
            budget["maximum_storage_cap_bytes"], 256 * 1024**3
        )
        self.assertEqual(budget["maximum_attempts_per_case"], 1)
        self.assertEqual(budget["maximum_retries_per_case"], 0)
        altered = copy.deepcopy(manifest)
        altered["case_count"] = 54
        with self.assertRaisesRegex(prep.PreparationError, "manifest drifted"):
            prep.validate_materialization(altered, files)

    def test_strict_validation_rejects_boolean_integer_alias(self) -> None:
        manifest, files = prep.build_materialization()
        altered = copy.deepcopy(manifest)
        altered["case_count"] = True
        with self.assertRaisesRegex(prep.PreparationError, "manifest drifted"):
            prep.validate_materialization(altered, files)

    def test_final_materialization_rejects_pending_q043_deck_binding(self) -> None:
        version = "a" * 64
        candidate = prep.AUTHORIZED_ORION_ROOT / "clean_candidates/q023-final"
        controller = prep.AUTHORIZED_ORION_ROOT / "control_plane" / version
        project = prep.CANONICAL_PROJECT_HOME_ROOT / "control_plane" / version
        final = {
            "record_type": prep.FINAL_BINDING_RECORD_TYPE,
            "schema_version": 1,
            "source_commit": "b" * 40,
            "source_bundle_sha256": "c" * 64,
            "source_archive_path": str(candidate / "source.tar"),
            "source_archive_sha256": "d" * 64,
            "clean_candidate_manifest_path": str(
                candidate / "clean_candidate_manifest.json"
            ),
            "clean_candidate_manifest_sha256": "e" * 64,
            "executable_path": str(candidate / "athena"),
            "executable_sha256": "f" * 64,
            "installed_control_plane_version": version,
            "orion_installed_control_plane_root": str(controller),
            "project_home_installed_control_plane_root": str(project),
            "environment_profile_path": str(
                controller / "frontier_pic_environment.sh"
            ),
            "environment_profile_sha256": "1" * 64,
            "job_script_path": str(controller / "frontier_job.sh"),
            "job_script_sha256": "2" * 64,
            "analysis_script_paths": [
                str(controller / "reconcile_q023_registered_execution.py")
            ],
            "analysis_script_sha256": ["3" * 64],
            "reconcile_q023_registered_execution_path": str(
                controller / "reconcile_q023_registered_execution.py"
            ),
            "reconcile_q023_registered_execution_sha256": "3" * 64,
            "q043_registered_matrix_path": str(
                prep.AUTHORIZED_ORION_ROOT / "analysis/q043/matrix.json"
            ),
            "q043_registered_matrix_sha256": "4" * 64,
            "q043_registered_matrix_record_type": (
                "q043_registered_execution_raw_oracle_matrix_qualification"
            ),
            "q043_registered_matrix_case_bindings_sha256": "5" * 64,
            "q043_registered_dependency_sha256": "6" * 64,
        }
        prep.validate_final_bindings(final)
        with self.assertRaisesRegex(
            prep.PreparationError, "decks rebound to the exact Q043 matrix"
        ):
            prep.build_materialization(final)

    def test_final_binding_materializer_derives_exact_immutable_identities(
        self,
    ) -> None:
        version = "a" * 64
        commit = "b" * 40
        candidate = (
            prep.AUTHORIZED_ORION_ROOT
            / "clean_candidates/11111111-2222-3333-4444-555555555555"
        )
        manifest_path = candidate / "clean_candidate_manifest.json"
        source_archive = candidate / "source.tar"
        executable = candidate / "athena"
        source_payload = b"source archive"
        executable_payload = b"athena executable"
        controller_payloads = {
            "frontier_pic_environment.sh": b"environment",
            "frontier_job.sh": b"job",
            "reconcile_q023_registered_execution.py": b"reconciler",
        }
        candidate_value = {
            "source": {
                "archive_path": str(source_archive),
                "archive_sha256": prep._sha256_bytes(source_payload),
                "source_bundle_sha256": "c" * 64,
                "git_commit": commit,
            },
            "build": {
                "executable_path": str(executable),
                "executable_sha256": prep._sha256_bytes(executable_payload),
            },
        }
        manifest_payload = json.dumps(
            candidate_value, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
        manifest_sha256 = prep._sha256_bytes(manifest_payload)
        matrix_path = (
            prep.AUTHORIZED_ORION_ROOT
            / "analysis/q043_registered_execution_raw_oracle_successor_v1"
            / "q043_registered_matrix_qualification.json"
        )
        dependency = {
            "binding_kind": "registered_matrix_qualification",
            "registered_matrix_sha256": "d" * 64,
            "registered_matrix_record_type": (
                "q043_registered_execution_raw_oracle_matrix_qualification"
            ),
            "registered_matrix_case_bindings_sha256": "e" * 64,
        }

        def stable_bytes(path: Path, **_: object) -> tuple[object, bytes]:
            if path == manifest_path:
                return object(), manifest_payload
            if path == source_archive:
                return object(), source_payload
            if path == executable:
                return object(), executable_payload
            if path.name in controller_payloads:
                return object(), controller_payloads[path.name]
            raise AssertionError(f"unexpected stable read: {path}")

        report = {
            "status": "passed",
            "current_control_plane_version": version,
            "build": {"receipt_control_plane_version": version},
            "source": {"git_commit": commit},
        }
        with patch.object(
            prep, "revalidate_clean_candidate", return_value=report
        ), patch.object(
            prep.q043, "_stable_regular_bytes", side_effect=stable_bytes
        ), patch.object(
            prep.bell,
            "registered_q043_raw_oracle_dependency",
            return_value=dependency,
        ), patch.object(
            prep.bell, "_dependency_digest", return_value="f" * 64
        ), patch.object(
            prep,
            "validate_final_binding_files",
            side_effect=prep.validate_final_bindings,
        ):
            final = prep.materialize_q023_final_bindings(
                clean_candidate_manifest=manifest_path,
                clean_candidate_manifest_sha256=manifest_sha256,
                expected_source_commit=commit,
                installed_control_plane_version=version,
                q043_registered_matrix=matrix_path,
            )

        self.assertEqual(final["source_commit"], commit)
        self.assertEqual(
            final["source_archive_sha256"],
            prep._sha256_bytes(source_payload),
        )
        self.assertEqual(
            final["executable_sha256"],
            prep._sha256_bytes(executable_payload),
        )
        self.assertEqual(final["q043_registered_dependency_sha256"], "f" * 64)
        self.assertEqual(
            final["reconcile_q023_registered_execution_sha256"],
            prep._sha256_bytes(
                controller_payloads["reconcile_q023_registered_execution.py"]
            ),
        )

    def test_final_binding_materializer_rejects_revalidation_identity_drift(
        self,
    ) -> None:
        with patch.object(
            prep,
            "revalidate_clean_candidate",
            return_value={
                "status": "passed",
                "current_control_plane_version": "0" * 64,
                "build": {"receipt_control_plane_version": "a" * 64},
                "source": {"git_commit": "b" * 40},
            },
        ):
            with self.assertRaisesRegex(
                prep.PreparationError,
                "revalidation report drifted",
            ):
                prep.materialize_q023_final_bindings(
                    clean_candidate_manifest=(
                        prep.AUTHORIZED_ORION_ROOT
                        / "clean_candidates/11111111-2222-3333-4444-555555555555"
                        / "clean_candidate_manifest.json"
                    ),
                    clean_candidate_manifest_sha256="c" * 64,
                    expected_source_commit="b" * 40,
                    installed_control_plane_version="a" * 64,
                    q043_registered_matrix=(
                        prep.AUTHORIZED_ORION_ROOT
                        / "analysis/q043/matrix.json"
                    ),
                )

    def test_final_binding_output_is_exclusive_and_read_only(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bindings.json"
            value = {"record_type": "test", "schema_version": 1}
            prep._write_json_exclusive(path, value)
            self.assertEqual(json.loads(path.read_text(encoding="utf-8")), value)
            self.assertEqual(path.stat().st_mode & 0o777, 0o444)
            with self.assertRaises(FileExistsError):
                prep._write_json_exclusive(path, value)

    def test_timeout_artifact_is_strict_and_fresh(self) -> None:
        final = {
            key: value
            for key, value in prep.PENDING_FINAL_BINDINGS.items()
        }
        final["environment_profile_sha256"] = "a" * 64
        now = datetime(2026, 6, 8, 12, 0, tzinfo=timezone.utc)
        with patch.object(prep, "validate_final_bindings", return_value=final):
            value = prep.materialize_q023_timeout_margin(
                final_bindings=final,
                measured_utc="2026-06-08T11:59:00Z",
                expires_utc="2026-06-08T12:59:00Z",
            )
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "timeout.json"
                path.write_text(json.dumps(value), encoding="utf-8")
                path.chmod(0o444)
                self.assertEqual(
                    prep.validate_q023_timeout_margin_artifact(
                        path, final_bindings=final, now=now
                    ),
                    str(path),
                )
                with self.assertRaisesRegex(prep.PreparationError, "stale"):
                    prep.validate_q023_timeout_margin_artifact(
                        path,
                        final_bindings=final,
                        now=now + timedelta(hours=2),
                    )


if __name__ == "__main__":
    unittest.main()
