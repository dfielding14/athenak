#!/usr/bin/env python3
"""Focused fail-closed tests for the Q011 resource-scaling pilot materializer."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import stat
import sys
import tempfile
import unittest
import uuid

from tst.publication import q011_resource_scaling_preproduction_pilot_materializer as materializer

CONTROL_PLANE = Path(__file__).resolve().parent / "frontier_control_plane"
sys.path.insert(0, str(CONTROL_PLANE))
import control_plane_common  # noqa: E402


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json(payload: bytes) -> dict[str, object]:
    value = json.loads(payload.decode("utf-8"))
    assert isinstance(value, dict)
    return value


class Q011ResourceScalingPilotMaterializerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest, cls.files = materializer.build_materialization()
        cls.configs = [
            _json(cls.files[record["path"]])
            for record in cls.manifest["config_fragments"]
        ]

    def test_source_bindings_and_selected_pressure_are_exact(self) -> None:
        self.assertEqual(
            self.manifest["qualification_effect"], materializer.QUALIFICATION_EFFECT
        )
        self.assertEqual(
            self.manifest["authorization_state"], materializer.AUTHORIZATION_STATE
        )
        self.assertEqual(
            self.manifest["selected_pressure"],
            {"case_id": "ps_p0_1p00", "problem_ps_p0": 1.0},
        )
        sources = self.manifest["source_bindings"]
        self.assertEqual(
            sources["resource_plan"]["sha256"], materializer.RESOURCE_PLAN_SHA256
        )
        self.assertEqual(
            sources["selected_pressure_receipt"]["sha256"],
            materializer.SELECTED_PRESSURE_RECEIPT_SHA256,
        )
        self.assertEqual(
            sources["production_deck"]["sha256"],
            materializer.PRODUCTION_DECK_SHA256,
        )

    def test_versioned_readiness_record_binds_materializer_and_gates(self) -> None:
        path = (
            materializer.READINESS_ROOT
            / "q011_resource_scaling_preproduction_pilot_materializer_2026-06-06.json"
        )
        readiness = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(
            readiness["status"], "source_local_materializer_ready_non_authorizing"
        )
        bindings = readiness["source_bindings"]
        for name in ["resource_scaling_plan", "production_deck", "materializer_source"]:
            binding = bindings[name]
            self.assertEqual(
                _sha256((materializer.REPO_ROOT / binding["path"]).read_bytes()),
                binding["sha256"],
            )
        self.assertEqual(
            bindings["selected_pressure_receipt"]["sha256"],
            materializer.SELECTED_PRESSURE_RECEIPT_SHA256,
        )
        compatibility = readiness["registered_control_plane_compatibility"]
        self.assertEqual(
            compatibility["phase_1_concrete_launch_contract_schema_validation"],
            "passed",
        )
        self.assertFalse(compatibility["planner_retention_materialized_by_this_artifact"])
        self.assertFalse(readiness["execution_boundary"]["scheduler_submission_authorized"])
        self.assertEqual(
            readiness["materialized_output_contract"]["deterministic_manifest_sha256"],
            _sha256(materializer._json_bytes(self.manifest)),
        )

    def test_complete_matrix_and_bytes_are_deterministic(self) -> None:
        self.assertEqual(self.manifest["deck_count"], 12)
        self.assertEqual(self.manifest["config_fragment_count"], 24)
        self.assertEqual(len(self.manifest["decks"]), 12)
        self.assertEqual(len(self.configs), 24)
        second_manifest, second_files = materializer.build_materialization()
        self.assertEqual(second_manifest, self.manifest)
        self.assertEqual(second_files, self.files)
        for record in [*self.manifest["decks"], *self.manifest["config_fragments"]]:
            self.assertEqual(_sha256(self.files[record["path"]]), record["sha256"])

    def test_phase_1_ladders_and_conditional_repeats_are_exact(self) -> None:
        ladder = [value for value in self.configs if value["case_id"].startswith("p1-ladder")]
        repeats = [value for value in self.configs if value["case_id"].startswith("p1-repeat")]
        self.assertEqual(len(ladder), 9)
        self.assertEqual(len(repeats), 9)
        observed: dict[str, list[int]] = {}
        for value in ladder:
            observed.setdefault(value["variant"], []).append(value["node_binding"]["nodes"])
            self.assertEqual(value["engineering_seed"], 24060601)
            self.assertEqual(
                value["registered_policy_slice_fragment"]["maximum_attempts"], 1
            )
        self.assertEqual(
            {key: sorted(values) for key, values in observed.items()},
            {key: list(values) for key, values in materializer.NODE_LADDERS.items()},
        )
        for value in repeats:
            nodes = value["node_binding"]["nodes"]
            self.assertEqual(value["engineering_seed"], 24060602)
            self.assertEqual(
                value["execution_condition"]["run_only_if_phase_1_selected_node_count_equals"],
                nodes,
            )
            self.assertEqual(
                value["registered_policy_slice_fragment"]["maximum_attempts"], 3
            )

    def test_concrete_phase_1_launch_contracts_match_control_plane_schema(self) -> None:
        concrete = [
            value for value in self.configs if value["launch_contract_candidate"] is not None
        ]
        self.assertEqual(len(concrete), 18)
        for value in concrete:
            contract = value["launch_contract_candidate"]
            self.assertEqual(control_plane_common.validate_launch_contract(contract), contract)
            resources = contract["actions"][0]["resources"]
            self.assertEqual(resources["tasks"], resources["nodes"] * 8)
            expected = _sha256(
                json.dumps(
                    contract, sort_keys=True, separators=(",", ":"), allow_nan=False
                ).encode("utf-8")
            )
            self.assertEqual(
                value["registered_policy_slice_fragment"]["launch_contract_sha256"],
                expected,
            )

    def test_phase_2_and_phase_3_fail_closed_on_unresolved_node_selection(self) -> None:
        blocked = [
            value
            for value in self.configs
            if value["phase"]
            in {
                "phase_2_reduced_transverse_full_time_calibration",
                "phase_3_held_out_cost_model_validation",
            }
        ]
        self.assertEqual(len(blocked), 6)
        for value in blocked:
            self.assertIsNone(value["node_binding"]["nodes"])
            self.assertIsNone(value["launch_contract_candidate"])
            self.assertIn(
                "maximum_nodes",
                value["registered_policy_slice_fragment"]["unresolved_policy_bindings"],
            )
            self.assertEqual(
                value["launch_contract_template_gate"],
                "materialize_and_review_only_after_phase_1_node_selection",
            )

    def test_all_fragments_are_explicitly_non_authorizing_registered_science(self) -> None:
        engineering = set()
        for value in self.configs:
            engineering.add(value["engineering_seed"])
            self.assertEqual(value["submission_scope"], "registered_science")
            self.assertEqual(value["selected_qos"], "normal")
            self.assertFalse(value["qualifying_seed"])
            self.assertEqual(value["physical_output_inspection"], "forbidden")
            self.assertEqual(
                value["authorization_state"], materializer.AUTHORIZATION_STATE
            )
            self.assertFalse(value["execution_boundary"]["launch_authorized"])
            self.assertFalse(value["execution_boundary"]["frontier_execution_authorized"])
            self.assertFalse(
                value["execution_boundary"]["physical_output_inspection_authorized"]
            )
            self.assertTrue(
                value["planner_retention"]["required_by_installed_q011_registered_control_plane"]
            )
            self.assertFalse(value["planner_retention"]["materialized_by_this_artifact"])
            self.assertIn("qualifying_campaign_outputs", value["forbidden_inspection"])
        self.assertEqual(engineering, set(materializer.ENGINEERING_SEEDS))
        self.assertTrue(engineering.isdisjoint(materializer.QUALIFYING_SEEDS))

    def test_decks_bind_pressure_seeds_geometry_time_and_output_policy(self) -> None:
        for record in self.manifest["decks"]:
            blocks = materializer._parse_deck(self.files[record["path"]].decode("utf-8"))
            phase = record["phase"]
            variant = record["variant"]
            self.assertEqual(blocks["problem"]["ps_p0"], "1.0")
            expected_seed = str(materializer.PHASE_SEEDS[phase])
            self.assertEqual(blocks["particles"]["pic_random_seed"], expected_seed)
            self.assertEqual(blocks["problem"]["ps_inject_seed"], expected_seed)
            self.assertEqual(blocks["problem"]["ps_seed_noise_seed"], expected_seed)
            if phase == "phase_2_reduced_transverse":
                self.assertEqual(blocks["mesh"]["x2max"], "240.0")
                self.assertEqual(
                    blocks["mesh"]["nx2"],
                    "80" if variant == "fine_uniform_dx3" else "20",
                )
                self.assertEqual(blocks["time"]["tlim"], "1200.0")
                for index in range(1, 7):
                    self.assertEqual(blocks[f"output{index}"]["dt"], "100.0")
            else:
                self.assertEqual(blocks["mesh"]["x2max"], "3120.0")
                self.assertEqual(blocks["time"]["nlim"], "512")
                for index in range(1, 7):
                    self.assertLessEqual(float(blocks[f"output{index}"]["dt"]), 0.0)

    def test_drifted_plan_receipt_and_deck_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan = root / "plan.json"
            plan.write_bytes(materializer.RESOURCE_PLAN.read_bytes() + b"\n")
            with self.assertRaisesRegex(materializer.MaterializationError, "plan SHA-256"):
                materializer.build_materialization(resource_plan=plan)

            receipt = root / "receipt.json"
            value = json.loads(materializer.SELECTED_PRESSURE_RECEIPT.read_text(encoding="utf-8"))
            value["selected_case"]["problem_ps_p0"] = 0.2
            receipt.write_text(json.dumps(value), encoding="utf-8")
            with self.assertRaisesRegex(materializer.MaterializationError, "receipt SHA-256"):
                materializer.build_materialization(selected_pressure_receipt=receipt)

            deck = root / "deck.athinput"
            deck.write_bytes(materializer.PRODUCTION_DECK.read_bytes() + b"\n")
            with self.assertRaisesRegex(materializer.MaterializationError, "deck SHA-256"):
                materializer.build_materialization(production_deck=deck)

    def test_materialization_is_safe_complete_read_only_and_no_overwrite(self) -> None:
        output = materializer.AUTHORIZED_OUTPUT_PARENT / f"q011-rs-test-{uuid.uuid4().hex}"
        try:
            manifest = materializer.materialize_pilot_artifacts(output)
            self.assertEqual(manifest, self.manifest)
            expected_count = len(self.files) + 2
            written = [path for path in output.rglob("*") if path.is_file()]
            self.assertEqual(len(written), expected_count)
            self.assertTrue(all(not path.stat().st_mode & 0o222 for path in [output, *output.rglob("*")]))
            with self.assertRaisesRegex(materializer.MaterializationError, "already exists"):
                materializer.materialize_pilot_artifacts(output)
        finally:
            if output.exists():
                for path in [output, *output.rglob("*")]:
                    path.chmod(path.stat().st_mode | stat.S_IWUSR)
                shutil.rmtree(output)

    def test_output_outside_direct_codex_child_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(materializer.MaterializationError, "direct tst/.codex child"):
                materializer.materialize_pilot_artifacts(Path(directory) / "not-authorized")


if __name__ == "__main__":
    unittest.main()
