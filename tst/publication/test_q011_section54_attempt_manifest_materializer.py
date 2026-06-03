#!/usr/bin/env python3
"""Planner-to-attempt integration coverage for the Q-011 Section 5.4 bridge."""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from tst.publication import analyze_q011_section54_campaign as campaign
from tst.publication import publish_q011_section54_campaign_attempt as attempt_publisher
from tst.publication import q011_section54_attempt_manifest_materializer as bridge
from tst.publication import q011_section54_qualifying_campaign_execution as planner
from tst.publication import test_analyze_q011_section54_campaign as raw_fixtures
from tst.publication import (
    test_q011_section54_qualifying_campaign_execution as planner_fixtures,
)


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_completed_raw_outputs(root: Path) -> None:
    for cycle, time in enumerate(raw_fixtures._TIMES):
        suffix = f"{int(time):05d}"
        for kind in ("rho", "bmag", "prtcl_jx", "j2"):
            raw_fixtures._product(
                root,
                kind,
                f"bin/q011.{kind}.{suffix}.bin",
                raw_fixtures._mesh_bin(kind, time, cycle),
                time,
            )
        raw_fixtures._product(
            root,
            "prtcl_all",
            f"pvtk/q011.prtcl_all.{suffix}.part.vtk",
            raw_fixtures._particle_vtk(time, cycle),
            time,
        )
        products: list[dict[str, object]] = []
        raw_fixtures._append_restart_publication(root, products, time, cycle)
    raw_fixtures._product(
        root, "stdout", "stdout.txt", raw_fixtures._stdout_telemetry(), None
    )


def _write_registered_execution_receipt(
    attempt_root: Path,
    artifact_dir: Path,
    planner_root: Path,
    planner_inventory_sha256: str,
    planner_result: dict[str, object],
    plan: dict[str, object],
    contract: dict[str, object],
) -> Path:
    candidate = plan["candidate_binding"]
    source_bindings = plan["source_bindings"]
    assert isinstance(candidate, dict)
    assert isinstance(source_bindings, dict)
    environment = candidate["environment_profile"]
    assert isinstance(environment, dict)
    receipt = raw_fixtures._registered_execution_receipt(
        attempt_id=str(contract["attempt_id"]),
        source_commit=str(candidate["git_commit"]),
        executable_sha256=str(candidate["executable"]["sha256"]),
        deck_sha256=str(source_bindings["paper_deck"]["sha256"]),
        environment_sha256=str(environment["sha256"]),
        control_plane_version=str(environment["control_plane_version"]),
        argv=list(contract["argv"]),
        raw_output_root=f"{attempt_root}/raw",
        planner_result=planner_result,
        artifact_dir=str(artifact_dir),
    )
    receipt["planner_retention"] = planner.materialize_planner_retention(
        planner_root=planner_root,
        planner_inventory_sha256=planner_inventory_sha256,
        attempt_id=str(contract["attempt_id"]),
        authorized_pic_root=artifact_dir.parents[2],
    )
    path = artifact_dir / "analysis" / campaign.REGISTERED_EXECUTION_RECEIPT_NAME
    path.parent.mkdir()
    path.write_bytes(raw_fixtures._json_bytes(receipt))
    path.chmod(0o444)
    return path


@contextmanager
def _prepared_real_attempt() -> object:
    with planner_fixtures._fixture() as fixture:
        planner_parent = fixture["orion"] / "plans"
        planner_parent.mkdir()
        planner_result = planner_fixtures._materialize(
            fixture, output_parent=planner_parent
        )
        planner_root = Path(planner_result["plan_root"])
        plan = _json(planner_root / "campaign_plan.json")
        descriptor = _json(
            planner_root / plan["baseline_attempt_descriptors"][8]["path"]
        )
        attempt_root = Path(descriptor["authorized_orion_attempt_root"])
        artifact_dir = (
            fixture["orion"]
            / "runs"
            / "q011_section54_qualifying"
            / "22222222-2222-4222-8222-222222222222"
        )
        raw_root = artifact_dir / "raw"
        raw_root.mkdir(parents=True)
        _write_completed_raw_outputs(raw_root)
        contract = _json(planner_root / descriptor["launch_contract"]["path"])
        receipt = _write_registered_execution_receipt(
            attempt_root,
            artifact_dir,
            planner_root,
            str(planner_result["inventory_sha256"]),
            planner_result,
            plan,
            contract,
        )
        try:
            yield {
                "fixture": fixture,
                "planner_parent": planner_parent,
                "planner_result": planner_result,
                "planner_root": planner_root,
                "plan": plan,
                "descriptor": descriptor,
                "attempt_root": attempt_root,
                "artifact_dir": artifact_dir,
                "raw_root": raw_root,
                "registered_execution_receipt": receipt,
            }
        finally:
            raw_fixtures._make_writable_tree(planner_root)
            if planner_root.exists():
                shutil.rmtree(planner_root)
            if planner_parent.exists():
                planner_parent.rmdir()
            if artifact_dir.exists():
                shutil.rmtree(artifact_dir)


@contextmanager
def _fixture_campaign_context(fixture: dict[str, Path]) -> object:
    def validate_fixture_ledger(
        receipt: dict[str, object], *, authorized_pic_root: Path, **_kwargs: object
    ) -> dict[str, str]:
        if authorized_pic_root != fixture["orion"]:
            raise campaign.QualificationError(
                "fixture registered-execution ledger root drifted",
                code="registered_execution_receipt_ledger_drift",
            )
        return {
            "reservation_id": str(receipt["reservation_id"]),
            "submission_id": str(receipt["submission_id"]),
            "reconciliation_event_sha256": str(
                receipt["reconciliation_event_sha256"]
            ),
        }

    @contextmanager
    def validate_fixture_ledger_snapshot(
        receipt: dict[str, object], *, authorized_pic_root: Path, **kwargs: object
    ) -> object:
        yield validate_fixture_ledger(
            receipt, authorized_pic_root=authorized_pic_root, **kwargs
        )

    with (
        patch.object(campaign, "ORION_BULK_ROOT", fixture["orion"]),
        patch.object(
            campaign,
            "_validate_registered_execution_receipt_ledger_binding",
            side_effect=validate_fixture_ledger,
        ),
        patch.object(
            campaign,
            "_validated_registered_execution_receipt_ledger_snapshot",
            side_effect=validate_fixture_ledger_snapshot,
        ),
    ):
        yield


class Q011Section54AttemptManifestMaterializerTests(unittest.TestCase):
    def test_real_planner_tree_materializes_authorized_raw_attempt(self) -> None:
        with _prepared_real_attempt() as prepared:
            fixture = prepared["fixture"]
            with _fixture_campaign_context(fixture):
                bridged = bridge.materialize_completed_attempt_manifest(
                    planner_root=prepared["planner_root"],
                    planner_inventory_sha256=prepared["planner_result"][
                        "inventory_sha256"
                    ],
                    raw_output_root=prepared["raw_root"],
                    registered_execution_receipt=prepared[
                        "registered_execution_receipt"
                    ],
                    attempt_id=prepared["descriptor"]["attempt_id"],
                    authorized_pic_root=fixture["orion"],
                )
            manifest = _json(prepared["raw_root"] / campaign.MANIFEST_NAME)
            self.assertFalse(bridged["launch_authorized"])
            self.assertFalse(bridged["scheduler_submission_authorized"])
            self.assertFalse(bridged["claim_closure_authorized"])
            self.assertEqual(bridged["product_count"], 118)
            self.assertIn(
                "registered_execution_receipt", manifest["artifact_bindings"]
            )
            self.assertFalse(prepared["attempt_root"].exists())
            self.assertEqual(
                _json(prepared["registered_execution_receipt"])["artifact_dir"],
                str(prepared["artifact_dir"]),
            )

    def test_separated_namespace_materializes_and_publishes_end_to_end(self) -> None:
        with _prepared_real_attempt() as prepared:
            fixture = prepared["fixture"]
            attempt_root = prepared["attempt_root"]
            attempt_root.parent.mkdir(parents=True)
            try:
                with _fixture_campaign_context(fixture):
                    bridge.materialize_completed_attempt_manifest(
                        planner_root=prepared["planner_root"],
                        planner_inventory_sha256=prepared["planner_result"][
                            "inventory_sha256"
                        ],
                        raw_output_root=prepared["raw_root"],
                        registered_execution_receipt=prepared[
                            "registered_execution_receipt"
                        ],
                        attempt_id=prepared["descriptor"]["attempt_id"],
                        authorized_pic_root=fixture["orion"],
                    )
                    with patch.object(
                        campaign,
                        "_validate_external_clean_candidate_closure",
                        return_value={"validation": "fixture_external_candidate_closure"},
                    ):
                        retained = attempt_publisher.freeze_campaign_attempt(
                            prepared["raw_root"],
                            attempt_root.parent,
                            attempt_status="completed",
                            authorized_pic_root=fixture["orion"],
                        )
                self.assertEqual(Path(retained["campaign_root"]), attempt_root)
                self.assertTrue((attempt_root / campaign.MANIFEST_NAME).is_file())
                self.assertTrue(prepared["raw_root"].is_dir())
            finally:
                if attempt_root.exists():
                    raw_fixtures._make_writable_tree(attempt_root)
                    shutil.rmtree(attempt_root)

    def test_rejects_arbitrarily_relabeled_raw_output_root(self) -> None:
        with _prepared_real_attempt() as prepared:
            fixture = prepared["fixture"]
            relabeled_raw = fixture["root"] / "arbitrarily-relabeled-raw"
            relabeled_raw.mkdir()
            _write_completed_raw_outputs(relabeled_raw)
            with (
                _fixture_campaign_context(fixture),
                self.assertRaisesRegex(
                    bridge.AttemptManifestMaterializationError,
                    "must be runs/<campaign>/<submission-id>/raw",
                ),
            ):
                bridge.materialize_completed_attempt_manifest(
                    planner_root=prepared["planner_root"],
                    planner_inventory_sha256=prepared["planner_result"][
                        "inventory_sha256"
                    ],
                    raw_output_root=relabeled_raw,
                    registered_execution_receipt=prepared[
                        "registered_execution_receipt"
                    ],
                    attempt_id=prepared["descriptor"]["attempt_id"],
                    authorized_pic_root=fixture["orion"],
                )
            self.assertFalse((relabeled_raw / "bindings").exists())
            self.assertFalse((relabeled_raw / campaign.MANIFEST_NAME).exists())

    def test_write_failure_rolls_back_and_retry_succeeds(self) -> None:
        with _prepared_real_attempt() as prepared:
            fixture = prepared["fixture"]
            raw_root = prepared["raw_root"]
            original_write_member = bridge._write_member
            write_count = 0

            def fail_after_write(
                root_fd: int, relative: str, payload: bytes, mode: int = 0o644
            ) -> None:
                nonlocal write_count
                original_write_member(root_fd, relative, payload, mode)
                write_count += 1
                if write_count == 3:
                    raise bridge.AttemptManifestMaterializationError(
                        "injected staged write failure"
                    )

            arguments = {
                "planner_root": prepared["planner_root"],
                "planner_inventory_sha256": prepared["planner_result"][
                    "inventory_sha256"
                ],
                "raw_output_root": raw_root,
                "registered_execution_receipt": prepared[
                    "registered_execution_receipt"
                ],
                "attempt_id": prepared["descriptor"]["attempt_id"],
                "authorized_pic_root": fixture["orion"],
            }
            with (
                _fixture_campaign_context(fixture),
                patch.object(bridge, "_write_member", side_effect=fail_after_write),
                self.assertRaisesRegex(
                    bridge.AttemptManifestMaterializationError,
                    "injected staged write failure",
                ),
            ):
                bridge.materialize_completed_attempt_manifest(**arguments)
            self.assertFalse((raw_root / "bindings").exists())
            self.assertFalse((raw_root / campaign.MANIFEST_NAME).exists())
            self.assertFalse(
                any(
                    path.name.startswith(".q011-attempt-manifest-staging-")
                    for path in raw_root.iterdir()
                )
            )
            with _fixture_campaign_context(fixture):
                bridged = bridge.materialize_completed_attempt_manifest(**arguments)
            self.assertEqual(bridged["product_count"], 118)

    def test_raw_root_replacement_during_staging_is_rejected(self) -> None:
        with _prepared_real_attempt() as prepared:
            fixture = prepared["fixture"]
            raw_root = prepared["raw_root"]
            detached = raw_root.with_name("raw.detached")
            original_write = bridge._write_member
            substituted = False

            def substitute_then_write(*args: object, **kwargs: object) -> None:
                nonlocal substituted
                if not substituted:
                    raw_root.rename(detached)
                    raw_root.mkdir()
                    (raw_root / "sentinel.txt").write_text(
                        "replacement raw root\n", encoding="utf-8"
                    )
                    substituted = True
                original_write(*args, **kwargs)

            with (
                _fixture_campaign_context(fixture),
                patch.object(bridge, "_write_member", side_effect=substitute_then_write),
                self.assertRaisesRegex(
                    bridge.AttemptManifestMaterializationError,
                    "completed raw output root public pathname changed",
                ),
            ):
                bridge.materialize_completed_attempt_manifest(
                    planner_root=prepared["planner_root"],
                    planner_inventory_sha256=prepared["planner_result"][
                        "inventory_sha256"
                    ],
                    raw_output_root=raw_root,
                    registered_execution_receipt=prepared[
                        "registered_execution_receipt"
                    ],
                    attempt_id=prepared["descriptor"]["attempt_id"],
                    authorized_pic_root=fixture["orion"],
                )
            self.assertEqual(
                (raw_root / "sentinel.txt").read_text(encoding="utf-8"),
                "replacement raw root\n",
            )
            self.assertFalse((raw_root / "bindings").exists())
            self.assertFalse((raw_root / campaign.MANIFEST_NAME).exists())
            self.assertFalse((detached / "bindings").exists())
            self.assertFalse((detached / campaign.MANIFEST_NAME).exists())
            self.assertFalse(
                any(
                    path.name.startswith(".q011-attempt-manifest-staging-")
                    for path in detached.iterdir()
                )
            )

    def test_minimal_self_authored_planner_is_rejected_before_write(self) -> None:
        with _prepared_real_attempt() as prepared:
            fixture = prepared["fixture"]
            holder = fixture["orion"] / "minimal-holder"
            holder.mkdir()
            planner_receipt = raw_fixtures._minimal_planner_materialization_receipt(
                holder,
                (prepared["planner_root"] / "campaign_plan.json").read_bytes(),
                (
                    prepared["planner_root"] / "helper_source_closure.json"
                ).read_bytes(),
            )
            outer_receipt = _json(holder / planner_receipt["path"])
            try:
                with (
                    _fixture_campaign_context(fixture),
                    patch.object(bridge, "_write_member") as write_member,
                    self.assertRaisesRegex(
                        bridge.AttemptManifestMaterializationError,
                        "unavailable in immutable planner tree|"
                        "production planner graph validation failed",
                    ),
                ):
                    bridge.materialize_completed_attempt_manifest(
                        planner_root=outer_receipt["plan_root"],
                        planner_inventory_sha256=outer_receipt["inventory_sha256"],
                        raw_output_root=prepared["raw_root"],
                        registered_execution_receipt=prepared[
                            "registered_execution_receipt"
                        ],
                        attempt_id=prepared["descriptor"]["attempt_id"],
                        authorized_pic_root=fixture["orion"],
                    )
                write_member.assert_not_called()
            finally:
                raw_fixtures._make_writable_tree(fixture["orion"] / "minimal-plans")

    def test_self_authored_execution_receipt_is_rejected_before_write(self) -> None:
        with _prepared_real_attempt() as prepared:
            fixture = prepared["fixture"]
            error = campaign.QualificationError(
                "fixture receipt lacks mirrored ledger authority",
                code="registered_execution_receipt_ledger_drift",
            )
            with (
                patch.object(campaign, "ORION_BULK_ROOT", fixture["orion"]),
                patch.object(
                    campaign,
                    "_validate_registered_execution_receipt_ledger_binding",
                    side_effect=error,
                ),
                patch.object(bridge, "_write_member") as write_member,
                self.assertRaisesRegex(
                    bridge.AttemptManifestMaterializationError,
                    "lacks mirrored ledger authority",
                ),
            ):
                bridge.materialize_completed_attempt_manifest(
                    planner_root=prepared["planner_root"],
                    planner_inventory_sha256=prepared["planner_result"][
                        "inventory_sha256"
                    ],
                    raw_output_root=prepared["raw_root"],
                    registered_execution_receipt=prepared[
                        "registered_execution_receipt"
                    ],
                    attempt_id=prepared["descriptor"]["attempt_id"],
                    authorized_pic_root=fixture["orion"],
                )
            write_member.assert_not_called()

    def test_rejects_unrecognized_completed_raw_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            raw_root = Path(directory)
            (raw_root / "stdout.txt").write_text("fixture\n", encoding="utf-8")
            (raw_root / "unexpected.dat").write_text("unsupported\n", encoding="utf-8")
            descriptor = 0
            snapshot = None
            try:
                descriptor = __import__("os").open(raw_root, bridge._DIRECTORY_FLAGS)
                snapshot = bridge._raw_snapshot(descriptor)
                with self.assertRaisesRegex(
                    bridge.AttemptManifestMaterializationError,
                    "unsupported files",
                ):
                    bridge._classify_products(descriptor, snapshot)
            finally:
                if descriptor:
                    __import__("os").close(descriptor)

    def test_stable_bound_source_payload_rejects_checksum_drift(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "analyzer.py"
            source.write_text("print('bound')\n", encoding="utf-8")
            expected = hashlib.sha256(source.read_bytes()).hexdigest()
            self.assertEqual(
                bridge._stable_bound_source_payload(
                    source, expected_sha256=expected, label="fixture analyzer"
                ),
                b"print('bound')\n",
            )
            source.write_text("print('drifted')\n", encoding="utf-8")
            with self.assertRaisesRegex(
                bridge.AttemptManifestMaterializationError,
                "fixture analyzer SHA-256 drifted",
            ):
                bridge._stable_bound_source_payload(
                    source, expected_sha256=expected, label="fixture analyzer"
                )


if __name__ == "__main__":
    unittest.main()
