#!/usr/bin/env python3
"""Adversarial tests for Q-011 Section 5.4 final-evidence orchestration."""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path
import stat
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent / "frontier_control_plane"))

from . import analyze_q011_section54_numerical_qualification as numerical
from . import q011_section54_final_evidence_orchestrator as orchestrator
from . import test_analyze_q011_section54_numerical_qualification as numerical_fixture


def _sha256(index: int) -> str:
    return f"{index:064x}"


def _pending() -> dict[str, object]:
    return {
        "status": "pending_external_review",
        "reviewer": None,
        "reviewed_at_utc": None,
        "notes": "Terminal external review remains deliberately open.",
    }


def _accepted() -> dict[str, object]:
    return {
        "status": "accepted",
        "reviewer": "independent-reviewer",
        "reviewed_at_utc": "2026-06-06T12:00:00Z",
        "notes": "Independent raw-artifact recomputation reviewed and accepted.",
    }


def _production_restart_parity(aggregate: dict[str, object]) -> dict[str, object]:
    source = next(
        attempt["attempt_sha256"]
        for attempt in aggregate["attempt_gate_results"]
        if attempt["identity"]["variant"] == numerical.AMR_VARIANT
        and attempt["identity"]["seed"] == numerical.QUALIFYING_SEEDS[0]
    )
    return {
        "source_attempt_sha256": source,
        "restart_parity_binding_sha256": _sha256(900),
        "screen_scope": numerical.RESTART_SCREEN_SCOPE,
        "full_state_equivalence_claimed": False,
        "result": {
            "result": "pass_deterministic_continuation_parity",
            "checkpoint_nominal_slot_omega0_inverse": 500.0,
            "checkpoint_observed_committed_cycle": 5000,
            "checkpoint_observed_committed_time_omega0_inverse": 500.01,
            "retained_output_nominal_slots_after_checkpoint_omega0_inverse": list(
                orchestrator.RESTART_OUTPUT_SLOTS
            ),
            "paired_output_observed_commits": [
                {
                    "nominal_slot_omega0_inverse": slot,
                    "observed_committed_cycle": 5100 + index,
                    "observed_committed_time_omega0_inverse": slot + 0.01,
                }
                for index, slot in enumerate(orchestrator.RESTART_OUTPUT_SLOTS)
            ],
            "maximum_absolute_difference_by_field": {
                field: 0.0 for field in orchestrator.RESTART_FIELDS
            },
        },
    }


def _aggregate() -> dict[str, object]:
    aggregate = numerical_fixture._qualify(numerical_fixture._attempts())
    aggregate["restart_parity"] = _production_restart_parity(aggregate)
    return aggregate


def _bundles(aggregate: dict[str, object]) -> list[dict[str, object]]:
    return [
        {
            "attempt_id": attempt["identity"]["attempt_id"],
            "variant": attempt["identity"]["variant"],
            "seed": attempt["identity"]["seed"],
            "attempt_sha256": attempt["attempt_sha256"],
            "raw_inventory_sha256": attempt["raw_inventory_sha256"],
            "admission_result_sha256": _sha256(1000 + index),
            "admitted_for_follow_on_numerical_qualification": True,
        }
        for index, attempt in enumerate(aggregate["attempt_gate_results"])
    ]


def _metric_row(
    *,
    attempt_id: str,
    variant: str,
    seed: int,
    observable: str,
    production_metric: float,
) -> dict[str, object]:
    return {
        "attempt_id": attempt_id,
        "variant": variant,
        "qualifying_seed": seed,
        "observable": observable,
        "production_metric": production_metric,
        "independent_metric": production_metric,
        "absolute_difference": 0.0,
        "relative_difference": None if production_metric == 0.0 else 0.0,
        "declared_tolerance": 0.0,
        "disposition": "pass_within_declared_tolerance",
    }


def _metric_rows(aggregate: dict[str, object]) -> list[dict[str, object]]:
    rows = []
    by_cell = {
        (attempt["identity"]["variant"], attempt["identity"]["seed"]): attempt
        for attempt in aggregate["attempt_gate_results"]
    }
    for variant, seed in numerical.EXPECTED_MATRIX_CELLS:
        identity = by_cell[(variant, seed)]["identity"]
        for observable in orchestrator.ATTEMPT_PRIMARY_OBSERVABLES:
            rows.append(
                _metric_row(
                    attempt_id=identity["attempt_id"],
                    variant=variant,
                    seed=seed,
                    observable=observable,
                    production_metric=1.0,
                )
            )
    for pair in aggregate["paired_amr_fine_results"]:
        identity = by_cell[(numerical.AMR_VARIANT, pair["seed"])]["identity"]
        for residual in pair["residuals"]:
            rows.append(
                _metric_row(
                    attempt_id=identity["attempt_id"],
                    variant=numerical.AMR_VARIANT,
                    seed=pair["seed"],
                    observable=(
                        f"{orchestrator.PAIR_OBSERVABLE_PREFIX}{residual['observable']}"
                    ),
                    production_metric=residual["maximum_absolute_difference"],
                )
            )
    restart = aggregate["restart_parity"]
    source = next(
        attempt
        for attempt in aggregate["attempt_gate_results"]
        if attempt["attempt_sha256"] == restart["source_attempt_sha256"]
    )
    maximum = max(
        restart["result"]["maximum_absolute_difference_by_field"].values()
    )
    rows.append(
        _metric_row(
            attempt_id=source["identity"]["attempt_id"],
            variant=source["identity"]["variant"],
            seed=source["identity"]["seed"],
            observable=orchestrator.RESTART_OBSERVABLE,
            production_metric=maximum,
        )
    )
    return rows


def _independent_recompute(aggregate: dict[str, object]) -> dict[str, object]:
    return {
        "schema_version": 1,
        "record_type": orchestrator.INDEPENDENT_RECOMPUTE_RECORD_TYPE,
        "campaign_id": numerical.admission.CAMPAIGN_ID,
        "status": "pass_independent_recompute_reviewed",
        "production_aggregate_sha256": orchestrator._canonical_sha256(aggregate),
        "independent_script": {
            "path_or_archive_locator": "/archive/independent/q011_recompute.py",
            "sha256": _sha256(2001),
        },
        "independent_environment_lock": {
            "path_or_archive_locator": "/archive/independent/environment.lock",
            "sha256": _sha256(2002),
        },
        "production_helper_imports_authorized": False,
        "production_helper_imports_detected": False,
        "attempt_inventory_sha256_values": list(
            aggregate["ordered_raw_inventory_sha256_values"]
        ),
        "input_artifact_checksums": [
            {
                "path_or_archive_locator": f"/archive/attempt-{index:03d}/artifact_inventory",
                "sha256": inventory,
            }
            for index, inventory in enumerate(
                aggregate["ordered_raw_inventory_sha256_values"]
            )
        ],
        "metric_comparison_table": _metric_rows(aggregate),
        "reviewer_identity": "independent-reviewer",
        "reviewer_disposition": _accepted(),
    }


def _build(
    *,
    aggregate: dict[str, object] | None = None,
    bundles: list[dict[str, object]] | None = None,
    recompute: dict[str, object] | None = None,
    external_review: dict[str, object] | None = None,
) -> dict[str, object]:
    aggregate = aggregate or _aggregate()
    return orchestrator.build_claim_candidate_evidence_manifest(
        admitted_bundle_identities=bundles or _bundles(aggregate),
        numerical_aggregate=aggregate,
        independent_recompute=recompute or _independent_recompute(aggregate),
        external_reviewer_disposition=external_review or _pending(),
    )


class Q011Section54FinalEvidenceOrchestratorTests(unittest.TestCase):
    def test_complete_evidence_emits_pending_non_authorizing_claim_candidate(self) -> None:
        manifest = _build()
        self.assertEqual(manifest["status"], "claim_candidate_pending_external_review")
        self.assertTrue(manifest["claim_candidate"])
        self.assertFalse(manifest["final_claim_closure"])
        self.assertFalse(manifest["claim_closure_authorized"])
        self.assertFalse(manifest["frontier_execution_authorized"])
        self.assertFalse(manifest["policy_mutation_authorized"])
        self.assertEqual(len(manifest["admitted_bundle_identities"]), 24)
        self.assertEqual(len(manifest["amr_fine_pair_coverage"]), 8)
        self.assertEqual(
            manifest["independent_recomputation"]["metric_comparison_row_count"],
            169,
        )
        self.assertEqual(
            manifest["external_reviewer_disposition"]["status"],
            "pending_external_review",
        )
        self.assertEqual(
            orchestrator.validate_claim_candidate_evidence_manifest(manifest),
            manifest,
        )

    def test_admitted_bundle_matrix_and_identity_drift_fail_closed(self) -> None:
        aggregate = _aggregate()
        bundles = _bundles(aggregate)
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "exactly 24"):
            _build(aggregate=aggregate, bundles=bundles[:-1])

        duplicate = copy.deepcopy(bundles)
        duplicate[1]["admission_result_sha256"] = duplicate[0]["admission_result_sha256"]
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "duplicate admission"):
            _build(aggregate=aggregate, bundles=duplicate)

        drift = copy.deepcopy(bundles)
        drift[0]["attempt_sha256"] = _sha256(9999)
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "attempt_sha256 binding"):
            _build(aggregate=aggregate, bundles=drift)

        not_admitted = copy.deepcopy(bundles)
        not_admitted[0]["admitted_for_follow_on_numerical_qualification"] = False
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "was not admitted"):
            _build(aggregate=aggregate, bundles=not_admitted)

        numeric_alias = copy.deepcopy(bundles)
        numeric_alias[0]["seed"] = float(numeric_alias[0]["seed"])
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "seed must be an integer"):
            _build(aggregate=aggregate, bundles=numeric_alias)

    def test_aggregate_numerical_and_pair_failures_fail_closed(self) -> None:
        aggregate = _aggregate()
        aggregate["numerical_gates_passed"] = False
        aggregate["numerical_gate_status"] = "failed"
        with self.assertRaisesRegex(
            orchestrator.FinalEvidenceError, "numerical gates did not pass"
        ):
            _build(aggregate=aggregate)

        aggregate = _aggregate()
        aggregate["paired_amr_fine_results"][0]["gates_passed"] = False
        with self.assertRaisesRegex(
            orchestrator.FinalEvidenceError, "paired numerical gates failed"
        ):
            _build(aggregate=aggregate)

        aggregate = _aggregate()
        aggregate["paired_amr_fine_results"][0]["fine_uniform_attempt_sha256"] = (
            aggregate["paired_amr_fine_results"][1]["fine_uniform_attempt_sha256"]
        )
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "fine-uniform pair identity"):
            _build(aggregate=aggregate)

    def test_restart_parity_must_be_retained_first_seed_production_screen(self) -> None:
        aggregate = _aggregate()
        aggregate["restart_parity"]["screen_scope"] = "unit_only_synthetic_restart_comparator"
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "retained production screen"):
            _build(aggregate=aggregate)

        aggregate = _aggregate()
        second_amr = next(
            attempt["attempt_sha256"]
            for attempt in aggregate["attempt_gate_results"]
            if attempt["identity"]["variant"] == numerical.AMR_VARIANT
            and attempt["identity"]["seed"] == numerical.QUALIFYING_SEEDS[1]
        )
        aggregate["restart_parity"]["source_attempt_sha256"] = second_amr
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "first-seed AMR"):
            _build(aggregate=aggregate)

        aggregate = _aggregate()
        aggregate["restart_parity"]["result"]["result"] = "failed"
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "did not pass"):
            _build(aggregate=aggregate)

        aggregate = _aggregate()
        aggregate["restart_parity"]["result"]["paired_output_observed_commits"][1][
            "observed_committed_cycle"
        ] = aggregate["restart_parity"]["result"]["paired_output_observed_commits"][0][
            "observed_committed_cycle"
        ]
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "sequence drifted"):
            _build(aggregate=aggregate)

    def test_independent_recompute_binding_and_coverage_fail_closed(self) -> None:
        aggregate = _aggregate()
        recompute = _independent_recompute(aggregate)
        recompute["production_aggregate_sha256"] = _sha256(9999)
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "aggregate binding drifted"):
            _build(aggregate=aggregate, recompute=recompute)

        recompute = _independent_recompute(aggregate)
        recompute["production_helper_imports_detected"] = True
        with self.assertRaisesRegex(
            orchestrator.FinalEvidenceError, "helper imports are forbidden"
        ):
            _build(aggregate=aggregate, recompute=recompute)

        recompute = _independent_recompute(aggregate)
        recompute["input_artifact_checksums"] = recompute["input_artifact_checksums"][:-1]
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "inventories are incomplete"):
            _build(aggregate=aggregate, recompute=recompute)

        recompute = _independent_recompute(aggregate)
        recompute["metric_comparison_table"] = recompute["metric_comparison_table"][:-1]
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "coverage count drifted"):
            _build(aggregate=aggregate, recompute=recompute)

        recompute = _independent_recompute(aggregate)
        recompute["metric_comparison_table"][0]["independent_metric"] = 2.0
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "absolute difference"):
            _build(aggregate=aggregate, recompute=recompute)

        recompute = _independent_recompute(aggregate)
        recompute["metric_comparison_table"][0]["qualifying_seed"] = float(
            recompute["metric_comparison_table"][0]["qualifying_seed"]
        )
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "seed must be an integer"):
            _build(aggregate=aggregate, recompute=recompute)

    def test_independent_pair_metric_and_reviewer_disposition_are_bound(self) -> None:
        aggregate = _aggregate()
        recompute = _independent_recompute(aggregate)
        pair_row = len(numerical.EXPECTED_MATRIX_CELLS) * len(
            orchestrator.ATTEMPT_PRIMARY_OBSERVABLES
        )
        recompute["metric_comparison_table"][pair_row]["production_metric"] = 1.0
        recompute["metric_comparison_table"][pair_row]["independent_metric"] = 1.0
        recompute["metric_comparison_table"][pair_row]["relative_difference"] = 0.0
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "aggregate metric binding"):
            _build(aggregate=aggregate, recompute=recompute)

        recompute = _independent_recompute(aggregate)
        recompute["reviewer_identity"] = "different-reviewer"
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "reviewer identity binding"):
            _build(aggregate=aggregate, recompute=recompute)

        recompute = _independent_recompute(aggregate)
        recompute["reviewer_disposition"] = _pending()
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "must remain accepted"):
            _build(aggregate=aggregate, recompute=recompute)

    def test_terminal_external_review_must_remain_pending(self) -> None:
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "pending_external_review"):
            _build(external_review=_accepted())

    def test_compact_manifest_validator_rejects_authority_and_status_tamper(self) -> None:
        manifest = _build()
        authorized = copy.deepcopy(manifest)
        authorized["claim_closure_authorized"] = True
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "must be false"):
            orchestrator.validate_claim_candidate_evidence_manifest(authorized)

        promoted = copy.deepcopy(manifest)
        promoted["status"] = "qualified"
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "candidate status drifted"):
            orchestrator.validate_claim_candidate_evidence_manifest(promoted)

        closed_gate = copy.deepcopy(manifest)
        closed_gate["claim_registry_boundary"]["gates_closed_by_this_manifest"] = ["Q-011"]
        with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "registry boundary drifted"):
            orchestrator.validate_claim_candidate_evidence_manifest(closed_gate)

    def test_cli_emits_canonical_read_only_manifest_and_refuses_overwrite(self) -> None:
        aggregate = _aggregate()
        inputs = {
            "bundles.json": _bundles(aggregate),
            "aggregate.json": aggregate,
            "recompute.json": _independent_recompute(aggregate),
            "external.json": _pending(),
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, value in inputs.items():
                (root / name).write_text(
                    json.dumps(value, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
            output = root / "claim-candidate.json"
            argv = [
                "--admitted-bundle-identities",
                str(root / "bundles.json"),
                "--numerical-aggregate",
                str(root / "aggregate.json"),
                "--independent-recompute",
                str(root / "recompute.json"),
                "--external-reviewer-disposition",
                str(root / "external.json"),
                "--output",
                str(output),
            ]
            self.assertEqual(orchestrator.main(argv), 0)
            emitted = json.loads(output.read_text(encoding="utf-8"))
            orchestrator.validate_claim_candidate_evidence_manifest(emitted)
            self.assertEqual(
                output.read_bytes(),
                orchestrator._canonical_json_bytes(emitted),
            )
            self.assertEqual(stat.S_IMODE(output.stat().st_mode), 0o444)
            with self.assertRaisesRegex(orchestrator.FinalEvidenceError, "cannot emit"):
                orchestrator.main(argv)
            os.chmod(output, 0o600)


if __name__ == "__main__":
    unittest.main()
