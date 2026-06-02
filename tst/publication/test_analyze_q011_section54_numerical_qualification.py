#!/usr/bin/env python3
"""Synthetic tests for source-local Q-011 Section 5.4 numerical aggregation."""

from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "frontier_control_plane"))

from . import analyze_q011_section54_campaign as admission
from . import analyze_q011_section54_numerical_qualification as qualifier
from . import q011_section54_particles as particles


def _sha256(character: str) -> str:
    return character * 64


def _identity(index: int, variant: str, seed: int) -> dict[str, object]:
    return {
        "variant": variant,
        "seed": seed,
        "physical_mode": qualifier.PHYSICAL_MODE,
        "attempt_id": f"baseline-{index:03d}-{variant}-seed-{seed}",
    }


def _pairing(identity: dict[str, object]) -> dict[str, object]:
    counterpart = None
    if identity["variant"] == qualifier.AMR_VARIANT:
        counterpart = {"variant": qualifier.FINE_VARIANT, "seed": identity["seed"]}
    elif identity["variant"] == qualifier.FINE_VARIANT:
        counterpart = {"variant": qualifier.AMR_VARIANT, "seed": identity["seed"]}
    return {
        "paired_seed_rule": qualifier.PAIRED_SEED_RULE,
        "pair_key": {"seed": identity["seed"]},
        "amr_fine_uniform_counterpart": counterpart,
        "comparison_status": "schema_wired_not_evaluated_by_artifact_admission_slice",
    }


def _tail_spectrum() -> dict[str, object]:
    edges = np.asarray(particles.CHI_BIN_EDGES)
    centers = np.sqrt(edges[:-1] * edges[1:])
    selected = (centers >= 20.0) & (centers <= 160.0)
    chi = centers[selected]
    macro_weight = 7.0 * chi**-1.5 * np.diff(edges)[selected]
    return particles.weighted_spectrum_record(chi, macro_weight)


def _particle_reductions() -> dict[str, object]:
    early = particles.weighted_spectrum_record([2.0], [1.0])
    late = _tail_spectrum()
    return {
        "t500": {
            "schema_version": particles.SCHEMA_VERSION,
            "record_type": "q011_section54_particle_snapshot_reduction",
            "snapshot_time_omega0_inverse": 500.0,
            "weighted_spectrum": early,
        },
        "t1200": {
            "schema_version": particles.SCHEMA_VERSION,
            "record_type": "q011_section54_particle_snapshot_reduction",
            "snapshot_time_omega0_inverse": 1200.0,
            "weighted_spectrum": late,
            "late_slope": particles.late_slope_record(late["f_chi"]),
        },
    }


def _spatial_reduction(*, passes: bool = True) -> dict[str, object]:
    mean = 2.0 if passes else 1.0
    return {
        "schema_version": 1,
        "record_type": "q011_section54_t500_spatial_reduction",
        "time_omega0_inverse": 500.0,
        "upstream_b_amplification": {
            "time_omega0_inverse": 500.0,
            "x_ideal_c_over_omega_pi": 0.0,
            "upstream_window_c_over_omega_pi": [120.0, 1200.0],
            "selected_cell_count": 1,
            "selected_area": 1.0,
            "mean_magnetic_magnitude": mean,
            "reference_b0": 1.0,
            "amplification_over_b0": mean,
            "acceptance_range": [1.2, 3.5],
            "passes_gate": passes,
        },
    }


def _admission(identity: dict[str, object], inventory_sha256: str) -> dict[str, object]:
    return {
        "schema_version": 1,
        "record_type": admission.RESULT_RECORD_TYPE,
        "campaign_id": admission.CAMPAIGN_ID,
        "qualification_scope": admission.QUALIFICATION_SCOPE,
        "admitted_for_follow_on_numerical_qualification": True,
        "final_claim_closure": False,
        "status": "admitted_for_follow_on_numerical_qualification",
        "failure_reasons": [],
        "admission": {
            "run_identity": copy.deepcopy(identity),
            "preregistration_binding": {
                "sha256": admission.EXPECTED_PREREGISTRATION_SHA256,
                "expected_sha256": admission.EXPECTED_PREREGISTRATION_SHA256,
            },
            "immutable_tree": {"inventory_sha256": inventory_sha256},
            "amr_pairing": _pairing(identity),
            "numerical_qualification_status": "not_evaluated_by_artifact_admission_slice",
        },
    }


def _attempts() -> list[dict[str, object]]:
    attempts = []
    for index, (variant, seed) in enumerate(qualifier.EXPECTED_MATRIX_CELLS):
        inventory_sha256 = f"{index + 1:064x}"
        identity = _identity(index, variant, seed)
        attempts.append(
            qualifier.bind_canonical_attempt(
                {
                    "schema_version": 1,
                    "record_type": qualifier.ATTEMPT_RECORD_TYPE,
                    "run_identity": identity,
                    "raw_inventory_sha256": inventory_sha256,
                    "admission_result": _admission(identity, inventory_sha256),
                    "particle_reductions": _particle_reductions(),
                    "spatial_reduction": _spatial_reduction(),
                }
            )
        )
    return attempts


def _residuals() -> list[dict[str, object]]:
    return [
        {
            "observable": observable,
            "maximum_absolute_difference": 0.0,
            "relative_mean_absolute": None if mean_bound is None else 0.0,
            "relative_root_mean_square": None if rms_bound is None else 0.0,
        }
        for observable, _, mean_bound, rms_bound in qualifier.PAIRED_RESIDUAL_THRESHOLDS
    ]


def _pairs(attempts: list[dict[str, object]]) -> list[dict[str, object]]:
    by_cell = {
        (
            wrapper["attempt"]["run_identity"]["variant"],
            wrapper["attempt"]["run_identity"]["seed"],
        ): wrapper["attempt_sha256"]
        for wrapper in attempts
    }
    return [
        {
            "seed": seed,
            "amr_attempt_sha256": by_cell[(qualifier.AMR_VARIANT, seed)],
            "fine_uniform_attempt_sha256": by_cell[(qualifier.FINE_VARIANT, seed)],
            "residuals": _residuals(),
        }
        for seed in qualifier.QUALIFYING_SEEDS
    ]


def _restart_binding(attempts: list[dict[str, object]]) -> dict[str, object]:
    source = next(
        wrapper["attempt_sha256"]
        for wrapper in attempts
        if wrapper["attempt"]["run_identity"]["variant"] == qualifier.AMR_VARIANT
    )
    return {
        "source_attempt_sha256": source,
        "uninterrupted": {"synthetic": "uninterrupted"},
        "continued": {"synthetic": "continued"},
    }


def _pending() -> dict[str, object]:
    return {
        "status": "pending_external_review",
        "reviewer": None,
        "reviewed_at_utc": None,
        "notes": "Synthetic aggregate remains pending named external review.",
    }


def _parity_result() -> dict[str, object]:
    return {
        "result": "pass_deterministic_continuation_parity",
        "checkpoint_time_omega0_inverse": 500.0,
        "retained_output_schedule_after_checkpoint_omega0_inverse": [
            600.0,
            700.0,
            800.0,
            900.0,
            1000.0,
            1100.0,
            1200.0,
        ],
        "maximum_absolute_difference_by_field": {"rho_bin": 0.0},
    }


def _qualify(
    attempts: list[dict[str, object]],
    *,
    pairs: list[dict[str, object]] | None = None,
    inventories: list[str] | None = None,
) -> dict[str, object]:
    with mock.patch.object(
        qualifier.restart,
        "compare_deterministic_continuation_parity",
        return_value=_parity_result(),
    ):
        return qualifier.qualify_numerical_aggregate(
            attempts=attempts,
            ordered_raw_inventory_sha256_values=inventories
            or [wrapper["attempt"]["raw_inventory_sha256"] for wrapper in attempts],
            paired_amr_fine_results=pairs or _pairs(attempts),
            restart_parity_binding=_restart_binding(attempts),
            reviewer_disposition=_pending(),
        )


class Q011Section54NumericalQualificationTests(unittest.TestCase):
    def test_complete_matrix_passes_numerically_but_remains_pending_review(self) -> None:
        attempts = _attempts()
        result = _qualify(attempts)
        self.assertEqual(result["attempt_count"], 24)
        self.assertTrue(result["numerical_gates_passed"])
        self.assertEqual(result["numerical_gate_status"], "passed")
        self.assertEqual(result["status"], "pending_external_review")
        self.assertFalse(result["final_claim_closure"])
        self.assertEqual(
            result["ordered_raw_inventory_sha256_values"],
            [wrapper["attempt"]["raw_inventory_sha256"] for wrapper in attempts],
        )

    def test_missing_duplicate_and_matrix_cell_drift_fail_closed(self) -> None:
        attempts = _attempts()
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "exactly 24"):
            _qualify(attempts[:-1], pairs=_pairs(attempts))

        duplicate = copy.deepcopy(attempts)
        duplicate[1] = copy.deepcopy(duplicate[0])
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "duplicate canonical"):
            _qualify(duplicate)

        drift = copy.deepcopy(attempts)
        drift[1]["attempt"]["run_identity"]["seed"] = qualifier.QUALIFYING_SEEDS[0]
        drift[1]["attempt"]["admission_result"]["admission"]["run_identity"]["seed"] = (
            qualifier.QUALIFYING_SEEDS[0]
        )
        drift[1]["attempt"]["admission_result"]["admission"]["amr_pairing"]["pair_key"]["seed"] = (
            qualifier.QUALIFYING_SEEDS[0]
        )
        drift[1] = qualifier.bind_canonical_attempt(drift[1]["attempt"])
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "matrix order drifted"):
            _qualify(drift)

    def test_cross_seed_amr_fine_pair_fails_closed(self) -> None:
        attempts = _attempts()
        pairs = _pairs(attempts)
        pairs[0]["fine_uniform_attempt_sha256"] = pairs[1]["fine_uniform_attempt_sha256"]
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "matched seed"):
            _qualify(attempts, pairs=pairs)

    def test_attempt_and_pair_numerical_failures_aggregate_without_claim_closure(self) -> None:
        attempts = _attempts()
        attempts[0]["attempt"]["spatial_reduction"] = _spatial_reduction(passes=False)
        attempts[0] = qualifier.bind_canonical_attempt(attempts[0]["attempt"])
        result = _qualify(attempts)
        self.assertFalse(result["numerical_gates_passed"])
        self.assertEqual(result["numerical_gate_status"], "failed")
        self.assertEqual(result["status"], "pending_external_review")

        attempts = _attempts()
        pairs = _pairs(attempts)
        pairs[0]["residuals"][0]["maximum_absolute_difference"] = 241.0
        result = _qualify(attempts, pairs=pairs)
        self.assertFalse(result["numerical_gates_passed"])
        self.assertFalse(result["paired_amr_fine_results"][0]["gates_passed"])

    def test_canonical_attempt_tamper_and_inventory_order_tamper_fail_closed(self) -> None:
        attempts = _attempts()
        attempts[0]["attempt"]["raw_inventory_sha256"] = _sha256("f")
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "canonical attempt SHA"):
            _qualify(attempts)

        attempts = _attempts()
        inventories = [
            wrapper["attempt"]["raw_inventory_sha256"] for wrapper in attempts
        ]
        inventories[0], inventories[1] = inventories[1], inventories[0]
        with self.assertRaisesRegex(qualifier.NumericalQualificationError, "inventory digest binding"):
            _qualify(attempts, inventories=inventories)


if __name__ == "__main__":
    unittest.main()
