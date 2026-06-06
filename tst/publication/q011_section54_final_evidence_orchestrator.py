#!/usr/bin/env python3
"""Build a non-authorizing Q-011 Section 5.4 claim-candidate manifest.

This component cross-binds already admitted campaign identities, the source-local
numerical aggregate, and a separately reviewed independent recomputation.  It
does not inspect raw simulation artifacts, authorize execution, mutate policy,
promote a claim, or replace terminal external review.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
import sys
from typing import Any, Mapping, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent / "frontier_control_plane"))

if __package__:
    from . import analyze_q011_section54_campaign as admission
    from . import analyze_q011_section54_numerical_qualification as numerical
    from . import q011_section54_artifacts as artifacts
else:
    import analyze_q011_section54_campaign as admission
    import analyze_q011_section54_numerical_qualification as numerical
    import q011_section54_artifacts as artifacts


SCHEMA_VERSION = 1
RECORD_TYPE = "q011_section54_claim_candidate_evidence_manifest"
INDEPENDENT_RECOMPUTE_RECORD_TYPE = (
    "q011_section54_independent_raw_artifact_recompute_result"
)
CLAIM_ID = "CLAIM-PAPER-SHOCK-001"
EVIDENCE_CLASS = "sun_bai_2023_reproduction"
QUALIFICATION_SCOPE = (
    "q011_section54_claim_candidate_only_external_review_and_separate_claim_"
    "promotion_required"
)
CLAIM_REQUIRED_GATES = (
    "Q-003",
    "Q-004",
    "Q-009",
    "Q-011",
    "Q-016",
    "Q-023",
    "Q-025",
    "Q-026",
    "Q-027",
)
ATTEMPT_PRIMARY_OBSERVABLES = (
    "shock_front_position_at_t500",
    "upstream_magnetic_amplification_at_t500",
    "downstream_chi_f_chi_at_t500",
    "downstream_chi_f_chi_at_t1200",
    "late_energy_tail_slope_at_t1200",
)
PAIR_OBSERVABLE_PREFIX = "paired_amr_fine:"
RESTART_OBSERVABLE = "restart_continuation_parity_maximum_absolute_difference"
RESTART_OUTPUT_SLOTS = (600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0)
RESTART_FIELDS = (
    "rho_bin",
    "bmag_bin",
    "prtcl_jx_bin",
    "j2_bin",
    "prtcl_all_pvtk_integer_payload",
    "prtcl_all_pvtk_float_payload",
)
REMAINING_GATES = (
    "named_terminal_external_review",
    "separate_reviewed_claim_registry_promotion",
    "manuscript_claim_wording_and_figure_review",
)
MAX_JSON_BYTES = 64 * 1024 * 1024
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_ATTEMPT_ID_PATTERN = re.compile(r"[a-z0-9][a-z0-9._-]{0,127}")
_READ_FLAGS = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)


class FinalEvidenceError(ValueError):
    """Raised when Q-011 final-evidence orchestration cannot proceed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FinalEvidenceError(message)


def _object(value: object, keys: set[str], *, label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label}: expected object")
    _require(set(value) == keys, f"{label}: keys drifted")
    return value


def _array(value: object, *, label: str) -> list[Any]:
    _require(type(value) is list, f"{label}: expected array")
    return value


def _text(value: object, *, label: str) -> str:
    _require(type(value) is str and bool(value.strip()), f"{label}: expected nonempty text")
    return value


def _sha256(value: object, *, label: str) -> str:
    _require(
        type(value) is str and _SHA256_PATTERN.fullmatch(value) is not None,
        f"{label}: malformed SHA-256",
    )
    return value


def _finite_float(value: object, *, label: str, minimum: float | None = None) -> float:
    _require(type(value) is float and math.isfinite(value), f"{label}: expected finite float")
    if minimum is not None:
        _require(value >= minimum, f"{label}: expected value >= {minimum}")
    return value


def _positive_integer(value: object, *, label: str) -> int:
    _require(type(value) is int and value > 0, f"{label}: expected positive integer")
    return value


def _strict_equal(left: object, right: object) -> bool:
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return left.keys() == right.keys() and all(
            _strict_equal(left[key], right[key]) for key in left
        )
    if type(left) is list:
        return len(left) == len(right) and all(
            _strict_equal(left_item, right_item)
            for left_item, right_item in zip(left, right)
        )
    return left == right


def _require_strict_equal(left: object, right: object, *, label: str) -> None:
    _require(_strict_equal(left, right), f"{label}: value drifted")


def _canonical_json_bytes(value: object) -> bytes:
    try:
        return artifacts.canonical_json_bytes(value)
    except artifacts.DerivedArtifactError as error:
        raise FinalEvidenceError("record is not canonical finite JSON") from error


def _canonical_sha256(value: object) -> str:
    return hashlib.sha256(_canonical_json_bytes(value)).hexdigest()


def _reviewer_disposition(
    value: object, *, required_status: str, label: str
) -> dict[str, object]:
    try:
        disposition = artifacts.validate_reviewer_disposition(value)
    except artifacts.DerivedArtifactError as error:
        raise FinalEvidenceError(f"{label}: invalid reviewer disposition") from error
    _require(
        disposition["status"] == required_status,
        f"{label}: must remain {required_status}",
    )
    return disposition


def _validate_campaign_matrix(value: object) -> dict[str, object]:
    matrix = _object(
        value,
        {
            "physical_mode",
            "grid_variants",
            "qualifying_seeds",
            "expected_baseline_attempts",
            "paired_seed_rule",
        },
        label="numerical_aggregate/campaign_matrix",
    )
    expected = {
        "physical_mode": numerical.PHYSICAL_MODE,
        "grid_variants": list(numerical.GRID_VARIANTS),
        "qualifying_seeds": list(numerical.QUALIFYING_SEEDS),
        "expected_baseline_attempts": numerical.EXPECTED_BASELINE_ATTEMPTS,
        "paired_seed_rule": numerical.PAIRED_SEED_RULE,
    }
    _require_strict_equal(matrix, expected, label="numerical_aggregate/campaign_matrix")
    return expected


def _validate_identity(value: object, *, label: str) -> dict[str, object]:
    identity = _object(
        value,
        {"variant", "seed", "physical_mode", "attempt_id"},
        label=label,
    )
    _require(identity["variant"] in numerical.GRID_VARIANTS, f"{label}: variant drifted")
    _require(
        type(identity["seed"]) is int and identity["seed"] in numerical.QUALIFYING_SEEDS,
        f"{label}: qualifying seed drifted",
    )
    _require(identity["physical_mode"] == numerical.PHYSICAL_MODE, f"{label}: mode drifted")
    _require(
        type(identity["attempt_id"]) is str
        and _ATTEMPT_ID_PATTERN.fullmatch(identity["attempt_id"]) is not None,
        f"{label}: attempt ID is noncanonical",
    )
    return dict(identity)


def _validate_checkpoint_lineage(
    value: object, *, identity: Mapping[str, object], label: str
) -> dict[str, object]:
    lineage = _object(
        value,
        {
            "nominal_slot_time",
            "observed_committed_time",
            "retained_attempt_id",
            "restart_manifest_path",
            "restart_member_path",
            "retained_restart_member_absolute_path",
            "restart_member_sha256",
        },
        label=label,
    )
    nominal = _finite_float(
        lineage["nominal_slot_time"], label=f"{label}/nominal_slot_time", minimum=0.0
    )
    observed = _finite_float(
        lineage["observed_committed_time"],
        label=f"{label}/observed_committed_time",
        minimum=0.0,
    )
    _require(
        nominal == 500.0,
        f"{label}: checkpoint nominal slot time drifted",
    )
    _require(
        lineage["retained_attempt_id"] == identity["attempt_id"],
        f"{label}: retained attempt identity drifted",
    )
    _text(lineage["restart_manifest_path"], label=f"{label}/restart_manifest_path")
    _text(lineage["restart_member_path"], label=f"{label}/restart_member_path")
    absolute = _text(
        lineage["retained_restart_member_absolute_path"],
        label=f"{label}/retained_restart_member_absolute_path",
    )
    _require(absolute.startswith("/"), f"{label}: retained restart path must be absolute")
    _sha256(lineage["restart_member_sha256"], label=f"{label}/restart_member_sha256")
    parsed = dict(lineage)
    parsed["nominal_slot_time"] = nominal
    parsed["observed_committed_time"] = observed
    return parsed


def _validate_snapshot_time(
    value: object, *, expected_nominal_slot_time: float, label: str
) -> dict[str, float]:
    snapshot = _object(
        value,
        {"nominal_slot_time", "observed_committed_time"},
        label=label,
    )
    nominal = _finite_float(
        snapshot["nominal_slot_time"], label=f"{label}/nominal_slot_time", minimum=0.0
    )
    observed = _finite_float(
        snapshot["observed_committed_time"],
        label=f"{label}/observed_committed_time",
        minimum=0.0,
    )
    _require(
        nominal == expected_nominal_slot_time,
        f"{label}: nominal slot time drifted",
    )
    return {
        "nominal_slot_time": nominal,
        "observed_committed_time": observed,
    }


def _validate_attempts(
    aggregate: Mapping[str, object],
    admitted_bundle_identities: Sequence[object],
) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    attempts = _array(aggregate["attempt_gate_results"], label="numerical_aggregate/attempts")
    _require(
        len(attempts) == numerical.EXPECTED_BASELINE_ATTEMPTS,
        "numerical_aggregate/attempts: expected exactly 24 attempts",
    )
    parsed_attempts = []
    by_sha256: dict[str, dict[str, object]] = {}
    observed_cells = []
    for index, (value, expected_cell) in enumerate(zip(attempts, numerical.EXPECTED_MATRIX_CELLS)):
        label = f"numerical_aggregate/attempts[{index}]"
        attempt = _object(
            value,
            {
                "attempt_sha256",
                "identity",
                "raw_inventory_sha256",
                "source_checkpoint_lineage",
                "snapshot_times",
                "gates_passed",
            },
            label=label,
        )
        attempt_sha256 = _sha256(attempt["attempt_sha256"], label=f"{label}/attempt_sha256")
        _require(attempt_sha256 not in by_sha256, f"{label}: duplicate attempt SHA-256")
        identity = _validate_identity(attempt["identity"], label=f"{label}/identity")
        cell = (identity["variant"], identity["seed"])
        _require(cell == expected_cell, f"{label}: canonical campaign cell order drifted")
        observed_cells.append(cell)
        raw_inventory = _sha256(
            attempt["raw_inventory_sha256"], label=f"{label}/raw_inventory_sha256"
        )
        checkpoint_lineage = _validate_checkpoint_lineage(
            attempt["source_checkpoint_lineage"],
            identity=identity,
            label=f"{label}/source_checkpoint_lineage",
        )
        snapshot_times = _object(
            attempt["snapshot_times"], {"t500", "t1200"}, label=f"{label}/snapshot_times"
        )
        t500 = _validate_snapshot_time(
            snapshot_times["t500"],
            expected_nominal_slot_time=500.0,
            label=f"{label}/snapshot_times/t500",
        )
        t1200 = _validate_snapshot_time(
            snapshot_times["t1200"],
            expected_nominal_slot_time=1200.0,
            label=f"{label}/snapshot_times/t1200",
        )
        _require(
            checkpoint_lineage["nominal_slot_time"] == t500["nominal_slot_time"]
            and checkpoint_lineage["observed_committed_time"]
            == t500["observed_committed_time"],
            f"{label}: source checkpoint and t500 snapshot time drifted",
        )
        _require(
            t1200["observed_committed_time"] > t500["observed_committed_time"],
            f"{label}: observed snapshot time sequence drifted",
        )
        _require(attempt["gates_passed"] is True, f"{label}: attempt numerical gates failed")
        parsed = {
            "attempt_sha256": attempt_sha256,
            "identity": identity,
            "raw_inventory_sha256": raw_inventory,
            "source_checkpoint_lineage": checkpoint_lineage,
            "snapshot_times": {"t500": t500, "t1200": t1200},
        }
        parsed_attempts.append(parsed)
        by_sha256[attempt_sha256] = parsed
    _require(
        observed_cells == list(numerical.EXPECTED_MATRIX_CELLS),
        "numerical_aggregate/attempts: campaign matrix is incomplete",
    )
    _require(
        type(aggregate["attempt_count"]) is int
        and aggregate["attempt_count"] == numerical.EXPECTED_BASELINE_ATTEMPTS,
        "numerical_aggregate/attempt_count drifted",
    )
    _require_strict_equal(
        aggregate["attempt_sha256_values"],
        [item["attempt_sha256"] for item in parsed_attempts],
        label="numerical_aggregate/attempt_sha256_values",
    )
    _require_strict_equal(
        aggregate["ordered_raw_inventory_sha256_values"],
        [item["raw_inventory_sha256"] for item in parsed_attempts],
        label="numerical_aggregate/ordered_raw_inventory_sha256_values",
    )

    bundles = _array(admitted_bundle_identities, label="admitted_bundle_identities")
    _require(
        len(bundles) == numerical.EXPECTED_BASELINE_ATTEMPTS,
        "admitted_bundle_identities: expected exactly 24 identities",
    )
    normalized_bundles = []
    admission_sha256_values = set()
    for index, (value, attempt) in enumerate(zip(bundles, parsed_attempts)):
        label = f"admitted_bundle_identities[{index}]"
        bundle = _object(
            value,
            {
                "attempt_id",
                "variant",
                "seed",
                "attempt_sha256",
                "raw_inventory_sha256",
                "admission_result_sha256",
                "admitted_for_follow_on_numerical_qualification",
            },
            label=label,
        )
        normalized = {
            "attempt_id": _text(bundle["attempt_id"], label=f"{label}/attempt_id"),
            "variant": _text(bundle["variant"], label=f"{label}/variant"),
            "seed": bundle["seed"],
            "attempt_sha256": _sha256(bundle["attempt_sha256"], label=f"{label}/attempt_sha256"),
            "raw_inventory_sha256": _sha256(
                bundle["raw_inventory_sha256"], label=f"{label}/raw_inventory_sha256"
            ),
            "admission_result_sha256": _sha256(
                bundle["admission_result_sha256"], label=f"{label}/admission_result_sha256"
            ),
            "admitted_for_follow_on_numerical_qualification": bundle[
                "admitted_for_follow_on_numerical_qualification"
            ],
        }
        _require(type(normalized["seed"]) is int, f"{label}: seed must be an integer")
        _require(
            normalized["admitted_for_follow_on_numerical_qualification"] is True,
            f"{label}: bundle was not admitted for numerical qualification",
        )
        expected = {
            "attempt_id": attempt["identity"]["attempt_id"],
            "variant": attempt["identity"]["variant"],
            "seed": attempt["identity"]["seed"],
            "attempt_sha256": attempt["attempt_sha256"],
            "raw_inventory_sha256": attempt["raw_inventory_sha256"],
        }
        for key, expected_value in expected.items():
            _require(normalized[key] == expected_value, f"{label}: {key} binding drifted")
        _require(
            normalized["admission_result_sha256"] not in admission_sha256_values,
            f"{label}: duplicate admission-result SHA-256",
        )
        admission_sha256_values.add(normalized["admission_result_sha256"])
        normalized_bundles.append(normalized)
    return normalized_bundles, by_sha256


def _validate_pairs(
    value: object, *, by_sha256: Mapping[str, Mapping[str, object]]
) -> tuple[list[dict[str, object]], dict[tuple[int, str], float]]:
    pairs = _array(value, label="numerical_aggregate/paired_amr_fine_results")
    _require(
        len(pairs) == len(numerical.QUALIFYING_SEEDS),
        "numerical_aggregate/paired_amr_fine_results: coverage count drifted",
    )
    by_cell = {
        (item["identity"]["variant"], item["identity"]["seed"]): sha256
        for sha256, item in by_sha256.items()
    }
    coverage = []
    production_metrics: dict[tuple[int, str], float] = {}
    for index, (value, seed) in enumerate(zip(pairs, numerical.QUALIFYING_SEEDS)):
        label = f"numerical_aggregate/paired_amr_fine_results[{index}]"
        pair = _object(
            value,
            {
                "seed",
                "amr_attempt_sha256",
                "fine_uniform_attempt_sha256",
                "residuals",
                "gates_passed",
            },
            label=label,
        )
        _require(type(pair["seed"]) is int and pair["seed"] == seed, f"{label}: seed drifted")
        amr_sha256 = _sha256(pair["amr_attempt_sha256"], label=f"{label}/amr_attempt_sha256")
        fine_sha256 = _sha256(
            pair["fine_uniform_attempt_sha256"], label=f"{label}/fine_uniform_attempt_sha256"
        )
        _require(
            amr_sha256 == by_cell[(numerical.AMR_VARIANT, seed)],
            f"{label}: AMR pair identity drifted",
        )
        _require(
            fine_sha256 == by_cell[(numerical.FINE_VARIANT, seed)],
            f"{label}: fine-uniform pair identity drifted",
        )
        _require(pair["gates_passed"] is True, f"{label}: paired numerical gates failed")
        residuals = _array(pair["residuals"], label=f"{label}/residuals")
        _require(
            len(residuals) == len(numerical.PAIRED_RESIDUAL_THRESHOLDS),
            f"{label}: residual coverage drifted",
        )
        for residual_index, (residual_value, threshold) in enumerate(
            zip(residuals, numerical.PAIRED_RESIDUAL_THRESHOLDS)
        ):
            residual_label = f"{label}/residuals[{residual_index}]"
            residual = _object(
                residual_value,
                {
                    "observable",
                    "maximum_absolute_difference",
                    "relative_mean_absolute",
                    "relative_root_mean_square",
                },
                label=residual_label,
            )
            _require(
                residual["observable"] == threshold[0],
                f"{residual_label}: observable order drifted",
            )
            maximum = _finite_float(
                residual["maximum_absolute_difference"],
                label=f"{residual_label}/maximum_absolute_difference",
                minimum=0.0,
            )
            for field, bound in (
                ("relative_mean_absolute", threshold[2]),
                ("relative_root_mean_square", threshold[3]),
            ):
                if bound is None:
                    _require(residual[field] is None, f"{residual_label}/{field}: expected null")
                else:
                    _finite_float(residual[field], label=f"{residual_label}/{field}", minimum=0.0)
            production_metrics[(seed, threshold[0])] = maximum
        coverage.append(
            {
                "seed": seed,
                "amr_attempt_sha256": amr_sha256,
                "fine_uniform_attempt_sha256": fine_sha256,
                "pair_result_sha256": _canonical_sha256(pair),
            }
        )
    return coverage, production_metrics


def _validate_restart(
    value: object, *, by_sha256: Mapping[str, Mapping[str, object]]
) -> tuple[dict[str, object], float]:
    restart = _object(
        value,
        {
            "source_attempt_sha256",
            "restart_parity_binding_sha256",
            "screen_scope",
            "full_state_equivalence_claimed",
            "result",
        },
        label="numerical_aggregate/restart_parity",
    )
    source_sha256 = _sha256(
        restart["source_attempt_sha256"],
        label="numerical_aggregate/restart_parity/source_attempt_sha256",
    )
    _require(source_sha256 in by_sha256, "restart parity source attempt is unknown")
    source_identity = by_sha256[source_sha256]["identity"]
    _require(
        (source_identity["variant"], source_identity["seed"])
        == (numerical.AMR_VARIANT, numerical.QUALIFYING_SEEDS[0]),
        "restart parity source must be the preregistered first-seed AMR baseline",
    )
    binding_sha256 = _sha256(
        restart["restart_parity_binding_sha256"],
        label="numerical_aggregate/restart_parity/restart_parity_binding_sha256",
    )
    _require(
        restart["screen_scope"] == numerical.RESTART_SCREEN_SCOPE,
        "restart parity must use the retained production screen",
    )
    _require(
        restart["full_state_equivalence_claimed"] is False,
        "restart parity must not claim full-state equivalence",
    )
    result = _object(
        restart["result"],
        {
            "result",
            "checkpoint_nominal_slot_omega0_inverse",
            "checkpoint_observed_committed_cycle",
            "checkpoint_observed_committed_time_omega0_inverse",
            "retained_output_nominal_slots_after_checkpoint_omega0_inverse",
            "paired_output_observed_commits",
            "maximum_absolute_difference_by_field",
        },
        label="numerical_aggregate/restart_parity/result",
    )
    _require(
        result["result"] == "pass_deterministic_continuation_parity",
        "restart continuation parity did not pass",
    )
    _require(
        type(result["checkpoint_nominal_slot_omega0_inverse"]) is float
        and result["checkpoint_nominal_slot_omega0_inverse"] == 500.0,
        "restart checkpoint nominal slot drifted",
    )
    _positive_integer(
        result["checkpoint_observed_committed_cycle"],
        label="restart/checkpoint_observed_committed_cycle",
    )
    _finite_float(
        result["checkpoint_observed_committed_time_omega0_inverse"],
        label="restart/checkpoint_observed_committed_time_omega0_inverse",
        minimum=0.0,
    )
    _require_strict_equal(
        result["retained_output_nominal_slots_after_checkpoint_omega0_inverse"],
        list(RESTART_OUTPUT_SLOTS),
        label="restart/retained_output_nominal_slots_after_checkpoint_omega0_inverse",
    )
    observed = _array(
        result["paired_output_observed_commits"],
        label="restart/paired_output_observed_commits",
    )
    _require(len(observed) == len(RESTART_OUTPUT_SLOTS), "restart observed output count drifted")
    previous_cycle = result["checkpoint_observed_committed_cycle"]
    previous_time = result["checkpoint_observed_committed_time_omega0_inverse"]
    for index, (commit, slot) in enumerate(zip(observed, RESTART_OUTPUT_SLOTS)):
        label = f"restart/paired_output_observed_commits[{index}]"
        parsed = _object(
            commit,
            {
                "nominal_slot_omega0_inverse",
                "observed_committed_cycle",
                "observed_committed_time_omega0_inverse",
            },
            label=label,
        )
        _require(
            type(parsed["nominal_slot_omega0_inverse"]) is float
            and parsed["nominal_slot_omega0_inverse"] == slot,
            f"{label}: nominal slot drifted",
        )
        cycle = _positive_integer(parsed["observed_committed_cycle"], label=f"{label}/cycle")
        time = _finite_float(
            parsed["observed_committed_time_omega0_inverse"],
            label=f"{label}/time",
            minimum=0.0,
        )
        _require(cycle > previous_cycle and time > previous_time, f"{label}: sequence drifted")
        previous_cycle = cycle
        previous_time = time
    maxima = _object(
        result["maximum_absolute_difference_by_field"],
        set(RESTART_FIELDS),
        label="restart/maximum_absolute_difference_by_field",
    )
    parsed_maxima = [
        _finite_float(maxima[field], label=f"restart/maximum/{field}", minimum=0.0)
        for field in RESTART_FIELDS
    ]
    return (
        {
            "source_attempt_sha256": source_sha256,
            "restart_parity_binding_sha256": binding_sha256,
            "restart_result_sha256": _canonical_sha256(result),
            "screen_scope": numerical.RESTART_SCREEN_SCOPE,
            "full_state_equivalence_claimed": False,
            "result": "pass_deterministic_continuation_parity",
        },
        max(parsed_maxima),
    )


def _validate_numerical_aggregate(
    value: object, admitted_bundle_identities: Sequence[object]
) -> dict[str, object]:
    aggregate = _object(
        value,
        {
            "schema_version",
            "record_type",
            "campaign_id",
            "qualification_scope",
            "final_claim_closure",
            "status",
            "numerical_gate_status",
            "numerical_gates_passed",
            "campaign_matrix",
            "attempt_count",
            "attempt_sha256_values",
            "ordered_raw_inventory_sha256_values",
            "attempt_gate_results",
            "paired_amr_fine_results",
            "restart_parity",
            "reviewer_disposition",
        },
        label="numerical_aggregate",
    )
    _require(
        type(aggregate["schema_version"]) is int
        and aggregate["schema_version"] == numerical.SCHEMA_VERSION,
        "numerical_aggregate: schema version drifted",
    )
    _require(aggregate["record_type"] == numerical.RECORD_TYPE, "numerical_aggregate: type drifted")
    _require(
        aggregate["campaign_id"] == admission.CAMPAIGN_ID,
        "numerical_aggregate: campaign drifted",
    )
    _require(
        aggregate["qualification_scope"] == numerical.QUALIFICATION_SCOPE,
        "numerical_aggregate: qualification scope drifted",
    )
    _require(aggregate["final_claim_closure"] is False, "numerical_aggregate closed a claim")
    _require(
        aggregate["status"] == "pending_external_review",
        "numerical_aggregate: status must remain pending external review",
    )
    _require(
        aggregate["numerical_gate_status"] == "passed"
        and aggregate["numerical_gates_passed"] is True,
        "numerical_aggregate: numerical gates did not pass",
    )
    matrix = _validate_campaign_matrix(aggregate["campaign_matrix"])
    bundles, by_sha256 = _validate_attempts(aggregate, admitted_bundle_identities)
    pairs, pair_metrics = _validate_pairs(
        aggregate["paired_amr_fine_results"], by_sha256=by_sha256
    )
    restart, restart_metric = _validate_restart(aggregate["restart_parity"], by_sha256=by_sha256)
    aggregate_disposition = _reviewer_disposition(
        aggregate["reviewer_disposition"],
        required_status="pending_external_review",
        label="numerical_aggregate/reviewer_disposition",
    )
    return {
        "sha256": _canonical_sha256(aggregate),
        "matrix": matrix,
        "bundles": bundles,
        "by_sha256": by_sha256,
        "pairs": pairs,
        "pair_metrics": pair_metrics,
        "restart": restart,
        "restart_metric": restart_metric,
        "aggregate_reviewer_disposition": aggregate_disposition,
        "ordered_raw_inventory_sha256_values": list(
            aggregate["ordered_raw_inventory_sha256_values"]
        ),
    }


def _binding(value: object, *, label: str) -> dict[str, str]:
    binding = _object(value, {"path_or_archive_locator", "sha256"}, label=label)
    return {
        "path_or_archive_locator": _text(
            binding["path_or_archive_locator"], label=f"{label}/path_or_archive_locator"
        ),
        "sha256": _sha256(binding["sha256"], label=f"{label}/sha256"),
    }


def _expected_metric_rows(
    aggregate: Mapping[str, object],
) -> list[tuple[tuple[str, str, int, str], float | None]]:
    rows: list[tuple[tuple[str, str, int, str], float | None]] = []
    by_cell = {
        (item["identity"]["variant"], item["identity"]["seed"]): item
        for item in aggregate["by_sha256"].values()
    }
    for variant, seed in numerical.EXPECTED_MATRIX_CELLS:
        identity = by_cell[(variant, seed)]["identity"]
        for observable in ATTEMPT_PRIMARY_OBSERVABLES:
            rows.append(((identity["attempt_id"], variant, seed, observable), None))
    for seed in numerical.QUALIFYING_SEEDS:
        identity = by_cell[(numerical.AMR_VARIANT, seed)]["identity"]
        for threshold in numerical.PAIRED_RESIDUAL_THRESHOLDS:
            observable = threshold[0]
            rows.append(
                (
                    (
                        identity["attempt_id"],
                        numerical.AMR_VARIANT,
                        seed,
                        f"{PAIR_OBSERVABLE_PREFIX}{observable}",
                    ),
                    aggregate["pair_metrics"][(seed, observable)],
                )
            )
    restart_identity = aggregate["by_sha256"][
        aggregate["restart"]["source_attempt_sha256"]
    ]["identity"]
    rows.append(
        (
            (
                restart_identity["attempt_id"],
                restart_identity["variant"],
                restart_identity["seed"],
                RESTART_OBSERVABLE,
            ),
            aggregate["restart_metric"],
        )
    )
    return rows


def _validate_metric_row(
    value: object,
    *,
    expected_key: tuple[str, str, int, str],
    expected_production_metric: float | None,
    label: str,
) -> dict[str, object]:
    row = _object(
        value,
        {
            "attempt_id",
            "variant",
            "qualifying_seed",
            "observable",
            "production_metric",
            "independent_metric",
            "absolute_difference",
            "relative_difference",
            "declared_tolerance",
            "disposition",
        },
        label=label,
    )
    _text(row["attempt_id"], label=f"{label}/attempt_id")
    _text(row["variant"], label=f"{label}/variant")
    _require(type(row["qualifying_seed"]) is int, f"{label}: qualifying seed must be an integer")
    _text(row["observable"], label=f"{label}/observable")
    key = (row["attempt_id"], row["variant"], row["qualifying_seed"], row["observable"])
    _require(key == expected_key, f"{label}: canonical metric coverage order drifted")
    production = _finite_float(row["production_metric"], label=f"{label}/production_metric")
    independent = _finite_float(row["independent_metric"], label=f"{label}/independent_metric")
    absolute = _finite_float(
        row["absolute_difference"], label=f"{label}/absolute_difference", minimum=0.0
    )
    tolerance = _finite_float(
        row["declared_tolerance"], label=f"{label}/declared_tolerance", minimum=0.0
    )
    _require(
        math.isclose(absolute, abs(production - independent), rel_tol=1.0e-12, abs_tol=1.0e-15),
        f"{label}: absolute difference is inconsistent",
    )
    relative = row["relative_difference"]
    if production == 0.0:
        _require(
            relative is None,
            f"{label}: zero production metric requires null relative difference",
        )
    else:
        parsed_relative = _finite_float(
            relative, label=f"{label}/relative_difference", minimum=0.0
        )
        _require(
            math.isclose(
                parsed_relative,
                absolute / abs(production),
                rel_tol=1.0e-12,
                abs_tol=1.0e-15,
            ),
            f"{label}: relative difference is inconsistent",
        )
    _require(
        row["disposition"] == "pass_within_declared_tolerance",
        f"{label}: independent metric disposition did not pass",
    )
    _require(absolute <= tolerance, f"{label}: declared tolerance exceeded")
    if expected_production_metric is not None:
        _require(
            math.isclose(
                production,
                expected_production_metric,
                rel_tol=1.0e-12,
                abs_tol=1.0e-15,
            ),
            f"{label}: production aggregate metric binding drifted",
        )
    return dict(row)


def _validate_independent_recompute(
    value: object, *, aggregate: Mapping[str, object]
) -> dict[str, object]:
    recompute = _object(
        value,
        {
            "schema_version",
            "record_type",
            "campaign_id",
            "status",
            "production_aggregate_sha256",
            "independent_script",
            "independent_environment_lock",
            "production_helper_imports_authorized",
            "production_helper_imports_detected",
            "attempt_inventory_sha256_values",
            "input_artifact_checksums",
            "metric_comparison_table",
            "reviewer_identity",
            "reviewer_disposition",
        },
        label="independent_recompute",
    )
    _require(
        type(recompute["schema_version"]) is int
        and recompute["schema_version"] == SCHEMA_VERSION,
        "independent_recompute: schema version drifted",
    )
    _require(
        recompute["record_type"] == INDEPENDENT_RECOMPUTE_RECORD_TYPE,
        "independent_recompute: record type drifted",
    )
    _require(
        recompute["campaign_id"] == admission.CAMPAIGN_ID,
        "independent_recompute: campaign drifted",
    )
    _require(
        recompute["status"] == "pass_independent_recompute_reviewed",
        "independent_recompute: reviewed pass is required",
    )
    _require(
        recompute["production_aggregate_sha256"] == aggregate["sha256"],
        "independent_recompute: production aggregate binding drifted",
    )
    script = _binding(recompute["independent_script"], label="independent_recompute/script")
    environment = _binding(
        recompute["independent_environment_lock"],
        label="independent_recompute/environment_lock",
    )
    _require(
        recompute["production_helper_imports_authorized"] is False
        and recompute["production_helper_imports_detected"] is False,
        "independent_recompute: production helper imports are forbidden",
    )
    _require_strict_equal(
        recompute["attempt_inventory_sha256_values"],
        aggregate["ordered_raw_inventory_sha256_values"],
        label="independent_recompute/attempt_inventory_sha256_values",
    )
    checksums = _array(
        recompute["input_artifact_checksums"],
        label="independent_recompute/input_artifact_checksums",
    )
    parsed_checksums = [
        _binding(value, label=f"independent_recompute/input_artifact_checksums[{index}]")
        for index, value in enumerate(checksums)
    ]
    _require(
        len({item["path_or_archive_locator"] for item in parsed_checksums})
        == len(parsed_checksums),
        "independent_recompute/input_artifact_checksums: duplicate locator",
    )
    checksum_values = {item["sha256"] for item in parsed_checksums}
    _require(
        set(aggregate["ordered_raw_inventory_sha256_values"]) <= checksum_values,
        "independent_recompute/input_artifact_checksums: attempt inventories are incomplete",
    )
    expected_rows = _expected_metric_rows(aggregate)
    rows = _array(
        recompute["metric_comparison_table"],
        label="independent_recompute/metric_comparison_table",
    )
    _require(
        len(rows) == len(expected_rows),
        "independent_recompute/metric_comparison_table: coverage count drifted",
    )
    parsed_rows = [
        _validate_metric_row(
            row,
            expected_key=expected_key,
            expected_production_metric=expected_metric,
            label=f"independent_recompute/metric_comparison_table[{index}]",
        )
        for index, (row, (expected_key, expected_metric)) in enumerate(
            zip(rows, expected_rows)
        )
    ]
    reviewer_identity = _text(
        recompute["reviewer_identity"], label="independent_recompute/reviewer_identity"
    )
    disposition = _reviewer_disposition(
        recompute["reviewer_disposition"],
        required_status="accepted",
        label="independent_recompute/reviewer_disposition",
    )
    _require(
        disposition["reviewer"] == reviewer_identity,
        "independent_recompute: reviewer identity binding drifted",
    )
    return {
        "sha256": _canonical_sha256(recompute),
        "status": "pass_independent_recompute_reviewed",
        "independent_script": script,
        "independent_environment_lock": environment,
        "metric_comparison_row_count": len(parsed_rows),
        "reviewer_identity": reviewer_identity,
        "reviewer_disposition": disposition,
    }


def build_claim_candidate_evidence_manifest(
    *,
    admitted_bundle_identities: Sequence[object],
    numerical_aggregate: object,
    independent_recompute: object,
    external_reviewer_disposition: object,
) -> dict[str, object]:
    """Cross-bind complete Q-011 evidence without authorizing claim closure."""
    aggregate = _validate_numerical_aggregate(
        numerical_aggregate, admitted_bundle_identities
    )
    recompute = _validate_independent_recompute(independent_recompute, aggregate=aggregate)
    external_review = _reviewer_disposition(
        external_reviewer_disposition,
        required_status="pending_external_review",
        label="external_reviewer_disposition",
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "campaign_id": admission.CAMPAIGN_ID,
        "claim_id": CLAIM_ID,
        "evidence_class": EVIDENCE_CLASS,
        "qualification_scope": QUALIFICATION_SCOPE,
        "status": "claim_candidate_pending_external_review",
        "claim_candidate": True,
        "final_claim_closure": False,
        "claim_closure_authorized": False,
        "frontier_execution_authorized": False,
        "policy_mutation_authorized": False,
        "campaign_matrix": aggregate["matrix"],
        "admitted_bundle_identities": aggregate["bundles"],
        "amr_fine_pair_coverage": aggregate["pairs"],
        "restart_parity": aggregate["restart"],
        "aggregate_analysis": {
            "record_type": numerical.RECORD_TYPE,
            "sha256": aggregate["sha256"],
            "status": "pending_external_review",
            "numerical_gates_passed": True,
            "reviewer_disposition": aggregate["aggregate_reviewer_disposition"],
        },
        "independent_recomputation": recompute,
        "external_reviewer_disposition": external_review,
        "claim_registry_boundary": {
            "required_gates": list(CLAIM_REQUIRED_GATES),
            "gates_closed_by_this_manifest": [],
            "claim_registry_promotion_authorized": False,
        },
        "remaining_gates": list(REMAINING_GATES),
    }
    validate_claim_candidate_evidence_manifest(manifest)
    return manifest


def validate_claim_candidate_evidence_manifest(value: object) -> dict[str, object]:
    """Validate the emitted compact manifest's non-authorizing boundary."""
    manifest = _object(
        value,
        {
            "schema_version",
            "record_type",
            "campaign_id",
            "claim_id",
            "evidence_class",
            "qualification_scope",
            "status",
            "claim_candidate",
            "final_claim_closure",
            "claim_closure_authorized",
            "frontier_execution_authorized",
            "policy_mutation_authorized",
            "campaign_matrix",
            "admitted_bundle_identities",
            "amr_fine_pair_coverage",
            "restart_parity",
            "aggregate_analysis",
            "independent_recomputation",
            "external_reviewer_disposition",
            "claim_registry_boundary",
            "remaining_gates",
        },
        label="claim_candidate_manifest",
    )
    _require(
        type(manifest["schema_version"]) is int and manifest["schema_version"] == SCHEMA_VERSION,
        "claim_candidate_manifest: schema version drifted",
    )
    _require(manifest["record_type"] == RECORD_TYPE, "claim_candidate_manifest: type drifted")
    _require(
        manifest["campaign_id"] == admission.CAMPAIGN_ID,
        "claim_candidate_manifest: campaign drifted",
    )
    _require(manifest["claim_id"] == CLAIM_ID, "claim_candidate_manifest: claim drifted")
    _require(
        manifest["evidence_class"] == EVIDENCE_CLASS,
        "claim_candidate_manifest: class drifted",
    )
    _require(
        manifest["qualification_scope"] == QUALIFICATION_SCOPE,
        "claim_candidate_manifest: scope drifted",
    )
    _require(
        manifest["status"] == "claim_candidate_pending_external_review"
        and manifest["claim_candidate"] is True,
        "claim_candidate_manifest: candidate status drifted",
    )
    for field in (
        "final_claim_closure",
        "claim_closure_authorized",
        "frontier_execution_authorized",
        "policy_mutation_authorized",
    ):
        _require(manifest[field] is False, f"claim_candidate_manifest: {field} must be false")
    _validate_campaign_matrix(manifest["campaign_matrix"])
    _require(
        len(_array(manifest["admitted_bundle_identities"], label="claim_candidate/bundles"))
        == numerical.EXPECTED_BASELINE_ATTEMPTS,
        "claim_candidate_manifest: admitted bundle count drifted",
    )
    _require(
        len(_array(manifest["amr_fine_pair_coverage"], label="claim_candidate/pairs"))
        == len(numerical.QUALIFYING_SEEDS),
        "claim_candidate_manifest: pair coverage count drifted",
    )
    restart = _object(
        manifest["restart_parity"],
        {
            "source_attempt_sha256",
            "restart_parity_binding_sha256",
            "restart_result_sha256",
            "screen_scope",
            "full_state_equivalence_claimed",
            "result",
        },
        label="claim_candidate/restart_parity",
    )
    _require(
        restart["screen_scope"] == numerical.RESTART_SCREEN_SCOPE
        and restart["full_state_equivalence_claimed"] is False
        and restart["result"] == "pass_deterministic_continuation_parity",
        "claim_candidate_manifest: restart boundary drifted",
    )
    aggregate = _object(
        manifest["aggregate_analysis"],
        {
            "record_type",
            "sha256",
            "status",
            "numerical_gates_passed",
            "reviewer_disposition",
        },
        label="claim_candidate/aggregate_analysis",
    )
    _require(
        aggregate["record_type"] == numerical.RECORD_TYPE
        and aggregate["status"] == "pending_external_review"
        and aggregate["numerical_gates_passed"] is True,
        "claim_candidate_manifest: aggregate boundary drifted",
    )
    _sha256(aggregate["sha256"], label="claim_candidate/aggregate_analysis/sha256")
    _reviewer_disposition(
        aggregate["reviewer_disposition"],
        required_status="pending_external_review",
        label="claim_candidate/aggregate_analysis/reviewer_disposition",
    )
    recompute = _object(
        manifest["independent_recomputation"],
        {
            "sha256",
            "status",
            "independent_script",
            "independent_environment_lock",
            "metric_comparison_row_count",
            "reviewer_identity",
            "reviewer_disposition",
        },
        label="claim_candidate/independent_recomputation",
    )
    _sha256(recompute["sha256"], label="claim_candidate/independent_recomputation/sha256")
    _require(
        recompute["status"] == "pass_independent_recompute_reviewed",
        "claim_candidate_manifest: independent recompute status drifted",
    )
    _binding(recompute["independent_script"], label="claim_candidate/independent_script")
    _binding(
        recompute["independent_environment_lock"],
        label="claim_candidate/independent_environment_lock",
    )
    _positive_integer(
        recompute["metric_comparison_row_count"],
        label="claim_candidate/metric_comparison_row_count",
    )
    reviewer_identity = _text(
        recompute["reviewer_identity"], label="claim_candidate/reviewer_identity"
    )
    recompute_disposition = _reviewer_disposition(
        recompute["reviewer_disposition"],
        required_status="accepted",
        label="claim_candidate/independent_recomputation/reviewer_disposition",
    )
    _require(
        recompute_disposition["reviewer"] == reviewer_identity,
        "claim_candidate_manifest: independent reviewer binding drifted",
    )
    _reviewer_disposition(
        manifest["external_reviewer_disposition"],
        required_status="pending_external_review",
        label="claim_candidate/external_reviewer_disposition",
    )
    boundary = _object(
        manifest["claim_registry_boundary"],
        {
            "required_gates",
            "gates_closed_by_this_manifest",
            "claim_registry_promotion_authorized",
        },
        label="claim_candidate/claim_registry_boundary",
    )
    _require_strict_equal(
        boundary["required_gates"],
        list(CLAIM_REQUIRED_GATES),
        label="claim_candidate/claim_registry_boundary/required_gates",
    )
    _require(
        boundary["gates_closed_by_this_manifest"] == []
        and boundary["claim_registry_promotion_authorized"] is False,
        "claim_candidate_manifest: claim registry boundary drifted",
    )
    _require_strict_equal(
        manifest["remaining_gates"],
        list(REMAINING_GATES),
        label="claim_candidate/remaining_gates",
    )
    _canonical_json_bytes(manifest)
    return manifest


def _decode_json(payload: bytes, *, label: str) -> object:
    def reject_duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result = {}
        for key, value in pairs:
            _require(key not in result, f"{label}: duplicate JSON key {key!r}")
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise FinalEvidenceError(f"{label}: forbidden JSON constant {value}")

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise FinalEvidenceError(f"{label}: JSON is not UTF-8") from error
    try:
        return json.loads(
            text,
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except json.JSONDecodeError as error:
        raise FinalEvidenceError(f"{label}: malformed JSON") from error


def _stable_file_identity(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _read_stable_json(path: Path, *, label: str) -> object:
    _require(path.is_absolute(), f"{label}: input path must be absolute")
    descriptor = None
    try:
        descriptor = os.open(path, _READ_FLAGS)
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode) and before.st_nlink == 1,
            f"{label}: input must be one regular file",
        )
        _require(before.st_size <= MAX_JSON_BYTES, f"{label}: input exceeds size limit")
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        _require(
            _stable_file_identity(before) == _stable_file_identity(after),
            f"{label}: input changed while reading",
        )
    except OSError as error:
        raise FinalEvidenceError(f"{label}: cannot read input: {error}") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)
    return _decode_json(bytes(payload), label=label)


def _write_exclusive(path: Path, payload: bytes) -> None:
    _require(path.is_absolute(), "output path must be absolute")
    descriptor = None
    try:
        descriptor = os.open(
            path,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
            0o600,
        )
        offset = 0
        while offset < len(payload):
            written = os.write(descriptor, payload[offset:])
            _require(written > 0, "short write while emitting claim-candidate manifest")
            offset += written
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
    except OSError as error:
        raise FinalEvidenceError(f"cannot emit claim-candidate manifest: {error}") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--admitted-bundle-identities", type=Path, required=True)
    parser.add_argument("--numerical-aggregate", type=Path, required=True)
    parser.add_argument("--independent-recompute", type=Path, required=True)
    parser.add_argument("--external-reviewer-disposition", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    manifest = build_claim_candidate_evidence_manifest(
        admitted_bundle_identities=_read_stable_json(
            args.admitted_bundle_identities,
            label="admitted bundle identities",
        ),
        numerical_aggregate=_read_stable_json(
            args.numerical_aggregate,
            label="numerical aggregate",
        ),
        independent_recompute=_read_stable_json(
            args.independent_recompute,
            label="independent recompute",
        ),
        external_reviewer_disposition=_read_stable_json(
            args.external_reviewer_disposition,
            label="external reviewer disposition",
        ),
    )
    _write_exclusive(args.output, _canonical_json_bytes(manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "ATTEMPT_PRIMARY_OBSERVABLES",
    "CLAIM_ID",
    "FinalEvidenceError",
    "INDEPENDENT_RECOMPUTE_RECORD_TYPE",
    "PAIR_OBSERVABLE_PREFIX",
    "RECORD_TYPE",
    "RESTART_OBSERVABLE",
    "build_claim_candidate_evidence_manifest",
    "validate_claim_candidate_evidence_manifest",
]
