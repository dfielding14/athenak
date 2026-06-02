#!/usr/bin/env python3
"""Fail-closed source-local Q-011 Section 5.4 numerical aggregation.

This module consumes already-retained admission and reducer records.  It does
not discover artifacts, execute AthenaK, publish evidence, automate human
review, or close a manuscript claim.
"""

from __future__ import annotations

import copy
import math
import re
from typing import Any, Mapping, Sequence

if __package__:
    from . import analyze_q011_section54_campaign as admission
    from . import q011_section54_artifacts as artifacts
    from . import q011_section54_particles as particles
    from . import q011_section54_restart as restart
    from . import q011_section54_spatial as spatial
else:
    import analyze_q011_section54_campaign as admission
    import q011_section54_artifacts as artifacts
    import q011_section54_particles as particles
    import q011_section54_restart as restart
    import q011_section54_spatial as spatial


SCHEMA_VERSION = 1
RECORD_TYPE = "q011_section54_source_local_numerical_qualification"
ATTEMPT_RECORD_TYPE = "q011_section54_source_local_numerical_baseline_attempt"
QUALIFICATION_SCOPE = (
    "bounded_source_local_numerical_aggregation_only_external_review_required_"
    "no_claim_closure"
)
PHYSICAL_MODE = "paper_mhd_pic_vl2_tsc"
GRID_VARIANTS = (
    "coarse_uniform_dx12",
    "three_level_amr_root_dx12_finest_dx3",
    "fine_uniform_dx3",
)
AMR_VARIANT = "three_level_amr_root_dx12_finest_dx3"
FINE_VARIANT = "fine_uniform_dx3"
QUALIFYING_SEEDS = (
    23050101,
    23050102,
    23050103,
    23050104,
    23050105,
    23050106,
    23050107,
    23050108,
)
EXPECTED_BASELINE_ATTEMPTS = 24
EXPECTED_MATRIX_CELLS = tuple(
    (variant, seed) for variant in GRID_VARIANTS for seed in QUALIFYING_SEEDS
)
PAIRED_SEED_RULE = (
    "Use the same qualifying seed for coarse-uniform, AMR and fine-uniform variants."
)
PAIRED_RESIDUAL_THRESHOLDS = (
    ("shock_front_position_at_t500", 240.0, None, None),
    ("upstream_magnetic_amplification_at_t500", 0.35, None, None),
    ("rho_y_average_at_t500", None, 0.2, 0.3),
    ("bmag_y_average_at_t500", None, 0.25, 0.35),
    ("normalized_downstream_chi_f_chi_at_t500", None, 0.2, 0.3),
    ("normalized_downstream_chi_f_chi_at_t1200", None, 0.2, 0.3),
)
_SHA256 = re.compile(r"[0-9a-f]{64}")
_ATTEMPT_ID = re.compile(r"[a-z0-9][a-z0-9._-]{0,127}")


class NumericalQualificationError(ValueError):
    """Raised when source-local numerical aggregation cannot proceed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise NumericalQualificationError(message)


def _object(value: object, keys: set[str], *, label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label}: expected object")
    _require(set(value) == keys, f"{label}: keys drifted")
    return value


def _required_object(
    value: object, keys: set[str], *, label: str
) -> dict[str, Any]:
    _require(type(value) is dict, f"{label}: expected object")
    _require(keys <= set(value), f"{label}: required keys drifted")
    return value


def _list(value: object, *, label: str) -> list[Any]:
    _require(type(value) is list, f"{label}: expected list")
    return value


def _sha256(value: object, *, label: str) -> str:
    _require(
        type(value) is str and _SHA256.fullmatch(value) is not None,
        f"{label}: malformed SHA-256",
    )
    return value


def _finite_float(value: object, *, label: str, minimum: float = 0.0) -> float:
    _require(
        type(value) is float and math.isfinite(value) and value >= minimum,
        f"{label}: expected finite float >= {minimum}",
    )
    return value


def _canonical_sha256(value: object) -> str:
    try:
        return artifacts.sha256_bytes(artifacts.canonical_json_bytes(value))
    except artifacts.DerivedArtifactError as error:
        raise NumericalQualificationError("record is not canonical finite JSON") from error


def bind_canonical_attempt(attempt: Mapping[str, object]) -> dict[str, object]:
    """Bind one caller-supplied baseline attempt without creating evidence."""
    _require(isinstance(attempt, Mapping), "attempt binding requires a mapping")
    retained = copy.deepcopy(dict(attempt))
    return {
        "attempt_sha256": _canonical_sha256(retained),
        "attempt": retained,
    }


def _expected_amr_pairing(identity: Mapping[str, object]) -> dict[str, object]:
    variant = identity["variant"]
    seed = identity["seed"]
    counterpart = None
    if variant == AMR_VARIANT:
        counterpart = {"variant": FINE_VARIANT, "seed": seed}
    elif variant == FINE_VARIANT:
        counterpart = {"variant": AMR_VARIANT, "seed": seed}
    return {
        "paired_seed_rule": PAIRED_SEED_RULE,
        "pair_key": {"seed": seed},
        "amr_fine_uniform_counterpart": counterpart,
        "comparison_status": (
            "schema_wired_not_evaluated_by_artifact_admission_slice"
        ),
    }


def _validate_identity(value: object, *, label: str) -> dict[str, object]:
    identity = _object(
        value,
        {"variant", "seed", "physical_mode", "attempt_id"},
        label=label,
    )
    _require(identity["variant"] in GRID_VARIANTS, f"{label}: variant is not canonical")
    _require(
        type(identity["seed"]) is int and identity["seed"] in QUALIFYING_SEEDS,
        f"{label}: seed is not canonical",
    )
    _require(identity["physical_mode"] == PHYSICAL_MODE, f"{label}: physical mode drifted")
    _require(
        type(identity["attempt_id"]) is str
        and _ATTEMPT_ID.fullmatch(identity["attempt_id"]) is not None,
        f"{label}: attempt ID is not canonical",
    )
    return identity


def _validate_admission(
    value: object,
    *,
    identity: Mapping[str, object],
    raw_inventory_sha256: str,
    label: str,
) -> None:
    result = _object(
        value,
        {
            "schema_version",
            "record_type",
            "campaign_id",
            "qualification_scope",
            "admitted_for_follow_on_numerical_qualification",
            "final_claim_closure",
            "status",
            "failure_reasons",
            "admission",
        },
        label=label,
    )
    _require(
        type(result["schema_version"]) is int and result["schema_version"] == 1,
        f"{label}: schema version drifted",
    )
    _require(result["record_type"] == admission.RESULT_RECORD_TYPE, f"{label}: record type drifted")
    _require(result["campaign_id"] == admission.CAMPAIGN_ID, f"{label}: campaign ID drifted")
    _require(
        result["qualification_scope"] == admission.QUALIFICATION_SCOPE,
        f"{label}: qualification scope drifted",
    )
    _require(
        result["admitted_for_follow_on_numerical_qualification"] is True,
        f"{label}: attempt was not admitted",
    )
    _require(result["final_claim_closure"] is False, f"{label}: claim closure must remain false")
    _require(
        result["status"] == "admitted_for_follow_on_numerical_qualification",
        f"{label}: admission status drifted",
    )
    _require(result["failure_reasons"] == [], f"{label}: admitted result has failure reasons")
    admitted = _required_object(
        result["admission"],
        {
            "run_identity",
            "preregistration_binding",
            "immutable_tree",
            "amr_pairing",
            "numerical_qualification_status",
        },
        label=f"{label}/admission",
    )
    _require(admitted["run_identity"] == identity, f"{label}: run identity binding drifted")
    immutable_tree = _required_object(
        admitted["immutable_tree"], {"inventory_sha256"}, label=f"{label}/immutable_tree"
    )
    _require(
        immutable_tree["inventory_sha256"] == raw_inventory_sha256,
        f"{label}: raw inventory binding drifted",
    )
    preregistration = _required_object(
        admitted["preregistration_binding"],
        {"sha256", "expected_sha256"},
        label=f"{label}/preregistration_binding",
    )
    _require(
        preregistration["sha256"] == admission.EXPECTED_PREREGISTRATION_SHA256
        and preregistration["expected_sha256"] == admission.EXPECTED_PREREGISTRATION_SHA256,
        f"{label}: preregistration binding drifted",
    )
    _require(
        admitted["amr_pairing"] == _expected_amr_pairing(identity),
        f"{label}: admission AMR pairing drifted",
    )
    _require(
        admitted["numerical_qualification_status"]
        == "not_evaluated_by_artifact_admission_slice",
        f"{label}: admission slice overclaimed numerical qualification",
    )


def _validate_overflow_gate(spectrum: object, *, label: str) -> bool:
    record = _required_object(
        spectrum,
        {
            "f_chi",
            "overflow_macro_weight",
            "total_post_filter_macro_weight",
            "overflow_macro_weight_fraction",
            "overflow_gate",
        },
        label=label,
    )
    overflow = _finite_float(record["overflow_macro_weight"], label=f"{label}/overflow")
    total = _finite_float(
        record["total_post_filter_macro_weight"], label=f"{label}/total"
    )
    _require(total > 0.0, f"{label}: total admitted macro weight must be positive")
    fraction = _finite_float(
        record["overflow_macro_weight_fraction"], label=f"{label}/overflow fraction"
    )
    _require(fraction == overflow / total, f"{label}: overflow fraction is inconsistent")
    gate = _object(
        record["overflow_gate"],
        {"maximum_macro_weight_fraction", "passed"},
        label=f"{label}/overflow_gate",
    )
    _require(
        gate["maximum_macro_weight_fraction"]
        == particles.MAX_OVERFLOW_MACRO_WEIGHT_FRACTION,
        f"{label}: overflow threshold drifted",
    )
    _require(type(gate["passed"]) is bool, f"{label}: overflow gate result must be boolean")
    _require(
        gate["passed"] == (fraction <= particles.MAX_OVERFLOW_MACRO_WEIGHT_FRACTION),
        f"{label}: overflow gate result is inconsistent",
    )
    return gate["passed"]


def _validate_particle_reduction(
    value: object, *, expected_time: float, require_late_slope: bool, label: str
) -> list[bool]:
    record = _required_object(
        value,
        {"schema_version", "record_type", "snapshot_time_omega0_inverse", "weighted_spectrum"},
        label=label,
    )
    _require(
        type(record["schema_version"]) is int and record["schema_version"] == particles.SCHEMA_VERSION,
        f"{label}: schema version drifted",
    )
    _require(
        record["record_type"] == "q011_section54_particle_snapshot_reduction",
        f"{label}: record type drifted",
    )
    _require(record["snapshot_time_omega0_inverse"] == expected_time, f"{label}: time drifted")
    try:
        particles.canonical_record_bytes(record)
    except particles.ParticleReducerError as error:
        raise NumericalQualificationError(f"{label}: particle record is invalid") from error
    gates = [_validate_overflow_gate(record["weighted_spectrum"], label=f"{label}/weighted_spectrum")]
    if require_late_slope:
        _require("late_slope" in record, f"{label}: late slope record is missing")
        try:
            recomputed = particles.late_slope_record(record["weighted_spectrum"]["f_chi"])
        except particles.ParticleReducerError as error:
            raise NumericalQualificationError(f"{label}: late slope record is invalid") from error
        _require(record["late_slope"] == recomputed, f"{label}: late slope record drifted")
        gates.append(recomputed["slope_gate_passed"])
    else:
        _require("late_slope" not in record, f"{label}: unexpected early late-slope record")
    return gates


def _validate_spatial_reduction(value: object, *, label: str) -> bool:
    record = _required_object(
        value,
        {"schema_version", "record_type", "time_omega0_inverse", "upstream_b_amplification"},
        label=label,
    )
    _require(
        type(record["schema_version"]) is int and record["schema_version"] == 1,
        f"{label}: schema version drifted",
    )
    _require(
        record["record_type"] == "q011_section54_t500_spatial_reduction",
        f"{label}: record type drifted",
    )
    _require(record["time_omega0_inverse"] == spatial.T500_OMEGA0_INVERSE, f"{label}: time drifted")
    amplification = _object(
        record["upstream_b_amplification"],
        {
            "time_omega0_inverse",
            "x_ideal_c_over_omega_pi",
            "upstream_window_c_over_omega_pi",
            "selected_cell_count",
            "selected_area",
            "mean_magnetic_magnitude",
            "reference_b0",
            "amplification_over_b0",
            "acceptance_range",
            "passes_gate",
        },
        label=f"{label}/upstream_b_amplification",
    )
    try:
        parsed = spatial.UpstreamBAmplificationRecord(
            time_omega0_inverse=amplification["time_omega0_inverse"],
            x_ideal_c_over_omega_pi=amplification["x_ideal_c_over_omega_pi"],
            upstream_window_c_over_omega_pi=tuple(amplification["upstream_window_c_over_omega_pi"]),
            selected_cell_count=amplification["selected_cell_count"],
            selected_area=amplification["selected_area"],
            mean_magnetic_magnitude=amplification["mean_magnetic_magnitude"],
            reference_b0=amplification["reference_b0"],
            amplification_over_b0=amplification["amplification_over_b0"],
            acceptance_range=tuple(amplification["acceptance_range"]),
            passes_gate=amplification["passes_gate"],
        )
        normalized = spatial.upstream_b_amplification_record(parsed)
    except (TypeError, spatial.AnalysisError) as error:
        raise NumericalQualificationError(f"{label}: spatial amplification is invalid") from error
    _require(amplification == normalized, f"{label}: spatial amplification record drifted")
    return parsed.passes_gate


def _validate_attempt_wrapper(value: object, *, index: int) -> dict[str, object]:
    label = f"attempts[{index}]"
    wrapper = _object(value, {"attempt_sha256", "attempt"}, label=label)
    expected_sha256 = _sha256(wrapper["attempt_sha256"], label=f"{label}/attempt_sha256")
    _require(
        _canonical_sha256(wrapper["attempt"]) == expected_sha256,
        f"{label}: canonical attempt SHA-256 drifted",
    )
    attempt = _object(
        wrapper["attempt"],
        {
            "schema_version",
            "record_type",
            "run_identity",
            "raw_inventory_sha256",
            "admission_result",
            "particle_reductions",
            "spatial_reduction",
        },
        label=f"{label}/attempt",
    )
    _require(
        type(attempt["schema_version"]) is int and attempt["schema_version"] == SCHEMA_VERSION,
        f"{label}: schema version drifted",
    )
    _require(attempt["record_type"] == ATTEMPT_RECORD_TYPE, f"{label}: record type drifted")
    identity = _validate_identity(attempt["run_identity"], label=f"{label}/run_identity")
    inventory = _sha256(attempt["raw_inventory_sha256"], label=f"{label}/raw_inventory_sha256")
    _validate_admission(
        attempt["admission_result"],
        identity=identity,
        raw_inventory_sha256=inventory,
        label=f"{label}/admission_result",
    )
    reductions = _object(
        attempt["particle_reductions"], {"t500", "t1200"}, label=f"{label}/particle_reductions"
    )
    gates = _validate_particle_reduction(
        reductions["t500"], expected_time=500.0, require_late_slope=False, label=f"{label}/t500"
    )
    gates.extend(
        _validate_particle_reduction(
            reductions["t1200"],
            expected_time=particles.LATE_SLOPE_SNAPSHOT_TIME,
            require_late_slope=True,
            label=f"{label}/t1200",
        )
    )
    gates.append(_validate_spatial_reduction(attempt["spatial_reduction"], label=f"{label}/spatial"))
    return {
        "attempt_sha256": expected_sha256,
        "identity": identity,
        "raw_inventory_sha256": inventory,
        "gates_passed": all(gates),
    }


def _validate_residual_metric(
    value: object,
    *,
    expected: tuple[str, float | None, float | None, float | None],
    label: str,
) -> tuple[dict[str, object], bool]:
    metric = _object(
        value,
        {
            "observable",
            "maximum_absolute_difference",
            "relative_mean_absolute",
            "relative_root_mean_square",
        },
        label=label,
    )
    observable, maximum_bound, mean_bound, rms_bound = expected
    _require(metric["observable"] == observable, f"{label}: observable order drifted")
    maximum = _finite_float(metric["maximum_absolute_difference"], label=f"{label}/maximum")
    passed = maximum_bound is None or maximum <= maximum_bound
    for field, bound in (
        ("relative_mean_absolute", mean_bound),
        ("relative_root_mean_square", rms_bound),
    ):
        if bound is None:
            _require(metric[field] is None, f"{label}/{field}: expected null")
        else:
            parsed = _finite_float(metric[field], label=f"{label}/{field}")
            passed = passed and parsed <= bound
    return metric, passed


def _validate_pairs(
    value: object, *, by_sha256: Mapping[str, Mapping[str, object]]
) -> list[dict[str, object]]:
    pairs = _list(value, label="paired_amr_fine_results")
    _require(len(pairs) == len(QUALIFYING_SEEDS), "paired AMR/fine result count drifted")
    results = []
    for index, (pair_value, seed) in enumerate(zip(pairs, QUALIFYING_SEEDS)):
        label = f"paired_amr_fine_results[{index}]"
        pair = _object(
            pair_value,
            {"seed", "amr_attempt_sha256", "fine_uniform_attempt_sha256", "residuals"},
            label=label,
        )
        _require(type(pair["seed"]) is int and pair["seed"] == seed, f"{label}: seed order drifted")
        amr_sha256 = _sha256(pair["amr_attempt_sha256"], label=f"{label}/amr_attempt_sha256")
        fine_sha256 = _sha256(
            pair["fine_uniform_attempt_sha256"], label=f"{label}/fine_uniform_attempt_sha256"
        )
        _require(amr_sha256 in by_sha256 and fine_sha256 in by_sha256, f"{label}: attempt binding is unknown")
        _require(
            (by_sha256[amr_sha256]["identity"]["variant"], by_sha256[amr_sha256]["identity"]["seed"])
            == (AMR_VARIANT, seed),
            f"{label}: AMR attempt is not the matched seed",
        )
        _require(
            (by_sha256[fine_sha256]["identity"]["variant"], by_sha256[fine_sha256]["identity"]["seed"])
            == (FINE_VARIANT, seed),
            f"{label}: fine-uniform attempt is not the matched seed",
        )
        metrics = _list(pair["residuals"], label=f"{label}/residuals")
        _require(
            len(metrics) == len(PAIRED_RESIDUAL_THRESHOLDS),
            f"{label}: residual metric count drifted",
        )
        parsed_metrics = []
        gates = []
        for metric_index, (metric, expected) in enumerate(
            zip(metrics, PAIRED_RESIDUAL_THRESHOLDS)
        ):
            parsed, passed = _validate_residual_metric(
                metric, expected=expected, label=f"{label}/residuals[{metric_index}]"
            )
            parsed_metrics.append(parsed)
            gates.append(passed)
        results.append(
            {
                "seed": seed,
                "amr_attempt_sha256": amr_sha256,
                "fine_uniform_attempt_sha256": fine_sha256,
                "residuals": parsed_metrics,
                "gates_passed": all(gates),
            }
        )
    return results


def _validate_restart_binding(
    value: object, *, by_sha256: Mapping[str, Mapping[str, object]]
) -> dict[str, object]:
    binding = _object(
        value,
        {"source_attempt_sha256", "uninterrupted", "continued"},
        label="restart_parity_binding",
    )
    source_sha256 = _sha256(
        binding["source_attempt_sha256"], label="restart_parity_binding/source_attempt_sha256"
    )
    _require(source_sha256 in by_sha256, "restart parity source attempt is unknown")
    _require(
        by_sha256[source_sha256]["identity"]["variant"] == AMR_VARIANT,
        "restart parity source attempt must be an AMR baseline",
    )
    binding_sha256 = _canonical_sha256(binding)
    try:
        result = restart.compare_deterministic_continuation_parity(
            binding["uninterrupted"], binding["continued"]
        )
    except restart.RestartPolicyError as error:
        raise NumericalQualificationError("restart continuation parity failed") from error
    _require(
        result.get("result") == "pass_deterministic_continuation_parity",
        "restart continuation comparator did not return a passing result",
    )
    return {
        "source_attempt_sha256": source_sha256,
        "restart_parity_binding_sha256": binding_sha256,
        "result": result,
    }


def qualify_numerical_aggregate(
    *,
    attempts: Sequence[object],
    ordered_raw_inventory_sha256_values: Sequence[object],
    paired_amr_fine_results: Sequence[object],
    restart_parity_binding: object,
    reviewer_disposition: object,
) -> dict[str, object]:
    """Aggregate one exact 3x8 admitted matrix without closing external review."""
    _require(type(attempts) is list, "attempts: expected list")
    _require(
        len(attempts) == EXPECTED_BASELINE_ATTEMPTS,
        "attempts: expected exactly 24 canonical baseline attempts",
    )
    parsed_attempts = [
        _validate_attempt_wrapper(value, index=index)
        for index, value in enumerate(attempts)
    ]
    attempt_sha256_values = [item["attempt_sha256"] for item in parsed_attempts]
    _require(
        len(set(attempt_sha256_values)) == EXPECTED_BASELINE_ATTEMPTS,
        "attempts: duplicate canonical baseline attempt",
    )
    attempt_ids = [item["identity"]["attempt_id"] for item in parsed_attempts]
    _require(len(set(attempt_ids)) == EXPECTED_BASELINE_ATTEMPTS, "attempts: duplicate attempt ID")
    cells = [
        (item["identity"]["variant"], item["identity"]["seed"])
        for item in parsed_attempts
    ]
    _require(cells == list(EXPECTED_MATRIX_CELLS), "attempts: canonical 3x8 matrix order drifted")

    inventories = _list(
        ordered_raw_inventory_sha256_values,
        label="ordered_raw_inventory_sha256_values",
    )
    _require(
        len(inventories) == EXPECTED_BASELINE_ATTEMPTS,
        "ordered raw inventory digest count drifted",
    )
    parsed_inventories = [
        _sha256(value, label=f"ordered_raw_inventory_sha256_values[{index}]")
        for index, value in enumerate(inventories)
    ]
    _require(
        parsed_inventories
        == [item["raw_inventory_sha256"] for item in parsed_attempts],
        "ordered raw inventory digest binding drifted",
    )

    by_sha256 = {item["attempt_sha256"]: item for item in parsed_attempts}
    pair_results = _validate_pairs(paired_amr_fine_results, by_sha256=by_sha256)
    restart_result = _validate_restart_binding(restart_parity_binding, by_sha256=by_sha256)
    try:
        disposition = artifacts.validate_reviewer_disposition(reviewer_disposition)
    except artifacts.DerivedArtifactError as error:
        raise NumericalQualificationError("reviewer disposition is invalid") from error
    _require(
        disposition["status"] == "pending_external_review",
        "numerical aggregation must remain pending external review",
    )
    numerical_gates_passed = (
        all(item["gates_passed"] for item in parsed_attempts)
        and all(item["gates_passed"] for item in pair_results)
    )
    result = {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "campaign_id": admission.CAMPAIGN_ID,
        "qualification_scope": QUALIFICATION_SCOPE,
        "final_claim_closure": False,
        "status": "pending_external_review",
        "numerical_gate_status": "passed" if numerical_gates_passed else "failed",
        "numerical_gates_passed": numerical_gates_passed,
        "campaign_matrix": {
            "physical_mode": PHYSICAL_MODE,
            "grid_variants": list(GRID_VARIANTS),
            "qualifying_seeds": list(QUALIFYING_SEEDS),
            "expected_baseline_attempts": EXPECTED_BASELINE_ATTEMPTS,
            "paired_seed_rule": PAIRED_SEED_RULE,
        },
        "attempt_count": len(parsed_attempts),
        "attempt_sha256_values": attempt_sha256_values,
        "ordered_raw_inventory_sha256_values": parsed_inventories,
        "attempt_gate_results": parsed_attempts,
        "paired_amr_fine_results": pair_results,
        "restart_parity": restart_result,
        "reviewer_disposition": disposition,
    }
    _canonical_sha256(result)
    return result


__all__ = [
    "AMR_VARIANT",
    "ATTEMPT_RECORD_TYPE",
    "EXPECTED_BASELINE_ATTEMPTS",
    "EXPECTED_MATRIX_CELLS",
    "FINE_VARIANT",
    "GRID_VARIANTS",
    "NumericalQualificationError",
    "PAIRED_RESIDUAL_THRESHOLDS",
    "PHYSICAL_MODE",
    "QUALIFYING_SEEDS",
    "QUALIFICATION_SCOPE",
    "RECORD_TYPE",
    "SCHEMA_VERSION",
    "bind_canonical_attempt",
    "qualify_numerical_aggregate",
]
