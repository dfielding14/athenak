#!/usr/bin/env python3
"""Adversarial checks for the Bell current-normalization supersession."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import re
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
RECORD_PATH = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q043_bell_current_normalization_supersession_2026-06-06.json"
)

LEGACY_RELATION = "PPC*deposit_qscale*species_charge*v=2*B0*C*k0"
VOLUME_AWARE_CLOSURE = (
    "PPC*deposit_qscale*species_charge*v/V_root_cell=2*B0*k0"
)
APPLICABLE_CAMPAIGNS = {
    "Q023-PAPER-BELL-LINEAR",
    "Q023-PROD-BELL-NONLINEAR",
    "Q019-HR-FIXED-CURRENT-LIKE-NOHALL",
    "Q029-HALL-BELL-LINEAR",
}
REQUIRED_OUTPUTS = {"prtcl_rho", "prtcl_jx", "prtcl_jy", "prtcl_jz"}
LEGACY_GENERATORS = {"q023_paper_bell_linear", "q029_hall_bell_linear"}
EXPECTED_AFFECTED_PATHS = {
    "src/pgen/tests/q023_paper_bell_linear.cpp",
    "src/pgen/tests/q029_hall_bell_linear.cpp",
    "inputs/tests/pic_q023_paper_bell_linear_1d_candidate.athinput",
    "inputs/tests/pic_q023_paper_bell_linear_1d_candidate_vl2_tsc.athinput",
    "inputs/tests/pic_q023_paper_bell_linear_2d_candidate.athinput",
    "inputs/tests/pic_q023_paper_bell_linear_2d_candidate_vl2_tsc.athinput",
    "inputs/tests/pic_q023_paper_bell_linear_3d_candidate.athinput",
    "inputs/tests/pic_q023_paper_bell_linear_3d_candidate_vl2_tsc.athinput",
    "inputs/tests/pic_q023_prod_bell_nonlinear_foundation_vl2_tsc.athinput",
    "inputs/tests/pic_q029_hall_bell_linear_1d_candidate.athinput",
    "inputs/tests/pic_q029_hall_bell_linear_2d_candidate.athinput",
    "inputs/tests/pic_q029_hall_bell_linear_3d_candidate.athinput",
    "tst/publication/analyze_q023_paper_bell_linear.py",
    "tst/publication/analyze_q023_paper_bell_linear_vl2_tsc.py",
    "tst/publication/materialize_q023_paper_bell_linear_variants.py",
    "tst/publication/analyze_q023_prod_bell_nonlinear_foundation.py",
    "tst/publication/analyze_q029_hall_bell_linear.py",
    "tst/publication/analyze_q029_hall_bell_linear_raw.py",
    "tst/publication/q023_paper_bell_linear_host_harness.cpp",
    "tst/publication/q029_hall_bell_linear_host_harness.cpp",
    "tst/publication/test_analyze_q023_paper_bell_linear.py",
    "tst/publication/test_materialize_q023_paper_bell_linear_variants.py",
    "tst/publication/test_q023_paper_bell_linear_host_harness.py",
    "tst/publication/test_analyze_q023_prod_bell_nonlinear_foundation.py",
    "tst/publication/test_analyze_q029_hall_bell_linear.py",
    "tst/publication/test_analyze_q029_hall_bell_linear_raw.py",
    "tst/publication/test_q029_hall_bell_linear_host_harness.py",
    (
        "tst/publication/readiness/"
        "q023_paper_bell_linear_source_local_implementation_2026-05-30.json"
    ),
    (
        "tst/publication/readiness/"
        "q023_paper_bell_linear_source_local_implementation_successor_2026-05-31.json"
    ),
    (
        "tst/publication/readiness/"
        "q023_paper_bell_linear_source_local_implementation_successor_v2_2026-05-31.json"
    ),
    (
        "tst/publication/readiness/"
        "q023_paper_bell_linear_source_local_implementation_successor_v3_2026-06-01.json"
    ),
    (
        "tst/publication/readiness/"
        "q023_paper_bell_linear_materialized_variants_local_2026-05-30.json"
    ),
    (
        "tst/publication/readiness/"
        "q023_prod_bell_nonlinear_source_local_foundation_2026-06-06.json"
    ),
    "tst/publication/readiness/q023_campaign_drafts_2026-05-30.json",
    (
        "tst/publication/readiness/"
        "q019_nonlinear_bell_saturation_campaign_design_2026-06-06.json"
    ),
    (
        "tst/publication/readiness/"
        "q019_nonlinear_bell_saturation_campaign_design_2026-06-06.md"
    ),
    (
        "tst/publication/readiness/"
        "q029_experimental_hall_normalization_derivation_2026-05-30.md"
    ),
    (
        "tst/publication/readiness/"
        "q029_hall_bell_linear_source_local_preparation_2026-05-30.json"
    ),
    (
        "tst/publication/readiness/"
        "q029_hall_bell_linear_raw_extractor_source_local_2026-05-30.json"
    ),
    (
        "tst/publication/readiness/"
        "q029_hall_bell_q022_prerequisite_successor_2026-05-30.json"
    ),
}


class ContractError(ValueError):
    """The source-local supersession record was weakened or drifted."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _affected_artifacts(record: dict[str, object]) -> list[dict[str, str]]:
    groups = record["affected_artifacts"]
    _require(set(groups) == {"Q023", "Q019", "Q029"}, "affected groups drifted")
    artifacts = []
    for group_name, group in groups.items():
        _require("preserved" in group["status"], f"{group_name} is not preserved")
        _require("invalidated" in group["status"], f"{group_name} is not invalidated")
        artifacts.extend(group["artifacts"])
    return artifacts


def validate_supersession(record: dict[str, object]) -> None:
    _require(record["record_type"] == "bell_current_normalization_supersession",
             "record type drifted")
    _require(record["schema_version"] == 1, "schema version drifted")
    _require(record["qualification_effect"] ==
             "source_local_invalidation_only_no_execution_authorization_"
             "no_policy_mutation_no_claim_closure",
             "qualification boundary drifted")

    authority = record["authority"]
    _require(authority and not any(authority.values()), "authority must remain false")

    normalization = record["normalization_supersession"]
    _require(normalization["legacy_relation"] == LEGACY_RELATION,
             "legacy relation drifted")
    _require(normalization["legacy_relation_status"] ==
             "invalid_for_physical_bell_execution_or_qualification",
             "legacy relation was not invalidated")
    _require(normalization["exact_volume_aware_closure"] == VOLUME_AWARE_CLOSURE,
             "volume-aware closure drifted")
    _require(normalization["deposited_moment_semantics"] == "prtcl_j_is_j_CR_over_c",
             "deposited current semantics drifted")
    _require(
        normalization["physical_current_mapping"] ==
        "external_normalization_required_artificial_CR_light_speed_is_not_a_"
        "source_level_conversion_factor",
        "physical-current mapping was overstated",
    )
    _require(normalization["legacy_runtime_overcurrent_factor"] == "C/V_root_cell",
             "legacy overcurrent factor drifted")
    _require({item["id"] for item in normalization["independent_defects"]} ==
             {"extra_artificial_light_speed", "missing_root_cell_volume"},
             "independent defect inventory drifted")

    prohibition = record["source_local_execution_and_qualification_prohibition"]
    _require(set(prohibition["applies_to"]) == APPLICABLE_CAMPAIGNS,
             "campaign prohibition scope drifted")
    for field in (
        "legacy_relation_physical_execution",
        "legacy_relation_scientific_qualification",
        "legacy_relation_publication_evidence",
    ):
        _require(prohibition[field] == "prohibited", f"{field} was weakened")

    consequences = record["campaign_consequences"]
    _require(set(consequences) == {"Q023", "Q019", "Q029"},
             "campaign consequence inventory drifted")
    for campaign, consequence in consequences.items():
        _require("prohibited" in consequence["physical_status"],
                 f"{campaign} physical status was weakened")

    interim = record["interim_cycle_zero_oracle_consequence"]
    _require(
        interim["status"] ==
        "preserved_non_authorizing_chronology_invalid_as_raw_deposited_current_"
        "evidence",
        "interim cycle-zero oracle was not invalidated",
    )
    _require(
        "precedes_particle_moment_deposition" in interim["reason"],
        "interim cycle-zero oracle failure mode drifted",
    )
    for path, digest in interim["artifact_bindings"].items():
        _require(_sha256(REPO_ROOT / path) == digest,
                 f"interim cycle-zero oracle artifact drifted: {path}")

    successor = record["required_corrected_successor"]
    _require(
        successor["implementation_status"] ==
        "source_local_implementation_and_registration_ready_"
        "runtime_oracle_observation_pending",
        "successor source-local implementation status drifted",
    )
    _require(successor["campaign_id"] == "Q043-BELL-CURRENT-VOLUME-AWARE",
             "successor campaign identity drifted")
    _require(successor["generator_name"] == "q043_bell_current_volume_aware",
             "successor generator identity drifted")
    _require(successor["readiness_record_prefix"] ==
             "q043_bell_current_volume_aware_successor",
             "successor readiness identity drifted")
    _require(successor["must_be_uniquely_named"], "unique successor no longer required")
    _require(successor["must_not_silently_modify_or_relabel_legacy_artifacts"],
             "historical preservation was weakened")
    _require(set(successor["legacy_generator_names_forbidden_for_successor"]) ==
             LEGACY_GENERATORS, "legacy generator exclusion drifted")
    _require(successor["generator_name"] not in LEGACY_GENERATORS,
             "successor reuses a legacy generator identity")

    oracle = record["required_output_level_deposited_current_oracle"]
    _require(
        oracle["status"] ==
        "source_local_deck_matrix_and_raw_output_analyzer_ready_"
        "runtime_observation_pending",
        "oracle source-local readiness status drifted",
    )
    _require(oracle["source_of_truth"] ==
             "raw_output_level_deposited_current_not_deck_arithmetic_or_"
             "generator_preflight",
             "oracle source of truth was weakened")
    _require(oracle["required_observation_cycle"] == 1,
             "oracle must observe deposited moments after one complete cycle")
    _require(set(oracle["required_outputs"]) == REQUIRED_OUTPUTS,
             "required deposited-current outputs drifted")
    _require(oracle["output_semantics"] == "deposited_j_CR_over_c",
             "oracle output semantics drifted")
    _require(oracle["parallel_closure"] ==
             "dot(volume_average(prtcl_j),b_hat_stream)=2*B0*k0",
             "oracle target drifted")
    _require(
        oracle["required_matrix"]["decompositions"]
        == ["serial", "MPI_x1", "MPI_x2", "MPI_x3", "MPI_multiaxis"],
        "decomposition oracle was weakened",
    )
    _require(
        "require_registered_raw_output_coverage_of_x1_x2_x3_and_multiaxis_"
        "MPI_moment_exchange_paths"
        in successor["required_source_contracts"],
        "multidirectional MPI oracle was weakened",
    )
    _require(
        "require_positive_integer_PPC_because_AthenaK_realizes_a_discrete_global_"
        "particle_count"
        in successor["required_source_contracts"],
        "discrete PPC contract was weakened",
    )
    _require(oracle["required_matrix"]["minimum_resolution_values_per_geometry"] >= 2,
             "resolution oracle was weakened")
    _require(oracle["required_matrix"]["minimum_PPC_values_per_geometry"] >= 2,
             "PPC oracle was weakened")
    _require(
        oracle["required_matrix"]["artificial_C_over_v_CR"] == [100, 1000, 10000],
        "artificial-C runtime oracle was weakened",
    )

    artifacts = _affected_artifacts(record)
    paths = [artifact["path"] for artifact in artifacts]
    _require(len(paths) == len(set(paths)), "affected artifacts contain duplicates")
    _require(set(paths) == EXPECTED_AFFECTED_PATHS, "affected artifact inventory drifted")
    for artifact in artifacts:
        _require(re.fullmatch(r"[0-9a-f]{64}", artifact["sha256"]) is not None,
                 f"invalid sha256 for {artifact['path']}")
        _require(_sha256(REPO_ROOT / artifact["path"]) == artifact["sha256"],
                 f"historical artifact bytes drifted: {artifact['path']}")

    generated = record["affected_artifacts"]["Q023"]["generated_variant_scope"]
    _require(generated["affected_materialized_deck_count"] == 405,
             "materialized variant invalidation scope drifted")
    _require(record["control_documents_consulted_not_modified"] ==
             ["tst/publication/PIC_PRODUCTION_READINESS_PLAN.md"],
             "control-document non-mutation boundary drifted")


class Q043BellCurrentNormalizationSupersessionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.record = json.loads(RECORD_PATH.read_text(encoding="utf-8"))

    def test_record_is_fail_closed_and_preserves_affected_artifacts(self) -> None:
        validate_supersession(self.record)

    def test_source_audit_evidence_still_exposes_the_legacy_defect(self) -> None:
        q023 = (REPO_ROOT / "src/pgen/tests/q023_paper_bell_linear.cpp").read_text()
        q029 = (REPO_ROOT / "src/pgen/tests/q029_hall_bell_linear.cpp").read_text()
        moments = (REPO_ROOT / "src/particles/particles_moments.cpp").read_text()
        self.assertIn("2.0*b_g*light_speed*k0", q023)
        self.assertIn("2.0*b_g*light_speed*k0", q029)
        self.assertIn("q_density = q_macro*inv_cell_vol", moments)
        self.assertIn("weighted_q_density*vx", moments)

    def test_authoritative_source_uses_charge_closure_and_mode_scoped_mass(self) -> None:
        source = (
            REPO_ROOT / "src/pgen/tests/q043_bell_current_volume_aware.cpp"
        ).read_text()
        self.assertIn(
            "ppc*qscale*species_charge*stream_speed/root_cell_volume", source
        )
        self.assertIn(
            "species_charge/species_mass must equal omega/b_g",
            source,
        )
        self.assertIn(
            "uniform_current_oracle accepts any positive species_mass", source
        )
        self.assertIn(
            "corrected_linear_eigenmode requires species_mass=1", source
        )
        self.assertNotIn("2.0*b_g*light_speed*k0", source)
        self.assertNotIn("ppc*qscale*(species_charge/species_mass)", source)

    def test_rejects_non_volume_aware_or_C_weighted_successor_closure(self) -> None:
        bad_closures = [
            "PPC*deposit_qscale*species_charge*v=2*B0*k0",
            "PPC*deposit_qscale*species_charge*v/V_root_cell=2*B0*C*k0",
            "PPC*deposit_qscale*(q/mc)*v/V_root_cell=2*B0*k0",
            LEGACY_RELATION,
        ]
        for bad_closure in bad_closures:
            with self.subTest(bad_closure=bad_closure):
                candidate = copy.deepcopy(self.record)
                candidate["normalization_supersession"]["exact_volume_aware_closure"] = (
                    bad_closure
                )
                with self.assertRaises(ContractError):
                    validate_supersession(candidate)

    def test_rejects_weakened_Q023_Q019_or_Q029_prohibition(self) -> None:
        for campaign in ("Q023", "Q019", "Q029"):
            with self.subTest(campaign=campaign):
                candidate = copy.deepcopy(self.record)
                candidate["campaign_consequences"][campaign]["physical_status"] = (
                    "physical_execution_allowed"
                )
                with self.assertRaises(ContractError):
                    validate_supersession(candidate)
        candidate = copy.deepcopy(self.record)
        candidate["source_local_execution_and_qualification_prohibition"][
            "applies_to"
        ].remove("Q029-HALL-BELL-LINEAR")
        with self.assertRaises(ContractError):
            validate_supersession(candidate)

    def test_rejects_legacy_identity_reuse_or_non_output_oracle(self) -> None:
        candidates = []
        reused_name = copy.deepcopy(self.record)
        reused_name["required_corrected_successor"]["generator_name"] = (
            "q023_paper_bell_linear"
        )
        candidates.append(reused_name)

        missing_output = copy.deepcopy(self.record)
        missing_output["required_output_level_deposited_current_oracle"][
            "required_outputs"
        ].remove("prtcl_jz")
        candidates.append(missing_output)

        arithmetic_only = copy.deepcopy(self.record)
        arithmetic_only["required_output_level_deposited_current_oracle"][
            "source_of_truth"
        ] = "deck_arithmetic_only"
        candidates.append(arithmetic_only)

        for candidate in candidates:
            with self.assertRaises(ContractError):
                validate_supersession(candidate)

    def test_rejects_artifact_omission_hash_drift_or_side_effect_claim(self) -> None:
        omitted = copy.deepcopy(self.record)
        omitted["affected_artifacts"]["Q019"]["artifacts"].pop()
        with self.assertRaises(ContractError):
            validate_supersession(omitted)

        drifted = copy.deepcopy(self.record)
        drifted["affected_artifacts"]["Q023"]["artifacts"][0]["sha256"] = "0" * 64
        with self.assertRaises(ContractError):
            validate_supersession(drifted)

        side_effect = copy.deepcopy(self.record)
        side_effect["authority"]["active_policy_modified"] = True
        with self.assertRaises(ContractError):
            validate_supersession(side_effect)


if __name__ == "__main__":
    unittest.main()
